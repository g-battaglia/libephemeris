//! Rise/set anchored phenomena windows, independent of event-window refinement.
//!
//! Solar and object horizon searches are consumed before margin sampling.
//! Missing horizon/window times retain TJD_INVALID, distinct from the zero
//! duration convention for unsupported acronychal targets.

use super::{
    HeliacalPhenoWindow, PHENO_WINDOW_GOLDEN_STEPS, PHENO_WINDOW_GRID_STEPS,
    PHENO_WINDOW_LOOKBACK_DAYS, PHENO_WINDOW_TOL_DAYS,
};
use crate::constants::{MERCURY, MOON, VENUS};
use crate::errors::{Error, Result};

/// This crossing assumes a visible optimum and reads the outer edge first.
/// Keep its 22-iteration bound and tolerance separate from event windows.
fn pheno_window_crossing(t_in: f64, t_out: f64, margin_at: impl Fn(f64) -> f64) -> f64 {
    if margin_at(t_out) > 0.0 {
        return t_out;
    }
    let (mut inner, mut outer) = (t_in, t_out);
    for _ in 0..PHENO_WINDOW_GOLDEN_STEPS {
        let mid = 0.5 * (inner + outer);
        if margin_at(mid) > 0.0 {
            inner = mid;
        } else {
            outer = mid;
        }
        if (outer - inner).abs() < PHENO_WINDOW_TOL_DAYS {
            break;
        }
    }
    0.5 * (inner + outer)
}

/// Compose UT1 horizon/visibility times for event codes 1–4.
///
/// Horizon callbacks receive a UT1 Julian-day search start and return the
/// selected disc-center crossing, or `None`. The caller binds rise for codes
/// 1/4 and set for 2/3; this function does not pass an event selector to them.
/// The Sun callback runs before the object callback, each starting the fixed
/// lookback before `jd_ut`. Margin sampling uses UT1 and a positive value for
/// visibility, normally limiting magnitude minus object magnitude.
///
/// The solar crossing selects the bounded visibility interval. Missing window
/// times and duration use TJD_INVALID; missing horizon pairs give zero lag.
/// Unsupported acronychal targets retain horizon words but zero duration, as
/// described on [`HeliacalPhenoWindow`]. Nonfinite `jd_ut` returns
/// `NonFiniteArgument`; unsupported event codes return `InvalidBody`.
/// Callback misses have no typed error channel and are not proof of complete
/// coverage. Supplied callback values are not independently range-validated.
#[allow(clippy::too_many_arguments)]
pub fn pheno_rise_window(
    jd_ut: f64,
    body_id: i32,
    is_star: bool,
    event_type: i32,
    named_for_window: bool,
    rise_s_at: impl Fn(f64) -> Option<f64>,
    rise_o_at: impl Fn(f64) -> Option<f64>,
    margin_at: impl Fn(f64) -> f64,
) -> Result<HeliacalPhenoWindow> {
    use crate::constants::{HELIACAL_RISING, MORNING_LAST, TJD_INVALID};
    if !jd_ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_ut });
    }
    if event_type != HELIACAL_RISING
        && event_type != crate::constants::HELIACAL_SETTING
        && event_type != crate::constants::EVENING_FIRST
        && event_type != MORNING_LAST
    {
        return Err(Error::InvalidBody {
            body_id: event_type,
            operation: "heliacal_pheno_ut",
        });
    }
    let start = jd_ut - PHENO_WINDOW_LOOKBACK_DAYS;
    let sun_hit = rise_s_at(start);
    let obj_hit = rise_o_at(start);
    let rise_s = sun_hit.unwrap_or(TJD_INVALID);
    let sun_found = sun_hit.is_some();
    let rise_o = obj_hit.unwrap_or(TJD_INVALID);
    let norise = obj_hit.is_none();
    let lag = if norise || !sun_found {
        0.0
    } else {
        rise_o - rise_s
    };
    let mut t_yallop = TJD_INVALID;
    if !is_star && body_id == MOON && !norise && sun_found {
        t_yallop = (rise_o * 4.0 + rise_s * 5.0) / 9.0;
    }
    let invalid = HeliacalPhenoWindow {
        rise_o,
        rise_s,
        lag,
        t_first: TJD_INVALID,
        t_best: TJD_INVALID,
        t_last: TJD_INVALID,
        tvis: TJD_INVALID,
        t_yallop,
    };
    // Acronychal visibility windows apply only to the Moon and inner planets.
    // Other targets retain the horizon words but use zero duration.
    if (event_type == crate::constants::EVENING_FIRST || event_type == MORNING_LAST)
        && (is_star || ![MOON, MERCURY, VENUS].contains(&body_id))
    {
        return Ok(HeliacalPhenoWindow {
            tvis: 0.0,
            ..invalid
        });
    }
    if !sun_found {
        return Ok(invalid);
    }
    if !named_for_window {
        return Ok(invalid);
    }
    let win = PHENO_WINDOW_LOOKBACK_DAYS;
    let (lo, hi) = if event_type == HELIACAL_RISING || event_type == MORNING_LAST {
        (rise_s - win, rise_s)
    } else {
        (rise_s, rise_s + win)
    };
    // Find a narrow positive spike before refining the locally unimodal margin.
    let n = PHENO_WINDOW_GRID_STEPS;
    let step = (hi - lo) / n as f64;
    let (mut best_i, mut best_m) = (0_usize, -1e99_f64);
    for i in 0..=n {
        let m = margin_at(lo + i as f64 * step);
        if m > best_m {
            best_i = i;
            best_m = m;
        }
    }
    let mut g_lo = lo + best_i.saturating_sub(1) as f64 * step;
    let mut g_hi = lo + best_i.saturating_add(1).min(n) as f64 * step;
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let mut a = g_hi - (g_hi - g_lo) / phi;
    let mut b = g_lo + (g_hi - g_lo) / phi;
    let mut fa = margin_at(a);
    let mut fb = margin_at(b);
    for _ in 0..PHENO_WINDOW_GOLDEN_STEPS {
        if fa > fb {
            g_hi = b;
            b = a;
            fb = fa;
            a = g_hi - (g_hi - g_lo) / phi;
            fa = margin_at(a);
        } else {
            g_lo = a;
            fa = fb;
            a = b;
            b = g_lo + (g_hi - g_lo) / phi;
            fb = margin_at(b);
        }
        if g_hi - g_lo < PHENO_WINDOW_TOL_DAYS {
            break;
        }
    }
    let t_peak = 0.5 * (g_lo + g_hi);
    if margin_at(t_peak) <= 0.0 {
        return Ok(invalid);
    }
    let t_best = t_peak;
    let mut t_first = pheno_window_crossing(t_peak, lo, &margin_at);
    let mut t_last = pheno_window_crossing(t_peak, hi, &margin_at);
    if !norise {
        if (event_type == HELIACAL_RISING || event_type == MORNING_LAST) && t_first != TJD_INVALID {
            t_first = t_first.max(rise_o);
        } else if event_type != HELIACAL_RISING
            && event_type != MORNING_LAST
            && t_last != TJD_INVALID
        {
            t_last = t_last.min(rise_o);
        }
    }
    let tvis = if t_first != TJD_INVALID && t_last != TJD_INVALID {
        t_last - t_first
    } else {
        TJD_INVALID
    };
    Ok(HeliacalPhenoWindow {
        rise_o,
        rise_s,
        lag,
        t_first,
        t_best,
        t_last,
        tvis,
        t_yallop,
    })
}
