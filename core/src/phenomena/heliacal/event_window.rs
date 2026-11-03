//! Symmetric visibility window around an already found heliacal event.
//!
//! This stencil is not the phenomena rise-window search: it has its own
//! symmetric bracket, scan count, optimum retention and crossing semantics.

use super::{
    HELIACAL_WINDOW_CROSS_STEPS, HELIACAL_WINDOW_GOLDEN_STEPS, HELIACAL_WINDOW_HALF_DAYS,
    HELIACAL_WINDOW_SCAN_STEPS, HELIACAL_WINDOW_TOL_DAYS,
};
use crate::errors::{Error, Result};

/// Visible optimum to outer edge; unlike the phenomena crossing, an
/// invisible optimum immediately returns the degenerate inner edge.
pub fn heliacal_window_crossing(t_in: f64, t_out: f64, margin_at: impl Fn(f64) -> f64) -> f64 {
    let f_in = margin_at(t_in);
    let f_out = margin_at(t_out);
    if f_in <= 0.0 {
        return t_in;
    }
    if f_out > 0.0 {
        return t_out;
    }
    let (mut inner, mut outer) = (t_in, t_out);
    for _ in 0..HELIACAL_WINDOW_CROSS_STEPS {
        let mid = 0.5 * (inner + outer);
        if margin_at(mid) > 0.0 {
            inner = mid;
        } else {
            outer = mid;
        }
        if (outer - inner).abs() < HELIACAL_WINDOW_TOL_DAYS {
            break;
        }
    }
    0.5 * (inner + outer)
}

/// Return `(start, optimum, end)` in UT Julian days from the symmetric
/// 145-sample margin scan and golden-section refinement. The scan optimum
/// survives unless the refined midpoint is better. Rejects nonfinite anchors;
/// a flat miss returns a degenerate window at the event epoch.
pub fn heliacal_visibility_window(
    jd_event: f64,
    margin_at: impl Fn(f64) -> f64,
) -> Result<(f64, f64, f64)> {
    if !jd_event.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_event });
    }
    let win = HELIACAL_WINDOW_HALF_DAYS;
    let n = HELIACAL_WINDOW_SCAN_STEPS;
    let mut jd_opt = jd_event;
    let mut best = margin_at(jd_event);
    for i in 0..=n {
        let t = jd_event - win + 2.0 * win * i as f64 / n as f64;
        let m = margin_at(t);
        if m > best {
            best = m;
            jd_opt = t;
        }
    }
    let step = 2.0 * win / n as f64;
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    let (mut lo, mut hi) = (jd_opt - step, jd_opt + step);
    let mut a = hi - (hi - lo) / phi;
    let mut b = lo + (hi - lo) / phi;
    let mut fa = margin_at(a);
    let mut fb = margin_at(b);
    for _ in 0..HELIACAL_WINDOW_GOLDEN_STEPS {
        if fa > fb {
            hi = b;
            b = a;
            fb = fa;
            a = hi - (hi - lo) / phi;
            fa = margin_at(a);
        } else {
            lo = a;
            fa = fb;
            a = b;
            b = lo + (hi - lo) / phi;
            fb = margin_at(b);
        }
        if hi - lo < HELIACAL_WINDOW_TOL_DAYS {
            break;
        }
    }
    // Preserve the scan maximum unless the refined midpoint is better.
    if margin_at(0.5 * (lo + hi)) > best {
        jd_opt = 0.5 * (lo + hi);
    }
    let jd_start = heliacal_window_crossing(jd_opt, jd_opt - win, &margin_at);
    let jd_end = heliacal_window_crossing(jd_opt, jd_opt + win, &margin_at);
    Ok((jd_start.min(jd_opt), jd_opt, jd_end.max(jd_opt)))
}
