//! Rise/set scanning and refinement over borrowed event callbacks.
use crate::errors::Error;
use crate::events::{self, RiseSetOutcome};

/// First Native rise/set height crossing within two days. `event_height`
/// returns disc/refraction/horizon-adjusted degrees; `get_ra_dec` supplies
/// apparent RA hours/Dec degrees for the declination-speed precheck. Positive
/// height means above the horizon. Rise uses h_prev<0<=h, set h_prev>=0>h.
/// The latitude, horizon and lift inputs are degrees and all acquired values
/// must be finite. The analytic extreme-altitude margin uses |dec2-dec|*24*2.
/// Grazing margin selects 10/20/60-minute steps; the scan acquires exactly
/// floor(2/step) samples after the initial height, with no extra sample.
/// A missing bracket returns Circumpolar (no event), not an acquisition error.
/// A bracket gets at most 60 bisections, a 0.05-second symmetric slope probe
/// and the 0.03-second grazing exit; the ordinary 2..100 arcsec/sec slope band
/// admits one in-bracket Newton correction. Returns a UT1 JD or no event;
/// callback/derived nonfinite errors propagate. Callbacks remain synchronous.
/// This bounded solve makes no certified physical-error or wider-body claim.
#[allow(clippy::too_many_arguments)]
#[cfg(test)]
fn native_rise_set_search(
    jd_start: f64,
    lat_deg: f64,
    rising: bool,
    event_height: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    get_ra_dec: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    horizon_alt_deg: f64,
    lift_max_deg: f64,
) -> Result<RiseSetOutcome, Error> {
    native_rise_set_search_with_policy(
        jd_start,
        lat_deg,
        rising,
        event_height,
        get_ra_dec,
        horizon_alt_deg,
        lift_max_deg,
        true,
    )
}

/// Same sample/bracket/refinement arithmetic, with checked admission retained
/// for planetary callers. The star policy lets actual declination/altitude
/// callbacks consume their values first and preserves IEEE branch semantics;
/// it does not turn a failed provider callback into a no-event result.
#[allow(clippy::too_many_arguments)]
pub(super) fn native_rise_set_search_with_policy(
    jd_start: f64,
    lat_deg: f64,
    rising: bool,
    event_height: &mut impl FnMut(f64) -> Result<f64, Error>,
    get_ra_dec: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    horizon_alt_deg: f64,
    lift_max_deg: f64,
    checked: bool,
) -> Result<RiseSetOutcome, Error> {
    if checked && !jd_start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: jd_start });
    }
    if checked && !lat_deg.is_finite() {
        return Err(Error::InvalidObserver {
            field: "lat_deg",
            value: lat_deg,
        });
    }
    if checked && (!horizon_alt_deg.is_finite() || !lift_max_deg.is_finite()) {
        return Err(Error::NonFiniteSearchValue { x: jd_start });
    }
    let height = &mut |jd: f64| -> Result<f64, Error> {
        let h = event_height(jd)?;
        if checked && !h.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok(h)
    };
    let declination = &mut |jd: f64| -> Result<f64, Error> {
        let (_, dec) = get_ra_dec(jd)?;
        if checked && !dec.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok(dec)
    };

    // Circumpolar pre-check from the declination extremes, widened by a
    // finite-difference declination-speed margin for fast movers (the Moon
    // moves ~13 deg/day, so its declination shifts inside the window).
    let dec = declination(jd_start)?;
    let dec2 = declination(jd_start + 1.0 / 24.0)?;
    let dec_speed = (dec2 - dec).abs() * 24.0;
    let margin = dec_speed * events::RISE_SET_MAX_SEARCH_DAYS;

    let max_alt = 90.0 - (lat_deg - dec).abs();
    let min_alt = if lat_deg >= 0.0 {
        dec - (90.0 - lat_deg)
    } else {
        -dec - (90.0 + lat_deg)
    };
    if min_alt - margin > horizon_alt_deg || max_alt + margin + lift_max_deg < horizon_alt_deg {
        return Ok(RiseSetOutcome::Circumpolar);
    }

    // Step size: near-grazing geometries can stay above/below the horizon
    // for less than an hour, so sample more densely there.
    let grazing_margin = (max_alt - horizon_alt_deg)
        .abs()
        .min((min_alt - horizon_alt_deg).abs());
    let search_step = if grazing_margin < 2.0 {
        10.0 / 1440.0
    } else if grazing_margin < 5.0 {
        20.0 / 1440.0
    } else {
        1.0 / 24.0
    };

    let mut h_prev = height(jd_start)?;
    let mut bracket: Option<(f64, f64)> = None;
    let n_steps = (events::RISE_SET_MAX_SEARCH_DAYS / search_step) as usize;
    for i in 1..=n_steps {
        let jd = jd_start + i as f64 * search_step;
        let h = height(jd)?;
        let crossed = if rising {
            h_prev < 0.0 && h >= 0.0
        } else {
            h_prev >= 0.0 && h < 0.0
        };
        if crossed {
            bracket = Some((jd - search_step, jd));
            break;
        }
        h_prev = h;
    }
    let (mut lo, mut hi) = match bracket {
        Some(b) => b,
        None => return Ok(RiseSetOutcome::Circumpolar),
    };

    for _ in 0..events::RISE_SET_BISECT_ITER {
        let jd_mid = 0.5 * (lo + hi);
        let h_mid = height(jd_mid)?;
        if h_mid.abs() < events::RISE_SET_ALT_TOLERANCE_DEG {
            let span = hi - lo;
            let mut dt = events::RISE_SET_SLOPE_PROBE_S / 86400.0;
            if dt > 0.25 * span {
                dt = 0.25 * span;
            }
            let mut slope = 0.0;
            if dt > 0.0 {
                let slope_signed =
                    (height(jd_mid + dt)? - height(jd_mid - dt)?) / (2.0 * dt * 86400.0);
                slope = slope_signed.abs();
            }
            if slope <= 0.0
                || slope >= events::RISE_SET_GRAZE_SLOPE_DEG_S
                || h_mid.abs() < slope * events::RISE_SET_TIME_TOL_S
            {
                let jd_ut = if slope > 0.0 {
                    events::refine_horizon_crossing(
                        height,
                        lo,
                        hi,
                        rising,
                        slope >= events::RISE_SET_DIP_SLOPE_DEG_S,
                    )?
                } else {
                    jd_mid
                };
                return Ok(RiseSetOutcome::Found { jd_ut });
            }
        }
        if rising {
            if h_mid < 0.0 {
                lo = jd_mid;
            } else {
                hi = jd_mid;
            }
        } else if h_mid >= 0.0 {
            lo = jd_mid;
        } else {
            hi = jd_mid;
        }
    }

    Ok(RiseSetOutcome::Found {
        jd_ut: 0.5 * (lo + hi),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn native_scan_does_not_acquire_forbidden_next_sample() {
        let calls = RefCell::new(Vec::new());
        let mut height = |jd: f64| {
            calls.borrow_mut().push(jd);
            assert!(jd <= 1002.0);
            Ok(-10.0)
        };
        let result = native_rise_set_search(
            1000.0,
            0.0,
            true,
            &mut height,
            &mut |_| Ok((0.0, 0.0)),
            0.0,
            0.0,
        )
        .unwrap();
        assert!(matches!(result, RiseSetOutcome::Circumpolar));
        assert_eq!(calls.borrow().len(), 49);
        assert_eq!(
            calls.borrow().last().unwrap().to_bits(),
            1002.0f64.to_bits()
        );
    }
    #[test]
    fn event_only_after_last_admitted_sample_is_not_an_event() {
        let result = native_rise_set_search(
            1000.0,
            0.0,
            true,
            &mut |jd| Ok(if jd > 1002.0 { 10.0 } else { -10.0 }),
            &mut |_| Ok((0.0, 0.0)),
            0.0,
            0.0,
        )
        .unwrap();
        assert!(matches!(result, RiseSetOutcome::Circumpolar));
    }
    #[test]
    fn final_admitted_scan_crossing_is_retained() {
        let result = native_rise_set_search(
            1000.0,
            0.0,
            true,
            &mut |jd| Ok(jd - 1002.0),
            &mut |_| Ok((0.0, 0.0)),
            0.0,
            0.0,
        )
        .unwrap();
        match result {
            RiseSetOutcome::Found { jd_ut } => assert!((jd_ut - 1002.0).abs() < 1e-6),
            _ => panic!("last admitted bracket lost"),
        }
    }
}
