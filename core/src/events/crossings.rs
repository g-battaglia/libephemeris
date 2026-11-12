//! Crossings implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Signed wrapped residual `(lon - target)` in `(-180, 180]` degrees.
/// A genuine crossing passes smoothly through 0; the antipodal wrap
/// (`target + 180`) jumps by ~360 and is filtered by the bracket scans.
pub fn wrapped_residual_deg(lon_deg: f64, target_deg: f64) -> f64 {
    let mut diff = (lon_deg - target_deg) % 360.0;
    if diff > 180.0 {
        diff -= 360.0;
    } else if diff <= -180.0 {
        diff += 360.0;
    }
    diff
}

/// Signed Newton residual `(target - lon)` in `(-180, 180]` degrees.
///
/// This is the release Newton-`diff` spelling (`(x2cross - lon) % 360`
/// folded past 180 into the negative half, as in `solcross_ut`,
/// `mooncross_ut` and `cross_ut`): the correction step is
/// `jd += diff / speed`. It is the negation of [`wrapped_residual_deg`],
/// which keeps the `(lon - target)` convention for the sign-agnostic
/// bracket/Brent machinery. The Newton steps must not reuse the bracket
/// residual: an earlier revision stepped away from the root whenever
/// the linear seed was not already exact (masked by exact-seed linear
/// tests); the sealed-serving sweep exposed the divergence.
pub fn newton_residual_deg(lon_deg: f64, target_deg: f64) -> f64 {
    -wrapped_residual_deg(lon_deg, target_deg)
}

/// Positive-divisor modulo into `[0, 360)`: Rust's `%` keeps the sign of the
/// dividend while the release `(x2cross - lon) % 360.0` always lands in
/// `[0, 360)`. Every forward/backward seed distance must use this form.
pub(super) fn mod360_deg(x: f64) -> f64 {
    let r = x % 360.0;
    if r < 0.0 { r + 360.0 } else { r }
}

/// Normalize a crossing target into `[0, 360)` degrees.
pub fn normalize_target_deg(target_deg: f64) -> f64 {
    mod360_deg(target_deg)
}

/// Adaptive Newton cap from the instantaneous speed, degrees/day.
pub fn adaptive_max_iterations(speed_deg_day: f64) -> usize {
    let s = speed_deg_day.abs();
    if s >= 0.1 {
        NR_MAX_ITER_PLANET
    } else if s >= 0.001 {
        // Release spells the 0.01 and 0.001 bands as two caps
        // (`NR_MAX_ITER_HELIO` / `NR_MAX_ITER_VERY_SLOW`), both 80.
        NR_MAX_ITER_SLOW
    } else {
        NR_MAX_ITER_STATION
    }
}

/// Station predicate: Newton divides by the
/// speed, so below this floor the search must use Brent's method.
pub fn is_near_station(speed_deg_day: f64) -> bool {
    speed_deg_day.abs() < STATION_SPEED_THRESHOLD_DEG_DAY
}

/// Geocentric mean direct-motion speed seed per planet id, degrees/day
/// (Meeus ch. 31 order-of-magnitude
/// seeds, non-final). Off-table bodies fall back to `0.5`.
pub fn typical_speed_deg_day(planet_id: i32) -> f64 {
    match planet_id {
        2 => 1.4,
        3 => 1.2,
        4 => 0.5,
        5 => 0.08,
        6 => 0.03,
        7 => 0.01,
        8 => 0.006,
        9 => 0.004,
        _ => 0.5,
    }
}

/// Forward search horizon in days, sized to cover at least one orbital
/// period.
pub fn cross_max_range_days(planet_id: i32, speed_default_deg_day: f64) -> f64 {
    let s = speed_default_deg_day.abs();
    if s < 0.01 {
        100000.0
    } else if s < 0.05 {
        40000.0
    } else if s < 0.1 {
        5000.0
    } else if planet_id == 2 || planet_id == 3 {
        500.0
    } else {
        800.0
    }
}

/// Retrograde-loop arc seed per planet, degrees.
pub fn retro_arc_deg(planet_id: i32) -> f64 {
    match planet_id {
        2 => 16.0,
        3 => 18.0,
        4 => 22.0,
        5 => 12.0,
        6 => 9.0,
        7 => 6.0,
        8 => 5.0,
        9 => 4.0,
        _ => 22.0,
    }
}

/// Synodic period seed per planet, days.
pub fn synodic_days(planet_id: i32) -> f64 {
    match planet_id {
        2 => 116.0,
        3 => 584.0,
        4 => 780.0,
        5 => 399.0,
        6 => 378.0,
        7 => 370.0,
        8 => 367.0,
        9 => 367.0,
        _ => 800.0,
    }
}

/// Conservative upper bound on geocentric daily motion per planet,
/// degrees/day.
pub fn slow_max_rate_deg_day(planet_id: i32) -> f64 {
    match planet_id {
        2 => 2.5,
        3 => 1.35,
        4 => 0.85,
        5 => 0.30,
        6 => 0.16,
        7 => 0.08,
        8 => 0.05,
        9 => 0.05,
        _ => 0.16,
    }
}

/// Whether the guarded first-crossing re-validation applies.
pub fn needs_first_crossing_guard(planet_id: i32) -> bool {
    (2..=9).contains(&planet_id)
}

/// Days to exclude before the start epoch in a backward search so an
/// endpoint root is stepped over without hiding a genuine distinct
/// crossing. Kept for the future
/// heliocentric path; the geocentric backward seed below handles the
/// endpoint root through the signed `diff` mapping instead.
pub fn backward_start_exclusion_days(rate_deg_day: f64, tolerance_deg: f64) -> f64 {
    let rate = rate_deg_day.abs();
    if rate <= 1e-9 {
        1.0
    } else {
        (10.0 * tolerance_deg / rate).max(1e-6)
    }
}

/// Direction-aware angular distance from `lon_start` to `target`, in
/// degrees: forward searches return a value in `[0, 360)` (a just-passed
/// target maps near 360, i.e. next cycle; a body exactly on target yields
/// 0, i.e. now), backward searches return a value in `(-360, 0]` (an
/// at-target start maps to `-360`, i.e. a full cycle back, per the
/// strictly-past contract). Mirrors the seed-`diff` blocks of
/// `solcross_ut`/`mooncross_ut`, including the 1e-6 deg at-crossing guard
/// and the forward retrograde guard.
pub fn seed_diff_deg(
    lon_start_deg: f64,
    target_deg: f64,
    direction: SearchDirection,
    speed_deg_day: f64,
) -> f64 {
    let diff = mod360_deg(target_deg - lon_start_deg);
    match direction {
        SearchDirection::Backward => {
            let signed = if diff > 180.0 { diff - 360.0 } else { diff };
            if signed.abs() < 1e-6 {
                -360.0
            } else {
                signed - 360.0 * (signed > 0.0) as u8 as f64
            }
        }
        SearchDirection::Forward => {
            if speed_deg_day < 0.0 && diff > 0.0 {
                diff - 360.0
            } else {
                diff
            }
        }
    }
}

/// Bracket `[a, b]` (`a < b`) holding a longitude crossing, found by a
/// uniform scan with the antipodal-wrap filter: a sign change whose jump
/// reaches 180 deg is the `target + 180` wrap, never a crossing. `from_end` scans backward from `hi` so
/// the first hit is the crossing nearest `hi`.
pub fn find_crossing_bracket(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    lo: f64,
    hi: f64,
    samples: usize,
    from_end: bool,
) -> Result<(f64, f64)> {
    if !lo.is_finite() || !hi.is_finite() || lo >= hi {
        return Err(Error::InvalidSearchInterval { a: lo, b: hi });
    }
    if samples == 0 {
        return Err(Error::InvalidScanSamples { count: 0 });
    }
    let eval = &mut |jd: f64| -> Result<f64> {
        let (lon, _speed) = get_position(jd)?;
        if !lon.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok(wrapped_residual_deg(lon, target_deg))
    };
    let step = (hi - lo) / samples as f64;
    if from_end {
        let mut next_jd = hi;
        let mut next_diff = eval(hi)?;
        for i in (0..samples).rev() {
            let curr_jd = lo + i as f64 * step;
            let curr_diff = eval(curr_jd)?;
            if curr_diff * next_diff <= 0.0 && (curr_diff - next_diff).abs() < 180.0 {
                return Ok((curr_jd, next_jd));
            }
            next_jd = curr_jd;
            next_diff = curr_diff;
        }
    } else {
        let mut prev_jd = lo;
        let mut prev_diff = eval(lo)?;
        for i in 1..=samples {
            let curr_jd = if i == samples {
                hi
            } else {
                lo + i as f64 * step
            };
            let curr_diff = eval(curr_jd)?;
            if prev_diff * curr_diff <= 0.0 && (prev_diff - curr_diff).abs() < 180.0 {
                return Ok((prev_jd, curr_jd));
            }
            prev_jd = curr_jd;
            prev_diff = curr_diff;
        }
    }
    Err(Error::NoCrossingInInterval {
        lo,
        hi,
        target: target_deg,
    })
}

/// Refine a bracketed crossing with Brent's method over the wrapped
/// residual. Returns the refined crossing
/// time; an unbracketed interval is `NoCrossingInInterval`, a non-finite
/// interpolant yields the best estimate as `CrossingNotConverged`.
pub fn brent_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    a: f64,
    b: f64,
    tolerance_deg: f64,
    max_iter: usize,
) -> Result<f64> {
    if !tolerance_deg.is_finite() || tolerance_deg <= 0.0 {
        return Err(Error::InvalidSearchTolerance {
            value: tolerance_deg,
        });
    }
    if max_iter == 0 {
        return Err(Error::InvalidSearchIterations { max_iter });
    }
    let mut f = |jd: f64| -> Result<f64> {
        let (lon, _) = get_position(jd)?;
        if !lon.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok(wrapped_residual_deg(lon, target_deg))
    };
    let mut xa = a;
    let mut xb = b;
    let mut fa = f(xa)?;
    let mut fb = f(xb)?;
    if fa.abs() <= tolerance_deg {
        return Ok(xa);
    }
    if fb.abs() <= tolerance_deg {
        return Ok(xb);
    }
    if fa * fb > 0.0 {
        return Err(Error::NoCrossingInInterval {
            lo: a,
            hi: b,
            target: target_deg,
        });
    }
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut xa, &mut xb);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = xa;
    let mut fc = fa;
    let mut mflag = true;
    let mut d = 0.0;
    for _ in 0..max_iter {
        if fb.abs() < tolerance_deg {
            return Ok(xb);
        }
        let s = if fa != fc && fb != fc {
            xa * fb * fc / ((fa - fb) * (fa - fc))
                + xb * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            xb - fb * (xb - xa) / (fb - fa)
        };
        let in_cone = (3.0 * xa + xb) / 4.0 < s && s < xb || xb < s && s < (3.0 * xa + xb) / 4.0;
        let cond1 = !in_cone;
        let cond2 = mflag && (s - xb).abs() >= (xb - c).abs() / 2.0;
        let cond3 = !mflag && (s - xb).abs() >= (c - d).abs() / 2.0;
        let cond4 = mflag && (xb - c).abs() < BRENT_TIME_EPS_DAYS;
        let cond5 = !mflag && (c - d).abs() < BRENT_TIME_EPS_DAYS;
        let s = if cond1 || cond2 || cond3 || cond4 || cond5 {
            mflag = true;
            0.5 * (xa + xb)
        } else {
            mflag = false;
            s
        };
        if !s.is_finite() {
            return Err(Error::CrossingNotConverged {
                best: xb,
                residual: fb,
            });
        }
        let fs = f(s)?;
        d = c;
        c = xb;
        fc = fb;
        if fa * fs < 0.0 {
            xb = s;
            fb = fs;
        } else {
            xa = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            std::mem::swap(&mut xa, &mut xb);
            std::mem::swap(&mut fa, &mut fb);
        }
    }
    Ok(xb)
}

/// First crossing at or after `start` via a forward bracket scan plus
/// Brent refinement. Forward-only by
/// construction; raises `NoCrossingInInterval` when none lies within
/// `max_range_days`.
pub fn forward_first_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    start: f64,
    max_range_days: f64,
    tolerance_deg: f64,
) -> Result<f64> {
    if !max_range_days.is_finite() || max_range_days <= 0.0 {
        return Err(Error::InvalidSearchStep {
            value: max_range_days,
        });
    }
    let samples = (200).max((max_range_days / 0.5) as usize).min(4000);
    let (a, b) = find_crossing_bracket(
        get_position,
        target_deg,
        start,
        start + max_range_days,
        samples,
        false,
    )?;
    brent_crossing(get_position, target_deg, a, b, tolerance_deg, 100)
}

/// Bisect the wrapped residual on `[lo, hi]` (60 iterations), the shared
/// refinement tail of the adaptive forward scan.
pub(super) fn bisect_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    lo: f64,
    hi: f64,
    f_lo: f64,
) -> Result<f64> {
    let mut lo = lo;
    let mut hi = hi;
    let mut f_lo = f_lo;
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let (lon, _) = get_position(mid)?;
        if !lon.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: mid });
        }
        let f_mid = wrapped_residual_deg(lon, target_deg);
        if f_lo * f_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// Earliest crossing in `[lo, hi]`, or `Ok(None)` when there is none
/// a skip-free adaptive step
/// — at most `(gap - band) / max_rate` far from the target, `fine_step`
/// inside the band — plus exact station-graze repair (a speed sign flip
/// inside a fine interval refines the station by speed bisection and the
/// extremum's sign decides the hidden double crossing).
#[allow(clippy::too_many_arguments)]
pub fn first_forward_crossing_adaptive(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    lo: f64,
    hi: f64,
    max_rate_deg_day: f64,
    band_deg: f64,
    fine_step_days: f64,
    coarse_cap_days: f64,
) -> Result<Option<f64>> {
    if !lo.is_finite() || !hi.is_finite() || lo >= hi {
        return Err(Error::InvalidSearchInterval { a: lo, b: hi });
    }
    if !max_rate_deg_day.is_finite()
        || max_rate_deg_day <= 0.0
        || !band_deg.is_finite()
        || band_deg <= 0.0
        || !fine_step_days.is_finite()
        || fine_step_days <= 0.0
        || !coarse_cap_days.is_finite()
        || coarse_cap_days <= 0.0
    {
        return Err(Error::InvalidSearchStep {
            value: max_rate_deg_day,
        });
    }
    let sample = &mut |jd: f64| -> Result<(f64, f64)> {
        let (lon, speed) = get_position(jd)?;
        if !lon.is_finite() || !speed.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok((wrapped_residual_deg(lon, target_deg), speed))
    };
    let (mut d_prev, mut v_prev) = sample(lo)?;
    let mut t_prev = lo;
    let mut t = lo;
    while t < hi {
        let gap = d_prev.abs();
        let in_band = gap <= band_deg;
        let step = if in_band {
            fine_step_days
        } else {
            coarse_cap_days.min(fine_step_days.max((gap - band_deg) / max_rate_deg_day))
        };
        t = (t + step).min(hi);
        let (d_cur, v_cur) = sample(t)?;
        if d_prev * d_cur < 0.0 && (d_cur - d_prev).abs() < 180.0 {
            return Ok(Some(bisect_crossing(
                get_position,
                target_deg,
                t_prev,
                t,
                d_prev,
            )?));
        }
        if in_band && v_prev * v_cur < 0.0 {
            let (mut s_lo, mut s_hi) = (t_prev, t);
            let mut sv_lo = v_prev;
            for _ in 0..50 {
                let s_mid = 0.5 * (s_lo + s_hi);
                let v_mid = sample(s_mid)?.1;
                if sv_lo * v_mid <= 0.0 {
                    s_hi = s_mid;
                } else {
                    s_lo = s_mid;
                    sv_lo = v_mid;
                }
            }
            let t_station = 0.5 * (s_lo + s_hi);
            let d_station = sample(t_station)?.0;
            if d_station == 0.0 {
                return Ok(Some(t_station));
            }
            if d_prev * d_station < 0.0 {
                return Ok(Some(bisect_crossing(
                    get_position,
                    target_deg,
                    t_prev,
                    t_station,
                    d_prev,
                )?));
            }
        }
        t_prev = t;
        d_prev = d_cur;
        v_prev = v_cur;
    }
    Ok(None)
}

/// Earliest crossing at or after `start`, re-validating a far Newton
/// result against an earlier crossing it may have skipped. Engages only for the slow on-table planets
/// and candidates more than 5 days ahead; a more precise Newton value is
/// kept when no earlier crossing exists.
pub fn guarded_first_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    planet_id: i32,
    target_deg: f64,
    start: f64,
    candidate: f64,
) -> Result<f64> {
    if !needs_first_crossing_guard(planet_id) || candidate - start <= 5.0 {
        return Ok(candidate);
    }
    let earlier = first_forward_crossing_adaptive(
        get_position,
        target_deg,
        start,
        candidate - 0.5,
        slow_max_rate_deg_day(planet_id),
        retro_arc_deg(planet_id),
        2.0,
        300.0,
    )?;
    if let Some(t) = earlier {
        if start - 1e-6 <= t && t < candidate - 1.0 {
            return Ok(t);
        }
    }
    Ok(candidate)
}

/// Sun/Moon-style Newton crossing search over an abstract provider
/// direction-aware
/// linear seed, Newton refinement with the body speed floor, and the
/// divergence window. `seed_speed` replaces a zero starting speed,
/// `floor_speed` replaces near-zero speeds mid-iteration.
#[allow(clippy::too_many_arguments)]
pub fn newton_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    start: f64,
    direction: SearchDirection,
    tolerance_deg: f64,
    max_iter: usize,
    seed_speed_deg_day: f64,
    floor_speed_deg_day: f64,
    divergence_days: f64,
) -> Result<f64> {
    if !start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start });
    }
    if !tolerance_deg.is_finite() || tolerance_deg <= 0.0 {
        return Err(Error::InvalidSearchTolerance {
            value: tolerance_deg,
        });
    }
    if max_iter == 0 {
        return Err(Error::InvalidSearchIterations { max_iter });
    }
    let target = normalize_target_deg(target_deg);
    let (lon_start, mut speed) = checked_position(get_position, start)?;
    if speed == 0.0 {
        speed = seed_speed_deg_day;
    }
    let diff = seed_diff_deg(lon_start, target, direction, speed);
    let jd_guess = start + diff / speed;
    let mut jd = jd_guess;
    for _ in 0..max_iter {
        let (lon, spd) = checked_position(get_position, jd)?;
        let mut step_speed = spd;
        let diff = newton_residual_deg(lon, target);
        if diff.abs() < tolerance_deg {
            return Ok(jd);
        }
        if step_speed.abs() < floor_speed_deg_day {
            step_speed = seed_speed_deg_day;
        }
        jd += diff / step_speed;
        if (jd - jd_guess).abs() > divergence_days {
            return Err(Error::CrossingDiverged {
                start,
                estimate: jd,
            });
        }
    }
    Err(Error::CrossingNotConverged {
        best: jd,
        residual: newton_residual_deg(checked_position(get_position, jd)?.0, target),
    })
}

/// Evaluate the provider once, rejecting non-finite longitudes/speeds.
pub(super) fn checked_position(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    jd: f64,
) -> Result<(f64, f64)> {
    let (lon, speed) = get_position(jd)?;
    if !lon.is_finite() || !speed.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: jd });
    }
    Ok((lon, speed))
}

/// Sun longitude crossing over an abstract longitude/rate provider.
///
/// `get_position(jd)` returns longitude degrees and degrees/day in the target's
/// chosen frame. `start` and the result are Julian days in that provider's
/// common clock; this reducer performs no TT/UT1 conversion. `target_deg` is
/// degrees and `direction` selects the search direction. Numerical stopping,
/// endpoint and error behavior follow [`newton_crossing`].
pub fn sun_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    start: f64,
    direction: SearchDirection,
) -> Result<f64> {
    newton_crossing(
        get_position,
        target_deg,
        start,
        direction,
        CROSSING_TOLERANCE_DEG,
        NR_MAX_ITER_SUN,
        SUN_MEAN_RATE_DEG_DAY,
        STATION_SPEED_THRESHOLD_DEG_DAY,
        SUN_DIVERGENCE_DAYS,
    )
}

/// Moon longitude crossing over an abstract longitude/rate provider.
///
/// Provider units and clock ownership follow [`sun_crossing`]; Moon-specific
/// iteration/rate bounds are used. The result is a Julian day in the provider's
/// clock, not implicitly UTC or TT. Consumed failures follow [`newton_crossing`].
pub fn moon_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    start: f64,
    direction: SearchDirection,
) -> Result<f64> {
    newton_crossing(
        get_position,
        target_deg,
        start,
        direction,
        CROSSING_TOLERANCE_DEG,
        NR_MAX_ITER_MOON,
        MOON_MEAN_RATE_DEG_DAY,
        MOONCROSS_SPEED_FLOOR_DEG_DAY,
        MOON_DIVERGENCE_DAYS,
    )
}

/// First forward planet longitude crossing over an abstract provider.
///
/// The callback maps a Julian day to `(longitude_deg, speed_deg_per_day)` in
/// one frame. `start` and the returned date share that provider's clock;
/// `target_deg` is degrees. `planet_id` selects search heuristics, while the
/// callback owns the actual target, so the caller must keep them aligned.
/// No backward search or clock conversion is performed here. Nonfinite sampled
/// values and unresolved searches return typed errors; provider acquisition
/// and numerical refinement remain part of this bounded search.
pub fn planet_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    planet_id: i32,
    target_deg: f64,
    start: f64,
) -> Result<f64> {
    if !start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start });
    }
    let target = normalize_target_deg(target_deg);
    let (lon_start, speed_start) = checked_position(get_position, start)?;
    let mut speed_default = typical_speed_deg_day(planet_id);
    let mut max_range = cross_max_range_days(planet_id, speed_default);
    if !typical_speed_known(planet_id) && speed_start.abs() > 1e-6 && speed_start.abs() < 0.45 {
        // Slow off-table body (e.g. Mean Node): derive the horizon from
        // the actual instantaneous rate instead of the generic fallback.
        speed_default = speed_start.abs();
        let horizon = (1.2 * 360.0 / speed_default).min(40000.0);
        max_range = max_range.max(horizon);
    }

    let diff = mod360_deg(target - lon_start);
    let effective_speed = if speed_start < 0.0 {
        speed_default
    } else if speed_start.abs() > 0.001 {
        speed_start
    } else {
        speed_default
    };

    // Retro-dip pre-scan: currently retrograde, near a station, or the
    // target lies behind within one retrograde-loop arc — Newton would
    // project a full orbit ahead and return a crossing months too late.
    let diff_back = mod360_deg(lon_start - target);
    if speed_start < 0.0 || is_near_station(speed_start) || diff_back < retro_arc_deg(planet_id) {
        let mut scan_window = synodic_days(planet_id) * 1.2;
        if !typical_speed_known(planet_id) && speed_start.abs() > 1e-6 && speed_start.abs() < 0.45 {
            scan_window = scan_window.max((1.2 * 360.0 / speed_default).min(40000.0));
        }
        let d_prev = wrapped_residual_deg(lon_start, target);
        let mut t_prev = start;
        let mut prev = d_prev;
        let mut t_scan = start + 0.5;
        while t_scan <= start + scan_window {
            let (lon, _) = checked_position(get_position, t_scan)?;
            let cur = wrapped_residual_deg(lon, target);
            if prev == 0.0 {
                return guarded_first_crossing(get_position, planet_id, target, start, t_prev);
            }
            if prev * cur < 0.0 && (cur - prev).abs() < 180.0 {
                let refined = bisect_crossing(get_position, target, t_prev, t_scan, prev)?;
                return guarded_first_crossing(get_position, planet_id, target, start, refined);
            }
            t_prev = t_scan;
            prev = cur;
            t_scan += 0.5;
        }
    }

    let dt_guess = diff / effective_speed;
    let jd_guess = start + dt_guess;
    let mut max_iter = adaptive_max_iterations(speed_start);

    // Near-station start: Brent over a bracketed window, Newton only as
    // a last resort.
    if is_near_station(speed_start) {
        let window = station_search_window(speed_default, dt_guess, diff);
        let end = start + window;
        let samples = (40).max((window / 30.0) as usize);
        match find_crossing_bracket(get_position, target, start, end, samples, false) {
            Ok((a, b)) => {
                let hit =
                    brent_crossing(get_position, target, a, b, CROSSING_TOLERANCE_DEG, max_iter)?;
                return guarded_first_crossing(get_position, planet_id, target, start, hit);
            }
            Err(_) => {
                // Fall through to Newton, as the release does.
            }
        }
    }

    let mut jd = jd_guess;
    let mut station_fallback_done = false;
    for iteration in 0.. {
        if iteration >= max_iter {
            break;
        }
        let step = match checked_position(get_position, jd) {
            Ok(ok) => ok,
            Err(first) => {
                // Newton wandered out of provider coverage: the in-range
                // forward scan may still hold a genuine crossing; otherwise surface the
                // original provider error.
                return match forward_first_crossing(
                    get_position,
                    target,
                    start,
                    max_range,
                    CROSSING_TOLERANCE_DEG,
                ) {
                    Ok(hit) => guarded_first_crossing(get_position, planet_id, target, start, hit),
                    Err(_) => Err(first),
                };
            }
        };
        let (lon, speed) = step;
        let diff = newton_residual_deg(lon, target);
        if diff.abs() < CROSSING_TOLERANCE_DEG {
            if jd >= start - 1e-6 {
                return guarded_first_crossing(get_position, planet_id, target, start, jd);
            }
            // Converged to a PAST crossing: the forward-only contract
            // requires the first crossing at or after `start`.
            return match forward_first_crossing(
                get_position,
                target,
                start,
                max_range,
                CROSSING_TOLERANCE_DEG,
            ) {
                Ok(hit) => guarded_first_crossing(get_position, planet_id, target, start, hit),
                Err(_) => Err(Error::CrossingDiverged {
                    start,
                    estimate: jd,
                }),
            };
        }
        if (is_near_station(speed) || speed < 0.0) && !station_fallback_done {
            station_fallback_done = true;
            let window = station_search_window(speed_default, dt_guess, diff);
            let end = start + window;
            let samples = (40).max((window / 30.0) as usize);
            if let Ok((a, b)) =
                find_crossing_bracket(get_position, target, start, end, samples, false)
            {
                if let Ok(hit) = brent_crossing(
                    get_position,
                    target,
                    a,
                    b,
                    CROSSING_TOLERANCE_DEG,
                    max_iter.saturating_sub(iteration),
                ) {
                    return guarded_first_crossing(get_position, planet_id, target, start, hit);
                }
            }
        }
        max_iter = max_iter.max(adaptive_max_iterations(speed));
        let step_speed = if speed.abs() < 0.001 {
            speed_default
        } else {
            speed
        };
        jd += diff / step_speed;
        if (jd - start).abs() > max_range {
            return match forward_first_crossing(
                get_position,
                target,
                start,
                max_range,
                CROSSING_TOLERANCE_DEG,
            ) {
                Ok(hit) => guarded_first_crossing(get_position, planet_id, target, start, hit),
                Err(_) => Err(Error::CrossingDiverged {
                    start,
                    estimate: jd,
                }),
            };
        }
    }

    // Newton exhausted its cap inside the range without converging: the
    // forward bracket scan is the robust catch-all.
    match forward_first_crossing(
        get_position,
        target,
        start,
        max_range,
        CROSSING_TOLERANCE_DEG,
    ) {
        Ok(hit) => guarded_first_crossing(get_position, planet_id, target, start, hit),
        Err(_) => Err(Error::CrossingNotConverged {
            best: jd,
            residual: newton_residual_deg(checked_position(get_position, jd)?.0, target),
        }),
    }
}

/// Whether `planet_id` has a calibrated typical-speed entry.
pub(super) fn typical_speed_known(planet_id: i32) -> bool {
    (2..=9).contains(&planet_id)
}

/// Brent search window for station/retrograde starts.
pub(super) fn station_search_window(
    speed_default_deg_day: f64,
    dt_guess_days: f64,
    diff_deg: f64,
) -> f64 {
    if speed_default_deg_day > 0.0 {
        let min_window: f64 = if speed_default_deg_day < 0.1 {
            500.0
        } else {
            60.0
        };
        min_window
            .max(dt_guess_days.abs() * 2.0)
            .max((diff_deg / speed_default_deg_day).abs() * 3.0)
    } else {
        500.0_f64.max(dt_guess_days.abs() * 2.0)
    }
}

// H: heliocentric longitude crossing search layer.

/// Heliocentric mean direct-motion speed seed per planet id,
/// degrees/day. Off-table bodies fall back to the measured rate (see
/// [`helio_crossing`]).
pub fn helio_typical_speed_deg_day(planet_id: i32) -> f64 {
    match planet_id {
        2 => 4.09,
        3 => 1.60,
        14 => 0.986,
        4 => 0.524,
        5 => 0.083,
        6 => 0.034,
        7 => 0.012,
        8 => 0.006,
        9 => 0.004,
        _ => 0.0,
    }
}

/// Whether `planet_id` has a calibrated heliocentric speed entry
pub fn helio_speed_known(planet_id: i32) -> bool {
    matches!(planet_id, 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 14)
}

/// Nearest strictly-past heliocentric crossing, validating a Newton
/// result: scans the interior
/// of `(jd_newton, start)` — both endpoints guarded by `margin` so a
/// body sitting at the target at `start` (which must resolve a full
/// cycle back) and a just-passed crossing are both preserved — for a
/// genuine crossing nearer to `start`. Returns it when one exists,
/// otherwise keeps `jd_newton`. A failed scan keeps `jd_newton`, like
/// the release `except` tail.
pub fn nearest_past_helio_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    start: f64,
    jd_newton: f64,
    period_days: f64,
    tolerance_deg: f64,
) -> Result<f64> {
    let margin = 1.0_f64.max(0.01 * period_days);
    let lo = jd_newton.min(start) + margin;
    let hi = start - margin;
    if hi - lo >= margin {
        let samples = (40).max(((hi - lo) / 15.0) as usize);
        if let Ok((a, b)) = find_crossing_bracket(get_position, target_deg, lo, hi, samples, true) {
            if let Ok(hit) = helio_brent(get_position, target_deg, a, b, tolerance_deg, 100) {
                return Ok(hit);
            }
        }
    }
    Ok(jd_newton)
}

/// First heliocentric crossing strictly after `start`: forward repair used when Newton
/// converges into the basin of a past revolution. Scans one
/// revolution-plus ahead of `start` and refines the first bracket.
pub fn first_forward_helio_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    start: f64,
    period_days: f64,
    tolerance_deg: f64,
) -> Result<f64> {
    let window = 1.25 * period_days;
    let samples = (60).max((window / 30.0) as usize);
    let (a, b) = find_crossing_bracket(
        get_position,
        target_deg,
        start + 1e-6,
        start + window,
        samples,
        false,
    )?;
    helio_brent(get_position, target_deg, a, b, tolerance_deg, 100)
}

/// Bracket-and-refine heliocentric crossing for out-of-table slow
/// bodies: the heliocentric
/// longitude is strictly monotonic, so the first bracketed crossing
/// adjacent to `start` on the requested side is the correct one. The
/// far end of the window expands geometrically (up to 64x) because an
/// eccentric orbit sampled near perihelion overestimates the mean rate;
/// the near end stays anchored at `start` so an expanded scan still
/// returns the adjacent crossing. A root exactly at `start` resolves to
/// `start` itself going forward (the no-dead-band contract) and is
/// stepped over going backward (strictly-past semantics) — but only
/// when one is actually there.
pub fn helio_cross_bracketed(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    start: f64,
    mean_rate_deg_day: f64,
    direction: SearchDirection,
    tolerance_deg: f64,
) -> Result<f64> {
    let period_days = if mean_rate_deg_day > 1e-9 {
        360.0 / mean_rate_deg_day
    } else {
        720.0
    };
    let base_window = 3.0 * period_days;
    let (lon_start, _) = checked_position(get_position, start)?;
    // Release residual spelling (no antipodal fold): `% 360` into
    // `[0, 360)` minus 180.
    let start_residual = mod360_deg(lon_start - target_deg + 180.0) - 180.0;
    let start_is_root = start_residual.abs() <= tolerance_deg;
    if direction == SearchDirection::Forward && start_is_root {
        return Ok(start);
    }
    let mut jd_near = start;
    if direction == SearchDirection::Backward && start_is_root {
        jd_near -= backward_start_exclusion_days(mean_rate_deg_day, tolerance_deg);
    } else if direction == SearchDirection::Forward {
        jd_near += 1e-6;
    }
    let mut last_error: Option<Error> = None;
    for expansion in 0..=6 {
        let window = base_window * 2.0_f64.powi(expansion);
        let samples = (120).max((window / 15.0) as usize).min(4000);
        let bracket = if direction == SearchDirection::Backward {
            find_crossing_bracket(
                get_position,
                target_deg,
                start - window,
                jd_near,
                samples,
                true,
            )
        } else {
            find_crossing_bracket(
                get_position,
                target_deg,
                jd_near,
                start + window,
                samples,
                false,
            )
        };
        match bracket {
            Ok((a, b)) => {
                return helio_brent(get_position, target_deg, a, b, tolerance_deg, 200);
            }
            Err(exc) => {
                last_error = Some(exc);
            }
        }
    }
    Err(last_error.unwrap_or(Error::NoCrossingInInterval {
        lo: start,
        hi: start,
        target: target_deg,
    }))
}

/// Refine a bracketed heliocentric crossing with Brent's method
/// An unbracketed interval is `NoCrossingInInterval`, a
/// non-finite interpolant yields the best estimate as
/// `CrossingNotConverged`.
pub fn helio_brent(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    target_deg: f64,
    a: f64,
    b: f64,
    tolerance_deg: f64,
    max_iter: usize,
) -> Result<f64> {
    if !tolerance_deg.is_finite() || tolerance_deg <= 0.0 {
        return Err(Error::InvalidSearchTolerance {
            value: tolerance_deg,
        });
    }
    if max_iter == 0 {
        return Err(Error::InvalidSearchIterations { max_iter });
    }
    let mut f = |jd: f64| -> Result<f64> {
        let (lon, _) = get_position(jd)?;
        if !lon.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok(wrapped_residual_deg(lon, target_deg))
    };
    let mut xa = a;
    let mut xb = b;
    let mut fa = f(xa)?;
    let mut fb = f(xb)?;
    if fa * fb > 0.0 {
        return Err(Error::NoCrossingInInterval {
            lo: a,
            hi: b,
            target: target_deg,
        });
    }
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut xa, &mut xb);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = xa;
    let mut fc = fa;
    let mut mflag = true;
    let mut d = 0.0;
    for _ in 0..max_iter {
        if fb.abs() < tolerance_deg {
            return Ok(xb);
        }
        let s = if fa != fc && fb != fc {
            xa * fb * fc / ((fa - fb) * (fa - fc))
                + xb * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            xb - fb * (xb - xa) / (fb - fa)
        };
        let in_cone = (3.0 * xa + xb) / 4.0 < s && s < xb || xb < s && s < (3.0 * xa + xb) / 4.0;
        let cond1 = !in_cone;
        let cond2 = mflag && (s - xb).abs() >= (xb - c).abs() / 2.0;
        let cond3 = !mflag && (s - xb).abs() >= (c - d).abs() / 2.0;
        let cond4 = mflag && (xb - c).abs() < BRENT_TIME_EPS_DAYS;
        let cond5 = !mflag && (c - d).abs() < BRENT_TIME_EPS_DAYS;
        let s = if cond1 || cond2 || cond3 || cond4 || cond5 {
            mflag = true;
            0.5 * (xa + xb)
        } else {
            mflag = false;
            s
        };
        if !s.is_finite() {
            return Err(Error::CrossingNotConverged {
                best: xb,
                residual: fb,
            });
        }
        let fs = f(s)?;
        d = c;
        c = xb;
        fc = fb;
        if fa * fs < 0.0 {
            xb = s;
            fb = fs;
        } else {
            xa = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            std::mem::swap(&mut xa, &mut xb);
            std::mem::swap(&mut fa, &mut fb);
        }
    }
    Ok(xb)
}

/// Heliocentric longitude crossing over an abstract provider. The ephemeris binding
/// (UT through Delta-T vs TT) arrives with the caller's
/// longitude source; the `FLG_RADIANS` target normalization belongs to
/// the binding layer.
pub fn helio_crossing(
    get_position: &mut impl FnMut(f64) -> std::result::Result<(f64, f64), Error>,
    planet_id: i32,
    target_deg: f64,
    start: f64,
    direction: SearchDirection,
) -> Result<f64> {
    if !start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start });
    }
    let target = normalize_target_deg(target_deg);
    let (lon_start, mut speed) = checked_position(get_position, start)?;
    let mut speed_default = helio_typical_speed_deg_day(planet_id);
    if speed_default <= 0.0 {
        // Off-table body (asteroids, centaurs, TNOs, fictitious ids):
        // estimate the mean rate from the measured heliocentric state —
        // on a bound orbit the instantaneous rate stays within a bounded
        // factor of the mean rate, which the window scaling absorbs.
        speed_default = if speed.abs() > 1e-6 { speed.abs() } else { 0.5 };
    }
    let backwards = direction == SearchDirection::Backward;
    let mut diff = mod360_deg(target - lon_start);
    let dt_guess = if backwards {
        // The two release branches (`diff > 1e-5`, `diff > -1e-5`)
        // collapse to one: any diff above `-1e-5` maps into the past so
        // a target essentially at the current longitude resolves a full
        // cycle back (nearest strictly-past crossing).
        if diff > -1e-5 {
            diff -= 360.0;
        }
        // Size the backward guess with the MEAN rate: near aphelion the
        // instantaneous speed runs well below the mean and would place
        // the guess a half-revolution too far back.
        let guess_speed = if speed_default > 1e-9 {
            speed_default
        } else {
            1.0
        };
        diff / guess_speed
    } else {
        // No forward dead-band: a body at (or just before) the target
        // has its first crossing now, not one revolution later.
        if speed.abs() < 0.0001 {
            speed = speed_default;
        }
        diff / speed
    };
    let jd_guess = start + dt_guess;
    let max_iter = adaptive_max_iterations(speed);

    if !helio_speed_known(planet_id) {
        // Out-of-table slow bodies: deterministic bracket-plus-Brent on
        // the strictly monotonic heliocentric longitude.
        return helio_cross_bracketed(
            get_position,
            target,
            start,
            speed_default,
            direction,
            CROSSING_TOLERANCE_DEG,
        );
    }

    // Near-station start: Brent over a bracketed window, Newton only as
    // a last resort.
    if is_near_station(speed) {
        let search_window = if speed_default > 0.0 {
            60.0_f64
                .max(dt_guess.abs() * 2.0)
                .max((diff / speed_default).abs() * 1.5)
        } else {
            120.0_f64.max(dt_guess.abs() * 2.0)
        };
        let (bracket_start, bracket_end) = if backwards {
            let mut end = start;
            let start_residual = mod360_deg(lon_start - target + 180.0) - 180.0;
            if start_residual.abs() <= CROSSING_TOLERANCE_DEG {
                end -= backward_start_exclusion_days(speed_default, CROSSING_TOLERANCE_DEG);
            }
            (start - search_window, end)
        } else {
            (start, start + search_window)
        };
        let samples = (40).max((search_window / 30.0) as usize);
        if let Ok((a, b)) = find_crossing_bracket(
            get_position,
            target,
            bracket_start,
            bracket_end,
            samples,
            backwards,
        ) {
            if let Ok(hit) =
                helio_brent(get_position, target, a, b, CROSSING_TOLERANCE_DEG, max_iter)
            {
                return Ok(hit);
            }
        }
        // Fall through to Newton-Raphson, as the release does.
    }

    let period_days = if speed_default > 1e-9 {
        360.0 / speed_default
    } else {
        720.0
    };
    let mut jd = jd_guess;
    // Fixed cap: the release `for... in range(max_iter)` binds the
    // range once, so the in-loop `max(max_iter,...)` refresh is dead.
    for _ in 0..max_iter {
        let (lon, step_speed) = checked_position(get_position, jd)?;
        let diff = newton_residual_deg(lon, target);
        if diff.abs() < CROSSING_TOLERANCE_DEG {
            if backwards {
                return nearest_past_helio_crossing(
                    get_position,
                    target,
                    start,
                    jd,
                    period_days,
                    CROSSING_TOLERANCE_DEG,
                );
            }
            if jd < start - 1e-6 {
                // Newton slipped into a past revolution; a forward
                // search must return a crossing after the start epoch.
                return first_forward_helio_crossing(
                    get_position,
                    target,
                    start,
                    period_days,
                    CROSSING_TOLERANCE_DEG,
                );
            }
            return Ok(jd);
        }
        let step_speed = if step_speed.abs() < 0.0001 {
            speed_default
        } else {
            step_speed
        };
        jd += diff / step_speed;
        // Revolution-scaled divergence guard: a single crossing lies at
        // most one revolution away, and Newton transients on an
        // eccentric slow orbit legitimately overshoot past one full
        // period before settling — the failsafe allows three.
        let base_range: f64 = if speed_default.abs() < 0.05 {
            500.0
        } else {
            400.0
        };
        let max_range = base_range.max(dt_guess.abs() * 2.0).max(3.0 * period_days);
        if (jd - start).abs() > max_range {
            return Err(Error::CrossingDiverged {
                start,
                estimate: jd,
            });
        }
    }
    Err(Error::CrossingNotConverged {
        best: jd,
        residual: newton_residual_deg(checked_position(get_position, jd)?.0, target),
    })
}

/// Next Moon longitude crossing over the sealed serving.
pub fn mooncross_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    target_deg: f64,
    start_ut: f64,
    direction: SearchDirection,
) -> Result<f64> {
    check_angular_serving_flags(serving.calc_flags(), "mooncross_ut")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_position = |jd: f64| serving.lon_speed_at_ut(jd);
    moon_crossing(&mut get_position, target_deg, start_ut, direction)
}

/// Next heliocentric longitude crossing over the sealed serving
///
/// The serving must carry `FLG_HELCTR` (the kernel forces it, like the
/// release `helio_flag = flags | FLG_HELCTR | FLG_SPEED` composition —
/// `FLG_SPEED` is added by the longitude feed itself); a serving bound
/// without it is refused fail-closed rather than searched
/// geocentrically. The epoch is UT through the serving Delta-T leg for
/// the `_ut` kernel and TT for the TT kernel (the B2
/// identity-timescale shape, reusing this same arm).
pub fn helio_crossing_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    target_deg: f64,
    start: f64,
    direction: SearchDirection,
) -> Result<f64> {
    check_angular_serving_flags(serving.calc_flags(), "helio_cross")?;
    if (serving.calc_flags() & crate::constants::FLG_HELCTR) == 0 {
        return Err(Error::UnsupportedFlagCombination {
            entry: "helio_cross",
            flags: serving.calc_flags(),
        });
    }
    if !start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start });
    }
    let mut get_position = |jd: f64| serving.lon_speed_at_ut(jd);
    helio_crossing(&mut get_position, planet_id, target_deg, start, direction)
}

/// Next lunar-node crossing over the sealed serving.
pub fn mooncross_node_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    start_ut: f64,
    direction: SearchDirection,
) -> Result<MoonNodeCrossing> {
    check_angular_serving_flags(serving.calc_flags(), "mooncross_node_ut")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_moon = |jd: f64| -> Result<(f64, f64, f64)> {
        // Latitude speed by central difference over the sealed triple
        // (same stencil association as the velocity path: the
        // release node solver differences latitudes, not longitudes).
        const H: f64 = 0.5 / 86400.0 * 10.0;
        let (lon, lat, _) = serving.lonlatdist_at_ut(jd)?;
        let (_, lat_m, _) = serving.lonlatdist_at_ut(jd - H)?;
        let (_, lat_p, _) = serving.lonlatdist_at_ut(jd + H)?;
        if !lon.is_finite() || !lat.is_finite() || !lat_m.is_finite() || !lat_p.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok((lon, lat, (lat_p - lat_m) / (2.0 * H)))
    };
    moon_node_crossing(&mut get_moon, start_ut, direction)
}
