//! Public scalar-position crossing orchestration for NativeJpl sessions.
//!
//! A callback receives the exact flags for each acquisition. The initial
//! geocentric scan uses the caller's flags; Newton, bracket and guard samples
//! request SPEED explicitly. Position errors are wrapped only at the entry and
//! Newton stages that own a public crossing error. Recovery preserves a fresh
//! range error and does not reinterpret unrelated input or I/O failures.
use super::{Error, ErrorCategory, SearchDirection};
use crate::constants::{FLG_HELCTR, FLG_SPEED, SUN};

/// Angular stopping threshold in degrees; one milliarcsecond.
const TOLERANCE: f64 = 0.001 / 3600.0;
/// Time-width guard in days for Brent's bisection fallback.
const BRENT_TIME_EPS: f64 = 1e-8;

/// Positive-divisor floating modulo, including positive zero at exact multiples.
/// Negative remainders receive one turn; NaN remains NaN.
fn modulo360(value: f64) -> f64 {
    let remainder = value % 360.0;
    if remainder < 0.0 {
        remainder + 360.0
    } else if remainder == 0.0 {
        0.0
    } else {
        remainder
    }
}

/// Normalize a degree target to one positive turn, retaining the modulo's IEEE
/// zero and NaN semantics. Radians conversion belongs to the public adapter.
pub(super) fn normalize_target(target: f64) -> f64 {
    modulo360(target)
}

/// Signed residual in (-180,180], retaining the positive antipodal endpoint.
/// Callers choose operand order: Newton uses target minus longitude, while
/// bracket refinement uses longitude minus target.
fn signed(value: f64) -> f64 {
    let diff = modulo360(value);
    if diff > 180.0 { diff - 360.0 } else { diff }
}

/// Scan residual in [-180,180); the chronological scan owns this spelling.
fn scan_delta(longitude: f64, target: f64) -> f64 {
    modulo360(longitude - target + 180.0) - 180.0
}

/// Public base crossing failure outside an active exception handler.
fn failure(message: impl Into<String>) -> Error {
    Error::PublicCrossing {
        message: message.into(),
        cause: None,
    }
}

/// Report failed crossing recovery with its native cause.
fn failure_with_cause(message: impl Into<String>, cause: Error) -> Error {
    Error::PublicCrossing {
        message: message.into(),
        cause: Some(Box::new(cause)),
    }
}

/// Standard Rust binary64 diagnostic formatting, independent of providers.
pub(crate) fn float_text(value: f64) -> String {
    value.to_string()
}

/// Native input/domain failures that terminate a position iteration.
pub(in crate::engine_api) fn input_failure(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::NonFiniteArgument { .. }
            | Error::ZeroVectorNormalize
            | Error::ZeroVectorToSpherical
            | Error::SingularMatrix
            | Error::NormalizedTimeOutOfRange { .. }
            | Error::InvalidSegment { .. }
            | Error::InvalidDerivativeOrder { .. }
            | Error::MathDomain { .. }
            | Error::EphemerisLoaderValue { .. }
            | Error::UnsupportedFlagCombination { .. }
    )
}

/// Position-stage catch for calculation/range failures and numerical values.
/// A base crossing error is not a CalculationError, despite its coarse category.
fn position_error(error: &Error) -> bool {
    input_failure(error)
        || error.category() == ErrorCategory::Calculation
            && !matches!(error.primary_error(), Error::PublicCrossing { .. })
}

/// Whether a bracket recovery catches this public Error/RuntimeError branch.
/// Built-in value, arithmetic and loader I/O failures do not enter that catch.
/// Typed body/file lookup and configuration failures remain public Error branches.
fn recoverable(error: &Error) -> bool {
    !input_failure(error)
        && !matches!(
            error.primary_error(),
            Error::EphemerisLoaderIo { .. }
                | Error::IntegerOverflow
                | Error::RationalDivisionByZero
                | Error::NumericalOverflow { .. }
        )
}

/// A genuine provider coverage failure eligible for range recovery.
fn is_range_error(error: &Error) -> bool {
    matches!(
        error.primary_error(),
        Error::EphemerisRange { .. }
            | Error::LebOutOfRange { .. }
            | Error::SpkOutOfRange { .. }
            | Error::EpochOutOfSegment { .. }
    )
}

/// Wrap a selected position failure with the operation's stage and owned cause.
fn stage(error: Error, name: &str) -> Error {
    if position_error(&error) {
        Error::PublicCrossing {
            message: format!("Failed to calculate {name}: {error}"),
            cause: Some(Box::new(error)),
        }
    } else {
        error
    }
}

/// Capture the Newton loop's cap from its initial daily speed.
/// Later rates may enlarge fallback caps but do not extend this loop.
fn adaptive_cap(speed: f64) -> usize {
    let speed = speed.abs();
    if speed >= 0.1 {
        60
    } else if speed >= 0.001 {
        80
    } else {
        100
    }
}

/// Integer sample count, refusing non-finite or unaddressable values before cast.
/// All ordinary scan windows fit usize; this guard does not cap their size.
fn sample_count(value: f64) -> Result<usize, Error> {
    if value.is_nan() {
        Err(Error::MathDomain { value })
    } else if !value.is_finite() || value < 0.0 || value >= usize::MAX as f64 {
        Err(Error::IntegerOverflow)
    } else {
        Ok(value as usize)
    }
}

/// Sample component zero and its daily component-three derivative with SPEED.
fn rate(
    sample: &mut impl FnMut(f64, i32) -> Result<[f64; 6], Error>,
    jd: f64,
    flags: i32,
) -> Result<(f64, f64), Error> {
    let value = sample(jd, flags | FLG_SPEED)?;
    Ok((value[0], value[3]))
}

/// Uniform chronological bracket scan. Both endpoints are acquired through
/// the callback; no endpoint shortcut, finite-value filter or epoch substitution
/// is applied. Sign changes jumping 180 degrees or more are antipodal wraps.
fn bracket(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    target: f64,
    start: f64,
    end: f64,
    samples: usize,
    from_end: bool,
) -> Result<(f64, f64), Error> {
    let step = (end - start) / samples as f64;
    let mut delta = |jd| get(jd).map(|p| signed(p.0 - target));
    if from_end {
        let mut next = end;
        let mut next_delta = delta(end)?;
        for i in (0..samples).rev() {
            let current = start + i as f64 * step;
            let current_delta = delta(current)?;
            if current_delta * next_delta <= 0.0 && (current_delta - next_delta).abs() < 180.0 {
                return Ok((current, next));
            }
            next = current;
            next_delta = current_delta;
        }
    } else {
        let mut previous = start;
        let mut previous_delta = delta(start)?;
        for i in 1..=samples {
            let current = start + i as f64 * step;
            let current_delta = delta(current)?;
            if previous_delta * current_delta <= 0.0
                && (current_delta - previous_delta).abs() < 180.0
            {
                return Ok((previous, current));
            }
            previous = current;
            previous_delta = current_delta;
        }
    }
    Err(failure(format!(
        "No crossing found in interval [{}, {}] for target {}°",
        float_text(start),
        float_text(end),
        float_text(target)
    )))
}

/// Brent refinement with inverse quadratic, secant and bisection steps.
/// The best endpoint is returned after the literal cap; exhaustion performs
/// no additional acquisition. The residual has angular units and the width
/// guard has day units. Callback failures propagate to the owning outer catch.
fn brent(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    target: f64,
    mut a: f64,
    mut b: f64,
    tolerance: f64,
    max_iter: usize,
) -> Result<f64, Error> {
    let mut f = |jd| get(jd).map(|p| signed(p.0 - target));
    let mut fa = f(a)?;
    let mut fb = f(b)?;
    if fa * fb > 0.0 {
        return Err(failure(format!(
            "Brent's method: root not bracketed. f(a)={fa:.6}, f(b)={fb:.6}"
        )));
    }
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    let mut c = a;
    let mut fc = fa;
    let mut mflag = true;
    let mut d = 0.0;
    for _ in 0..max_iter {
        if fb.abs() < tolerance {
            return Ok(b);
        }
        let mut s = if fa != fc && fb != fc {
            a * fb * fc / ((fa - fb) * (fa - fc))
                + b * fa * fc / ((fb - fa) * (fb - fc))
                + c * fa * fb / ((fc - fa) * (fc - fb))
        } else {
            b - fb * (b - a) / (fb - fa)
        };
        let cond1 = !((3.0 * a + b) / 4.0 < s && s < b || b < s && s < (3.0 * a + b) / 4.0);
        let cond2 = mflag && (s - b).abs() >= (b - c).abs() / 2.0;
        let cond3 = !mflag && (s - b).abs() >= (c - d).abs() / 2.0;
        let cond4 = mflag && (b - c).abs() < BRENT_TIME_EPS;
        let cond5 = !mflag && (c - d).abs() < BRENT_TIME_EPS;
        if cond1 || cond2 || cond3 || cond4 || cond5 {
            s = (a + b) / 2.0;
            mflag = true;
        } else {
            mflag = false;
        }
        let fs = f(s)?;
        d = c;
        c = b;
        fc = fb;
        if fa * fs < 0.0 {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }
        if fa.abs() < fb.abs() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }
    }
    Ok(b)
}

/// Sixty midpoint acquisitions refine a chronological sign-change bracket.
fn bisect(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    target: f64,
    mut lo: f64,
    mut hi: f64,
    mut f_lo: f64,
) -> Result<f64, Error> {
    for _ in 0..60 {
        let mid = 0.5 * (lo + hi);
        let f_mid = scan_delta(get(mid)?.0, target);
        if f_lo * f_mid <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    Ok(0.5 * (lo + hi))
}

/// Scan for an earlier crossing, shrinking steps in the target's retrograde
/// band. A speed sign change triggers a fifty-step station bisection to detect
/// a grazing pair that two longitude samples would otherwise miss.
fn adaptive_first(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    target: f64,
    start: f64,
    end: f64,
    max_rate: f64,
    band: f64,
) -> Result<Option<f64>, Error> {
    let (lon, mut v_prev) = get(start)?;
    let mut d_prev = scan_delta(lon, target);
    let mut t_prev = start;
    let mut t = start;
    while t < end {
        let gap = d_prev.abs();
        let in_band = gap <= band;
        let step = if in_band {
            2.0
        } else {
            300.0_f64.min(2.0_f64.max((gap - band) / max_rate))
        };
        t = (t + step).min(end);
        let (lon, v_cur) = get(t)?;
        let d_cur = scan_delta(lon, target);
        if d_prev * d_cur < 0.0 && (d_cur - d_prev).abs() < 180.0 {
            return bisect(get, target, t_prev, t, d_prev).map(Some);
        }
        if in_band && v_prev * v_cur < 0.0 {
            let (mut s_lo, mut s_hi, mut sv_lo) = (t_prev, t, v_prev);
            for _ in 0..50 {
                let s_mid = 0.5 * (s_lo + s_hi);
                let v_mid = get(s_mid)?.1;
                if sv_lo * v_mid <= 0.0 {
                    s_hi = s_mid;
                } else {
                    s_lo = s_mid;
                    sv_lo = v_mid;
                }
            }
            let station = 0.5 * (s_lo + s_hi);
            let d_station = scan_delta(get(station)?.0, target);
            if d_station == 0.0 {
                return Ok(Some(station));
            }
            if d_prev * d_station < 0.0 {
                return bisect(get, target, t_prev, station, d_prev).map(Some);
            }
        }
        t_prev = t;
        d_prev = d_cur;
        v_prev = v_cur;
    }
    Ok(None)
}

/// Revalidate a far on-table result against the chronological first crossing.
/// Near results and off-table bodies do not acquire this extra scan.
fn first_guard(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    body: i32,
    target: f64,
    start: f64,
    candidate: f64,
) -> Result<f64, Error> {
    if !(2..=9).contains(&body) || candidate - start <= 5.0 {
        return Ok(candidate);
    }
    let earlier = adaptive_first(
        get,
        target,
        start,
        candidate - 0.5,
        crate::events::slow_max_rate_deg_day(body),
        crate::events::retro_arc_deg(body),
    )?;
    Ok(match earlier {
        Some(t) if start - 1e-6 <= t && t < candidate - 1.0 => t,
        _ => candidate,
    })
}

/// Forward range recovery, using at least 200 and at most 4000 intervals,
/// followed by a hundred-step Brent refinement and the optional first guard.
fn forward(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    body: i32,
    target: f64,
    start: f64,
    horizon: f64,
) -> Result<f64, Error> {
    let samples = sample_count(horizon / 0.5)?.clamp(200, 4000);
    let (a, b) = bracket(get, target, start, start + horizon, samples, false)?;
    let t = brent(get, target, a, b, TOLERANCE, 100)?;
    first_guard(get, body, target, start, t)
}

/// Sun or Moon longitude crossing in the callback's time scale.
/// `body` is SUN or MOON; the private caller fixes that choice.
/// Fifty Newton acquisitions follow one entry sample. The seed uses a mean
/// rate only at zero speed; the iteration floor is 0.01 deg/day for the Sun,
/// 0.1 deg/day for the Moon. Backward seeds exclude an at-start crossing.
pub(super) fn luminary(
    sample: &mut impl FnMut(f64, i32) -> Result<[f64; 6], Error>,
    body: i32,
    target: f64,
    start: f64,
    flags: i32,
    direction: SearchDirection,
) -> Result<f64, Error> {
    let (name, mean_speed, speed_floor, span, divergence, exhaustion) = if body == SUN {
        (
            "Sun",
            0.9856,
            0.01,
            366.0,
            "Solar crossing search diverged",
            "Maximum iterations reached in solar crossing calculation",
        )
    } else {
        (
            "Moon",
            13.176,
            0.1,
            31.0,
            "Moon crossing search diverged",
            "Maximum iterations reached in moon crossing calculation",
        )
    };
    let (lon, mut speed) =
        rate(sample, start, flags).map_err(|e| stage(e, &format!("{name} position")))?;
    if speed == 0.0 {
        speed = mean_speed;
    }
    let mut diff = modulo360(target - lon);
    if direction == SearchDirection::Backward {
        let signed = if diff > 180.0 { diff - 360.0 } else { diff };
        if signed.abs() < 1e-6 {
            diff = -360.0;
        } else {
            diff -= 360.0;
        }
    } else if speed < 0.0 && diff > 0.0 {
        diff -= 360.0;
    }
    let guess = start + diff / speed;
    let mut jd = guess;
    for _ in 0..50 {
        let (lon, mut speed) = rate(sample, jd, flags)
            .map_err(|e| stage(e, &format!("{name} position during iteration")))?;
        let diff = signed(target - lon);
        if diff.abs() < TOLERANCE {
            return Ok(jd);
        }
        if speed.abs() < speed_floor {
            speed = mean_speed;
        }
        jd += diff / speed;
        if (jd - guess).abs() > span {
            return Err(failure(divergence));
        }
    }
    Err(failure(exhaustion))
}

/// First forward geocentric longitude crossing. Raw scan acquisitions retain
/// the caller's flags and reacquire the start before sampling start+0.5 days.
/// A scan return is final; only later Newton/Brent results use the rate guard.
/// Newton range recovery keeps a fresh failure's already acquired channels;
/// an unchained fresh failure inherits the active range as context. Re-raising
/// the original range after a caught lookup instead replaces implicit context,
/// retaining any explicit cause and suppression. The callback owns acquisition
/// channels; this orchestration does not reconstruct a missing provider chain.
pub(super) fn planet(
    sample: &mut impl FnMut(f64, i32) -> Result<[f64; 6], Error>,
    body: i32,
    target: f64,
    start: f64,
    flags: i32,
) -> Result<f64, Error> {
    let (lon_start, initial_speed) =
        rate(sample, start, flags).map_err(|e| stage(e, "planet position"))?;
    let mut speed_default = crate::events::typical_speed_deg_day(body);
    let mut off_table_horizon = None;
    if !(2..=9).contains(&body) && 1e-6 < initial_speed.abs() && initial_speed.abs() < 0.45 {
        speed_default = initial_speed.abs();
        off_table_horizon = Some((1.2 * 360.0 / speed_default).min(40000.0));
    }
    let mut diff = modulo360(target - lon_start);
    let effective_speed = if initial_speed < 0.0 {
        speed_default
    } else if initial_speed.abs() > 0.001 {
        initial_speed
    } else {
        speed_default
    };
    let diff_back = modulo360(lon_start - target);
    if initial_speed < 0.0
        || initial_speed.abs() < 0.01
        || diff_back < crate::events::retro_arc_deg(body)
    {
        let mut window = crate::events::synodic_days(body) * 1.2;
        if let Some(horizon) = off_table_horizon {
            window = window.max(horizon);
        }
        let mut delta = |jd| sample(jd, flags).map(|v| scan_delta(v[0], target));
        let mut d_prev = delta(start)?;
        let mut t_prev = start;
        let mut t_scan = start + 0.5;
        while t_scan <= start + window {
            let d_cur = delta(t_scan)?;
            if d_prev == 0.0 {
                return Ok(t_prev);
            }
            if d_prev * d_cur < 0.0 && (d_cur - d_prev).abs() < 180.0 {
                let (mut lo, mut hi) = (t_prev, t_scan);
                let mut f_lo = delta(lo)?;
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    let f_mid = delta(mid)?;
                    if f_lo * f_mid <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        f_lo = f_mid;
                    }
                }
                return Ok(0.5 * (lo + hi));
            }
            t_prev = t_scan;
            d_prev = d_cur;
            t_scan += 0.5;
        }
    }
    let dt_guess = diff / effective_speed;
    let mut jd = start + dt_guess;
    let loop_cap = adaptive_cap(initial_speed);
    let mut fallback_cap = loop_cap;
    let mut get = |jd| rate(sample, jd, flags);
    if initial_speed.abs() < 0.01 {
        let attempt = (|| {
            let window = if speed_default > 0.0 {
                let minimum: f64 = if speed_default < 0.1 { 500.0 } else { 60.0 };
                minimum
                    .max(dt_guess.abs() * 2.0)
                    .max((diff / speed_default).abs() * 3.0)
            } else {
                500.0_f64.max(dt_guess.abs() * 2.0)
            };
            let samples = sample_count(window / 30.0)?.max(40);
            let (a, b) = bracket(&mut get, target, start, start + window, samples, false)?;
            let t = brent(&mut get, target, a, b, TOLERANCE, fallback_cap)?;
            first_guard(&mut get, body, target, start, t)
        })();
        match attempt {
            Ok(t) => return Ok(t),
            Err(error) if recoverable(&error) => {}
            Err(error) => return Err(error),
        }
    }
    let mut station_fallback_triggered = false;
    let mut horizon = crate::events::cross_max_range_days(body, speed_default);
    if let Some(off_table) = off_table_horizon {
        horizon = horizon.max(off_table);
    }
    for iteration in 0..loop_cap {
        let (lon, mut speed) = match get(jd) {
            Ok(value) => value,
            Err(range_error) if is_range_error(&range_error) => {
                return match forward(&mut get, body, target, start, horizon) {
                    Ok(t) => Ok(t),
                    Err(new_range) if is_range_error(&new_range) => Err(new_range),
                    Err(error) if recoverable(&error) => Err(range_error),
                    Err(error) => Err(error),
                };
            }
            Err(error) => return Err(stage(error, "planet position during iteration")),
        };
        diff = signed(target - lon);
        if diff.abs() < TOLERANCE {
            if jd >= start - 1e-6 {
                return first_guard(&mut get, body, target, start, jd);
            }
            return match forward(&mut get, body, target, start, horizon) {
                Ok(t) => Ok(t),
                Err(range_error) if is_range_error(&range_error) => Err(range_error),
                Err(error) if recoverable(&error) => {
                    Err(failure_with_cause("Planet crossing search diverged", error))
                }
                Err(error) => Err(error),
            };
        }
        if (speed.abs() < 0.01 || speed < 0.0) && !station_fallback_triggered {
            station_fallback_triggered = true;
            let attempt = (|| {
                let window = if speed_default > 0.0 {
                    let minimum: f64 = if speed_default < 0.1 { 500.0 } else { 60.0 };
                    minimum
                        .max(dt_guess.abs() * 2.0)
                        .max((diff / speed_default).abs() * 3.0)
                } else {
                    500.0_f64.max(dt_guess.abs() * 2.0)
                };
                let samples = sample_count(window / 30.0)?.max(40);
                let (a, b) = bracket(&mut get, target, start, start + window, samples, false)?;
                let t = brent(&mut get, target, a, b, TOLERANCE, fallback_cap - iteration)?;
                first_guard(&mut get, body, target, start, t)
            })();
            match attempt {
                Ok(t) => return Ok(t),
                Err(error) if recoverable(&error) => {}
                Err(error) => return Err(error),
            }
        }
        fallback_cap = fallback_cap.max(adaptive_cap(speed));
        if speed.abs() < 0.001 {
            speed = speed_default;
        }
        jd += diff / speed;
        if (jd - start).abs() > horizon {
            return match forward(&mut get, body, target, start, horizon) {
                Ok(t) => Ok(t),
                Err(range_error) if is_range_error(&range_error) => Err(range_error),
                Err(error) if recoverable(&error) => {
                    Err(failure_with_cause("Planet crossing search diverged", error))
                }
                Err(error) => Err(error),
            };
        }
    }
    match forward(&mut get, body, target, start, horizon) {
        Ok(t) => Ok(t),
        Err(range_error) if is_range_error(&range_error) => Err(range_error),
        Err(error) if recoverable(&error) => Err(failure_with_cause(
            "Maximum iterations reached in planet crossing calculation",
            error,
        )),
        Err(error) => Err(error),
    }
}

/// Latitude-zero result in TT, retaining actual components zero and one.
pub(super) struct Node {
    /// Scalar TT Julian Day at the accepted crossing.
    pub(super) time: f64,
    /// Component zero in the caller's selected frame and units.
    pub(super) longitude: f64,
    /// Component one in the caller's selected frame and units.
    pub(super) latitude: f64,
}

/// First direction-aware lunar component-one zero, over scalar TT samples.
/// Component four supplies Newton's derivative. A two-day chronological scan
/// and sixty midpoint samples select the first root; fifty Newton samples
/// refine it. The at-start epsilon is 1e-6 days, reduced to 2.5e-7 days when
/// the entry sample estimates an imminent distinct root. No stored stencil or
/// configured Delta-T replaces these position components.
pub(super) fn node(
    sample: &mut impl FnMut(f64, i32) -> Result<[f64; 6], Error>,
    start: f64,
    flags: i32,
    direction: SearchDirection,
) -> Result<Node, Error> {
    let dir = if direction == SearchDirection::Backward {
        -1.0
    } else {
        1.0
    };
    let initial = sample(start, flags | FLG_SPEED).map_err(|e| stage(e, "Moon position"))?;
    let lat = initial[1];
    let mut accept_eps = 1e-6;
    let mut guess = start + 13.6 * dir;
    if lat.abs() < 10.0 * TOLERANCE {
        let t_est = if initial[4].abs() > 0.1 {
            -lat / initial[4]
        } else {
            0.0
        };
        if t_est * dir > 2.5e-7 {
            guess = start + t_est;
            accept_eps = 2.5e-7;
        }
    } else {
        let step = 2.0 * dir;
        let mut lat_prev = lat;
        let mut t_prev = start;
        let mut moon_lat = |jd| {
            sample(jd, flags | FLG_SPEED)
                .map(|v| v[1])
                .map_err(|e| stage(e, "Moon position"))
        };
        for i in 1..=9 {
            let current = start + i as f64 * step;
            let lat_current = moon_lat(current)?;
            if lat_prev * lat_current < 0.0 {
                let (mut lo, mut hi, mut f_lo) = (t_prev, current, lat_prev);
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    let f_mid = moon_lat(mid)?;
                    if f_lo * f_mid <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        f_lo = f_mid;
                    }
                }
                guess = 0.5 * (lo + hi);
                break;
            }
            t_prev = current;
            lat_prev = lat_current;
        }
    }
    let mut jd = guess;
    for _ in 0..50 {
        let value = sample(jd, flags | FLG_SPEED)
            .map_err(|e| stage(e, "Moon position during iteration"))?;
        let lat = value[1];
        let mut speed = value[4];
        if lat.abs() < TOLERANCE {
            if direction == SearchDirection::Forward && jd > start + accept_eps
                || direction == SearchDirection::Backward && jd < start - accept_eps
            {
                return Ok(Node {
                    time: jd,
                    longitude: value[0],
                    latitude: lat,
                });
            }
            jd += (13.6 / 2.0) * dir;
            continue;
        }
        if speed.abs() < 0.1 {
            speed = if lat >= 0.0 { 1.0 } else { -1.0 };
        }
        jd -= lat / speed;
        if direction == SearchDirection::Forward && jd < start {
            jd = start + 13.6 / 2.0;
        } else if direction == SearchDirection::Backward && jd > start {
            jd = start - 13.6 / 2.0;
        } else if (jd - start).abs() > 30.0 {
            return Err(failure("Moon node crossing search diverged"));
        }
    }
    Err(failure(
        "Maximum iterations reached in moon node crossing calculation",
    ))
}

/// Days excluded before an at-start backward root, scaled by mean daily rate.
fn backward_exclusion(rate: f64) -> f64 {
    let rate = rate.abs();
    if rate <= 1e-9 {
        1.0
    } else {
        1e-6_f64.max(10.0 * TOLERANCE / rate)
    }
}

/// Nearest strictly-past repair. Error/RuntimeError may keep Newton's result;
/// arithmetic, value and filesystem failures remain visible to the caller.
fn nearest_past(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    target: f64,
    start: f64,
    newton: f64,
    period: f64,
) -> Result<f64, Error> {
    let margin = 1.0_f64.max(0.01 * period);
    let lo = newton.min(start) + margin;
    let hi = start - margin;
    if hi - lo >= margin {
        let attempt = (|| {
            let samples = sample_count((hi - lo) / 15.0)?.max(40);
            let (a, b) = bracket(get, target, lo, hi, samples, true)?;
            brent(get, target, a, b, TOLERANCE, 100)
        })();
        match attempt {
            Ok(t) => return Ok(t),
            Err(error) if recoverable(&error) => {}
            Err(error) => return Err(error),
        }
    }
    Ok(newton)
}

/// Repair a past Newton root with a forward scan beginning start+1e-6 days.
fn helio_forward(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    target: f64,
    start: f64,
    period: f64,
) -> Result<f64, Error> {
    let window = 1.25 * period;
    let samples = sample_count(window / 30.0)?.max(60);
    let (a, b) = bracket(get, target, start + 1e-6, start + window, samples, false)?;
    brent(get, target, a, b, TOLERANCE, 100)
}

/// Out-of-table heliocentric search, expanding only the far end of a bracket
/// by powers of two through six expansions. A forward at-start root is included;
/// a backward at-start root is excluded using the measured mean-rate seed.
/// Only public Error failures expand the window; a refinement error propagates.
fn helio_bracketed(
    get: &mut impl FnMut(f64) -> Result<(f64, f64), Error>,
    target: f64,
    start: f64,
    mean_rate: f64,
    direction: SearchDirection,
) -> Result<f64, Error> {
    let period = if mean_rate > 1e-9 {
        360.0 / mean_rate
    } else {
        720.0
    };
    let base_window = 3.0 * period;
    let lon_start = get(start)?.0;
    let start_is_root = scan_delta(lon_start, target).abs() <= TOLERANCE;
    if direction == SearchDirection::Forward && start_is_root {
        return Ok(start);
    }
    let mut near = start;
    if direction == SearchDirection::Backward && start_is_root {
        near -= backward_exclusion(mean_rate);
    } else if direction == SearchDirection::Forward {
        near += 1e-6;
    }
    let mut last_error = None;
    for expansion in 0..=6 {
        let window = base_window * 2.0_f64.powi(expansion);
        let samples = sample_count(window / 15.0)?.clamp(120, 4000);
        let bracket_result = if direction == SearchDirection::Backward {
            bracket(get, target, start - window, near, samples, true)
        } else {
            bracket(get, target, near, start + window, samples, false)
        };
        let (a, b) = match bracket_result {
            Ok(pair) => pair,
            Err(error) if recoverable(&error) => {
                last_error = Some(error);
                continue;
            }
            Err(error) => return Err(error),
        };
        return brent(get, target, a, b, TOLERANCE, 200);
    }
    Err(last_error.unwrap_or_else(|| failure("No crossing found")))
}

/// Heliocentric longitude crossing in the callback's input time scale.
/// Every sample adds HELCTR and SPEED. Mean-rate backward seeds, nearest-past
/// validation and out-of-table bracket expansion preserve the chosen side of
/// the start epoch. The initial adaptive cap is frozen for Newton; later rates
/// update only the fallback-cap variable. Exhaustion reads no extra position.
pub(super) fn helio(
    sample: &mut impl FnMut(f64, i32) -> Result<[f64; 6], Error>,
    body: i32,
    target: f64,
    start: f64,
    flags: i32,
    direction: SearchDirection,
) -> Result<f64, Error> {
    let flags = flags | FLG_HELCTR | FLG_SPEED;
    let (lon_start, mut speed) =
        rate(sample, start, flags).map_err(|e| stage(e, "heliocentric planet position"))?;
    let mut speed_default = crate::events::helio_typical_speed_deg_day(body);
    if speed_default <= 0.0 {
        speed_default = if speed.abs() > 1e-6 { speed.abs() } else { 0.5 };
    }
    let mut diff = modulo360(target - lon_start);
    // Keep both ordered backward-exclusion comparisons explicit.
    #[allow(clippy::if_same_then_else)]
    let dt_guess = if direction == SearchDirection::Backward {
        if diff > 1e-5 {
            diff -= 360.0;
        } else if diff > -1e-5 {
            diff -= 360.0;
        }
        let guess_speed = if speed_default > 1e-9 {
            speed_default
        } else {
            1.0
        };
        diff / guess_speed
    } else {
        if speed.abs() < 0.0001 {
            speed = speed_default;
        }
        diff / speed
    };
    let mut jd = start + dt_guess;
    let loop_cap = adaptive_cap(speed);
    let mut get = |jd| rate(sample, jd, flags);
    if !crate::events::helio_speed_known(body) {
        return helio_bracketed(&mut get, target, start, speed_default, direction);
    }
    if speed.abs() < 0.01 {
        let attempt = (|| {
            let window = if speed_default > 0.0 {
                60.0_f64
                    .max(dt_guess.abs() * 2.0)
                    .max((diff / speed_default).abs() * 1.5)
            } else {
                120.0_f64.max(dt_guess.abs() * 2.0)
            };
            let (lo, mut hi) = if direction == SearchDirection::Backward {
                (start - window, start)
            } else {
                (start, start + window)
            };
            if direction == SearchDirection::Backward
                && scan_delta(lon_start, target).abs() <= TOLERANCE
            {
                hi -= backward_exclusion(speed_default);
            }
            let samples = sample_count(window / 30.0)?.max(40);
            let (a, b) = bracket(
                &mut get,
                target,
                lo,
                hi,
                samples,
                direction == SearchDirection::Backward,
            )?;
            brent(&mut get, target, a, b, TOLERANCE, loop_cap)
        })();
        match attempt {
            Ok(t) => return Ok(t),
            Err(error) if recoverable(&error) => {}
            Err(error) => return Err(error),
        }
    }
    // The mutable fallback cap has no later consumer on this branch, but its
    // update does not change the already captured Newton iteration count.
    let mut _fallback_cap = loop_cap;
    for _ in 0..loop_cap {
        let (lon, mut speed) =
            get(jd).map_err(|e| stage(e, "heliocentric position during iteration"))?;
        diff = signed(target - lon);
        if diff.abs() < TOLERANCE {
            let period = if speed_default > 1e-9 {
                360.0 / speed_default
            } else {
                720.0
            };
            if direction == SearchDirection::Backward {
                return nearest_past(&mut get, target, start, jd, period);
            }
            if jd < start - 1e-6 {
                return helio_forward(&mut get, target, start, period);
            }
            return Ok(jd);
        }
        _fallback_cap = _fallback_cap.max(adaptive_cap(speed));
        if speed.abs() < 0.0001 {
            speed = speed_default;
        }
        jd += diff / speed;
        let base_range: f64 = if speed_default.abs() < 0.05 {
            500.0
        } else {
            400.0
        };
        let period = if speed_default > 1e-9 {
            360.0 / speed_default
        } else {
            720.0
        };
        let horizon = base_range.max(dt_guess.abs() * 2.0).max(3.0 * period);
        if (jd - start).abs() > horizon {
            return Err(failure("Heliocentric crossing search diverged"));
        }
    }
    Err(failure(
        "Maximum iterations reached in heliocentric crossing calculation",
    ))
}

#[cfg(test)]
#[path = "native_crossing_tests.rs"]
mod tests;
