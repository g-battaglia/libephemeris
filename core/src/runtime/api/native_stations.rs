//! Scalar-UT station orchestration. Every callback is the public position leg;
//! component three and the caller's coordinate flags retain their actual units.
//! No vector, model epoch, finite preflight or provider substitute is inserted.
use super::{Error, StationInformation, StationKind, StationMotionKind, StationRequest};
use crate::constants::{FLG_SPEED, MOON, SUN};

/// One original scalar UT1 acquisition, with the exact per-stage flags.
pub(super) type Sample<'a> = dyn FnMut(f64, i32) -> Result<[f64; 6], Error> + 'a;

/// Sun/Moon have no stations; the summary owns a distinct message. Unused time
/// and flags are not consumed by this identity check.
pub(super) fn check_body(body: i32, _summary: bool) -> Result<(), Error> {
    if matches!(body, SUN | MOON) {
        return Err(Error::BodyHasNoStations { body });
    }
    Ok(())
}

/// Station periods in days; other body IDs use 400 days. These size scans and
/// are not orbital solutions, domain checks or numerical fitting parameters.
fn period(body: i32) -> f64 {
    match body {
        2 => 116.0,
        3 => 584.0,
        4 => 780.0,
        5 => 399.0,
        6 => 378.0,
        7 => 370.0,
        8 | 9 => 367.0,
        _ => 400.0,
    }
}

/// Acquire the public component-three daily rate, adding SPEED only here.
fn speed(sample: &mut Sample<'_>, jd: f64, flags: i32) -> Result<f64, Error> {
    sample(jd, flags | FLG_SPEED).map(|p| p[3])
}

/// Sun/Moon have no retrograde phase. Numerical/range failures produce false;
/// lookup, configuration and data-acquisition errors propagate.
pub(super) fn retrograde(
    sample: &mut Sample<'_>,
    body: i32,
    jd: f64,
    flags: i32,
) -> Result<bool, Error> {
    if matches!(body, SUN | MOON) {
        return Ok(false);
    }
    match speed(sample, jd, flags) {
        Ok(v) => Ok(v < 0.0),
        Err(e) if super::super::crossings::station_position_error(&e) => Ok(false),
        Err(e) => Err(e),
    }
}

/// Motion kind at the requested UT1. Near means |speed| < 0.01; the previous
/// and following one-day samples determine decreasing versus increasing speed.
/// This trend is distinct from the sign-change finder classification.
pub(super) fn motion(
    sample: &mut Sample<'_>,
    body: i32,
    jd: f64,
    flags: i32,
) -> Result<StationMotionKind, Error> {
    if matches!(body, SUN | MOON) {
        return Ok(StationMotionKind::Direct);
    }
    let result: Result<StationMotionKind, Error> = (|| {
        let v = speed(sample, jd, flags)?;
        if v.abs() < 0.01 {
            let before = speed(sample, jd - 1.0, flags)?;
            let after = speed(sample, jd + 1.0, flags)?;
            Ok(if before > after {
                StationMotionKind::StationaryRetrograde
            } else {
                StationMotionKind::StationaryDirect
            })
        } else {
            Ok(if v < 0.0 {
                StationMotionKind::Retrograde
            } else {
                StationMotionKind::Direct
            })
        }
    })();
    match result {
        Err(e) if super::super::crossings::station_position_error(&e) => {
            Ok(StationMotionKind::Direct)
        }
        other => other,
    }
}

/// Floating division at the consumed arithmetic stage. A zero divisor refuses
/// after its operands have been acquired, including signed zero. Other IEEE
/// values follow binary64 division. A zero denominator is a typed error.
fn divide(numerator: f64, denominator: f64) -> Result<f64, Error> {
    if denominator == 0.0 {
        Err(Error::RationalDivisionByZero)
    } else {
        Ok(numerator / denominator)
    }
}

/// Refine the rate-zero bracket with inverse quadratic interpolation, secant
/// steps and conditional bisection. The stopping rate is 1e-6 in the selected
/// component units/day; a 1e-8-day time guard controls bisection. Return the
/// latest b after 100 iterations; no extra provider call or cap refusal occurs.
/// Endpoint reads and subsequent arithmetic failures propagate to the caller.
fn brent(
    sample: &mut Sample<'_>,
    body: i32,
    mut a: f64,
    mut b: f64,
    flags: i32,
) -> Result<f64, Error> {
    let mut fa = speed(sample, a, flags)?;
    let mut fb = speed(sample, b, flags)?;
    if fa * fb > 0.0 {
        return Err(Error::NoStationInInterval {
            lo: a,
            hi: b,
            planet: body,
        });
    }
    if fa.abs() < fb.abs() {
        std::mem::swap(&mut a, &mut b);
        std::mem::swap(&mut fa, &mut fb);
    }
    let (mut c, mut fc, mut mflag, mut d) = (a, fa, true, 0.0);
    for _ in 0..100 {
        if fb.abs() < 1e-6 {
            return Ok(b);
        }
        let mut s = if fa != fc && fb != fc {
            divide(a * fb * fc, (fa - fb) * (fa - fc))?
                + divide(b * fa * fc, (fb - fa) * (fb - fc))?
                + divide(c * fa * fb, (fc - fa) * (fc - fb))?
        } else {
            b - divide(fb * (b - a), fb - fa)?
        };
        let cone = (3.0 * a + b) / 4.0;
        let cond1 = !(cone < s && s < b || b < s && s < cone);
        let cond2 = mflag && (s - b).abs() >= (b - c).abs() / 2.0;
        let cond3 = !mflag && (s - b).abs() >= (c - d).abs() / 2.0;
        let cond4 = mflag && (b - c).abs() < 1e-8;
        let cond5 = !mflag && (c - d).abs() < 1e-8;
        if cond1 || cond2 || cond3 || cond4 || cond5 {
            s = (a + b) / 2.0;
            mflag = true;
        } else {
            mflag = false;
        }
        let fs = speed(sample, s, flags)?;
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

/// First rate sign-change in 50 equal forward intervals. The last epoch is
/// lo + 50*step, retaining that expression instead of replacing it by hi.
fn bracket(
    sample: &mut Sample<'_>,
    body: i32,
    lo: f64,
    hi: f64,
    flags: i32,
) -> Result<(f64, f64), Error> {
    let step = (hi - lo) / 50.0;
    let (mut previous, mut v) = (lo, speed(sample, lo, flags)?);
    for i in 1..=50 {
        let jd = lo + f64::from(i) * step;
        let next = speed(sample, jd, flags)?;
        if v * next <= 0.0 {
            return Ok((previous, jd));
        }
        previous = jd;
        v = next;
    }
    Err(Error::NoStationInInterval {
        lo,
        hi,
        planet: body,
    })
}

/// Finder classification reads the minus-one-day sample before plus one day;
/// only positive-before/negative-after is SR, every other trend is SD.
fn classify(sample: &mut Sample<'_>, jd: f64, flags: i32) -> Result<StationKind, Error> {
    let before = speed(sample, jd - 1.0, flags)?;
    let after = speed(sample, jd + 1.0, flags)?;
    Ok(if before > 0.0 && after < 0.0 {
        StationKind::StationaryRetrograde
    } else {
        StationKind::StationaryDirect
    })
}

/// Forward search: one initial speed read, then four half-period attempts.
/// Public Error/RuntimeError recovery spans scan, refinement and classification;
/// builtin ValueError, arithmetic and loader I/O remain terminal. A wrong-kind
/// root advances by one day; a caught attempt advances by its window. No
/// numerical preflight or rewritten endpoint changes the consumed samples.
pub(super) fn find(
    sample: &mut Sample<'_>,
    body: i32,
    jd: f64,
    request: StationRequest,
    flags: i32,
) -> Result<(f64, StationKind), Error> {
    check_body(body, false)?;
    let window = period(body) / 2.0;
    speed(sample, jd, flags)?;
    let mut start = jd;
    for _ in 0..4 {
        let result = (|| {
            let (a, b) = bracket(sample, body, start, start + window, flags)?;
            let root = brent(sample, body, a, b, flags)?;
            Ok((root, classify(sample, root, flags)?))
        })();
        match result {
            Ok((root, kind)) => {
                if crate::events::station_matches_request(kind, request) {
                    return Ok((root, kind));
                }
                start = root + 1.0;
            }
            Err(e) if super::super::crossings::station_recoverable(&e) => start += window,
            Err(e) => return Err(e),
        }
    }
    Err(Error::NoStationInInterval {
        lo: jd,
        hi: start + window,
        planet: body,
    })
}

/// Most recent station, scanning one period backward with max(60,int(P/2))
/// intervals. Read the original epoch before the descending scan; there is no
/// local recovery. A station at the initial endpoint is retained when bracketed.
fn previous(
    sample: &mut Sample<'_>,
    body: i32,
    jd: f64,
    flags: i32,
) -> Result<(f64, StationKind), Error> {
    check_body(body, false)?;
    let window = period(body);
    let count = 60.max((window / 2.0) as i32);
    let step = window / f64::from(count);
    let (mut next_jd, mut next_v) = (jd, speed(sample, jd, flags)?);
    for i in 1..=count {
        let t = jd - f64::from(i) * step;
        let v = speed(sample, t, flags)?;
        if next_v * v <= 0.0 {
            let root = brent(sample, body, t, next_jd, flags)?;
            return Ok((root, classify(sample, root, flags)?));
        }
        next_jd = t;
        next_v = v;
    }
    Err(Error::NoStationInInterval {
        lo: jd - window,
        hi: jd,
        planet: body,
    })
}

/// Next complete retrograde period. The predicate is consumed before the
/// station-body guard; an already retrograde body first finds the current SD.
pub(super) fn next_period(
    sample: &mut Sample<'_>,
    body: i32,
    jd: f64,
    flags: i32,
) -> Result<(f64, f64), Error> {
    let start = if retrograde(sample, body, jd, flags)? {
        find(sample, body, jd, StationRequest::Direct, flags)?.0 + 1.0
    } else {
        jd
    };
    let sr = find(sample, body, start, StationRequest::Retrograde, flags)?.0;
    let sd = find(sample, body, sr + 1.0, StationRequest::Direct, flags)?.0;
    Ok((sr, sd))
}

/// Central difference of actual public components with the original flags.
/// No SPEED bit is inserted. Acquire jd-dt then jd+dt before dividing by 2*dt;
/// negative dt is accepted, zero division is a consumed arithmetic refusal.
/// Component zero receives one +/-360 fold; other coordinate shapes still
/// supply their actual components and are not converted to physical longitude.
pub(super) fn velocity(
    sample: &mut Sample<'_>,
    jd: f64,
    flags: i32,
    dt: f64,
) -> Result<(f64, f64, f64), Error> {
    let before = sample(jd - dt, flags)?;
    let after = sample(jd + dt, flags)?;
    let mut lon = after[0] - before[0];
    if lon > 180.0 {
        lon -= 360.0;
    } else if lon < -180.0 {
        lon += 360.0;
    }
    Ok((
        divide(lon, 2.0 * dt)?,
        divide(after[1] - before[1], 2.0 * dt)?,
        divide(after[2] - before[2], 2.0 * dt)?,
    ))
}

/// Current rate, forward attempt, backward attempt, then selected longitude.
/// Keep the nearest root with forward winning a tie. Double caught failures
/// return a partial success owning the first error. Final longitude uses the
/// caller's original flags, not the SPEED flags of rate samples.
pub(super) fn info(
    sample: &mut Sample<'_>,
    body: i32,
    jd: f64,
    flags: i32,
) -> Result<StationInformation, Error> {
    check_body(body, true)?;
    let current = speed(sample, jd, flags)?;
    let mut first_error = None;
    let mut candidate = match find(sample, body, jd, StationRequest::Any, flags) {
        Ok(v) => Some(v),
        Err(e) if super::super::crossings::station_recoverable(&e) => {
            first_error = Some(e);
            None
        }
        Err(e) => return Err(e),
    };
    match previous(sample, body, jd, flags) {
        Ok(v) => {
            if candidate.is_none_or(|old| (v.0 - jd).abs() < (old.0 - jd).abs()) {
                candidate = Some(v);
            }
        }
        Err(e) if super::super::crossings::station_recoverable(&e) => {
            if first_error.is_none() {
                first_error = Some(e);
            }
        }
        Err(e) => return Err(e),
    }
    match candidate {
        Some((root, kind)) => {
            let longitude = sample(root, flags)?[0];
            Ok(StationInformation::Complete(crate::events::StationInfo {
                jd_station: root,
                kind,
                days_to_station: root - jd,
                longitude_at_station_deg: longitude,
                currently_retrograde: current < 0.0,
                speed_now_deg_day: current,
            }))
        }
        None => Ok(StationInformation::Unavailable {
            currently_retrograde: current < 0.0,
            velocity: current,
            error: first_error.expect("both unsuccessful searches own a caught error"),
        }),
    }
}
