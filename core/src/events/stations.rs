//! Stations implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Retrograde predicate over a longitude speed. Sun/Moon identity lives with
/// the caller — see [`check_station_body`].
pub fn is_retrograde_speed(speed_deg_day: f64) -> bool {
    speed_deg_day < 0.0
}

/// Refuse station searches for bodies that never station.
pub fn check_station_body(planet_id: i32) -> Result<()> {
    if planet_id == 0 || planet_id == 1 {
        return Err(Error::BodyHasNoStations { body: planet_id });
    }
    Ok(())
}

/// Classify a station from the velocity trend across it.
pub fn classify_station(speed_before_deg_day: f64, speed_after_deg_day: f64) -> StationKind {
    if speed_before_deg_day > 0.0 && speed_after_deg_day < 0.0 {
        StationKind::StationaryRetrograde
    } else {
        StationKind::StationaryDirect
    }
}

/// Whether a found station satisfies the request.
pub fn station_matches_request(kind: StationKind, request: StationRequest) -> bool {
    match request {
        StationRequest::Any => true,
        StationRequest::Retrograde => kind == StationKind::StationaryRetrograde,
        StationRequest::Direct => kind == StationKind::StationaryDirect,
    }
}

/// Evaluate the speed provider once, rejecting non-finite velocities.
pub(super) fn checked_speed(
    get_speed: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    jd: f64,
) -> Result<f64> {
    let speed = get_speed(jd)?;
    if !speed.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: jd });
    }
    Ok(speed)
}

/// Refine a bracketed station with Brent's method over the speed residual
/// An unbracketed interval is
/// `NoStationInInterval`; a non-finite interpolant or a non-finite speed
/// evaluation is an error. Cap exhaustion returns the best estimate,
/// mirroring the release tail (and the `brent_crossing` root-finder precedent).
pub fn brent_station(
    get_speed: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    planet_id: i32,
    a: f64,
    b: f64,
    tolerance_deg_day: f64,
    max_iter: usize,
) -> Result<f64> {
    if !a.is_finite() || !b.is_finite() || a == b {
        return Err(Error::InvalidSearchInterval { a, b });
    }
    if !tolerance_deg_day.is_finite() || tolerance_deg_day <= 0.0 {
        return Err(Error::InvalidSearchTolerance {
            value: tolerance_deg_day,
        });
    }
    if max_iter == 0 {
        return Err(Error::InvalidSearchIterations { max_iter });
    }
    let mut xa = a;
    let mut xb = b;
    let mut fa = checked_speed(get_speed, xa)?;
    let mut fb = checked_speed(get_speed, xb)?;
    if fa * fb > 0.0 {
        return Err(Error::NoStationInInterval {
            lo: a.min(b),
            hi: a.max(b),
            planet: planet_id,
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
        if fb.abs() < tolerance_deg_day {
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
            return Err(Error::NonFiniteSearchValue { x: s });
        }
        let fs = checked_speed(get_speed, s)?;
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

/// Isolate a station bracket `[jd_a, jd_b]` by uniform scan: the first velocity sign change, endpoints
/// included (`<=`), wins. No sign change is `NoStationInInterval`.
pub fn find_station_bracket(
    get_speed: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    planet_id: i32,
    lo: f64,
    hi: f64,
    samples: usize,
) -> Result<(f64, f64)> {
    if !lo.is_finite() || !hi.is_finite() || lo >= hi {
        return Err(Error::InvalidSearchInterval { a: lo, b: hi });
    }
    if samples == 0 {
        return Err(Error::InvalidScanSamples { count: 0 });
    }
    let step = (hi - lo) / samples as f64;
    let mut jd_prev = lo;
    let mut speed_prev = checked_speed(get_speed, lo)?;
    for i in 1..=samples {
        let jd_curr = if i == samples {
            hi
        } else {
            lo + i as f64 * step
        };
        let speed_curr = checked_speed(get_speed, jd_curr)?;
        if speed_prev * speed_curr <= 0.0 {
            return Ok((jd_prev, jd_curr));
        }
        jd_prev = jd_curr;
        speed_prev = speed_curr;
    }
    Err(Error::NoStationInInterval {
        lo,
        hi,
        planet: planet_id,
    })
}

/// Next station at or after `start` matching `request`.
///
/// `get_speed(jd)` returns longitude speed in degrees/day for `planet_id` in
/// one selected frame; input/output dates share its Julian-day clock. This
/// reducer performs no TT/UT1 conversion. Returns `(jd_station, kind)` after
/// sign-change refinement and before/after classification. Unsupported station
/// bodies and a nonfinite start are typed errors.
///
/// Every bracket-search error currently advances to the next search window,
/// including a callback error; exhausting those attempts returns
/// `NoStationInInterval`. Once a bracket is accepted, refinement and
/// classification errors propagate. Thus exhaustion is not proof that every
/// requested epoch had usable provider data. This error policy is distinct
/// from the backward scanner's immediate propagation of sample errors.
pub fn find_station(
    get_speed: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    planet_id: i32,
    start: f64,
    request: StationRequest,
) -> Result<(f64, StationKind)> {
    check_station_body(planet_id)?;
    if !start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start });
    }
    let search_window = synodic_days(planet_id) / 2.0;
    let mut jd_search_start = start;
    for _ in 0..STATION_SEARCH_ATTEMPTS {
        let jd_search_end = jd_search_start + search_window;
        match find_station_bracket(
            get_speed,
            planet_id,
            jd_search_start,
            jd_search_end,
            STATION_BRACKET_SAMPLES,
        ) {
            Ok((jd_a, jd_b)) => {
                let jd_station = brent_station(
                    get_speed,
                    planet_id,
                    jd_a,
                    jd_b,
                    STATION_VELOCITY_TOLERANCE_DEG_DAY,
                    100,
                )?;
                let speed_before = checked_speed(get_speed, jd_station - STATION_TYPE_OFFSET_DAYS)?;
                let speed_after = checked_speed(get_speed, jd_station + STATION_TYPE_OFFSET_DAYS)?;
                let found = classify_station(speed_before, speed_after);
                if station_matches_request(found, request) {
                    return Ok((jd_station, found));
                }
                jd_search_start = jd_station + STATION_TYPE_OFFSET_DAYS;
            }
            Err(_) => {
                jd_search_start += search_window;
            }
        }
    }
    Err(Error::NoStationInInterval {
        lo: start,
        hi: start + STATION_SEARCH_ATTEMPTS as f64 * search_window,
        planet: planet_id,
    })
}

/// Search backward from `start` for the nearest station.
///
/// Provider units and clock ownership follow [`find_station`]. A full-synodic
/// backward scan at roughly two-day resolution refines its first sign change
/// with Brent, then classifies samples on either side. An exact zero at the
/// starting endpoint can return `start`; this low-level search does not enforce
/// strict exclusion of that endpoint. Callback/refinement errors propagate.
pub fn find_previous_station(
    get_speed: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    planet_id: i32,
    start: f64,
) -> Result<(f64, StationKind)> {
    check_station_body(planet_id)?;
    if !start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start });
    }
    let search_window = synodic_days(planet_id);
    let num_samples = 60.max((search_window / 2.0) as usize);
    let step = search_window / num_samples as f64;
    let mut jd_next = start;
    let mut speed_next = checked_speed(get_speed, start)?;
    for i in 1..=num_samples {
        let jd_curr = start - i as f64 * step;
        let speed_curr = checked_speed(get_speed, jd_curr)?;
        if speed_next * speed_curr <= 0.0 {
            let jd_station = brent_station(
                get_speed,
                planet_id,
                jd_curr,
                jd_next,
                STATION_VELOCITY_TOLERANCE_DEG_DAY,
                100,
            )?;
            let speed_before = checked_speed(get_speed, jd_station - STATION_TYPE_OFFSET_DAYS)?;
            let speed_after = checked_speed(get_speed, jd_station + STATION_TYPE_OFFSET_DAYS)?;
            return Ok((jd_station, classify_station(speed_before, speed_after)));
        }
        jd_next = jd_curr;
        speed_next = speed_curr;
    }
    Err(Error::NoStationInInterval {
        lo: start - search_window,
        hi: start,
        planet: planet_id,
    })
}

/// Nearest station among the forward and backward search candidates.
///
/// Dates and callback speed units follow [`find_station`]. An endpoint station
/// can be selected. Absolute time distance chooses the winner, with a forward
/// candidate retained on an exact tie; `days_to_station` is signed. One failed
/// direction does not discard a successful candidate from the other. If both
/// fail, return the first error rather than a partial record.
pub fn nearest_station(
    get_speed: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    planet_id: i32,
    start: f64,
) -> Result<NearestStation> {
    check_station_body(planet_id)?;
    let speed_now = checked_speed(get_speed, start)?;
    let mut best: Option<(f64, StationKind)> = None;
    let mut first_err: Option<Error> = None;
    match find_station(get_speed, planet_id, start, StationRequest::Any) {
        Ok(candidate) => best = Some(candidate),
        Err(e) => first_err = Some(e),
    }
    match find_previous_station(get_speed, planet_id, start) {
        Ok(candidate) => {
            let replace = match best {
                None => true,
                Some((jd_best, _)) => (candidate.0 - start).abs() < (jd_best - start).abs(),
            };
            if replace {
                best = Some(candidate);
            }
        }
        Err(e) => {
            if first_err.is_none() {
                first_err = Some(e);
            }
        }
    }
    match best {
        Some((jd_station, kind)) => Ok(NearestStation {
            jd_station,
            kind,
            days_to_station: jd_station - start,
            currently_retrograde: is_retrograde_speed(speed_now),
            speed_now_deg_day: speed_now,
        }),
        None => Err(first_err.unwrap_or(Error::NoStationInInterval {
            lo: start,
            hi: start,
            planet: planet_id,
        })),
    }
}

/// Next complete retrograde interval `(station_retrograde, station_direct)`.
///
/// Dates and callback speed units follow [`find_station`]. If the body is
/// already retrograde at `start`, skip the remainder of that interval and
/// return the next complete one. Consumed station-search errors retain their
/// underlying forward-search retry/propagation policy.
pub fn next_retrograde(
    get_speed: &mut impl FnMut(f64) -> std::result::Result<f64, Error>,
    planet_id: i32,
    start: f64,
) -> Result<(f64, f64)> {
    check_station_body(planet_id)?;
    let speed_now = checked_speed(get_speed, start)?;
    if is_retrograde_speed(speed_now) {
        let (jd_sd, _) = find_station(get_speed, planet_id, start, StationRequest::Direct)?;
        let (jd_sr, _) = find_station(
            get_speed,
            planet_id,
            jd_sd + STATION_TYPE_OFFSET_DAYS,
            StationRequest::Retrograde,
        )?;
        let (jd_sd_next, _) = find_station(
            get_speed,
            planet_id,
            jd_sr + STATION_TYPE_OFFSET_DAYS,
            StationRequest::Direct,
        )?;
        Ok((jd_sr, jd_sd_next))
    } else {
        let (jd_sr, _) = find_station(get_speed, planet_id, start, StationRequest::Retrograde)?;
        let (jd_sd, _) = find_station(
            get_speed,
            planet_id,
            jd_sr + STATION_TYPE_OFFSET_DAYS,
            StationRequest::Direct,
        )?;
        Ok((jd_sr, jd_sd))
    }
}

/// Stable velocity near a station by central difference: positions at `jd ± dt`, longitude
/// difference wrapped into `(-180, 180]`, divided by `2·dt`. Returns
/// `(velocity_lon, velocity_lat, velocity_dist)` in degrees/day,
/// degrees/day and distance/day. The provider's distance unit is its own
/// (AU in the release binding).
pub fn station_central_velocity(
    get_state: &mut impl FnMut(f64) -> std::result::Result<(f64, f64, f64), Error>,
    jd_station: f64,
    dt_days: f64,
) -> Result<(f64, f64, f64)> {
    if !jd_station.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: jd_station });
    }
    if !dt_days.is_finite() || dt_days <= 0.0 {
        return Err(Error::InvalidSearchStep { value: dt_days });
    }
    let checked = &mut |jd: f64| -> Result<(f64, f64, f64)> {
        let (lon, lat, dist) = get_state(jd)?;
        if !lon.is_finite() || !lat.is_finite() || !dist.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok((lon, lat, dist))
    };
    let (lon_before, lat_before, dist_before) = checked(jd_station - dt_days)?;
    let (lon_after, lat_after, dist_after) = checked(jd_station + dt_days)?;
    let mut dlon = lon_after - lon_before;
    if dlon > 180.0 {
        dlon -= 360.0;
    } else if dlon < -180.0 {
        dlon += 360.0;
    }
    let span = 2.0 * dt_days;
    Ok((
        dlon / span,
        (lat_after - lat_before) / span,
        (dist_after - dist_before) / span,
    ))
}

/// Next (or, for [`SearchDirection::Backward`], strictly previous) lunar
/// node crossing over an abstract `(lon, lat, lat_speed)` provider
/// Returns the crossing time with the longitude and
/// latitude there. The `calc` provider binding lives at the serving layer; the
/// UT↔TT wrapping around it is binding-layer work.
pub fn moon_node_crossing(
    get_moon: &mut impl FnMut(f64) -> std::result::Result<(f64, f64, f64), Error>,
    start: f64,
    direction: SearchDirection,
) -> Result<MoonNodeCrossing> {
    if !start.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start });
    }
    let checked = &mut |jd: f64| -> Result<(f64, f64, f64)> {
        let (lon, lat, lat_speed) = get_moon(jd)?;
        if !lon.is_finite() || !lat.is_finite() || !lat_speed.is_finite() {
            return Err(Error::NonFiniteSearchValue { x: jd });
        }
        Ok((lon, lat, lat_speed))
    };
    let dir: f64 = match direction {
        SearchDirection::Forward => 1.0,
        SearchDirection::Backward => -1.0,
    };
    let (_, lat_start, lat_speed_start) = checked(start)?;

    let mut accept_eps = 1e-6;
    let mut jd_guess = start + HALF_NODAL_MONTH_DAYS * dir;
    if lat_start.abs() < 10.0 * MOON_NODE_TOLERANCE_DEG {
        let t_est = if lat_speed_start.abs() > 0.1 {
            -lat_start / lat_speed_start
        } else {
            0.0
        };
        if t_est * dir > MOON_NODE_AT_NODE_GATE_DAYS {
            jd_guess = start + t_est;
            accept_eps = MOON_NODE_AT_NODE_GATE_DAYS;
        }
    } else {
        // Skip-free scan: successive node crossings are >= ~6.8 days
        // apart, so a 2-day step cannot jump a root. The first sign
        // change brackets the first crossing in the requested direction;
        // bisection refines the bracket before Newton so an edge-sited
        // root cannot seed the wrong basin.
        let scan_step = 2.0;
        let scan_span = HALF_NODAL_MONTH_DAYS + 4.0;
        let n_steps = (scan_span / scan_step) as usize + 1;
        let mut lat_prev = lat_start;
        let mut t_prev = start;
        for i in 1..=n_steps {
            let t_cur = start + i as f64 * scan_step * dir;
            let lat_cur = checked(t_cur)?.1;
            if lat_prev * lat_cur < 0.0 {
                let (mut lo, mut hi) = (t_prev, t_cur);
                let mut f_lo = lat_prev;
                for _ in 0..60 {
                    let mid = 0.5 * (lo + hi);
                    let f_mid = checked(mid)?.1;
                    if f_lo * f_mid <= 0.0 {
                        hi = mid;
                    } else {
                        lo = mid;
                        f_lo = f_mid;
                    }
                }
                jd_guess = 0.5 * (lo + hi);
                break;
            }
            t_prev = t_cur;
            lat_prev = lat_cur;
        }
    }

    let mut jd = jd_guess;
    let mut best = (jd_guess, lat_start);
    for _ in 0..MOON_NODE_MAX_ITER {
        let (lon, lat, mut lat_speed) = checked(jd)?;
        best = (jd, lat);
        if lat.abs() < MOON_NODE_TOLERANCE_DEG {
            // Strictly past the start in the requested direction; the
            // epsilon accepts legitimate seconds-away results yet rejects
            // the numerically identical "same event" at the start itself.
            let past = match direction {
                SearchDirection::Forward => jd > start + accept_eps,
                SearchDirection::Backward => jd < start - accept_eps,
            };
            if past {
                return Ok(MoonNodeCrossing {
                    jd,
                    longitude_deg: lon,
                    latitude_deg: lat,
                });
            }
            jd += (HALF_NODAL_MONTH_DAYS / 2.0) * dir;
            continue;
        }
        if lat_speed.abs() < 0.1 {
            lat_speed = if lat >= 0.0 { 1.0 } else { -1.0 };
        }
        jd -= lat / lat_speed;
        match direction {
            SearchDirection::Forward => {
                if jd < start {
                    jd = start + HALF_NODAL_MONTH_DAYS / 2.0;
                }
            }
            SearchDirection::Backward => {
                if jd > start {
                    jd = start - HALF_NODAL_MONTH_DAYS / 2.0;
                }
            }
        }
        if (jd - start).abs() > MOON_NODE_DIVERGENCE_DAYS {
            return Err(Error::MoonNodeDiverged {
                start,
                estimate: jd,
            });
        }
    }
    Err(Error::MoonNodeNotConverged {
        best: best.0,
        latitude: best.1,
    })
}

/// Next station of the bound body, returning `(jd_ut1, station_kind)`.
///
/// `planet_id` controls admission/search heuristics and must match the body's
/// ID retained in `serving`; that equality is not checked here. XYZ/RADIANS are
/// refused. Other frame flags and [`find_station`]'s bracket-error retry policy
/// remain unchanged; a no-station error is not a complete-coverage certificate.
pub fn find_station_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    start_ut: f64,
    request: StationRequest,
) -> Result<(f64, StationKind)> {
    check_angular_serving_flags(serving.calc_flags(), "find_station_ut")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_speed = |jd: f64| serving.speed_at_ut(jd);
    find_station(&mut get_speed, planet_id, start_ut, request)
}

/// Backward station candidate at a UT1 Julian-day start.
///
/// `planet_id` must match the bound serving body. XYZ/RADIANS are refused;
/// endpoint inclusion and errors follow [`find_previous_station`], including
/// the possibility of returning an exact station at `start_ut`.
pub fn find_previous_station_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    start_ut: f64,
) -> Result<(f64, StationKind)> {
    check_angular_serving_flags(serving.calc_flags(), "find_previous_station_ut")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_speed = |jd: f64| serving.speed_at_ut(jd);
    find_previous_station(&mut get_speed, planet_id, start_ut)
}

/// Nearest station in either direction over the sealed serving.
pub fn nearest_station_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    start_ut: f64,
) -> Result<NearestStation> {
    check_angular_serving_flags(serving.calc_flags(), "nearest_station_ut")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_speed = |jd: f64| serving.speed_at_ut(jd);
    nearest_station(&mut get_speed, planet_id, start_ut)
}

/// Next retrograde `(start_SR, end_SD)` pair over the sealed serving
pub fn next_retrograde_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    start_ut: f64,
) -> Result<(f64, f64)> {
    check_angular_serving_flags(serving.calc_flags(), "next_retrograde_ut")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_speed = |jd: f64| serving.speed_at_ut(jd);
    next_retrograde(&mut get_speed, planet_id, start_ut)
}

impl StationMotionKind {
    /// Sweep wire code: 0 = direct, 1 = retrograde,
    /// 2 = stationary-retrograde, 3 = stationary-direct.
    pub fn code(self) -> i32 {
        match self {
            StationMotionKind::Direct => 0,
            StationMotionKind::Retrograde => 1,
            StationMotionKind::StationaryRetrograde => 2,
            StationMotionKind::StationaryDirect => 3,
        }
    }

    /// Release spelling of the kind.
    pub fn as_str(self) -> &'static str {
        match self {
            StationMotionKind::Direct => "direct",
            StationMotionKind::Retrograde => "retrograde",
            StationMotionKind::StationaryRetrograde => "stationary_retrograde",
            StationMotionKind::StationaryDirect => "stationary_direct",
        }
    }
}

/// Whether the bound body is retrograde at a UT epoch over the sealed
/// serving. Retrograde means a
/// negative longitude speed (deg/day) from the sealed triple.
pub fn is_retrograde_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    jd_ut: f64,
) -> Result<bool> {
    check_angular_serving_flags(serving.calc_flags(), "is_retrograde")?;
    if planet_id == crate::constants::SUN || planet_id == crate::constants::MOON {
        return Ok(false);
    }
    Ok(serving
        .speed_at_ut(jd_ut)
        .map(is_retrograde_speed)
        .unwrap_or(false))
}

/// Station-motion kind at a UT epoch over the sealed serving (same
/// acquisition as release `get_station_type`: Sun/Moon read `"direct"`
/// without consulting the ephemeris, and every provider failure folds to
/// `"direct"`. The near-station trend rule compares the ±1-day speeds —
/// decreasing reads stationary-retrograde, anything else
/// stationary-direct. Note this trend rule differs from
/// [`classify_station`] (the `find_station_ut` sign rule) by design: it
/// answers "which way is the speed heading", not "which transition was
/// crossed").
pub fn station_type_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    jd_ut: f64,
) -> Result<StationMotionKind> {
    check_angular_serving_flags(serving.calc_flags(), "station_type")?;
    if planet_id == crate::constants::SUN || planet_id == crate::constants::MOON {
        return Ok(StationMotionKind::Direct);
    }
    let speed = match serving.speed_at_ut(jd_ut) {
        Ok(v) => v,
        Err(_) => return Ok(StationMotionKind::Direct),
    };
    if is_near_station(speed) {
        let before = serving.speed_at_ut(jd_ut - STATION_TYPE_OFFSET_DAYS);
        let after = serving.speed_at_ut(jd_ut + STATION_TYPE_OFFSET_DAYS);
        match (before, after) {
            (Ok(b), Ok(a)) => Ok(if b > a {
                StationMotionKind::StationaryRetrograde
            } else {
                StationMotionKind::StationaryDirect
            }),
            _ => Ok(StationMotionKind::Direct),
        }
    } else if is_retrograde_speed(speed) {
        Ok(StationMotionKind::Retrograde)
    } else {
        Ok(StationMotionKind::Direct)
    }
}

/// Stable velocity near a station by central difference over the sealed
/// serving (same acquisition as release `calc_velocity_at_station`:
/// positions at `jd_station ± dt_days` through the plain `calc_ut` path —
/// mirrored here by the position triple, not the speed stencil — with
/// the longitude difference wrapped into `(-180, 180]`; the
/// [`station_central_velocity`] station core). The returned triple is
/// `(velocity_lon, velocity_lat, velocity_dist)` in deg/day, deg/day and
/// AU/day; `dt_days` is the half-stencil in days.
/// Like the release entry there is no Sun/Moon guard: any served body
/// differentiates.
pub fn velocity_at_station_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    jd_station: f64,
    dt_days: f64,
) -> Result<(f64, f64, f64)> {
    check_angular_serving_flags(serving.calc_flags(), "velocity_at_station")?;
    let mut get_state = |jd: f64| serving.lonlatdist_at_ut(jd);
    station_central_velocity(&mut get_state, jd_station, dt_days)
}

/// Nearest-station summary over the sealed serving.
pub fn station_info_ut_serving<P: StateProvider>(
    serving: &SealedLongitudeServing<'_, P>,
    planet_id: i32,
    start_ut: f64,
) -> Result<StationInfo> {
    check_angular_serving_flags(serving.calc_flags(), "station_info")?;
    if !start_ut.is_finite() {
        return Err(Error::NonFiniteSearchValue { x: start_ut });
    }
    let mut get_speed = |jd: f64| serving.speed_at_ut(jd);
    let near = nearest_station(&mut get_speed, planet_id, start_ut)?;
    let (lon, _, _) = serving.lonlatdist_at_ut(near.jd_station)?;
    Ok(StationInfo {
        jd_station: near.jd_station,
        kind: near.kind,
        days_to_station: near.days_to_station,
        longitude_at_station_deg: lon,
        currently_retrograde: near.currently_retrograde,
        speed_now_deg_day: near.speed_now_deg_day,
    })
}
