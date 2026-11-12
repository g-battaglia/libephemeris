//! Frames implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Resolve center flags with topocentric, barycentric, then heliocentric
/// priority. The returned bits select the physical observer; callers retain
/// their independent returned-flag policy.
pub fn resolve_center_flags(flags: i32) -> i32 {
    if flags & FLG_TOPOCTR != 0 {
        flags & !(FLG_HELCTR | FLG_BARYCTR)
    } else if flags & FLG_BARYCTR != 0 {
        flags & !FLG_HELCTR
    } else {
        flags
    }
}

/// Map (already priority-resolved) center flags to the geometric
/// observer. Resolution is idempotent: resolving again changes
/// nothing.
pub fn observation_center(resolved_flags: i32) -> ObservationCenter {
    if resolved_flags & FLG_TOPOCTR != 0 {
        ObservationCenter::Topocentric
    } else if resolved_flags & FLG_BARYCTR != 0 {
        ObservationCenter::Barycentric
    } else if resolved_flags & FLG_HELCTR != 0 {
        ObservationCenter::Heliocentric
    } else {
        ObservationCenter::Geocentric
    }
}

/// Whether a geocentrically defined lunar point has no heliocentric or
/// barycentric place. Such requests return the public zero sentinel.
pub fn is_helio_bary_sentinel_point(body_id: i32) -> bool {
    body_id == MEAN_NODE
        || body_id == TRUE_NODE
        || body_id == MEAN_APOG
        || body_id == OSCU_APOG
        || body_id == INTP_APOG
        || body_id == INTP_PERG
}

impl ObserverRelativeState {
    /// Instantaneous target-minus-observer position (AU).
    pub fn geo_pos(self) -> [f64; 3] {
        self.geo_pos
    }

    /// Instantaneous target-minus-observer velocity (AU/day).
    pub fn geo_vel(self) -> [f64; 3] {
        self.geo_vel
    }

    /// Earth velocity driving the later aberration stage (AU/day):
    /// the observer velocity for geocentric/topocentric work, exactly
    /// zero for heliocentric/barycentric work (aberration skipped
    /// there). The aberration/deflection layer consumes this, never a re-derivation.
    pub fn aberration_vel(self) -> [f64; 3] {
        self.aberration_vel
    }

    /// True when the relative position is exactly the zero vector
    /// (degenerate self-observation or sentinel point). The entry layer uses
    /// this distance check to skip the sidereal ayanamsha subtraction
    /// on an undefined longitude.
    pub fn is_origin(self) -> bool {
        self.geo_pos == [0.0, 0.0, 0.0]
    }
}

/// Zero position and velocity for an undefined place at the observation
/// origin. All output components are zero, including the undefined longitude.
pub fn degenerate_origin() -> ObserverRelativeState {
    ObserverRelativeState {
        geo_pos: [0.0, 0.0, 0.0],
        geo_vel: [0.0, 0.0, 0.0],
        aberration_vel: [0.0, 0.0, 0.0],
    }
}

pub(super) fn check_finite_vec3(v: [f64; 3]) -> Result<()> {
    for &c in &v {
        if !c.is_finite() {
            return Err(Error::NonFiniteArgument { value: c });
        }
    }
    Ok(())
}

/// Select the observer's barycentric ICRS position and velocity.
///
/// Heliocentric work uses the Sun; barycentric work uses the SSB origin;
/// geocentric work uses Earth. Topocentric work adds the supplied site offset
/// to Earth. Inputs and outputs are AU and AU/day. The serving layer requires
/// a site for a real topocentric request.
pub fn observer_state_for_center(
    center: ObservationCenter,
    sun_pos: [f64; 3],
    sun_vel: [f64; 3],
    earth_pos: [f64; 3],
    earth_vel: [f64; 3],
    topo: Option<([f64; 3], [f64; 3])>,
) -> ([f64; 3], [f64; 3]) {
    match center {
        ObservationCenter::Heliocentric => (sun_pos, sun_vel),
        ObservationCenter::Barycentric => ([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]),
        ObservationCenter::Geocentric | ObservationCenter::Topocentric => {
            let (mut op, mut ov) = (earth_pos, earth_vel);
            if center == ObservationCenter::Topocentric {
                if let Some((tp, tv)) = topo {
                    op = [op[0] + tp[0], op[1] + tp[1], op[2] + tp[2]];
                    ov = [ov[0] + tv[0], ov[1] + tv[1], ov[2] + tv[2]];
                }
            }
            (op, ov)
        }
    }
}
/// Compose an observer-relative geometric state from barycentric ICRS states.
///
/// Target, Sun and Earth positions are AU, and velocities are AU/day. The
/// optional site offset is added only for a topocentric center. Every consumed
/// component must be finite; invalid inputs return [`Error::NonFiniteArgument`].
/// Identical target and observer states naturally subtract to zero.
#[allow(clippy::too_many_arguments)]
pub fn compose_relative_state(
    target_pos: [f64; 3],
    target_vel: [f64; 3],
    sun_pos: [f64; 3],
    sun_vel: [f64; 3],
    earth_pos: [f64; 3],
    earth_vel: [f64; 3],
    center: ObservationCenter,
    topo: Option<([f64; 3], [f64; 3])>,
) -> Result<ObserverRelativeState> {
    check_finite_vec3(target_pos)?;
    check_finite_vec3(target_vel)?;
    check_finite_vec3(sun_pos)?;
    check_finite_vec3(sun_vel)?;
    check_finite_vec3(earth_pos)?;
    check_finite_vec3(earth_vel)?;
    if let Some((tp, tv)) = topo {
        check_finite_vec3(tp)?;
        check_finite_vec3(tv)?;
    }
    let (obs_pos, obs_vel) =
        observer_state_for_center(center, sun_pos, sun_vel, earth_pos, earth_vel, topo);
    // The complete geocentric or topocentric observer velocity drives aberration.
    // Heliocentric and barycentric requests suppress this correction.
    let aberration_vel = match center {
        ObservationCenter::Heliocentric | ObservationCenter::Barycentric => [0.0, 0.0, 0.0],
        ObservationCenter::Geocentric | ObservationCenter::Topocentric => obs_vel,
    };
    // Barycentric observer is the SSB origin with zero velocity;
    // every other center differences against its observer state.
    let geo_pos = [
        target_pos[0] - obs_pos[0],
        target_pos[1] - obs_pos[1],
        target_pos[2] - obs_pos[2],
    ];
    // Subtract the observer velocity at the same epoch. Light-time derivatives
    // are composed separately after the target's retarded state is known.
    let geo_vel = [
        target_vel[0] - obs_vel[0],
        target_vel[1] - obs_vel[1],
        target_vel[2] - obs_vel[2],
    ];
    Ok(ObserverRelativeState {
        geo_pos,
        geo_vel,
        aberration_vel,
    })
}

/// Compose a body's geometric state with resolved center flags.
///
/// Sun-heliocentric, topocentric Earth and unsupported centers for abstract
/// lunar points return [`degenerate_origin`]. Other inputs enter
/// [`compose_relative_state`], which checks their finiteness. Positions and
/// velocities are AU and AU/day; no apparent or frame correction occurs here.
#[allow(clippy::too_many_arguments)]
pub fn compose_body_geometric(
    body_id: i32,
    flags: i32,
    target_pos: [f64; 3],
    target_vel: [f64; 3],
    sun_pos: [f64; 3],
    sun_vel: [f64; 3],
    earth_pos: [f64; 3],
    earth_vel: [f64; 3],
    topo: Option<([f64; 3], [f64; 3])>,
) -> Result<ObserverRelativeState> {
    let center = observation_center(resolve_center_flags(flags));
    match center {
        ObservationCenter::Heliocentric if body_id == SUN => Ok(degenerate_origin()),
        ObservationCenter::Topocentric if body_id == EARTH => Ok(degenerate_origin()),
        ObservationCenter::Heliocentric | ObservationCenter::Barycentric
            if is_helio_bary_sentinel_point(body_id) =>
        {
            Ok(degenerate_origin())
        }
        _ => compose_relative_state(
            target_pos, target_vel, sun_pos, sun_vel, earth_pos, earth_vel, center, topo,
        ),
    }
}

/// Compose the Moon's instantaneous geocentric state by subtracting Earth's
/// barycentric position and velocity. Units are AU and AU/day.
pub fn compose_moon_geocentric(
    moon_pos: [f64; 3],
    moon_vel: [f64; 3],
    earth_pos: [f64; 3],
    earth_vel: [f64; 3],
) -> Result<ObserverRelativeState> {
    compose_relative_state(
        moon_pos,
        moon_vel,
        [0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0],
        earth_pos,
        earth_vel,
        ObservationCenter::Geocentric,
        None,
    )
}

/// Whether `FLG_TRUEPOS` leaves light-time correction enabled.
pub fn wants_light_time(flags: i32) -> bool {
    flags & FLG_TRUEPOS == 0
}

/// Euclidean norm of a Cartesian vector, in the vector's length unit.
pub(crate) fn vec3_dist(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

impl LightTimeSolution {
    /// Retarded target-minus-observer position (AU): `target(jd_tt -
    /// lt) - observer(jd_tt)`, or the instantaneous difference when
    /// light-time is skipped.
    pub fn geo_pos(self) -> [f64; 3] {
        self.geo_pos
    }

    /// Light-time-corrected geometric velocity in AU per TT day: retarded target
    /// velocity times `1 − light_time_rate`, minus the observer velocity.
    pub fn geo_vel(self) -> [f64; 3] {
        self.geo_vel
    }

    /// Retarded target state (AU, AU/day) at [`LightTimeSolution::retarded_jd_tt`].
    pub fn retarded_state(self) -> ([f64; 3], [f64; 3]) {
        (self.retarded_pos, self.retarded_vel)
    }

    /// One-way light travel time in days (0 when skipped or degenerate).
    pub fn light_time_days(self) -> f64 {
        self.light_time_days
    }

    /// Derivative of the retained one-way delay, in days per TT day.
    /// The stored-polynomial solver propagates its actual finite iterations;
    /// the converged solvers use their implicit light-path velocity.
    pub fn light_time_rate(self) -> f64 {
        self.light_time_rate
    }

    /// Epoch the target was evaluated at (`jd_tt - lt`, or `jd_tt`
    /// when skipped).
    pub fn retarded_jd_tt(self) -> f64 {
        self.retarded_jd_tt
    }
}

/// Solve the one-way propagation equation by three fixed-point iterations.
///
/// `target_at` supplies the physical target center at a TT Julian Day in AU
/// and AU/day. The observer remains at the reception epoch. With `apply = false`
/// only the unretarded target is read, and the returned light time is zero.
/// Otherwise each iteration updates `lt = |target(t−lt)−observer(t)| / c`;
/// a zero ray ends the iteration. Non-finite states and provider errors
/// propagate, including coverage misses at an emitted epoch.
pub fn solve_light_time(
    target_at: impl Fn(f64) -> Result<([f64; 3], [f64; 3])>,
    jd_tt: f64,
    observer_pos: [f64; 3],
    observer_vel: [f64; 3],
    apply: bool,
) -> Result<LightTimeSolution> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    check_finite_vec3(observer_pos)?;
    check_finite_vec3(observer_vel)?;
    let (target_pos, target_vel) = target_at(jd_tt)?;
    check_finite_vec3(target_pos)?;
    check_finite_vec3(target_vel)?;
    let mut geo = [
        target_pos[0] - observer_pos[0],
        target_pos[1] - observer_pos[1],
        target_pos[2] - observer_pos[2],
    ];
    let mut retarded_pos = target_pos;
    let mut retarded_vel = target_vel;
    let mut lt = 0.0;
    let mut retarded_jd = jd_tt;
    if apply {
        for _ in 0..LIGHT_TIME_ITERATIONS {
            let dist = vec3_dist(geo);
            if dist == 0.0 {
                break;
            }
            lt = dist / C_LIGHT_AU_DAY;
            retarded_jd = jd_tt - lt;
            let (pos, vel) = target_at(retarded_jd)?;
            check_finite_vec3(pos)?;
            check_finite_vec3(vel)?;
            retarded_pos = pos;
            retarded_vel = vel;
            geo = [
                pos[0] - observer_pos[0],
                pos[1] - observer_pos[1],
                pos[2] - observer_pos[2],
            ];
        }
    }
    let geo_vel = light_time_corrected_velocity(retarded_vel, observer_vel, geo, lt, apply);
    Ok(LightTimeSolution {
        geo_pos: geo,
        geo_vel,
        retarded_pos,
        retarded_vel,
        light_time_days: lt,
        light_time_rate: if apply {
            super::planet_range::RangeState::from_cartesian(geo, geo_vel).rate / C_LIGHT_AU_DAY
        } else {
            0.0
        },
        retarded_jd_tt: retarded_jd,
    })
}

/// Solve one-way light time while retaining the observation's day/fraction split.
///
/// `whole` and `tt_fraction` are TT Julian-Day parts; `jd_tt` is their scalar
/// central epoch. Target states are AU and AU/day. The observer is fixed at
/// reception. Retarded reads pass through [`split_retarded_epoch_tt`](
/// super::reference_frames::split_retarded_epoch_tt), reducing cancellation
/// when subtracting a short delay from a large Julian Day.
///
/// Iteration uses [`SPLIT_LIGHT_TIME_ITERATIONS`] and
/// [`SPLIT_LIGHT_TIME_TOLERANCE_DAYS`], with [`super::FRAME_LIGHT_SPEED_AU_DAY`].
/// A geometric request reads only the unretarded state. With debug assertions
/// disabled, exhaustion returns the final iterate rather than a convergence
/// error. The returned velocity differentiates the propagation equation; this
/// helper does not establish an exact solution of that equation.
///
/// # Errors
///
/// Returns nonfinite-input/state errors and propagates callback failures,
/// including missing coverage at an emitted epoch. Finite components alone do
/// not guarantee that every intermediate norm or quotient remains finite.
///
/// # Panics
///
/// With debug assertions enabled, panics if an enabled light-time iteration
/// exhausts its limit without satisfying the settling criterion. A panic from
/// `target_at` also propagates; callbacks are invoked synchronously.
pub fn solve_light_time_split(
    target_at: impl Fn(f64) -> Result<([f64; 3], [f64; 3])>,
    whole: f64,
    tt_fraction: f64,
    jd_tt: f64,
    observer_pos: [f64; 3],
    observer_vel: [f64; 3],
    apply: bool,
) -> Result<LightTimeSolution> {
    if !whole.is_finite() {
        return Err(Error::NonFiniteArgument { value: whole });
    }
    if !tt_fraction.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt_fraction });
    }
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    check_finite_vec3(observer_pos)?;
    check_finite_vec3(observer_vel)?;
    let (target_pos, target_vel) = target_at(jd_tt)?;
    check_finite_vec3(target_pos)?;
    check_finite_vec3(target_vel)?;
    let mut geo = [
        target_pos[0] - observer_pos[0],
        target_pos[1] - observer_pos[1],
        target_pos[2] - observer_pos[2],
    ];
    let mut retarded_pos = target_pos;
    let mut retarded_vel = target_vel;
    let mut lt = 0.0;
    let mut lt_prev = 0.0;
    let mut retarded_jd = jd_tt;
    let mut settled = !apply;
    if apply {
        for _ in 0..super::SPLIT_LIGHT_TIME_ITERATIONS {
            let dist = vec3_dist(geo);
            lt = dist / super::FRAME_LIGHT_SPEED_AU_DAY;
            if (lt - lt_prev).abs() < super::SPLIT_LIGHT_TIME_TOLERANCE_DAYS {
                settled = true;
                break;
            }
            retarded_jd = super::reference_frames::split_retarded_epoch_tt(whole, tt_fraction, lt)?;
            let (pos, vel) = target_at(retarded_jd)?;
            check_finite_vec3(pos)?;
            check_finite_vec3(vel)?;
            retarded_pos = pos;
            retarded_vel = vel;
            geo = [
                pos[0] - observer_pos[0],
                pos[1] - observer_pos[1],
                pos[2] - observer_pos[2],
            ];
            lt_prev = lt;
        }
        debug_assert!(settled, "split light-time loop exhausted without settling");
    }
    let geo_vel = light_time_corrected_velocity(retarded_vel, observer_vel, geo, lt, apply);
    Ok(LightTimeSolution {
        geo_pos: geo,
        geo_vel,
        retarded_pos,
        retarded_vel,
        light_time_days: lt,
        light_time_rate: if apply {
            super::planet_range::RangeState::from_cartesian(geo, geo_vel).rate / C_LIGHT_AU_DAY
        } else {
            0.0
        },
        retarded_jd_tt: retarded_jd,
    })
}

/// Differentiate the converged one-way propagation equation.
///
/// For the retarded ray direction `n`, target velocity `r′` and observer
/// velocity `o′`, `lt′ = n·(r′−o′)/(c+n·r′)`. The returned geometric velocity
/// is `(1−lt′)r′−o′`, in AU per TT day. A disabled correction, nonpositive delay or
/// zero ray returns the instantaneous velocity difference. A nonfinite ray
/// norm also selects that difference.
///
/// This is an unchecked floating-point helper, not a `Result`-returning
/// validator. Supply finite positions in AU, velocities in AU/TT-day and a
/// delay in days. A zero denominator or other nonfinite arithmetic can produce
/// IEEE infinity/NaN; no typed error or convergence check is performed here.
pub fn light_time_corrected_velocity(
    retarded_vel: [f64; 3],
    observer_vel: [f64; 3],
    geo_pos: [f64; 3],
    light_time_days: f64,
    apply: bool,
) -> [f64; 3] {
    let plain = [
        retarded_vel[0] - observer_vel[0],
        retarded_vel[1] - observer_vel[1],
        retarded_vel[2] - observer_vel[2],
    ];
    if !apply || light_time_days <= 0.0 {
        return plain;
    }
    let dist = vec3_dist(geo_pos);
    if dist <= 0.0 || !dist.is_finite() {
        return plain;
    }
    let unit = [geo_pos[0] / dist, geo_pos[1] / dist, geo_pos[2] / dist];
    let radial_target =
        unit[0] * retarded_vel[0] + unit[1] * retarded_vel[1] + unit[2] * retarded_vel[2];
    let radial_rel = unit[0] * plain[0] + unit[1] * plain[1] + unit[2] * plain[2];
    let ltdot = radial_rel / (C_LIGHT_AU_DAY + radial_target);
    [
        retarded_vel[0] * (1.0 - ltdot) - observer_vel[0],
        retarded_vel[1] * (1.0 - ltdot) - observer_vel[1],
        retarded_vel[2] * (1.0 - ltdot) - observer_vel[2],
    ]
}

/// Evaluate a target at `jd_tt − light_time_days` through the same stored
/// provider and center-offset composition as its observation state.
///
/// The result is a geometric state in the channel's declared frame and units.
/// Coverage and missing-channel errors propagate without changing providers.
pub fn retarded_geometric_state<P: StateProvider>(
    provider: &P,
    body_id: i32,
    jd_tt: f64,
    light_time_days: f64,
    offset_at: impl Fn(f64) -> Option<[f64; 3]>,
) -> Result<GeometricState> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !light_time_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: light_time_days,
        });
    }
    provider.geometric_state_with_offset(body_id, jd_tt - light_time_days, offset_at)
}

/// Whether stellar aberration is enabled by the center and correction flags.
/// `FLG_NOABERR`, heliocentric, barycentric and geometric requests disable it;
/// `FLG_NOGDEFL` controls deflection independently.
pub fn wants_aberration(flags: i32) -> bool {
    flags & (FLG_NOABERR | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS) == 0
}

/// Whether the scalar-place model applies planetary light deflection.
///
/// Geometric, heliocentric, barycentric, no-deflection and nonpositive-delay
/// requests omit it. The scalar lunar model omits planetary deflection; the
/// surface-vector model can include that small term separately in
/// [`apply_deflection_and_aberration`].
pub fn wants_deflection(flags: i32, body_id: i32, light_time_days: f64) -> bool {
    if flags & (FLG_NOGDEFL | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS) != 0 {
        return false;
    }
    if body_id == MOON {
        return false;
    }
    light_time_days > 0.0
}

/// Apply stellar aberration to an observer-relative ICRS ray in AU.
///
/// Observer velocity is AU/day and light time is days. Positive light time
/// selects the Lorentz transformation, normalized by `light_time_days * c`;
/// nonpositive delay uses normalized first-order Bradley aberration. A zero
/// ray or zero observer velocity is unchanged. Consumed inputs must be finite;
/// non-finite transformations return a typed numerical error.
pub fn apply_aberration(
    geo: [f64; 3],
    earth_vel: [f64; 3],
    light_time_days: f64,
) -> Result<[f64; 3]> {
    apply_aberration_with_c(geo, earth_vel, light_time_days, C_LIGHT_AU_DAY)
}

/// Apply stellar aberration with an explicit light speed in AU/day.
///
/// Ray and velocity units are AU and AU/day; delay is days. The same Lorentz
/// or first-order transformation as [`apply_aberration`] is used, with the
/// supplied physical speed constant. Invalid numerical inputs return errors.
pub fn apply_aberration_with_c(
    geo: [f64; 3],
    earth_vel: [f64; 3],
    light_time_days: f64,
    c_au_day: f64,
) -> Result<[f64; 3]> {
    apply_aberration_with_correction(geo, earth_vel, light_time_days, c_au_day)
        .map(|value| value.position)
}

/// An apparent Cartesian position and its separately accumulated correction.
/// Both vectors use the input length unit. Keeping the small increment avoids
/// recovering a differential by subtracting two nearly equal astronomical rays.
pub(super) struct PlaceCorrection {
    pub(super) position: [f64; 3],
    pub(super) increment: [f64; 3],
}

/// Apply the same Lorentz map while retaining its small physical increment.
/// Position arithmetic, flags and guards match [`apply_aberration_with_c`].
/// The equivalent increment uses `sqrt(1-beta²)-1 = -beta²/(1+sqrt(1-beta²))`
/// to avoid cancellation; it supplies a well-conditioned correction derivative.
pub(super) fn apply_aberration_with_correction(
    geo: [f64; 3],
    earth_vel: [f64; 3],
    light_time_days: f64,
    c_au_day: f64,
) -> Result<PlaceCorrection> {
    check_finite_vec3(geo)?;
    check_finite_vec3(earth_vel)?;
    if !light_time_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: light_time_days,
        });
    }
    if !c_au_day.is_finite() {
        return Err(Error::NonFiniteArgument { value: c_au_day });
    }
    let dist = vec3_dist(geo);
    if dist == 0.0 {
        return Ok(PlaceCorrection {
            position: geo,
            increment: [0.0; 3],
        });
    }
    if light_time_days > 0.0 {
        let p1mag = light_time_days * c_au_day;
        let vemag = vec3_dist(earth_vel);
        if vemag == 0.0 || p1mag == 0.0 {
            return Ok(PlaceCorrection {
                position: geo,
                increment: [0.0; 3],
            });
        }
        let beta = vemag / c_au_day;
        let dot = geo[0] * earth_vel[0] + geo[1] * earth_vel[1] + geo[2] * earth_vel[2];
        let cosd = dot / (p1mag * vemag);
        let gammai = (1.0 - beta * beta).sqrt();
        let p = beta * cosd;
        let q = (1.0 + p / (1.0 + gammai)) * light_time_days;
        let mut r = 1.0 + p;
        let mut scalar_increment = -(beta * beta) / (1.0 + gammai) - p;
        if r.abs() < ABERRATION_R_FLOOR {
            r = 1.0;
            scalar_increment = -(beta * beta) / (1.0 + gammai);
        }
        let position = [
            (gammai * geo[0] + q * earth_vel[0]) / r,
            (gammai * geo[1] + q * earth_vel[1]) / r,
            (gammai * geo[2] + q * earth_vel[2]) / r,
        ];
        let increment = std::array::from_fn(|axis| {
            compensated_dot3(
                [scalar_increment, q, 0.0],
                [geo[axis], earth_vel[axis], 0.0],
            ) / r
        });
        Ok(PlaceCorrection {
            position,
            increment,
        })
    } else {
        let ux = geo[0] / dist;
        let uy = geo[1] / dist;
        let uz = geo[2] / dist;
        let vx = earth_vel[0] / c_au_day;
        let vy = earth_vel[1] / c_au_day;
        let vz = earth_vel[2] / c_au_day;
        let dot = ux * vx + uy * vy + uz * vz;
        let ax = ux + vx - ux * dot;
        let ay = uy + vy - uy * dot;
        let az = uz + vz - uz * dot;
        let a_dist = (ax * ax + ay * ay + az * az).sqrt();
        if a_dist == 0.0 {
            return Ok(PlaceCorrection {
                position: geo,
                increment: [0.0; 3],
            });
        }
        let unit = [ux, uy, uz];
        let delta = [vx - ux * dot, vy - uy * dot, vz - uz * dot];
        let scale_increment = -(2.0 * compensated_dot3(unit, delta)
            + compensated_dot3(delta, delta))
            / (1.0 + a_dist);
        let increment = std::array::from_fn(|axis| {
            dist / a_dist
                * compensated_dot3([1.0, scale_increment, 0.0], [delta[axis], unit[axis], 0.0])
        });
        Ok(PlaceCorrection {
            position: [ax / a_dist * dist, ay / a_dist * dist, az / a_dist * dist],
            increment,
        })
    }
}

/// Acquire one deflector's barycentric ICRS position in AU at a TT Julian Day.
///
/// Missing channels or uncovered LEB/SPK epochs omit an optional deflector;
/// corrupt data and other errors propagate. The strict-lookup variant of the
/// shared PPN kernel lets callers require all configured deflectors.
pub fn deflector_position(
    lookup: &impl Fn(i32, f64) -> Result<[f64; 3]>,
    body_id: i32,
    jd_tt: f64,
) -> Result<Option<[f64; 3]>> {
    match lookup(body_id, jd_tt) {
        Ok(pos) => {
            check_finite_vec3(pos)?;
            Ok(Some(pos))
        }
        Err(
            Error::LebBodyNotFound { .. }
            | Error::LebOutOfRange { .. }
            | Error::SpkOutOfRange { .. }
            | Error::HorizonsApiError { .. }
            | Error::HorizonsMalformed { .. }
            | Error::UnknownBody { .. },
        ) => Ok(None),
        Err(other) => Err(other),
    }
}

/// Translate an acquired SSB deflector to a fixed reception origin, then
/// rotate it into the ray/observer axes. Units remain AU. The caller supplies
/// its actual rotation; this helper acquires no states. Rotation errors propagate.
pub(crate) fn deflector_in_frame(
    position: [f64; 3],
    origin: [f64; 3],
    rotate: impl FnOnce([f64; 3]) -> Result<[f64; 3]>,
) -> Result<[f64; 3]> {
    rotate(std::array::from_fn(|axis| position[axis] - origin[axis]))
}

/// Apply finite-source PPN light deflection by the Sun, Jupiter and Saturn.
///
/// The observer-relative ray and barycentric observer/deflector positions
/// are normally ICRS AU. A common orthonormal rotation and fixed origin shift
/// are also valid, provided the ray and every observer/deflector position use
/// that same frame and origin at both query epochs. The epoch is TT JD and
/// the one-way delay is days. Each
/// deflector is sampled at closest approach, clamped between emission and
/// reception. The gamma-one term uses [`GRAV_GS`], the deflector mass ratio
/// and the SOFA near-conjunction denominator limit. Interior or degenerate
/// geometries are omitted. A zero ray is unchanged. Optional missing support
/// follows [`deflector_position`]; invalid states and other errors propagate.
pub fn apply_gravitational_deflection(
    geo: [f64; 3],
    observer_pos: [f64; 3],
    jd_tt: f64,
    light_time_days: f64,
    deflector_pos: impl Fn(i32, f64) -> Result<[f64; 3]>,
) -> Result<[f64; 3]> {
    apply_gravitational_deflection_with_lookup(
        geo,
        observer_pos,
        jd_tt,
        light_time_days,
        |body, jd| deflector_position(&deflector_pos, body, jd),
    )
}

/// The same optional-deflector PPN map with its small accumulated increment.
/// Positions, AU units, closest-approach clocks and missing/range policy match
/// [`apply_gravitational_deflection`]; no large corrected rays are subtracted.
pub(super) fn apply_gravitational_deflection_with_correction(
    geo: [f64; 3],
    observer_pos: [f64; 3],
    jd_tt: f64,
    light_time_days: f64,
    deflector_pos: impl Fn(i32, f64) -> Result<[f64; 3]>,
) -> Result<PlaceCorrection> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    gravitational_deflection_correction_at_epoch(
        geo,
        observer_pos,
        &jd_tt,
        light_time_days,
        |epoch, offset| Ok(*epoch - offset),
        |body, epoch| deflector_position(&deflector_pos, body, *epoch),
    )
}

/// Apply the same PPN deflection formula with a caller-owned acquisition
/// policy. Arguments and output are ICRS AU vectors; observation TT and
/// one-way light time are days. The lookup receives each deflector body
/// and its observation or closest-approach TT. `None` omits a deflector;
/// errors propagate unchanged. This permits a public vector reader's
/// strict alias/coverage errors without changing the raw LEB reader's
/// optional-deflector policy in [`apply_gravitational_deflection`].
///
/// Numerical guards, deflector order and closest-approach clipping are
/// identical to [`apply_gravitational_deflection`]. Non-finite inputs or
/// acquired vectors return [`Error::NonFiniteArgument`]; a zero target
/// vector returns unchanged before acquiring any deflector.
pub fn apply_gravitational_deflection_with_lookup(
    geo: [f64; 3],
    observer_pos: [f64; 3],
    jd_tt: f64,
    light_time_days: f64,
    lookup: impl Fn(i32, f64) -> Result<Option<[f64; 3]>>,
) -> Result<[f64; 3]> {
    apply_gravitational_deflection_with_lookup_correction(
        geo,
        observer_pos,
        jd_tt,
        light_time_days,
        lookup,
    )
    .map(|value| value.position)
}

/// PPN position/increment with the caller's actual required/optional lookups.
pub(super) fn apply_gravitational_deflection_with_lookup_correction(
    geo: [f64; 3],
    observer_pos: [f64; 3],
    jd_tt: f64,
    light_time_days: f64,
    lookup: impl Fn(i32, f64) -> Result<Option<[f64; 3]>>,
) -> Result<PlaceCorrection> {
    check_finite_vec3(geo)?;
    check_finite_vec3(observer_pos)?;
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    gravitational_deflection_correction_at_epoch(
        geo,
        observer_pos,
        &jd_tt,
        light_time_days,
        |epoch, offset| Ok(*epoch - offset),
        |body, epoch| lookup(body, *epoch),
    )
}

/// PPN deflection over an actual acquired observation epoch. Numerical rays,
/// light time and returned vectors retain the existing finite checked domain.
/// The epoch remains opaque: original deflector queries receive it directly,
/// and closest-approach subtraction occurs only after its geometry/clipping.
/// The ordinary scalar entry above keeps its independent finite-JD admission.
pub fn apply_gravitational_deflection_with_epoch_lookup<E>(
    geo: [f64; 3],
    observer_pos: [f64; 3],
    jd_tt: &E,
    light_time_days: f64,
    subtract: impl Fn(&E, f64) -> Result<E>,
    lookup: impl Fn(i32, &E) -> Result<Option<[f64; 3]>>,
) -> Result<[f64; 3]> {
    gravitational_deflection_correction_at_epoch(
        geo,
        observer_pos,
        jd_tt,
        light_time_days,
        subtract,
        lookup,
    )
    .map(|value| value.position)
}

/// Shared PPN arithmetic with a separately accumulated finite-source increment.
fn gravitational_deflection_correction_at_epoch<E>(
    geo: [f64; 3],
    observer_pos: [f64; 3],
    jd_tt: &E,
    light_time_days: f64,
    subtract: impl Fn(&E, f64) -> Result<E>,
    lookup: impl Fn(i32, &E) -> Result<Option<[f64; 3]>>,
) -> Result<PlaceCorrection> {
    check_finite_vec3(geo)?;
    check_finite_vec3(observer_pos)?;
    if !light_time_days.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: light_time_days,
        });
    }
    let pmag = vec3_dist(geo);
    if pmag == 0.0 {
        return Ok(PlaceCorrection {
            position: geo,
            increment: [0.0; 3],
        });
    }
    let mut result = geo;
    let mut increment = [0.0; 3];
    for (defl_body_id, rmass) in DEFLECTORS {
        let Some(defl_pos) = lookup(defl_body_id, jd_tt)? else {
            continue;
        };
        check_finite_vec3(defl_pos)?;
        let gpv = [
            defl_pos[0] - observer_pos[0],
            defl_pos[1] - observer_pos[1],
            defl_pos[2] - observer_pos[2],
        ];
        let phat = [result[0] / pmag, result[1] / pmag, result[2] / pmag];
        let dlt = (phat[0] * gpv[0] + phat[1] * gpv[1] + phat[2] * gpv[2]) / C_LIGHT_AU_DAY;
        let tclose_offset = dlt.max(0.0).min(light_time_days);
        let tclose_jd = subtract(jd_tt, tclose_offset)?;
        let Some(defl_close) = lookup(defl_body_id, &tclose_jd)? else {
            continue;
        };
        check_finite_vec3(defl_close)?;
        let pe = [
            observer_pos[0] - defl_close[0],
            observer_pos[1] - defl_close[1],
            observer_pos[2] - defl_close[2],
        ];
        let pq = [result[0] + pe[0], result[1] + pe[1], result[2] + pe[2]];
        let qmag = vec3_dist(pq);
        let emag = vec3_dist(pe);
        if qmag == 0.0 || emag < DEFLECTION_EMAG_FLOOR_AU {
            continue;
        }
        let qhat = [pq[0] / qmag, pq[1] / qmag, pq[2] / qmag];
        let ehat = [pe[0] / emag, pe[1] / emag, pe[2] / emag];
        let pdotq = phat[0] * qhat[0] + phat[1] * qhat[1] + phat[2] * qhat[2];
        let qdote = qhat[0] * ehat[0] + qhat[1] * ehat[1] + qhat[2] * ehat[2];
        let edotp = ehat[0] * phat[0] + ehat[1] * phat[1] + ehat[2] * phat[2];
        if edotp.abs() > DEFLECTION_EDOTP_LIMIT {
            continue;
        }
        let fac1 =
            2.0 * GRAV_GS / (SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S * emag * AU_IN_M * rmass);
        let dlim = DEFLECTION_DLIM_SCALE / (if emag * emag > 1.0 { emag * emag } else { 1.0 });
        let mut fac2 = 1.0 + qdote;
        if fac2 < dlim {
            fac2 = dlim;
        }
        let coeff = fac1 / fac2 * pmag;
        for axis in 0..3 {
            let correction = coeff * (pdotq * ehat[axis] - edotp * qhat[axis]);
            result[axis] += correction;
            increment[axis] += correction;
        }
    }
    Ok(PlaceCorrection {
        position: result,
        increment,
    })
}

/// Apply Earth's finite-source PPN deflection for a surface observer.
///
/// `target_au` is an observer-relative ICRS ray in AU and `observer_gcrs_au`
/// is the geocentric site vector in AU. The gamma-one finite-deflection term
/// uses Earth's mass and is suppressed for a nadir-to-apparent-limb-radius
/// ratio below 0.8, where the far-field model is unsuitable. Zero rays and a
/// geocentric observer are unchanged. Non-finite inputs are errors.
pub fn apply_earth_deflection(target_au: [f64; 3], observer_gcrs_au: [f64; 3]) -> Result<[f64; 3]> {
    apply_earth_deflection_with_correction(target_au, observer_gcrs_au).map(|value| value.position)
}

/// Surface-Earth PPN map retaining its actual small limb-gated increment.
fn apply_earth_deflection_with_correction(
    target_au: [f64; 3],
    observer_gcrs_au: [f64; 3],
) -> Result<PlaceCorrection> {
    check_finite_vec3(target_au)?;
    check_finite_vec3(observer_gcrs_au)?;
    let pmag = vec3_dist(target_au);
    if pmag == 0.0 {
        return Ok(PlaceCorrection {
            position: target_au,
            increment: [0.0; 3],
        });
    }
    // `pe`: observer relative to the deflector. The Earth *is* the GCRS
    // origin, so the geocentric observer vector is used directly.
    let pe = observer_gcrs_au;
    let emag = vec3_dist(pe);
    if emag == 0.0 {
        return Ok(PlaceCorrection {
            position: target_au,
            increment: [0.0; 3],
        });
    }
    let pq = [
        target_au[0] + pe[0],
        target_au[1] + pe[1],
        target_au[2] + pe[2],
    ];
    let qmag = vec3_dist(pq);
    // Use a unit denominator for a degenerate source-to-deflector vector.
    let q_denom = if qmag != 0.0 { qmag } else { 1.0 };
    let phat = [
        target_au[0] / pmag,
        target_au[1] / pmag,
        target_au[2] / pmag,
    ];
    let qhat = [pq[0] / q_denom, pq[1] / q_denom, pq[2] / q_denom];
    let ehat = [pe[0] / emag, pe[1] / emag, pe[2] / emag];
    let pdotq = phat[0] * qhat[0] + phat[1] * qhat[1] + phat[2] * qhat[2];
    let qdote = qhat[0] * ehat[0] + qhat[1] * ehat[1] + qhat[2] * ehat[2];
    let edotp = ehat[0] * phat[0] + ehat[1] * phat[1] + ehat[2] * phat[2];
    // Omit the far-field term where the ray's nadir angle is below 0.8 of the
    // apparent Earth-limb radius.
    let earth_radius_au = IERS_EARTH_RADIUS_M / AU_IN_M;
    let aprad = (earth_radius_au / emag).min(1.0).asin();
    let coszd =
        (target_au[0] * pe[0] + target_au[1] * pe[1] + target_au[2] * pe[2]) / (pmag * emag);
    let zdobj = coszd.clamp(-1.0, 1.0).acos();
    let nadir_angle = (std::f64::consts::PI - zdobj) / aprad;
    if nadir_angle < EARTH_DEFLECTION_LIMB_GATE {
        return Ok(PlaceCorrection {
            position: target_au,
            increment: [0.0; 3],
        });
    }
    let flag = if edotp.abs() <= DEFLECTION_EDOTP_LIMIT {
        1.0
    } else {
        0.0
    };
    let fac1 = 2.0 * GRAV_GS
        / (SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S * emag * AU_IN_M * EARTH_DEFLECTION_RMASS);
    let fac2 = 1.0 + qdote;
    // Scale the finite-source angular deflection by the ray length.
    let d = [
        ((flag * fac1 * (pdotq * ehat[0] - edotp * qhat[0])) / fac2) * pmag,
        ((flag * fac1 * (pdotq * ehat[1] - edotp * qhat[1])) / fac2) * pmag,
        ((flag * fac1 * (pdotq * ehat[2] - edotp * qhat[2])) / fac2) * pmag,
    ];
    Ok(PlaceCorrection {
        position: [
            target_au[0] + d[0],
            target_au[1] + d[1],
            target_au[2] + d[2],
        ],
        increment: d,
    })
}

impl AberrationDeflection {
    /// Deflection- and aberration-corrected position (ICRS, AU).
    pub fn pos(self) -> [f64; 3] {
        self.pos
    }

    /// Deflection increment (`deflected - geo`, AU): zero when deflection
    /// was gated off. The speed stencil freezes this term for Sun targets instead of
    /// re-evaluating the degenerate PPN geometry on extrapolated vectors.
    pub fn deflection_delta(self) -> [f64; 3] {
        self.deflection_delta
    }

    /// Algebraically retained total PPN/Earth/Lorentz correction in ICRS AU.
    /// It is accumulated before adding corrections to the large ray, so its
    /// derivative is well-conditioned on short actual provider intervals.
    pub(super) fn correction_increment(self) -> [f64; 3] {
        self.correction_increment
    }
}

/// Apply planetary deflection, optional surface-Earth deflection, then aberration.
///
/// The relative ray and barycentric observer position are ICRS AU; observer
/// velocity is AU/day, epoch is TT JD and delay is days. `deflector_pos`
/// acquires barycentric positions at reception and closest approach.
/// `earth_site_au` optionally enables Earth's surface-observer term.
/// Geometric requests suppress all apparent corrections. Scalar lunar
/// requests omit planetary deflection; surface-vector requests include it.
///
/// The result retains the corrected position and algebraic correction
/// increment for stable derivatives. Optional support follows
/// [`deflector_position`]; strict-vector acquisition, non-finite and corrupt
/// state errors propagate unchanged.
#[allow(clippy::too_many_arguments)]
pub fn apply_deflection_and_aberration(
    geo: [f64; 3],
    aberration_vel: [f64; 3],
    observer_pos: [f64; 3],
    jd_tt: f64,
    light_time_days: f64,
    body_id: i32,
    flags: i32,
    deflector_pos: impl Fn(i32, f64) -> Result<[f64; 3]>,
    earth_site_au: Option<[f64; 3]>,
) -> Result<AberrationDeflection> {
    apply_deflection_and_aberration_with_vector_route(
        geo,
        aberration_vel,
        observer_pos,
        jd_tt,
        light_time_days,
        body_id,
        flags,
        deflector_pos,
        earth_site_au,
        earth_site_au.is_some(),
    )
}

/// Apply apparent-place corrections with an independent vector-route marker.
/// The input vectors are ICRS AU/AU-day, TT is JD, and light time is days.
/// `vector_route` includes lunar planetary deflection even for a geocenter;
/// `earth_site_au` independently adds the surface observer's Earth-mass term.
/// Geometric requests retain neither correction, and scalar Moon requests
/// retain their established planetary-deflection skip. Correction order,
/// finite checks and returned increments follow the public correction helper.
/// Vector-route deflector lookup failures propagate; scalar routes preserve
/// their optional raw lookup policy. Selected correction paths acquire their
/// required deflectors without a separate state-prefetch stage.
#[allow(clippy::too_many_arguments)]
pub(super) fn apply_deflection_and_aberration_with_vector_route(
    geo: [f64; 3],
    aberration_vel: [f64; 3],
    observer_pos: [f64; 3],
    jd_tt: f64,
    light_time_days: f64,
    body_id: i32,
    flags: i32,
    deflector_pos: impl Fn(i32, f64) -> Result<[f64; 3]>,
    earth_site_au: Option<[f64; 3]>,
    vector_route: bool,
) -> Result<AberrationDeflection> {
    let mut pos = geo;
    let mut deflection_delta = [0.0, 0.0, 0.0];
    let mut correction_increment = [0.0; 3];
    // Surface-vector lunar places include the planetary deflection term when
    // apparent corrections and a positive light time are enabled.
    let moon_vector_defl = vector_route
        && body_id == MOON
        && (flags & (FLG_NOGDEFL | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS)) == 0
        && light_time_days > 0.0;
    if wants_deflection(flags, body_id, light_time_days) || moon_vector_defl {
        let pre = pos;
        let corrected = if vector_route {
            apply_gravitational_deflection_with_lookup_correction(
                pos,
                observer_pos,
                jd_tt,
                light_time_days,
                |body, jd| deflector_pos(body, jd).map(Some),
            )?
        } else {
            apply_gravitational_deflection_with_correction(
                pos,
                observer_pos,
                jd_tt,
                light_time_days,
                deflector_pos,
            )?
        };
        pos = corrected.position;
        correction_increment = corrected.increment;
        deflection_delta = [pos[0] - pre[0], pos[1] - pre[1], pos[2] - pre[2]];
    }
    if let Some(site) = earth_site_au {
        // Geometric requests omit the Earth surface-deflection term.
        if (flags & FLG_TRUEPOS) == 0 {
            let pre = pos;
            let corrected = apply_earth_deflection_with_correction(pos, site)?;
            pos = corrected.position;
            for (total, correction) in correction_increment.iter_mut().zip(corrected.increment) {
                *total += correction;
            }
            deflection_delta = [
                deflection_delta[0] + (pos[0] - pre[0]),
                deflection_delta[1] + (pos[1] - pre[1]),
                deflection_delta[2] + (pos[2] - pre[2]),
            ];
        }
    }
    if wants_aberration(flags) {
        let corrected =
            apply_aberration_with_correction(pos, aberration_vel, light_time_days, C_LIGHT_AU_DAY)?;
        pos = corrected.position;
        for (total, correction) in correction_increment.iter_mut().zip(corrected.increment) {
            *total += correction;
        }
    }
    Ok(AberrationDeflection {
        pos,
        deflection_delta,
        correction_increment,
    })
}

#[cfg(test)]
mod acquired_deflector_epoch_tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn deflector_queries_retain_original_and_actual_closest_approach_epochs() {
        struct Epoch(u8);
        let original = Epoch(0);
        let calls = RefCell::new(Vec::new());
        let result = apply_gravitational_deflection_with_epoch_lookup(
            [1.0, 2.0, 3.0],
            [0.0; 3],
            &original,
            0.02,
            |epoch, _| {
                assert!(std::ptr::eq(epoch, &original));
                calls.borrow_mut().push((-1, epoch.0));
                Ok(Epoch(1))
            },
            |body, epoch| {
                calls.borrow_mut().push((body, epoch.0));
                Ok(Some([2.0, 3.0, 4.0]))
            },
        )
        .unwrap();
        assert!(result.into_iter().all(f64::is_finite));
        assert_eq!(
            *calls.borrow(),
            [
                (SUN, 0),
                (-1, 0),
                (SUN, 1),
                (JUPITER, 0),
                (-1, 0),
                (JUPITER, 1),
                (SATURN, 0),
                (-1, 0),
                (SATURN, 1)
            ]
        );
    }
}

#[cfg(test)]
#[path = "correction_increment_tests.rs"]
mod correction_increment_tests;

#[cfg(test)]
mod deflection_frame_tests {
    use super::*;
    use std::cell::{Cell, RefCell};

    #[test]
    fn fixed_origin_adapter_is_covariant_under_cyclic_axes() {
        let tt = 2451545.0;
        let origin = [8.0, -4.0, 2.0];
        let ray = [0.0, 4.0, 0.0];
        let observer = [1.0, 0.0, 0.0];
        let run = |rotate: fn([f64; 3]) -> [f64; 3]| {
            let queries = RefCell::new(Vec::new());
            let converted = RefCell::new(Vec::new());
            let output = apply_gravitational_deflection(
                rotate(ray),
                rotate(observer),
                tt,
                0.1,
                |body, epoch| {
                    queries.borrow_mut().push((body, epoch));
                    if body != SUN {
                        return Err(Error::UnknownBody { body_id: body });
                    }
                    // A moving deflector: the second acquisition must keep the
                    // reception origin, not move that origin with the deflector.
                    let y = 0.5 + 8.0 * (epoch - tt);
                    let source = [origin[0], origin[1] + y, origin[2]];
                    let value = deflector_in_frame(source, origin, |v| Ok(rotate(v)))?;
                    assert_eq!(value, rotate([0.0, y, 0.0]));
                    converted.borrow_mut().push(value);
                    Ok(value)
                },
            )
            .unwrap();
            assert_ne!(output, rotate(ray));
            let queries = queries.into_inner();
            assert_eq!(queries.len(), 4);
            assert_eq!(
                queries.iter().map(|v| v.0).collect::<Vec<_>>(),
                [SUN, SUN, JUPITER, SATURN]
            );
            assert_eq!(queries[0].1, tt);
            assert!(queries[1].1 < tt && queries[1].1 >= tt - 0.1);
            let converted = converted.into_inner();
            assert_ne!(converted[0], converted[1]);
            (output, queries)
        };
        let plain = run(|v| v);
        let rotated = run(|v| [v[1], v[2], v[0]]);
        assert_eq!(rotated.0, [plain.0[1], plain.0[2], plain.0[0]]);
        assert_eq!(rotated.1, plain.1);
    }

    #[test]
    fn adapter_keeps_zero_ray_and_optional_source_errors_lazy() {
        assert_eq!(
            apply_gravitational_deflection([0.0; 3], [1.0; 3], 2451545.0, 0.1, |_, _| panic!(
                "zero ray must not acquire or transform a deflector"
            ))
            .unwrap(),
            [0.0; 3]
        );
        let errors = [
            Error::UnknownBody { body_id: SUN },
            Error::SpkOutOfRange {
                target: 10,
                et: 2.0,
                start: 0.0,
                end: 1.0,
            },
            Error::SpkBodyNotFound { target: 10 },
            Error::SpkCorrupt {
                what: "synthetic invalid record".into(),
                offset: 8,
            },
        ];
        for (index, error) in errors.into_iter().enumerate() {
            let reads = RefCell::new(Vec::new());
            let rotations = Cell::new(0);
            let result = apply_gravitational_deflection(
                [1.0, 2.0, 3.0],
                [0.0; 3],
                2451545.0,
                0.1,
                |body, _| {
                    reads.borrow_mut().push(body);
                    let acquired: Result<[f64; 3]> = Err(error.clone());
                    acquired.and_then(|position| {
                        deflector_in_frame(position, [0.0; 3], |v| {
                            rotations.set(rotations.get() + 1);
                            Ok(v)
                        })
                    })
                },
            );
            assert_eq!(rotations.get(), 0);
            if index < 2 {
                assert_eq!(result.unwrap(), [1.0, 2.0, 3.0]);
                assert_eq!(*reads.borrow(), [SUN, JUPITER, SATURN]);
            } else {
                assert_eq!(result, Err(error));
                assert_eq!(*reads.borrow(), [SUN]);
            }
        }
    }
}
