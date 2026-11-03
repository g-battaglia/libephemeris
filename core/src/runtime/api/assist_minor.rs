//! Borrowed ASSIST minor-body fallback and its scalar-TT projection.
//!
//! Public suppliers return one checked heliocentric ecliptic-J2000 state in
//! AU and AU/day; their projection retains that held state for speed samples.
//! The concrete native owner instead evaluates its numerical trajectory at the
//! finite emission epochs and speed samples, retaining the simulation and its
//! force ephemeris together. Geocentric observers use the retained main SPK
//! image. Neither route substitutes analytical elements for an unavailable
//! numerical backend.

use super::Session;
use crate::constants::{
    EARTH, FLG_EQUATORIAL, FLG_HELCTR, FLG_J2000, FLG_NOABERR, FLG_NONUT, FLG_SIDEREAL, FLG_SPEED,
    FLG_TRUEPOS, J2000, SUN,
};
use crate::errors::{Error, Result};
use crate::nbody::{CartesianState, PropagationResult};
use crate::positions::{C_LIGHT_AU_DAY, LocalSpkStateProvider};
use crate::types::TwoPartEpoch;

/// Native minor-body numerical provider consulted after registered-SPK coverage declines.
/// Availability false preserves the configured analytical or strict-precision policy.
pub trait MinorAssistAcquisition: Send + Sync {
    /// Retain a newly provisioned minor-body SPK after a registered-source miss.
    /// The caller merges only these resources into its immutable configuration
    /// snapshot. This hook runs only on the native branch after a LEB miss;
    /// the default performs no download or registry mutation.
    fn refreshed_spk(
        &self,
        _body: i32,
        _tt: f64,
        _context: &crate::context::EphemerisContext,
    ) -> Result<Option<crate::spk::SpkRegistry>> {
        Ok(None)
    }

    /// Whether the body has an orbital-elements record for numerical integration.
    fn contains(&self, body: i32) -> Result<bool>;
    /// Whether the numerical backend and required ephemeris data are available.
    fn available(&self) -> Result<bool>;
    /// Acquire the curated body's checked heliocentric state at scalar TT.
    fn state_at(&self, body: i32, tt: f64) -> Result<PropagationResult>;

    /// Acquire the same central PV plus an optional heliocentric acceleration
    /// in fixed ecliptic J2000 AU/day². Apparent SPEED differentiates both the
    /// position and velocity parameters of its light-time model. None explicitly
    /// retains a supplier's existing constant-velocity model; it is not an error
    /// retry. For integrated trajectories the native N-body supplier returns
    /// its actual retained force acceleration. A selected stored trajectory
    /// retains its original PV without inventing that force channel.
    /// Provider/state/force errors propagate without substituting another model.
    fn state_with_acceleration_at(
        &self,
        body: i32,
        tt: f64,
    ) -> Result<(PropagationResult, Option<[f64; 3]>)> {
        self.state_at(body, tt).map(|state| (state, None))
    }
}

/// Geometric observer states for one TT sample, in SSB ICRS AU and AU/day.
pub trait AssistMinorEphemeris {
    /// Geometric Earth position at the TT Julian day.
    fn earth_position(&self, tt: f64) -> Result<[f64; 3]>;
    /// Geometric Sun position at the TT Julian day.
    fn sun_position(&self, tt: f64) -> Result<[f64; 3]>;
    /// Barycentric Earth velocity for the aberration calculation, AU/day.
    fn earth_velocity(&self, tt: f64) -> Result<[f64; 3]>;
}

struct Main<'a> {
    session: &'a Session,
}

fn scalar_time(tt: f64) -> TwoPartEpoch {
    let whole = tt.floor();
    TwoPartEpoch {
        hi: whole,
        lo: tt - whole,
    }
}

impl AssistMinorEphemeris for Main<'_> {
    fn earth_position(&self, tt: f64) -> Result<[f64; 3]> {
        #[cfg(test)]
        let observation = crate::positions::context_serving::hypothetical_support_observations::ObservationGuard::enter(None, Vec::new());
        let result =
            LocalSpkStateProvider::new(self.session.ctx.spk_registry(), self.session.ctx.mode())?
                .nodaps_state_two_part(EARTH, scalar_time(tt));
        #[cfg(test)]
        crate::nbody::owned_input_observations::note(|| {
            serde_json::json!({
                "kind":"main_read", "role":"earth_position", "tt":tt,
                "spk_reads":observation.finish(),
                "result":match &result { Ok(s) => serde_json::json!({"Ok":crate::nbody::owned_input_observations::pv(CartesianState{position:s.pos(),velocity:s.vel()})}), Err(e) => serde_json::json!({"Err":format!("{e:?}")}) }
            })
        });
        result.map(|s| s.pos())
    }
    fn sun_position(&self, tt: f64) -> Result<[f64; 3]> {
        #[cfg(test)]
        let observation = crate::positions::context_serving::hypothetical_support_observations::ObservationGuard::enter(None, Vec::new());
        let result =
            LocalSpkStateProvider::new(self.session.ctx.spk_registry(), self.session.ctx.mode())?
                .nodaps_state_two_part(SUN, scalar_time(tt));
        #[cfg(test)]
        crate::nbody::owned_input_observations::note(|| {
            serde_json::json!({
                "kind":"main_read", "role":"sun_position", "tt":tt,
                "spk_reads":observation.finish(),
                "result":match &result { Ok(s) => serde_json::json!({"Ok":crate::nbody::owned_input_observations::pv(CartesianState{position:s.pos(),velocity:s.vel()})}), Err(e) => serde_json::json!({"Err":format!("{e:?}")}) }
            })
        });
        result.map(|s| s.pos())
    }
    fn earth_velocity(&self, tt: f64) -> Result<[f64; 3]> {
        #[cfg(test)]
        let observation = crate::positions::context_serving::hypothetical_support_observations::ObservationGuard::enter(None, Vec::new());
        let result =
            LocalSpkStateProvider::new(self.session.ctx.spk_registry(), self.session.ctx.mode())?
                .nodaps_state_two_part(EARTH, scalar_time(tt));
        #[cfg(test)]
        crate::nbody::owned_input_observations::note(|| {
            serde_json::json!({
                "kind":"main_read", "role":"earth_velocity", "tt":tt,
                "spk_reads":observation.finish(),
                "result":match &result { Ok(s) => serde_json::json!({"Ok":crate::nbody::owned_input_observations::pv(CartesianState{position:s.pos(),velocity:s.vel()})}), Err(e) => serde_json::json!({"Err":format!("{e:?}")}) }
            })
        });
        result.map(|s| s.vel())
    }
}

impl Session {
    /// Project one acquired ASSIST state at scalar TT JD using this Session's
    /// already acquired main image. Returns true/mean ecliptic-of-date longitude
    /// and latitude in degrees and distance AU. Two constant-velocity light-time
    /// iterations precede Bradley aberration, Vondrák precession and conditional
    /// longitude nutation. HELCTR consumes no main provider; TRUEPOS suppresses
    /// light time and aberration. SIDEREAL+EQUATORIAL suppresses this nutation.
    ///
    /// This method does not apply the final public output-frame/sidereal tail or
    /// topocentric parallax. The
    /// numerical state, SPK, time and frame errors propagate; no file
    /// acquisition or fallback is performed here.
    pub fn assist_minor_position_from_state(
        &self,
        state: &PropagationResult,
        tt: f64,
        flags: i32,
    ) -> Result<[f64; 3]> {
        let main = Main { session: self };
        assist_minor_sample(state, tt, flags, &main)
    }

    /// Acquire one central numerical state and its local derivative samples.
    /// TT is JD; six output components are degrees/degrees/AU and daily rates. SPEED
    /// uses a fourth-order centered position derivative over local kinematic
    /// states, with a dyadic TT spacing. Native apparent samples include the
    /// acquired force acceleration in the changing velocity parameter of light
    /// time; suppliers without acceleration retain their declared linear model.
    /// Geocentric range uses the same
    /// already acquired samples, avoiding one-second norm cancellation and
    /// inexact full-JD spacing. Heliocentric radial velocity is the analytic
    /// derivative of the Cartesian norm, including both light-time iterations
    /// and the velocity parameter's derivative when acceleration is supplied.
    /// Only the central state invokes the integrator. Longitude increments take the short arc;
    /// absent SPEED retains three positive zeros.
    /// Availability, exception fallback and final output frames belong to the
    /// enclosing public branch and are deliberately not repeated by this method.
    pub fn assist_minor_ecliptic(
        &self,
        supplier: &dyn MinorAssistAcquisition,
        body: i32,
        tt: f64,
        flags: i32,
    ) -> Result<[f64; 6]> {
        self.assist_minor_in_frame(supplier, body, tt, flags, false)
    }

    /// ASSIST position and rates in the final equatorial output frame.
    ///
    /// Arguments and units match [`Self::assist_minor_ecliptic`], with
    /// `FLG_EQUATORIAL` selected. One central integration is retained; four
    /// shifted samples pass through the complete time-dependent position
    /// transformation before the fourth-order angular derivative is formed.
    /// This includes the obliquity derivative, J2000 precession and the
    /// mean-obliquity choice under NONUT/SIDEREAL. Radial rates use the same
    /// rotation-invariant model as the ecliptic entry. The caller must not
    /// apply a second equatorial/precession conversion to this result.
    /// Underlying provider, epoch and coordinate errors propagate unchanged.
    pub fn assist_minor_equatorial(
        &self,
        supplier: &dyn MinorAssistAcquisition,
        body: i32,
        tt: f64,
        flags: i32,
    ) -> Result<[f64; 6]> {
        self.assist_minor_in_frame(supplier, body, tt, flags | FLG_EQUATORIAL, true)
    }

    /// Final nonsidereal J2000-ecliptic state and rates from one held PV.
    ///
    /// Apply the complete requested frame to each position sample before its
    /// derivative. Differentiating an of-date rate and then using a separate
    /// half-day frame stencil loses the mixed motion/frame terms for fast arcs.
    /// The public of-date helper retains its existing contract; the caller must
    /// not precess this completed J2000 state a second time.
    pub(crate) fn assist_minor_inertial_ecliptic(
        &self,
        supplier: &dyn MinorAssistAcquisition,
        body: i32,
        tt: f64,
        flags: i32,
    ) -> Result<[f64; 6]> {
        debug_assert!(flags & FLG_J2000 != 0 && flags & (FLG_EQUATORIAL | FLG_SIDEREAL) == 0);
        self.assist_minor_in_frame(supplier, body, tt, flags, true)
    }

    /// Shared held-state calculation; the selected frame is applied to every
    /// position sample, so angular rates differentiate that actual frame.
    fn assist_minor_in_frame(
        &self,
        supplier: &dyn MinorAssistAcquisition,
        body: i32,
        tt: f64,
        flags: i32,
        apply_output_frame: bool,
    ) -> Result<[f64; 6]> {
        let (central, acceleration) = if flags & FLG_SPEED != 0 && flags & FLG_TRUEPOS == 0 {
            supplier.state_with_acceleration_at(body, tt)?
        } else {
            (supplier.state_at(body, tt)?, None)
        };
        if let Some(acceleration) = acceleration {
            for value in acceleration {
                if !value.is_finite() {
                    return Err(Error::NonFiniteArgument { value });
                }
            }
        }
        let project = |state: &PropagationResult, epoch: f64| -> Result<[f64; 3]> {
            let position = self.assist_minor_position_from_state(state, epoch, flags)?;
            if !apply_output_frame {
                return Ok(position);
            }
            let value = crate::positions::horizons_output::local_convert(
                [position[0], position[1], position[2], 0.0, 0.0, 0.0],
                epoch,
                flags & !FLG_SPEED,
                false,
            )?;
            Ok([value[0], value[1], value[2]])
        };
        let [lon, lat, distance] = project(&central, tt)?;
        let mut rates = [0.0; 3];
        if flags & FLG_SPEED != 0 {
            let sample = |offset| {
                let epoch = tt + offset;
                project(
                    &shifted_kinematic_state(central, acceleration, epoch)?,
                    epoch,
                )
            };
            let h = POSITION_RATE_STEP_DAYS;
            rates = fourth_order_position_rates(
                sample(-h)?,
                sample(h)?,
                sample(-2.0 * h)?,
                sample(2.0 * h)?,
            );
            if flags & FLG_HELCTR != 0 {
                rates[2] = heliocentric_radial_velocity_with_acceleration(
                    &central,
                    acceleration,
                    flags & FLG_TRUEPOS == 0,
                )?;
            }
        }
        Ok([lon, lat, distance, rates[0], rates[1], rates[2]])
    }

    /// Actual owned numerical trajectory at the finite emission epochs.
    /// The opaque public supplier entry above deliberately retains its held PV.
    pub(crate) fn assist_minor_owned_position(
        &self,
        supplier: &super::NativeMinorAssist<'_>,
        body: i32,
        tt: f64,
        flags: i32,
    ) -> std::result::Result<[f64; 6], crate::nbody::OwnedTrajectoryError> {
        if flags & FLG_TRUEPOS != 0 {
            let result = if flags & FLG_EQUATORIAL != 0 {
                self.assist_minor_equatorial(supplier, body, tt, flags)
            } else {
                self.assist_minor_ecliptic(supplier, body, tt, flags)
            };
            return result.map_err(crate::nbody::OwnedTrajectoryError::Initial);
        }
        supplier.with_native_trajectory(body, tt, |central, target, sun| {
            // Original state validation precedes every observer acquisition.
            crate::positions::cartesian_to_spherical_deg([central.x, central.y, central.z])?;
            let main = Main { session: self };
            let mut project = |offset: f64, want_rate: bool| -> Result<([f64; 3], f64)> {
                let epoch = tt + offset;
                #[cfg(test)]
                crate::nbody::owned_input_observations::note(|| {
                    serde_json::json!({
                        "kind":"project_start","offset":offset,"reception_tt":epoch,
                        "want_rate":want_rate,"flags":flags
                    })
                });
                if offset != 0.0 && (epoch == tt || epoch - tt != offset) {
                    return Err(Error::MathDomain { value: offset });
                }
                let central = target(offset)?;
                let mut earth_velocity = None;
                let observer = if flags & FLG_HELCTR != 0 {
                    // Same force-file Sun and TT-as-TDB clock as the orbit.
                    sun(offset)?
                } else {
                    let position = main.earth_position(epoch)?;
                    // Preserve this original acquisition/error contract. The
                    // barycentric emitter subtracts Earth SSB directly: mixing
                    // a force-Sun clock with this main Sun would shift origins.
                    main.sun_position(epoch)?;
                    if flags & FLG_NOABERR == 0 {
                        earth_velocity = Some(main.earth_velocity(epoch)?);
                    }
                    let velocity = if want_rate {
                        super::NativeMinorAssist::earth_position_tt_rate(
                            &self.ctx, epoch, position,
                        )?
                    } else {
                        [0.0; 3]
                    };
                    CartesianState { position, velocity }
                };
                #[cfg(test)]
                crate::nbody::owned_input_observations::note(|| {
                    serde_json::json!({
                        "kind":"observer","reception_tt":epoch,
                        "ssb_icrf":crate::nbody::owned_input_observations::pv(observer),
                        "physical_earth_velocity":earth_velocity
                    })
                });
                let mut emission = |delay_offset| target(offset + delay_offset);
                let ray = finite_emission_ray(central, observer, &mut emission, want_rate)?;
                let radius = magnitude(ray.position);
                let radial_rate = if want_rate {
                    ray.position
                        .iter()
                        .zip(ray.velocity)
                        .map(|(p, v)| p / radius * v)
                        .sum()
                } else {
                    0.0
                };
                let apparent = if let Some(velocity) = earth_velocity {
                    bradley_aberration(ray.position, velocity)
                } else {
                    ray.position
                };
                #[cfg(test)]
                crate::nbody::owned_input_observations::note(|| {
                    serde_json::json!({
                        "kind":"apparent_ray","reception_tt":epoch,
                        "geometric_icrf":crate::nbody::owned_input_observations::pv(ray),
                        "apparent_icrf":apparent,"physical_earth_velocity":earth_velocity
                    })
                });
                let orbital = crate::spk::icrs_to_ecliptic_j2000(apparent);
                let mut position = orbital_ray_place(orbital, epoch, flags)?;
                if flags & FLG_EQUATORIAL != 0 {
                    let value = crate::positions::horizons_output::local_convert(
                        [position[0], position[1], position[2], 0.0, 0.0, 0.0],
                        epoch,
                        flags & !FLG_SPEED,
                        false,
                    )?;
                    position = [value[0], value[1], value[2]];
                }
                #[cfg(test)]
                crate::nbody::owned_input_observations::note(|| {
                    serde_json::json!({
                        "kind":"project_finish","reception_tt":epoch,
                        "position":position,"radial_rate":radial_rate
                    })
                });
                Ok((position, radial_rate))
            };
            let (position, radial) = project(0.0, flags & FLG_SPEED != 0)?;
            let mut rates = [0.0; 3];
            if flags & FLG_SPEED != 0 {
                let h = POSITION_RATE_STEP_DAYS;
                let angular = fourth_order_position_rates(
                    project(-h, false)?.0,
                    project(h, false)?.0,
                    project(-2.0 * h, false)?.0,
                    project(2.0 * h, false)?.0,
                );
                rates[..2].copy_from_slice(&angular[..2]);
                rates[2] = radial;
            }
            Ok([
                position[0],
                position[1],
                position[2],
                rates[0],
                rates[1],
                rates[2],
            ])
        })
    }
}

/// Local Taylor state whose first derivatives are the acquired v and a.
/// The public constant-velocity shift remains unchanged for linear suppliers.
fn shifted_kinematic_state(
    state: PropagationResult,
    acceleration: Option<[f64; 3]>,
    epoch: f64,
) -> Result<PropagationResult> {
    let Some(acceleration) = acceleration else {
        return state.shifted(epoch);
    };
    let dt = epoch - state.jd_tt;
    let position = [state.x, state.y, state.z];
    let velocity = [state.vx, state.vy, state.vz];
    PropagationResult::from_state(
        crate::nbody::CartesianState {
            position: std::array::from_fn(|axis| {
                position[axis] + velocity[axis] * dt + 0.5 * acceleration[axis] * dt * dt
            }),
            velocity: std::array::from_fn(|axis| velocity[axis] + acceleration[axis] * dt),
        },
        epoch,
    )
}

/// TT-day spacing shared by the fourth-order local position derivatives.
/// The dyadic offset is exact throughout the supported Julian-day interval.
const POSITION_RATE_STEP_DAYS: f64 = 1.0 / 64.0;

/// Differentiate the same local position model without amplifying one-second
/// rounding or full-JD sample-spacing error. Units are deg/day, deg/day, AU/day.
fn fourth_order_position_rates(
    before: [f64; 3],
    after: [f64; 3],
    before2: [f64; 3],
    after2: [f64; 3],
) -> [f64; 3] {
    let longitude_delta = |next: f64, previous: f64| {
        let delta = next - previous;
        if delta > 180.0 {
            delta - 360.0
        } else if delta < -180.0 {
            delta + 360.0
        } else {
            delta
        }
    };
    [
        (8.0 * longitude_delta(after[0], before[0]) - longitude_delta(after2[0], before2[0]))
            / (12.0 * POSITION_RATE_STEP_DAYS),
        (8.0 * (after[1] - before[1]) - (after2[1] - before2[1]))
            / (12.0 * POSITION_RATE_STEP_DAYS),
        (8.0 * (after[2] - before[2]) - (after2[2] - before2[2]))
            / (12.0 * POSITION_RATE_STEP_DAYS),
    ]
}

/// Project an owned numerical state through a geometric observer.
/// See [`Session::assist_minor_position_from_state`] for units, corrections and limits.
pub fn assist_minor_sample(
    state: &PropagationResult,
    tt: f64,
    flags: i32,
    ephemeris: &dyn AssistMinorEphemeris,
) -> Result<[f64; 3]> {
    let initial = [state.x, state.y, state.z];
    // Preserve state validation before any observer acquisition. All geometry
    // below stays in the orbital frame until the complete ray is available.
    crate::positions::cartesian_to_spherical_deg(initial)?;
    let mut ray = initial;
    if flags & FLG_HELCTR != 0 {
        if flags & FLG_TRUEPOS == 0 {
            let velocity = [state.vx, state.vy, state.vz];
            let mut light_time = 0.0;
            for _ in 0..2 {
                let delayed = std::array::from_fn(|i| initial[i] - velocity[i] * light_time);
                light_time = magnitude(delayed) / C_LIGHT_AU_DAY;
            }
            ray = std::array::from_fn(|i| initial[i] - velocity[i] * light_time);
        }
    } else {
        let velocity = [state.vx, state.vy, state.vz];
        let earth = ephemeris.earth_position(tt)?;
        let sun = ephemeris.sun_position(tt)?;
        let observer =
            crate::spk::icrs_to_ecliptic_j2000(std::array::from_fn(|i| earth[i] - sun[i]));
        ray = std::array::from_fn(|i| initial[i] - observer[i]);
        if flags & FLG_TRUEPOS == 0 {
            for _ in 0..2 {
                let light_time = magnitude(ray) / C_LIGHT_AU_DAY;
                ray =
                    std::array::from_fn(|i| (initial[i] - velocity[i] * light_time) - observer[i]);
            }
        }
        if flags & (FLG_TRUEPOS | FLG_NOABERR) == 0 {
            let earth_velocity = crate::spk::icrs_to_ecliptic_j2000(ephemeris.earth_velocity(tt)?);
            ray = bradley_aberration(ray, earth_velocity);
        }
    }
    orbital_ray_place(ray, tt, flags)
}

/// Shared numerical-source frame boundary; the geometry has already finished.
fn orbital_ray_place(ray: [f64; 3], tt: f64, flags: i32) -> Result<[f64; 3]> {
    let [longitude, latitude, distance] =
        crate::positions::cartesian_to_spherical_deg(orbital_ray_in_engine_j2000(ray)?)?;
    let (mut longitude, latitude) =
        crate::hypothetical::precess_ecliptic_source_angles(longitude, latitude, J2000, tt)?;
    let sidereal_equatorial = flags & FLG_SIDEREAL != 0 && flags & FLG_EQUATORIAL != 0;
    if flags & (FLG_NONUT | FLG_J2000) == 0 && !sidereal_equatorial {
        longitude = (longitude + crate::erfa::nut06a(tt)?.0.to_degrees()).rem_euclid(360.0);
    }
    Ok([longitude, latitude, distance])
}

/// Re-express the complete orbital-frame ray in the output engine's J2000 basis.
///
/// PropagationResult uses the fixed orbital obliquity, while the frame tail
/// uses IAU 2006 J2000 obliquity. Their connecting rotation is mandatory for
/// every output frame. Applying it after geometry preserves the exact light-
/// time and aberration acquisitions; rotating only the target would mix bases.
fn orbital_ray_in_engine_j2000([x, y, z]: [f64; 3]) -> Result<[f64; 3]> {
    let orbital_obliquity = crate::nbody::ECLIPTIC_J2000_OBLIQUITY_DEG;
    let engine_obliquity = crate::positions::mean_obliquity_iau2006_deg(J2000)?;
    let (sin_delta, cos_delta) = (orbital_obliquity - engine_obliquity)
        .to_radians()
        .sin_cos();
    Ok([
        x,
        y * cos_delta - z * sin_delta,
        y * sin_delta + z * cos_delta,
    ])
}

/// Two finite light-time updates on a barycentric trajectory, in AU/TT day.
/// Source offsets are relative to reception, avoiding a rounded full-JD delay.
/// The observer stays at reception; its velocity is its actual TT derivative.
fn finite_emission_ray(
    central: CartesianState,
    observer: CartesianState,
    emission: &mut dyn FnMut(f64) -> Result<CartesianState>,
    want_rate: bool,
) -> Result<CartesianState> {
    let mut delay = 0.0;
    let mut delay_rate = 0.0;
    let mut ray = CartesianState {
        position: [0.0; 3],
        velocity: [0.0; 3],
    };
    for step in 0..=2 {
        let target = if step == 0 {
            central
        } else {
            emission(-delay)?
        };
        ray.position = std::array::from_fn(|i| target.position[i] - observer.position[i]);
        for value in ray.position {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        if want_rate {
            ray.velocity = std::array::from_fn(|i| {
                target.velocity[i] * (1.0 - delay_rate) - observer.velocity[i]
            });
            for value in ray.velocity {
                if !value.is_finite() {
                    return Err(Error::NonFiniteArgument { value });
                }
            }
        }
        #[cfg(test)]
        crate::nbody::owned_input_observations::note(|| {
            serde_json::json!({
                "kind":"finite_iterate","step":step,"delay":delay,"delay_rate":delay_rate,
                "target":crate::nbody::owned_input_observations::pv(target),
                "observer":crate::nbody::owned_input_observations::pv(observer),
                "ray":crate::nbody::owned_input_observations::pv(ray),"want_rate":want_rate
            })
        });
        let radius = magnitude(ray.position);
        if want_rate && radius == 0.0 {
            return Err(Error::MathDomain { value: radius });
        }
        if step < 2 {
            delay = radius / C_LIGHT_AU_DAY;
            if want_rate {
                delay_rate = ray
                    .position
                    .iter()
                    .zip(ray.velocity)
                    .map(|(p, v)| p / radius * v)
                    .sum::<f64>()
                    / C_LIGHT_AU_DAY;
            }
        }
    }
    Ok(ray)
}

fn magnitude([x, y, z]: [f64; 3]) -> f64 {
    (x * x + y * y + z * z).sqrt()
}

/// Derivative of the heliocentric range, AU/day, under the local linear state.
///
/// For a light-time iterate `L`, the ray is `q = p - v L` and its derivative
/// is `q' = v (1 - L')`. Thus the next delay satisfies
/// `L_next' = (q · q') / (c |q|)`. The same two delay updates used by the
/// position calculation are differentiated here. Precession and nutation are
/// rotations and do not change range. A moving state at the origin has no
/// defined norm derivative; the stationary origin has derivative zero.
#[cfg(test)]
fn heliocentric_radial_velocity(state: &PropagationResult, apparent: bool) -> Result<f64> {
    heliocentric_radial_velocity_with_acceleration(state, None, apparent)
}

/// Differentiate q=p-vL as q'=v-aL-vL', including the changing target velocity.
/// Zero acceleration retains the declared linear model and its old arithmetic.
fn heliocentric_radial_velocity_with_acceleration(
    state: &PropagationResult,
    acceleration: Option<[f64; 3]>,
    apparent: bool,
) -> Result<f64> {
    let position = [state.x, state.y, state.z];
    let velocity = [state.vx, state.vy, state.vz];
    let mut delay = 0.0;
    let mut delay_rate = 0.0;
    let iterations = if apparent { 2 } else { 0 };
    for step in 0..=iterations {
        let ray: [f64; 3] = std::array::from_fn(|i| position[i] - velocity[i] * delay);
        let radius = magnitude(ray);
        if radius == 0.0 {
            return if velocity == [0.0; 3] {
                Ok(0.0)
            } else {
                Err(Error::ZeroVectorToSpherical)
            };
        }
        let numerator = (ray[0] * velocity[0] + ray[1] * velocity[1] + ray[2] * velocity[2])
            * (1.0 - delay_rate);
        let numerator = if let Some(acceleration) = acceleration {
            numerator
                - (ray[0] * acceleration[0] + ray[1] * acceleration[1] + ray[2] * acceleration[2])
                    * delay
        } else {
            numerator
        };
        let radial = numerator / radius;
        if step == iterations {
            return Ok(radial);
        }
        delay = radius / C_LIGHT_AU_DAY;
        delay_rate = radial / C_LIGHT_AU_DAY;
    }
    unreachable!("the last delay iteration returns its derivative")
}

/// First-order Bradley aberration, with velocities in AU/day.
fn bradley_aberration(position: [f64; 3], velocity: [f64; 3]) -> [f64; 3] {
    let distance = magnitude(position);
    if distance < 1e-30 {
        return position;
    }
    let direction = position.map(|v| v / distance);
    let norm = magnitude(direction);
    let displacement = if norm < 1e-30 {
        [0.0; 3]
    } else {
        let unit = direction.map(|v| v / norm);
        let beta = velocity.map(|v| v / C_LIGHT_AU_DAY);
        let dot = unit[0] * beta[0] + unit[1] * beta[1] + unit[2] * beta[2];
        std::array::from_fn(|i| beta[i] - dot * unit[i])
    };
    let adjusted: [f64; 3] = std::array::from_fn(|i| direction[i] + displacement[i]);
    let norm = magnitude(adjusted);
    if norm < 1e-30 {
        return position;
    }
    adjusted.map(|v| v / norm * distance)
}

#[cfg(test)]
#[path = "assist_minor_data_tests.rs"]
mod data_tests;

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nbody::CartesianState;
    use std::cell::RefCell;
    use std::sync::atomic::{AtomicUsize, Ordering};
    // This optional capture records Unix inode/time fingerprints. Portable
    // integration and numerical contract tests remain enabled on every OS.
    #[cfg(unix)]
    #[test]
    #[ignore = "three actual owned ASSIST requests; explicit immutable inputs/output required"]
    fn capture_runtime_n16_owned_inputs() {
        use serde_json::{Value, json};
        use sha2::{Digest, Sha256};
        use std::collections::BTreeMap;
        use std::io::{Read, Write};
        use std::path::{Path, PathBuf};

        fn fingerprint(path: &Path) -> Value {
            use std::os::unix::fs::MetadataExt;
            let stat = || {
                let m = std::fs::metadata(path).unwrap();
                json!([
                    m.dev(),
                    m.ino(),
                    m.len(),
                    m.mtime(),
                    m.mtime_nsec(),
                    m.ctime(),
                    m.ctime_nsec()
                ])
            };
            let before = stat();
            let mut file = std::fs::File::open(path).unwrap();
            let mut digest = Sha256::new();
            let mut buffer = [0_u8; 65536];
            loop {
                let count = file.read(&mut buffer).unwrap();
                if count == 0 {
                    break;
                }
                digest.update(&buffer[..count]);
            }
            assert_eq!(before, stat());
            json!({"sha256":digest.finalize().iter().map(|byte|format!("{byte:02x}")).collect::<String>(),"stat":before})
        }
        fn word(value: &Value) -> String {
            format!("{:016x}", value.as_f64().unwrap().to_bits())
        }
        let current_input = std::env::var_os("LIBEPHEMERIS_ERIS_CURRENT_INPUT")
            .map(|input| (input, 0))
            .or_else(|| {
                std::env::var_os("LIBEPHEMERIS_CURATED_CURRENT_INPUT").map(|input| (input, 1))
            })
            .or_else(|| {
                std::env::var_os("LIBEPHEMERIS_BENNU_CURRENT_INPUT").map(|input| (input, 2))
            });
        if let Some((current_input, scope_kind)) = current_input {
            let curated = scope_kind == 1;
            let bennu = scope_kind == 2;
            // A separate current four-request scope. The historical three-call
            // input and every original word guard below remain unchanged.
            let input = PathBuf::from(current_input);
            let identity = fingerprint(&input);
            assert_eq!(
                identity["sha256"].as_str().unwrap(),
                std::env::var(match scope_kind {
                    1 => "LIBEPHEMERIS_CURATED_CURRENT_SHA256",
                    2 => "LIBEPHEMERIS_BENNU_CURRENT_SHA256",
                    _ => "LIBEPHEMERIS_ERIS_CURRENT_SHA256",
                })
                .unwrap()
            );
            let packet: Value =
                serde_json::from_reader(std::fs::File::open(&input).unwrap()).unwrap();
            let selected = packet["rows"].as_array().unwrap();
            let expected_calls = match scope_kind {
                1 => 25,
                2 => 2,
                _ => 4,
            };
            assert_eq!(selected.len(), expected_calls);
            let config = &selected[0]["request"]["config"];
            assert!(
                selected
                    .iter()
                    .all(|row| row["request"]["config"] == *config)
            );
            assert_eq!(
                std::env::var_os("ASSIST_DIR").map(PathBuf::from),
                Some(PathBuf::from(config["assist_dir"].as_str().unwrap()))
            );
            assert_eq!(
                std::env::var_os("HOME").map(PathBuf::from),
                Some(PathBuf::from(config["assist_home"].as_str().unwrap()))
            );
            let scope: std::collections::BTreeSet<_> = selected
                .iter()
                .map(|row| {
                    let request = &row["request"];
                    assert_eq!(request["kwargs"], json!({}));
                    assert_eq!(request["args"][0], 2451545.0);
                    let body = request["args"][1].as_i64().unwrap();
                    if curated {
                        assert!(matches!(
                            body,
                            18 | 19
                                | 146472
                                | 100482
                                | 60000
                                | 109942
                                | 10704
                                | 10052
                                | 10016
                                | 11221
                                | 11566
                                | 11685
                                | 10055
                                | 35143
                                | 172173
                        ));
                        assert_eq!(request["function"], "calc");
                        assert!(matches!(request["args"][2].as_i64().unwrap(), 258 | 266));
                    } else if bennu {
                        assert_eq!(body, 111955);
                        assert_eq!(request["function"], "calc");
                        assert!(matches!(request["args"][2].as_i64().unwrap(), 258 | 266));
                    } else {
                        assert_eq!(body, 146199);
                    }
                    (
                        request["function"].as_str().unwrap(),
                        body,
                        request["args"][2].as_i64().unwrap(),
                    )
                })
                .collect();
            assert_eq!(scope.len(), expected_calls);
            if !curated && !bennu {
                assert_eq!(
                    scope,
                    [
                        ("calc", 146199, 258),
                        ("calc", 146199, 266),
                        ("calc_ut", 146199, 258),
                        ("calc_ut", 146199, 266)
                    ]
                    .into_iter()
                    .collect()
                );
            }
            let mut paths = BTreeMap::new();
            for source in packet["inputs"].as_object().unwrap().values() {
                let path = PathBuf::from(source["path"].as_str().unwrap());
                let actual = fingerprint(&path);
                assert_eq!(actual["sha256"], source["sha256"]);
                assert_eq!(actual["stat"][2], source["bytes"]);
                paths.insert(path, actual);
            }
            let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
            let mut source_paths: Vec<_> = [
                "src/runtime/api/assist_minor.rs",
                "src/runtime/api/nbody_minor.rs",
                "src/bodies/nbody.rs",
                "src/bodies/nbody/source.rs",
                "src/positions/main_owner.rs",
                "src/providers/spk/reader.rs",
                "src/providers/spk/positional.rs",
                "src/providers/spk/evaluation.rs",
                "src/providers/spk/mod.rs",
            ]
            .map(|path| root.join(path))
            .into_iter()
            .collect();
            source_paths.push(std::env::current_exe().unwrap());
            let sources: BTreeMap<_, _> = source_paths
                .iter()
                .map(|path| (path.clone(), fingerprint(path)))
                .collect();
            let runtime = crate::config::Runtime::default();
            runtime.set_network_policy(Some("sealed")).unwrap();
            runtime.set_auto_spk_download(Some(false));
            runtime.set_iers_auto_download(Some(false));
            runtime
                .set_precision_tier(config["tier"].as_str().unwrap())
                .unwrap();
            runtime.set_calc_mode(config["mode"].as_str()).unwrap();
            runtime.set_iers_delta_t_enabled(config["iers_delta_t_enabled"].as_bool());
            runtime.set_delta_t_userdef(config["delta_t_userdef"].as_f64());
            runtime.set_strict_precision(config["strict_precision"].as_bool());
            runtime.set_ephe_path(Some(PathBuf::from(config["ephe_path"].as_str().unwrap())));
            runtime.set_ephemeris_file(config["ephemeris_file"].as_str().unwrap());
            let mut rows = Vec::new();
            for saved in selected {
                let request = &saved["request"];
                let args = request["args"].as_array().unwrap();
                let observation = crate::nbody::owned_input_observations::Scope::enter();
                let spk_scope = crate::positions::context_serving::hypothetical_support_observations::ObservationGuard::enter(None, Vec::new());
                let trace = crate::trace::start_session();
                let call = if request["function"] == "calc" {
                    runtime.calc(
                        args[0].as_f64().unwrap(),
                        args[1].as_i64().unwrap() as i32,
                        args[2].as_i64().unwrap() as i32,
                    )
                } else {
                    runtime.calc_ut(
                        args[0].as_f64().unwrap(),
                        args[1].as_i64().unwrap() as i32,
                        args[2].as_i64().unwrap() as i32,
                    )
                };
                let providers = crate::trace::get_trace_results();
                let trace_reset = trace.reset();
                let events = observation.finish();
                let auxiliary_spk_reads = spk_scope.finish();
                let mut row = saved.clone();
                row["events"] = json!(events);
                row["auxiliary_spk_reads"] = json!(auxiliary_spk_reads);
                row["providers"] = json!(providers);
                row["trace_reset"] = json!(trace_reset);
                match call {
                    Ok(value) => {
                        let actual = value.values.map(|v| format!("{:016x}", v.to_bits()));
                        let expected: Vec<_> = saved["native_result"][0]
                            .as_array()
                            .unwrap()
                            .iter()
                            .map(word)
                            .collect();
                        row["result"] = json!([value.values, value.returned_flags]);
                        row["actual_words"] = json!(actual);
                        row["expected_words"] = json!(expected);
                        row["exact_public_result"] = json!(
                            actual.as_slice() == expected.as_slice()
                                && value.returned_flags as i64
                                    == saved["native_result"][1].as_i64().unwrap()
                        );
                    }
                    Err(error) => {
                        row["error"] = json!(format!("{error:?}"));
                        row["exact_public_result"] = json!(false);
                    }
                }
                rows.push(row);
            }
            let after: BTreeMap<_, _> = paths
                .keys()
                .map(|path| (path.clone(), fingerprint(path)))
                .collect();
            let sources_after: BTreeMap<_, _> = source_paths
                .iter()
                .map(|path| (path.clone(), fingerprint(path)))
                .collect();
            let complete = paths == after
                && sources == sources_after
                && fingerprint(&input) == identity
                && rows
                    .iter()
                    .all(|row| row["exact_public_result"] == true && row["trace_reset"] == true);
            let result = json!({"complete":complete,"actual_calls":expected_calls,"rows":rows,
                "input":input,"input_identity":identity,"inputs_before":paths,"inputs_after":after,
                "sources_before":sources,"sources_after":sources_after,
                "ordinary_native_binary_sha256":packet["ordinary_native_binary_sha256"],
                "clock_probe_binary_sha256":packet["clock_probe_binary_sha256"],
                "public_qualified":false});
            let output = std::env::var_os("LIBEPHEMERIS_N16_OWNED_CAPTURE").unwrap();
            let mut file = std::fs::OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(output)
                .unwrap();
            serde_json::to_writer_pretty(&mut file, &result).unwrap();
            file.write_all(b"\n").unwrap();
            assert!(
                complete,
                "four current Eris requests failed original full-word/input guards"
            );
            return;
        }
        let input = PathBuf::from(
            std::env::var_os("LIBEPHEMERIS_N16_OWNED_SERIAL").expect("serial capture"),
        );
        let output =
            PathBuf::from(std::env::var_os("LIBEPHEMERIS_N16_OWNED_CAPTURE").expect("new output"));
        let serial_identity = fingerprint(&input);
        assert_eq!(
            serial_identity["sha256"],
            "dae2e1f86a0c6113d28e488ea977234f74fcc05b7e310c12e8f1872cba81e9b2"
        );
        let serial: Value = serde_json::from_reader(std::fs::File::open(&input).unwrap()).unwrap();
        assert_eq!(serial["complete"], true);
        assert_eq!(serial["calls"], 7);
        assert_eq!(serial["identities_before"], serial["identities_after"]);
        assert_eq!(
            serial["extension_after"]["sha256"],
            "38770e88e60eab058c890a3a18fe9a40aab5996b221635eee6ad197ac5ba232a"
        );
        let selected: Vec<_> = serial["rows"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|row| {
                let args = row["args"].as_array().unwrap();
                matches!(
                    (args[1].as_i64().unwrap(), args[2].as_i64().unwrap()),
                    (19, 266) | (10016, 258) | (10016, 266)
                )
            })
            .collect();
        assert_eq!(selected.len(), 3);
        let config = &selected[0]["config"];
        assert!(selected.iter().all(|r| r["config"] == *config
            && r["function"] == "calc"
            && r["kwargs"] == json!({})));
        assert_eq!(
            std::env::var_os("ASSIST_DIR").map(PathBuf::from),
            Some(PathBuf::from(config["assist_dir"].as_str().unwrap()))
        );
        let mut paths = BTreeMap::new();
        for source in selected[0]["inputs"].as_object().unwrap().values() {
            let path = PathBuf::from(source["path"].as_str().unwrap());
            let identity = fingerprint(&path);
            assert_eq!(identity["sha256"], source["sha256"]);
            paths.insert(path, identity);
        }
        let runtime = crate::config::Runtime::default();
        runtime.set_network_policy(Some("sealed")).unwrap();
        runtime.set_auto_spk_download(Some(false));
        runtime.set_iers_auto_download(Some(false));
        runtime
            .set_precision_tier(config["tier"].as_str().unwrap())
            .unwrap();
        runtime.set_calc_mode(config["mode"].as_str()).unwrap();
        runtime.set_iers_delta_t_enabled(config["iers_delta_t_enabled"].as_bool());
        runtime.set_delta_t_userdef(config["delta_t_userdef"].as_f64());
        runtime.set_strict_precision(config["strict_precision"].as_bool());
        runtime.set_ephe_path(Some(PathBuf::from(config["ephe_path"].as_str().unwrap())));
        runtime.set_ephemeris_file(config["ephemeris_file"].as_str().unwrap());
        let mut rows = Vec::new();
        for saved in selected {
            let args = saved["args"].as_array().unwrap();
            let scope = crate::nbody::owned_input_observations::Scope::enter();
            let spk_scope = crate::positions::context_serving::hypothetical_support_observations::ObservationGuard::enter(None, Vec::new());
            let trace = crate::trace::start_session();
            let result = runtime.calc(
                args[0].as_f64().unwrap(),
                args[1].as_i64().unwrap() as i32,
                args[2].as_i64().unwrap() as i32,
            );
            let providers = crate::trace::get_trace_results();
            let trace_reset = trace.reset();
            let events = scope.finish();
            let auxiliary_spk_reads = spk_scope.finish();
            let mut row = json!({"id":saved["id"],"function":"calc","args":args,"events":events,
                             "auxiliary_spk_reads":auxiliary_spk_reads,"providers":providers,"trace_reset":trace_reset});
            match result {
                Ok(value) => {
                    let actual = value.values.map(|v| format!("{:016x}", v.to_bits()));
                    let expected: Vec<_> = saved["result"][0]
                        .as_array()
                        .unwrap()
                        .iter()
                        .map(word)
                        .collect();
                    row["result"] = json!([value.values, value.returned_flags]);
                    row["actual_words"] = json!(actual);
                    row["expected_words"] = json!(expected);
                    row["exact_public_result"] = json!(
                        actual.as_slice() == expected.as_slice()
                            && value.returned_flags as i64 == saved["result"][1].as_i64().unwrap()
                    );
                }
                Err(error) => {
                    row["error"] = json!(format!("{error:?}"));
                    row["exact_public_result"] = json!(false);
                }
            }
            rows.push(row);
        }
        let after: BTreeMap<_, _> = paths
            .keys()
            .map(|path| (path.clone(), fingerprint(path)))
            .collect();
        let complete = paths == after
            && fingerprint(&input) == serial_identity
            && rows
                .iter()
                .all(|r| r["exact_public_result"] == true && r["trace_reset"] == true);
        let result = json!({"complete":complete,"actual_calls":3,"rows":rows,
                        "serial_capture":input,"serial_identity":serial_identity,
                        "inputs_before":paths,"inputs_after":after,
                        "public_qualified":false,"rigorous_error_bound":null});
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(output)
            .unwrap();
        serde_json::to_writer_pretty(&mut file, &result).unwrap();
        file.write_all(b"\n").unwrap();
        assert!(
            complete,
            "preserved output does not match original 387 serial capture"
        );
    }

    fn state() -> PropagationResult {
        PropagationResult::from_state(
            CartesianState {
                position: [2.0, 0.0, 0.0],
                velocity: [0.0, 0.01, 0.0],
            },
            J2000,
        )
        .unwrap()
    }
    #[test]
    fn complete_j2000_frame_preserves_fast_arc_kinematic_derivative() {
        struct Supplier(PropagationResult);
        impl MinorAssistAcquisition for Supplier {
            fn contains(&self, _: i32) -> Result<bool> {
                Ok(true)
            }
            fn available(&self) -> Result<bool> {
                Ok(true)
            }
            fn state_at(&self, _: i32, _: f64) -> Result<PropagationResult> {
                Ok(self.0)
            }
        }
        // Raw heliocentric PV from the observed fast-arc failure. Its fixed
        // orbital ecliptic differs from the output's IAU 2006 J2000 ecliptic.
        // Apply the constant basis change to both P and V before evaluating
        // the exact Cartesian-to-spherical kinematic identity. No neighbouring
        // epoch or fitted correction participates in this reference.
        let epoch = 2460311.5;
        let position = [-1.0517507657155696, 0.40478519300293, 0.03189628032597022];
        let velocity = [
            -0.0021599040146951477,
            -0.01722746538710307,
            -0.00011703341743518007,
        ];
        let supplier = Supplier(
            PropagationResult::from_state(CartesianState { position, velocity }, epoch).unwrap(),
        );
        let session = Session::from_context(crate::context::EphemerisContext::new());
        let flags = FLG_HELCTR | FLG_TRUEPOS | FLG_J2000 | FLG_NONUT | FLG_SPEED;
        let output = session
            .assist_minor_inertial_ecliptic(&supplier, 35143, epoch, flags)
            .unwrap();
        let (sin_orbital, cos_orbital) = crate::nbody::ECLIPTIC_J2000_OBLIQUITY_DEG
            .to_radians()
            .sin_cos();
        let (sin_output, cos_output) = crate::positions::mean_obliquity_iau2006_deg(J2000)
            .unwrap()
            .to_radians()
            .sin_cos();
        // Independently compose orbital ecliptic -> ICRF -> output ecliptic,
        // rather than reusing the production complete-ray bridge.
        let in_output_basis = |[x, y, z]: [f64; 3]| {
            let equatorial_y = y * cos_orbital - z * sin_orbital;
            let equatorial_z = y * sin_orbital + z * cos_orbital;
            [
                x,
                equatorial_y * cos_output + equatorial_z * sin_output,
                -equatorial_y * sin_output + equatorial_z * cos_output,
            ]
        };
        let [x, y, z] = in_output_basis(position);
        let [vx, vy, vz] = in_output_basis(velocity);
        let rho = x.hypot(y);
        let latitude_rate =
            ((rho * vz - z * (x * vx + y * vy) / rho) / (rho * rho + z * z)).to_degrees();
        let bound = latitude_rate.abs() * 1e-6;
        assert!(
            (output[4] - latitude_rate).abs() <= bound,
            "output latitude rate={} deg/day, invariant={latitude_rate} deg/day, bound={bound}",
            output[4]
        );
        let old = session
            .assist_minor_ecliptic(&supplier, 35143, epoch, flags)
            .unwrap();
        let old =
            crate::positions::horizons_output::local_convert(old, epoch, flags, false).unwrap();
        assert!((old[4] - latitude_rate).abs() > bound);
        assert_eq!(output[..3], old[..3]);
        let no_speed = session
            .assist_minor_inertial_ecliptic(&supplier, 35143, epoch, flags & !FLG_SPEED)
            .unwrap();
        assert_eq!(no_speed[..3], output[..3]);
        assert!(no_speed[3..].iter().all(|rate| rate.to_bits() == 0));
    }

    #[test]
    fn apparent_kinematics_are_acquired_once_and_not_requested_for_geometric_positions() {
        struct Supplier {
            pv: AtomicUsize,
            kinematics: AtomicUsize,
        }
        impl MinorAssistAcquisition for Supplier {
            fn contains(&self, _: i32) -> Result<bool> {
                Ok(true)
            }
            fn available(&self) -> Result<bool> {
                Ok(true)
            }
            fn state_at(&self, _: i32, _: f64) -> Result<PropagationResult> {
                self.pv.fetch_add(1, Ordering::Relaxed);
                Ok(state())
            }
            fn state_with_acceleration_at(
                &self,
                _: i32,
                _: f64,
            ) -> Result<(PropagationResult, Option<[f64; 3]>)> {
                self.kinematics.fetch_add(1, Ordering::Relaxed);
                Ok((state(), Some([0.0; 3])))
            }
        }
        let supplier = Supplier {
            pv: AtomicUsize::new(0),
            kinematics: AtomicUsize::new(0),
        };
        let session = Session::from_context(crate::context::EphemerisContext::new());
        for flags in [
            FLG_HELCTR | FLG_SPEED,
            FLG_HELCTR,
            FLG_HELCTR | FLG_TRUEPOS | FLG_SPEED,
        ] {
            session
                .assist_minor_ecliptic(&supplier, 17, J2000, flags)
                .unwrap();
        }
        assert_eq!(supplier.kinematics.load(Ordering::Relaxed), 1);
        assert_eq!(supplier.pv.load(Ordering::Relaxed), 2);
    }

    #[test]
    fn invalid_supplied_acceleration_is_an_error_not_a_linear_model_retry() {
        struct Supplier;
        impl MinorAssistAcquisition for Supplier {
            fn contains(&self, _: i32) -> Result<bool> {
                Ok(true)
            }
            fn available(&self) -> Result<bool> {
                Ok(true)
            }
            fn state_at(&self, _: i32, _: f64) -> Result<PropagationResult> {
                Ok(state())
            }
            fn state_with_acceleration_at(
                &self,
                _: i32,
                _: f64,
            ) -> Result<(PropagationResult, Option<[f64; 3]>)> {
                Ok((state(), Some([f64::NAN, 0.0, 0.0])))
            }
        }
        let session = Session::from_context(crate::context::EphemerisContext::new());
        assert!(matches!(
            session.assist_minor_ecliptic(&Supplier, 17, J2000, FLG_HELCTR | FLG_SPEED),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn circular_apparent_range_has_exactly_zero_derivative_with_its_acceleration() {
        let epoch = J2000;
        let state = PropagationResult::from_state(
            CartesianState {
                position: [1.0, 0.0, 0.0],
                velocity: [0.0, 1.0, 0.0],
            },
            epoch,
        )
        .unwrap();
        // p(t)=(cos t,sin t,0), v(t)=(-sin t,cos t,0). Both
        // light-time updates have constant norm, so the apparent range is
        // constant too. Holding v fixed instead introduces a false range rate.
        let acceleration = Some([-1.0, 0.0, 0.0]);
        assert_eq!(
            heliocentric_radial_velocity_with_acceleration(&state, acceleration, true).unwrap(),
            0.0
        );
        assert_ne!(heliocentric_radial_velocity(&state, true).unwrap(), 0.0);
        let local = shifted_kinematic_state(state, acceleration, epoch + 0.125).unwrap();
        assert_eq!(
            [local.x, local.y, local.z],
            [1.0 - 0.5 * 0.125 * 0.125, 0.125, 0.0]
        );
        assert_eq!([local.vx, local.vy, local.vz], [-0.125, 1.0, 0.0]);
        assert_eq!(
            shifted_kinematic_state(state, None, epoch + 0.125).unwrap(),
            state.shifted(epoch + 0.125).unwrap()
        );
    }

    #[test]
    fn analytic_range_derivative_has_radial_and_tangent_invariants() {
        let radial = PropagationResult::from_state(
            CartesianState {
                position: [3.0, 4.0, 12.0],
                velocity: [0.3, 0.4, 1.2],
            },
            J2000,
        )
        .unwrap();
        assert!((heliocentric_radial_velocity(&radial, false).unwrap() - 1.3).abs() < 1e-14);
        let tangent = state();
        assert_eq!(heliocentric_radial_velocity(&tangent, false).unwrap(), 0.0);
        let speed = 1.3;
        let beta = speed / C_LIGHT_AU_DAY;
        let delayed = speed * (1.0 - beta + beta * beta);
        assert!((heliocentric_radial_velocity(&radial, true).unwrap() - delayed).abs() < 1e-14);
    }
    struct Observer;
    impl AssistMinorEphemeris for Observer {
        fn earth_position(&self, _: f64) -> Result<[f64; 3]> {
            Ok([1.0, 0.0, 0.0])
        }
        fn sun_position(&self, _: f64) -> Result<[f64; 3]> {
            Ok([0.0; 3])
        }
        fn earth_velocity(&self, _: f64) -> Result<[f64; 3]> {
            Ok([0.0, 0.01, 0.0])
        }
    }

    #[derive(Default)]
    struct RecordingObserver {
        reads: RefCell<Vec<(&'static str, f64)>>,
        fail: Option<&'static str>,
    }

    impl RecordingObserver {
        fn read(&self, kind: &'static str, tt: f64, value: [f64; 3]) -> Result<[f64; 3]> {
            self.reads.borrow_mut().push((kind, tt));
            if self.fail == Some(kind) {
                return Err(Error::SpkBodyNotFound { target: 399 });
            }
            Ok(value)
        }
    }

    impl AssistMinorEphemeris for RecordingObserver {
        fn earth_position(&self, tt: f64) -> Result<[f64; 3]> {
            self.read("earth", tt, [0.0, 1.0, 0.0])
        }
        fn sun_position(&self, tt: f64) -> Result<[f64; 3]> {
            self.read("sun", tt, [0.0; 3])
        }
        fn earth_velocity(&self, tt: f64) -> Result<[f64; 3]> {
            self.read("velocity", tt, [0.0; 3])
        }
    }

    #[test]
    fn true_position_uses_the_geometric_earth_minus_sun_observer() {
        let value =
            assist_minor_sample(&state(), J2000, FLG_TRUEPOS | FLG_NONUT, &Observer).unwrap();
        assert!(
            value[0].abs() < 1e-10 && value[1].abs() < 1e-10,
            "{value:?}"
        );
        assert!((value[2] - 1.0).abs() < 1e-14, "{value:?}");
    }
    #[test]
    fn heliocentric_true_position_keeps_the_initial_geometric_radius() {
        let value = assist_minor_sample(
            &state(),
            J2000,
            FLG_HELCTR | FLG_TRUEPOS | FLG_NONUT,
            &Observer,
        )
        .unwrap();
        assert!(
            value[0].abs() < 1e-10 && value[1].abs() < 1e-10,
            "{value:?}"
        );
        assert!((value[2] - 2.0).abs() < 1e-14, "{value:?}");
    }

    #[test]
    fn orbital_basis_roundtrip_recovers_closed_icrf_position_and_velocity() {
        struct Supplier(AtomicUsize);
        impl MinorAssistAcquisition for Supplier {
            fn contains(&self, _: i32) -> Result<bool> {
                Ok(true)
            }
            fn available(&self) -> Result<bool> {
                Ok(true)
            }
            fn state_at(&self, _: i32, tt: f64) -> Result<PropagationResult> {
                self.0.fetch_add(1, Ordering::Relaxed);
                // A known ICRF state expressed in the documented public
                // orbital frame, independently of the output-frame helper.
                PropagationResult::from_state(
                    CartesianState {
                        position: crate::spk::icrs_to_ecliptic_j2000([0.0, 2.0, 0.0]),
                        velocity: crate::spk::icrs_to_ecliptic_j2000([0.01, 0.02, 0.03]),
                    },
                    tt,
                )
            }
        }
        let session = Session::from_context(crate::context::EphemerisContext::new());
        let supplier = Supplier(AtomicUsize::new(0));
        let value = session
            .assist_minor_equatorial(
                &supplier,
                crate::constants::CERES,
                J2000,
                FLG_HELCTR | FLG_TRUEPOS | FLG_J2000 | FLG_SPEED | crate::constants::FLG_ICRS,
            )
            .unwrap();
        assert_eq!(supplier.0.load(Ordering::Relaxed), 1);
        assert!((value[0] - 90.0).abs() < 1e-10, "{value:?}");
        assert!(value[1].abs() < 1e-10, "{value:?}");
        assert!((value[2] - 2.0).abs() < 1e-14, "{value:?}");
        // At p=(0,2,0), dRA=-vx/2, dDec=vz/2, dr=vy.
        for (actual, expected) in
            value[3..]
                .iter()
                .zip([(-0.005_f64).to_degrees(), 0.015_f64.to_degrees(), 0.02])
        {
            assert!((actual - expected).abs() <= expected.abs() * 1e-6);
        }
    }

    #[test]
    fn complete_ray_bridge_covers_true_apparent_ecliptic_and_equatorial_frames() {
        let central = PropagationResult::from_state(
            CartesianState {
                position: crate::spk::icrs_to_ecliptic_j2000([0.0, 2.0, 0.0]),
                velocity: crate::spk::icrs_to_ecliptic_j2000([0.0, 0.02, 0.0]),
            },
            J2000,
        )
        .unwrap();
        let engine_obliquity = crate::positions::mean_obliquity_iau2006_deg(J2000).unwrap();
        let angular_limit = 0.001 / 3600.0;
        for modifiers in [
            FLG_HELCTR | FLG_TRUEPOS,
            FLG_HELCTR,
            FLG_TRUEPOS,
            0,
            FLG_NOABERR,
        ] {
            for frame in [0, FLG_J2000, FLG_EQUATORIAL, FLG_J2000 | FLG_EQUATORIAL] {
                let observer = RecordingObserver::default();
                let flags = modifiers | frame | FLG_NONUT;
                let sample = assist_minor_sample(&central, J2000, flags, &observer).unwrap();
                assert!(
                    (sample[0] - 90.0).abs() <= angular_limit,
                    "{flags}: {sample:?}"
                );
                assert!(
                    (sample[1] + engine_obliquity).abs() <= angular_limit,
                    "{flags}: {sample:?}"
                );
                let radius = if modifiers & FLG_HELCTR != 0 {
                    2.0
                } else {
                    1.0
                };
                // For this radial constant-velocity state, the two held-state
                // delay updates have the closed result r*(1-beta+beta²).
                let beta = 0.02 / C_LIGHT_AU_DAY;
                let expected_radius = if modifiers & FLG_TRUEPOS != 0 {
                    radius
                } else {
                    radius * (1.0 - beta + beta * beta)
                };
                assert!(
                    (sample[2] - expected_radius).abs() < 1e-13,
                    "{flags}: {sample:?}"
                );
                let output = crate::positions::horizons_output::local_convert(
                    [sample[0], sample[1], sample[2], 0.0, 0.0, 0.0],
                    J2000,
                    flags,
                    false,
                )
                .unwrap();
                if frame & FLG_EQUATORIAL != 0 {
                    assert!(
                        (output[0] - 90.0).abs() <= angular_limit,
                        "{flags}: {output:?}"
                    );
                    assert!(output[1].abs() <= angular_limit, "{flags}: {output:?}");
                }
                let expected_reads = if modifiers & FLG_HELCTR != 0 {
                    vec![]
                } else if modifiers & (FLG_TRUEPOS | FLG_NOABERR) != 0 {
                    vec![("earth", J2000), ("sun", J2000)]
                } else {
                    vec![("earth", J2000), ("sun", J2000), ("velocity", J2000)]
                };
                assert_eq!(*observer.reads.borrow(), expected_reads);
            }
        }
    }

    #[test]
    fn basis_bridge_preserves_initial_validation_and_observer_error_order() {
        let observer = RecordingObserver {
            fail: Some("earth"),
            ..Default::default()
        };
        let mut invalid = state();
        invalid.x = f64::NAN;
        assert!(matches!(
            assist_minor_sample(&invalid, J2000, FLG_TRUEPOS, &observer),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(observer.reads.borrow().is_empty());
        for (failure, expected_reads) in [
            ("earth", vec![("earth", J2000)]),
            ("sun", vec![("earth", J2000), ("sun", J2000)]),
            (
                "velocity",
                vec![("earth", J2000), ("sun", J2000), ("velocity", J2000)],
            ),
        ] {
            let observer = RecordingObserver {
                fail: Some(failure),
                ..Default::default()
            };
            assert!(matches!(
                assist_minor_sample(&state(), J2000, 0, &observer),
                Err(Error::SpkBodyNotFound { target: 399 })
            ));
            assert_eq!(*observer.reads.borrow(), expected_reads);
        }
    }

    #[test]
    fn owned_emission_samples_curvature_and_moving_sun_finite_updates() {
        let radius = 2.0;
        let acceleration = 0.001;
        let sun_speed = 0.005;
        let central = CartesianState {
            position: [radius, 0.0, 0.0],
            velocity: [0.0; 3],
        };
        let observer = CartesianState {
            position: [0.0; 3],
            velocity: [sun_speed, 0.0, 0.0],
        };
        let mut offsets = Vec::new();
        let mut source = |offset: f64| {
            offsets.push(offset);
            Ok(CartesianState {
                position: [radius, 0.0, acceleration * offset * offset / 2.0],
                velocity: [0.0, 0.0, acceleration * offset],
            })
        };
        let output = finite_emission_ray(central, observer, &mut source, true).unwrap();
        let first_delay = radius / C_LIGHT_AU_DAY;
        let first_z = acceleration * first_delay * first_delay / 2.0;
        let first_radius = radius.hypot(first_z);
        let second_delay = first_radius / C_LIGHT_AU_DAY;
        let first_delay_rate = -sun_speed / C_LIGHT_AU_DAY;
        let first_z_rate = -acceleration * first_delay * (1.0 - first_delay_rate);
        let second_delay_rate =
            (-radius * sun_speed + first_z * first_z_rate) / first_radius / C_LIGHT_AU_DAY;
        let expected_z = acceleration * second_delay * second_delay / 2.0;
        let expected_z_rate = -acceleration * second_delay * (1.0 - second_delay_rate);
        assert_eq!(output.position, [radius, 0.0, expected_z]);
        assert_eq!(offsets, [-first_delay, -second_delay]);
        assert!((output.velocity[0] + sun_speed).abs() < 1e-16);
        assert!((output.velocity[2] - expected_z_rate).abs() < 1e-16);
    }

    #[test]
    fn owned_emission_domain_failure_stops_before_another_source_read() {
        let central = CartesianState {
            position: [2.0, 0.0, 0.0],
            velocity: [0.0; 3],
        };
        let observer = CartesianState {
            position: [0.0; 3],
            velocity: [0.0; 3],
        };
        let first_delay = 2.0 / C_LIGHT_AU_DAY;
        for fail_at_first in [true, false] {
            let mut reads = 0;
            let start = if fail_at_first {
                -first_delay / 2.0
            } else {
                -first_delay
            };
            let mut source = |offset: f64| {
                reads += 1;
                if offset < start {
                    return Err(Error::SpkOutOfRange {
                        target: 999,
                        et: offset * 86400.0,
                        start: start * 86400.0,
                        end: 0.0,
                    });
                }
                Ok(CartesianState {
                    position: [2.0, 0.0, offset * offset / 2.0],
                    velocity: [0.0, 0.0, offset],
                })
            };
            assert!(matches!(
                finite_emission_ray(central, observer, &mut source, true),
                Err(Error::SpkOutOfRange { target: 999, .. })
            ));
            assert_eq!(reads, if fail_at_first { 1 } else { 2 });
        }
    }

    #[test]
    fn owned_finite_delay_derivative_carries_source_and_observer_clock_rates() {
        let source_rate = 1.002;
        let observer_rate = 0.999;
        let speed = 0.02;
        let observer_speed = 0.005;
        let central = CartesianState {
            position: [2.0, 0.0, 0.0],
            velocity: [speed * source_rate, 0.0, 0.0],
        };
        let observer = CartesianState {
            position: [0.0; 3],
            velocity: [observer_speed * observer_rate, 0.0, 0.0],
        };
        let mut source = |offset: f64| {
            Ok(CartesianState {
                position: [2.0 + speed * source_rate * offset, 0.0, 0.0],
                velocity: central.velocity,
            })
        };
        let value = finite_emission_ray(central, observer, &mut source, true).unwrap();
        let beta = speed * source_rate / C_LIGHT_AU_DAY;
        let relative = speed * source_rate - observer_speed * observer_rate;
        let expected = relative * (1.0 - beta + beta * beta);
        assert!((value.velocity[0] - expected).abs() < 1e-16);
        let implicit = relative / (1.0 + beta);
        assert_ne!(value.velocity[0].to_bits(), implicit.to_bits());
    }

    #[test]
    fn one_integration_supplies_the_tangent_velocity_stencil() {
        struct Supplier(AtomicUsize);
        impl MinorAssistAcquisition for Supplier {
            fn contains(&self, _: i32) -> Result<bool> {
                Ok(true)
            }
            fn available(&self) -> Result<bool> {
                Ok(true)
            }
            fn state_at(&self, _: i32, _: f64) -> Result<PropagationResult> {
                self.0.fetch_add(1, Ordering::Relaxed);
                Ok(state())
            }
        }
        let supplier = Supplier(AtomicUsize::new(0));
        let session = Session::from_context(crate::context::EphemerisContext::new());
        let value = session
            .assist_minor_ecliptic(
                &supplier,
                crate::constants::CERES,
                J2000,
                FLG_HELCTR | FLG_TRUEPOS | FLG_NONUT | FLG_SPEED,
            )
            .unwrap();
        // The result uses the ecliptic of date, whose precession contributes
        // to angular speed independently of the body's Cartesian motion.
        let dt = 0.5;
        let before = J2000 - dt;
        let after = J2000 + dt;
        let frame_before =
            crate::hypothetical::precess_ecliptic_source_angles(0.0, 0.0, J2000, before)
                .unwrap()
                .0;
        let frame_after =
            crate::hypothetical::precess_ecliptic_source_angles(0.0, 0.0, J2000, after)
                .unwrap()
                .0;
        let frame_rate =
            ((frame_after - frame_before + 180.0).rem_euclid(360.0) - 180.0) / (2.0 * dt);
        let angular_speed = (0.01_f64 / 2.0).to_degrees() + frame_rate;
        assert!((value[3] - angular_speed).abs() < 1e-6, "{value:?}");
        assert!(value[4].abs() < 1e-6 && value[5].abs() < 1e-10, "{value:?}");
        assert_eq!(supplier.0.load(Ordering::Relaxed), 1);
    }

    #[test]
    fn equatorial_rates_include_the_time_dependent_obliquity() {
        struct Supplier(AtomicUsize);
        impl MinorAssistAcquisition for Supplier {
            fn contains(&self, _: i32) -> Result<bool> {
                Ok(true)
            }
            fn available(&self) -> Result<bool> {
                Ok(true)
            }
            fn state_at(&self, _: i32, tt: f64) -> Result<PropagationResult> {
                self.0.fetch_add(1, Ordering::Relaxed);
                PropagationResult::from_state(
                    CartesianState {
                        position: [1.5, 1.0, 0.7],
                        velocity: [0.002, -0.003, 0.004],
                    },
                    tt,
                )
            }
        }
        let session = Session::from_context(crate::context::EphemerisContext::new());
        let epoch = 2460310.5;
        for modifiers in [0, FLG_NONUT, FLG_J2000, FLG_SIDEREAL] {
            let supplier = Supplier(AtomicUsize::new(0));
            let flags = FLG_HELCTR | FLG_TRUEPOS | FLG_EQUATORIAL | FLG_SPEED | modifiers;
            let output = session
                .assist_minor_equatorial(&supplier, 17, epoch, flags)
                .unwrap();
            assert_eq!(supplier.0.load(Ordering::Relaxed), 1);
            let central = supplier.state_at(17, epoch).unwrap();
            let position = |offset| {
                let tt = epoch + offset;
                let state = central.shifted(tt).unwrap();
                let value = session
                    .assist_minor_position_from_state(&state, tt, flags)
                    .unwrap();
                crate::positions::horizons_output::local_convert(
                    [value[0], value[1], value[2], 0.0, 0.0, 0.0],
                    tt,
                    flags & !FLG_SPEED,
                    false,
                )
                .unwrap()
            };
            let h = 1.0 / 32.0;
            let before = position(-h);
            let after = position(h);
            let before2 = position(-2.0 * h);
            let after2 = position(2.0 * h);
            let before3 = position(-3.0 * h);
            let after3 = position(3.0 * h);
            for component in 0..2 {
                let sixth = (45.0 * (after[component] - before[component])
                    - 9.0 * (after2[component] - before2[component])
                    + after3[component]
                    - before3[component])
                    / (60.0 * h);
                // Fixed PLAN speed bound against an independent higher-order
                // position derivative, rather than a self-derived expected value.
                let limit = if sixth.abs() <= 1e-6 {
                    1e-6
                } else {
                    sixth.abs() * 1e-6
                };
                assert!(
                    (output[component + 3] - sixth).abs() <= limit,
                    "flags={flags}, slot={component}: {} vs {sixth}",
                    output[component + 3]
                );
            }
        }
    }

    #[test]
    fn projected_range_retains_a_small_radial_component_without_one_second_cancellation() {
        struct LinearObserver;
        impl AssistMinorEphemeris for LinearObserver {
            fn earth_position(&self, tt: f64) -> Result<[f64; 3]> {
                // ICRF and ecliptic-J2000 share the x axis, so the
                // independently stated relative-state invariant needs no
                // obliquity conversion in its expected vector.
                Ok([1.0 + 0.01 * (tt - J2000), 0.0, 0.0])
            }
            fn sun_position(&self, _: f64) -> Result<[f64; 3]> {
                Ok([0.0; 3])
            }
            fn earth_velocity(&self, _: f64) -> Result<[f64; 3]> {
                Ok([0.01, 0.0, 0.0])
            }
        }
        let central = PropagationResult::from_state(
            CartesianState {
                position: [31.0, 40.0, 0.0],
                velocity: [0.00206, 0.00608, 0.0],
            },
            J2000,
        )
        .unwrap();
        let sample = |offset| {
            let tt = J2000 + offset;
            assist_minor_sample(
                &central.shifted(tt).unwrap(),
                tt,
                FLG_TRUEPOS | FLG_NONUT,
                &LinearObserver,
            )
            .unwrap()
        };
        let h = POSITION_RATE_STEP_DAYS;
        let radial =
            fourth_order_position_rates(sample(-h), sample(h), sample(-2.0 * h), sample(2.0 * h))
                [2];
        // Independently, q=(30,40,0), q'=(-0.00794,0.00608,0),
        // so d|q|/dt=q·q'/|q|=0.0001 AU/day. Keep the PLAN 1e-6
        // relative speed budget, not a norm-only substitute for this component.
        assert!((radial - 1e-4).abs() <= 1e-10, "{radial:.17e}");
    }

    #[test]
    fn fourth_order_position_rates_recover_cubic_drift_and_exact_constant_range() {
        let sample = |t: f64| {
            [
                (0.005 + 0.5 * t + 0.1 * t.powi(2) + 0.01 * t.powi(3)).rem_euclid(360.0),
                0.0002 - 4e-6 * t + 0.1 * t.powi(2) - 0.2 * t.powi(3),
                1.0,
            ]
        };
        let h = POSITION_RATE_STEP_DAYS;
        let rates =
            fourth_order_position_rates(sample(-h), sample(h), sample(-2.0 * h), sample(2.0 * h));
        assert!((rates[0] - 0.5).abs() < 1e-11);
        assert!((rates[1] - -4e-6).abs() < 1e-17);
        assert_eq!(rates[2], 0.0);
        for jd in [-1_000_000.5_f64, J2000, 8_000_000.5] {
            assert_eq!((jd + h) - (jd - h), 2.0 * h);
        }
    }
}
