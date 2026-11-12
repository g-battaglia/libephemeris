//! Rate calculations for the stored core without a surface observer.
//!
//! The caller selects analytic eligibility and acquires acceleration before
//! entering either route. Projection, motion caching and final trace publication
//! remain in that caller; the distinct topocentric model is not handled here.

use super::{
    AngularModelSample, CentralAngularPlace, CentralAngularQuery, CentralStateOrder,
    LightTimeSolution, OwnedApparentMotion, PpnScalarSamples, SealedServing, StateProvider,
    StencilCenter, TwoPartEpoch, angular_geometry_support, central_angular_sample,
    compensated_dot3, local_model_velocity, serving_position_core_recorded, shift_state_word_parts,
};
use crate::constants::{FLG_SPEED, FLG_SPEED3};
use crate::errors::Result;
use crate::positions::correction_velocity::{ExpandedCenter, apparent_motion};
use crate::positions::derivative::AngularWord as Word;
use crate::positions::state_reads::ApparentRateReads;

/// Borrow the solved propagation and its matching central observation geometry.
/// Epochs are TT Julian Days; stencil positions/velocities are ICRS AU/AU per
/// TT day. Original split parts remain separate from the scalar sampling clock.
/// No provider reads or epoch preparation occur when this record is assembled.
pub(super) struct RateObservation<'a> {
    pub(super) jd_tt: f64,
    pub(super) obs_split: Option<(f64, f64)>,
    pub(super) observer_body: Option<i32>,
    pub(super) light_time: &'a LightTimeSolution,
    pub(super) stencil: &'a StencilCenter,
}

/// Differentiate the solved propagation and enabled apparent corrections.
///
/// Acceleration has already passed the caller's lazy analytic-kinematics gate.
/// Retain relative-velocity lookup before Word reads, including its errors,
/// and use the original emission Jacobian without another propagation solve.
pub(super) fn analytic_apparent_motion<S: ApparentRateReads + ?Sized>(
    source: &S,
    rate_observation: &RateObservation<'_>,
    acceleration: [f64; 3],
    ppn_scalar_samples: Option<&PpnScalarSamples>,
) -> Result<OwnedApparentMotion> {
    let jd_tt = rate_observation.jd_tt;
    let obs_split = rate_observation.obs_split;
    let observer_body = rate_observation.observer_body;
    let lt = rate_observation.light_time;
    let stencil_center = rate_observation.stencil;
    let body_id = stencil_center.target_body_id;
    let obs_pos = stencil_center.observer;
    let obs_vel = stencil_center.observer_vel;

    // Keep the small emission Jacobian separate from unity. Rounding
    // `1 - delay_rate` first loses barycentric velocity bits before
    // the observer subtraction, especially for close relative rays.
    let mut analytic_center = *stencil_center;
    let target_velocity = lt.retarded_state().1;
    let observation = match obs_split {
        Some((hi, lo)) => TwoPartEpoch { hi, lo },
        None => TwoPartEpoch::split(jd_tt)?,
    };
    let emission = observation.add_offset(-lt.light_time_days())?;
    let paired = match observer_body {
        Some(observer) => {
            source.eval_relative_velocity_two_part(body_id, emission, observer, observation)?
        }
        None => None,
    };
    analytic_center.geo_vel = std::array::from_fn(|i| match paired {
        Some(relative) => (-lt.light_time_rate()).mul_add(target_velocity[i], relative[i]),
        None => compensated_dot3(
            [1.0, -1.0, -lt.light_time_rate()],
            [
                target_velocity[i],
                stencil_center.observer_vel[i],
                target_velocity[i],
            ],
        ),
    });
    let expand = |pair: [f64; 2]| Word {
        hi: pair[0],
        lo: pair[1],
    };
    let (target_position, target_rate) =
        source.geometric_state_word_parts_two_part(body_id, emission)?;
    let (observer_position, observer_rate) = match observer_body {
        Some(observer) => {
            let words = source.geometric_state_word_parts_two_part(observer, observation)?;
            let central = source.geometric_state_two_part(observer, observation)?;
            shift_state_word_parts(words, (central.pos(), central.vel()), (obs_pos, obs_vel))
        }
        None => (obs_pos.map(|v| [v, 0.0]), obs_vel.map(|v| [v, 0.0])),
    };
    let observer_position = observer_position.map(expand);
    let observer_rate = observer_rate.map(expand);
    let target_position = target_position.map(expand);
    let target_rate = target_rate.map(expand);
    // The exact existing observation/emission parts and final light-time
    // rate remain authoritative. No new propagation iteration is added.
    let expanded_center = ExpandedCenter {
        position: std::array::from_fn(|i| target_position[i] - observer_position[i]),
        velocity: std::array::from_fn(|i| {
            target_rate[i] * (Word::from(1.0) - Word::from(lt.light_time_rate())) - observer_rate[i]
        }),
        observer: observer_position,
        observer_velocity: observer_rate,
    };
    apparent_motion(
        source,
        &analytic_center,
        jd_tt,
        lt.light_time_days(),
        lt.light_time_rate(),
        acceleration,
        obs_split.is_some(),
        expanded_center,
        ppn_scalar_samples,
    )
}

/// Differentiate actual sampled places on their owned source domains.
///
/// Immutable central reuse still occurs only inside the zero-offset callback;
/// domain validation and trace publication retain that callback's timing. Opaque
/// providers and nonidentical split epochs execute the complete recursive place.
pub(super) fn sampled_apparent_rates<P: StateProvider>(
    serving: &SealedServing<'_, P>,
    observation: &RateObservation<'_>,
    calc_flags: i32,
    want_xyz: bool,
    order: CentralStateOrder,
    central: Option<&CentralAngularPlace<'_>>,
) -> Result<(f64, f64, f64)> {
    let jd_tt = observation.jd_tt;
    let obs_split = observation.obs_split;
    let observer_body = observation.observer_body;
    let lt = observation.light_time;
    let body_id = observation.stencil.target_body_id;
    let angular_at = |offset: f64| {
        let epoch = jd_tt + offset;
        let split = obs_split.map(|_| {
            let whole = epoch.floor();
            (whole, epoch - whole)
        });
        // A supplied split can contain a sub-ulp tail lost by the sample's
        // canonical reconstruction. Such a central request must still run
        // the original calculation; equality of the summed JD is not enough.
        if offset == 0.0 {
            if let Some(central) = central {
                let query = CentralAngularQuery::new(epoch, split, body_id, calc_flags, order);
                if let Some(sample) = central_angular_sample(serving, central, query)? {
                    return Ok(sample);
                }
            }
        }
        let capture = std::cell::Cell::new(AngularModelSample::default());
        serving_position_core_recorded(
            serving,
            epoch,
            body_id,
            calc_flags & !(FLG_SPEED | FLG_SPEED3),
            true,
            split,
            order,
            Some(&capture),
        )
        .map(|_| capture.get())
    };
    let support = angular_geometry_support(
        serving.provider,
        body_id,
        observer_body,
        jd_tt,
        lt.retarded_jd_tt(),
    )?;
    local_model_velocity(
        serving,
        observation.stencil,
        lt.geo_vel(),
        jd_tt,
        calc_flags,
        want_xyz,
        obs_split.is_some(),
        None,
        &angular_at,
        support,
    )
}

#[cfg(test)]
#[path = "rate_read_tests.rs"]
mod read_contract_tests;
