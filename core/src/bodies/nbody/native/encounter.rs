// SPDX-License-Identifier: AGPL-3.0
//! Local sampling time of prescribed Newtonian sources, in common units.
//!
//! A source can change direction on the separation/relative-speed crossing
//! time even when its weak acceleration is hidden by a large background force.
//! Sampling no less often than that time complements a total-force derivative
//! controller. This geometric cap neither changes forces nor guarantees a
//! global trajectory error. It uses no body-specific or fitted constant.

use super::forces::Source;
use super::state::{NumericalError, State};

/// Minimum r/|v_rel| across positive-GM, physically nonself sources. Stationary
/// relative motion imposes no cap. All supplied states must be finite, including
/// excluded rows; selected ephemeris ownership guarantees unique physical IDs.
pub(super) fn maximum_step(
    target: State,
    sources: &[Source],
    exclude: Option<i32>,
) -> Result<f64, NumericalError> {
    target.checked()?;
    let mut limit = f64::INFINITY;
    for source in sources {
        let point = source.point;
        if !point.gravitational_parameter.is_finite()
            || !point
                .position
                .into_iter()
                .chain(source.velocity)
                .all(f64::is_finite)
        {
            return Err(NumericalError::NonFinite);
        }
        if point.gravitational_parameter < 0.0 {
            return Err(NumericalError::InvalidMass { body: point.body });
        }
        if Some(point.body) == exclude || point.gravitational_parameter == 0.0 {
            continue;
        }
        let offset: [f64; 3] = std::array::from_fn(|i| point.position[i] - target.position[i]);
        let velocity: [f64; 3] = std::array::from_fn(|i| source.velocity[i] - target.velocity[i]);
        let separation = offset[0].hypot(offset[1]).hypot(offset[2]);
        let speed = velocity[0].hypot(velocity[1]).hypot(velocity[2]);
        if !separation.is_finite() || !speed.is_finite() {
            return Err(NumericalError::NonFinite);
        }
        if separation == 0.0 {
            return Err(NumericalError::SingularForce { body: point.body });
        }
        if speed > 0.0 {
            limit = limit.min(separation / speed);
        }
    }
    Ok(limit)
}

#[cfg(test)]
#[path = "encounter_tests.rs"]
mod tests;
