// SPDX-License-Identifier: AGPL-3.0
//! Newtonian acceleration from prescribed point masses in a common frame.
//!
//! Newton's inverse-square law is summed as mu/r^2 times a unit direction.
//! This avoids cubing a very small or large radius. The caller supplies all
//! states and gravitational parameters in consistent units; no ephemeris,
//! identifier order, astronomical constant or third-party table is embedded.

use super::state::{NumericalError, Sum};

/// One physical source with position and GM, not an integrator particle.
#[derive(Clone, Copy, Debug)]
pub(super) struct PointMass {
    pub(super) body: i32,
    pub(super) position: [f64; 3],
    pub(super) gravitational_parameter: f64,
}

/// Sum direct accelerations, excluding only the explicitly identified self body.
///
/// This kernel does not decide identity from proximity, disable an entire force
/// family, soften singularities or manufacture a missing source. Validate the
/// complete source selection, including an excluded row, before evaluation.
pub(super) fn acceleration(
    target: [f64; 3],
    sources: &[PointMass],
    exclude: Option<i32>,
) -> Result<[f64; 3], NumericalError> {
    acceleration_from(target, sources.iter().copied(), exclude)
}

/// The same Newtonian kernel over a borrowed frame projection. Cloning this
/// iterator copies only its cursor; no projected source vector is allocated.
pub(super) fn acceleration_from(
    target: [f64; 3],
    sources: impl Iterator<Item = PointMass> + Clone,
    exclude: Option<i32>,
) -> Result<[f64; 3], NumericalError> {
    if !target.into_iter().all(f64::is_finite) {
        return Err(NumericalError::NonFinite);
    }
    for (index, source) in sources.clone().enumerate() {
        if !source.gravitational_parameter.is_finite()
            || !source.position.into_iter().all(f64::is_finite)
        {
            return Err(NumericalError::NonFinite);
        }
        if source.gravitational_parameter < 0.0 {
            return Err(NumericalError::InvalidMass { body: source.body });
        }
        if sources
            .clone()
            .take(index)
            .any(|other| other.body == source.body)
        {
            return Err(NumericalError::DuplicateSource { body: source.body });
        }
    }
    let mut total = [Sum::default(); 3];
    for source in sources {
        if Some(source.body) == exclude || source.gravitational_parameter == 0.0 {
            continue;
        }
        let offset: [f64; 3] = std::array::from_fn(|axis| source.position[axis] - target[axis]);
        let radius = offset[0].hypot(offset[1]).hypot(offset[2]);
        if !radius.is_finite() {
            return Err(NumericalError::NonFinite);
        }
        if radius == 0.0 {
            return Err(NumericalError::SingularForce { body: source.body });
        }
        let magnitude = (source.gravitational_parameter / radius) / radius;
        for axis in 0..3 {
            total[axis].add(magnitude * (offset[axis] / radius));
        }
    }
    let result = total.map(Sum::value);
    if result.into_iter().all(f64::is_finite) {
        Ok(result)
    } else {
        Err(NumericalError::NonFinite)
    }
}
