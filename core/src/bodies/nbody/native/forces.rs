// SPDX-License-Identifier: AGPL-3.0
//! Independently derived prescribed-ephemeris force equations.
//!
//! All states share SSB ICRF axes, AU, AU/day and TDB days. GM is AU^3/day^2.
//! Newtonian sources remain distinct from the selected post-Newtonian sources;
//! a small-body perturber is never silently given a relativistic acceleration.
//! Constants, orientations and default model selections belong to the retained
//! ephemeris/configuration boundary, not this mathematical kernel.
//!
//! PPN: Will & Nordtvedt (1972), Moyer (2003), and the published test-particle
//! equation in Holman et al. (2023), section III.2, beta=gamma=1. Zonal terms
//! are minus the gradient of the Legendre gravitational potential. Radial,
//! transverse and normal nongravitational terms follow Marsden et al. (1973).
//! No software implementation, coefficient table or fitted output is an input.

use super::newtonian::{self, PointMass};
use super::state::{NumericalError, State, Sum};

/// One immutable massive-body snapshot. `post_newtonian_acceleration` is Some
/// only for selected major-body sources with an actual ephemeris acceleration.
#[derive(Clone, Copy, Debug)]
pub(super) struct Source {
    pub(super) point: PointMass,
    pub(super) velocity: [f64; 3],
    pub(super) post_newtonian_acceleration: Option<[f64; 3]>,
}

impl Source {
    fn checked(self) -> Result<Self, NumericalError> {
        if !self.velocity.into_iter().all(f64::is_finite)
            || self
                .post_newtonian_acceleration
                .is_some_and(|value| !value.into_iter().all(f64::is_finite))
        {
            return Err(NumericalError::NonFinite);
        }
        Ok(self)
    }
}

/// Which sources appear in the PPN outer sum. Both choices retain all selected
/// major bodies in the potentials; asteroid Newtonian sources stay separate.
#[derive(Clone, Copy, Debug)]
pub(super) enum RelativisticSources {
    Sun,
    /// Independent wider-source witness; the public API selects Sun only.
    #[cfg(test)]
    MajorBodies,
}

/// Parameter-free general relativity (beta=gamma=1), with c in AU/day.
#[derive(Clone, Copy, Debug)]
pub(super) struct Eih {
    pub(super) speed_of_light: f64,
    pub(super) outer: RelativisticSources,
}

/// One axisymmetric potential, degree2 through degree4. Coefficients are
/// unnormalized Jn, not normalized Cn0. The pole is normalized on construction.
#[derive(Clone, Copy, Debug)]
pub(super) struct Zonals {
    body: i32,
    reference_radius: f64,
    pole: [f64; 3],
    coefficients: [f64; 3],
}

impl Zonals {
    pub(super) fn new(
        body: i32,
        reference_radius: f64,
        pole: [f64; 3],
        coefficients: [f64; 3],
    ) -> Result<Self, NumericalError> {
        if !reference_radius.is_finite()
            || !pole.into_iter().chain(coefficients).all(f64::is_finite)
        {
            return Err(NumericalError::NonFinite);
        }
        let norm = length(pole)?;
        if reference_radius <= 0.0 || norm == 0.0 {
            return Err(NumericalError::InvalidSettings);
        }
        Ok(Self {
            body,
            reference_radius,
            pole: pole.map(|value| value / norm),
            coefficients,
        })
    }

    /// Change only the independently supplied inertial pole at a force epoch.
    /// Body identity, selected-file radius and all Jn coefficients are retained.
    #[cfg(feature = "nbody")]
    pub(super) fn with_pole(self, pole: [f64; 3]) -> Result<Self, NumericalError> {
        Self::new(self.body, self.reference_radius, pole, self.coefficients)
    }

    fn acceleration(self, position: [f64; 3], source: Source) -> Result<[f64; 3], NumericalError> {
        let displacement = subtract(position, source.point.position);
        let radius = distance(displacement, self.body)?;
        let direction = displacement.map(|value| value / radius);
        let sine_latitude = dot(direction, self.pole);
        let s = sine_latitude;
        let square = s * s;
        // Pn and Pn' are elementary Legendre polynomials, not an upstream table.
        let polynomial = [
            ((3.0 * square - 1.0) / 2.0, 3.0 * s),
            (s * (5.0 * square - 3.0) / 2.0, (15.0 * square - 3.0) / 2.0),
            (
                (35.0 * square * square - 30.0 * square + 3.0) / 8.0,
                s * (35.0 * square - 15.0) / 2.0,
            ),
        ];
        let ratio = self.reference_radius / radius;
        let monopole = (source.point.gravitational_parameter / radius) / radius;
        let mut power = ratio;
        let mut total = [Sum::default(); 3];
        for (index, (value, derivative)) in polynomial.into_iter().enumerate() {
            power *= ratio;
            let coefficient = self.coefficients[index];
            if coefficient == 0.0 {
                continue;
            }
            let radial = (index + 3) as f64 * value + s * derivative;
            let magnitude = monopole * coefficient * power;
            for axis in 0..3 {
                total[axis]
                    .add(magnitude * (radial * direction[axis] - derivative * self.pole[axis]));
            }
        }
        finite_vector(total.map(Sum::value))
    }
}

/// Distance law for the given A1/A2/A3 accelerations. Inverse-square uses AU;
/// the Marsden scale distance is also in AU. All other parameters are unitless.
#[derive(Clone, Copy, Debug)]
pub(super) enum DistanceLaw {
    InverseSquare,
    /// General scientific-law witness, not an unexposed runtime option.
    #[cfg(test)]
    Marsden {
        alpha: f64,
        scale: f64,
        m: f64,
        n: f64,
        k: f64,
    },
}

impl DistanceLaw {
    fn value(self, radius: f64) -> Result<f64, NumericalError> {
        let value = match self {
            Self::InverseSquare => (1.0 / radius) / radius,
            #[cfg(test)]
            Self::Marsden {
                alpha,
                scale,
                m,
                n,
                k,
            } => {
                if ![alpha, scale, m, n, k].into_iter().all(f64::is_finite) {
                    return Err(NumericalError::NonFinite);
                }
                if alpha <= 0.0 || scale <= 0.0 {
                    return Err(NumericalError::InvalidSettings);
                }
                // Evaluate through logs to avoid an overflowing (r/r0)^n.
                // This is algebraically the same distance law, not a changed
                // tail or saturation. Ordinary IEEE underflow to zero is kept.
                let logarithm = radius.ln() - scale.ln();
                let y = n * logarithm;
                if !y.is_finite() {
                    return Err(NumericalError::NonFinite);
                }
                let log_one_plus = if y > 0.0 {
                    y + (-y).exp().ln_1p()
                } else {
                    y.exp().ln_1p()
                };
                (alpha.ln() - m * logarithm - k * log_one_plus).exp()
            }
        };
        if value.is_finite() {
            Ok(value)
        } else {
            Err(NumericalError::NonFinite)
        }
    }
}

/// Physical radial/transverse/normal accelerations in AU/day^2.
#[derive(Clone, Copy, Debug)]
pub(super) struct Nongravitational {
    pub(super) parameters: [f64; 3],
    pub(super) law: DistanceLaw,
}

impl Nongravitational {
    fn acceleration(self, state: State, sun: Source) -> Result<[f64; 3], NumericalError> {
        if !self.parameters.into_iter().all(f64::is_finite) {
            return Err(NumericalError::NonFinite);
        }
        if self.parameters == [0.0; 3] {
            return Ok([0.0; 3]);
        }
        let position = subtract(state.position, sun.point.position);
        let velocity = subtract(state.velocity, sun.velocity);
        let radius = distance(position, sun.point.body)?;
        let radial = position.map(|value| value / radius);
        let mut directions = [[0.0; 3]; 3];
        directions[0] = radial;
        if self.parameters[1] != 0.0 || self.parameters[2] != 0.0 {
            // r-hat x v has the same normal as r x v without multiplying
            // by a potentially large orbital radius. A zero normal has no
            // transverse/normal frame; do not invent an arbitrary one.
            let angular = cross(radial, velocity);
            let magnitude = length(angular)?;
            if magnitude == 0.0 {
                return Err(NumericalError::InvalidOrbit(
                    "undefined nongravitational orbital plane",
                ));
            }
            let normal = angular.map(|value| value / magnitude);
            directions[1] = cross(normal, radial);
            directions[2] = normal;
        }
        let factor = self.law.value(radius)?;
        let mut total = [Sum::default(); 3];
        for (direction, parameter) in directions.into_iter().zip(self.parameters) {
            for axis in 0..3 {
                total[axis].add(factor * parameter * direction[axis]);
            }
        }
        finite_vector(total.map(Sum::value))
    }
}

/// Fully explicit force selection. An absent family is a configuration choice,
/// never an automatic fallback after an evaluation error. The caller establishes
/// the actual defaults and retained data provenance once per ephemeris owner.
#[derive(Clone, Debug)]
pub(super) struct ForceModel {
    pub(super) relativity: Option<Eih>,
    pub(super) zonals: Vec<Zonals>,
    pub(super) nongravitational: Option<Nongravitational>,
}

impl ForceModel {
    /// Evaluate every selected force in one common frame. Only the identified
    /// self source is excluded; every other Newtonian source remains active.
    pub(super) fn acceleration(
        &self,
        state: State,
        sources: &[Source],
        exclude: Option<i32>,
    ) -> Result<[f64; 3], NumericalError> {
        state.checked()?;
        for source in sources {
            source.checked()?;
        }
        let newtonian = newtonian::acceleration_from(
            state.position,
            sources.iter().map(|source| source.point),
            exclude,
        )?;
        let mut result = newtonian.map(Sum::new);
        if let Some(model) = self.relativity {
            add(&mut result, model.acceleration(state, sources, exclude)?);
        }
        for zonals in &self.zonals {
            if Some(zonals.body) == exclude {
                continue;
            }
            let source = find(sources, zonals.body)?;
            add(&mut result, zonals.acceleration(state.position, source)?);
        }
        if let Some(model) = self.nongravitational {
            if !model.parameters.into_iter().all(f64::is_finite) {
                return Err(NumericalError::NonFinite);
            }
            // A1=A2=A3=0 requires neither a Sun nor an orbital-plane choice.
            if model.parameters != [0.0; 3] {
                add(&mut result, model.acceleration(state, find(sources, 10)?)?);
            }
        }
        finite_vector(result.map(Sum::value))
    }
}

impl Eih {
    fn acceleration(
        self,
        state: State,
        sources: &[Source],
        exclude: Option<i32>,
    ) -> Result<[f64; 3], NumericalError> {
        if !self.speed_of_light.is_finite() {
            return Err(NumericalError::NonFinite);
        }
        if self.speed_of_light <= 0.0 {
            return Err(NumericalError::InvalidSettings);
        }
        if matches!(self.outer, RelativisticSources::Sun)
            && find(sources, 10)?.post_newtonian_acceleration.is_none()
        {
            return Err(NumericalError::MissingSource { body: 10 });
        }
        let inverse_c_squared = (1.0 / self.speed_of_light) / self.speed_of_light;
        let selected = |source: &&Source| {
            source.post_newtonian_acceleration.is_some()
                && source.point.gravitational_parameter != 0.0
                && Some(source.point.body) != exclude
        };
        let mut potential = Sum::default();
        for source in sources.iter().filter(selected) {
            let radius = distance(
                subtract(state.position, source.point.position),
                source.point.body,
            )?;
            potential.add(source.point.gravitational_parameter / radius);
        }
        let target_potential = potential.value();
        let target_speed_squared = dot(state.velocity, state.velocity);
        let mut total = [Sum::default(); 3];
        for source in sources.iter().filter(selected) {
            if matches!(self.outer, RelativisticSources::Sun) && source.point.body != 10 {
                continue;
            }
            let displacement = subtract(state.position, source.point.position);
            let radius = distance(displacement, source.point.body)?;
            let unit = displacement.map(|value| value / radius);
            let mut source_potential = Sum::default();
            for other in sources.iter().filter(selected) {
                if other.point.body == source.point.body {
                    continue;
                }
                let separation = distance(
                    subtract(source.point.position, other.point.position),
                    other.point.body,
                )?;
                source_potential.add(other.point.gravitational_parameter / separation);
            }
            let source_acceleration = source
                .post_newtonian_acceleration
                .expect("selected physical acceleration");
            // ri-rj is displacement. The published radial force is in the
            // opposite direction, rj-ri. Keep those signs explicit here.
            let projection = dot(unit, source.velocity);
            let radial = -4.0 * target_potential - source_potential.value()
                + target_speed_squared
                + 2.0 * dot(source.velocity, source.velocity)
                - 4.0 * dot(state.velocity, source.velocity)
                - 1.5 * projection * projection
                - 0.5 * dot(displacement, source_acceleration);
            // Contract with the unit separation before scaling, avoiding an
            // unnecessary r^-3 intermediate in the velocity-dependent term.
            let velocity_factor = dot(
                unit,
                std::array::from_fn(|axis| {
                    4.0 * state.velocity[axis] - 3.0 * source.velocity[axis]
                }),
            );
            let relative_velocity = subtract(state.velocity, source.velocity);
            let monopole = (source.point.gravitational_parameter / radius) / radius;
            let velocity_scale = monopole * inverse_c_squared;
            let acceleration_scale =
                3.5 * (source.point.gravitational_parameter / radius) * inverse_c_squared;
            for axis in 0..3 {
                total[axis].add(-monopole * inverse_c_squared * radial * unit[axis]);
                total[axis].add(velocity_scale * velocity_factor * relative_velocity[axis]);
                total[axis].add(acceleration_scale * source_acceleration[axis]);
            }
        }
        finite_vector(total.map(Sum::value))
    }
}

fn find(sources: &[Source], body: i32) -> Result<Source, NumericalError> {
    sources
        .iter()
        .find(|source| source.point.body == body)
        .copied()
        .ok_or(NumericalError::MissingSource { body })
}
fn subtract(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|axis| first[axis] - second[axis])
}
fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    let mut value = Sum::default();
    for axis in 0..3 {
        value.add(first[axis] * second[axis]);
    }
    value.value()
}
fn cross(first: [f64; 3], second: [f64; 3]) -> [f64; 3] {
    [
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    ]
}
fn length(vector: [f64; 3]) -> Result<f64, NumericalError> {
    let result = vector[0].hypot(vector[1]).hypot(vector[2]);
    if result.is_finite() {
        Ok(result)
    } else {
        Err(NumericalError::NonFinite)
    }
}
fn distance(vector: [f64; 3], body: i32) -> Result<f64, NumericalError> {
    let result = length(vector)?;
    if result == 0.0 {
        Err(NumericalError::SingularForce { body })
    } else {
        Ok(result)
    }
}
fn finite_vector(value: [f64; 3]) -> Result<[f64; 3], NumericalError> {
    if value.into_iter().all(f64::is_finite) {
        Ok(value)
    } else {
        Err(NumericalError::NonFinite)
    }
}
fn add(total: &mut [Sum; 3], value: [f64; 3]) {
    for axis in 0..3 {
        total[axis].add(value[axis]);
    }
}

#[cfg(test)]
#[path = "forces_tests.rs"]
mod tests;
