// SPDX-License-Identifier: AGPL-3.0
//! Original conversion of finite elliptic/hyperbolic elements to Cartesian PV.
//!
//! Kepler's equations are M=E-e*sin(E) and M=e*sinh(H)-H. Their derivatives
//! are strictly positive on the admitted domains, permitting a safeguarded
//! Newton solve within a physical bracket. Cartesian coordinates follow the
//! elementary conic parameterizations and an active periapsis/node/inclination
//! rotation. No third-party orbit constructor or numerical table is an input.
//!
//! Elements use a caller-selected length unit, radians, and a compatible GM.
//! A finite semimajor axis cannot describe the parabolic e=1 case; reject it
//! explicitly rather than substituting a neighboring elliptic orbit.

use super::state::{NumericalError, State, Sum};

/// Osculating elements in one declared inertial plane, with angles in radians.
#[derive(Clone, Copy, Debug)]
pub(super) struct Elements {
    pub(super) semimajor_axis: f64,
    pub(super) eccentricity: f64,
    pub(super) inclination: f64,
    pub(super) ascending_node: f64,
    pub(super) periapsis_argument: f64,
    pub(super) mean_anomaly: f64,
}

/// Compute x-sin(x) or sinh(x)-x without losing its cubic leading term.
///
/// For |x|<=1 the elementary Taylor series converges rapidly. Termination is
/// based on binary64 resolution of the sum, with a fixed work bound, not on a
/// selected orbital input. Outside this interval direct subtraction is stable.
fn cubic_remainder(argument: f64, hyperbolic: bool) -> f64 {
    if argument.abs() > 1.0 {
        return if hyperbolic {
            argument.sinh() - argument
        } else {
            argument - argument.sin()
        };
    }
    let square = argument * argument;
    let mut term = argument * (square / 6.0);
    let mut sum = Sum::new(term);
    for order in 1..32 {
        let denominator = f64::from((2 * order + 2) * (2 * order + 3));
        term *= if hyperbolic {
            square / denominator
        } else {
            -square / denominator
        };
        let before = sum.value();
        sum.add(term);
        if sum.value() == before {
            break;
        }
    }
    sum.value()
}

/// Safeguarded monotone root solve; no unbounded iteration or solver fallback.
///
/// The callback returns a residual and its strictly positive derivative. A
/// floating-point fixed point is the numerical resolution limit, not a license
/// to return an unconverged iterate after exhausting the iteration bound.
fn monotone_root(
    mut low: f64,
    mut high: f64,
    initial: f64,
    equation: impl Fn(f64) -> (f64, f64),
) -> Result<f64, NumericalError> {
    let mut value = initial.clamp(low, high);
    for _ in 0..128 {
        let (residual, derivative) = equation(value);
        if !residual.is_finite() || !derivative.is_finite() || derivative <= 0.0 {
            return Err(NumericalError::NonFinite);
        }
        if residual == 0.0 {
            return Ok(value);
        }
        if residual < 0.0 {
            low = value;
        } else {
            high = value;
        }
        let proposed = value - residual / derivative;
        if proposed == value {
            return Ok(value);
        }
        let next = if low < proposed && proposed < high && proposed.is_finite() {
            proposed
        } else {
            low + (high - low) / 2.0
        };
        if next == low || next == high {
            return Ok(next);
        }
        value = next;
    }
    Err(NumericalError::IterationLimit)
}

/// Eccentric anomaly for an ellipse, preserving small signed input anomalies.
fn elliptic_anomaly(mean: f64, eccentricity: f64) -> Result<f64, NumericalError> {
    // Trigonometric argument reduction uses mathematical 2*pi rather than a
    // division by its rounded binary64 spelling. It also avoids adding pi to
    // a tiny signed mean anomaly before a remainder operation.
    let (sine_mean, cosine_mean) = mean.sin_cos();
    let reduced = sine_mean.atan2(cosine_mean);
    let magnitude = reduced.abs();
    if magnitude == 0.0 {
        return Ok(reduced);
    }
    let deficit = 1.0 - eccentricity;
    let high = (magnitude + eccentricity).min(magnitude / deficit);
    let low = (magnitude - eccentricity).max(0.0);
    let root = monotone_root(low, high, high, |anomaly| {
        let mut residual = Sum::new(deficit * anomaly);
        residual.add(eccentricity * cubic_remainder(anomaly, false));
        residual.add(-magnitude);
        let sine_half = (anomaly / 2.0).sin();
        (
            residual.value(),
            deficit + eccentricity * (2.0 * sine_half * sine_half),
        )
    })?;
    Ok(root.copysign(reduced))
}

/// Hyperbolic anomaly on its unwrapped time domain; unlike E it is not periodic.
fn hyperbolic_anomaly(mean: f64, eccentricity: f64) -> Result<f64, NumericalError> {
    let magnitude = mean.abs();
    if magnitude == 0.0 {
        return Ok(mean);
    }
    let excess = eccentricity - 1.0;
    let equation = |anomaly: f64| {
        let mut residual = Sum::new(excess * anomaly.sinh());
        residual.add(cubic_remainder(anomaly, true));
        residual.add(-magnitude);
        let sine_half = (anomaly / 2.0).sinh();
        (
            residual.value(),
            excess + eccentricity * (2.0 * sine_half * sine_half),
        )
    };
    // Linear and cubic lower bounds on e*sinh(H)-H give upper bounds on H.
    // The logarithmic bound keeps large mean anomalies away from huge H values.
    let mut high = ((magnitude / eccentricity).asinh() + 1.0)
        .min(magnitude / excess)
        .min(((magnitude / eccentricity) * 6.0).cbrt());
    for _ in 0..64 {
        let (residual, _) = equation(high);
        if !residual.is_finite() {
            return Err(NumericalError::NonFinite);
        }
        if residual >= 0.0 {
            return monotone_root(
                0.0,
                high,
                (magnitude / eccentricity).asinh().min(high),
                equation,
            )
            .map(|root| root.copysign(mean));
        }
        high *= 2.0;
        if !high.is_finite() {
            return Err(NumericalError::NonFinite);
        }
    }
    Err(NumericalError::IterationLimit)
}

/// Position/velocity from a finite nonparabolic conic, without propagation.
///
/// Ellipses require a>0 and 0<=e<1; hyperbolas require a<0 and e>1. All
/// angles may be any finite radians. GM must be finite and positive. Invalid
/// domains, nonfinite derived states and unresolved roots return explicit errors.
pub(super) fn initial_state(elements: Elements, gm: f64) -> Result<State, NumericalError> {
    if ![
        elements.semimajor_axis,
        elements.eccentricity,
        elements.inclination,
        elements.ascending_node,
        elements.periapsis_argument,
        elements.mean_anomaly,
        gm,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(NumericalError::NonFinite);
    }
    if gm <= 0.0 {
        return Err(NumericalError::InvalidOrbit("GM must be positive"));
    }
    let a = elements.semimajor_axis;
    let e = elements.eccentricity;
    let speed = gm.sqrt() / a.abs().sqrt();
    let (x, y, vx, vy) = if a > 0.0 && (0.0..1.0).contains(&e) {
        let anomaly = elliptic_anomaly(elements.mean_anomaly, e)?;
        let (sine, cosine) = anomaly.sin_cos();
        let sine_half = (anomaly / 2.0).sin();
        let deficit = 1.0 - e;
        let shape = (deficit * (1.0 + e)).sqrt();
        let radius_over_a = deficit + e * (2.0 * sine_half * sine_half);
        (
            a * (deficit - 2.0 * sine_half * sine_half),
            a * shape * sine,
            -speed * sine / radius_over_a,
            speed * shape * cosine / radius_over_a,
        )
    } else if a < 0.0 && e > 1.0 {
        let anomaly = hyperbolic_anomaly(elements.mean_anomaly, e)?;
        let sine = anomaly.sinh();
        let cosine = anomaly.cosh();
        let sine_half = (anomaly / 2.0).sinh();
        let excess = e - 1.0;
        let shape = excess.sqrt() * (e + 1.0).sqrt();
        let radius_over_a = excess + e * (2.0 * sine_half * sine_half);
        let length = -a;
        (
            length * (excess - 2.0 * sine_half * sine_half),
            length * shape * sine,
            -speed * sine / radius_over_a,
            speed * shape * cosine / radius_over_a,
        )
    } else {
        return Err(NumericalError::InvalidOrbit(
            "inconsistent semimajor axis and eccentricity",
        ));
    };
    let (sin_node, cos_node) = elements.ascending_node.sin_cos();
    let (sin_periapsis, cos_periapsis) = elements.periapsis_argument.sin_cos();
    let (sin_inclination, cos_inclination) = elements.inclination.sin_cos();
    let periapsis = [
        cos_periapsis * cos_node - sin_periapsis * sin_node * cos_inclination,
        cos_periapsis * sin_node + sin_periapsis * cos_node * cos_inclination,
        sin_periapsis * sin_inclination,
    ];
    let transverse = [
        -sin_periapsis * cos_node - cos_periapsis * sin_node * cos_inclination,
        -sin_periapsis * sin_node + cos_periapsis * cos_node * cos_inclination,
        cos_periapsis * sin_inclination,
    ];
    let state = State {
        position: std::array::from_fn(|axis| x * periapsis[axis] + y * transverse[axis]),
        velocity: std::array::from_fn(|axis| vx * periapsis[axis] + vy * transverse[axis]),
    }
    .checked()?;
    if state.position.into_iter().all(|value| value == 0.0) {
        return Err(NumericalError::InvalidOrbit(
            "periapsis is not representable",
        ));
    }
    Ok(state)
}

#[cfg(test)]
#[path = "conic_tests.rs"]
mod tests;
