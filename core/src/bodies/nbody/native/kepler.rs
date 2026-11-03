// SPDX-License-Identifier: AGPL-3.0
//! Original compensated Cartesian Kepler flow, in caller-selected units/frame.
//!
//! Sundman time gives G0'=-b G1 and Gn'=G(n-1), with b=2-r0|v0|²/GM.
//! The normalized Kepler equation is G1+eta G2+G3=tau and its derivative is
//! radius/r0. Power series and exact doubling identities evaluate those same
//! functions through zero energy; no element/anomaly retiming is involved.
//! See Wisdom & Hernandez (2015), the universal initial-value formulation:
//! <https://arxiv.org/html/1508.02699v1>. No program is an implementation input.
//!
//! The root, parameters and f/g products retain two binary64 words until the
//! final state is rounded. Dropping their residuals earlier loses near-zero
//! components or periapsis-return timing despite an apparently small residual.
//! The domain is a nonradial orbit. Radial collision/ejection continuation,
//! nonfinite/unrepresentable arithmetic and unfinished solves are explicit
//! errors, not alternate analytic/backend dispatches.

use super::state::{NumericalError, State};
#[path = "kepler_words.rs"]
mod words;
use words::Words;

#[derive(Clone, Copy)]
struct Functions {
    zero: Words,
    one: Words,
    two: Words,
    three: Words,
}

/// Propagate the actual supplied Cartesian initial-value problem. All fields
/// use one compensated route; only the final exposed state rounds to binary64.
pub(super) fn advance(initial: State, gm: f64, interval: f64) -> Result<State, NumericalError> {
    initial.checked()?;
    if !gm.is_finite() || !interval.is_finite() {
        return Err(NumericalError::NonFinite);
    }
    if gm <= 0.0 {
        return Err(NumericalError::InvalidMass { body: 10 });
    }
    if interval == 0.0 {
        return Ok(initial);
    }
    let position = initial.position.map(Words::new);
    let original_velocity = initial.velocity.map(Words::new);
    let radius = vector_norm(position)?;
    if radius.sign() == 0.0 {
        return Err(NumericalError::SingularForce { body: 10 });
    }
    let speed = Words::new(gm).square_root().divided(radius.square_root());
    if !speed.finite() || speed.sign() <= 0.0 {
        return Err(NumericalError::NonFinite);
    }
    let direction = interval.signum();
    let radial = position.map(|value| value.divided(radius));
    let velocity = original_velocity.map(|value| value.scaled(direction).divided(speed));
    let angular = [
        radial[1]
            .times(velocity[2])
            .minus(radial[2].times(velocity[1])),
        radial[2]
            .times(velocity[0])
            .minus(radial[0].times(velocity[2])),
        radial[0]
            .times(velocity[1])
            .minus(radial[1].times(velocity[0])),
    ];
    if vector_norm(angular)?.sign() == 0.0 {
        return Err(NumericalError::InvalidOrbit(
            "radial Kepler flow is outside the admitted domain",
        ));
    }
    let eta = product_sum(radial, velocity);
    let b = normalized_energy(initial, gm, radius)?;
    let elapsed = Words::new(interval.abs());
    let tau = elapsed.divided(radius).times(speed);
    if !eta.finite() || !tau.finite() {
        return Err(NumericalError::NonFinite);
    }
    let values = if tau.sign() == 0.0 {
        universal(Words::new(0.0), b)?
    } else {
        solve(tau, b, eta)?
    };
    // f=G0+(b-1)G2 and g=r0(G1+eta G2)/v_c. Retaining both
    // words preserves the cancellation through the Cartesian products too.
    let f = values.zero.plus(b.minus(Words::new(1.0)).times(values.two));
    let g = if tau.sign() == 0.0 {
        // Removable scaled-time underflow: gravity has no represented
        // increment, but v0*dt can still be a nonzero Cartesian displacement.
        elapsed
    } else {
        radius
            .times(values.one.plus(eta.times(values.two)))
            .divided(speed)
    };
    let endpoint = std::array::from_fn(|axis| {
        f.times(position[axis])
            .plus(g.scaled(direction).times(original_velocity[axis]))
    });
    let final_radius = vector_norm(endpoint)?;
    if final_radius.sign() == 0.0 {
        return Err(NumericalError::SingularForce { body: 10 });
    }
    let ratio = final_radius.divided(radius);
    if !ratio.finite() || ratio.sign() <= 0.0 {
        return Err(NumericalError::NonFinite);
    }
    let gd = values.zero.plus(eta.times(values.one)).divided(ratio);
    let radial_speed = speed.times(values.one).divided(ratio).scaled(-direction);
    let final_velocity = std::array::from_fn(|axis| {
        radial_speed
            .times(radial[axis])
            .plus(gd.times(original_velocity[axis]))
            .value()
    });
    State {
        position: endpoint.map(Words::value),
        velocity: final_velocity,
    }
    .checked()
}

fn solve(tau: Words, b: Words, eta: Words) -> Result<Functions, NumericalError> {
    let mut lower = Words::new(0.0);
    let mut upper = if tau.value() < 1.0 {
        tau
    } else {
        Words::new(1.0)
    };
    if upper.sign() == 0.0 {
        return Err(NumericalError::StepUnderflow);
    }
    let mut bracketed = false;
    for _ in 0..1024 {
        let (value, _, values) = equation(upper, b, eta, tau)?;
        if value.sign() == 0.0 {
            return Ok(values);
        }
        if value.sign() > 0.0 {
            bracketed = true;
            break;
        }
        lower = upper;
        upper = upper.scaled(2.0);
        if !upper.finite() {
            return Err(NumericalError::NonFinite);
        }
    }
    if !bracketed {
        return Err(NumericalError::IterationLimit);
    }
    let mut root = lower.plus(upper.minus(lower).scaled(0.5));
    for _ in 0..512 {
        let (value, derivative, values) = equation(root, b, eta, tau)?;
        if value.sign() == 0.0 {
            return Ok(values);
        }
        if value.sign() < 0.0 {
            lower = root;
        } else {
            upper = root;
        }
        let newton = root.minus(value.divided(derivative));
        let next = if newton.finite()
            && newton.minus(lower).sign() > 0.0
            && upper.minus(newton).sign() > 0.0
        {
            newton
        } else {
            lower.plus(upper.minus(lower).scaled(0.5))
        };
        if next == root || next == lower || next == upper {
            let a = equation(lower, b, eta, tau)?;
            let c = equation(upper, b, eta, tau)?;
            return Ok(if a.0.magnitude() < c.0.magnitude() {
                a.2
            } else {
                c.2
            });
        }
        root = next;
    }
    Err(NumericalError::IterationLimit)
}

fn equation(
    x: Words,
    b: Words,
    eta: Words,
    tau: Words,
) -> Result<(Words, Words, Functions), NumericalError> {
    let values = universal(x, b)?;
    let residual = values
        .one
        .plus(eta.times(values.two))
        .plus(values.three)
        .minus(tau);
    let derivative = values.zero.plus(eta.times(values.one)).plus(values.two);
    if !residual.finite() || !derivative.finite() {
        return Err(NumericalError::NonFinite);
    }
    if derivative.sign() < 0.0 {
        return Err(NumericalError::InvalidOrbit("negative Kepler radius"));
    }
    Ok((residual, derivative, values))
}

fn universal(x: Words, b: Words) -> Result<Functions, NumericalError> {
    // Halving limits |b*x²| to1/2, so the elementary defining series
    // converge rapidly. Doubling then restores the same universal functions,
    // without a rounded trigonometric phase or a fitted circular period.
    let mut reduced = x;
    let mut angle = b.magnitude().sqrt() * x.magnitude();
    if !angle.is_finite() {
        return Err(NumericalError::NonFinite);
    }
    let mut doublings = 0;
    while angle > 0.5_f64.sqrt() {
        reduced = reduced.scaled(0.5);
        angle *= 0.5;
        doublings += 1;
    }
    let z = b.times(reduced).times(reduced).negative();
    let mut terms = [
        Words::new(1.0),
        Words::new(1.0),
        Words::new(0.5),
        Words::new(1.0).divided(Words::new(6.0)),
    ];
    let mut sums = terms;
    let mut converged = false;
    for k in 1..=64 {
        for i in 0..4 {
            let degree = 2 * k + i;
            terms[i] = terms[i]
                .times(z)
                .divided(Words::new(((degree - 1) * degree) as f64));
            sums[i] = sums[i].plus(terms[i]);
        }
        if (0..4).all(|i| {
            terms[i].magnitude() <= (f64::EPSILON * f64::EPSILON / 16.0) * sums[i].magnitude()
        }) {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(NumericalError::IterationLimit);
    }
    let mut values = Functions {
        zero: sums[0],
        one: sums[1].times(reduced),
        two: sums[2].times(reduced).times(reduced),
        three: sums[3].times(reduced).times(reduced).times(reduced),
    };
    for _ in 0..doublings {
        let square = values.one.times(values.one);
        values = Functions {
            zero: values.zero.times(values.zero).minus(b.times(square)),
            one: values.zero.times(values.one).scaled(2.0),
            two: square.scaled(2.0),
            three: values.three.plus(values.one.times(values.two)).scaled(2.0),
        };
    }
    if [values.zero, values.one, values.two, values.three]
        .into_iter()
        .all(Words::finite)
    {
        Ok(values)
    } else {
        Err(NumericalError::NonFinite)
    }
}

/// b=2-r0|v0|²/GM, retaining products/divisions through binding-energy
/// cancellation rather than rounding the normalized speed before subtracting.
fn normalized_energy(initial: State, gm: f64, radius: Words) -> Result<Words, NumericalError> {
    let mut energy = Words::new(2.0);
    for value in initial.velocity.map(Words::new) {
        let kinetic = radius.times(value).divided(Words::new(gm)).times(value);
        energy = energy.minus(kinetic);
    }
    if energy.finite() {
        Ok(energy)
    } else {
        Err(NumericalError::NonFinite)
    }
}

fn product_sum(first: [Words; 3], second: [Words; 3]) -> Words {
    first
        .into_iter()
        .zip(second)
        .fold(Words::new(0.0), |sum, (a, b)| sum.plus(a.times(b)))
}

fn vector_norm(vector: [Words; 3]) -> Result<Words, NumericalError> {
    let estimate = norm(vector.map(Words::value))?;
    if estimate == 0.0 {
        return Ok(Words::new(0.0));
    }
    // Binary exponent scaling is exact and bounds squared coordinates. It
    // avoids overflowing an otherwise representable physical distance.
    let exponent = estimate.to_bits() & 0x7ff0_0000_0000_0000;
    let scale = if exponent == 0 {
        f64::MIN_POSITIVE
    } else {
        f64::from_bits(exponent)
    };
    let scaled = vector.map(|value| value.divided(Words::new(scale)));
    let result = product_sum(scaled, scaled).square_root().scaled(scale);
    if result.finite() {
        Ok(result)
    } else {
        Err(NumericalError::NonFinite)
    }
}

fn norm(vector: [f64; 3]) -> Result<f64, NumericalError> {
    let value = vector[0].hypot(vector[1]).hypot(vector[2]);
    if value.is_finite() {
        Ok(value)
    } else {
        Err(NumericalError::NonFinite)
    }
}
#[cfg(test)]
fn dot(first: [f64; 3], second: [f64; 3]) -> f64 {
    product_sum(first.map(Words::new), second.map(Words::new)).value()
}

#[cfg(test)]
#[path = "kepler_tests.rs"]
mod tests;
