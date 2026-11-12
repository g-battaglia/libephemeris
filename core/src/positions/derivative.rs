//! Adaptive fourth-order derivatives of a continuous apparent place.
//!
//! Exact dyadic day offsets keep the geometry, observer and frame clocks on
//! the same TT spacing. Richardson estimates select the coarsest resolved
//! derivative before binary64 evaluation noise grows. Longitude is periodic;
//! latitude is linear. This is a differentiation method, not an acceptance
//! tolerance or an estimate of ephemeris-model accuracy. Exhaustion returns
//! the best available estimate, without a universal accuracy guarantee.

use crate::errors::{Error, Result};

/// Differentiate longitude/latitude in degrees per TT day.
///
/// `sample(offset)` evaluates the same local physical model at a TT-day
/// displacement. `support` is its continuous available interval relative to
/// the central epoch and must contain zero. The initial positive dyadic step
/// is reduced to fit support. A fourth-order one-sided stencil is used when
/// a centered stencil would cross an owned polynomial boundary. Every actual
/// sample/provider error propagates. Undefined support or a nonrepresentable
/// initial stencil is [`Error::MathDomain`]. At most eight refinements are
/// attempted. A third scale confirms early resolution; exhaustion, a refused
/// refinement or an underflowed step retains the best available estimate.
pub(crate) fn angular_rates(
    sample: &dyn Fn(f64) -> Result<[f64; 2]>,
    initial_half_step: f64,
    support: (f64, f64),
) -> Result<[f64; 2]> {
    place_rates(sample, initial_half_step, support, true, true)
}

/// Angular derivative restricted by actual sample-level record ownership.
///
/// `None` rejects a TT offset solely because a consumed physical provider read
/// would leave its central record. Such a step is halved before differentiating;
/// all actual calculation errors remain errors. This additionally handles
/// nonlinear retarded-time mappings without treating a coverage envelope or a
/// linearized emission clock as a continuous polynomial domain.
pub(crate) fn angular_rates_on_owned_domain(
    sample: &dyn Fn(f64) -> Result<Option<[f64; 2]>>,
    initial_half_step: f64,
    support: (f64, f64),
) -> Result<[f64; 2]> {
    admitted_place_rates(sample, initial_half_step, support, true, true)
}

/// Angular velocity from actual Cartesian places in the central tangent plane.
///
/// The fixed central longitude/latitude directions are the spherical Jacobian.
/// Projecting AU displacements before differencing avoids subtraction of
/// rounded angles on short owned intervals. Every sample remains the actual
/// physical output ray; no trajectory or record is extrapolated. `None` has
/// the same owned-record admission meaning as [`angular_rates_on_owned_domain`].
/// Rates are degrees per TT day. The established spherical convention remains
/// in effect for a pole or origin where the longitude Jacobian is undefined.
#[cfg(test)]
pub(crate) fn cartesian_angular_rates_on_owned_domain(
    sample: &dyn Fn(f64) -> Result<Option<[f64; 3]>>,
    initial_half_step: f64,
    support: (f64, f64),
) -> Result<[f64; 2]> {
    let central = sample(0.0)?.ok_or(Error::MathDomain { value: 0.0 })?;
    crate::positions::check_finite_vec3(central)?;
    let horizontal = central[0].hypot(central[1]);
    let radius = horizontal.hypot(central[2]);
    if horizontal == 0.0 || radius == 0.0 {
        return angular_rates_on_owned_domain(
            &|offset| {
                sample(offset)?
                    .map(|value| {
                        crate::positions::cartesian_to_spherical_deg(value)
                            .map(|angles| [angles[0], angles[1]])
                    })
                    .transpose()
            },
            initial_half_step,
            support,
        );
    }
    let longitude = [-central[1] / horizontal, central[0] / horizontal, 0.0];
    let latitude = [
        -(central[0] / horizontal) * (central[2] / radius),
        -(central[1] / horizontal) * (central[2] / radius),
        horizontal / radius,
    ];
    admitted_place_rates(
        &|offset| {
            let Some(value) = sample(offset)? else {
                return Ok(None);
            };
            crate::positions::check_finite_vec3(value)?;
            let displacement = std::array::from_fn(|i| value[i] - central[i]);
            Ok(Some([
                (crate::positions::compensated_dot3(longitude, displacement) / horizontal)
                    .to_degrees(),
                (crate::positions::compensated_dot3(latitude, displacement) / radius).to_degrees(),
            ]))
        },
        initial_half_step,
        support,
        false,
        true,
    )
}

/// Differentiate actual Cartesian places in AU per TT day on owned support.
/// The same fourth-order domain/refinement rule as angular rates applies;
/// Cartesian components use relative resolution and never wrap or receive an
/// angular absolute floor. Every sample and domain failure propagates.
pub(crate) fn cartesian_rates(
    sample: &dyn Fn(f64) -> Result<[f64; 3]>,
    initial_half_step: f64,
    support: (f64, f64),
) -> Result<[f64; 3]> {
    place_rates(sample, initial_half_step, support, false, false)
}

/// Differentiate small Cartesian corrections on their actual owned domains.
/// Units are input units per TT day. `None` rejects a side sample solely on
/// provider/clock ownership, with the same admission and genuine-error rule as
/// [`angular_rates_on_owned_domain`]. No component is wrapped or angular-scaled.
pub(crate) fn cartesian_rates_on_owned_domain(
    sample: &dyn Fn(f64) -> Result<Option<[f64; 3]>>,
    initial_half_step: f64,
    support: (f64, f64),
) -> Result<[f64; 3]> {
    admitted_place_rates(sample, initial_half_step, support, false, false)
}

/// Exact spherical angular Jacobian of a Cartesian ray and its TT velocity.
/// Inputs use a common length unit and that unit per TT day; rates are degrees
/// per TT day. Compensated tangent products retain nearly stationary residuals.
/// The established origin/pole convention returns zero angular rates where
/// longitude is undefined. Nonfinite vectors are rejected.
pub(crate) fn spherical_angular_velocity(
    position: [f64; 3],
    velocity: [f64; 3],
) -> Result<[f64; 2]> {
    crate::positions::check_finite_vec3(position)?;
    crate::positions::check_finite_vec3(velocity)?;
    let horizontal = position[0].hypot(position[1]);
    let radius = horizontal.hypot(position[2]);
    if horizontal == 0.0 || radius == 0.0 {
        return Ok([0.0; 2]);
    }
    let longitude = [-position[1] / horizontal, position[0] / horizontal, 0.0];
    let latitude = [
        -(position[0] / horizontal) * (position[2] / radius),
        -(position[1] / horizontal) * (position[2] / radius),
        horizontal / radius,
    ];
    Ok([
        (crate::positions::compensated_dot3(longitude, velocity) / horizontal).to_degrees(),
        (crate::positions::compensated_dot3(latitude, velocity) / radius).to_degrees(),
    ])
}

/// Analytic-path Jacobian with one compensated determinant per numerator.
pub(crate) fn spherical_angular_velocity_direct(
    position: [f64; 3],
    velocity: [f64; 3],
) -> Result<[f64; 2]> {
    crate::positions::check_finite_vec3(position)?;
    crate::positions::check_finite_vec3(velocity)?;
    let horizontal = position[0].hypot(position[1]);
    let radius = horizontal.hypot(position[2]);
    if horizontal == 0.0 || radius == 0.0 {
        return Ok([0.0; 2]);
    }
    // Exact power-of-two scaling keeps quadratic/cubic Jacobian products
    // representable for finite near-origin and large Cartesian inputs.
    let largest = position.into_iter().map(f64::abs).fold(0.0, f64::max);
    let exponent = ((largest.to_bits() >> 52) & 0x7ff).max(1);
    let scale = f64::from_bits(exponent << 52);
    let [x, y, z] = position.map(|value| AngularWord::from(value / scale));
    let [vx, vy, vz] = velocity.map(|value| AngularWord::from(value / scale));
    let horizontal_squared = x * x + y * y;
    let radius_squared = horizontal_squared + z * z;
    let longitude = (x * vy - y * vx) / horizontal_squared;
    let latitude = (horizontal_squared * vz - z * (x * vx + y * vy))
        / (radius_squared * horizontal_squared.sqrt());
    // The two words are the nearest binary64 decomposition of mathematical
    // 180/pi, computed independently to 100 decimal places. They are unit
    // conversion constants, independent of ephemeris inputs or measurements.
    let degrees_per_radian = AngularWord {
        hi: 57.29577951308232,
        lo: -1.9878495670576283e-15,
    };
    Ok([
        (longitude * degrees_per_radian).rounded(),
        (latitude * degrees_per_radian).rounded(),
    ])
}

/// Two binary64 words retaining product and sum residuals in the angular
/// Jacobian only. Physical scalar place arithmetic and models stay unchanged.
#[derive(Clone, Copy)]
pub(crate) struct AngularWord {
    pub(crate) hi: f64,
    pub(crate) lo: f64,
}
impl From<f64> for AngularWord {
    fn from(hi: f64) -> Self {
        Self { hi, lo: 0.0 }
    }
}
impl AngularWord {
    fn normalize(a: f64, b: f64) -> Self {
        let hi = a + b;
        let bp = hi - a;
        let lo = (a - (hi - bp)) + (b - bp);
        Self { hi, lo }
    }
    pub(crate) fn rounded(self) -> f64 {
        self.hi + self.lo
    }
    pub(crate) fn sqrt(self) -> Self {
        let q = self.hi.sqrt();
        let error = self - Self::from(q) * Self::from(q);
        Self::normalize(q, error.rounded() / (2.0 * q))
    }
}
impl std::ops::Add for AngularWord {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let sum = Self::normalize(self.hi, other.hi);
        Self::normalize(sum.hi, sum.lo + self.lo + other.lo)
    }
}
impl std::ops::Sub for AngularWord {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + Self {
            hi: -other.hi,
            lo: -other.lo,
        }
    }
}
impl std::ops::Mul for AngularWord {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let hi = self.hi * other.hi;
        let lo = self.hi.mul_add(other.hi, -hi)
            + self.hi * other.lo
            + self.lo * other.hi
            + self.lo * other.lo;
        Self::normalize(hi, lo)
    }
}
impl std::ops::Div for AngularWord {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let q = self.hi / other.hi;
        let residual = self - other * Self::from(q);
        Self::normalize(q, residual.rounded() / other.hi)
    }
}

fn place_rates<const N: usize>(
    sample: &dyn Fn(f64) -> Result<[f64; N]>,
    initial_half_step: f64,
    support: (f64, f64),
    wrap_longitude: bool,
    angular_resolution: bool,
) -> Result<[f64; N]> {
    admitted_place_rates(
        &|offset| sample(offset).map(Some),
        initial_half_step,
        support,
        wrap_longitude,
        angular_resolution,
    )
}

fn admitted_place_rates<const N: usize>(
    sample: &dyn Fn(f64) -> Result<Option<[f64; N]>>,
    initial_half_step: f64,
    support: (f64, f64),
    wrap_longitude: bool,
    angular_resolution: bool,
) -> Result<[f64; N]> {
    // One local physical model: equal offset bits imply the same place and
    // record admission. Cache refusals too, but propagate every actual error.
    // Admission can halve arbitrarily often; saturation only stops retention,
    // never refinement. No borrow survives the call into the sampler.
    let samples = std::cell::RefCell::new([None::<(u64, Option<[f64; N]>)>; 64]);
    let cached = |offset: f64| {
        let key = offset.to_bits();
        if let Some((_, value)) = samples
            .borrow()
            .iter()
            .flatten()
            .find(|(bits, _)| *bits == key)
        {
            return Ok(*value);
        }
        let value = sample(offset)?;
        if let Some(slot) = samples.borrow_mut().iter_mut().find(|slot| slot.is_none()) {
            *slot = Some((key, value));
        }
        Ok(value)
    };
    admitted_place_rates_impl(
        &cached,
        initial_half_step,
        support,
        wrap_longitude,
        angular_resolution,
    )
}

/// Fourth-order stencil and refinement, independent of sample retention.
fn admitted_place_rates_impl<const N: usize>(
    sample: &dyn Fn(f64) -> Result<Option<[f64; N]>>,
    initial_half_step: f64,
    support: (f64, f64),
    wrap_longitude: bool,
    angular_resolution: bool,
) -> Result<[f64; N]> {
    if !initial_half_step.is_finite()
        || initial_half_step <= 0.0
        || support.0.is_nan()
        || support.1.is_nan()
        || support.0 > 0.0
        || support.1 < 0.0
        || support.0 == support.1
    {
        return Err(Error::MathDomain {
            value: initial_half_step,
        });
    }
    let direction = if support.0 <= -2.0 * initial_half_step && support.1 >= 2.0 * initial_half_step
    {
        0.0
    } else if support.1 >= -support.0 {
        1.0
    } else {
        -1.0
    };
    let mut h = initial_half_step;
    let fits = |step: f64| {
        if direction == 0.0 {
            support.0 <= -2.0 * step && 2.0 * step <= support.1
        } else {
            let end = direction * 4.0 * step;
            support.0 <= end && end <= support.1
        }
    };
    while !fits(h) {
        h *= 0.5;
        if h == 0.0 {
            return Err(Error::MathDomain { value: h });
        }
    }
    let derivative = |step: f64| -> Result<Option<[f64; N]>> {
        if direction == 0.0 {
            let Some(minus) = sample(-step)? else {
                return Ok(None);
            };
            let Some(plus) = sample(step)? else {
                return Ok(None);
            };
            let Some(outer_minus) = sample(-2.0 * step)? else {
                return Ok(None);
            };
            let Some(outer_plus) = sample(2.0 * step)? else {
                return Ok(None);
            };
            let inner = [minus, plus];
            let outer = [outer_minus, outer_plus];
            Ok(Some(std::array::from_fn(|component| {
                let delta = |pair: [[f64; N]; 2]| {
                    if wrap_longitude && component == 0 {
                        short_arc(pair[1][component], pair[0][component])
                    } else {
                        pair[1][component] - pair[0][component]
                    }
                };
                (8.0 * delta(inner) - delta(outer)) / (12.0 * step)
            })))
        } else {
            let Some(center) = sample(0.0)? else {
                return Ok(None);
            };
            let mut sum = [0.0; N];
            for (offset, weight) in [(1.0, 48.0), (2.0, -36.0), (3.0, 16.0), (4.0, -3.0)] {
                let Some(value) = sample(direction * offset * step)? else {
                    return Ok(None);
                };
                for component in 0..N {
                    let delta = if wrap_longitude && component == 0 {
                        short_arc(value[component], center[component])
                    } else {
                        value[component] - center[component]
                    };
                    sum[component] += weight * delta;
                }
            }
            Ok(Some(sum.map(|value| value / (12.0 * direction * step))))
        }
    };
    let mut coarse = loop {
        if let Some(value) = derivative(h)? {
            break value;
        }
        h *= 0.5;
        if h == 0.0 {
            return Err(Error::MathDomain { value: h });
        }
    };
    let mut best = coarse;
    let mut best_score = f64::INFINITY;
    let mut provisional_resolution = false;
    for _ in 0..8 {
        let half = h / 2.0;
        if half == 0.0 {
            break;
        }
        let Some(fine) = derivative(half)? else { break };
        let mut score = 0.0_f64;
        for component in 0..N {
            if !coarse[component].is_finite() || !fine[component].is_finite() {
                return Err(Error::NonFiniteArgument {
                    value: fine[component],
                });
            }
            // The numerical resolution goal is two orders finer than PLAN's
            // relative angular-speed criterion, including its near-zero domain.
            let goal = if angular_resolution && fine[component].abs() <= 1e-6 {
                1e-8
            } else {
                fine[component].abs().max(f64::MIN_POSITIVE) * 1e-8
            };
            let estimate = (coarse[component] - fine[component]).abs() * (16.0 / 15.0);
            score = score.max(estimate / goal);
        }
        if score <= 1.0 {
            if provisional_resolution {
                // Two consecutive refinement pairs agree: retain the
                // coarsest resolved derivative before evaluation noise grows.
                return Ok(best);
            }
            // A single agreement can cancel competing truncation terms,
            // especially in a one-sided stencil. Confirm at a third scale.
            best = coarse;
            best_score = score;
            provisional_resolution = true;
        } else {
            if provisional_resolution {
                // The next scale contradicted the provisional agreement.
                // Discard that candidate rather than treating it as noise.
                best_score = f64::INFINITY;
                provisional_resolution = false;
            }
            if score < best_score {
                best = coarse;
                best_score = score;
            } else {
                break;
            }
        }
        h = half;
        coarse = fine;
    }
    // Confirmation controls early termination, not the valid-input domain.
    // Preserve the available estimate on exhaustion, refusal or underflow,
    // without claiming that an unconfirmed pair proves its accuracy.
    Ok(best)
}

fn short_arc(next: f64, previous: f64) -> f64 {
    let delta = next - previous;
    if delta > 180.0 {
        delta - 360.0
    } else if delta < -180.0 {
        delta + 360.0
    } else {
        delta
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn compensated_angular_jacobian_retains_near_stationary_products() {
        // Direct symbolic Cartesian Jacobians evaluated independently to
        // 100 decimal places, including exact binary64 input conversion.
        let cases: [([f64; 3], [f64; 3], [f64; 2]); 3] = [
            (
                [1.0, 1.0, 0.25],
                [1.0, 1.0 + 2.0_f64.powi(-30), -0.125],
                [2.668042644536231e-08, -14.732448814417108],
            ),
            (
                [1e-10, -2e-10, 3e-11],
                [-3e-7, -2e-7, 1e-7],
                [-91673.2472209317, 23660.161561136807],
            ),
            (
                [0.0005, 0.001, -0.0001],
                [0.02, 0.04 + 2.0_f64.powi(-40), -0.004],
                [2.08440831604393e-08, 3.6991155815798067e-09],
            ),
        ];
        for (position, velocity, expected) in cases {
            let actual = spherical_angular_velocity_direct(position, velocity).unwrap();
            for (actual, expected) in actual.into_iter().zip(expected) {
                assert!(
                    actual.to_bits().abs_diff(expected.to_bits()) <= 1,
                    "{actual} != {expected}"
                );
            }
        }
        for scale in [1e-300, 1e300] {
            let actual = spherical_angular_velocity_direct(
                [scale, scale, 0.25 * scale],
                [scale, 0.5 * scale, -0.125 * scale],
            )
            .unwrap();
            assert!(actual.into_iter().all(f64::is_finite));
        }
    }

    #[test]
    fn exact_sample_reuse_preserves_values_refusals_and_saturation() {
        use std::cell::{Cell, RefCell};
        for support in [(-1.0, 1.0), (0.0, 1.0), (-1.0, 0.0), (-0.0001, 0.0002)] {
            for admission in [1.0, 2.0_f64.powi(-80)] {
                let model = |t: f64| {
                    Ok((t.abs() <= admission)
                        .then_some([(359.99 + t + t.powi(5)).rem_euclid(360.0), t.sin()]))
                };
                let uncached = Cell::new(0);
                let expected = admitted_place_rates_impl(
                    &|t| {
                        uncached.set(uncached.get() + 1);
                        model(t)
                    },
                    1.0 / 64.0,
                    support,
                    true,
                    true,
                )
                .unwrap();
                let seen = RefCell::new(Vec::new());
                let actual = admitted_place_rates(
                    &|t| {
                        seen.borrow_mut().push(t.to_bits());
                        model(t)
                    },
                    1.0 / 64.0,
                    support,
                    true,
                    true,
                )
                .unwrap();
                assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
                assert!(seen.borrow().len() <= uncached.get());
                if admission == 1.0 {
                    let mut unique = seen.borrow().clone();
                    unique.sort_unstable();
                    unique.dedup();
                    assert_eq!(unique.len(), seen.borrow().len());
                    assert!(seen.borrow().len() < uncached.get());
                } else {
                    assert!(seen.borrow().len() > 64, "exercise a saturated cache");
                }
            }
            let model = |t: f64| Ok(Some([400.0 + t.sin(), t.powi(5), t.cos()]));
            let expected =
                admitted_place_rates_impl(&model, 1.0 / 64.0, support, false, false).unwrap();
            let actual = cartesian_rates_on_owned_domain(&model, 1.0 / 64.0, support).unwrap();
            assert_eq!(actual.map(f64::to_bits), expected.map(f64::to_bits));
        }
    }

    #[test]
    fn differentiates_wrapped_cubic_places_on_centered_and_owned_boundary_domains() {
        let sample = |t: f64| {
            Ok([
                (359.99 + 2.0 * t + t * t * t).rem_euclid(360.0),
                0.001 + 0.000002 * t - t * t * t,
            ])
        };
        for support in [(-1.0, 1.0), (0.0, 1.0), (-1.0, 0.0)] {
            let rate = angular_rates(&sample, 1.0 / 64.0, support).unwrap();
            assert!((rate[0] - 2.0).abs() < 1e-10);
            assert!((rate[1] - 0.000002).abs() < 1e-16);
        }
    }

    #[test]
    fn cartesian_cubic_derivatives_use_linear_components_on_every_owned_domain() {
        let sample = |t: f64| {
            Ok([
                400.0 + 3.0 * t + t * t * t,
                -500.0 - 2.0 * t - t * t,
                0.001 + 0.000002 * t + t * t * t,
            ])
        };
        for support in [(-1.0, 1.0), (0.0, 1.0), (-1.0, 0.0)] {
            let rate = cartesian_rates(&sample, 1.0 / 64.0, support).unwrap();
            for (actual, exact) in rate.into_iter().zip([3.0, -2.0, 0.000002]) {
                // A cubic has zero fourth-order truncation error; the bound
                // covers binary64 evaluation and the weighted subtraction.
                assert!((actual - exact).abs() <= exact.abs() * 1e-8);
            }
        }
    }

    #[test]
    fn no_sample_leaves_the_continuous_owned_support() {
        for support in [(-0.0001, 0.0002), (0.0, 0.0002), (-0.0001, 0.0)] {
            let sample = |t: f64| {
                assert!(support.0 <= t && t <= support.1);
                Ok([t, t * t])
            };
            let rate = angular_rates(&sample, 1.0 / 64.0, support).unwrap();
            assert!((rate[0] - 1.0).abs() < 1e-14);
            assert!(rate[1].abs() < 1e-16);
        }
    }

    #[test]
    fn owned_record_admission_halves_only_for_nonlinear_retarded_domain() {
        let domain = (0.0, 1.0 / 64.0);
        let seen = std::cell::RefCell::new(Vec::new());
        let sample = |offset: f64| {
            // A varying emission clock reaches the record end sooner than
            // a linear observation-time support estimate predicts.
            let emission = offset + offset * offset;
            if emission < domain.0 || emission > domain.1 {
                return Ok(None);
            }
            seen.borrow_mut().push(offset);
            Ok(Some([10.0 + 0.000003 * offset, 7.0 - 0.000004 * offset]))
        };
        let rate = angular_rates_on_owned_domain(&sample, 1.0 / 64.0, domain).unwrap();
        assert!((rate[0] - 0.000003).abs() < 0.000003 * 1e-6);
        assert!((rate[1] + 0.000004).abs() < 0.000004 * 1e-6);
        assert!(seen.borrow().iter().all(|&offset| {
            let emission = offset + offset * offset;
            domain.0 <= emission && emission <= domain.1
        }));
    }

    #[test]
    fn admitted_in_domain_calculation_failure_is_not_a_step_retry() {
        let calls = std::cell::Cell::new(0);
        let sample = |_: f64| {
            calls.set(calls.get() + 1);
            Err(Error::SpkBodyNotFound { target: 499 })
        };
        assert!(matches!(
            angular_rates_on_owned_domain(&sample, 1.0 / 64.0, (-1.0, 1.0)),
            Err(Error::SpkBodyNotFound { target: 499 })
        ));
        assert_eq!(calls.get(), 1);
    }

    #[test]
    fn cartesian_tangent_jacobian_retains_small_latitude_rate_on_short_owned_record() {
        let central = [1.0, 1.0, 2.0_f64.powi(-20)];
        let velocity = [0.001, 0.002, 2.0_f64.powi(-24)];
        let support = (0.0, 2.0_f64.powi(-11));
        let sample = |offset: f64| {
            assert!(support.0 <= offset && offset <= support.1);
            Ok(Some(std::array::from_fn(|i| {
                velocity[i].mul_add(offset, central[i])
            })))
        };
        let rates = cartesian_angular_rates_on_owned_domain(&sample, 1.0 / 64.0, support).unwrap();
        let horizontal = central[0].hypot(central[1]);
        let radius_squared = horizontal * horizontal + central[2] * central[2];
        let horizontal_rate = (central[0] * velocity[0] + central[1] * velocity[1]) / horizontal;
        let longitude = ((central[0] * velocity[1] - central[1] * velocity[0])
            / (horizontal * horizontal))
            .to_degrees();
        let latitude = ((horizontal * velocity[2] - central[2] * horizontal_rate) / radius_squared)
            .to_degrees();
        assert!((rates[0] - longitude).abs() <= longitude.abs() * 1e-6);
        assert!((rates[1] - latitude).abs() <= latitude.abs() * 1e-6);
    }
}

impl std::ops::Neg for AngularWord {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}

/// Compensated projection and spherical/radial Jacobian. The frame coefficients
/// arrive from the actual analytic construction before coefficient rounding;
/// available value/rate residuals survive the projection and Jacobian. Fixed
/// J2000 output still has an exactly zero matrix derivative.
#[allow(clippy::too_many_arguments)]
pub(crate) fn expanded_output_rates(
    matrix: [[AngularWord; 3]; 3],
    matrix_rate: [[AngularWord; 3]; 3],
    apparent: [AngularWord; 3],
    velocity: [AngularWord; 3],
    want_xyz: bool,
) -> Result<(f64, f64, f64)> {
    expanded_motion_rates(
        expanded_output_motion(matrix, matrix_rate, apparent, velocity)?,
        want_xyz,
    )
}

/// Word motion after the actual output-frame projection, before any tuple reduction.
#[derive(Clone, Copy)]
pub(crate) struct OutputWordMotion {
    pub(crate) position: [AngularWord; 3],
    pub(crate) velocity: [AngularWord; 3],
}

/// Project the same motion once, retaining its fixed-size arithmetic residuals.
pub(super) fn expanded_output_motion(
    matrix: [[AngularWord; 3]; 3],
    matrix_rate: [[AngularWord; 3]; 3],
    apparent: [AngularWord; 3],
    velocity: [AngularWord; 3],
) -> Result<OutputWordMotion> {
    let dot = |row: [AngularWord; 3], vector: [AngularWord; 3]| {
        row.into_iter()
            .zip(vector)
            .fold(AngularWord::from(0.0), |sum, (a, b)| sum + a * b)
    };
    let position = matrix.map(|row| dot(row, apparent));
    let rate = std::array::from_fn(|i| dot(matrix[i], velocity) + dot(matrix_rate[i], apparent));
    crate::positions::check_finite_vec3(position.map(AngularWord::rounded))?;
    crate::positions::check_finite_vec3(rate.map(AngularWord::rounded))?;
    Ok(OutputWordMotion {
        position,
        velocity: rate,
    })
}

/// Reduce a projected word motion with the existing angular/radial Jacobian.
pub(super) fn expanded_motion_rates(
    motion: OutputWordMotion,
    want_xyz: bool,
) -> Result<(f64, f64, f64)> {
    let position = motion.position;
    let rate = motion.velocity;
    if want_xyz {
        let rate = rate.map(AngularWord::rounded);
        return Ok((rate[0], rate[1], rate[2]));
    }
    let largest = position
        .map(AngularWord::rounded)
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max);
    if largest == 0.0 {
        return Ok((0.0, 0.0, 0.0));
    }
    let exponent = ((largest.to_bits() >> 52) & 0x7ff).max(1);
    let scale = AngularWord::from(f64::from_bits(exponent << 52));
    let [x, y, z] = position.map(|value| value / scale);
    let [vx, vy, vz] = rate.map(|value| value / scale);
    let h2 = x * x + y * y;
    let r2 = h2 + z * z;
    if h2.rounded() == 0.0 || r2.rounded() == 0.0 {
        return Ok((0.0, 0.0, 0.0));
    }
    let longitude = (x * vy - y * vx) / h2;
    let latitude = (h2 * vz - z * (x * vx + y * vy)) / (r2 * h2.sqrt());
    let degrees = AngularWord {
        hi: 57.29577951308232,
        lo: -1.9878495670576283e-15,
    };
    let radial = (x * vx + y * vy + z * vz) / r2.sqrt() * scale;
    let result = (
        (longitude * degrees).rounded(),
        (latitude * degrees).rounded(),
        radial.rounded(),
    );
    crate::positions::check_finite_vec3([result.0, result.1, result.2])?;
    Ok(result)
}

#[cfg(test)]
mod expanded_output_tests {
    use super::*;

    #[test]
    fn output_jacobian_receives_the_central_value_residual() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let epsilon = 2.0_f64.powi(-54);
        let result = expanded_output_rates(
            identity.map(|row| row.map(AngularWord::from)),
            [[AngularWord::from(0.0); 3]; 3],
            [
                AngularWord {
                    hi: 1.0,
                    lo: epsilon,
                },
                AngularWord::from(1.0),
                AngularWord::from(0.0),
            ],
            [
                AngularWord::from(1.0),
                AngularWord::from(1.0),
                AngularWord::from(0.0),
            ],
            false,
        )
        .unwrap();
        // The determinant is epsilon, whereas discarding the value residual
        // before this Jacobian produces an exactly zero longitude derivative.
        let expected = epsilon * 0.5 * (180.0 / std::f64::consts::PI);
        assert!(result.0 > 0.0);
        assert!((result.0 - expected).abs() <= 4.0 * f64::EPSILON * expected);
    }
}

#[cfg(test)]
mod owned_record_seam_tests {
    use super::*;

    #[test]
    fn last_iteration_preserves_available_derivative_without_confirmation() {
        // The fifth-order term first reaches the unchanged resolution goal
        // on the last permitted comparison. There is no third scale left.
        let result = cartesian_rates_on_owned_domain(
            &|offset| Ok(Some([40.0 + offset + 524288.0 * offset.powi(5), 0.0, 0.0])),
            1.0 / 64.0,
            (0.0, 1.0),
        )
        .unwrap();
        // The closed derivative is 1 AU/day. Reaching the iteration cap
        // cannot turn an available finite estimate into a domain error.
        assert!(result.into_iter().all(f64::is_finite));
        assert!((result[0] - 1.0).abs() <= 1e-8);
        assert_eq!([result[1], result[2]], [0.0, 0.0]);
    }

    #[test]
    fn refused_confirmation_preserves_available_derivative() {
        let result = cartesian_rates_on_owned_domain(
            &|offset| {
                if offset != 0.0 && offset <= 1.0 / 256.0 {
                    Ok(None)
                } else {
                    Ok(Some([40.0 + offset, 0.0, 0.0]))
                }
            },
            1.0 / 64.0,
            (0.0, 1.0),
        )
        .unwrap();
        assert_eq!(
            result.map(f64::to_bits),
            [1.0_f64, 0.0, 0.0].map(f64::to_bits)
        );
    }

    #[test]
    fn underflowed_confirmation_preserves_available_derivative() {
        let result = cartesian_rates_on_owned_domain(
            &|offset| Ok(Some([offset, 0.0, 0.0])),
            f64::from_bits(2),
            (0.0, 1.0),
        )
        .unwrap();
        assert_eq!(
            result.map(f64::to_bits),
            [1.0_f64, 0.0, 0.0].map(f64::to_bits)
        );
    }

    #[test]
    fn one_sided_radius_requires_a_third_scale_to_confirm_resolution() {
        // r'(0)=1. For the exact polynomial with a=-496/3, the forward
        // fourth-order stencils at h=1/64 and h/2 both equal 1+2^-17:
        // their t^5 and t^6 truncation terms cancel in the first comparison.
        // Binary64 coefficient rounding is far below the unchanged goal.
        let support = (-2.0_f64.powi(-30), 1.0);
        let rate = cartesian_rates_on_owned_domain(
            &|offset| {
                assert!(support.0 < offset && offset < support.1);
                Ok(Some([
                    40.0 + offset - (496.0 / 3.0) * offset.powi(5) + 1024.0 * offset.powi(6),
                    0.0,
                    0.0,
                ]))
            },
            1.0 / 64.0,
            support,
        )
        .unwrap()[0];
        assert!(
            (rate - 1.0).abs() <= 1e-8,
            "unconfirmed radius slope: {rate}"
        );
    }

    #[test]
    fn stable_cross_record_radius_secants_do_not_certify_the_local_derivative() {
        // Two exact affine pieces meet continuously. The central epoch is
        // strictly inside the right record, whose radius slope is 1 AU/day.
        // The left record's slope is 3 AU/day. A centered stencil crossing
        // both records estimates their mean, even as its three estimates agree.
        let seam = -2.0_f64.powi(-30);
        let radius = |offset: f64| 40.0 + offset + 2.0 * (offset - seam).min(0.0);
        let support = (seam, 1.0);
        let seen = std::cell::RefCell::new(Vec::new());
        let local = cartesian_rates_on_owned_domain(
            &|offset| {
                assert!(support.0 < offset && offset < support.1);
                seen.borrow_mut().push(offset);
                Ok(Some([radius(offset), 0.0, 0.0]))
            },
            1.0 / 64.0,
            support,
        )
        .unwrap()[0];
        assert_eq!(local, 1.0);
        assert!(seen.borrow().iter().all(|&offset| offset >= 0.0));

        // This is the existing diagnostic's independent eight-point stencil,
        // not an alternative production derivative. No step is fitted to data.
        let public_secant = |denominator: f64| {
            let sum = [
                (1.0, 4.0 / 5.0),
                (2.0, -1.0 / 5.0),
                (3.0, 4.0 / 105.0),
                (4.0, -1.0 / 280.0),
            ]
            .into_iter()
            .map(|(n, weight)| weight * (radius(n / denominator) - radius(-n / denominator)))
            .sum::<f64>();
            sum * denominator
        };
        let secants = [32.0, 64.0, 16.0].map(public_secant);
        let low = secants.into_iter().fold(f64::INFINITY, f64::min);
        let high = secants.into_iter().fold(f64::NEG_INFINITY, f64::max);
        assert!(high - low < 1e-6);
        assert!(secants.into_iter().all(|rate| (rate - 2.0).abs() < 1e-6));
        assert!(secants.into_iter().all(|rate| (rate - local).abs() > 0.9));
    }
}
