//! First-order kinematics of the same finite-source PPN and aberration models.
//!
//! Stack-only tangents propagate TT-day derivatives through central geometry.
//! Closest-approach clipping differentiates its actual active branch; deflector
//! velocities include that acquisition epoch's rate. No neighboring apparent
//! places, polynomial extrapolation or reduced astronomical model is used.

use super::derivative::AngularWord as Word;
use super::state_reads::CorrectionStateReads;
use super::*;
use std::ops::{Add, Div, Mul, Neg, Sub};

#[derive(Clone, Copy)]
struct Tangent {
    value: Word,
    rate: Word,
    // Branch decisions retain the existing scalar central arithmetic. The
    // expanded value cannot silently select a different clipping branch.
    scalar: f64,
}

impl Tangent {
    fn fixed(value: f64) -> Self {
        Self {
            value: Word::from(value),
            rate: Word::from(0.0),
            scalar: value,
        }
    }

    /// Apply one already computed Word inverse without regrouping the rate.
    /// Shared vector quotients keep their original scalar divisions and branches.
    fn divide_with_inverse(self, denominator: Self, inverse: Word) -> Self {
        let value = self.value * inverse;
        Self {
            value,
            rate: (self.rate - value * denominator.rate) * inverse,
            scalar: self.scalar / denominator.scalar,
        }
    }

    fn sqrt(self) -> Self {
        let value = self.value.sqrt();
        Self {
            value,
            rate: self.rate / (Word::from(2.0) * value),
            scalar: self.scalar.sqrt(),
        }
    }
}
impl Add for Tangent {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self {
            value: self.value + rhs.value,
            rate: self.rate + rhs.rate,
            scalar: self.scalar + rhs.scalar,
        }
    }
}
impl Sub for Tangent {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self {
            value: self.value - rhs.value,
            rate: self.rate - rhs.rate,
            scalar: self.scalar - rhs.scalar,
        }
    }
}
impl Mul for Tangent {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        Self {
            value: self.value * rhs.value,
            rate: self.rate * rhs.value + self.value * rhs.rate,
            scalar: self.scalar * rhs.scalar,
        }
    }
}
impl Div for Tangent {
    type Output = Self;
    fn div(self, rhs: Self) -> Self {
        let value = self.value / rhs.value;
        Self {
            value,
            rate: (self.rate - value * rhs.rate) / rhs.value,
            scalar: self.scalar / rhs.scalar,
        }
    }
}
impl Neg for Tangent {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            value: -self.value,
            rate: -self.rate,
            scalar: -self.scalar,
        }
    }
}

type Motion = [Tangent; 3];
type State = ([f64; 3], [f64; 3], [Word; 3], [Word; 3]);

#[cfg(test)]
fn moving(position: [f64; 3], rate: [f64; 3]) -> Motion {
    moving_words(position.map(Word::from), rate.map(Word::from), position)
}
fn moving_words(position: [Word; 3], rate: [Word; 3], scalar: [f64; 3]) -> Motion {
    std::array::from_fn(|i| Tangent {
        value: position[i],
        rate: rate[i],
        scalar: scalar[i],
    })
}
fn dot(a: Motion, b: Motion) -> Tangent {
    let value = (0..3).fold(Word::from(0.0), |sum, i| sum + a[i].value * b[i].value);
    let rate = (0..3).fold(Word::from(0.0), |sum, i| {
        sum + a[i].rate * b[i].value + a[i].value * b[i].rate
    });
    Tangent {
        value,
        rate,
        scalar: compensated_dot3(a.map(|x| x.scalar), b.map(|x| x.scalar)),
    }
}

fn norm(a: Motion) -> Tangent {
    dot(a, a).sqrt()
}

fn aberration(ray: Motion, velocity: Motion, delay: Tangent) -> Motion {
    aberration_with_c(ray, velocity, delay, C_LIGHT_AU_DAY)
}

fn aberration_with_c(ray: Motion, velocity: Motion, delay: Tangent, light_speed: f64) -> Motion {
    let c = Tangent::fixed(light_speed);
    let one = Tangent::fixed(1.0);
    let distance = norm(ray);
    if distance.scalar == 0.0 {
        return ray;
    }
    if delay.scalar > 0.0 {
        // Algebraically cancel |v| from beta*cos(angle). This differentiates
        // the Lorentz formula continuously even at zero observer velocity.
        let beta2 = dot(velocity, velocity) / (c * c);
        let p = dot(ray, velocity) / (delay * c * c);
        let gamma_inverse = (one - beta2).sqrt();
        let q = (one + p / (one + gamma_inverse)) * delay;
        let mut denominator = one + p;
        if denominator.scalar.abs() < ABERRATION_R_FLOOR {
            denominator = one;
        }
        let inverse = Word::from(1.0) / denominator.value;
        std::array::from_fn(|i| {
            (gamma_inverse * ray[i] + q * velocity[i]).divide_with_inverse(denominator, inverse)
        })
    } else {
        let unit = ray.map(|x| x / distance);
        let beta = velocity.map(|x| x / c);
        let projection = dot(unit, beta);
        let shifted = std::array::from_fn(|i| unit[i] + beta[i] - unit[i] * projection);
        let length = norm(shifted);
        if length.scalar == 0.0 {
            return ray;
        }
        shifted.map(|x| x / length * distance)
    }
}

fn deflection(
    ray: Motion,
    observer: Motion,
    epoch: f64,
    delay: Tangent,
    lookup: &dyn Fn(i32, f64) -> Result<Option<State>>,
) -> Result<Motion> {
    deflection_with_scalar_dot(ray, observer, epoch, delay, lookup, dot)
}

/// Share the physical tangent formula while retaining each scalar acquisition's
/// written branch arithmetic. Stellar legacy dot products are ordinary sums;
/// the existing planetary analytic path keeps its compensated scalar dots.
fn deflection_with_scalar_dot(
    mut ray: Motion,
    observer: Motion,
    epoch: f64,
    delay: Tangent,
    lookup: &dyn Fn(i32, f64) -> Result<Option<State>>,
    scalar_dot: fn(Motion, Motion) -> Tangent,
) -> Result<Motion> {
    let norm = |a| scalar_dot(a, a).sqrt();
    let dot = scalar_dot;
    let distance = norm(ray);
    if distance.scalar == 0.0 {
        return Ok(ray);
    }
    let one = Tangent::fixed(1.0);
    let c = Tangent::fixed(C_LIGHT_AU_DAY);
    for (body, mass_ratio) in DEFLECTORS {
        let Some((position, _, position_words, velocity_words)) = lookup(body, epoch)? else {
            continue;
        };
        let at_reception = moving_words(position_words, velocity_words, position);
        let observer_to_body = std::array::from_fn(|i| at_reception[i] - observer[i]);
        // Keep this reciprocal after the original successful reception lookup.
        let distance_inverse = Word::from(1.0) / distance.value;
        let unit = ray.map(|x| x.divide_with_inverse(distance, distance_inverse));
        let projected_delay = dot(unit, observer_to_body) / c;
        let closest = if projected_delay.scalar < 0.0 {
            Tangent::fixed(0.0)
        } else if projected_delay.scalar > delay.scalar {
            delay
        } else {
            projected_delay
        };
        let Some((position, _, position_words, velocity_words)) =
            lookup(body, epoch - closest.scalar)?
        else {
            continue;
        };
        let at_closest = moving_words(
            position_words,
            velocity_words.map(|v| v * (Word::from(1.0) - closest.rate)),
            position,
        );
        let body_to_observer = std::array::from_fn(|i| observer[i] - at_closest[i]);
        let body_to_target = std::array::from_fn(|i| ray[i] + body_to_observer[i]);
        let source_range = norm(body_to_target);
        let observer_range = norm(body_to_observer);
        if source_range.scalar == 0.0 || observer_range.scalar < DEFLECTION_EMAG_FLOOR_AU {
            continue;
        }
        let source_inverse = Word::from(1.0) / source_range.value;
        let source_direction =
            body_to_target.map(|x| x.divide_with_inverse(source_range, source_inverse));
        let observer_inverse = Word::from(1.0) / observer_range.value;
        let observer_direction =
            body_to_observer.map(|x| x.divide_with_inverse(observer_range, observer_inverse));
        let e_dot_p = dot(observer_direction, unit);
        if e_dot_p.scalar.abs() > DEFLECTION_EDOTP_LIMIT {
            continue;
        }
        let p_dot_q = dot(unit, source_direction);
        let q_dot_e = dot(source_direction, observer_direction);
        let strength = Tangent::fixed(2.0 * GRAV_GS)
            / (Tangent::fixed(SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S)
                * observer_range
                * Tangent::fixed(AU_IN_M)
                * Tangent::fixed(mass_ratio));
        let squared = observer_range * observer_range;
        let limit = Tangent::fixed(DEFLECTION_DLIM_SCALE)
            / if squared.scalar > 1.0 { squared } else { one };
        let mut denominator = one + q_dot_e;
        if denominator.scalar < limit.scalar {
            denominator = limit;
        }
        let coefficient = strength / denominator * distance;
        ray = std::array::from_fn(|i| {
            ray[i] + coefficient * (p_dot_q * observer_direction[i] - e_dot_p * source_direction[i])
        });
    }
    Ok(ray)
}

/// Central source expansions for the existing apparent kinematic formulas.
/// Positions are ICRS AU; velocities are AU per actual TT observation day.
pub(super) struct ExpandedCenter {
    pub position: [Word; 3],
    pub velocity: [Word; 3],
    pub observer: [Word; 3],
    pub observer_velocity: [Word; 3],
}

/// Carry source value/rate residuals through the same PPN and Lorentz chain.
/// Lookup epochs and clipping gates remain the original scalar branches.
#[allow(clippy::too_many_arguments)]
pub(super) fn apparent_motion<S: CorrectionStateReads + ?Sized>(
    provider: &S,
    center: &StencilCenter,
    epoch: f64,
    delay: f64,
    delay_rate: f64,
    observer_acceleration: [f64; 3],
    required_lookup: bool,
    seed: ExpandedCenter,
    scalar_samples: Option<&super::serving::PpnScalarSamples>,
) -> Result<([Word; 3], [Word; 3])> {
    let mut ray = moving_words(seed.position, seed.velocity, center.geo_geom);
    let delay = Tangent {
        value: Word::from(delay),
        rate: Word::from(delay_rate),
        scalar: delay,
    };
    if center.do_defl {
        let scalar_at = |body, tt| -> Result<GeometricState> {
            match scalar_samples.and_then(|samples| samples.state_at(body, tt)) {
                Some(state) => Ok(state),
                None => provider.geometric_state(body, tt),
            }
        };
        let lookup = |body, tt| {
            let scalar = if required_lookup {
                let state = scalar_at(body, tt)?;
                Some((state.pos(), state.vel()))
            } else {
                let captured = std::cell::Cell::new(None);
                super::frames::deflector_position(
                    &|id, time| {
                        let state = scalar_at(id, time)?;
                        captured.set(Some((state.pos(), state.vel())));
                        Ok(state.pos())
                    },
                    body,
                    tt,
                )?
                .and(captured.get())
            };
            let Some((position, velocity)) = scalar else {
                return Ok(None);
            };
            let (p, v) =
                provider.geometric_state_word_parts_two_part(body, TwoPartEpoch::split(tt)?)?;
            let expand = |pair: [f64; 2]| Word {
                hi: pair[0],
                lo: pair[1],
            };
            Ok(Some((position, velocity, p.map(expand), v.map(expand))))
        };
        ray = deflection(
            ray,
            moving_words(seed.observer, seed.observer_velocity, center.observer),
            epoch,
            delay,
            &lookup,
        )?;
    }
    if center.do_aber {
        ray = aberration(
            ray,
            moving_words(
                seed.observer_velocity,
                observer_acceleration.map(Word::from),
                center.observer_vel,
            ),
            delay,
        );
    }
    let position = ray.map(|x| x.value);
    let velocity = ray.map(|x| x.rate);
    check_finite_vec3(position.map(Word::rounded))?;
    check_finite_vec3(velocity.map(Word::rounded))?;
    Ok((position, velocity))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn excluded_deflector_preserves_motion_and_lookups_without_unused_dots() {
        std::thread_local! {
            static DOT_CALLS: std::cell::Cell<usize> = const { std::cell::Cell::new(0) };
        }
        fn counted_dot(a: Motion, b: Motion) -> Tangent {
            DOT_CALLS.with(|count| count.set(count.get() + 1));
            dot(a, b)
        }

        DOT_CALLS.with(|count| count.set(0));
        let mut ray = moving([1.0, 0.0, -0.0], [0.125, -0.0625, 0.0]);
        ray[0].value.lo = 2.0_f64.powi(-54);
        ray[0].rate.lo = 2.0_f64.powi(-58);
        let epoch: f64 = 2451545.0;
        let lookups = std::cell::RefCell::new(Vec::new());
        // The aligned directions reach the existing exclusion after both Sun reads.
        let actual = deflection_with_scalar_dot(
            ray,
            moving([0.0; 3], [0.0; 3]),
            epoch,
            Tangent::fixed(1.0 / C_LIGHT_AU_DAY),
            &|body, tt| {
                lookups.borrow_mut().push((body, tt.to_bits()));
                Ok((body == SUN).then(|| {
                    let position = [-1.0, 0.0, 0.0];
                    let velocity = [0.0; 3];
                    (
                        position,
                        velocity,
                        position.map(Word::from),
                        velocity.map(Word::from),
                    )
                }))
            },
            counted_dot,
        )
        .unwrap();
        let bits = |motion: Motion| {
            motion.map(|value| {
                [
                    value.value.hi.to_bits(),
                    value.value.lo.to_bits(),
                    value.rate.hi.to_bits(),
                    value.rate.lo.to_bits(),
                    value.scalar.to_bits(),
                ]
            })
        };
        assert_eq!(bits(actual), bits(ray));
        assert_eq!(
            lookups.into_inner(),
            vec![
                (SUN, epoch.to_bits()),
                (SUN, epoch.to_bits()),
                (JUPITER, epoch.to_bits()),
                (SATURN, epoch.to_bits()),
            ]
        );
        // Initial range, delay projection, two ranges and e_dot_p; the other two dots are unused.
        assert_eq!(DOT_CALLS.with(|count| count.get()), 5);
    }

    #[test]
    fn moving_deflector_and_observer_acceleration_match_high_precision_places() {
        // Independent 100-decimal symmetric place differences at 1e-20 and
        // 1e-30 TT day agree to these words. The reference uses the scalar
        // finite-source geometry and norm/cosine Lorentz formula, not tangents.
        let cases = [
            (
                false,
                true,
                [
                    0.019998416326847842,
                    -0.030002707139256863,
                    0.009999788371848211,
                ],
            ),
            (
                true,
                false,
                [
                    0.01999999974957686,
                    -0.029999999974511457,
                    0.009999999939786699,
                ],
            ),
            (
                true,
                true,
                [
                    0.01999841607641417,
                    -0.03000270711372821,
                    0.009999788311626113,
                ],
            ),
        ];
        let position = [1.0, 0.2, 0.3];
        let velocity = [0.02, -0.03, 0.01];
        let range = vec3_dist(position);
        let delay = Tangent {
            value: Word::from(range / C_LIGHT_AU_DAY),
            rate: Word::from(compensated_dot3(position, velocity) / (range * C_LIGHT_AU_DAY)),
            scalar: range / C_LIGHT_AU_DAY,
        };
        for (ppn, aber, expected) in cases {
            let mut ray = moving(position, velocity);
            if ppn {
                ray = deflection(
                    ray,
                    moving([0.1, -0.2, 0.3], [0.002, -0.01, 0.003]),
                    0.0,
                    delay,
                    &|body, t| {
                        Ok((body == SUN).then(|| {
                            (
                                [0.3 + 0.001 * t, 1.0 + 0.002 * t, 0.2 - 0.003 * t],
                                [0.001, 0.002, -0.003],
                                [0.3 + 0.001 * t, 1.0 + 0.002 * t, 0.2 - 0.003 * t].map(Word::from),
                                [0.001, 0.002, -0.003].map(Word::from),
                            )
                        }))
                    },
                )
                .unwrap();
            }
            if aber {
                ray = aberration(
                    ray,
                    moving([0.002, -0.01, 0.003], [0.0003, -0.0002, 0.0001]),
                    delay,
                );
            }
            for (actual, expected) in ray.into_iter().zip(expected) {
                assert!(
                    (actual.rate.rounded() - expected).abs()
                        <= 16.0 * f64::EPSILON * expected.abs()
                );
            }
        }
    }

    #[test]
    fn transverse_observer_acceleration_has_the_exact_zero_velocity_limit() {
        let actual = aberration(
            moving([1.0, 0.0, 0.0], [0.0; 3]),
            moving([0.0; 3], [0.0, 0.03, 0.0]),
            Tangent::fixed(1.0 / C_LIGHT_AU_DAY),
        );
        assert_eq!(actual.map(|x| x.value.rounded()), [1.0, 0.0, 0.0]);
        let expected = 0.03 / C_LIGHT_AU_DAY;
        assert!((actual[1].rate.rounded() - expected).abs() <= 2.0 * f64::EPSILON * expected);
        assert_eq!(actual[0].rate.rounded(), 0.0);
        assert_eq!(actual[2].rate.rounded(), 0.0);
    }
}

#[cfg(test)]
mod expanded_correction_tests {
    use super::*;

    #[test]
    fn quotient_value_residual_contributes_to_the_final_derivative() {
        let epsilon = 2.0_f64.powi(-54);
        let numerator = Tangent {
            value: Word {
                hi: 1.0,
                lo: epsilon,
            },
            rate: Word::from(0.0),
            scalar: 1.0,
        };
        let denominator = Tangent {
            value: Word::from(1.0),
            rate: Word::from(1.0),
            scalar: 1.0,
        };
        let cancel = Tangent {
            value: Word::from(0.0),
            rate: Word::from(1.0),
            scalar: 0.0,
        };
        let result = numerator / denominator + cancel;
        assert_eq!(result.scalar, 1.0);
        assert_eq!(result.rate.rounded(), -epsilon);
    }
}

/// Complete native static-catalog anchor value and rate, degrees per TT day.
/// This differentiates the actual legacy PPN and stellar-speed aberration law;
/// it does not replace scalar positions or opaque catalog/ray callbacks.
pub(crate) fn static_anchor_motion<P: StateProvider>(
    provider: &P,
    data: &crate::stars::StarData,
    tt: f64,
    noaberr: bool,
    nogdefl: bool,
    wilhelm: bool,
) -> Result<(f64, f64)> {
    let w = Word::from;
    let fixed = Tangent::fixed;
    let zero = fixed(0.0);
    let (initial, velocity) = catalog_words(data)?;
    let earth = provider.geometric_state(EARTH, tt)?;
    let (ep, ev) = provider.geometric_state_word_parts_two_part(EARTH, TwoPartEpoch::split(tt)?)?;
    let expand = |pair: [f64; 2]| Word {
        hi: pair[0],
        lo: pair[1],
    };
    let observer = moving_words(ep.map(expand), ev.map(expand), earth.pos());
    let initial_motion = initial.map(|value| Tangent {
        value,
        rate: w(0.0),
        scalar: value.rounded(),
    });
    let direction = initial_motion.map(|v| v / norm(initial_motion));
    // The source clock is the complete seven-term TT->TDB realization. Its
    // value residual matters for proper motion, so do not round the JD first.
    let t = Tangent {
        value: (w(tt) - w(FRAME_J2000)) / w(36525.0),
        rate: w(1.0) / w(36525.0),
        scalar: (tt - FRAME_J2000) / 36525.0,
    };
    let mut seconds = zero;
    for (a, b, c) in [
        (0.001657, 628.3076, 6.2401),
        (0.000022, 575.3385, 4.2970),
        (0.000014, 1256.6152, 6.1969),
        (0.000005, 606.9777, 4.0212),
        (0.000005, 52.9691, 0.4444),
        (0.000002, 21.3299, 5.5431),
    ] {
        seconds = seconds + fixed(a) * tangent_sin_cos(fixed(b) * t + fixed(c))?.0;
    }
    seconds =
        seconds + fixed(0.000010) * t * tangent_sin_cos(fixed(628.3076) * t + fixed(4.2490))?.0;
    let span = Tangent {
        value: w(tt) - w(FRAME_J2000),
        rate: w(1.0),
        scalar: tt - FRAME_J2000,
    } + seconds / fixed(86400.0)
        + dot(direction, observer) / fixed(FRAME_LIGHT_SPEED_AU_DAY);
    let mut astro: Motion = std::array::from_fn(|i| {
        initial_motion[i]
            + Tangent {
                value: velocity[i],
                rate: w(0.0),
                scalar: velocity[i].rounded(),
            } * span
            - observer[i]
    });
    // Scalar branch operands remain the actual catalog propagation, including
    // the original rounded scalar JD addition before its J2000 subtraction.
    let scalar_epoch = crate::stars::topo_star_epoch_days(data, tt, earth.pos())?;
    let scalar_position = crate::stars::star_icrs_state_catalog_angles(data, scalar_epoch).0;
    for i in 0..3 {
        astro[i].scalar = scalar_position[i] - earth.pos()[i];
    }
    let legacy_dot = |a: Motion, b: Motion| {
        let mut value = dot(a, b);
        value.scalar =
            a[0].scalar * b[0].scalar + a[1].scalar * b[1].scalar + a[2].scalar * b[2].scalar;
        value
    };
    let delay = legacy_dot(astro, astro).sqrt() / fixed(FRAME_LIGHT_SPEED_AU_DAY);
    let aberrate = |ray| -> Result<Motion> {
        let acceleration =
            provider
                .geometric_acceleration_tt(EARTH, tt)?
                .ok_or(Error::ConfigurationMissing {
                    missing_config: "anchor_analytic_acceleration",
                    suggestion: "native analytic anchors require complete source acceleration",
                })?;
        Ok(aberration_with_c(
            ray,
            moving_words(ev.map(expand), acceleration.map(w), earth.vel()),
            delay,
            FRAME_LIGHT_SPEED_AU_DAY,
        ))
    };
    let ray = if noaberr && nogdefl {
        astro
    } else if nogdefl {
        aberrate(astro)?
    } else {
        let lookup = |body, epoch| {
            let captured = std::cell::Cell::new(None);
            let present = super::frames::deflector_position(
                &|id, time| {
                    let state = provider.geometric_state(id, time)?;
                    captured.set(Some((state.pos(), state.vel())));
                    Ok(state.pos())
                },
                body,
                epoch,
            )?
            .and(captured.get());
            let Some((position, velocity)) = present else {
                return Ok(None);
            };
            let (p, v) =
                provider.geometric_state_word_parts_two_part(body, TwoPartEpoch::split(epoch)?)?;
            Ok(Some((position, velocity, p.map(expand), v.map(expand))))
        };
        let full = aberrate(deflection_with_scalar_dot(
            astro, observer, tt, delay, &lookup, legacy_dot,
        )?)?;
        if noaberr {
            let a = aberrate(astro)?;
            std::array::from_fn(|i| full[i] - (a[i] - astro[i]))
        } else {
            full
        }
    };
    let (frame, cos_eps, sin_eps) = anchor_frame_motion(tt)?;
    let eq: Motion = std::array::from_fn(|i| dot(frame[i], ray));
    // The hour-circle formula cancels the positive equatorial horizontal
    // norm algebraically; this is the same atan2(sin(alpha),cos(alpha)*cos(eps)).
    let (x, y) = if wilhelm {
        (eq[0] * cos_eps, eq[1])
    } else {
        (eq[0], eq[1] * cos_eps + eq[2] * sin_eps)
    };
    let denominator = x.value * x.value + y.value * y.value;
    if denominator.rounded() == 0.0 {
        return Err(Error::MathDomain { value: 0.0 });
    }
    let rate = (x.value * y.rate - y.value * x.rate) / denominator;
    let degrees = Word {
        hi: 57.29577951308232,
        lo: -1.9878495670576283e-15,
    };
    Ok((
        y.value.rounded().atan2(x.value.rounded()).to_degrees(),
        (rate * degrees).rounded(),
    ))
}

/// The immutable catalog constructor before intermediate elementary rounding.
fn catalog_words(data: &crate::stars::StarData) -> Result<([Word; 3], [Word; 3])> {
    let w = Word::from;
    let tau = w(std::f64::consts::TAU);
    let (sr, cr) = super::frame_velocity::frame_sin_cos_words(
        w(data.ra_j2000_deg()) / w(15.0) / w(24.0) * tau,
    )?;
    let (sd, cd) =
        super::frame_velocity::frame_sin_cos_words(w(data.dec_j2000_deg()) / w(360.0) * tau)?;
    let px = if data.parallax_mas() > 0.0 {
        data.parallax_mas()
    } else {
        1e-6
    };
    let distance = w(1.0)
        / super::frame_velocity::frame_sin_cos_words(
            w(px) * w(1e-3) * w(std::f64::consts::PI / 180.0 / 3600.0),
        )?
        .0;
    let initial = [distance * cd * cr, distance * cd * sr, distance * sd];
    let k = w(1.0) / (w(1.0) - w(data.radial_velocity_km_s()) / w(299792458.0) * w(1000.0));
    let ra = w(data.pm_ra_arcsec_yr()) * w(1000.0) / (w(px) * w(365.25)) * k;
    let dec = w(data.pm_dec_arcsec_yr()) * w(1000.0) / (w(px) * w(365.25)) * k;
    let radial = w(data.radial_velocity_km_s()) * w(86400.0) / w(149597870.7) * k;
    Ok((
        initial,
        [
            -ra * sr - dec * sd * cr + radial * cd * cr,
            ra * cr - dec * sd * sr + radial * cd * sr,
            dec * cd + radial * sd,
        ],
    ))
}

fn tangent_sin_cos(angle: Tangent) -> Result<(Tangent, Tangent)> {
    let (s, c) = super::frame_velocity::frame_sin_cos_words(angle.value)?;
    Ok((
        Tangent {
            value: s,
            rate: c * angle.rate,
            scalar: angle.scalar.sin(),
        },
        Tangent {
            value: c,
            rate: -s * angle.rate,
            scalar: angle.scalar.cos(),
        },
    ))
}

/// Full Vondrak pole pair, with an explicit angular scale for the distinct
/// ERFA mean-equator and direct-TT obliquity realizations.
#[allow(clippy::too_many_arguments)]
fn anchor_pole_pair(
    t: Tangent,
    periods: &[f64],
    ac: &[f64],
    as_: &[f64],
    bc: &[f64],
    bs: &[f64],
    poly: &[(f64, f64)],
    scale: f64,
) -> Result<(Tangent, Tangent)> {
    let fixed = Tangent::fixed;
    let mut a = fixed(0.0);
    let mut b = a;
    for i in 0..periods.len() {
        let (s, c) = tangent_sin_cos(fixed(std::f64::consts::TAU) * t / fixed(periods[i]))?;
        a = a + c * fixed(ac[i]) + s * fixed(as_[i]);
        b = b + c * fixed(bc[i]) + s * fixed(bs[i]);
    }
    let mut power = fixed(1.0);
    for &(x, y) in poly {
        a = a + fixed(x) * power;
        b = b + fixed(y) * power;
        power = power * t;
    }
    Ok((a * fixed(scale), b * fixed(scale)))
}
fn anchor_poles(t: Tangent, scale: f64, eps: Word) -> Result<(Motion, Motion)> {
    let (p, q) = anchor_pole_pair(
        t,
        &VONDRAK_PQPER_PERIODS,
        &VONDRAK_PQPER_P_COS,
        &VONDRAK_PQPER_P_SIN,
        &VONDRAK_PQPER_Q_COS,
        &VONDRAK_PQPER_Q_SIN,
        &VONDRAK_PQPOL,
        scale,
    )?;
    let (x, y) = anchor_pole_pair(
        t,
        &VONDRAK_XYPER_PERIODS,
        &VONDRAK_XYPER_X_COS,
        &VONDRAK_XYPER_X_SIN,
        &VONDRAK_XYPER_Y_COS,
        &VONDRAK_XYPER_Y_SIN,
        &VONDRAK_XYPOL,
        scale,
    )?;
    let one = Tangent::fixed(1.0);
    let zero = Tangent::fixed(0.0);
    let pole_z = |a: Tangent, b: Tangent| {
        let squared = one - a * a - b * b;
        if squared.scalar > 0.0 {
            squared.sqrt()
        } else {
            zero
        }
    };
    let z = pole_z(p, q);
    let (s, c) = super::frame_velocity::frame_sin_cos_words(eps)?;
    let constant = |value: Word| Tangent {
        value,
        rate: Word::from(0.0),
        scalar: value.rounded(),
    };
    Ok((
        [x, y, pole_z(x, y)],
        [
            p,
            -q * constant(c) - z * constant(s),
            -q * constant(s) + z * constant(c),
        ],
    ))
}
fn motion_cross(a: Motion, b: Motion) -> Motion {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        a[j] * b[k] - a[k] * b[j]
    })
}
fn anchor_frame_motion(tt: f64) -> Result<([Motion; 3], Tangent, Tangent)> {
    let w = Word::from;
    let t = Tangent {
        value: (w(tt) - w(FRAME_J2000)) / w(36525.0),
        rate: w(1.0) / w(36525.0),
        scalar: (tt - FRAME_J2000) / 36525.0,
    };
    let das2r = 4.848_136_811_095_36e-6;
    let (q, e) = anchor_poles(t, das2r, w(84381.406) * w(das2r))?;
    let raw = motion_cross(q, e);
    let length = norm(raw);
    let x = raw.map(|v| v / length);
    let rows = [x, motion_cross(q, x), q];
    let fixed = Tangent::fixed;
    let dx = fixed(-0.016617) * fixed(das2r);
    let de = fixed(-0.0068192) * fixed(das2r);
    let dr = fixed(-0.0146) * fixed(das2r);
    let frame = rows.map(|[a, b, c]| {
        [
            a - b * dr + c * dx,
            a * dr + b + c * de,
            -a * dx - b * de + c,
        ]
    });
    let (direct_q, direct_e) = anchor_poles(t, ARCSEC_TO_RAD, w(VONDRAK_EPS0_RAD))?;
    let cosine = dot(direct_q, direct_e);
    let sine = (fixed(1.0) - cosine * cosine).sqrt();
    Ok((frame, cosine, sine))
}

#[cfg(test)]
mod live_anchor_toggle_tests {
    use super::*;
    struct Linear;
    impl StateProvider for Linear {
        fn provider_name(&self) -> &'static str {
            "independent-linear-anchor-test"
        }
        fn has_body(&self, id: i32) -> bool {
            [EARTH, SUN, JUPITER, SATURN].contains(&id)
        }
        fn has_analytic_kinematics(&self) -> bool {
            true
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            None
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            Some(COORD_ICRS_BARY)
        }
        fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("explicit test state")
        }
        fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("explicit test state")
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
            panic!("metadata is not consumed")
        }
        fn geometric_acceleration_tt(&self, _: i32, _: f64) -> Result<Option<[f64; 3]>> {
            Ok(Some([0.0; 3]))
        }
        fn geometric_state(&self, id: i32, tt: f64) -> Result<GeometricState> {
            let (p, v) = match id {
                EARTH => ([0.1, -0.2, 0.3], [0.001, -0.01, 0.003]),
                SUN => ([0.3, 1.0, 0.2], [0.001, 0.002, -0.003]),
                JUPITER => ([0.5, 1.2, -0.7], [0.0; 3]),
                SATURN => ([4.0, 0.3, -0.5], [0.0; 3]),
                _ => return Err(Error::UnknownBody { body_id: id }),
            };
            Ok(GeometricState {
                body_id: id,
                jd_tt: tt,
                pos: std::array::from_fn(|i| p[i] + v[i] * (tt - FRAME_J2000)),
                vel: v,
                frame: StateFrame::IcrsBarycentric,
                planet_center: false,
                tier: None,
            })
        }
        fn geometric_state_two_part(&self, id: i32, epoch: TwoPartEpoch) -> Result<GeometricState> {
            self.geometric_state(id, epoch.value())
        }
    }
    #[test]
    fn live_anchor_all_four_toggle_compositions_match_independent_mp100() {
        // Full catalog propagation, moving Earth/deflectors, both frame models,
        // and finite-source PPN were independently differentiated at 100 digits
        // with 1e-20/1e-30-day steps. In particular NOABERR alone removes the
        // astrometric aberration displacement from the full apparent result.
        let cases = [
            (false, false, 3.818135765664671e-05),
            (true, false, 3.818134433062274e-05),
            (false, true, 3.818577675438103e-05),
            (true, true, 3.8185763435052386e-05),
        ];
        let data = crate::sidereal::anchor_star("SPICA").unwrap();
        for (noaberr, nogdefl, expected) in cases {
            let (_, rate) =
                static_anchor_motion(&Linear, &data, FRAME_J2000, noaberr, nogdefl, false).unwrap();
            assert!(
                (rate - expected).abs() < 1e-18,
                "{noaberr} {nogdefl}: {rate:.20e} {expected:.20e}"
            );
        }
    }
    #[test]
    fn live_anchor_geometric_pole_ignores_both_correction_toggles() {
        use crate::sidereal::{DynamicAnchor, DynamicAnchorProvider};
        let provider = crate::stars::SealedAnchorProvider::new(&Linear);
        let expected = provider
            .anchor_rate_deg_day(DynamicAnchor::GalNorthPole, FRAME_J2000, false, false)
            .unwrap();
        for a in [false, true] {
            for b in [false, true] {
                assert_eq!(
                    provider
                        .anchor_rate_deg_day(DynamicAnchor::GalNorthPole, FRAME_J2000, a, b)
                        .unwrap(),
                    expected
                );
            }
        }
    }
}

#[cfg(test)]
mod shared_aberration_inverse_tests {
    use super::*;

    fn assert_float_same(actual: f64, expected: f64) {
        if expected.is_nan() {
            // Arithmetic NaN payload/sign is not a portable Rust contract.
            assert!(actual.is_nan());
        } else {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
    }

    fn assert_word_same(actual: Word, expected: Word) {
        assert_float_same(actual.hi, expected.hi);
        assert_float_same(actual.lo, expected.lo);
    }

    /// The original positive-delay expression, using its serial quotients.
    /// This detects arithmetic/refusal changes; it is not a physical oracle.
    fn serial_positive_delay(
        ray: Motion,
        velocity: Motion,
        delay: Tangent,
        light_speed: f64,
    ) -> (Motion, bool) {
        assert!(delay.scalar > 0.0);
        assert_ne!(norm(ray).scalar, 0.0);
        let c = Tangent::fixed(light_speed);
        let one = Tangent::fixed(1.0);
        let beta2 = dot(velocity, velocity) / (c * c);
        let p = dot(ray, velocity) / (delay * c * c);
        let gamma_inverse = (one - beta2).sqrt();
        let q = (one + p / (one + gamma_inverse)) * delay;
        let mut denominator = one + p;
        let clamped = denominator.scalar.abs() < ABERRATION_R_FLOOR;
        if clamped {
            denominator = one;
        }
        (
            std::array::from_fn(|i| (gamma_inverse * ray[i] + q * velocity[i]) / denominator),
            clamped,
        )
    }

    #[test]
    fn shared_inverse_retains_exact_dyadic_value_and_moving_rate_residuals() {
        let epsilon = 2.0_f64.powi(-54);
        let values = [
            Word {
                hi: 1.0,
                lo: epsilon,
            },
            Word {
                hi: -2.0,
                lo: 2.0 * epsilon,
            },
            Word {
                hi: 0.25,
                lo: epsilon / 4.0,
            },
        ];
        let rates = [
            Word {
                hi: 0.625,
                lo: -epsilon / 8.0,
            },
            Word {
                hi: -0.125,
                lo: epsilon / 8.0,
            },
            Word {
                hi: 1.0,
                lo: epsilon / 2.0,
            },
        ];
        let expected_values = [
            Word {
                hi: 0.5,
                lo: epsilon / 2.0,
            },
            Word {
                hi: -1.0,
                lo: epsilon,
            },
            Word {
                hi: 0.125,
                lo: epsilon / 8.0,
            },
        ];
        // Exact dyadic evaluation of (rate - (value/2)*(1/2))/2.
        let expected_rates = [
            Word {
                hi: 0.1875,
                lo: -3.0 * epsilon / 16.0,
            },
            Word {
                hi: 0.1875,
                lo: -3.0 * epsilon / 16.0,
            },
            Word {
                hi: 0.46875,
                lo: 7.0 * epsilon / 32.0,
            },
        ];
        let denominator = Tangent {
            value: Word::from(2.0),
            rate: Word::from(0.5),
            scalar: 2.0,
        };
        let inverse = Word::from(1.0) / denominator.value;
        let numerators = moving_words(values, rates, values.map(Word::rounded));
        for i in 0..3 {
            let actual = numerators[i].divide_with_inverse(denominator, inverse);
            assert_word_same(actual.value, expected_values[i]);
            assert_word_same(actual.rate, expected_rates[i]);
            assert_float_same(actual.scalar, numerators[i].scalar / denominator.scalar);
        }

        let residual = Tangent {
            value: Word {
                hi: 1.0,
                lo: epsilon,
            },
            rate: Word::from(0.0),
            scalar: 1.0,
        }
        .divide_with_inverse(
            Tangent {
                value: Word::from(1.0),
                rate: Word::from(1.0),
                scalar: 1.0,
            },
            Word::from(1.0),
        ) + Tangent {
            value: Word::from(0.0),
            rate: Word::from(1.0),
            scalar: 0.0,
        };
        assert_eq!(residual.rate.rounded(), -epsilon);
    }

    #[test]
    fn positive_delay_identity_preserves_mixed_source_and_delay_words() {
        let epsilon = 2.0_f64.powi(-54);
        let values = [
            Word {
                hi: 1.0,
                lo: epsilon,
            },
            Word {
                hi: 0.5,
                lo: -epsilon / 8.0,
            },
            Word {
                hi: -0.25,
                lo: epsilon / 16.0,
            },
        ];
        let rates = [
            Word {
                hi: 0.125,
                lo: -epsilon / 16.0,
            },
            Word {
                hi: -0.5,
                lo: epsilon / 8.0,
            },
            Word {
                hi: 0.25,
                lo: epsilon / 16.0,
            },
        ];
        let ray = moving_words(values, rates, values.map(Word::rounded));
        let actual = aberration(
            ray,
            moving([0.0; 3], [0.0; 3]),
            Tangent {
                value: Word {
                    hi: 2.0,
                    lo: epsilon,
                },
                rate: Word {
                    hi: 0.25,
                    lo: -epsilon / 16.0,
                },
                scalar: 2.0,
            },
        );
        for i in 0..3 {
            assert_word_same(actual[i].value, values[i]);
            assert_word_same(actual[i].rate, rates[i]);
            assert_float_same(actual[i].scalar, ray[i].scalar);
        }
    }

    #[test]
    fn positive_delay_floor_and_extreme_refusal_match_serial_quotients() {
        let compare = |ray: Motion, velocity: Motion, clamped: bool, finite: bool| {
            let delay = Tangent::fixed(1.0);
            let (expected, was_clamped) = serial_positive_delay(ray, velocity, delay, 1.0);
            assert_eq!(was_clamped, clamped);
            let actual = aberration_with_c(ray, velocity, delay, 1.0);
            let all_finite = |motion: Motion| {
                motion.iter().all(|value| {
                    value.value.rounded().is_finite() && value.rate.rounded().is_finite()
                })
            };
            assert_eq!(all_finite(expected), finite);
            assert_eq!(all_finite(actual), finite);
            for i in 0..3 {
                assert_float_same(actual[i].value.rounded(), expected[i].value.rounded());
                assert_float_same(actual[i].rate.rounded(), expected[i].rate.rounded());
                assert_float_same(actual[i].scalar, expected[i].scalar);
            }
        };
        // The original scalar denominator is zero, even though the retained
        // ray residual makes the Word denominator nonzero. The same floor wins.
        let epsilon = 2.0_f64.powi(-54);
        compare(
            moving_words(
                [
                    Word {
                        hi: -2.0,
                        lo: epsilon,
                    },
                    Word::from(0.0),
                    Word::from(0.0),
                ],
                [Word::from(0.0); 3],
                [-2.0, 0.0, 0.0],
            ),
            moving([0.5, 0.0, 0.0], [0.0; 3]),
            true,
            true,
        );
        compare(
            moving([2.0, 0.0, 0.0], [0.0; 3]),
            moving([0.5, 0.0, 0.0], [0.0; 3]),
            false,
            true,
        );
        // Outside the Lorentz square-root domain both arithmetic routes must
        // remain nonfinite; no new clamp or successful finite answer is allowed.
        compare(
            moving([1.0, 0.0, 0.0], [0.0; 3]),
            moving([2.0, 0.0, 0.0], [0.0; 3]),
            false,
            false,
        );
    }

    #[test]
    fn shared_inverse_exposes_subnormal_divisor_range_limit() {
        // This is a deliberate adverse helper case, not a physical assertion
        // about the final aberration denominator after its existing floor.
        // It prevents describing the helper as a generic Word::Div identity.
        let smallest = f64::from_bits(1);
        let numerator = Tangent::fixed(smallest);
        let denominator = Tangent::fixed(smallest);
        let serial = numerator / denominator;
        let inverse = Word::from(1.0) / denominator.value;
        let shared = numerator.divide_with_inverse(denominator, inverse);
        assert_eq!(serial.value.rounded(), 1.0);
        assert!(shared.value.rounded().is_nan());
        assert_float_same(shared.scalar, serial.scalar);
    }

    #[test]
    fn actual_positive_delay_branch_emits_fixed_multiprecision_witnesses() {
        fn encode(value: Tangent) -> serde_json::Value {
            let bits = |value: f64| format!("{:016x}", value.to_bits());
            serde_json::json!({
                "value_bits": [bits(value.value.hi), bits(value.value.lo)],
                "rate_bits": [bits(value.rate.hi), bits(value.rate.lo)],
                "scalar_bits": bits(value.scalar),
                "value_rounded_bits": bits(value.value.rounded()),
                "rate_rounded_bits": bits(value.rate.rounded()),
            })
        }
        let emit =
            |id: &str, ray: Motion, observer: Motion, delay: Tangent, c: f64, floor: bool| {
                let one = Tangent::fixed(1.0);
                let speed = Tangent::fixed(c);
                let raw = one + dot(ray, observer) / (delay * speed * speed);
                let active = raw.scalar.abs() < ABERRATION_R_FLOOR;
                assert_eq!(active, floor);
                let selected = if active { one } else { raw };
                let candidate = aberration_with_c(ray, observer, delay, c);
                let (serial, serial_floor) = serial_positive_delay(ray, observer, delay, c);
                assert_eq!(active, serial_floor);
                assert!(delay.scalar > 0.0);
                for motion in [candidate, serial] {
                    assert!(motion.iter().all(|v| {
                        v.value.rounded().is_finite() && v.rate.rounded().is_finite()
                    }));
                }
                println!(
                    "SHARED_WORD_ABERRATION_BRANCH {}",
                    serde_json::json!({
                        "schema": 1,
                        "case": id,
                        "law": "positive_delay_linear_inputs_active_scalar_branch",
                        "light_speed_bits": format!("{:016x}", c.to_bits()),
                        "floor_bits": format!("{:016x}", ABERRATION_R_FLOOR.to_bits()),
                        "floor_active": active,
                        "raw_denominator": encode(raw),
                        "selected_denominator": encode(selected),
                        "ray": ray.map(encode),
                        "observer_velocity": observer.map(encode),
                        "delay": encode(delay),
                        "candidate": candidate.map(encode),
                        "serial": serial.map(encode),
                    })
                );
            };
        let epsilon = 2.0_f64.powi(-54);
        let ray = moving_words(
            [
                Word {
                    hi: 1.0,
                    lo: epsilon,
                },
                Word {
                    hi: 0.25,
                    lo: -epsilon / 8.0,
                },
                Word {
                    hi: 0.125,
                    lo: epsilon / 16.0,
                },
            ],
            [
                Word {
                    hi: 0.015625,
                    lo: epsilon / 16.0,
                },
                Word {
                    hi: -0.03125,
                    lo: epsilon / 32.0,
                },
                Word {
                    hi: 0.0078125,
                    lo: -epsilon / 64.0,
                },
            ],
            [1.0, 0.25, 0.125],
        );
        let observer = moving_words(
            [
                Word {
                    hi: 0.015625,
                    lo: 2.0_f64.powi(-70),
                },
                Word {
                    hi: -0.0078125,
                    lo: -2.0_f64.powi(-72),
                },
                Word {
                    hi: 0.00390625,
                    lo: 2.0_f64.powi(-74),
                },
            ],
            [
                Word {
                    hi: 2.0_f64.powi(-16),
                    lo: 2.0_f64.powi(-80),
                },
                Word {
                    hi: -2.0_f64.powi(-17),
                    lo: 2.0_f64.powi(-82),
                },
                Word {
                    hi: 2.0_f64.powi(-18),
                    lo: -2.0_f64.powi(-84),
                },
            ],
            [0.015625, -0.0078125, 0.00390625],
        );
        let moving_delay = Tangent {
            value: Word {
                hi: 0.0078125,
                lo: 2.0_f64.powi(-80),
            },
            rate: Word {
                hi: 2.0_f64.powi(-14),
                lo: 2.0_f64.powi(-90),
            },
            scalar: 0.0078125,
        };
        let c = Tangent::fixed(C_LIGHT_AU_DAY);
        let moving_p = dot(ray, observer) / (moving_delay * c * c);
        assert_ne!(moving_p.rate.rounded(), 0.0);
        emit(
            "smooth_moving_mixed_words",
            ray,
            observer,
            moving_delay,
            C_LIGHT_AU_DAY,
            false,
        );
        // Exact radial geometry with ray = c*delay: changing radial observing
        // velocity cancels from the position law. Inputs are dyadic, not fitted.
        emit(
            "radial_stationary_mixed_words",
            moving_words(
                [
                    Word {
                        hi: 2.0,
                        lo: epsilon,
                    },
                    Word::from(0.0),
                    Word::from(0.0),
                ],
                [Word::from(0.0); 3],
                [2.0, 0.0, 0.0],
            ),
            moving_words(
                [
                    Word {
                        hi: 0.25,
                        lo: epsilon / 16.0,
                    },
                    Word::from(0.0),
                    Word::from(0.0),
                ],
                [
                    Word {
                        hi: 0.0625,
                        lo: epsilon / 64.0,
                    },
                    Word::from(0.0),
                    Word::from(0.0),
                ],
                [0.25, 0.0, 0.0],
            ),
            Tangent {
                value: Word {
                    hi: 1.0,
                    lo: epsilon / 2.0,
                },
                rate: Word::from(0.0),
                scalar: 1.0,
            },
            2.0,
            false,
        );
        // Near-floor continuations use r(t)=r(0)*(1+t/8) and delay(t)=1+t/8.
        // Their p is constant, keeping the mandated MP spacings resolvable.
        // The consumed scalar floor is held even at its discontinuous zero.
        for (id, central, floor) in [
            ("scalar_floor_inside", -2.0, true),
            (
                "scalar_floor_positive_side",
                -2.0 + 2.0_f64.powi(-51),
                false,
            ),
            (
                "scalar_floor_negative_side",
                -2.0 - 2.0_f64.powi(-51),
                false,
            ),
        ] {
            let low = 2.0_f64.powi(-80);
            emit(
                id,
                moving_words(
                    [
                        Word {
                            hi: central,
                            lo: low,
                        },
                        Word::from(0.0),
                        Word::from(0.0),
                    ],
                    [
                        Word {
                            hi: central / 8.0,
                            lo: low / 8.0,
                        },
                        Word::from(0.0),
                        Word::from(0.0),
                    ],
                    [central, 0.0, 0.0],
                ),
                moving([0.5, 0.125, 0.0], [0.0; 3]),
                Tangent {
                    value: Word::from(1.0),
                    rate: Word::from(0.125),
                    scalar: 1.0,
                },
                1.0,
                floor,
            );
        }
    }
}

#[cfg(test)]
mod shared_ppn_inverse_tests {
    use super::*;
    use std::cell::RefCell;

    fn bits(value: f64) -> String {
        format!("{:016x}", value.to_bits())
    }

    fn encode(value: Tangent) -> serde_json::Value {
        serde_json::json!({
            "value_bits": [bits(value.value.hi), bits(value.value.lo)],
            "rate_bits": [bits(value.rate.hi), bits(value.rate.lo)],
            "scalar_bits": bits(value.scalar),
            "value_rounded_bits": bits(value.value.rounded()),
            "rate_rounded_bits": bits(value.rate.rounded()),
        })
    }

    #[derive(Clone, Copy)]
    struct LinearDeflector {
        body: i32,
        position: [Word; 3],
        velocity: [Word; 3],
        scalar_position: [f64; 3],
        scalar_velocity: [f64; 3],
    }

    fn lookup(
        sources: &[LinearDeflector],
        calls: &RefCell<Vec<serde_json::Value>>,
        body: i32,
        tt: f64,
    ) -> Result<Option<State>> {
        let state = sources
            .iter()
            .find(|source| source.body == body)
            .map(|source| {
                let position = std::array::from_fn(|i| {
                    source.position[i] + source.velocity[i] * Word::from(tt)
                });
                let scalar = std::array::from_fn(|i| {
                    source.scalar_position[i] + source.scalar_velocity[i] * tt
                });
                (scalar, source.scalar_velocity, position, source.velocity)
            });
        calls.borrow_mut().push(serde_json::json!({
            "body": body,
            "tt_bits": bits(tt),
            "state": state.map(|(position, velocity, p, v)| serde_json::json!({
                "position_scalar_bits": position.map(bits),
                "velocity_scalar_bits": velocity.map(bits),
                "position_words": p.map(|x| [bits(x.hi), bits(x.lo)]),
                "velocity_words": v.map(|x| [bits(x.hi), bits(x.lo)]),
            })),
        }));
        Ok(state)
    }

    // Accepted serial PPN arithmetic, with test-only branch observations.
    // This is a rounding/control-flow comparison, not a physical oracle.
    fn serial_deflection_with_branches(
        mut ray: Motion,
        observer: Motion,
        epoch: f64,
        delay: Tangent,
        lookup: &dyn Fn(i32, f64) -> Result<Option<State>>,
        scalar_dot: fn(Motion, Motion) -> Tangent,
        branches: &RefCell<Vec<serde_json::Value>>,
    ) -> Result<Motion> {
        let norm = |a| scalar_dot(a, a).sqrt();
        let dot = scalar_dot;
        let distance = norm(ray);
        if distance.scalar == 0.0 {
            return Ok(ray);
        }
        let one = Tangent::fixed(1.0);
        let c = Tangent::fixed(C_LIGHT_AU_DAY);
        for (body, mass_ratio) in DEFLECTORS {
            let branch_index = branches.borrow().len();
            branches.borrow_mut().push(serde_json::json!({
                "body": body,
                "reception_available": false,
            }));
            let Some((position, _, position_words, velocity_words)) = lookup(body, epoch)? else {
                continue;
            };
            branches.borrow_mut()[branch_index]["reception_available"] = true.into();
            let at_reception = moving_words(position_words, velocity_words, position);
            let observer_to_body = std::array::from_fn(|i| at_reception[i] - observer[i]);
            let unit = ray.map(|x| x / distance);
            let projected_delay = dot(unit, observer_to_body) / c;
            let closest = if projected_delay.scalar < 0.0 {
                Tangent::fixed(0.0)
            } else if projected_delay.scalar > delay.scalar {
                delay
            } else {
                projected_delay
            };
            branches.borrow_mut()[branch_index]["closest_branch"] = if projected_delay.scalar < 0.0
            {
                "zero"
            } else if projected_delay.scalar > delay.scalar {
                "delay"
            } else {
                "projected"
            }
            .into();
            branches.borrow_mut()[branch_index]["projected_delay"] = encode(projected_delay);
            branches.borrow_mut()[branch_index]["closest"] = encode(closest);
            branches.borrow_mut()[branch_index]["closest_available"] = false.into();
            let Some((position, _, position_words, velocity_words)) =
                lookup(body, epoch - closest.scalar)?
            else {
                continue;
            };
            branches.borrow_mut()[branch_index]["closest_available"] = true.into();
            let at_closest = moving_words(
                position_words,
                velocity_words.map(|v| v * (Word::from(1.0) - closest.rate)),
                position,
            );
            let body_to_observer = std::array::from_fn(|i| observer[i] - at_closest[i]);
            let body_to_target = std::array::from_fn(|i| ray[i] + body_to_observer[i]);
            let source_range = norm(body_to_target);
            let observer_range = norm(body_to_observer);
            branches.borrow_mut()[branch_index]["range_skipped"] = (source_range.scalar == 0.0
                || observer_range.scalar < DEFLECTION_EMAG_FLOOR_AU)
                .into();
            if source_range.scalar == 0.0 || observer_range.scalar < DEFLECTION_EMAG_FLOOR_AU {
                continue;
            }
            let source_direction = body_to_target.map(|x| x / source_range);
            let observer_direction = body_to_observer.map(|x| x / observer_range);
            let p_dot_q = dot(unit, source_direction);
            let q_dot_e = dot(source_direction, observer_direction);
            let e_dot_p = dot(observer_direction, unit);
            branches.borrow_mut()[branch_index]["angle_skipped"] =
                (e_dot_p.scalar.abs() > DEFLECTION_EDOTP_LIMIT).into();
            if e_dot_p.scalar.abs() > DEFLECTION_EDOTP_LIMIT {
                continue;
            }
            let strength = Tangent::fixed(2.0 * GRAV_GS)
                / (Tangent::fixed(SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S)
                    * observer_range
                    * Tangent::fixed(AU_IN_M)
                    * Tangent::fixed(mass_ratio));
            let squared = observer_range * observer_range;
            let limit = Tangent::fixed(DEFLECTION_DLIM_SCALE)
                / if squared.scalar > 1.0 { squared } else { one };
            branches.borrow_mut()[branch_index]["squared_above_one"] =
                (squared.scalar > 1.0).into();
            let mut denominator = one + q_dot_e;
            branches.borrow_mut()[branch_index]["denominator_limited"] =
                (denominator.scalar < limit.scalar).into();
            if denominator.scalar < limit.scalar {
                denominator = limit;
            }
            let coefficient = strength / denominator * distance;
            ray = std::array::from_fn(|i| {
                ray[i]
                    + coefficient
                        * (p_dot_q * observer_direction[i] - e_dot_p * source_direction[i])
            });
        }
        Ok(ray)
    }

    #[test]
    fn ppn_normalization_retains_exact_radial_cancellation_and_transverse_low_words() {
        let epsilon = 2.0_f64.powi(-54);
        let vector = moving_words(
            [Word::from(2.0), Word::from(0.0), Word::from(0.0)],
            [
                Word {
                    hi: 0.125,
                    lo: epsilon,
                },
                Word {
                    hi: 0.25,
                    lo: epsilon / 8.0,
                },
                Word {
                    hi: -0.5,
                    lo: epsilon / 16.0,
                },
            ],
            [2.0, 0.0, 0.0],
        );
        let length = norm(vector);
        let inverse = Word::from(1.0) / length.value;
        let normalized = vector.map(|value| value.divide_with_inverse(length, inverse));
        // |(2,0,0)| is exactly two. Longitudinal motion cancels; the two
        // transverse rates retain the exact dyadic low words after /2.
        let expected_rates = [
            Word::from(0.0),
            Word {
                hi: 0.125,
                lo: epsilon / 16.0,
            },
            Word {
                hi: -0.25,
                lo: epsilon / 32.0,
            },
        ];
        for i in 0..3 {
            assert_eq!(
                normalized[i].value.rounded(),
                if i == 0 { 1.0 } else { 0.0 }
            );
            assert_eq!(
                normalized[i].rate.hi.to_bits(),
                expected_rates[i].hi.to_bits()
            );
            assert_eq!(
                normalized[i].rate.lo.to_bits(),
                expected_rates[i].lo.to_bits()
            );
            assert_eq!(
                normalized[i].scalar.to_bits(),
                (vector[i].scalar / length.scalar).to_bits()
            );
        }
    }

    #[test]
    fn ppn_second_acquisition_failure_preserves_serial_call_order() {
        let calls = RefCell::new(Vec::new());
        let fail_second = |body, tt: f64| -> Result<Option<State>> {
            calls.borrow_mut().push((body, tt.to_bits()));
            if calls.borrow().len() == 2 {
                return Err(Error::ConfigurationMissing {
                    missing_config: "ppn_second_lookup",
                    suggestion: "test acquisition must fail at the original closest lookup",
                });
            }
            Ok(Some((
                [1.0, 0.25, 0.125],
                [0.0; 3],
                [1.0, 0.25, 0.125].map(Word::from),
                [Word::from(0.0); 3],
            )))
        };
        let ray = moving([2.0, 0.0, 0.0], [0.0, 0.015625, 0.0]);
        let observer = moving([0.0; 3], [0.0; 3]);
        let delay = Tangent::fixed(0.0625);
        let Err(candidate) = deflection(ray, observer, 0.0, delay, &fail_second) else {
            panic!("the second lookup must fail");
        };
        let candidate_calls = calls.replace(Vec::new());
        let branches = RefCell::new(Vec::new());
        let Err(serial) = serial_deflection_with_branches(
            ray,
            observer,
            0.0,
            delay,
            &fail_second,
            dot,
            &branches,
        ) else {
            panic!("the second lookup must fail");
        };
        assert_eq!(candidate.to_string(), serial.to_string());
        assert_eq!(candidate_calls, *calls.borrow());
        assert_eq!(candidate_calls.len(), 2);
        assert_eq!(candidate_calls[0], (SUN, 0.0_f64.to_bits()));
        assert_eq!(candidate_calls[1].0, SUN);
        assert_ne!(candidate_calls[1].1, 0.0_f64.to_bits());
    }

    #[test]
    fn actual_moving_ppn_branch_emits_fixed_multiprecision_witnesses() {
        let epsilon = 2.0_f64.powi(-54);
        let make_source = |body, position: [f64; 3], velocity: [f64; 3]| LinearDeflector {
            body,
            position: std::array::from_fn(|i| Word {
                hi: position[i],
                lo: if position[i] == 0.0 {
                    0.0
                } else {
                    epsilon * position[i] / 8.0
                },
            }),
            velocity: std::array::from_fn(|i| Word {
                hi: velocity[i],
                lo: if velocity[i] == 0.0 {
                    0.0
                } else {
                    epsilon * velocity[i] / 16.0
                },
            }),
            scalar_position: position,
            scalar_velocity: velocity,
        };
        let emit = |id: &str,
                    ray: Motion,
                    observer: Motion,
                    delay: Tangent,
                    sources: &[LinearDeflector],
                    expected_sun_branch: &str,
                    limited: bool| {
            let candidate_calls = RefCell::new(Vec::new());
            let serial_calls = RefCell::new(Vec::new());
            let branches = RefCell::new(Vec::new());
            let candidate = deflection(ray, observer, 0.0, delay, &|body, tt| {
                lookup(sources, &candidate_calls, body, tt)
            })
            .unwrap();
            let serial = serial_deflection_with_branches(
                ray,
                observer,
                0.0,
                delay,
                &|body, tt| lookup(sources, &serial_calls, body, tt),
                dot,
                &branches,
            )
            .unwrap();
            assert_eq!(*candidate_calls.borrow(), *serial_calls.borrow());
            assert_eq!(branches.borrow()[0]["closest_branch"], expected_sun_branch);
            assert_eq!(branches.borrow()[0]["range_skipped"], false);
            assert_eq!(branches.borrow()[0]["angle_skipped"], false);
            assert_eq!(branches.borrow()[0]["denominator_limited"], limited);
            for i in 0..3 {
                assert_eq!(candidate[i].scalar.to_bits(), serial[i].scalar.to_bits());
                for output in [candidate[i], serial[i]] {
                    assert!(output.value.rounded().is_finite());
                    assert!(output.rate.rounded().is_finite());
                }
            }
            println!(
                "SHARED_WORD_PPN_BRANCH {}",
                serde_json::json!({
                    "schema": 1,
                    "case": id,
                    "law": "finite_source_ppn_linear_sources_consumed_scalar_branches",
                    "epoch_bits": bits(0.0),
                    "constants": {
                        "light_speed_au_day_bits": bits(C_LIGHT_AU_DAY),
                        "grav_gs_bits": bits(GRAV_GS),
                        "light_speed_m_s_bits": bits(SPEED_OF_LIGHT_M_S),
                        "light_speed_m_s_squared_bits": bits(SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S),
                        "au_m_bits": bits(AU_IN_M),
                        "observer_range_floor_bits": bits(DEFLECTION_EMAG_FLOOR_AU),
                        "edotp_limit_bits": bits(DEFLECTION_EDOTP_LIMIT),
                        "dlim_scale_bits": bits(DEFLECTION_DLIM_SCALE),
                        "deflectors": DEFLECTORS.map(|(body, ratio)| serde_json::json!({
                            "body": body, "mass_ratio_bits": bits(ratio),
                        })),
                    },
                    "ray": ray.map(encode),
                    "observer": observer.map(encode),
                    "delay": encode(delay),
                    "source_laws": sources.iter().map(|source| serde_json::json!({
                        "body": source.body,
                        "position_words": source.position.map(|x| [bits(x.hi), bits(x.lo)]),
                        "velocity_words": source.velocity.map(|x| [bits(x.hi), bits(x.lo)]),
                        "position_scalar_bits": source.scalar_position.map(bits),
                        "velocity_scalar_bits": source.scalar_velocity.map(bits),
                    })).collect::<Vec<_>>(),
                    "acquisitions": &*candidate_calls.borrow(),
                    "serial_branches": &*branches.borrow(),
                    "candidate": candidate.map(encode),
                    "serial": serial.map(encode),
                })
            );
        };
        let ray = moving_words(
            [
                Word {
                    hi: 2.0,
                    lo: epsilon,
                },
                Word::from(0.0),
                Word::from(0.0),
            ],
            [
                Word {
                    hi: 0.015625,
                    lo: epsilon / 16.0,
                },
                Word {
                    hi: -0.0078125,
                    lo: epsilon / 32.0,
                },
                Word {
                    hi: 0.00390625,
                    lo: -epsilon / 64.0,
                },
            ],
            [2.0, 0.0, 0.0],
        );
        let observer = moving_words(
            [Word::from(0.0); 3],
            [
                Word {
                    hi: 0.000244140625,
                    lo: 2.0_f64.powi(-80),
                },
                Word {
                    hi: -0.0001220703125,
                    lo: -2.0_f64.powi(-82),
                },
                Word {
                    hi: 0.00006103515625,
                    lo: 2.0_f64.powi(-84),
                },
            ],
            [0.0; 3],
        );
        let delay = Tangent {
            value: Word {
                hi: 0.0625,
                lo: 2.0_f64.powi(-80),
            },
            rate: Word {
                hi: 0.0001220703125,
                lo: 2.0_f64.powi(-90),
            },
            scalar: 0.0625,
        };
        emit(
            "moving_three_deflectors_projected",
            ray,
            observer,
            delay,
            &[
                make_source(
                    SUN,
                    [1.0, 0.25, 0.125],
                    [0.0009765625, 0.00048828125, -0.000244140625],
                ),
                make_source(
                    5,
                    [4.0, 1.0, -0.5],
                    [-0.000244140625, 0.0001220703125, 0.00006103515625],
                ),
                make_source(
                    6,
                    [8.0, -2.0, 1.0],
                    [0.0001220703125, -0.00006103515625, 0.000030517578125],
                ),
            ],
            "projected",
            false,
        );
        emit(
            "moving_closest_zero",
            ray,
            observer,
            delay,
            &[make_source(
                SUN,
                [-1.0, 0.25, 0.125],
                [0.0009765625, 0.00048828125, -0.000244140625],
            )],
            "zero",
            false,
        );
        emit(
            "moving_closest_delay",
            ray,
            observer,
            Tangent {
                value: Word::from(0.0009765625),
                rate: delay.rate,
                scalar: 0.0009765625,
            },
            &[make_source(
                SUN,
                [4.0, 0.25, 0.125],
                [0.0009765625, 0.00048828125, -0.000244140625],
            )],
            "delay",
            false,
        );
        let conjunction_ray = moving_words(
            [
                Word {
                    hi: 2.0,
                    lo: epsilon,
                },
                Word::from(2.0_f64.powi(-16)),
                Word::from(0.0),
            ],
            ray.map(|x| x.rate),
            [2.0, 2.0_f64.powi(-16), 0.0],
        );
        emit(
            "moving_conjunction_limited",
            conjunction_ray,
            observer,
            delay,
            &[make_source(SUN, [1.0, 0.0, 0.0], [0.0009765625, 0.0, 0.0])],
            "projected",
            true,
        );
    }
}
