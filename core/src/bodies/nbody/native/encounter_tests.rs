// SPDX-License-Identifier: AGPL-3.0
//! Geometric sampling identities and a weak moving-source Kepler invariant.
use super::super::newtonian::{PointMass, acceleration};
use super::super::radau::{ForceProvider, IntegrationError, Radau, Settings};
use super::*;

fn source(position: [f64; 3], velocity: [f64; 3]) -> Source {
    Source {
        point: PointMass {
            body: 42,
            position,
            gravitational_parameter: 1.0,
        },
        velocity,
        post_newtonian_acceleration: None,
    }
}

fn target() -> State {
    State {
        position: [0.0; 3],
        velocity: [0.0; 3],
    }
}

#[test]
fn crossing_time_is_translation_boost_rotation_and_direction_invariant() {
    let base = source([3.0, 4.0, 0.0], [0.0, 0.0, 2.0]);
    assert_eq!(maximum_step(target(), &[base], None), Ok(2.5));
    for direction in [-1.0, 1.0] {
        let transform = |v: [f64; 3]| [v[2], -v[0], v[1]];
        let origin = [128.0, -32.0, 16.0];
        let boost = [4.0, -2.0, 8.0];
        let state = State {
            position: origin,
            velocity: boost,
        };
        let mut shifted = base;
        shifted.point.position =
            std::array::from_fn(|i| origin[i] + transform(base.point.position)[i]);
        shifted.velocity =
            std::array::from_fn(|i| boost[i] + direction * transform(base.velocity)[i]);
        assert_eq!(maximum_step(state, &[shifted], None), Ok(2.5));
    }
}

#[test]
fn only_explicit_self_and_zero_mass_are_excluded_and_stationary_motion_is_unlimited() {
    let mut self_source = source([0.0; 3], [1.0; 3]);
    assert!(matches!(
        maximum_step(target(), &[self_source], None),
        Err(NumericalError::SingularForce { body: 42 })
    ));
    assert_eq!(
        maximum_step(target(), &[self_source], Some(42)),
        Ok(f64::INFINITY)
    );
    self_source.point.gravitational_parameter = 0.0;
    assert_eq!(
        maximum_step(target(), &[self_source], None),
        Ok(f64::INFINITY)
    );
    assert_eq!(
        maximum_step(target(), &[source([1.0; 3], [0.0; 3])], None),
        Ok(f64::INFINITY)
    );
    let distant = source([10.0, 0.0, 0.0], [1.0, 0.0, 0.0]);
    let mut near = source([1.0, 0.0, 0.0], [2.0, 0.0, 0.0]);
    near.point.body = 43;
    assert_eq!(maximum_step(target(), &[distant, near], None), Ok(0.5));
    assert_eq!(maximum_step(target(), &[near, distant], Some(43)), Ok(10.0));
}

#[test]
fn invalid_excluded_sources_are_not_hidden_by_the_step_guard() {
    let mut row = source([1.0; 3], [0.0; 3]);
    row.velocity[0] = f64::NAN;
    assert_eq!(
        maximum_step(target(), &[row], Some(42)),
        Err(NumericalError::NonFinite)
    );
    row.velocity = [0.0; 3];
    row.point.gravitational_parameter = -1.0;
    assert_eq!(
        maximum_step(target(), &[row], Some(42)),
        Err(NumericalError::InvalidMass { body: 42 })
    );
}

struct Flyby {
    direction: f64,
    guarded: bool,
}
impl Flyby {
    fn source(&self, time: f64) -> Source {
        let mut row = source(
            [
                -1.0 + self.direction * 0.01 * time + 0.5e-4 * time * time,
                1e-4,
                0.0,
            ],
            [self.direction * 0.01 + 1e-4 * time, 0.0, 0.0],
        );
        row.point.gravitational_parameter = 1e-15;
        row
    }
}
impl ForceProvider for Flyby {
    type Error = NumericalError;
    fn acceleration(&mut self, time: f64, state: &State) -> Result<[f64; 3], NumericalError> {
        let mut value = acceleration(state.position, &[self.source(time).point], None)?;
        value[0] += 1e-4;
        Ok(value)
    }
    fn maximum_step(&mut self, time: f64, state: &State) -> Result<f64, NumericalError> {
        if self.guarded {
            maximum_step(*state, &[self.source(time)], None)
        } else {
            Ok(f64::INFINITY)
        }
    }
}

#[test]
fn weak_fast_flyby_in_a_uniform_acceleration_field_resolves_the_exact_relative_kepler_orbit() {
    // In the uniformly accelerating source frame the common acceleration
    // cancels exactly: r''=-GM r/|r|^3. This is an independent conic invariant,
    // not a trajectory obtained by tightening this integrator's settings.
    for direction in [-1.0, 1.0] {
        let endpoint = direction * 200.0;
        let relative = State {
            position: [1.0, -1e-4, 0.0],
            velocity: [-direction * 0.01, 0.0, 0.0],
        };
        let conic = super::super::kepler::advance(relative, 1e-15, endpoint).unwrap();
        let prescribed = Flyby {
            direction,
            guarded: true,
        }
        .source(endpoint);
        let expected = State {
            position: std::array::from_fn(|i| conic.position[i] + prescribed.point.position[i]),
            velocity: std::array::from_fn(|i| conic.velocity[i] + prescribed.velocity[i]),
        };
        let mut errors = Vec::new();
        for guarded in [false, true] {
            let mut solver = Radau::new(
                target(),
                0.0,
                Settings {
                    coefficient_tolerance: 1e-9,
                    initial_step: 0.001,
                    maximum_step: 32.0,
                    maximum_attempts: 100_000,
                    maximum_iterations: 64,
                },
            )
            .unwrap();
            let statistics = solver
                .advance_with_provider(endpoint, &mut Flyby { direction, guarded })
                .unwrap();
            let actual = solver.state();
            let error = actual
                .position
                .into_iter()
                .zip(expected.position)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0_f64, f64::max);
            eprintln!(
                "direction={direction} guard={guarded} position error={error:e}, {statistics:?}"
            );
            errors.push(error);
            if guarded {
                for (a, b) in actual.velocity.into_iter().zip(expected.velocity) {
                    assert!((a - b).abs() <= 1e-13, "velocity {a:e} vs {b:e}");
                }
            }
        }
        assert!(errors[1] <= 1e-12, "guarded position error: {}", errors[1]);
        assert!(
            errors[0] > 1e-10,
            "unguarded positive control: {}",
            errors[0]
        );
    }
}

#[test]
fn sampling_failure_propagates_once_before_any_force_evaluation_or_state_change() {
    struct Failed {
        calls: usize,
    }
    impl ForceProvider for Failed {
        type Error = &'static str;
        fn acceleration(&mut self, _: f64, _: &State) -> Result<[f64; 3], Self::Error> {
            panic!("failed sampling must not call acceleration")
        }
        fn maximum_step(&mut self, _: f64, _: &State) -> Result<f64, Self::Error> {
            self.calls += 1;
            Err("original record I/O error")
        }
    }
    let mut solver = Radau::new(
        target(),
        0.0,
        Settings {
            coefficient_tolerance: 1e-9,
            initial_step: 0.001,
            maximum_step: 32.0,
            maximum_attempts: 100,
            maximum_iterations: 64,
        },
    )
    .unwrap();
    let mut provider = Failed { calls: 0 };
    assert!(matches!(
        solver.advance_with_provider(1.0, &mut provider),
        Err(IntegrationError::Force("original record I/O error"))
    ));
    assert_eq!(provider.calls, 1);
    assert_eq!(solver.state(), target());
    assert_eq!(solver.time(), 0.0);
}

#[test]
fn malformed_sampling_limits_are_explicit_and_positive_infinity_is_unconstrained() {
    struct Limited(f64);
    impl ForceProvider for Limited {
        type Error = NumericalError;
        fn acceleration(&mut self, _: f64, _: &State) -> Result<[f64; 3], NumericalError> {
            Ok([0.0; 3])
        }
        fn maximum_step(&mut self, _: f64, _: &State) -> Result<f64, NumericalError> {
            Ok(self.0)
        }
    }
    for limit in [f64::NAN, -1.0, f64::NEG_INFINITY, 0.0, f64::INFINITY] {
        let mut solver = Radau::new(
            target(),
            0.0,
            Settings {
                coefficient_tolerance: 1e-9,
                initial_step: 0.001,
                maximum_step: 32.0,
                maximum_attempts: 100,
                maximum_iterations: 64,
            },
        )
        .unwrap();
        let result = solver.advance_with_provider(1.0, &mut Limited(limit));
        if limit == f64::INFINITY {
            assert!(result.is_ok());
            assert_eq!(solver.time(), 1.0);
        } else {
            let expected = if limit == 0.0 {
                NumericalError::StepUnderflow
            } else {
                NumericalError::InvalidSettings
            };
            assert!(matches!(result, Err(IntegrationError::Numerical(error)) if error == expected));
            assert_eq!(solver.time(), 0.0);
        }
        assert_eq!(solver.state(), target());
    }
}
