// SPDX-License-Identifier: AGPL-3.0
//! Analytic and manufactured-solution checks independent of the legacy backend.

use super::newtonian::{PointMass, acceleration};
use super::radau::{Radau, Settings, Statistics, trial};
use super::radau_tables::{NODES, POSITION_INTEGRALS, VELOCITY_INTEGRALS};
use super::state::{AccumulatedState, NumericalError, State, Sum};

fn settings() -> Settings {
    Settings {
        coefficient_tolerance: 1e-9,
        initial_step: 0.125,
        maximum_step: 1.0,
        maximum_attempts: 100_000,
        maximum_iterations: 64,
    }
}

fn initial() -> State {
    State {
        position: [1.0, 0.0, 0.0],
        velocity: [0.0, 1.0, 0.0],
    }
}

fn roundoff_close(actual: f64, expected: f64, operations: f64) {
    let bound = operations * f64::EPSILON * expected.abs().max(1.0);
    assert!(
        (actual - expected).abs() <= bound,
        "{actual:.17e} != {expected:.17e}, bound {bound:.3e}"
    );
}

fn state_close(actual: State, expected: State, operations: f64) {
    for (a, b) in actual
        .position
        .into_iter()
        .chain(actual.velocity)
        .zip(expected.position.into_iter().chain(expected.velocity))
    {
        roundoff_close(a, b, operations);
    }
}

#[test]
fn compensated_sum_retains_both_large_and_small_addition_orders() {
    for terms in [[1e16, 1.0, -1e16], [1.0, 1e16, -1e16]] {
        let mut sum = Sum::default();
        for term in terms {
            sum.add(term);
        }
        assert_eq!(sum.value(), 1.0);
    }
    let mut state = AccumulatedState::new(initial()).unwrap();
    let increment = State {
        position: [f64::EPSILON / 8.0, 0.0, 0.0],
        velocity: [0.0; 3],
    };
    for _ in 0..8 {
        state = state.advanced(increment).unwrap();
    }
    assert_eq!(state.state().position[0], 1.0 + f64::EPSILON);
}

#[test]
fn generated_basis_has_prescribed_endpoint_and_integral_constants() {
    assert_eq!(NODES[0], 0.0);
    assert!(NODES.windows(2).all(|pair| pair[0] < pair[1]));
    const { assert!(NODES[7] < 1.0) };
    assert_eq!(VELOCITY_INTEGRALS[0], [0.0; 8]);
    assert_eq!(POSITION_INTEGRALS[0], [0.0; 8]);
    for row in 0..=8 {
        let argument = if row == 8 { 1.0 } else { NODES[row] };
        assert_eq!(VELOCITY_INTEGRALS[row][0], argument);
        roundoff_close(POSITION_INTEGRALS[row][0], argument * argument / 2.0, 2.0);
    }
}

#[test]
fn zero_duration_performs_no_force_evaluation() {
    let mut integrator = Radau::new(initial(), 0.0, settings()).unwrap();
    let mut force = |_: f64, _: &State| -> Result<[f64; 3], NumericalError> {
        panic!("zero duration must not evaluate the force")
    };
    assert_eq!(
        integrator.advance(0.0, &mut force).unwrap(),
        Statistics::default()
    );
    assert_eq!(integrator.state(), initial());
}

#[test]
fn constant_acceleration_is_integrated_without_interpolation_error() {
    let initial = State {
        position: [3.0, -2.0, 4.0],
        velocity: [0.5, -1.0, 2.0],
    };
    for endpoint in [0.5, -0.5, 3.0, -3.0] {
        let acceleration = [2.0, -0.25, 0.0];
        let mut integrator = Radau::new(initial, 0.0, settings()).unwrap();
        let mut force = |_: f64, _: &State| Ok(acceleration);
        let statistics = integrator.advance(endpoint, &mut force).unwrap();
        let expected = State {
            position: std::array::from_fn(|axis| {
                initial.position[axis]
                    + endpoint * initial.velocity[axis]
                    + 0.5 * acceleration[axis] * endpoint * endpoint
            }),
            velocity: std::array::from_fn(|axis| {
                initial.velocity[axis] + endpoint * acceleration[axis]
            }),
        };
        state_close(integrator.state(), expected, 16.0);
        assert_eq!(integrator.time(), endpoint);
        assert!(statistics.accepted_steps > 0);
        assert_eq!(statistics.rejected_steps, 0);
    }
}

#[test]
fn prescribed_polynomial_accelerations_have_exact_radau_moments() {
    // Eight-node left-Radau integrates acceleration moments up to degree 14;
    // the second integral adds a weight (1-s), hence degree <=13 for position.
    let initial = State {
        position: [0.0; 3],
        velocity: [0.0; 3],
    };
    for degree in 0..=14 {
        let mut force = |time: f64, _: &State| Ok([time.powi(degree), 0.0, 0.0]);
        let value = trial(
            AccumulatedState::new(initial).unwrap(),
            0.0,
            1.0,
            64,
            &mut force,
            &mut Statistics::default(),
        )
        .unwrap()
        .increment;
        roundoff_close(value.velocity[0], 1.0 / f64::from(degree + 1), 128.0);
        if degree <= 13 {
            roundoff_close(
                value.position[0],
                1.0 / f64::from((degree + 1) * (degree + 2)),
                128.0,
            );
        }
    }
}

#[test]
fn fixed_step_has_fifteenth_order_on_an_analytic_oscillator() {
    let mut errors = Vec::new();
    // A unit step already rounds to the analytic binary64 answer, so it
    // cannot measure order. These two coarse steps retain resolvable defects.
    for step in [4.0, 2.0] {
        let mut force = |_: f64, state: &State| Ok(state.position.map(|value| -value));
        let increment = trial(
            AccumulatedState::new(initial()).unwrap(),
            0.0,
            step,
            64,
            &mut force,
            &mut Statistics::default(),
        )
        .unwrap()
        .increment;
        let actual = AccumulatedState::new(initial())
            .unwrap()
            .advanced(increment)
            .unwrap()
            .state();
        let expected = State {
            position: [step.cos(), step.sin(), 0.0],
            velocity: [-step.sin(), step.cos(), 0.0],
        };
        let error = actual
            .position
            .into_iter()
            .chain(actual.velocity)
            .zip(expected.position.into_iter().chain(expected.velocity))
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        errors.push(error);
    }
    // Step-local error scales as h^16. A one-order margin permits mixed
    // leading terms and binary64 rounding, but rejects a low-order substitute.
    assert!(errors[0] > 0.0);
    assert!(
        errors[0] / errors[1].max(f64::EPSILON) > 2.0_f64.powi(14),
        "{errors:?}"
    );
}

#[test]
fn adaptive_oscillator_supports_retained_forward_backward_endpoints() {
    let mut integrator = Radau::new(initial(), 0.0, settings()).unwrap();
    let mut force = |_: f64, state: &State| Ok(state.position.map(|value| -value));
    for time in [0.75, 1.0, 7.0, -1.0, 0.0] {
        integrator.advance(time, &mut force).unwrap();
        let expected = State {
            position: [time.cos(), time.sin(), 0.0],
            velocity: [-time.sin(), time.cos(), 0.0],
        };
        state_close(integrator.state(), expected, 512.0);
        assert_eq!(integrator.time(), time);
    }
}

#[test]
fn velocity_dependent_damped_oscillator_matches_closed_solution() {
    // x=e^(-t/4) cos(t), v=e^(-t/4)(-sin(t)-cos(t)/4).
    // Therefore x'' = -x' / 2 - 17 x / 16, including nonconservative force.
    let initial = State {
        position: [1.0, 0.0, 0.0],
        velocity: [-0.25, 0.0, 0.0],
    };
    let mut integrator = Radau::new(initial, 0.0, settings()).unwrap();
    let mut force = |_: f64, state: &State| {
        Ok(std::array::from_fn(|axis| {
            -0.5 * state.velocity[axis] - (17.0 / 16.0) * state.position[axis]
        }))
    };
    for time in [1.0_f64, 4.0, -1.0] {
        integrator.advance(time, &mut force).unwrap();
        let decay = (-time / 4.0).exp();
        let expected = State {
            position: [decay * time.cos(), 0.0, 0.0],
            velocity: [decay * (-time.sin() - time.cos() / 4.0), 0.0, 0.0],
        };
        state_close(integrator.state(), expected, 512.0);
    }
}

#[test]
fn coefficient_controller_rejects_overlarge_steps_without_changing_models() {
    let mut options = settings();
    options.initial_step = 2.0;
    options.maximum_step = 2.0;
    let mut integrator = Radau::new(initial(), 0.0, options).unwrap();
    let mut force = |_: f64, state: &State| Ok(state.position.map(|value| -value));
    let counts = integrator.advance(2.0, &mut force).unwrap();
    assert!(counts.rejected_steps > 0);
    assert!(counts.accepted_steps > 0);
    state_close(
        integrator.state(),
        State {
            position: [2.0_f64.cos(), 2.0_f64.sin(), 0.0],
            velocity: [-2.0_f64.sin(), 2.0_f64.cos(), 0.0],
        },
        128.0,
    );
}

#[test]
fn iteration_and_attempt_bounds_refuse_unfinished_propagation() {
    let mut force = |_: f64, state: &State| Ok(state.position.map(|value| -value));
    assert_eq!(
        trial(
            AccumulatedState::new(initial()).unwrap(),
            0.0,
            1.0,
            1,
            &mut force,
            &mut Statistics::default()
        )
        .unwrap_err(),
        NumericalError::IterationLimit
    );
    let mut options = settings();
    options.maximum_attempts = 1;
    options.initial_step = 0.125;
    let mut integrator = Radau::new(initial(), 0.0, options).unwrap();
    assert_eq!(
        integrator.advance(1.0, &mut force),
        Err(NumericalError::StepLimit)
    );
    assert_eq!(integrator.time(), 0.125);
}

#[test]
fn nonfinite_inputs_settings_outputs_and_unrepresentable_steps_are_errors() {
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut bad = initial();
        bad.velocity[0] = invalid;
        assert_eq!(bad.checked(), Err(NumericalError::NonFinite));
        let mut integrator = Radau::new(initial(), 0.0, settings()).unwrap();
        assert_eq!(
            integrator.advance(invalid, &mut |_, _| Ok([0.0; 3])),
            Err(NumericalError::NonFinite)
        );
        assert_eq!(
            integrator.advance(1.0, &mut |_, _| Ok([invalid; 3])),
            Err(NumericalError::NonFinite)
        );
        assert_eq!(integrator.state(), initial());
    }
    let mut options = settings();
    options.coefficient_tolerance = 0.0;
    assert!(matches!(
        Radau::new(initial(), 0.0, options),
        Err(NumericalError::InvalidSettings)
    ));
    assert_eq!(
        trial(
            AccumulatedState::new(initial()).unwrap(),
            0.0,
            0.0,
            64,
            &mut |_, _| Ok([0.0; 3]),
            &mut Statistics::default()
        )
        .unwrap_err(),
        NumericalError::InvalidStep
    );
    let mut options = settings();
    options.initial_step = f64::MIN_POSITIVE;
    let mut integrator = Radau::new(initial(), 1.0, options).unwrap();
    assert_eq!(
        integrator.advance(2.0, &mut |_, _| Ok([0.0; 3])),
        Err(NumericalError::StepUnderflow)
    );
}

#[test]
fn force_error_propagates_without_retries_or_analytic_success() {
    let mut calls = 0;
    let mut integrator = Radau::new(initial(), 0.0, settings()).unwrap();
    let error = NumericalError::SingularForce { body: 10 };
    let result = integrator.advance(1.0, &mut |_, _| {
        calls += 1;
        Err(error)
    });
    assert_eq!(result, Err(error));
    assert_eq!(calls, 1);
    assert_eq!(integrator.time(), 0.0);
    assert_eq!(integrator.state(), initial());
}

#[test]
fn newtonian_force_has_inverse_square_scaling_and_rotation_covariance() {
    let sun = PointMass {
        body: 10,
        position: [0.0; 3],
        gravitational_parameter: 1.0,
    };
    assert_eq!(
        acceleration([2.0, 0.0, 0.0], &[sun], None).unwrap(),
        [-0.25, 0.0, 0.0]
    );
    assert_eq!(
        acceleration([0.0, -2.0, 0.0], &[sun], None).unwrap(),
        [0.0, 0.25, 0.0]
    );
    for radius in [1e-100_f64, 1e100] {
        let result = acceleration([radius, 0.0, 0.0], &[sun], None).unwrap();
        roundoff_close(result[0] * radius * radius, -1.0, 8.0);
    }
}

#[test]
fn self_exclusion_preserves_every_other_configured_mass() {
    let sources = [
        PointMass {
            body: 1,
            position: [0.0; 3],
            gravitational_parameter: 5.0,
        },
        PointMass {
            body: 2,
            position: [2.0, 0.0, 0.0],
            gravitational_parameter: 4.0,
        },
        PointMass {
            body: 3,
            position: [0.0, -3.0, 0.0],
            gravitational_parameter: 9.0,
        },
    ];
    assert_eq!(
        acceleration([0.0; 3], &sources, Some(1)).unwrap(),
        [1.0, -1.0, 0.0]
    );
    assert_eq!(
        acceleration([0.0; 3], &sources, None),
        Err(NumericalError::SingularForce { body: 1 })
    );
    let duplicate = [sources[1], sources[1]];
    assert_eq!(
        acceleration([0.0; 3], &duplicate, None),
        Err(NumericalError::DuplicateSource { body: 2 })
    );
    let invalid = [PointMass {
        gravitational_parameter: -1.0,
        ..sources[1]
    }];
    assert_eq!(
        acceleration([0.0; 3], &invalid, Some(2)),
        Err(NumericalError::InvalidMass { body: 2 })
    );
}

#[test]
fn circular_point_mass_orbit_preserves_analytic_state_energy_and_momentum() {
    let source = PointMass {
        body: 10,
        position: [0.0; 3],
        gravitational_parameter: 1.0,
    };
    let mut force = |_: f64, state: &State| acceleration(state.position, &[source], None);
    let mut integrator = Radau::new(initial(), 0.0, settings()).unwrap();
    for time in [1.0, std::f64::consts::TAU, -std::f64::consts::TAU, 0.0] {
        integrator.advance(time, &mut force).unwrap();
        let state = integrator.state();
        state_close(
            state,
            State {
                position: [time.cos(), time.sin(), 0.0],
                velocity: [-time.sin(), time.cos(), 0.0],
            },
            4096.0,
        );
        let radius = state.position[0]
            .hypot(state.position[1])
            .hypot(state.position[2]);
        let speed2: f64 = state.velocity.into_iter().map(|v| v * v).sum();
        roundoff_close(speed2 / 2.0 - 1.0 / radius, -0.5, 4096.0);
        roundoff_close(
            state.position[0] * state.velocity[1] - state.position[1] * state.velocity[0],
            1.0,
            4096.0,
        );
    }
}

#[test]
fn independent_threads_share_no_mutable_solver_state() {
    let evaluate = |endpoint: f64| {
        let mut integrator = Radau::new(initial(), 0.0, settings()).unwrap();
        let mut force = |_: f64, state: &State| Ok(state.position.map(|value| -value));
        integrator.advance(endpoint, &mut force).unwrap();
        integrator.state()
    };
    let endpoints = [1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 0.5, -0.5];
    let sequential = endpoints.map(evaluate);
    let threads = endpoints.map(|endpoint| std::thread::spawn(move || evaluate(endpoint)));
    for (thread, expected) in threads.into_iter().zip(sequential) {
        assert_eq!(thread.join().unwrap(), expected);
    }
}
