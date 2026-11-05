// SPDX-License-Identifier: AGPL-3.0
//! Prepared typed-force-error checks; execute in the final qualification phase.
use super::*;

fn initial() -> State {
    State {
        position: [1.0, 0.0, 0.0],
        velocity: [0.0, 1.0, 0.0],
    }
}
fn trajectory() -> Radau {
    Radau::new(
        initial(),
        0.0,
        Settings {
            coefficient_tolerance: 1e-9,
            initial_step: 0.1,
            maximum_step: 1.0,
            maximum_attempts: 100,
            maximum_iterations: 32,
        },
    )
    .unwrap()
}

#[test]
fn typed_io_failure_propagates_once_without_a_different_step_or_endpoint() {
    let mut trajectory = trajectory();
    let mut calls = 0;
    let result = trajectory.advance_with_error(1.0, &mut |_, _| {
        calls += 1;
        Err(std::io::Error::new(
            std::io::ErrorKind::UnexpectedEof,
            "synthetic coefficient record",
        ))
    });
    match result {
        Err(IntegrationError::Force(error)) => {
            assert_eq!(error.kind(), std::io::ErrorKind::UnexpectedEof)
        }
        other => panic!("original typed I/O failure was not retained: {other:?}"),
    }
    assert_eq!(calls, 1);
    assert_eq!(trajectory.time(), 0.0);
    assert_eq!(trajectory.state(), initial());
}

#[test]
fn callback_named_iteration_limit_is_not_an_implicit_solver_rejection() {
    let mut trajectory = trajectory();
    let mut calls = 0;
    let error = trajectory
        .advance_with_error(1.0, &mut |_, _| {
            calls += 1;
            Err(NumericalError::IterationLimit)
        })
        .unwrap_err();
    assert!(matches!(
        error,
        IntegrationError::Force(NumericalError::IterationLimit)
    ));
    assert_eq!(calls, 1);
    assert_eq!(trajectory.time(), 0.0);
}

#[test]
fn numerical_convenience_api_uses_the_same_nonretrying_force_route() {
    let mut trajectory = trajectory();
    let mut calls = 0;
    let result = trajectory.advance(1.0, &mut |_, _| {
        calls += 1;
        Err(NumericalError::IterationLimit)
    });
    assert_eq!(result, Err(NumericalError::IterationLimit));
    assert_eq!(calls, 1);
    assert_eq!(trajectory.state(), initial());
}

#[test]
fn a_real_solver_iteration_limit_remains_distinct_from_the_callback_error_type() {
    let mut statistics = Statistics::default();
    let result = trial_with_error(
        AccumulatedState::new(initial()).unwrap(),
        0.0,
        0.1,
        0,
        &mut |_, _| Ok::<_, std::io::Error>([0.0; 3]),
        &mut statistics,
    );
    assert!(matches!(
        result,
        Err(IntegrationError::Numerical(NumericalError::IterationLimit))
    ));
    assert_eq!(statistics.force_evaluations, 1);
}

#[test]
fn successful_callback_with_nonfinite_acceleration_is_a_numerical_failure() {
    let mut trajectory = trajectory();
    let mut calls = 0;
    let result = trajectory.advance_with_error(1.0, &mut |_, _| {
        calls += 1;
        Ok::<_, std::io::Error>([f64::NAN, 0.0, 0.0])
    });
    assert!(matches!(
        result,
        Err(IntegrationError::Numerical(NumericalError::NonFinite))
    ));
    assert_eq!(calls, 1);
    assert_eq!(trajectory.time(), 0.0);
}
