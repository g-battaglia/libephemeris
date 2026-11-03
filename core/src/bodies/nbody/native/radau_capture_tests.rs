// SPDX-License-Identifier: AGPL-3.0
//! Diagnostic capture must preserve the actual numerical and force-call route.
use super::*;

fn settings() -> Settings {
    Settings {
        coefficient_tolerance: 1e-9,
        initial_step: 0.001,
        maximum_step: 1.0,
        maximum_attempts: 20_000,
        maximum_iterations: 64,
    }
}

#[test]
fn accepted_capture_preserves_words_force_counts_and_both_time_directions() {
    let initial = State {
        position: [1.0, 0.25, -0.5],
        velocity: [0.1, 0.5, 0.0],
    };
    for endpoint in [-8.0, 8.0] {
        let mut ordinary = Radau::new(initial, 0.0, settings()).unwrap();
        let mut captured = ordinary.clone();
        captured.enable_step_capture();
        let mut force = |_: f64, state: &State| -> Result<[f64; 3], NumericalError> {
            Ok(state.position.map(|coordinate| -coordinate))
        };
        let a = ordinary.advance(endpoint, &mut force).unwrap();
        let b = captured.advance(endpoint, &mut force).unwrap();
        assert_eq!(a, b);
        assert_eq!(ordinary.state(), captured.state());
        assert_eq!(ordinary.time(), captured.time());
        assert!(ordinary.step_capture().is_none());
        let (steps, complete) = captured.step_capture().unwrap();
        assert!(complete);
        assert_eq!(steps.len(), b.accepted_steps);
        assert_eq!(steps.first().unwrap().initial, initial);
        assert_eq!(steps.first().unwrap().start_time, 0.0);
        assert_eq!(steps.last().unwrap().endpoint, captured.state());
        assert_eq!(steps.last().unwrap().end_time, endpoint);
        for step in steps {
            assert_eq!(step.step.signum(), endpoint.signum());
            assert!(
                step.acceleration_coefficients
                    .iter()
                    .flatten()
                    .all(|x| x.is_finite())
            );
        }
        for pair in steps.windows(2) {
            assert_eq!(pair[0].end_time, pair[1].start_time);
            assert_eq!(pair[0].endpoint, pair[1].initial);
        }
    }
}

#[test]
fn capture_overflow_is_explicit_and_does_not_stop_or_replace_integration() {
    let initial = State {
        position: [0.0; 3],
        velocity: [1.0, 0.0, 0.0],
    };
    let mut config = settings();
    config.initial_step = 1.0;
    let mut ordinary = Radau::new(initial, 0.0, config).unwrap();
    let mut captured = ordinary.clone();
    captured.enable_step_capture();
    let endpoint = (MAX_CAPTURED_STEPS + 1) as f64;
    let mut force = |_: f64, _: &State| Ok([0.0; 3]);
    let a = ordinary.advance(endpoint, &mut force).unwrap();
    let b = captured.advance(endpoint, &mut force).unwrap();
    assert_eq!(a, b);
    assert_eq!(ordinary.state(), captured.state());
    assert_eq!(captured.time(), endpoint);
    let (steps, complete) = captured.step_capture().unwrap();
    assert!(!complete);
    assert_eq!(steps.len(), MAX_CAPTURED_STEPS);
    assert_eq!(b.accepted_steps, MAX_CAPTURED_STEPS + 1);
}
