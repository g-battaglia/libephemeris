// SPDX-License-Identifier: AGPL-3.0
//! Analytic conic identities; high-precision/reference qualification is separate.

use super::*;

#[path = "conic_reference.rs"]
mod reference;

fn elements(a: f64, eccentricity: f64, mean_anomaly: f64) -> Elements {
    Elements {
        semimajor_axis: a,
        eccentricity,
        inclination: 0.0,
        ascending_node: 0.0,
        periapsis_argument: 0.0,
        mean_anomaly,
    }
}

fn close(actual: f64, expected: f64, scale: f64) {
    let bound = 512.0 * f64::EPSILON * scale;
    assert!(
        (actual - expected).abs() <= bound,
        "actual={actual:.17e}, expected={expected:.17e}, bound={bound:.3e}"
    );
}

fn invariants(state: State, orbit: Elements, gm: f64) {
    let r = state.position[0]
        .hypot(state.position[1])
        .hypot(state.position[2]);
    let kinetic = state.velocity.into_iter().map(|v| v * v).sum::<f64>() / 2.0;
    let potential = gm / r;
    close(
        kinetic - potential,
        -gm / (2.0 * orbit.semimajor_axis),
        kinetic + potential,
    );
    let cross = [
        state.position[1] * state.velocity[2] - state.position[2] * state.velocity[1],
        state.position[2] * state.velocity[0] - state.position[0] * state.velocity[2],
        state.position[0] * state.velocity[1] - state.position[1] * state.velocity[0],
    ];
    let momentum = cross[0].hypot(cross[1]).hypot(cross[2]);
    let shape = ((1.0 - orbit.eccentricity).abs() * (1.0 + orbit.eccentricity)).sqrt();
    let expected = (gm * orbit.semimajor_axis.abs()).sqrt() * shape;
    close(momentum, expected, expected);
}

#[test]
fn exact_input_states_match_independent_decimal120_and_decimal160_references() {
    for (index, (parameters, expected)) in reference::CASES.iter().enumerate() {
        let [
            a,
            e,
            mean,
            gm,
            inclination,
            ascending_node,
            periapsis_argument,
        ] = *parameters;
        let orbit = Elements {
            semimajor_axis: a,
            eccentricity: e,
            mean_anomaly: mean,
            inclination,
            ascending_node,
            periapsis_argument,
        };
        let actual = initial_state(orbit, gm).unwrap();
        let position_scale = expected[0].hypot(expected[1]).hypot(expected[2]);
        let velocity_scale = expected[3].hypot(expected[4]).hypot(expected[5]);
        for (axis, (value, truth)) in actual
            .position
            .into_iter()
            .chain(actual.velocity)
            .zip(*expected)
            .enumerate()
        {
            let scale = if axis < 3 {
                position_scale
            } else {
                velocity_scale
            };
            let bound = 512.0 * f64::EPSILON * scale;
            assert!(
                (value - truth).abs() <= bound,
                "case={index}, component={axis}, actual={value:.17e}, reference={truth:.17e}, bound={bound:.3e}"
            );
        }
    }
}

#[test]
fn circular_orbits_reduce_to_sine_and_cosine_even_for_large_mean_anomaly() {
    for mean in [0.0_f64, 0.2, -0.2, 3.0, 1e12, -1e12] {
        let orbit = elements(2.0, 0.0, mean);
        let state = initial_state(orbit, 8.0).unwrap();
        close(state.position[0], 2.0 * mean.cos(), 2.0);
        close(state.position[1], 2.0 * mean.sin(), 2.0);
        close(state.velocity[0], -2.0 * mean.sin(), 2.0);
        close(state.velocity[1], 2.0 * mean.cos(), 2.0);
        invariants(state, orbit, 8.0);
    }
}

#[test]
fn elliptic_states_match_eccentric_anomaly_parameterization() {
    for eccentricity in [0.1_f64, 0.8, 0.99] {
        for anomaly in [0.2_f64, 1.0, 3.0] {
            let mean = anomaly - eccentricity * anomaly.sin();
            let orbit = elements(2.0, eccentricity, mean);
            let state = initial_state(orbit, 8.0).unwrap();
            let shape = (1.0 - eccentricity * eccentricity).sqrt();
            let radius = 1.0 - eccentricity * anomaly.cos();
            let expected = [
                2.0 * (anomaly.cos() - eccentricity),
                2.0 * shape * anomaly.sin(),
                -2.0 * anomaly.sin() / radius,
                2.0 * shape * anomaly.cos() / radius,
            ];
            for (actual, expected) in [
                state.position[0],
                state.position[1],
                state.velocity[0],
                state.velocity[1],
            ]
            .into_iter()
            .zip(expected)
            {
                close(actual, expected, expected.abs().max(1.0));
            }
            invariants(state, orbit, 8.0);
        }
    }
}

#[test]
fn hyperbolic_states_match_hyperbolic_anomaly_parameterization() {
    for eccentricity in [1.01_f64, 1.5, 3.0] {
        for anomaly in [0.2_f64, 1.0, 3.0] {
            let mean = eccentricity * anomaly.sinh() - anomaly;
            let orbit = elements(-2.0, eccentricity, mean);
            let state = initial_state(orbit, 8.0).unwrap();
            let shape = (eccentricity * eccentricity - 1.0).sqrt();
            let radius = eccentricity * anomaly.cosh() - 1.0;
            let expected = [
                2.0 * (eccentricity - anomaly.cosh()),
                2.0 * shape * anomaly.sinh(),
                -2.0 * anomaly.sinh() / radius,
                2.0 * shape * anomaly.cosh() / radius,
            ];
            for (actual, expected) in [
                state.position[0],
                state.position[1],
                state.velocity[0],
                state.velocity[1],
            ]
            .into_iter()
            .zip(expected)
            {
                close(actual, expected, expected.abs().max(1.0));
            }
            invariants(state, orbit, 8.0);
        }
    }
}

#[test]
fn signed_small_anomalies_and_near_parabolic_states_are_not_rounded_to_periapsis() {
    for (a, eccentricity) in [
        (1.0, 1.0 - 2.0_f64.powi(-52)),
        (-1.0, 1.0 + 2.0_f64.powi(-52)),
    ] {
        for mean in [1e-30_f64, 1e-12, 0.1] {
            let positive = initial_state(elements(a, eccentricity, mean), 1.0).unwrap();
            let negative = initial_state(elements(a, eccentricity, -mean), 1.0).unwrap();
            assert!(positive.position[1] > 0.0);
            assert!(positive.velocity[0] < 0.0);
            assert_eq!(positive.position[0], negative.position[0]);
            assert_eq!(positive.position[1], -negative.position[1]);
            assert_eq!(positive.velocity[0], -negative.velocity[0]);
            assert_eq!(positive.velocity[1], negative.velocity[1]);
            invariants(positive, elements(a, eccentricity, mean), 1.0);
        }
    }
}

#[test]
fn three_elementary_rotations_match_the_periapsis_basis_and_preserve_invariants() {
    let mut orbit = elements(2.0, 0.3, 0.5);
    let planar = initial_state(orbit, 1.0).unwrap();
    orbit.inclination = 2.0; // retrograde, not folded into a prograde orbit
    orbit.ascending_node = 0.7;
    orbit.periapsis_argument = -0.4;
    let transformed = initial_state(orbit, 1.0).unwrap();
    let rotate = |mut vector: [f64; 3]| {
        for (axis, angle) in [
            (2, orbit.periapsis_argument),
            (0, orbit.inclination),
            (2, orbit.ascending_node),
        ] {
            let (s, c) = angle.sin_cos();
            vector = if axis == 2 {
                [
                    c * vector[0] - s * vector[1],
                    s * vector[0] + c * vector[1],
                    vector[2],
                ]
            } else {
                [
                    vector[0],
                    c * vector[1] - s * vector[2],
                    s * vector[1] + c * vector[2],
                ]
            };
        }
        vector
    };
    for (actual, expected) in transformed
        .position
        .into_iter()
        .chain(transformed.velocity)
        .zip(
            rotate(planar.position)
                .into_iter()
                .chain(rotate(planar.velocity)),
        )
    {
        close(actual, expected, expected.abs().max(1.0));
    }
    invariants(transformed, orbit, 1.0);
}

#[test]
fn invalid_domains_and_nonfinite_elements_do_not_construct_a_replacement_orbit() {
    for (a, e) in [(0.0, 0.1), (1.0, -0.1), (1.0, 1.0), (1.0, 1.1), (-1.0, 0.9)] {
        assert!(matches!(
            initial_state(elements(a, e, 0.0), 1.0),
            Err(NumericalError::InvalidOrbit(_))
        ));
    }
    assert!(matches!(
        initial_state(elements(1.0, 0.1, 0.0), 0.0),
        Err(NumericalError::InvalidOrbit(_))
    ));
    for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(
            initial_state(elements(1.0, 0.1, value), 1.0),
            Err(NumericalError::NonFinite)
        );
    }
}
