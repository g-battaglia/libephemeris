// SPDX-License-Identifier: AGPL-3.0
//! Exact-flow, cancellation and physical-invariant regression witnesses.

use super::super::conic::{Elements, initial_state};
use super::*;

fn circle() -> State {
    State {
        position: [1.0, 0.0, 0.0],
        velocity: [0.0, 1.0, 0.0],
    }
}
fn close(actual: [f64; 3], expected: [f64; 3]) {
    let scale = norm(expected).unwrap();
    for axis in 0..3 {
        assert!(
            (actual[axis] - expected[axis]).abs() <= 512.0 * f64::EPSILON * scale,
            "axis {axis}: {actual:?} != {expected:?}"
        );
    }
}
fn energy(state: State) -> f64 {
    dot(state.velocity, state.velocity) / 2.0 - 1.0 / norm(state.position).unwrap()
}
fn angular(state: State) -> [f64; 3] {
    let p = state.position;
    let v = state.velocity;
    [
        p[1] * v[2] - p[2] * v[1],
        p[2] * v[0] - p[0] * v[2],
        p[0] * v[1] - p[1] * v[0],
    ]
}

#[test]
fn universal_circular_flow_retains_quarter_period_cosines_and_both_directions() {
    for interval in [
        0.0_f64,
        0.1,
        -0.2,
        std::f64::consts::FRAC_PI_2,
        std::f64::consts::TAU,
        10.0,
        100_000.0,
    ] {
        let actual = advance(circle(), 1.0, interval).unwrap();
        close(actual.position, [interval.cos(), interval.sin(), 0.0]);
        close(actual.velocity, [-interval.sin(), interval.cos(), 0.0]);
    }
    let quarter = advance(circle(), 1.0, std::f64::consts::FRAC_PI_2).unwrap();
    assert!(quarter.position[0] > 0.0);
    assert!(quarter.velocity[1] > 0.0);
}

#[test]
fn ellipse_and_hyperbola_retain_energy_angular_momentum_and_reverse_flow() {
    for (a, e) in [(1.7, 0.2), (-2.0, 1.3)] {
        let initial = initial_state(
            Elements {
                semimajor_axis: a,
                eccentricity: e,
                inclination: 0.4,
                ascending_node: 0.7,
                periapsis_argument: 0.3,
                mean_anomaly: 0.2,
            },
            1.0,
        )
        .unwrap();
        let evolved = advance(initial, 1.0, 3.1).unwrap();
        let reverse = advance(evolved, 1.0, -3.1).unwrap();
        close(reverse.position, initial.position);
        close(reverse.velocity, initial.velocity);
        close(angular(evolved), angular(initial));
        let reference = energy(initial);
        assert!((energy(evolved) - reference).abs() <= 512.0 * f64::EPSILON * reference.abs());
    }
}

#[test]
fn parabolic_cartesian_flow_has_the_barker_limit_without_an_eccentricity_switch() {
    let initial = State {
        position: [1.0, 0.0, 0.0],
        velocity: [0.0, 2.0_f64.sqrt(), 0.0],
    };
    let d = 0.5_f64;
    let interval = 2.0_f64.sqrt() * (d + d * d * d / 3.0);
    let actual = advance(initial, 1.0, interval).unwrap();
    close(actual.position, [1.0 - d * d, 2.0 * d, 0.0]);
    close(
        actual.velocity,
        [
            -2.0_f64.sqrt() * d / (1.0 + d * d),
            2.0_f64.sqrt() / (1.0 + d * d),
            0.0,
        ],
    );
}

#[test]
fn universal_functions_satisfy_their_defining_recurrences_through_zero_energy() {
    for b in [-0.1_f64, 0.0, 0.1, 1.0] {
        for x in [0.0_f64, 0.1, 1.0, 3.0] {
            let values = universal(Words::new(x), Words::new(b)).unwrap();
            assert!(
                (values.zero.value() + b * values.two.value() - 1.0).abs() <= 64.0 * f64::EPSILON
            );
            assert!(
                (values.one.value() + b * values.three.value() - x).abs()
                    <= 64.0 * f64::EPSILON * x.max(1.0)
            );
        }
    }
    let zero = universal(Words::new(2.0), Words::new(0.0)).unwrap();
    assert_eq!(
        [zero.zero, zero.one, zero.two, zero.three].map(Words::value),
        [1.0, 2.0, 2.0, 4.0 / 3.0]
    );
}

#[test]
fn normalized_time_does_not_erase_representable_transverse_motion() {
    // These are scaling identities, not astronomical force fixtures.
    // In each case the tiny dimensionless phase still gives nonzero y=v*dt.
    for (radius, interval) in [(1e200_f64, 1.0_f64), (1e308, 1.0), (1e300, 1e300)] {
        let speed = 1.0 / radius.sqrt();
        let state = State {
            position: [radius, 0.0, 0.0],
            velocity: [0.0, speed, 0.0],
        };
        let actual = advance(state, 1.0, interval).unwrap();
        assert!(actual.position[1] > 0.0);
        let expected = speed * interval;
        assert!((actual.position[1] - expected).abs() <= 64.0 * f64::EPSILON * expected);
    }
}

#[test]
fn binding_energy_keeps_raw_cartesian_product_residuals_near_the_parabolic_boundary() {
    // Exact binary64 input words; references independently evaluated with
    // Decimal100 from b=2-r|v|²/GM, not inferred from a trajectory residual.
    for (speed, expected) in [
        (94_906_265.624_251_54, 3.546_046_371_670_331e-16_f64),
        (94_906_265.624_251_56, -2.734_323_463_064_769_3e-16),
    ] {
        for axis in 0..3 {
            let mut position = [0.0; 3];
            let mut velocity = [0.0; 3];
            position[axis] = 2.0_f64.powi(-52);
            velocity[(axis + 1) % 3] = speed;
            let state = State { position, velocity };
            let actual =
                normalized_energy(state, 1.0, vector_norm(position.map(Words::new)).unwrap())
                    .unwrap()
                    .value();
            assert!(
                (actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs(),
                "{actual} != {expected}"
            );
        }
    }
}

#[test]
fn compensated_root_and_endpoint_preserve_tiny_components_and_periapsis_return() {
    // Decimal120/160 osculating-conic solutions of the exact input PV words.
    // These previously failed despite a small rounded time residual. Keep the
    // unchanged component/exact-zero budgets, including components below1ulp
    // of an ordinary coordinate. No nominal element time is a reference.
    for (initial, interval, expected) in [
        (
            State {
                position: [0.9950041652780258, 0.09983341664682815, 0.0],
                velocity: [-0.09983341664682815, 0.9950041652780258, 0.0],
            },
            -0.1,
            [
                1.0_f64,
                -8.503152428150948e-18,
                0.0,
                -7.801237259850591e-18,
                1.0,
                0.0,
            ],
        ),
        (
            State {
                position: [-104.35520818950674, 2.220094986138211e-6, 0.0],
                velocity: [-1.0095371765387784, 2.127536361833131e-8, 0.0],
            },
            -100.0,
            [
                -1.4551449814069393e-10,
                -3.5950387233893144e-13,
                0.0,
                117235.96724259364,
                144.81965862460234,
                0.0,
            ],
        ),
    ] {
        let result = advance(initial, 1.0, interval).unwrap();
        for (axis, (actual, reference)) in result
            .position
            .into_iter()
            .chain(result.velocity)
            .zip(expected)
            .enumerate()
        {
            let limit = reference.abs() * if axis < 3 { 1e-8 } else { 1e-6 };
            assert!(
                (actual - reference).abs() <= limit,
                "axis{axis}: {actual} != {reference}"
            );
        }
    }
}

#[test]
fn zero_time_is_exact_and_invalid_physical_or_nonfinite_states_fail() {
    assert_eq!(advance(circle(), 1.0, 0.0).unwrap(), circle());
    for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert_eq!(advance(circle(), 1.0, bad), Err(NumericalError::NonFinite));
        assert_eq!(advance(circle(), bad, 1.0), Err(NumericalError::NonFinite));
    }
    for gm in [0.0, -1.0] {
        assert_eq!(
            advance(circle(), gm, 1.0),
            Err(NumericalError::InvalidMass { body: 10 })
        );
    }
    let origin = State {
        position: [0.0; 3],
        ..circle()
    };
    assert_eq!(
        advance(origin, 1.0, 1.0),
        Err(NumericalError::SingularForce { body: 10 })
    );
    let radial = State {
        velocity: [1.0, 0.0, 0.0],
        ..circle()
    };
    assert!(matches!(
        advance(radial, 1.0, 1.0),
        Err(NumericalError::InvalidOrbit(_))
    ));
}
