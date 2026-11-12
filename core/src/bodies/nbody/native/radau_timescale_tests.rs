// SPDX-License-Identifier: AGPL-3.0
//! Independent polynomial identities for the lower-derivative time criterion.
use super::*;

// Preserve the rejected PRS23 candidate's algebraic witness too. This helper
// is test-only, not a runtime alternate controller or a fallback selection.
fn prs23_ratio(polynomial: &[[f64; 3]; STAGES]) -> f64 {
    let [acceleration, jerk, snap, _] = derivative_norms(polynomial);
    if acceleration == 0.0 {
        return 0.0;
    }
    let phase_squared = 0.5 * ((jerk / acceleration).powi(2) + snap / acceleration);
    phase_squared.sqrt().powi(7) / 5040.0
}

#[test]
fn nested_newton_differentiation_has_the_exact_quadratic_derivatives() {
    let mut polynomial = [[0.0; 3]; STAGES];
    polynomial[0] = [2.0, 3.0, 5.0];
    polynomial[1] = [7.0, 11.0, 13.0];
    polynomial[2] = [17.0, 19.0, 23.0];
    for fraction in [0.0, 0.25, 1.0] {
        let actual = scaled_derivatives(&polynomial, fraction);
        for axis in 0..3 {
            let expected = [
                polynomial[0][axis]
                    + polynomial[1][axis] * fraction
                    + polynomial[2][axis] * fraction * (fraction - NODES[1]),
                polynomial[1][axis] + polynomial[2][axis] * (2.0 * fraction - NODES[1]),
                2.0 * polynomial[2][axis],
            ];
            assert_eq!(actual[3][axis], 0.0);
            for derivative in 0..3 {
                assert!(
                    (actual[derivative][axis] - expected[derivative]).abs()
                        <= 32.0 * f64::EPSILON * expected[derivative].abs()
                );
            }
        }
    }
}

#[test]
fn prs23_anchor_witness_is_rotation_and_acceleration_scale_invariant() {
    // a(s)=(1,s,0): the maximum norm is at s=1, |a|=sqrt(2),
    // |da/ds|=1 and d²a/ds²=0. Hence (h/tau)^7/7! = 1/(128*5040).
    let mut polynomial = [[0.0; 3]; STAGES];
    polynomial[0] = [1.0, 0.0, 0.0];
    polynomial[1] = [0.0, 1.0, 0.0];
    let expected = 1.0 / (128.0 * 5040.0);
    for scale in [0.125, 1.0, 1024.0] {
        for rotated in [false, true] {
            let transformed = polynomial.map(|value| {
                if rotated {
                    [-value[1] * scale, 0.0, value[0] * scale]
                } else {
                    value.map(|entry| entry * scale)
                }
            });
            let actual = prs23_ratio(&transformed);
            assert!((actual - expected).abs() <= 64.0 * f64::EPSILON * expected);
        }
    }
    assert_eq!(prs23_ratio(&[[0.0; 3]; STAGES]), 0.0);
    polynomial[1] = [0.0; 3];
    assert_eq!(prs23_ratio(&polynomial), 0.0);
}

#[test]
fn a85_retains_circular_calibration_and_a_genuine_acceleration_zero() {
    let expected = 1.0 / 5040.0;
    for norms in [[1.0, 1.0, 1.0, 1.0], [0.0, 1.0, 0.0, 1.0]] {
        for scale in [0.125, 1.0, 1024.0] {
            let actual = a85_ratio(norms.map(|value| value * scale));
            assert!((actual - expected).abs() <= 32.0 * f64::EPSILON * expected);
        }
    }
    assert_eq!(a85_ratio([1.0, 0.0, 0.0, 0.0]), 0.0);
    assert_eq!(a85_ratio([0.0, 1.0, 0.0, 0.0]), 0.0);
}

#[test]
fn small_force_sample_roundoff_does_not_become_a_seventh_derivative_timescale() {
    // A small linear physical change plus explicit sample quantization. Every
    // sample still contributes to the collocation polynomial/integral. Only
    // the estimator uses lower derivatives instead of amplifying a top-degree
    // difference into a fictitious short physical timescale.
    let samples = std::array::from_fn(|stage| {
        [
            1.0 + if stage % 2 == 0 { 1e-13 } else { -1e-13 },
            1e-4 * NODES[stage],
            0.0,
        ]
    });
    let polynomial = coefficients(&samples);
    assert!(timestep_ratio(&polynomial) < 1e-9);
    assert!(scaled_derivatives(&polynomial, 1.0)[1][1] > 0.0);
}
