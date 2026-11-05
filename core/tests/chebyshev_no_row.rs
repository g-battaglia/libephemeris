//! The allocation-free derivative retains the existing coefficient arithmetic.
use libephemeris_core::math::chebyshev::{
    chebyshev_deriv_coeffs, chebyshev_dx_derivatives, chebyshev_value,
    chebyshev_value_and_derivative, chebyshev_value_and_two_derivatives,
};

#[test]
fn rolling_higher_derivatives_preserve_the_full_basis_arithmetic() {
    for degree in [0, 1, 2, 3, 8, 32, 128] {
        let coefficients: Vec<_> = (0..=degree).map(|k| (k as f64 - 3.0) / 1024.0).collect();
        for x in [-1.0, -0.75, -0.0, 0.0, 0.125, 0.75, 1.0] {
            for order in 0..=4 {
                let mut basis = vec![vec![0.0; order + 1]; degree + 1];
                basis[0][0] = 1.0;
                if degree > 0 {
                    basis[1][0] = x;
                    if order > 0 {
                        basis[1][1] = 1.0;
                    }
                }
                for k in 2..=degree {
                    basis[k][0] = (2.0 * basis[k - 1][0]) * x - basis[k - 2][0];
                    for j in 1..=order {
                        basis[k][j] = 2.0 * j as f64 * basis[k - 1][j - 1]
                            + (2.0 * basis[k - 1][j]) * x
                            - basis[k - 2][j];
                    }
                }
                let mut expected = vec![0.0; order + 1];
                for (k, coefficient) in coefficients.iter().enumerate() {
                    for j in 0..=order {
                        expected[j] += coefficient * basis[k][j];
                    }
                }
                let actual = chebyshev_dx_derivatives(&coefficients, x, order).unwrap();
                assert_eq!(
                    actual.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    expected.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
                );
            }
        }
    }
}

#[test]
fn no_row_derivative_preserves_binary64_words() {
    let mut seed = 123456789_u64;
    for degree in 0..=128 {
        let coefficients: Vec<_> = (0..=degree)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((seed >> 11) as f64 / (1_u64 << 53) as f64 - 0.5) * 0.003
            })
            .collect();
        let row = chebyshev_deriv_coeffs(&coefficients);
        for i in -128..=128 {
            let x = i as f64 / 128.0;
            let (value, rate) = chebyshev_value_and_derivative(&coefficients, x).unwrap();
            assert_eq!(
                value.to_bits(),
                chebyshev_value(&coefficients, x).unwrap().to_bits()
            );
            assert_eq!(
                rate.to_bits(),
                chebyshev_value(&row, x).unwrap().to_bits(),
                "degree {degree}, x={x}"
            );
        }
    }
    for coefficients in [vec![0.0], vec![-0.0], vec![0.0, -0.0], vec![1.0, -0.0, 0.0]] {
        for x in [-1.0, -0.0, 0.0, 1.0] {
            let (_, rate) = chebyshev_value_and_derivative(&coefficients, x).unwrap();
            assert_eq!(
                rate.to_bits(),
                chebyshev_value(&chebyshev_deriv_coeffs(&coefficients), x)
                    .unwrap()
                    .to_bits()
            );
        }
    }
}

#[test]
fn second_derivative_preserves_two_explicit_coefficient_rows() {
    let mut seed = 987654321_u64;
    for degree in 0..=128 {
        let coefficients: Vec<_> = (0..=degree)
            .map(|_| {
                seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                ((seed >> 11) as f64 / (1_u64 << 53) as f64 - 0.5) * 0.003
            })
            .collect();
        let first = chebyshev_deriv_coeffs(&coefficients);
        let second = chebyshev_deriv_coeffs(&first);
        for i in -128..=128 {
            let x = i as f64 / 128.0;
            let values = chebyshev_value_and_two_derivatives(&coefficients, x).unwrap();
            for (actual, row) in values.into_iter().zip([&coefficients, &first, &second]) {
                assert_eq!(
                    actual.to_bits(),
                    chebyshev_value(row, x).unwrap().to_bits(),
                    "degree={degree}, x={x}"
                );
            }
        }
    }
}

#[test]
fn second_derivative_matches_an_exact_cubic_and_segment_scaling() {
    for x in [-1.0_f64, -0.5, 0.0, 0.5, 1.0] {
        let [value, first, second] =
            chebyshev_value_and_two_derivatives(&[2.0, -3.0, 4.0, 5.0], x).unwrap();
        assert_eq!(value, 20.0 * x * x * x + 8.0 * x * x - 18.0 * x - 2.0);
        assert_eq!(first, 60.0 * x * x + 16.0 * x - 18.0);
        assert_eq!(second, 120.0 * x + 16.0);
    }
    // A record spanning four TT days has dx/dTT=1/2. Its acceleration
    // therefore scales the second normalized derivative by 1/4.
    assert_eq!(
        chebyshev_value_and_two_derivatives(&[0.0, 0.0, 1.0], 0.25).unwrap()[2] / 4.0,
        1.0
    );
}

#[test]
fn second_derivative_retains_signed_zero_and_checked_failure_contracts() {
    for coefficients in [vec![-0.0], vec![1.0, -0.0], vec![1.0, 0.0, -0.0]] {
        let expected = chebyshev_deriv_coeffs(&chebyshev_deriv_coeffs(&coefficients));
        for x in [-1.25, -0.0, 0.0, 1.25] {
            assert_eq!(
                chebyshev_value_and_two_derivatives(&coefficients, x).unwrap()[2].to_bits(),
                chebyshev_value(&expected, x).unwrap().to_bits()
            );
        }
    }
    for (coefficients, x) in [
        (vec![], 0.0),
        (vec![1.0], f64::NAN),
        (vec![f64::INFINITY], 0.0),
        // First-derivative coefficients fit; the second row overflows.
        (vec![0.0, 0.0, 0.0, f64::MAX / 8.0], 0.0),
    ] {
        assert!(chebyshev_value_and_two_derivatives(&coefficients, x).is_err());
    }
}

#[test]
fn checked_inputs_are_still_rejected() {
    for (coefficients, x) in [
        (vec![], 0.0),
        (vec![1.0], f64::NAN),
        (vec![f64::INFINITY], 0.0),
        (vec![1.0, 2.0, f64::MAX], 0.0),
    ] {
        assert!(chebyshev_value_and_derivative(&coefficients, x).is_err());
    }
}

#[test]
fn streamed_derivative_satisfies_low_degree_polynomial_identity() {
    // T3 = 4x^3 - 3x; c = [2,-3,4,5] gives f = 20x^3+8x^2-18x-2.
    for x in [-1.0_f64, -0.5, 0.0, 0.5, 1.0] {
        let (value, rate) = chebyshev_value_and_derivative(&[2.0, -3.0, 4.0, 5.0], x).unwrap();
        assert_eq!(value, 20.0 * x * x * x + 8.0 * x * x - 18.0 * x - 2.0);
        assert_eq!(rate, 60.0 * x * x + 16.0 * x - 18.0);
    }
}
