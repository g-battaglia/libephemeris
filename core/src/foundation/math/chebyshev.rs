//! Chebyshev polynomials and their derivatives on ephemeris segments.
//!
//! For `f(x) = sum c[k] T_k(x)`, Clenshaw's backward recurrence evaluates
//! the value. First-derivative coefficients satisfy
//! `d[k] = d[k+2] + 2(k+1)c[k+1]`, with the constant term halved.
//! [`chebyshev_value_and_derivative`] streams that recurrence with constant
//! auxiliary space; [`chebyshev_deriv_coeffs`] also exposes the explicit row.
//! Higher derivatives use the differentiated polynomial basis recurrence.
//!
//! The argument `x` is dimensionless. Ephemeris readers enforce their owned
//! segment coverage; the raw finite polynomial kernel also accepts arguments
//! outside `[-1,1]`. Segment-time derivatives multiply by `(dx/dt)^k` through
//! [`scale_derivatives_to_segment`] and [`segment_derivatives`]. No coordinate
//! frame or astronomical clock is selected here.
//!
//! Empty coefficient slices, nonfinite admitted inputs, invalid derivative
//! orders and unrepresentable segment-time scales return typed errors. Very
//! large finite polynomial arguments can overflow ordinary binary64 arithmetic.
//! Constant values avoid an unnecessary overflowing `2*x` factor.

use crate::errors::{Error, Result};
use crate::types::Segment;

/// Maximum derivative order supported by [`chebyshev_dx_derivatives`].
pub const MAX_DERIVATIVE_ORDER: usize = 4;

/// Evaluate `sum_k coeffs[k] * T_k(x)` with the Clenshaw recurrence.
///
/// With `b[n+1] = b[n+2] = 0`, evaluate backward as
/// `b[k] = c[k] + 2*x*b[k+1] - b[k+2]` and return `c[0] + x*b[1] - b[2]`.
/// All arithmetic is IEEE binary64.
///
/// Returns an error for an empty slice, a non-finite `x`, or any
/// non-finite coefficient (such coefficients previously produced NaN
/// against the documented rejection policy); `x` outside
/// `[-1, 1]` is accepted (pure polynomial evaluation) and must be range
/// checked by the caller when segment semantics apply.
pub fn chebyshev_value(coeffs: &[f64], x: f64) -> Result<f64> {
    validate_series_inputs(coeffs, x)?;
    chebyshev_value_ieee(coeffs, x)
}

/// Check input admission before either checked polynomial evaluation route.
fn validate_series_inputs(coeffs: &[f64], x: f64) -> Result<()> {
    if coeffs.is_empty() {
        return Err(Error::EmptyCoefficients);
    }
    if !x.is_finite() {
        return Err(Error::NonFiniteArgument { value: x });
    }
    if let Some(&bad) = coeffs.iter().find(|c| !c.is_finite()) {
        return Err(Error::NonFiniteArgument { value: bad });
    }
    Ok(())
}

/// Evaluate the same Clenshaw recurrence with IEEE operands retained.
///
/// This arithmetic kernel propagates nonfinite operands through Clenshaw's
/// recurrence. Its caller owns range, coefficient and metadata admission.
/// An empty slice returns an error before indexing. Use [`chebyshev_value`]
/// when finite-input checks belong to this call.
pub fn chebyshev_value_ieee(coeffs: &[f64], x: f64) -> Result<f64> {
    if coeffs.is_empty() {
        return Err(Error::EmptyCoefficients);
    }
    let n = coeffs.len() - 1;
    if n == 0 {
        return Ok(coeffs[0]);
    }
    // Release order: hoisted `two_x`, loop over k = n..=1, tail outside.
    let two_x = 2.0 * x;
    let mut b1 = 0.0_f64; // b_{k+1}
    let mut b2 = 0.0_f64; // b_{k+2}
    for k in (1..=n).rev() {
        let b = coeffs[k] + two_x * b1 - b2;
        b2 = b1;
        b1 = b;
    }
    Ok(coeffs[0] + x * b1 - b2)
}

/// Derivative-series coefficients in release `_deriv_coeffs` order.
///
/// Given `f(x) = sum_{k=0}^{n} c_k T_k(x)`, returns
/// `(d_0,..., d_{n-1})` with `f'(x) = sum d_k T_k(x)` via
/// `d_{n-1} = 2*n*c_n`, `d_k = d_{k+2} + 2*(k+1)*c_{k+1}` for
/// `k = n-2..=1`, and `d_0 = d_2/2 + c_1` (with `d_2 = 0` when
/// `n < 3`). Degree-0 input yields the single zero coefficient;
/// degree-1 input yields `(c_1,)`. The caller evaluates the returned
/// row with [`chebyshev_value`], which also enforces the finite
/// coefficient policy; a non-finite input coefficient therefore
/// surfaces there, never as a silent NaN here.
pub fn chebyshev_deriv_coeffs(coeffs: &[f64]) -> Vec<f64> {
    let length = coeffs.len().saturating_sub(1).max(1);
    let mut derivative = vec![0.0; length];
    chebyshev_deriv_coeffs_into(coeffs, &mut derivative);
    derivative
}

/// Write the same derivative coefficients into caller-owned storage.
/// The output must contain `max(coeffs.len() - 1, 1)` words. Every word is
/// assigned; input coefficient validation remains with the subsequent value
/// evaluator, including the unused constant coefficient and signed zeros.
pub(crate) fn chebyshev_deriv_coeffs_into(coeffs: &[f64], d: &mut [f64]) {
    let n = coeffs.len().saturating_sub(1);
    assert_eq!(d.len(), n.max(1));
    if n == 0 {
        d[0] = 0.0;
        return;
    }
    if n == 1 {
        d[0] = coeffs[1];
        return;
    }
    d[n - 1] = (2.0 * n as f64) * coeffs[n];
    for k in (1..n - 1).rev() {
        let d_k2 = if k + 2 < n { d[k + 2] } else { 0.0 };
        d[k] = d_k2 + (2.0 * (k + 1) as f64) * coeffs[k + 1];
    }
    // d[0] special case (T_0 coefficient halving).
    let d_2 = if n >= 3 { d[2] } else { 0.0 };
    d[0] = d_2 / 2.0 + coeffs[1];
}

/// Evaluate a polynomial and its first x-derivative without allocating a row.
///
/// For `f(x) = sum c[k] T_k(x)`, derivative coefficients satisfy
/// `d[k] = d[k+2] + 2(k+1)c[k+1]`, with the constant term halved.
/// Streaming this recurrence into Clenshaw retains two derivative coefficients
/// and two accumulators: auxiliary space is constant in polynomial degree.
/// Value and derivative updates occupy independent lanes with ordinary unfused
/// multiply/add/subtract operations. Inputs are dimensionless; empty/nonfinite
/// inputs and overflowing derivative
/// coefficients return the same typed failures as the separate checked kernels.
/// Segment callers apply their own time scale and coverage validation.
pub fn chebyshev_value_and_derivative(coeffs: &[f64], x: f64) -> Result<(f64, f64)> {
    use wide::f64x2;

    validate_series_inputs(coeffs, x)?;
    let n = coeffs.len() - 1;
    if n == 0 {
        return Ok((chebyshev_value_ieee(coeffs, x)?, 0.0));
    }
    if n == 1 {
        return Ok((chebyshev_value_ieee(coeffs, x)?, coeffs[1]));
    }
    let two_x = 2.0 * x;
    // The value stream starts one coefficient before the derivative stream.
    // Keep its zero products: finite x can still make two_x overflow.
    let leading_value = coeffs[n] + two_x * 0.0 - 0.0;
    let factor = f64x2::new([two_x, two_x]);
    let (mut next, mut next2) = (0.0, 0.0);
    let (mut b1, mut b2) = (f64x2::new([leading_value, 0.0]), f64x2::new([0.0, 0.0]));
    for k in (1..n).rev() {
        let d = if k == n - 1 {
            (2.0 * n as f64) * coeffs[n]
        } else {
            next2 + (2.0 * (k + 1) as f64) * coeffs[k + 1]
        };
        if !d.is_finite() {
            return Err(Error::NonFiniteArgument { value: d });
        }
        next2 = next;
        next = d;
        // Ordinary multiply, add, subtract in each lane; no fused contraction.
        let b = f64x2::new([coeffs[k], d]) + factor * b1 - b2;
        b2 = b1;
        b1 = b;
    }
    let d0 = next2 / 2.0 + coeffs[1];
    if !d0.is_finite() {
        return Err(Error::NonFiniteArgument { value: d0 });
    }
    let result = f64x2::new([coeffs[0], d0]) + f64x2::new([x, x]) * b1 - b2;
    let [value, derivative] = result.to_array();
    Ok((value, derivative))
}

/// Evaluate a series and its first two derivatives with constant auxiliary space.
///
/// Returns `[f, df/dx, d²f/dx²]` for dimensionless `x`. Values and first
/// derivatives retain [`chebyshev_value_and_derivative`]'s arithmetic. The
/// second derivative streams the derivative-coefficient recurrence twice into
/// Clenshaw, producing the same words as differentiating two explicit rows.
/// No row or basis matrix is allocated. Ephemeris callers enforce record
/// coverage and multiply derivatives by their own segment-time scales.
/// Empty/nonfinite inputs and overflowing derivative coefficients fail as in
/// the first-derivative kernel; large finite extrapolation can overflow values.
pub fn chebyshev_value_and_two_derivatives(coeffs: &[f64], x: f64) -> Result<[f64; 3]> {
    let (value, first) = chebyshev_value_and_derivative(coeffs, x)?;
    let n = coeffs.len() - 1;
    if n < 2 {
        return Ok([value, first, 0.0]);
    }
    if n == 2 {
        let second = 4.0 * coeffs[2];
        if !second.is_finite() {
            return Err(Error::NonFiniteArgument { value: second });
        }
        return Ok([value, first, second]);
    }
    let two_x = 2.0 * x;
    let (mut d_next, mut d_next2) = (0.0, 0.0);
    let (mut e_next, mut e_next2) = (0.0, 0.0);
    let (mut b1, mut b2) = (0.0, 0.0);
    for k in (1..n).rev() {
        let d = if k == n - 1 {
            (2.0 * n as f64) * coeffs[n]
        } else {
            d_next2 + (2.0 * (k + 1) as f64) * coeffs[k + 1]
        };
        d_next2 = d_next;
        d_next = d;
        if k == 1 {
            break;
        }
        let e = if k == n - 1 {
            (2.0 * k as f64) * d
        } else {
            e_next2 + (2.0 * k as f64) * d
        };
        if !e.is_finite() {
            return Err(Error::NonFiniteArgument { value: e });
        }
        e_next2 = e_next;
        e_next = e;
        let b = e + two_x * b1 - b2;
        b2 = b1;
        b1 = b;
    }
    let e0 = e_next2 / 2.0 + d_next;
    if !e0.is_finite() {
        return Err(Error::NonFiniteArgument { value: e0 });
    }
    Ok([value, first, e0 + x * b1 - b2])
}

/// Evaluate a Chebyshev series and its first `order` derivatives.
///
/// Returns `d[0..=order]`, where `d[k]` is the kth derivative with respect
/// to normalized coordinate `x`. Differentiated basis recurrences avoid
/// finite differences. Physical time scaling is applied separately by
/// [`scale_derivatives_to_segment`] or [`segment_derivatives`].
pub fn chebyshev_dx_derivatives(coeffs: &[f64], x: f64, order: usize) -> Result<Vec<f64>> {
    if coeffs.is_empty() {
        return Err(Error::EmptyCoefficients);
    }
    if !x.is_finite() {
        return Err(Error::NonFiniteArgument { value: x });
    }
    if let Some(&bad) = coeffs.iter().find(|c| !c.is_finite()) {
        return Err(Error::NonFiniteArgument { value: bad });
    }
    if order > MAX_DERIVATIVE_ORDER {
        return Err(Error::InvalidDerivativeOrder {
            order,
            max: MAX_DERIVATIVE_ORDER,
        });
    }
    let mut out = vec![0.0; order + 1];
    // Only T[k-1] and T[k-2] and their first four derivatives are needed.
    // Accumulating each degree immediately keeps the same sum order while
    // eliminating storage proportional to polynomial degree.
    let mut older = [0.0; MAX_DERIVATIVE_ORDER + 1];
    let mut previous = [0.0; MAX_DERIVATIVE_ORDER + 1];
    for (k, &c) in coeffs.iter().enumerate() {
        let mut current = [0.0; MAX_DERIVATIVE_ORDER + 1];
        match k {
            0 => current[0] = 1.0,
            1 => {
                current[0] = x;
                current[1] = 1.0;
            }
            _ => {
                current[0] = (2.0 * previous[0]) * x - older[0];
                for j in 1..=order {
                    current[j] =
                        2.0 * j as f64 * previous[j - 1] + (2.0 * previous[j]) * x - older[j];
                }
            }
        }
        for j in 0..=order {
            out[j] += c * current[j];
        }
        older = previous;
        previous = current;
    }
    Ok(out)
}

/// Scale `x`-derivatives to segment-time derivatives: `d[k]/dt^k =
/// d[k]/dx^k * (dx/dt)^k` for a linear segment map.
///
/// A non-finite scale propagates through the products; callers that must
/// reject unrepresentable time derivatives use [`segment_derivatives`],
/// which errors with [`Error::DerivativeScaleOverflow`] instead.
pub fn scale_derivatives_to_segment(dx: &[f64], segment: &Segment) -> Vec<f64> {
    let scale = segment.dx_dt();
    let mut power = 1.0;
    dx.iter()
        .map(|&v| {
            let s = v * power;
            power *= scale;
            s
        })
        .collect()
}

/// Evaluate value and time derivatives at epoch `t` on `segment`.
///
/// Returns `d[0..=order]` with `d[k] = d^k/dt^k`; range/seam handling comes
/// from [`Segment::normalize`]. Time derivatives of order >= 1 need a finite
/// `dx/dt`; on degenerate-scale segments this errors instead of
/// producing infinities or NaN.
pub fn segment_derivatives(
    coeffs: &[f64],
    segment: &Segment,
    t: f64,
    order: usize,
) -> Result<Vec<f64>> {
    let x = segment.normalize(t)?.value();
    let dx = chebyshev_dx_derivatives(coeffs, x, order)?;
    let scale = segment.dx_dt();
    if order > 0 && !scale.is_finite() {
        return Err(Error::DerivativeScaleOverflow { scale });
    }
    Ok(scale_derivatives_to_segment(&dx, segment))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::{NormalizedTime, Segment};

    const TOL: f64 = 1e-12;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= TOL * b.abs().max(1.0)
    }

    #[test]
    fn empty_and_non_finite_inputs_error() {
        assert!(matches!(
            chebyshev_value(&[], 0.0),
            Err(Error::EmptyCoefficients)
        ));
        assert!(matches!(
            chebyshev_value(&[1.0], f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        // FIX-03: non-finite coefficients are rejected like any other
        // non-finite input instead of yielding NaN.
        assert!(matches!(
            chebyshev_value(&[1.0, f64::NAN], 0.0),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            chebyshev_dx_derivatives(&[f64::INFINITY], 0.0, 1),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            chebyshev_dx_derivatives(&[], 0.0, 1),
            Err(Error::EmptyCoefficients)
        ));
        assert!(matches!(
            chebyshev_dx_derivatives(&[1.0], 0.0, MAX_DERIVATIVE_ORDER + 1),
            Err(Error::InvalidDerivativeOrder { .. })
        ));
    }

    #[test]
    fn zero_and_constant_polynomials() {
        let zero = [0.0, 0.0, 0.0];
        for x in [-1.0, -0.3, 0.0, 0.7, 1.0] {
            assert_eq!(chebyshev_value(&zero, x).unwrap(), 0.0);
        }
        for x in [-1.0, 0.0, 0.42, 1.0] {
            assert_eq!(chebyshev_value(&[2.5], x).unwrap(), 2.5);
        }
        let d = chebyshev_dx_derivatives(&[2.5], 0.3, 2).unwrap();
        assert_eq!(d, vec![2.5, 0.0, 0.0]);
    }

    #[test]
    fn low_degree_identities_at_interior_point() {
        // T_2(x) = 2x^2 - 1, T_3(x) = 4x^3 - 3x.
        let x = 0.5;
        assert!(close(
            chebyshev_value(&[0.0, 0.0, 1.0], x).unwrap(),
            2.0 * x * x - 1.0
        ));
        assert!(close(
            chebyshev_value(&[0.0, 0.0, 0.0, 1.0], x).unwrap(),
            4.0 * x * x * x - 3.0 * x
        ));
        // T_3(0.5) = -1 exactly.
        assert!(close(
            chebyshev_value(&[0.0, 0.0, 0.0, 1.0], x).unwrap(),
            -1.0
        ));
    }

    #[test]
    fn both_normalized_endpoints() {
        // T_n(1) = 1, T_n(-1) = (-1)^n.
        for n in 0..6 {
            let mut c = vec![0.0; n + 1];
            c[n] = 1.0;
            assert!(
                close(chebyshev_value(&c, 1.0).unwrap(), 1.0),
                "T_{n}(1) != 1"
            );
            let expected = if n % 2 == 0 { 1.0 } else { -1.0 };
            assert!(
                close(chebyshev_value(&c, -1.0).unwrap(), expected),
                "T_{n}(-1) wrong"
            );
        }
        // Endpoint derivatives of T_3: T_3' = 12x^2 - 3 -> 9 at x = 1.
        let d = chebyshev_dx_derivatives(&[0.0, 0.0, 0.0, 1.0], 1.0, 1).unwrap();
        assert!(close(d[0], 1.0));
        assert!(close(d[1], 9.0));
    }

    #[test]
    fn derivatives_match_analytic_forms() {
        // p(x) = T_3(x): p' = 12x^2 - 3, p'' = 24x, p''' = 24, p'''' = 0.
        let c = [0.0, 0.0, 0.0, 1.0];
        let x = 0.5;
        let d = chebyshev_dx_derivatives(&c, x, 4).unwrap();
        assert!(close(d[0], -1.0));
        assert!(close(d[1], 12.0 * x * x - 3.0));
        assert!(close(d[2], 24.0 * x));
        assert!(close(d[3], 24.0));
        assert!(close(d[4], 0.0));
        // Linear combination: 1.5*T_0 + 2*T_1 + 0.5*T_2 at x = -0.25.
        let c2 = [1.5, 2.0, 0.5];
        let y = -0.25;
        let d2 = chebyshev_dx_derivatives(&c2, y, 2).unwrap();
        assert!(close(d2[0], 1.5 + 2.0 * y + 0.5 * (2.0 * y * y - 1.0)));
        assert!(close(d2[1], 2.0 + 0.5 * 4.0 * y));
        assert!(close(d2[2], 0.5 * 4.0));
    }

    #[test]
    fn deriv_coeff_rows_match_analytic_forms() {
        // Release `_deriv_coeffs` identities: T_3' = 3*T_0 + 6*T_2,
        // T_2' = 4*T_1, T_1' = 1, constant' = 0.
        assert_eq!(
            chebyshev_deriv_coeffs(&[0.0, 0.0, 0.0, 1.0]),
            vec![3.0, 0.0, 6.0]
        );
        assert_eq!(chebyshev_deriv_coeffs(&[0.0, 0.0, 1.0]), vec![0.0, 4.0]);
        assert_eq!(chebyshev_deriv_coeffs(&[0.0, 1.0]), vec![1.0]);
        assert_eq!(chebyshev_deriv_coeffs(&[2.5]), vec![0.0]);
        // Evaluated derivative rows recover the analytic slope:
        // T_3'(0.5) = 12*0.25 - 3 = 0 via the release row, like the
        // basis kernel (both constructions agree mathematically).
        let row = chebyshev_deriv_coeffs(&[0.0, 0.0, 0.0, 1.0]);
        assert!(close(chebyshev_value(&row, 0.5).unwrap(), 0.0));
        let d = chebyshev_dx_derivatives(&[0.0, 0.0, 0.0, 1.0], 0.5, 1).unwrap();
        assert!(close(d[1], 0.0));
    }

    #[test]
    fn derivative_row_storage_preserves_bits_and_checked_error_policy() {
        let nan = f64::from_bits(0x7ff8_0000_0000_1234);
        let mut cases = vec![
            vec![],
            vec![nan],
            vec![nan, -0.0],
            vec![f64::INFINITY, 0.0, -0.0],
            vec![0.0, f64::INFINITY, 1.0, f64::NEG_INFINITY],
            vec![0.0, nan, 2.0],
        ];
        for degree in [0, 1, 2, 3, 63, 64, 65, 129] {
            cases.push(
                (0..=degree)
                    .map(|k| (k as f64 - 3.0) * 2.0_f64.powi(-20))
                    .collect(),
            );
        }
        for coefficients in cases {
            let expected = chebyshev_deriv_coeffs(&coefficients);
            let mut actual = vec![nan; expected.len()];
            chebyshev_deriv_coeffs_into(&coefficients, &mut actual);
            assert_eq!(
                actual.iter().copied().map(f64::to_bits).collect::<Vec<_>>(),
                expected
                    .iter()
                    .copied()
                    .map(f64::to_bits)
                    .collect::<Vec<_>>()
            );
            for tau in [-1.0, -0.0, 0.0, 0.25, 1.0, nan, f64::INFINITY] {
                match (
                    chebyshev_value(&actual, tau),
                    chebyshev_value(&expected, tau),
                ) {
                    (Ok(a), Ok(b)) => assert_eq!(a.to_bits(), b.to_bits()),
                    (
                        Err(Error::NonFiniteArgument { value: a }),
                        Err(Error::NonFiniteArgument { value: b }),
                    ) => assert_eq!(a.to_bits(), b.to_bits()),
                    other => panic!("different derivative-row policy: {other:?}"),
                }
            }
        }
        // A constant's value never enters its derivative coefficient row.
        assert_eq!(chebyshev_deriv_coeffs(&[nan]), vec![0.0]);
        assert_eq!(
            chebyshev_deriv_coeffs(&[nan, -0.0])[0].to_bits(),
            (-0.0_f64).to_bits()
        );
    }

    #[test]
    fn non_unit_segment_scaling() {
        // Segment [10, 14]: dx/dt = 0.5; T_1(x(t)) = x, so d/dt = 0.5,
        // d2/dt2 = 0. T_2 gives d2/dt2 = 4 * 0.25 = 1.
        let seg = Segment::new(10.0, 14.0).unwrap();
        let d1 = segment_derivatives(&[0.0, 1.0], &seg, 12.0, 2).unwrap();
        assert!(close(d1[0], 0.0));
        assert!(close(d1[1], 0.5));
        assert!(close(d1[2], 0.0));
        let d2 = segment_derivatives(&[0.0, 0.0, 1.0], &seg, 12.0, 2).unwrap();
        // T_2(0) = -1, T_2'(0) = 0, T_2''(0) = 4.
        assert!(close(d2[0], -1.0));
        assert!(close(d2[1], 0.0));
        assert!(close(d2[2], 4.0 * 0.25));
        // Unit segment [0, 2]: dx/dt = 1, derivatives pass through unchanged.
        let unit = Segment::new(0.0, 2.0).unwrap();
        let du = segment_derivatives(&[0.0, 0.0, 0.0, 1.0], &unit, 1.5, 1).unwrap();
        let x = NormalizedTime::new(0.5).unwrap().value();
        assert!(close(du[0], 4.0 * x * x * x - 3.0 * x));
        assert!(close(du[1], 12.0 * x * x - 3.0));
    }

    #[test]
    fn segment_endpoints_evaluate() {
        let seg = Segment::new(0.0, 4.0).unwrap();
        let c = [1.0, 1.0]; // 1 + x
        let lo = segment_derivatives(&c, &seg, 0.0, 0).unwrap();
        let hi = segment_derivatives(&c, &seg, 4.0, 0).unwrap();
        assert!(close(lo[0], 0.0));
        assert!(close(hi[0], 2.0));
    }

    #[test]
    fn finite_domain_extremes() {
        // FIX-03 review probes: a constant polynomial is exactly 1 even at
        // f64::MAX (no spurious NaN from an overflowing 2*x factor), and a
        // near-MAX span normalizes to the mathematically correct ~0.5.
        assert_eq!(chebyshev_value(&[1.0], f64::MAX).unwrap(), 1.0);
        let d = chebyshev_dx_derivatives(&[1.0], f64::MAX, 1).unwrap();
        assert_eq!(d, vec![1.0, 0.0]);
        let seg = Segment::new(0.0, f64::MAX).unwrap();
        assert_eq!(seg.normalize(f64::MAX * 0.75).unwrap().value(), 0.5);
        // A subnormal-length segment has no finite dx/dt: order-0 values
        // still evaluate, but time derivatives are an explicit error.
        // (Rust's MIN_POSITIVE is the smallest *normal*; the true
        // subnormal floor is from_bits(1).)
        let tiny = Segment::new(0.0, f64::from_bits(1)).unwrap();
        assert!(!tiny.dx_dt().is_finite());
        assert!(segment_derivatives(&[1.0], &tiny, 0.0, 0).is_ok());
        assert!(matches!(
            segment_derivatives(&[1.0], &tiny, 0.0, 1),
            Err(Error::DerivativeScaleOverflow { .. })
        ));
    }

    fn assert_same_ieee(actual: f64, expected: f64) {
        if expected.is_nan() {
            // Arithmetic NaN payload/sign is not a portable Rust contract.
            assert!(actual.is_nan());
        } else {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
    }

    #[test]
    fn paired_value_rate_matches_independent_explicit_series_words() {
        let widths = [1, 2, 3, 4, 65, 66, 128, 257, 1025];
        let arguments = [
            -1.0,
            -0.0,
            0.0,
            0.25,
            1.0,
            f64::from_bits(1),
            -f64::from_bits(1),
            1.0 + f64::EPSILON,
            -1.0 - f64::EPSILON,
        ];
        let mut compared = 0;
        for width in widths {
            for pattern in 0..4 {
                let coefficients: Vec<_> = (0..width)
                    .map(|k| {
                        let sign = if k % 2 == 0 { 1.0 } else { -1.0 };
                        match pattern {
                            0 => sign * (1.0 + (k % 7) as f64 / 16.0) * 2.0_f64.powi(-(k * 7 % 48)),
                            1 => sign * 2.0_f64.powi(900 - (k % 16)),
                            2 => {
                                sign * if k % 2 == 0 {
                                    f64::from_bits(1)
                                } else {
                                    2.0_f64.powi(-1020)
                                }
                            }
                            _ => {
                                if k % 2 == 0 {
                                    -0.0
                                } else {
                                    0.0
                                }
                            }
                        }
                    })
                    .collect();
                let derivative = chebyshev_deriv_coeffs(&coefficients);
                for x in arguments {
                    // Separate checked value kernels and the explicit derivative
                    // row are independent of the paired recurrence under test.
                    let expected_value = chebyshev_value(&coefficients, x).unwrap();
                    let expected_rate = chebyshev_value(&derivative, x).unwrap();
                    let actual = chebyshev_value_and_derivative(&coefficients, x).unwrap();
                    assert_same_ieee(actual.0, expected_value);
                    assert_same_ieee(actual.1, expected_rate);
                    compared += 1;
                }
            }
        }
        assert_eq!(compared, widths.len() * 4 * arguments.len());

        // An ordinary rounded product is exactly one here. Fusing the final
        // multiply/add would instead retain a nonzero -epsilon² residual.
        let epsilon = f64::EPSILON;
        let coefficients = [-1.0, 1.0 - epsilon, 0.0];
        let x = 1.0 + epsilon;
        let actual = chebyshev_value_and_derivative(&coefficients, x).unwrap();
        assert_same_ieee(actual.0, chebyshev_value(&coefficients, x).unwrap());
        assert_eq!(actual.0.to_bits(), 0.0_f64.to_bits());
        assert_same_ieee(
            actual.1,
            chebyshev_value(&chebyshev_deriv_coeffs(&coefficients), x).unwrap(),
        );
    }

    #[test]
    fn paired_value_rate_preserves_admission_and_descending_overflow_errors() {
        let argument_nan = f64::from_bits(0x7ff8_0000_0000_1234);
        let coefficient_nan = f64::from_bits(0xfff8_0000_0000_5678);
        assert!(matches!(
            chebyshev_value_and_derivative(&[], argument_nan),
            Err(Error::EmptyCoefficients)
        ));
        assert!(matches!(
            chebyshev_value(&[], argument_nan),
            Err(Error::EmptyCoefficients)
        ));
        let cases = [
            (
                vec![coefficient_nan, f64::INFINITY],
                argument_nan,
                argument_nan,
            ),
            (
                vec![coefficient_nan, f64::NEG_INFINITY],
                0.0,
                coefficient_nan,
            ),
            (
                vec![0.0, f64::INFINITY, coefficient_nan],
                0.0,
                f64::INFINITY,
            ),
            (
                vec![0.0, f64::NEG_INFINITY, coefficient_nan],
                0.0,
                f64::NEG_INFINITY,
            ),
        ];
        for (coefficients, x, expected) in cases {
            for result in [
                chebyshev_value(&coefficients, x).map(|value| (value, 0.0)),
                chebyshev_value_and_derivative(&coefficients, x),
            ] {
                match result {
                    Err(Error::NonFiniteArgument { value }) => {
                        // These are original input payloads, not arithmetic NaNs.
                        assert_eq!(value.to_bits(), expected.to_bits());
                    }
                    other => panic!("input admission or precedence changed: {other:?}"),
                }
            }
        }

        // Highest derivative coefficient first; then a lower negative
        // overflow; finally d0 overflow after all streamed coefficients fit.
        for (coefficients, expected) in [
            ([0.0, 1.0, -f64::MAX, f64::MAX], f64::INFINITY),
            ([0.0, 1.0, -f64::MAX, 0.0], f64::NEG_INFINITY),
            ([0.0, f64::MAX, 0.0, f64::MAX / 8.0], f64::INFINITY),
        ] {
            assert!(chebyshev_value(&coefficients, 0.0).is_ok());
            match chebyshev_value_and_derivative(&coefficients, 0.0) {
                Err(Error::NonFiniteArgument { value }) => {
                    assert_eq!(value.to_bits(), expected.to_bits());
                }
                other => panic!("derivative overflow order changed: {other:?}"),
            }
        }
    }

    #[test]
    fn paired_value_rate_keeps_low_degree_signs_and_huge_argument_results() {
        for coefficient in [0.0, -0.0, 1.0, -1.0, f64::MAX] {
            for x in [f64::MAX, -f64::MAX, 0.0, -0.0, f64::from_bits(1)] {
                let (value, rate) = chebyshev_value_and_derivative(&[coefficient], x).unwrap();
                assert_eq!(value.to_bits(), coefficient.to_bits());
                assert_eq!(rate.to_bits(), 0.0_f64.to_bits());
            }
        }
        for slope in [-0.0, 0.0, 1.0, f64::MAX] {
            for x in [f64::MAX, -f64::MAX, -0.0, 0.0, 2.0, -2.0] {
                let coefficients = [-0.0, slope];
                let (value, rate) = chebyshev_value_and_derivative(&coefficients, x).unwrap();
                assert_same_ieee(value, chebyshev_value(&coefficients, x).unwrap());
                assert_eq!(rate.to_bits(), slope.to_bits());
            }
        }
        for coefficients in [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 1.0],
            [f64::MAX / 16.0, 0.0, f64::MAX / 32.0],
        ] {
            for x in [f64::MAX, -f64::MAX] {
                let actual = chebyshev_value_and_derivative(&coefficients, x).unwrap();
                assert_same_ieee(actual.0, chebyshev_value(&coefficients, x).unwrap());
                assert_same_ieee(
                    actual.1,
                    chebyshev_value(&chebyshev_deriv_coeffs(&coefficients), x).unwrap(),
                );
                // Preserve successful IEEE propagation, including the leading
                // infinity-times-zero products; do not add a finite-result gate.
                assert!(actual.0.is_nan() && actual.1.is_nan());
            }
        }
    }
}
