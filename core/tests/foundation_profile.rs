//! Foundation profile groundwork (the source-independent kernels part).
//!
//! Seeded deterministic sweep over the Chebyshev kernels, checking
//! mathematical identities that must hold in every build profile:
//!
//! 1. Clenshaw evaluation agrees with direct `T_n`-recurrence summation.
//!    Sanity bound `1e-12` relative to the coefficient scale, fixed before
//!    the first run, not fitted to it. This is a numerical sanity check in
//!    the sense of the project rules — agreement between two evaluation
//!    orders — not a certified rounding bound for the Clenshaw operation
//!    sequence.
//! 2. Exact low-degree polynomials (`x^3`, `x^4` in Chebyshev basis) match
//!    `powi` values and analytic derivatives.
//! 3. Analytic first derivatives match central finite differences on the
//!    exact low-degree cases with step `h = 2^-13`: rounding contributes at
//!    most ~`eps/h ~= 2e-12` in value scale and truncation `h^2/6 * max|f'''|`
//!    stays below `1e-6` for these polynomials; the asserted bound is `1e-6`.
//!
//! The profile comparison itself is the recorded dual run: this file passes
//! under both `cargo test` (debug) and `cargo test --release`. A single test
//! binary only ever sees one profile, so no in-test profile detection is
//! attempted.

use libephemeris_core::math::chebyshev::{chebyshev_dx_derivatives, chebyshev_value};

/// Deterministic xorshift64* generator (fixed seed: reproducible by
/// construction, no stored vectors needed).
struct Rng(u64);

impl Rng {
    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_F491_4F6C_DD1D)
    }

    fn next_unit(&mut self) -> f64 {
        // 53-bit mantissa mapped to [0, 1).
        ((self.next_u64() >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.next_unit()
    }
}

/// Direct evaluation: explicit `T_n` recurrence plus summation, an
/// algorithmically independent path from Clenshaw.
fn direct_value(coeffs: &[f64], x: f64) -> f64 {
    let mut sum = 0.0;
    let (mut t_prev2, mut t_prev1) = (1.0, x);
    for (n, &c) in coeffs.iter().enumerate() {
        let t = match n {
            0 => t_prev2,
            1 => t_prev1,
            _ => {
                let t = 2.0 * x * t_prev1 - t_prev2;
                t_prev2 = t_prev1;
                t_prev1 = t;
                t
            }
        };
        sum += c * t;
    }
    sum
}

#[test]
fn seeded_sweep_clenshaw_agrees_with_direct_summation() {
    let mut rng = Rng(0x1E9E5EED);
    for _ in 0..300 {
        let degree = (rng.next_u64() % 13) as usize;
        let coeffs: Vec<f64> = (0..=degree).map(|_| rng.range(-2.0, 2.0)).collect();
        let x = rng.range(-1.0, 1.0);
        let a = chebyshev_value(&coeffs, x).unwrap();
        let b = direct_value(&coeffs, x);
        let scale = coeffs.iter().map(|c| c.abs()).sum::<f64>().max(1.0);
        assert!(
            (a - b).abs() <= 1e-12 * scale,
            "Clenshaw/direct disagree: {a} vs {b} (coeffs={coeffs:?}, x={x})"
        );
    }
}

#[test]
fn exact_low_degree_values_and_derivatives() {
    // x^3 = (T_3 + 3 T_1) / 4; x^4 = (T_4 + 4 T_2 + 3 T_0) / 8.
    let cube = [0.0, 0.75, 0.0, 0.25];
    let fourth = [3.0 / 8.0, 0.0, 0.5, 0.0, 1.0 / 8.0];
    let mut rng = Rng(0xC0FFEE);
    for _ in 0..100 {
        let x = rng.range(-1.0, 1.0);
        let v3 = chebyshev_value(&cube, x).unwrap();
        assert!((v3 - x.powi(3)).abs() <= 4e-16 * x.abs().powi(3).max(1.0));
        let d3 = chebyshev_dx_derivatives(&cube, x, 2).unwrap();
        assert!((d3[1] - 3.0 * x * x).abs() <= 4e-16 * (3.0 * x * x).max(1.0));
        assert!((d3[2] - 6.0 * x).abs() <= 1e-15 * (6.0 * x.abs()).max(1.0));
        let d4 = chebyshev_dx_derivatives(&fourth, x, 4).unwrap();
        assert!((d4[0] - x.powi(4)).abs() <= 8e-16);
        assert!((d4[1] - 4.0 * x.powi(3)).abs() <= 2e-15);
        assert!((d4[4] - 24.0).abs() <= 1e-13);
    }
}

#[test]
fn analytic_derivatives_match_central_differences() {
    let fourth = [3.0 / 8.0, 0.0, 0.5, 0.0, 1.0 / 8.0];
    let h = 2f64.powi(-13);
    let mut rng = Rng(0xDECADE);
    for _ in 0..100 {
        let x = rng.range(-0.9, 0.9);
        let analytic = chebyshev_dx_derivatives(&fourth, x, 1).unwrap()[1];
        let numeric = (chebyshev_value(&fourth, x + h).unwrap()
            - chebyshev_value(&fourth, x - h).unwrap())
            / (2.0 * h);
        assert!(
            (analytic - numeric).abs() <= 1e-6 * analytic.abs().max(1.0),
            "derivative mismatch at x={x}: {analytic} vs {numeric}"
        );
    }
}
