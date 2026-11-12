// SPDX-License-Identifier: AGPL-3.0
//! Algebraic witnesses for the concrete two-word Kepler compensation.

use super::Words;

#[test]
fn increments_and_product_residuals_survive_cancellation() {
    let small = 2.0_f64.powi(-60);
    let one = Words::new(1.0);
    let positive = one.plus(Words::new(small));
    assert_eq!(positive.value(), 1.0);
    assert_eq!(positive.minus(one).value(), small);
    let negative = one.minus(Words::new(small));
    let product = positive.times(negative);
    assert_eq!(product.minus(one).value(), -small * small);
}

#[test]
fn quotient_and_square_root_retain_their_physical_residuals() {
    let limit = 32.0 * f64::EPSILON * f64::EPSILON;
    for value in [0.125, 0.9, 2.0, 17.0, 1e100] {
        let original = Words::new(value);
        let root = original.square_root();
        assert!(root.times(root).minus(original).magnitude() <= limit * value);
        let denominator = Words::new(3.0).plus(Words::new(2.0_f64.powi(-60)));
        let quotient = original.divided(denominator);
        assert!(quotient.times(denominator).minus(original).magnitude() <= limit * value);
    }
}

#[test]
fn zero_and_sign_are_explicit_and_unrepresentable_arithmetic_is_not_success() {
    let zero = Words::new(0.0);
    assert_eq!(zero.sign(), 0.0);
    assert_eq!(zero.square_root(), zero);
    assert_eq!(Words::new(-0.25).sign(), -1.0);
    assert!(!Words::new(1e308).times(Words::new(1e308)).finite());
    assert!(!Words::new(1.0).divided(zero).finite());
    assert!(!Words::new(-1.0).square_root().finite());
}
