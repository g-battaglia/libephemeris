// SPDX-License-Identifier: AGPL-3.0
//! Two-term binary64 compensation for cancellation-sensitive Kepler parameters.
//!
//! This concrete, private expansion carries the low word through Sundman time,
//! the root residual and f/g cancellation. There is no alternate operand type,
//! arbitrary precision or backend dispatch. Addition uses the existing Neumaier
//! sum; multiplication recovers a rounded product's residual with IEEE FMA.
//! Division and square root correct their binary64 leading approximation using
//! those same residuals. Unrepresentable operations remain nonfinite so the
//! physical caller can return its explicit numerical failure.

use super::super::state::Sum;

#[derive(Clone, Copy, Debug, PartialEq)]
pub(super) struct Words {
    high: f64,
    low: f64,
}

impl Words {
    pub(super) const fn new(high: f64) -> Self {
        Self { high, low: 0.0 }
    }

    fn accumulated(terms: &[f64]) -> Self {
        let mut sum = Sum::default();
        for &term in terms {
            sum.add(term);
        }
        let [high, low] = sum.parts();
        Self { high, low }
    }

    pub(super) fn value(self) -> f64 {
        self.high + self.low
    }

    pub(super) fn finite(self) -> bool {
        self.high.is_finite() && self.low.is_finite()
    }

    pub(super) fn negative(self) -> Self {
        Self {
            high: -self.high,
            low: -self.low,
        }
    }

    pub(super) fn sign(self) -> f64 {
        if self.high != 0.0 {
            self.high.signum()
        } else if self.low != 0.0 {
            self.low.signum()
        } else {
            0.0
        }
    }

    pub(super) fn magnitude(self) -> f64 {
        self.high.abs() + self.low.abs()
    }

    pub(super) fn plus(self, other: Self) -> Self {
        Self::accumulated(&[self.high, other.high, self.low, other.low])
    }

    pub(super) fn minus(self, other: Self) -> Self {
        self.plus(other.negative())
    }

    pub(super) fn times(self, other: Self) -> Self {
        let high = self.high * other.high;
        Self::accumulated(&[
            high,
            self.high.mul_add(other.high, -high),
            self.high * other.low,
            self.low * other.high,
            self.low * other.low,
        ])
    }

    pub(super) fn scaled(self, factor: f64) -> Self {
        self.times(Self::new(factor))
    }

    pub(super) fn divided(self, other: Self) -> Self {
        let leading = self.high / other.high;
        let residual = self.minus(other.scaled(leading));
        let correction = residual.high / other.high;
        let remaining = residual.minus(other.scaled(correction));
        Self::accumulated(&[leading, correction, remaining.high / other.high])
    }

    pub(super) fn square_root(self) -> Self {
        if self.high == 0.0 && self.low == 0.0 {
            return Self::new(0.0);
        }
        let mut root = Self::new(self.high.sqrt());
        // Two Newton residual corrections retain the low root word, including
        // the quadratic error of the leading binary64 square-root rounding.
        for _ in 0..2 {
            let residual = self.minus(root.times(root));
            root = root.plus(residual.divided(root.scaled(2.0)));
        }
        root
    }
}

#[cfg(test)]
#[path = "kepler_words_tests.rs"]
mod tests;
