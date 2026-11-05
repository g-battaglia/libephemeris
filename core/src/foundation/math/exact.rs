//! Exact rational arithmetic and outward-rounded intervals.
//!
//! Ordered layers: [`Rational`] (reduced `i128` exact algebra with a
//! correctly rounded binary64 conversion via integer long division) and
//! [`Interval`] (outward rounding of `+`, `-`, `*`, `/`, `sqrt` over
//! neighboring floats).
//!
//! Frames, units, epochs: rationals and intervals are dimensionless;
//! there are no frames or epochs in this module.
//!
//! These are the primitives that accepted proof components build on. They
//! are deliberately narrow:
//!
//! - [`Rational`] is exact `i128` numerator/denominator arithmetic with
//!   overflow-checked operations. It is exact algebra inside its range, not
//!   a substitute for arbitrary precision: overflow is an error.
//! - [`Interval`] encloses real results of `+`, `-`, `*`, `/`, and `sqrt`
//!   only. Outward rounding rests on one stated fact: IEEE binary64 rounds
//!   those five operations correctly, so the true result is within half an
//!   ulp of the computed one, and one step to the neighboring float on each
//!   side always encloses it (proof sketch in the numerical contract).
//!   Transcendental interval extensions do **not** exist yet; one-ULP
//!   padding of `sin`/`cos`/`exp` would not be a proof and is not offered.
//!
//! No B-16 proof component is implemented here — only the primitives, with
//! arithmetic identities and containment checks.
//!
//! Refusals: typed errors, never silent substitution — zero rational
//! denominators, integer overflow, non-finite or inverted interval bounds,
//! zero-spanning interval divisors, and negative interval square roots.
//!
//! Domains and limits: rationals are exact inside the checked `i128` range
//! and report overflow outside it; the correctly rounded conversion keeps
//! a (proved-unreachable under the reduced `den >= 1` invariant) subnormal
//! arm so a future invariant change fails safe. Intervals enclose the five
//! stated operations only and offer no transcendental extensions.
//!

use crate::errors::{Error, Result};

/// Exact rational number, always stored reduced with a positive denominator.
///
/// Fields are private so the reduced positive-denominator invariant cannot
/// be bypassed; construction goes through [`Rational::new`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rational {
    num: i128,
    den: i128,
}

fn gcd(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let r = a % b;
        a = b;
        b = r;
    }
    a
}

impl Rational {
    /// Construct `num / den`, reduced; zero denominator and `i128::MIN`
    /// negation overflow are errors.
    pub fn new(num: i128, den: i128) -> Result<Self> {
        if den == 0 {
            return Err(Error::RationalDivisionByZero);
        }
        if num == 0 {
            return Ok(Self { num: 0, den: 1 });
        }
        let sign = if den < 0 { -1 } else { 1 };
        let num = num.checked_mul(sign).ok_or(Error::IntegerOverflow)?;
        let den = den.checked_mul(sign).ok_or(Error::IntegerOverflow)?;
        let g = gcd(num.unsigned_abs(), den as u128);
        Ok(Self {
            num: num / g as i128,
            den: den / g as i128,
        })
    }

    /// Integer literal.
    pub fn from_int(n: i128) -> Self {
        Self { num: n, den: 1 }
    }

    /// Numerator of the reduced fraction.
    pub fn num(self) -> i128 {
        self.num
    }

    /// Positive denominator of the reduced fraction.
    pub fn den(self) -> i128 {
        self.den
    }

    /// Correctly rounded binary64 value (round-to-nearest-even).
    ///
    /// Computed by integer long division: 55 significant bits plus a
    /// sticky bit, rounded once with ties-to-even, including a (provably
    /// unreachable — see below) subnormal arm. No intermediate `as f64`
    /// cast enters, so the old cast-numerator/cast-denominator-then-divide
    /// hazard (small inexactness in each cast could compound
    /// past the true nearest float) cannot occur. All `u128`
    /// intermediates are overflow-safe by construction (remainders stay
    /// below the denominator).
    ///
    /// Reachability: the reduced invariant gives `den >= 1`, so
    /// `|value| <= 2^127` (the `e > 1023` arm is defensive) and the
    /// smallest positive value `1 / (2^127 - 1) > 2^-127` stays far above
    /// the subnormal range (kept as a proved-unreachable arm so a future
    /// invariant change fails safe toward analysis, not silent error).
    pub fn to_f64(self) -> f64 {
        if self.num == 0 {
            return 0.0;
        }
        let neg = self.num < 0;
        let n = self.num.unsigned_abs();
        let d = self.den as u128;
        // (sig, sticky, exp): top 55 significant bits with the leading 1
        // at bit 54, sticky folding every lower bit, exp = binary exponent
        // of the leading bit. 55 (not 54) bits so the subnormal guard bit
        // below always lands on a retained bit (see below).
        let (sig, sticky, exp) = sig_bits(n, d);
        let sign = if neg { 1u64 << 63 } else { 0 };
        if exp >= -1022 {
            let (m, overflow) = round_sig(sig, 55, sticky, 53);
            let e = exp + if overflow { 1 } else { 0 };
            if e > 1023 {
                // Unreachable under the den >= 1 invariant (|value| <=
                // 2^127), kept as a defensive arm.
                return f64::from_bits(sign | 0x7FF0_0000_0000_0000);
            }
            return f64::from_bits(sign | ((e + 1023) as u64) << 52 | (m & 0xF_FFFF_FFFF_FFFF));
        }
        // Subnormal range: round W = value * 2^1074 = (sig + f) / 2^shift
        // to the nearest integer (ties to even), with 0 <= f < 1 folded in
        // `sticky`. From value = (sig+f) * 2^(exp-54): shift = -(exp+1020),
        // always >= 3 here (exp <= -1023), so the unknown fraction f sits
        // strictly below the guard bit and only ever feeds `sticky`; the
        // guard bit always lands on a retained sig bit. For shift > 55 the
        // whole value is strictly below half a unit, hence exactly 0.
        let shift = (-(exp + 1020)) as u32;
        if shift > 55 {
            return f64::from_bits(sign);
        }
        let i = sig >> shift;
        let guard = shift >= 1 && ((sig >> (shift - 1)) & 1) == 1;
        let lower_mask = if shift <= 1 {
            0
        } else {
            (1u64 << (shift - 1)) - 1
        };
        let sticky = sticky || (sig & lower_mask) != 0;
        let mut m = i;
        if guard && (sticky || m & 1 == 1) {
            m += 1;
        }
        if m >= 1 << 52 {
            // Rounded up into the smallest normal: 2^-1022.
            return f64::from_bits(sign | 1u64 << 52);
        }
        f64::from_bits(sign | m)
    }

    /// Exact sum; overflow is an error.
    pub fn checked_add(self, other: Self) -> Result<Self> {
        let n = self
            .num
            .checked_mul(other.den)
            .and_then(|a| {
                other
                    .num
                    .checked_mul(self.den)
                    .and_then(|b| a.checked_add(b))
            })
            .ok_or(Error::IntegerOverflow)?;
        let d = self
            .den
            .checked_mul(other.den)
            .ok_or(Error::IntegerOverflow)?;
        Self::new(n, d)
    }

    /// Exact difference.
    pub fn checked_sub(self, other: Self) -> Result<Self> {
        self.checked_add(Self {
            num: other.num.checked_neg().ok_or(Error::IntegerOverflow)?,
            den: other.den,
        })
    }

    /// Exact product.
    pub fn checked_mul(self, other: Self) -> Result<Self> {
        let n = self
            .num
            .checked_mul(other.num)
            .ok_or(Error::IntegerOverflow)?;
        let d = self
            .den
            .checked_mul(other.den)
            .ok_or(Error::IntegerOverflow)?;
        Self::new(n, d)
    }

    /// Exact quotient; division by zero is an error.
    pub fn checked_div(self, other: Self) -> Result<Self> {
        if other.num == 0 {
            return Err(Error::RationalDivisionByZero);
        }
        self.checked_mul(Self {
            num: other.den,
            den: other.num,
        })
    }
}

/// Top 55 significant bits of `n/d` (both > 0) with sticky and exponent.
///
/// Returns `(sig, sticky, exp)`: `sig` in `[2^54, 2^55)` with the leading 1
/// at bit 54, `sticky` true iff any discarded lower bit (or the exact
/// remainder) is nonzero, `exp` the binary exponent of the leading bit, so
/// `n/d = (sig + f) * 2^(exp - 54)` for some `f` in `[0, 1)` with `f != 0`
/// exactly when `sticky`. Remainder discipline keeps every `u128`
/// intermediate overflow-safe: remainders stay below `d <= 2^127 - 1`, so
/// doubling them stays below `2^128 - 1`.
fn sig_bits(n: u128, d: u128) -> (u64, bool, i32) {
    let q = n / d;
    let mut r = n % d;
    if q > 0 {
        let bits = (128 - q.leading_zeros()) as i32;
        let exp = bits - 1;
        if bits >= 55 {
            let shift = (bits - 55) as u32;
            let v = (q >> shift) as u64;
            let sticky = r != 0 || (q & ((1u128 << shift) - 1)) != 0;
            (v, sticky, exp)
        } else {
            // `q` already holds `bits` significant bits; append fraction
            // bits to reach exactly 55 (no pre-shift: shifting first and
            // then appending overfills the field — FIX-02 debugging catch).
            let mut v = q as u64;
            for _ in 0..(55 - bits) {
                r *= 2;
                v <<= 1;
                if r >= d {
                    v |= 1;
                    r -= d;
                }
            }
            (v, r != 0, exp)
        }
    } else {
        // Pure fraction: skip leading zeros, then collect 55 significant
        // bits. n/d >= 1/(2^127 - 1) > 2^-127, so the leading 1 arrives
        // within 127 steps; the cap below is unreachable-safe.
        let mut exp = 0i32;
        let mut v = 0u64;
        let mut count = 0;
        for j in 1..=400u32 {
            r *= 2;
            let bit = r >= d;
            if bit {
                r -= d;
            }
            if count == 0 {
                if !bit {
                    continue;
                }
                exp = -(j as i32);
            }
            v = (v << 1) | (bit as u64);
            count += 1;
            if count == 55 {
                break;
            }
        }
        debug_assert_eq!(count, 55);
        (v, r != 0, exp)
    }
}

/// Round a `bits`-bit significand with sticky to `target` bits,
/// ties-to-even. Returns `(rounded, overflowed)`, where overflow means the
/// result carried past `2^target` (caller bumps the exponent).
fn round_sig(mut v: u64, mut bits: u32, mut sticky: bool, target: u32) -> (u64, bool) {
    debug_assert!(bits > target && bits <= 64 && target >= 1);
    while bits > target + 1 {
        sticky |= v & 1 == 1;
        v >>= 1;
        bits -= 1;
    }
    // Exactly target+1 bits remain: the low bit is the guard.
    let guard = v & 1 == 1;
    let mut m = v >> 1;
    if guard && (sticky || m & 1 == 1) {
        m += 1;
    }
    if m == 1 << target {
        (m >> 1, true)
    } else {
        (m, false)
    }
}

/// Closed interval `[lo, hi]` with finite bounds, `lo <= hi`.
///
/// Fields are private so the ordering invariant cannot be bypassed;
/// construction goes through [`Interval::new`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Interval {
    lo: f64,
    hi: f64,
}

impl Interval {
    /// Construct from bounds; rejects non-finite and inverted bounds.
    pub fn new(lo: f64, hi: f64) -> Result<Self> {
        if !lo.is_finite() || !hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: lo + hi });
        }
        if lo > hi {
            return Err(Error::IntervalInvertedBounds);
        }
        Ok(Self { lo, hi })
    }

    /// Point interval (a single float; encloses exactly that float).
    pub fn point(x: f64) -> Result<Self> {
        Self::new(x, x)
    }

    /// Lower bound.
    pub fn lo(self) -> f64 {
        self.lo
    }

    /// Upper bound.
    pub fn hi(self) -> f64 {
        self.hi
    }

    /// True when `x` lies within the bounds.
    pub fn contains(self, x: f64) -> bool {
        self.lo <= x && x <= self.hi
    }

    /// Width `hi - lo` (non-negative by construction).
    pub fn width(self) -> f64 {
        self.hi - self.lo
    }

    fn enclose(lo: f64, hi: f64) -> Result<Self> {
        if !lo.is_finite() || !hi.is_finite() {
            return Err(Error::NonFiniteArgument { value: lo + hi });
        }
        Self::new(lo.next_down(), hi.next_up())
    }

    /// Enclosing sum.
    pub fn checked_add(self, other: Self) -> Result<Self> {
        Self::enclose(self.lo + other.lo, self.hi + other.hi)
    }

    /// Enclosing difference.
    pub fn checked_sub(self, other: Self) -> Result<Self> {
        Self::enclose(self.lo - other.hi, self.hi - other.lo)
    }

    /// Enclosing product (outward rounding of all four corner products).
    pub fn checked_mul(self, other: Self) -> Result<Self> {
        let corners = [
            self.lo * other.lo,
            self.lo * other.hi,
            self.hi * other.lo,
            self.hi * other.hi,
        ];
        let mut lo = corners[0];
        let mut hi = corners[0];
        for &c in &corners[1..] {
            if !c.is_finite() {
                return Err(Error::NonFiniteArgument { value: c });
            }
            lo = lo.min(c);
            hi = hi.max(c);
        }
        Self::enclose(lo, hi)
    }

    /// Enclosing quotient; a divisor spanning zero is an error (the true
    /// quotient is unbounded there, not representable here).
    pub fn checked_div(self, other: Self) -> Result<Self> {
        if other.lo <= 0.0 && other.hi >= 0.0 {
            return Err(Error::IntervalSpansZero);
        }
        let corners = [
            self.lo / other.lo,
            self.lo / other.hi,
            self.hi / other.lo,
            self.hi / other.hi,
        ];
        let mut lo = corners[0];
        let mut hi = corners[0];
        for &c in &corners[1..] {
            if !c.is_finite() {
                return Err(Error::NonFiniteArgument { value: c });
            }
            lo = lo.min(c);
            hi = hi.max(c);
        }
        Self::enclose(lo, hi)
    }

    /// Enclosing square root; negative lower bound is a domain error.
    pub fn sqrt(self) -> Result<Self> {
        if self.lo < 0.0 {
            return Err(Error::IntervalDomainViolation);
        }
        Self::enclose(self.lo.sqrt(), self.hi.sqrt())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(num: i128, den: i128) -> Rational {
        Rational::new(num, den).unwrap()
    }

    #[test]
    fn rational_exact_identities() {
        // 1/3 + 1/6 = 1/2 exactly; (2/3)(3/4) = 1/2; 1/2 - 1/3 = 1/6.
        assert_eq!(r(1, 3).checked_add(r(1, 6)).unwrap(), r(1, 2));
        assert_eq!(r(2, 3).checked_mul(r(3, 4)).unwrap(), r(1, 2));
        assert_eq!(r(1, 2).checked_sub(r(1, 3)).unwrap(), r(1, 6));
        assert_eq!(r(3, 7).checked_div(r(6, 7)).unwrap(), r(1, 2));
    }

    #[test]
    fn rational_normalization_and_zero() {
        assert_eq!(Rational::new(-2, -4).unwrap(), r(1, 2));
        assert_eq!(Rational::new(3, -9).unwrap(), r(-1, 3));
        assert_eq!(Rational::new(0, -5).unwrap(), Rational::from_int(0));
        assert_eq!(r(6, 8), r(3, 4));
    }

    #[test]
    fn rational_to_f64_is_correctly_rounded() {
        // Review probe (FIX-02): independent Python Fraction arithmetic
        // puts the nearest float to 9007199254740993/7 at these bits; the
        // old cast-and-divide produced...492.
        assert_eq!(
            r(9007199254740993, 7).to_f64().to_bits(),
            0x4312_4924_9249_2493
        );
        // Exact values and zeros.
        assert_eq!(r(42, 1).to_f64(), 42.0);
        assert_eq!(r(-7, 1).to_f64(), -7.0);
        assert_eq!(r(0, 5).to_f64().to_bits(), 0.0_f64.to_bits());
        // Exact midpoints round ties-to-even: (2^53+1) and (2^53+3) sit
        // exactly halfway between adjacent floats.
        assert_eq!(
            r(9007199254740993, 1).to_f64().to_bits(),
            0x4340_0000_0000_0000
        );
        assert_eq!(
            r(9007199254740995, 1).to_f64().to_bits(),
            0x4340_0000_0000_0002
        );
        assert_eq!(
            r(9007199254740993, 2).to_f64().to_bits(),
            0x4330_0000_0000_0000
        );
        // Ordinary fractions agree with correctly-rounded hardware division.
        assert_eq!(r(1, 3).to_f64(), 1.0 / 3.0);
        assert_eq!(r(-1, 3).to_f64(), -1.0 / 3.0);
        assert_eq!(r(22, 7).to_f64(), 22.0 / 7.0);
        // Extremes: exact powers and hardware cross-check (-2^127 is
        // 0xC7E0..., verified against Python struct packing, not hand-set).
        assert_eq!(
            Rational::from_int(i128::MIN).to_f64().to_bits(),
            0xC7E0_0000_0000_0000
        );
        assert_eq!(Rational::from_int(i128::MAX).to_f64(), i128::MAX as f64);
        // Accessors expose the reduced invariant.
        assert_eq!((r(6, 8).num(), r(6, 8).den()), (3, 4));
    }

    #[test]
    fn rational_to_f64_matches_independent_oracle() {
        // Oracle vectors from Python Fraction arithmetic (seed 20260922;
        // correctly rounded by construction). The first four are
        // double-rounding witnesses: naive `(n as f64) / (d as f64)` gives
        // a different float there, so this table genuinely discriminates
        // correct rounding from cast-and-divide.
        let cases: &[(i128, i128, u64)] = &[
            (
                77286021465917791728876059106508869789,
                1310364285758773085537290346944367863,
                0x404d7d830ad391ee,
            ),
            (
                -38513332450145456663089138216238380882,
                161901408864662527297456008813075182735,
                0xbfce72e5beb90f5d,
            ),
            (
                -162997372889253240634252603048763145073,
                10312145515319871565058952417929139181,
                0xc02f9cd9ccfc6ec9,
            ),
            (
                115440357504663945277177456024086378149,
                108576228604466501590060350389352040146,
                0x3ff102f264841c74,
            ),
            (
                57900870944242203967641322758558883473,
                10265302066415789607988383010485483740,
                0x40168fd0cb5f4876,
            ),
            (
                45210845522179827542146108606280681739,
                19542994004973576047343417636001941996,
                0x400281da174bde7d,
            ),
            (
                -923327437311243264063250615715546977,
                53250801987345030039362303117142810792,
                0xbf91c15f70cc9d6d,
            ),
            (
                -164142565447789360666396725049488128702,
                131456056576294938716719591161823775670,
                0xbff3fa781b473e85,
            ),
            (
                16144684579046641305913836462477365670,
                8747880437152394236856690415181314123,
                0x3ffd87637e2c84aa,
            ),
            (
                8901049831292063560287548019216879797,
                20341986315675905814732287913505319257,
                0x3fdc0127050a31c2,
            ),
            (
                59569541711201930690496616286632903026,
                148947242418045520823694261410430781633,
                0x3fd9989228be3cf9,
            ),
            (
                -138351089386602724441782654274179746659,
                35368947911524840141531266364468720723,
                0xc00f4b114b46942d,
            ),
        ];
        for &(n, d, bits) in cases {
            assert_eq!(r(n, d).to_f64().to_bits(), bits, "n={n} d={d}");
        }
    }

    #[test]
    fn rational_errors() {
        assert!(matches!(
            Rational::new(1, 0),
            Err(Error::RationalDivisionByZero)
        ));
        assert!(matches!(
            r(1, 2).checked_div(Rational::from_int(0)),
            Err(Error::RationalDivisionByZero)
        ));
        assert!(matches!(
            Rational::from_int(i128::MAX).checked_add(Rational::from_int(1)),
            Err(Error::IntegerOverflow)
        ));
        assert!(matches!(
            Rational::from_int(i128::MAX).checked_mul(Rational::from_int(2)),
            Err(Error::IntegerOverflow)
        ));
    }

    #[test]
    fn intervals_contain_exact_results() {
        // Singleton arithmetic encloses the exact integer result.
        let three = Interval::point(1.0)
            .unwrap()
            .checked_add(Interval::point(2.0).unwrap())
            .unwrap();
        assert!(three.contains(3.0));
        assert!(three.width() > 0.0 && three.width() < 1e-14);
        // 0.1 + 0.2 encloses the real 0.3 even though neither endpoint is
        // exact in binary64.
        let s = Interval::point(0.1)
            .unwrap()
            .checked_add(Interval::point(0.2).unwrap())
            .unwrap();
        assert!(s.lo() <= 0.3 && 0.3 <= s.hi());
        // x - x always contains 0.
        let x = Interval::new(1.23456789, 9.87654321).unwrap();
        let z = x.checked_sub(x).unwrap();
        assert!(z.contains(0.0));
    }

    #[test]
    fn interval_product_quotient_sqrt_containment() {
        let a = Interval::new(1.0, 2.0).unwrap();
        let b = Interval::new(3.0, 4.0).unwrap();
        let p = a.checked_mul(b).unwrap();
        assert!(p.lo() <= 3.0 && 8.0 <= p.hi());
        let q = b.checked_div(a).unwrap();
        assert!(q.lo() <= 1.5 && 4.0 <= q.hi());
        // sqrt([1, 4]) encloses [1, 2].
        let s = Interval::new(1.0, 4.0).unwrap().sqrt().unwrap();
        assert!(s.lo() <= 1.0 && 2.0 <= s.hi());
        // Negative products: enclosure covers the full corner range.
        let n = Interval::new(-2.0, 3.0)
            .unwrap()
            .checked_mul(Interval::new(-4.0, -1.0).unwrap())
            .unwrap();
        assert!(n.lo() <= -12.0 && 8.0 <= n.hi());
    }

    #[test]
    fn interval_errors() {
        assert!(matches!(
            Interval::new(2.0, 1.0),
            Err(Error::IntervalInvertedBounds)
        ));
        assert!(matches!(
            Interval::new(1.0, 2.0)
                .unwrap()
                .checked_div(Interval::new(-1.0, 1.0).unwrap()),
            Err(Error::IntervalSpansZero)
        ));
        assert!(matches!(
            Interval::new(-4.0, 1.0).unwrap().sqrt(),
            Err(Error::IntervalDomainViolation)
        ));
    }
}
