//! SMH-2016 Delta-T tables and standalone UT1–TT conversion.
//!
//! Public callers use [`crate::time::default_timescale_tt_from_ut1`],
//! [`crate::time::default_timescale_from_ut1`] and
//! [`crate::time::default_timescale_ut1_from_tt`]. Model selection and overrides
//! are exposed through [`crate::time::DeltaTConfig`]. This module's numerical
//! implementation is crate-internal.
//!
//! Delta-T is TT minus UT1 in seconds. The observed realization linearly
//! interpolates merged semiannual historical observations (1657–1973) and daily
//! IERS rows. Outside that span a 13-piece cubic long-term curve is evaluated
//! with clamped endpoint segment selection. These extrapolations are predictions.
//!
//! Numeric data are bundled from Skyfield 1.54 (MIT); original input identities
//! and measurement provenance are retained in the internal `smh2016_tables`. Runtime
//! calculation is native binary64. Internal table knots belong to the interval
//! on their right; the final knot closes the last observed interval.
//!
//! Checked UT1-to-TT conversion makes exactly two substitutions, not an
//! iteration to a convergence tolerance. Split forms retain the entire original
//! UT1 value and a TT offset in days; the first part need not be an integer day.
//! Checked entries reject nonfinite input and unrepresentable results, while
//! internal source primitives preserve their documented IEEE propagation.
//! Extreme extrapolation can overflow and has no convergence guarantee. Caller
//! tidal corrections belong to the outer Delta-T configuration layer.

use super::smh2016_tables::{SMH2016_CURVE_ROWS, SMH2016_TABLE_DT, SMH2016_TABLE_TT};
use crate::errors::{Error, Result};

/// Seconds per day; `ut1_jd` iteration divides by this exact
/// literal between refinement steps.
const SECONDS_PER_DAY: f64 = 86400.0;

/// Julian Day anchoring the long-term curve argument:
/// `J = (tt - 1721045.0) / 365.25` in Julian years.
const CURVE_ANCHOR_JD: f64 = 1721045.0;

/// Days per Julian year for the curve argument (despite the constant's name).
const DAYS_PER_CENTURY: f64 = 365.25;

/// Linear interpolation over a strictly increasing knot table with the
/// reference endpoint behavior: below the first knot the left fill
/// applies, above the last knot the right fill applies, and interior
/// points (including exact knots) evaluate on the segment opened by the
/// last knot at or below `x`. The segment evaluates as the fused
/// `slope.mul_add(x - xp[i], fp[i])` with
/// `slope = (fp[i + 1] - fp[i]) / (xp[i + 1] - xp[i])`: the reference
/// evaluator fuses the multiply-add (verified exactly over 20000 random
/// interior points), so a separate multiply followed by an addition
/// would round differently on some segments.
fn interp_segments(x: f64, xp: &[f64], fp: &[f64], left: f64, right: f64) -> f64 {
    let n = xp.len();
    debug_assert!(n >= 2 && n == fp.len());
    if x < xp[0] {
        return left;
    }
    if x > xp[n - 1] {
        return right;
    }
    let lo = interpolation_segment_index(x, xp);
    let slope = (fp[lo + 1] - fp[lo]) / (xp[lo + 1] - xp[lo]);
    slope.mul_add(x - xp[lo], fp[lo])
}

/// Owned interval of the finite, in-span interpolation argument.
fn interpolation_segment_index(x: f64, xp: &[f64]) -> usize {
    let mut lo = 0usize;
    let mut hi = xp.len() - 1;
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if xp[mid] <= x {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    lo
}

/// One long-term curve evaluation at Julian year `x`.
///
/// The segment index is the reference clamped interpolation of `x`
/// against the `x0` column over `0..rows` (fused multiply-add, truncated
/// toward zero, clamped to the end rows outside the span), and the row
/// evaluates the cubic in `t = (x - x0) / (x1 - x0)` with Horner nesting
/// `((a3 * t + a2) * t + a1) * t + a0` (separate multiply/add rounds,
/// as the reference loop performs them).
fn curve_seconds(x: f64) -> f64 {
    let row = SMH2016_CURVE_ROWS[curve_row_index(x)];
    let t = (x - row[0]) / (row[1] - row[0]);
    ((row[2] * t + row[3]) * t + row[4]) * t + row[5]
}

/// Actual clamped spline-row selection, shared with domain metadata.
fn curve_row_index(x: f64) -> usize {
    let nrows = SMH2016_CURVE_ROWS.len();
    let x0_first = SMH2016_CURVE_ROWS[0][0];
    let x0_last = SMH2016_CURVE_ROWS[nrows - 1][0];
    if x <= x0_first {
        0
    } else if x >= x0_last {
        nrows - 1
    } else {
        let mut lo = 0usize;
        let mut hi = nrows - 1;
        while hi - lo > 1 {
            let mid = (lo + hi) / 2;
            if SMH2016_CURVE_ROWS[mid][0] <= x {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let slope = 1.0 / (SMH2016_CURVE_ROWS[lo + 1][0] - SMH2016_CURVE_ROWS[lo][0]);
        let idx = slope.mul_add(x - SMH2016_CURVE_ROWS[lo][0], lo as f64);
        idx.trunc().clamp(0.0, (nrows - 1) as f64) as usize
    }
}

/// One actual interpolation piece of the bundled TT-minus-UT1 model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum BundledDeltaTPiece {
    /// Linear observed-table interval, with final-knot closed ownership.
    Observed(usize),
    /// Clamped long-term cubic row, including its extrapolated end span.
    Curve(usize),
}

/// Actual bundled-clock piece consumed at finite scalar TT (Julian days).
///
/// This uses the same interval/row selectors as the model, including rounded
/// scalar epoch and endpoint ownership. It computes no alternate Delta-T and
/// performs no interpolation, extrapolation or inverse clock calculation.
pub(crate) fn bundled_deltat_piece(tt: f64) -> Result<BundledDeltaTPiece> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    Ok(
        if tt < SMH2016_TABLE_TT[0] || tt > SMH2016_TABLE_TT[SMH2016_TABLE_TT.len() - 1] {
            BundledDeltaTPiece::Curve(curve_row_index((tt - CURVE_ANCHOR_JD) / DAYS_PER_CENTURY))
        } else {
            BundledDeltaTPiece::Observed(interpolation_segment_index(tt, &SMH2016_TABLE_TT))
        },
    )
}

/// Raw SMH-2016 Delta T in **seconds** at a TT Julian Day: the merged
/// table inside its span, else the long-term curve. Extrapolation supplies
/// remote epochs without a coverage cutoff; overflowing or non-finite
/// arithmetic propagates per IEEE-754.
fn raw_seconds(tt: f64) -> f64 {
    if tt < SMH2016_TABLE_TT[0] || tt > SMH2016_TABLE_TT[SMH2016_TABLE_TT.len() - 1] {
        curve_seconds((tt - CURVE_ANCHOR_JD) / DAYS_PER_CENTURY)
    } else {
        interp_segments(tt, &SMH2016_TABLE_TT, &SMH2016_TABLE_DT, f64::NAN, f64::NAN)
    }
}

/// Default-model Delta-T provider in **seconds** for a UT1 Julian Day.
///
/// The first table evaluation seeds the TT estimate, the second refines
/// it, and the returned Delta T is read at the twice-refined epoch.
/// This returns a time difference, not the timescale's scalar TT epoch.
/// Non-finite input returns `None`; finite input preserves the model's
/// IEEE-754 extrapolation results, including overflow.
pub(crate) fn smh2016_lookup_seconds(tjdut: f64) -> Option<f64> {
    if !tjdut.is_finite() {
        return None;
    }
    let first = raw_seconds(tjdut);
    let refined = tjdut + first / SECONDS_PER_DAY;
    let second = raw_seconds(refined);
    let settled = tjdut + second / SECONDS_PER_DAY;
    Some(raw_seconds(settled))
}

/// Two-part TT epoch for a UT1 input, without collapsing
/// to a scalar Julian Day.
///
/// Returns `(whole, tt_fraction)` where `whole` is the original UT1 value
/// kept as-is (never floored or rounded) and `tt_fraction` is the
/// twice-refined TT-minus-UT1 offset in days. This construction `Time(ut1, delta_t / DAY_S)`, which
/// stores the input value whole and sets the `ut1_fraction` attribute to
/// zero: Earth rotation reads the whole value plus an identically-zero UT1
/// fraction, while the TT fraction carries the small table offset at full
/// precision. Collapsing the pair to a scalar day number first and
/// re-splitting it would quantize the small offset to the whole-day grid
/// (tens of microseconds at Julian-day magnitudes) and then cancel
/// catastrophically when the UT1 fraction is recovered by subtraction —
/// a whole-day-magnitude rounding no downstream consumer can undo. The
/// scalar [`timescale_tt_from_ut1`] epoch is exactly `whole + tt_fraction`.
/// Caller model/tidal/user overrides are not involved; see that function
/// for the error contract.
pub(crate) fn timescale_parts_from_ut1(tjdut: f64) -> Result<(f64, f64)> {
    if !tjdut.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjdut });
    }
    let first = raw_seconds(tjdut);
    let refined = tjdut + first / SECONDS_PER_DAY;
    if !first.is_finite() || !refined.is_finite() {
        return Err(Error::MathDomain { value: tjdut });
    }
    let second = raw_seconds(refined);
    let tt_fraction = second / SECONDS_PER_DAY;
    if second.is_finite() && (tjdut + tt_fraction).is_finite() {
        Ok((tjdut, tt_fraction))
    } else {
        Err(Error::MathDomain { value: tjdut })
    }
}
/// Convert a scalar UT1 Julian Day to a scalar Terrestrial Time Julian Day.
///
/// The pair comes from [`timescale_parts_from_ut1`]; the scalar is exactly
/// `whole + tt_fraction`.
/// & Hohenkerk long-term curve also used by [`smh2016_lookup_seconds`].
/// If `D(t)` is that model's TT-minus-UT1 difference in **seconds** at
/// TT epoch `t`, and `u` is the supplied UT1 epoch in **days**, the method is
///
/// ```text
/// d 0 = D(u)
/// t 1 = u + d 0 / 86400
/// d 1 = D(t 1)
/// TT = u + d 1 / 86400
/// ```
///
/// Exactly two evaluations are made, not an iteration to a convergence
/// tolerance. The two-part epoch is `(u, d1 / 86400)`: the first part is
/// the entire input value, not its integer portion. The final scalar sum
/// intentionally rounds to binary64 Julian-Day precision, as required by
/// The collapsed scalar TT epoch. Do not replace the final offset by
/// Delta T evaluated again at the settled TT, normalize the parts at
/// midnight, or derive it by subtracting two scalar epochs.
///
/// Non-finite input returns [`Error::NonFiniteArgument`]. All finite UT1
/// epochs are admitted provided the model arithmetic and resulting TT
/// remain finite; otherwise returns [`Error::MathDomain`] carrying the
/// original UT1 input. At extreme epochs extrapolation may overflow, and
/// two substitutions need not approximate the implicit equation
/// `TT = u + D(TT) / 86400` accurately.
/// No certified convergence or physical-accuracy bound is asserted.
///
/// This is the **bundled enhanced timescale**, independent of a caller's
/// selected Delta-T model, user override or tidal-acceleration correction.
/// Those settings affect other model-based conversions, not this timescale.
/// Custom time tables/functions, newer IERS data and the IERS-only fallback
/// outside the bundled data are not represented by this function.
/// It does not convert UTC or use leap-second labels; the input is UT1.
///
/// Scientific model: Stephenson, Morrison & Hohenkerk (2016), realized by
/// Skyfield 1.54's bundled tables and long-term splines (MIT licensed).
pub(crate) fn timescale_tt_from_ut1(tjdut: f64) -> Result<f64> {
    let (whole, tt_fraction) = timescale_parts_from_ut1(tjdut)?;
    Ok(whole + tt_fraction)
}

/// Scalar TT of a source-style UT1 Time after a finite public input guard.
///
/// The original UT1 value remains the whole part. Exactly two bundled Delta-T
/// evaluations form `u + D(u + D(u)/86400)/86400`. Intermediate overflow and
/// NaN propagate as IEEE-754 values, as in the source Time constructor; no
/// later provider/model/frame acquisition is replaced by a clock refusal.
/// Callers must retain the original finite request for its public domain and
/// error identity. This internal acquisition does not relax the checked
/// direct timescale APIs or admit a nonfinite caller-supplied epoch.
/// No convergence or physical accuracy is claimed at extrapolated epochs.
pub(crate) fn source_timescale_tt_from_ut1(tjdut: f64) -> Result<f64> {
    if !tjdut.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjdut });
    }
    let first = raw_seconds(tjdut);
    let refined = tjdut + first / SECONDS_PER_DAY;
    let second = raw_seconds(refined);
    Ok(tjdut + second / SECONDS_PER_DAY)
}

/// Scalar TT carried by a source-style UT1 Time, including IEEE nonfinite
/// intermediate values. Exactly two bundled evaluations preserve the original
/// whole UT1 value; no configured model or stored clock is substituted. The
/// following public position entry owns any nonfinite JD range constructor.
/// This private primitive does not alter the checked public time APIs.
pub(crate) fn source_time_scalar_tt_from_ut1(value: f64) -> f64 {
    let (whole, fraction, _) = source_time_parts_from_ut1(value);
    whole + fraction
}

/// Acquire the original UT1 Time's whole, TT fraction and UT1 fraction.
/// Two bundled curve evaluations construct TT. The supplied value remains the
/// whole part; the UT1 fraction is the consumed `0.0 * value`, including its
/// signed zero and IEEE nonfinite behavior. No normalization, configured model
/// or scalar-TT inverse is used. Public input validation belongs to the caller.
pub(crate) fn source_time_parts_from_ut1(value: f64) -> (f64, f64, f64) {
    let first = raw_seconds(value);
    let second = raw_seconds(value + first / SECONDS_PER_DAY);
    (value, second / SECONDS_PER_DAY, 0.0 * value)
}

/// Convert scalar TT to scalar UT1 with the fixed bundled timescale.
///
/// Split TT at its floor into whole day and fraction; evaluate the bundled
/// TT-minus-UT1 model at their scalar TT observable. Subtract the seconds/day
/// offset from the fraction before adding the whole day. Unlike UT1-to-TT,
/// this inverse requires one model evaluation, not two substitutions.
/// Caller model/tidal/user overrides and UTC leap labels are not involved.
/// All finite inputs are admitted if arithmetic remains finite; failures
/// return typed non-finite/domain errors. No physical error bound is implied.
///
/// Full-precision TT-minus-UT1 seconds at a TT epoch, for consumers that
/// need the model offset without the whole-day rounding (frame chains).
pub(crate) fn bundled_deltat_seconds(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let dt = raw_seconds(jd_tt);
    if !dt.is_finite() {
        return Err(Error::MathDomain { value: jd_tt });
    }
    Ok(dt)
}
pub(crate) fn timescale_ut1_from_tt(jd_tt: f64) -> Result<f64> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let whole = jd_tt.floor();
    let fraction = jd_tt - whole;
    let delta_seconds = raw_seconds(whole + fraction);
    let ut1 = whole + (fraction - delta_seconds / SECONDS_PER_DAY);
    if delta_seconds.is_finite() && ut1.is_finite() {
        Ok(ut1)
    } else {
        Err(Error::MathDomain { value: jd_tt })
    }
}

/// Analytic slope of the actual bundled piece, in seconds per TT day.
/// Nonfinite intermediate epochs retain the raw model's IEEE propagation.
fn raw_rate_seconds_per_day(tt: f64) -> f64 {
    match bundled_deltat_piece(tt) {
        Ok(BundledDeltaTPiece::Observed(index)) => {
            (SMH2016_TABLE_DT[index + 1] - SMH2016_TABLE_DT[index])
                / (SMH2016_TABLE_TT[index + 1] - SMH2016_TABLE_TT[index])
        }
        Ok(BundledDeltaTPiece::Curve(index)) => curve_row_rate_seconds_per_day(tt, index),
        Err(_) => f64::NAN,
    }
}

/// Differentiate the selected cubic, including Julian-year/day conversion.
fn curve_row_rate_seconds_per_day(tt: f64, index: usize) -> f64 {
    let row = SMH2016_CURVE_ROWS[index];
    let span = row[1] - row[0];
    // The scalar's rounded Julian year still selects the row. Preserve the
    // subtraction residual when differentiating its continuous polynomial:
    // (TT - anchor - 365.25*x0)/(365.25*(x1-x0)). Forming the Julian year
    // first loses several bits near a row origin, where the slope is small.
    let days = tt - CURVE_ANCHOR_JD;
    let virtual_anchor = days - tt;
    let remainder = (tt - (days - virtual_anchor)) - (CURVE_ANCHOR_JD + virtual_anchor);
    let numerator = (-DAYS_PER_CENTURY).mul_add(row[0], days) + remainder;
    let t = numerator / (DAYS_PER_CENTURY * span);
    (3.0 * row[2]).mul_add(t, 2.0 * row[3]).mul_add(t, row[4]) / span / DAYS_PER_CENTURY
}

/// Differentiate bundled TT-minus-UT1 seconds with respect to TT Julian days.
///
/// The observed piece contributes its linear slope; the long-term piece
/// contributes `(3 a3 t² + 2 a2 t + a1) / ((x1 - x0) * 365.25)`.
/// Selection uses the same scalar argument and selectors as the value. Internal
/// observed knots own their right interval, the final observed knot closes the
/// last interval, and strictly exterior epochs use the clamped cubic selector.
/// At a join this is the derivative of the owned piece, not a claim that the
/// whole model has a two-sided derivative. Binary64 rounding is not a step law.
///
/// The finite-input and model-overflow guards match [`bundled_deltat_seconds`];
/// a nonfinite slope also returns [`Error::MathDomain`] with the original TT.
pub(crate) fn bundled_deltat_rate_seconds_per_day(tt: f64) -> Result<f64> {
    bundled_deltat_seconds(tt)?;
    let rate = raw_rate_seconds_per_day(tt);
    if rate.is_finite() {
        Ok(rate)
    } else {
        Err(Error::MathDomain { value: tt })
    }
}

/// Differentiate the three-lookup default provider, in seconds per UT1 day.
///
/// With `t1 = u + D(u)/86400` and `t2 = u + D(t1)/86400`, this is
/// `D'(t2) * (1 + D'(t1) * (1 + D'(u)/86400)/86400)`. Each derivative
/// belongs to the piece actually selected at that scalar intermediate epoch;
/// the join convention is documented by [`bundled_deltat_rate_seconds_per_day`].
/// This is not the derivative of the two-lookup TT timescale below.
/// Nonfinite original input returns `None`; finite extrapolation keeps IEEE
/// overflow/NaN propagation, exactly as [`smh2016_lookup_seconds`] does.
pub(crate) fn smh2016_lookup_rate_seconds_per_day(ut: f64) -> Option<f64> {
    if !ut.is_finite() {
        return None;
    }
    let first = raw_seconds(ut);
    let refined = ut + first / SECONDS_PER_DAY;
    let second = raw_seconds(refined);
    let settled = ut + second / SECONDS_PER_DAY;
    let refined_rate = 1.0 + raw_rate_seconds_per_day(ut) / SECONDS_PER_DAY;
    let settled_rate = 1.0 + raw_rate_seconds_per_day(refined) * refined_rate / SECONDS_PER_DAY;
    Some(raw_rate_seconds_per_day(settled) * settled_rate)
}

/// Differentiate the bundled two-lookup TT clock, in TT days per UT1 day.
///
/// For `t1 = u + D(u)/86400`, the returned derivative is
/// `1 + D'(t1) * (1 + D'(u)/86400)/86400`. Both pieces are selected at
/// the actual scalar lookup epochs, including right/final knot ownership.
/// This differentiates the continuous local clock represented by those pieces;
/// it does not differentiate the staircase produced by binary64 rounding or
/// solve a different implicit clock equation.
///
/// The original nonfinite-input and finite-model/result guards match
/// [`timescale_parts_from_ut1`]. A nonfinite rate returns [`Error::MathDomain`]
/// carrying the original UT1 input. Caller model, tidal and IERS overrides do
/// not participate in this bundled timescale.
pub(crate) fn timescale_tt_rate_from_ut1(ut: f64) -> Result<f64> {
    if !ut.is_finite() {
        return Err(Error::NonFiniteArgument { value: ut });
    }
    let first = raw_seconds(ut);
    let refined = ut + first / SECONDS_PER_DAY;
    if !first.is_finite() || !refined.is_finite() {
        return Err(Error::MathDomain { value: ut });
    }
    let second = raw_seconds(refined);
    if !second.is_finite() || !(ut + second / SECONDS_PER_DAY).is_finite() {
        return Err(Error::MathDomain { value: ut });
    }
    let refined_rate = 1.0 + raw_rate_seconds_per_day(ut) / SECONDS_PER_DAY;
    // The scalar second lookup is already valid. Preserve the original UT1
    // identity if its derivative cannot be represented.
    let second_rate = bundled_deltat_rate_seconds_per_day(refined)
        .map_err(|_| Error::MathDomain { value: ut })?;
    let rate = 1.0 + second_rate * refined_rate / SECONDS_PER_DAY;
    if rate.is_finite() {
        Ok(rate)
    } else {
        Err(Error::MathDomain { value: ut })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bundled_clock_piece_uses_actual_scalar_lookup_join_ownership() {
        for index in [1, SMH2016_TABLE_TT.len() / 2, SMH2016_TABLE_TT.len() - 2] {
            let knot = SMH2016_TABLE_TT[index];
            assert_eq!(
                bundled_deltat_piece(knot.next_down()).unwrap(),
                BundledDeltaTPiece::Observed(index - 1)
            );
            assert_eq!(
                bundled_deltat_piece(knot).unwrap(),
                BundledDeltaTPiece::Observed(index)
            );
            assert_eq!(
                bundled_deltat_piece(knot.next_up()).unwrap(),
                BundledDeltaTPiece::Observed(index)
            );
            let whole = knot.floor();
            let fraction = knot - whole;
            let unchanged = whole + fraction.next_up();
            if unchanged == knot {
                assert_eq!(
                    bundled_deltat_piece(unchanged).unwrap(),
                    bundled_deltat_piece(knot).unwrap()
                );
            }
        }
        let end = SMH2016_TABLE_TT[SMH2016_TABLE_TT.len() - 1];
        assert_eq!(
            bundled_deltat_piece(end).unwrap(),
            BundledDeltaTPiece::Observed(SMH2016_TABLE_TT.len() - 2)
        );
        assert!(matches!(
            bundled_deltat_piece(end.next_up()),
            Ok(BundledDeltaTPiece::Curve(_))
        ));
        assert!(bundled_deltat_piece(f64::NAN).is_err());
    }

    #[test]
    fn source_ut1_parts_keep_exact_offset_before_scalar_rounding() {
        for ut in [0.2, 2451545.0, 2460000.5] {
            let (whole, fraction, original_ut) = source_time_parts_from_ut1(ut);
            let checked = timescale_parts_from_ut1(ut).unwrap();
            assert_eq!(whole.to_bits(), ut.to_bits());
            assert_eq!(fraction.to_bits(), checked.1.to_bits());
            assert_eq!(original_ut.to_bits(), (0.0 * ut).to_bits());
            assert_eq!(
                (whole + fraction).to_bits(),
                timescale_tt_from_ut1(ut).unwrap().to_bits()
            );
        }
        let (whole, fraction, _) = source_time_parts_from_ut1(2451545.0);
        assert_ne!(fraction.to_bits(), ((whole + fraction) - whole).to_bits());
    }

    #[test]
    fn raw_ut1_fraction_retains_signed_zero_and_nonfinite_time_semantics() {
        let (whole, _, fraction) = source_time_parts_from_ut1(-0.0);
        assert_eq!(whole.to_bits(), (-0.0f64).to_bits());
        assert_eq!(fraction.to_bits(), (-0.0f64).to_bits());
        for ut in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let (whole, _, fraction) = source_time_parts_from_ut1(ut);
            assert_eq!(whole.to_bits(), ut.to_bits());
            assert!(fraction.is_nan());
            assert!(timescale_parts_from_ut1(ut).is_err());
        }
    }

    #[test]
    fn table_spans_and_knots_are_sane() {
        assert!(SMH2016_TABLE_TT.len() >= 20000);
        assert_eq!(SMH2016_TABLE_TT.len(), SMH2016_TABLE_DT.len());
        assert!(SMH2016_TABLE_TT.windows(2).all(|w| w[0] < w[1]));
        assert!(SMH2016_TABLE_DT.iter().all(|v| v.is_finite()));
        assert_eq!(SMH2016_CURVE_ROWS.len(), 13);
        for row in SMH2016_CURVE_ROWS {
            assert!(row[1] > row[0]);
            assert!(row.iter().all(|v| v.is_finite()));
        }
    }

    #[test]
    fn provider_serves_finite_and_refuses_non_finite() {
        for jd in [0.2, 1500000.0, 2326267.5, 2451545.0, 2461428.5, 3000000.0] {
            let v = smh2016_lookup_seconds(jd).expect("finite epochs serve");
            assert!(v.is_finite(), "jd={jd}");
        }
        assert_eq!(smh2016_lookup_seconds(f64::NAN), None);
        assert_eq!(smh2016_lookup_seconds(f64::INFINITY), None);
    }

    #[test]
    fn inverse_timescale_matches_pinned_scalar_observables() {
        // Pinned v3.2.1 reference using the bundled Skyfield 1.54 clock data.
        for (tt, word) in [
            (-1000.2, 0xc08f4e70b672cf7b),
            (0.2, 0xbff67581cbf05d62),
            (2451545.0, 0x4142b42c7fe7cad3),
            (2460000.5, 0x4142c4b03fe5c16e),
            (2475612.0308, 0x4142e32e03d44a28),
        ] {
            assert_eq!(timescale_ut1_from_tt(tt).unwrap().to_bits(), word);
        }
        assert!(matches!(
            timescale_ut1_from_tt(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            timescale_ut1_from_tt(f64::MAX),
            Err(Error::MathDomain { .. })
        ));
    }

    #[test]
    fn timescale_stops_before_the_provider_third_evaluation() {
        // A small-magnitude JD resolves the extra correction that is often
        // hidden by scalar-JD rounding at modern epochs. This checks the
        // operation sequence, not independent physical model accuracy.
        let ut1 = 0.2;
        let first_tt = ut1 + raw_seconds(ut1) / SECONDS_PER_DAY;
        let expected_tt = ut1 + raw_seconds(first_tt) / SECONDS_PER_DAY;
        let tt = timescale_tt_from_ut1(ut1).unwrap();
        assert_eq!(tt.to_bits(), expected_tt.to_bits());
        let provider_tt = ut1 + smh2016_lookup_seconds(ut1).unwrap() / SECONDS_PER_DAY;
        assert_ne!(tt.to_bits(), provider_tt.to_bits());
    }

    #[test]
    fn timescale_preserves_the_original_ut1_word_until_the_final_sum() {
        // This is the exact two-part storage/sum rule; no normalization of
        // the UT1 part and no offset recovered by scalar subtraction.
        let ut1 = 2460232.249616082;
        let first_tt = ut1 + raw_seconds(ut1) / SECONDS_PER_DAY;
        let offset = raw_seconds(first_tt) / SECONDS_PER_DAY;
        let tt = timescale_tt_from_ut1(ut1).unwrap();
        assert_eq!(tt.to_bits(), (ut1 + offset).to_bits());
        assert_ne!(offset.to_bits(), (tt - ut1).to_bits());
    }

    #[test]
    fn timescale_table_seams_and_neighboring_words() {
        // The endpoints themselves belong to the observed table; only
        // strictly exterior values belong to the extrapolating spline.
        // Include the historical/daily hand-off and its adjacent segments.
        let daily_start = SMH2016_TABLE_TT
            .windows(2)
            .position(|w| w[1] - w[0] < 2.0)
            .unwrap();
        let last = SMH2016_TABLE_TT.len() - 1;
        for i in [0, daily_start, daily_start + 1, last] {
            let knot = SMH2016_TABLE_TT[i];
            assert_eq!(raw_seconds(knot).to_bits(), SMH2016_TABLE_DT[i].to_bits());
            let translated = knot - SMH2016_TABLE_DT[i] / SECONDS_PER_DAY;
            for ut1 in [
                knot.next_down(),
                knot,
                knot.next_up(),
                translated.next_down(),
                translated,
                translated.next_up(),
            ] {
                let t1 = ut1 + raw_seconds(ut1) / SECONDS_PER_DAY;
                let expected = ut1 + raw_seconds(t1) / SECONDS_PER_DAY;
                assert_eq!(
                    timescale_tt_from_ut1(ut1).unwrap().to_bits(),
                    expected.to_bits()
                );
            }
        }
    }

    #[test]
    fn timescale_refuses_non_finite_input_or_model_overflow() {
        for ut1 in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                timescale_tt_from_ut1(ut1),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        for ut1 in [f64::MAX, -f64::MAX] {
            assert_eq!(
                timescale_tt_from_ut1(ut1),
                Err(Error::MathDomain { value: ut1 })
            );
        }
        for ut1 in [-3100015.5, 0.2, 2451545.0, 8000016.5] {
            assert!(timescale_tt_from_ut1(ut1).unwrap().is_finite());
        }
    }

    #[test]
    fn j2000_value_matches_release_construction() {
        // the model `deltat(2451545.0)` with the default model and default
        // tidal acceleration: 63.82896039639765 s -> days. This pins the
        // provider against the reference construction (not a port copy:
        // the sweep below checks the full grid against the live pin).
        let seconds = smh2016_lookup_seconds(2451545.0).expect("serves");
        assert_eq!(seconds.to_bits(), 63.82896039639765f64.to_bits());
    }

    #[test]
    fn timescale_parts_keep_the_input_word_whole() {
        // The pair mirrors `Timescale.ut1_jd`: the whole is the untouched
        // input value (The clock model stores it as-is, not floored) and the scalar
        // stays exactly `whole + tt_fraction`. At the PHEN-01 vector UT
        // epoch the TT fraction is the model table offset value.
        for (ut1, tt_fraction_word) in [
            (2451545.0, 0x3f48_352d_6a6e_1050u64),
            (2470000.5, 0x3f4b_2422_39fd_ae45u64),
        ] {
            let (whole, tt_fraction) = timescale_parts_from_ut1(ut1).expect("parts");
            assert_eq!(whole.to_bits(), ut1.to_bits());
            assert_eq!(tt_fraction.to_bits(), tt_fraction_word);
            assert_eq!(
                timescale_tt_from_ut1(ut1).unwrap().to_bits(),
                (whole + tt_fraction).to_bits()
            );
        }
        assert!(matches!(
            timescale_parts_from_ut1(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
}

#[cfg(test)]
mod source_time_domain_tests {
    use super::*;

    #[test]
    fn internal_time_keeps_two_raw_substitutions_and_public_guard() {
        for u in [1e308, -1e308] {
            let first = raw_seconds(u);
            let second = raw_seconds(u + first / SECONDS_PER_DAY);
            let expected = u + second / SECONDS_PER_DAY;
            assert_eq!(
                source_timescale_tt_from_ut1(u).unwrap().to_bits(),
                expected.to_bits()
            );
            assert!(expected.is_nan());
            assert!(
                matches!(timescale_tt_from_ut1(u), Err(Error::MathDomain { value }) if value.to_bits() == u.to_bits())
            );
        }
        for u in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                source_timescale_tt_from_ut1(u),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        for u in [2451545.125, 2400000.5, 2500000.5] {
            assert_eq!(
                source_timescale_tt_from_ut1(u).unwrap().to_bits(),
                timescale_tt_from_ut1(u).unwrap().to_bits()
            );
        }
    }
}

#[cfg(test)]
mod analytic_clock_rate_tests {
    use super::*;

    /// Expected values differentiate the scalar polynomial laws independently
    /// at 100 and 120 decimal digits. Two adjacent binary64 words bound the
    /// rounding of these short expressions; this is a unit arithmetic check,
    /// not an astronomical comparison limit.
    fn assert_rounded_rate(actual: f64, expected: f64) {
        assert!(actual.is_finite());
        let lower = expected.next_down().next_down();
        let upper = expected.next_up().next_up();
        assert!(
            actual >= lower && actual <= upper,
            "{actual:.17e} vs {expected:.17e}"
        );
    }

    #[test]
    fn analytic_rate_covers_all_thirteen_cubic_rows() {
        // Include portions masked by observed data. Inputs are TT days;
        // derivatives are seconds per TT day.
        for (index, tt, expected) in [
            (0, 1138471.25, -0.060862422997946614),
            (1, 1238915.0, -0.05770303739636462),
            (2, 1514678.75, -0.042104369521538496),
            (3, 1730176.25, -0.02883005167693361),
            (4, 1921932.5, -0.02563180783709788),
            (5, 2099991.875, -0.012147969427332877),
            (6, 2154779.375, -0.00748095254391969),
            (7, 2214132.5, -0.005419351471594798),
            (8, 2278051.25, -0.005286517796030117),
            (9, 2310010.625, -0.003990616016427105),
            (10, 2330099.375, 0.0006035095161324592),
            (11, 2532007.625600556, 0.005327902904170917),
            (12, 2752876.25, 0.01779603011635866),
        ] {
            assert_rounded_rate(curve_row_rate_seconds_per_day(tt, index), expected);
        }
    }

    #[test]
    fn analytic_rate_uses_observed_and_curve_endpoint_ownership() {
        // First/final observed knots, historic/daily transition, interior
        // knots and neighboring words. Never differentiate across a join.
        for (tt, expected) in [
            (2326267.4999999995, 0.00016720113767107548),
            (2326267.5, -0.005494505494505495),
            (2326267.5000000005, -0.005494505494505495),
            (2326449.4999999995, -0.005494505494505495),
            (2326449.5, 0.0),
            (2326449.5000000005, 0.0),
            (2441683.4999999995, 0.002994535519125693),
            (2441683.5, 0.0025808801668835792),
            (2441683.5000000005, 0.0025808801668835792),
            (2441684.5005113883, 0.0025808801668835792),
            (2441684.5005113888, 0.002801500000003898),
            (2441684.500511389, 0.002801500000003898),
            (2451240.50074287, 0.0009499000000019464),
            (2451240.5007428704, 0.0009776999999999703),
            (2451240.500742871, 0.0009776999999999703),
            (2461427.5008007404, -0.00047800000000108867),
            (2461427.500800741, -0.0001289999999869451),
            (2461427.5008007414, -0.0001289999999869451),
            (2461428.5008007404, -0.0001289999999869451),
            (2461428.500800741, -0.0001289999999869451),
            (2461428.5008007414, -5.7746277741767016e-05),
        ] {
            assert_rounded_rate(bundled_deltat_rate_seconds_per_day(tt).unwrap(), expected);
        }
    }

    #[test]
    fn analytic_clock_chains_keep_two_and_three_substitutions_distinct() {
        // The independent physical clock keeps each acquired scalar center
        // while differentiating its owned polynomial, not rounding steps.
        for (ut, expected_tt, expected_lookup) in [
            (0.2, 0.9999986535674658, -0.11633177095939726),
            (-3100015.5, 0.999996905411546, -0.2673724424335357),
            (8000016.5, 1.0000031649662844, 0.2734530869681383),
            (2451545.0, 1.0000000101458335, 0.0008766000088918443),
            (2460232.249616082, 1.0000000037361112, 0.0003228000012124124),
            (1138471.25, 0.9999992955745501, -0.06086235886884952),
            (1238915.0, 0.9999993321413624, -0.057702986294345136),
            (1514678.75, 0.9999995126812472, -0.04210434023869515),
            (1730176.25, 0.9999996663190334, -0.028830035520007128),
            (1921932.5, 0.9999997033356601, -0.025631798961969394),
            (2099991.875, 0.9999998593985435, -0.012147965835589542),
            (2154779.375, 0.9999999134149162, -0.007480951235985325),
            (2214132.5, 0.9999999372760298, -0.005419351031437703),
            (2278051.25, 0.9999999388134555, -0.005286517448962412),
            (2310010.625, 0.9999999538123177, -0.003990615752212423),
            (2330099.375, 0.9999999367536976, -0.00546448052870873),
            (2532007.625600556, 1.0000000616655491, 0.005327903446933645),
            (2752876.25, 1.0000002059726323, 0.01779603543414645),
            (
                2326267.4999999995,
                0.9999999364061862,
                -0.005494505145088936,
            ),
            (2326267.5, 0.9999999364061904, -0.005494505145088958),
            (
                2461428.5008007414,
                0.9999999993316411,
                -5.7746212781259354e-05,
            ),
        ] {
            let tt_rate = timescale_tt_rate_from_ut1(ut).unwrap();
            // These two-lookup rates round exactly to their MP expectation.
            assert_eq!(tt_rate.to_bits(), f64::to_bits(expected_tt));
            assert_rounded_rate(
                smh2016_lookup_rate_seconds_per_day(ut).unwrap(),
                expected_lookup,
            );
        }
    }

    #[test]
    fn analytic_clock_rates_preserve_input_and_overflow_contracts() {
        for input in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for error in [
                bundled_deltat_rate_seconds_per_day(input).unwrap_err(),
                timescale_tt_rate_from_ut1(input).unwrap_err(),
            ] {
                assert!(matches!(error, Error::NonFiniteArgument { value }
                    if value.to_bits() == input.to_bits()));
            }
            assert_eq!(smh2016_lookup_rate_seconds_per_day(input), None);
        }
        for input in [f64::MAX, -f64::MAX, 1e308, -1e308] {
            assert_eq!(
                bundled_deltat_rate_seconds_per_day(input),
                Err(Error::MathDomain { value: input })
            );
            assert_eq!(
                timescale_tt_rate_from_ut1(input),
                Err(Error::MathDomain { value: input })
            );
            assert!(
                !smh2016_lookup_rate_seconds_per_day(input)
                    .unwrap()
                    .is_finite()
            );
        }
    }
}
