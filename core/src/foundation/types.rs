//! Small explicit types for units, epochs and segments.
//!
//! These types describe numerical shapes, not time scales or physical units.
//! Callers must keep TT, UT1, TDB and UTC distinct and supply matching units to
//! each operation. A Julian day number by itself does not identify its scale.
//!
//! Ordered layers:
//! - [`NormalizedTime`]: the checked Chebyshev argument on the closed
//!   interval `[-1, 1]`, with exact `-1`/`1` endpoints at the segment ends.
//! - [`Segment`]: a validated `[begin, end]` epoch span with its linear map
//!   onto the Chebyshev interval (`normalize`, `dx_dt`); seam ownership
//!   between adjacent segments belongs to the reader layer, never here.
//! - [`TwoPartEpoch`]: additive `hi + lo` parts, with a checked split and
//!   compensated offsets to retain small corrections beside large epochs.
//!   Public fields also permit caller-created pairs; they enforce no invariant.
//! - [`segment_index`]: the owning segment of an epoch in a uniform grid,
//!   with half-open seam ownership verified by exact error-free comparison.
//! - [`round_half_away`]: nearest-integer rounding with ties away from zero.
//!
//! Frames, units, epochs: this module owns no frame computation. Units are
//! the caller's (Julian days for the segment bounds, dimensionless `x` for
//! the normalized argument, split days for the two-part epoch); the types
//! carry no unit tags. Epochs are consumed as supplied binary64 parts;
//! conversion between scales belongs to the time layer.
//!
//! Checked constructors and evaluators return typed errors for their stated
//! invalid inputs. [`Segment::normalize`] checks the epoch before clamping its
//! computed argument to `[-1, 1]`; it never clamps an out-of-range epoch.
//! [`Segment::new`] checks endpoint ordering, not span representability, so
//! finite endpoints can still have an overflowing difference. See each method
//! before using extreme values or unchecked public epoch pairs.
//!
//! [`segment_index`] additionally rejects nonfinite grid spans and grids with
//! more than 2^53 segments. It verifies seam ownership using compensated
//! comparisons; the global end belongs to the last segment.
//!

use crate::errors::{Error, Result};

/// Chebyshev argument restricted to the closed interval `[-1, 1]`.
///
/// LEB-style segments own epochs mapped onto this interval, with exact
/// endpoints `x = -1` (segment begin) and `x = 1` (segment end).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NormalizedTime(f64);

impl NormalizedTime {
    /// Validate `x`; rejects NaN, infinities and values outside `[-1, 1]`.
    pub fn new(x: f64) -> Result<Self> {
        if !x.is_finite() {
            return Err(Error::NonFiniteArgument { value: x });
        }
        if x < -1.0 || x > 1.0 {
            return Err(Error::NormalizedTimeOutOfRange { x });
        }
        Ok(Self(x))
    }

    /// Value without re-validation (for internally derived arguments).
    pub fn value(self) -> f64 {
        self.0
    }
}

/// Finite ordered endpoints with an inclusive map onto `[-1, 1]`.
///
/// [`Self::normalize`] admits both `begin` and `end`. A reader separately
/// decides which adjacent polynomial owns a shared seam; [`segment_index`]
/// provides half-open ownership with a closed global endpoint.
///
/// Construction does not guarantee a finite span or derivative scale. Callers
/// needing a valid linear map must ensure `end - begin` is finite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Segment {
    begin: f64,
    end: f64,
}

impl Segment {
    /// Validate finite bounds with `end > begin`, otherwise [`Error::InvalidSegment`].
    ///
    /// The subtraction is not checked: widely separated finite endpoints can
    /// have an infinite [`Self::length`].
    pub fn new(begin: f64, end: f64) -> Result<Self> {
        if !begin.is_finite() || !end.is_finite() || end <= begin {
            return Err(Error::InvalidSegment { begin, end });
        }
        Ok(Self { begin, end })
    }

    /// Segment begin, in the caller's epoch unit (e.g. Julian day).
    pub fn begin(self) -> f64 {
        self.begin
    }

    /// Segment end, same unit as [`Segment::begin`]; always above `begin`.
    pub fn end(self) -> f64 {
        self.end
    }

    /// Rounded `end - begin` in the caller's epoch units.
    /// This unchecked subtraction can overflow to positive infinity.
    pub fn length(self) -> f64 {
        self.end - self.begin
    }

    /// `dx/dt = 2 / (end - begin)`, in inverse epoch units.
    /// No finiteness check is made: a tiny span can give infinity and an
    /// overflowing span gives zero. These are not usable finite map scales.
    pub fn dx_dt(self) -> f64 {
        2.0 / self.length()
    }

    /// Map an epoch onto `[-1, 1]`; rejects epochs outside `[begin, end]`.
    ///
    /// For a finite span, endpoints map exactly (`begin -> -1`, `end -> 1`).
    /// The computed argument is clamped to absorb rounding spillover only after
    /// the inclusive epoch-range check. Nonfinite epochs return
    /// [`Error::NonFiniteArgument`]; uncovered epochs return
    /// [`Error::EpochOutOfSegment`].
    ///
    /// An overflowing span admitted by [`Self::new`] can produce a degenerate
    /// argument or NaN (which is refused); this method does not independently
    /// validate span representability. Seam ownership remains the reader's job.
    pub fn normalize(self, t: f64) -> Result<NormalizedTime> {
        if !t.is_finite() {
            return Err(Error::NonFiniteArgument { value: t });
        }
        if t < self.begin || t > self.end {
            return Err(Error::EpochOutOfSegment {
                t,
                begin: self.begin,
                end: self.end,
            });
        }
        // For a finite span, divide before multiplying to avoid overflowing
        // twice the offset. Clamp rounding spillover, not the input epoch.
        // Endpoint exactness requires a representable span; new checks only
        // endpoint finiteness and ordering.
        let x = ((t - self.begin) / self.length()) * 2.0 - 1.0;
        NormalizedTime::new(x.clamp(-1.0, 1.0))
    }
}

/// Round-half-away-from-zero, matching Rust `f64::round` semantics.
///
/// Exact half-integer ties round away from zero. Other finite values use the
/// nearest integer; signed zero, infinities and NaNs follow `f64::round`.
/// This operation does not apply the split-epoch arithmetic of [`TwoPartEpoch`].
pub fn round_half_away(x: f64) -> f64 {
    x.round()
}

/// Additive epoch parts `hi + lo` for retaining sub-ulp corrections.
///
/// The parts use the same units and time scale; neither field carries a unit
/// tag. They need not be an integer day and a fraction in `[0, 1)`.
/// [`Self::split`] truncates the leading significand and retains its remainder;
/// [`Self::add_offset`] uses Knuth TwoSum compensation. The latter can change
/// the split layout and does not promise exact arbitrary sequences of sums.
/// Keep parts separate until a consuming API explicitly needs [`Self::value`].
///
/// Public fields allow arbitrary caller-created pairs, including nonfinite or
/// overlapping parts. Such pairs have not passed constructor validation.
/// Splitting an already rounded scalar cannot recover information lost earlier.
///
/// ```
/// use libephemeris_core::types::TwoPartEpoch;
///
/// let base = TwoPartEpoch::split(2_451_545.0)?;
/// let shifted = base.add_offset(1e-10)?; // days, smaller than one JD ulp
/// assert_eq!(shifted.value(), base.value()); // scalar collapse loses the offset
/// assert!((shifted.diff(base) - 1e-10).abs() < 1e-15);
/// # Ok::<(), libephemeris_core::errors::Error>(())
/// ```
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TwoPartEpoch {
    /// Leading additive part; a 26-bit significand after [`Self::split`]
    /// for normal finite inputs, but not necessarily after offset arithmetic.
    pub hi: f64,
    /// Trailing additive part; [`Self::split`] supplies the exact remainder.
    /// Caller construction and later arithmetic need not retain that layout.
    pub lo: f64,
}

impl TwoPartEpoch {
    /// Exact split of a finite value into leading part + remainder.
    ///
    /// The leading part keeps the sign, exponent, and top 25 explicit
    /// mantissa bits (26-bit significand); the remainder is the exact
    /// difference. Nonzero inputs are reproduced by the collapsed sum.
    /// For either signed zero, `hi` keeps the input sign and `lo` is positive
    /// zero; collapsing `-0.0 + 0.0` need not preserve the negative-zero bit.
    /// Bit truncation cannot overflow, unlike multiply-based Dekker splitting
    /// near `f64::MAX`. Nonfinite inputs return [`Error::NonFiniteArgument`].
    pub fn split(value: f64) -> Result<Self> {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
        if value == 0.0 {
            return Ok(Self { hi: value, lo: 0.0 });
        }
        const KEEP: u64 = !((1 << 27) - 1);
        let hi = f64::from_bits(value.to_bits() & KEEP);
        // Exact: hi shares value's top 26 significand bits, so the
        // difference needs at most 27 bits and is always representable.
        let lo = value - hi;
        Ok(Self { hi, lo })
    }

    /// Zero epoch.
    pub fn zero() -> Self {
        Self { hi: 0.0, lo: 0.0 }
    }

    /// Rounded scalar `hi + lo`, discarding any sub-ulp remainder.
    /// This does not validate caller-created parts and can return infinity or
    /// NaN; it can also lose a stored negative-zero sign.
    pub fn value(self) -> f64 {
        self.hi + self.lo
    }

    /// Error-free sum of two parts: returns `(sum, error)` with
    /// `sum + error == a + b` exactly (Knuth TwoSum, branchless form valid
    /// for all finite inputs without ordering assumptions, provided no
    /// intermediate overflows; callers working near `f64::MAX` must check
    /// finiteness of the result).
    fn two_sum(a: f64, b: f64) -> (f64, f64) {
        let s = a + b;
        let bp = s - a;
        let ap = s - bp;
        let da = a - ap;
        let db = b - bp;
        (s, da + db)
    }

    /// Add a finite offset, keeping the correction in `lo` where possible.
    ///
    /// Uses TwoSum compensation without first collapsing the input pair.
    /// Returns [`Error::NonFiniteArgument`] for a nonfinite offset or any
    /// nonfinite intermediate, including invalid caller-created parts. The
    /// trailing-part addition still rounds; no universal exact-round-trip
    /// guarantee is made for arbitrary pairs or repeated offsets.
    pub fn add_offset(self, offset: f64) -> Result<Self> {
        if !offset.is_finite() {
            return Err(Error::NonFiniteArgument { value: offset });
        }
        let (s, e) = Self::two_sum(self.hi, offset);
        let lo = self.lo + e;
        let (hi, rem) = Self::two_sum(s, lo);
        if !s.is_finite() || !lo.is_finite() || !hi.is_finite() || !rem.is_finite() {
            return Err(Error::NonFiniteArgument {
                value: s + lo + hi + rem,
            });
        }
        Ok(Self { hi, lo: rem })
    }

    /// Rounded difference `self - other`, in the parts' common epoch units.
    /// Leading differences precede trailing corrections. This is not an
    /// error-free subtraction; nonfinite inputs and overflow propagate as IEEE
    /// values rather than returning an error.
    pub fn diff(self, other: Self) -> f64 {
        (self.hi - other.hi) + (self.lo - other.lo)
    }

    /// Error-free product: returns `(product, error)` with
    /// `product + error == a * b` exactly (Dekker), or `None` if any
    /// intermediate overflows. Inputs must be finite.
    fn two_product(a: f64, b: f64) -> Option<(f64, f64)> {
        let p = a * b;
        // Bit-truncation split is infallible on finite inputs.
        let ah = f64::from_bits(a.to_bits() & !((1 << 27) - 1));
        let bh = f64::from_bits(b.to_bits() & !((1 << 27) - 1));
        let al = a - ah;
        let bl = b - bh;
        let e = ((ah * bh - p) + ah * bl + al * bh) + al * bl;
        if p.is_finite() && e.is_finite() {
            Some((p, e))
        } else {
            None
        }
    }
}

/// Exact sign of `(t - begin) - q * length` for finite inputs.
///
/// Both differences are evaluated with error-free transformations, so the
/// returned ordering is the true real ordering, not a rounded one. Returns
/// `None` when an intermediate overflows: only possible for grids spanning
/// a sizable fraction of the `f64` range near its extremes, far outside
/// any astronomical domain.
fn cmp_to_seam(t: f64, begin: f64, q: f64, length: f64) -> Option<std::cmp::Ordering> {
    let (ah, al) = TwoPartEpoch::two_sum(t, -begin);
    let (ph, pl) = TwoPartEpoch::two_product(q, length)?;
    if !ah.is_finite() || !al.is_finite() || !ph.is_finite() || !pl.is_finite() {
        return None;
    }
    // Both expansions are non-overlapping (error term below half an ulp of
    // the leader), so the leaders decide unless they tie; float comparison
    // of the tied case preserves the exact sign.
    match ah.partial_cmp(&ph) {
        Some(std::cmp::Ordering::Equal) => al.partial_cmp(&pl),
        ord => ord,
    }
}

/// Locate the owning segment of `t` in a uniform grid.
///
/// Grid: `n` segments of `length` covering `[begin, begin + n * length]`.
/// Ownership is half-open per segment: an epoch exactly on an interior seam
/// belongs to the segment it starts; the global end belongs to the last
/// segment. Epochs outside the global span are an error, never silently
/// clamped. Adjacent-`f64`-word neighbors of a seam therefore resolve to
/// different segments, which the tests pin down.
///
/// This matches the release contract (`int` truncation of the quotient
/// for in-range epochs, global end clamped to the last segment). The float
/// quotient only seeds the search: ownership is then verified and adjusted
/// with exact error-free comparisons, so translated grids (nonzero `begin`,
/// non-power-of-two `length`) resolve seams exactly. Grids whose
/// span or end is not finite, or with more than 2^53 segments, are
/// rejected: ownership is unrepresentable there.
pub fn segment_index(t: f64, begin: f64, length: f64, n: usize) -> Result<usize> {
    use std::cmp::Ordering;
    if !t.is_finite() {
        return Err(Error::NonFiniteArgument { value: t });
    }
    if !begin.is_finite() || !length.is_finite() || length <= 0.0 || n == 0 {
        return Err(Error::InvalidSegment { begin, end: length });
    }
    if n > (1u64 << 53) as usize {
        return Err(Error::InvalidSegment { begin, end: length });
    }
    let span = length * n as f64;
    let end = begin + span;
    if !span.is_finite() || !end.is_finite() {
        return Err(Error::InvalidSegment { begin, end: length });
    }
    if t < begin || t > end {
        return Err(Error::EpochOutOfSegment { t, begin, end });
    }
    if t == end {
        return Ok(n - 1);
    }
    let unresolvable = || Error::SegmentOwnershipUnresolved { t, begin, length };
    // Seed from the float quotient, clamped exactly as the release clamps.
    let mut q = ((t - begin) / length).floor() as usize;
    q = q.min(n - 1);
    // Step down while `t` is strictly below its seam. Terminates: `q`
    // strictly decreases toward 0.
    while let Ordering::Less = cmp_to_seam(t, begin, q as f64, length).ok_or_else(unresolvable)? {
        if q == 0 {
            // Unreachable: t >= begin == seam(0) exactly, so t
            // cannot be strictly below it. Kept as an error rather
            // than a silent assignment.
            return Err(unresolvable());
        }
        q -= 1;
    }
    // Step up while `t` reaches the next seam. Terminates: seams are
    // strictly increasing (`length > 0`), bounded by `n`.
    loop {
        match cmp_to_seam(t, begin, (q + 1) as f64, length).ok_or_else(unresolvable)? {
            Ordering::Less => return Ok(q),
            _ => {
                q += 1;
                if q >= n {
                    // t >= seam(n): float `end` and seam(n) can differ by
                    // rounding; the release clamps such epochs to the last
                    // segment, and so do we.
                    return Ok(n - 1);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normalized_time_accepts_endpoints() {
        assert_eq!(NormalizedTime::new(-1.0).unwrap().value(), -1.0);
        assert_eq!(NormalizedTime::new(1.0).unwrap().value(), 1.0);
        assert_eq!(NormalizedTime::new(0.0).unwrap().value(), 0.0);
    }

    #[test]
    fn normalized_time_rejects_out_of_range_and_non_finite() {
        assert!(matches!(
            NormalizedTime::new(1.5),
            Err(Error::NormalizedTimeOutOfRange { .. })
        ));
        assert!(matches!(
            NormalizedTime::new(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(matches!(
            NormalizedTime::new(f64::INFINITY),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn segment_endpoints_map_exactly() {
        let seg = Segment::new(10.0, 14.0).unwrap();
        assert_eq!(seg.normalize(10.0).unwrap().value(), -1.0);
        assert_eq!(seg.normalize(14.0).unwrap().value(), 1.0);
        assert_eq!(seg.normalize(12.0).unwrap().value(), 0.0);
        assert_eq!(seg.dx_dt(), 0.5);
    }

    #[test]
    fn segment_rejects_degenerate_bounds_and_outside_epochs() {
        assert!(matches!(
            Segment::new(2.0, 2.0),
            Err(Error::InvalidSegment { .. })
        ));
        assert!(matches!(
            Segment::new(2.0, 1.0),
            Err(Error::InvalidSegment { .. })
        ));
        let seg = Segment::new(0.0, 2.0).unwrap();
        assert!(matches!(
            seg.normalize(2.5),
            Err(Error::EpochOutOfSegment { .. })
        ));
    }

    #[test]
    fn split_is_bitwise_exact() {
        // hi + lo reproduces the input with identical bits, including a
        // large Julian-day-like epoch with a fractional tail.
        for t in [0.0, 1.5, 2451545.0, 2451545.0000001, -123.456e10, 1e-300] {
            let e = TwoPartEpoch::split(t).unwrap();
            assert_eq!(
                (e.hi + e.lo).to_bits(),
                t.to_bits(),
                "split not exact for {t}"
            );
            // Leading part carries the magnitude; truncation keeps the top
            // 26 significand bits, so the remainder holds at most the low
            // 27 bits (~2^-25 of |t| for nonzero t; rounding-based splits
            // achieve 2^-26, truncation honestly needs one bit more).
            if t != 0.0 {
                assert!(
                    e.lo.abs() <= 2.0_f64.powi(-25) * t.abs(),
                    "split remainder too large for {t}"
                );
            }
        }
        // Extreme magnitudes no longer overflow (FIX-03): the split is
        // pure bit truncation, exact for every finite input.
        for t in [f64::MAX, f64::MIN, f64::MIN_POSITIVE, 1e300] {
            let e = TwoPartEpoch::split(t).unwrap();
            assert!(e.hi.is_finite() && e.lo.is_finite());
            assert_eq!((e.hi + e.lo).to_bits(), t.to_bits());
        }
        assert!(matches!(
            TwoPartEpoch::split(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn small_offsets_survive_large_epochs() {
        // ulp(2451545) ~= 9.3e-10, so a 1e-10 offset is lost in plain f64
        // addition but kept in two parts.
        let base = TwoPartEpoch::split(2451545.0).unwrap();
        assert_eq!(2451545.0_f64 + 1e-10, 2451545.0);
        let shifted = base.add_offset(1e-10).unwrap();
        assert!((shifted.value() - 2451545.0 - 1e-10).abs() < 1e-9);
        // diff recovers the offset to ~1 ulp of the correction scale.
        assert!((shifted.diff(base) - 1e-10).abs() < 1e-15);
    }

    #[test]
    fn segment_index_translated_grid_seams_are_exact() {
        // Review reproducer (FIX-01): on grid [-2, 4) with length 3 the
        // epoch just below the seam at 1 belongs to segment 0. The old
        // float quotient rounded the displacement 3 - eps up to 3.
        assert_eq!(segment_index(1.0_f64.next_down(), -2.0, 3.0, 2).unwrap(), 0);
        assert_eq!(segment_index(1.0, -2.0, 3.0, 2).unwrap(), 1);
        assert_eq!(segment_index(1.0_f64.next_up(), -2.0, 3.0, 2).unwrap(), 1);
        // Translated Julian-day-like grid: seams at non-representable
        // positions resolve exactly on both sides.
        let (begin, length, n) = (2451545.0, 0.125, 8);
        for k in 1..n {
            // Nearest floats straddling the mathematical seam.
            let seam_approx = begin + k as f64 * length;
            for t in [seam_approx.next_down(), seam_approx.next_up()] {
                if t < begin || t >= begin + n as f64 * length {
                    continue;
                }
                let got = segment_index(t, begin, length, n).unwrap();
                // Verify against exact seam comparisons, not the
                // implementation under test.
                let below =
                    cmp_to_seam(t, begin, k as f64, length) == Some(std::cmp::Ordering::Less);
                assert_eq!(got, if below { k - 1 } else { k }, "t={t}");
            }
        }
        // Agreement with the release's truncation rule away from hazards:
        // where the float quotient is more than 4 ulps from an integer,
        // exact ownership must equal the naive floor.
        let mut x = 1234567usize;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x as f64 / u64::MAX as f64
        };
        for _ in 0..500 {
            let b = (next() - 0.5) * 2e7;
            let l = 0.01 + next() * 100.0;
            let nn = 1 + (next() * 20.0) as usize;
            let t = b + next() * l * nn as f64;
            if t >= b + l * nn as f64 {
                continue;
            }
            let q = (t - b) / l;
            if (q - q.floor()).min(q.ceil() - q) * l.max(1.0)
                < 4.0 * f64::EPSILON * (t - b).abs().max(1.0)
            {
                continue; // too close to a seam to call naive exact
            }
            assert_eq!(
                segment_index(t, b, l, nn).unwrap(),
                q.floor() as usize,
                "t={t} b={b} l={l} n={nn}"
            );
        }
    }

    #[test]
    fn segment_index_seam_ownership_and_adjacent_words() {
        // Grid [0, 30), 3 segments of length 10.
        assert_eq!(segment_index(0.0, 0.0, 10.0, 3).unwrap(), 0);
        assert_eq!(segment_index(9.999, 0.0, 10.0, 3).unwrap(), 0);
        // Exact seam belongs to the segment it starts.
        assert_eq!(segment_index(10.0, 0.0, 10.0, 3).unwrap(), 1);
        assert_eq!(segment_index(20.0, 0.0, 10.0, 3).unwrap(), 2);
        // Adjacent f64 words resolve to different segments.
        assert_eq!(
            segment_index(10.0_f64.next_down(), 0.0, 10.0, 3).unwrap(),
            0
        );
        assert_eq!(segment_index(10.0_f64.next_up(), 0.0, 10.0, 3).unwrap(), 1);
        // Global end belongs to the last segment; outside is an error.
        assert_eq!(segment_index(30.0, 0.0, 10.0, 3).unwrap(), 2);
        assert!(matches!(
            segment_index(30.5, 0.0, 10.0, 3),
            Err(Error::EpochOutOfSegment { .. })
        ));
        assert!(matches!(
            segment_index(-0.5, 0.0, 10.0, 3),
            Err(Error::EpochOutOfSegment { .. })
        ));
    }

    #[test]
    fn rounding_documents_half_away_rule() {
        assert_eq!(round_half_away(2.5), 3.0);
        assert_eq!(round_half_away(-2.5), -3.0);
        assert_eq!(round_half_away(2.4), 2.0);
    }
}
