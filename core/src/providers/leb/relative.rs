//! Relative Cartesian states from aligned, owned Chebyshev polynomials.

use super::*;

/// ICRS Cartesian position in AU and its analytical TT-day derivative in AU/day.
pub type CartesianState = ([f64; 3], [f64; 3]);

/// Inline capacity is an output-storage choice, never a polynomial degree limit.
const INLINE_TARGET_COEFFICIENTS: usize = 64;

enum RelativeTargetRow<'a> {
    Inline(&'a [f64]),
    Owned(CoefficientWords),
}

impl RelativeTargetRow<'_> {
    fn as_slice(&self) -> &[f64] {
        match self {
            Self::Inline(words) => words,
            Self::Owned(words) => words,
        }
    }
}

impl FileIndex {
    /// Consume one selected row while retaining its actual coefficient owner.
    /// Cached LEB2 words are borrowed without a copied-row allocation. LEB1
    /// keeps its ordinary owned row and resource reservation. No cache lock is
    /// held across the callback; an evicted decoded owner stays live until its
    /// last callback returns, with its real working-memory lease intact.
    #[allow(clippy::too_many_arguments)]
    fn with_relative_segment_coefficients<R>(
        &self,
        bytes: LebSource<'_>,
        body: i32,
        entry: &BodyEntry,
        segment: usize,
        jd: f64,
        cache: &LebChunkCache,
        consume: impl FnOnce(&[f64]) -> Result<R>,
    ) -> Result<R> {
        if self.header().format() == LebFormat::Leb1 {
            let coefficients = read_segment_coefficients(bytes, entry, segment)?;
            consume(&coefficients)
        } else {
            self.with_leb2_segment_coefficients(
                bytes,
                body,
                entry,
                segment,
                jd,
                Some(cache),
                consume,
            )
        }
    }

    /// Release the target chunk before acquiring the center. Small selected
    /// rows copy only their words to caller-owned stack storage; larger rows
    /// retain the ordinary copied-row heap allocation and actual live lease.
    /// LEB1 already owns its selected row and needs no additional copy. No
    /// decoded chunk owner or cache borrow escapes this synchronous operation.
    #[allow(clippy::too_many_arguments)]
    fn relative_target_row<'a>(
        &self,
        bytes: LebSource<'_>,
        body: i32,
        entry: &BodyEntry,
        segment: usize,
        jd: f64,
        cache: &LebChunkCache,
        inline: &'a mut [f64; INLINE_TARGET_COEFFICIENTS],
    ) -> Result<RelativeTargetRow<'a>> {
        if self.header().format() == LebFormat::Leb1 {
            return read_segment_coefficients(bytes, entry, segment).map(RelativeTargetRow::Owned);
        }
        let fits_inline = (entry.degree() as usize)
            .checked_add(1)
            .and_then(|width| width.checked_mul(entry.components() as usize))
            .is_some_and(|words| words <= inline.len());
        if !fits_inline {
            return self
                .leb2_segment_coefficients(bytes, body, entry, segment, jd, Some(cache))
                .map(RelativeTargetRow::Owned);
        }
        let length = self.with_leb2_segment_coefficients(
            bytes,
            body,
            entry,
            segment,
            jd,
            Some(cache),
            |words| {
                inline[..words.len()].copy_from_slice(words);
                Ok(words.len())
            },
        )?;
        Ok(RelativeTargetRow::Inline(&inline[..length]))
    }

    /// Each epoch is validated just before its own row is acquired, retaining
    /// target-payload failures before any later center epoch or payload error.
    fn relative_record(body: i32, entry: &BodyEntry, epoch: TwoPartEpoch) -> Result<(usize, f64)> {
        let jd = epoch.value();
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start() || jd > entry.jd_end() {
            return Err(Error::LebOutOfRange {
                body_id: body,
                jd,
                start: entry.jd_start(),
                end: entry.jd_end(),
            });
        }
        let (segment, tau) = select_segment_two_part(
            entry.jd_start(),
            entry.interval_days(),
            entry.segment_count(),
            epoch,
        )?;
        Ok((segment, tau.value()))
    }

    /// Evaluate target minus center before rounding their barycentric values.
    ///
    /// `jd_tt` is TT. Both selected records must be Cartesian ICRS body-center
    /// channels with identical polynomial degree, interval, origin and segment
    /// ownership. Their coefficient difference is evaluated by Clenshaw, and
    /// differentiated analytically. This avoids cancellation of two rounded
    /// barycentric positions. No polynomial is refitted or extrapolated.
    ///
    /// Returns `None` when those alignment conditions do not hold, allowing
    /// the caller's ordinary independently selected state calculation. Missing
    /// bodies, invalid metadata and coefficient decoding errors propagate.
    /// `cache` belongs to this image; it uses the ordinary body/chunk keys.
    pub(crate) fn eval_relative_icrs_cached_source(
        &self,
        bytes: LebSource<'_>,
        target: i32,
        center: i32,
        jd_tt: f64,
        cache: &LebChunkCache,
    ) -> Result<Option<CartesianState>> {
        let target_entry = self.find_body(target)?;
        let center_entry = self.find_body(center)?;
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        if target_entry.coord_type() != COORD_ICRS_BARY
            || center_entry.coord_type() != COORD_ICRS_BARY
            || target_entry.components() != 3
            || center_entry.components() != 3
            || target_entry.degree() != center_entry.degree()
            || target_entry.interval_days() != center_entry.interval_days()
            || target_entry.jd_start() != center_entry.jd_start()
        {
            return Ok(None);
        }
        let domain = target_entry.smooth_interval(jd_tt)?;
        if domain.is_none() || domain != center_entry.smooth_interval(jd_tt)? {
            return Ok(None);
        }
        let (segment, tau) = select_segment(
            target_entry.jd_start(),
            target_entry.interval_days(),
            target_entry.segment_count(),
            jd_tt,
        )?;
        let mut inline = [0.0; INLINE_TARGET_COEFFICIENTS];
        let target_row = self.relative_target_row(
            bytes,
            target,
            target_entry,
            segment,
            jd_tt,
            cache,
            &mut inline,
        )?;
        let target_coefficients = target_row.as_slice();
        self.with_relative_segment_coefficients(
            bytes,
            center,
            center_entry,
            segment,
            jd_tt,
            cache,
            |center_coefficients| {
                let difference_bytes = target_coefficients
                    .len()
                    .checked_mul(std::mem::size_of::<f64>())
                    .ok_or_else(|| Error::LebCorrupt {
                        what: "relative coefficient row size".into(),
                        offset: 0,
                    })?;
                let _difference_allocation = bytes.reserve_work(difference_bytes)?;
                let difference: Vec<_> = target_coefficients
                    .iter()
                    .copied()
                    .zip(center_coefficients.iter().copied())
                    .map(|(target, center)| target - center)
                    .collect();
                eval_one_segment_inline(target_entry, &difference, tau.value()).map(Some)
            },
        )
    }
}

/// Evaluate two Cartesian polynomials jointly without rounding their large
/// common origin first. The two normalized arguments and degrees may differ;
/// each polynomial remains in its own selected record and domain. Paired
/// Clenshaw recurrences propagate their difference, including the argument
/// displacement, while a second recurrence carries the center polynomial.
fn polynomial_difference(
    target: &[f64],
    target_tau: f64,
    center: &[f64],
    center_tau: f64,
) -> Result<f64> {
    let count = target.len().max(center.len());
    let (mut c1, mut c2, mut d1, mut d2) = (0.0, 0.0, 0.0, 0.0);
    let coefficient = |values: &[f64], index| values.get(index).copied().unwrap_or(0.0);
    for index in (1..count).rev() {
        let a = coefficient(target, index);
        let b = coefficient(center, index);
        if !a.is_finite() || !b.is_finite() {
            return Err(Error::NonFiniteArgument { value: a + b });
        }
        let difference = (2.0 * target_tau).mul_add(
            d1,
            (2.0 * (target_tau - center_tau)).mul_add(c1, (a - b) - d2),
        );
        let value = (2.0 * center_tau).mul_add(c1, b - c2);
        (c2, c1, d2, d1) = (c1, value, d1, difference);
    }
    let a = coefficient(target, 0);
    let b = coefficient(center, 0);
    let difference = target_tau.mul_add(d1, (target_tau - center_tau).mul_add(c1, (a - b) - d2));
    if difference.is_finite() {
        Ok(difference)
    } else {
        Err(Error::NonFiniteArgument { value: difference })
    }
}

impl FileIndex {
    /// Target position at `target_epoch` minus center at `center_epoch`.
    ///
    /// Both epochs are two-part TT Julian days and the result is ICRS AU.
    /// The selected channels must be Cartesian body centers; system-barycenter
    /// channels return `None` so their actual center corrections remain owned
    /// by the ordinary provider. Each polynomial is decoded at its own epoch
    /// and evaluated by a paired Clenshaw recurrence. Neither polynomial is
    /// refitted or extrapolated, including when the epochs straddle a join.
    /// Missing bodies, out-of-range/nonfinite epochs and malformed data fail
    /// with the ordinary reader errors. The cache belongs to this image.
    pub(crate) fn eval_relative_position_two_part_cached_source(
        &self,
        bytes: LebSource<'_>,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
        cache: &LebChunkCache,
    ) -> Result<Option<[f64; 3]>> {
        let target_entry = self.find_body(target)?;
        let center_entry = self.find_body(center)?;
        if [target_entry, center_entry]
            .iter()
            .any(|entry| entry.coord_type() != COORD_ICRS_BARY || entry.components() != 3)
        {
            return Ok(None);
        }
        let (target_segment, target_tau) =
            Self::relative_record(target, target_entry, target_epoch)?;
        let mut inline = [0.0; INLINE_TARGET_COEFFICIENTS];
        let target_row = self.relative_target_row(
            bytes,
            target,
            target_entry,
            target_segment,
            target_epoch.value(),
            cache,
            &mut inline,
        )?;
        let target_coefficients = target_row.as_slice();
        let (center_segment, center_tau) =
            Self::relative_record(center, center_entry, center_epoch)?;
        self.with_relative_segment_coefficients(
            bytes,
            center,
            center_entry,
            center_segment,
            center_epoch.value(),
            cache,
            |center_coefficients| {
                let target_count = target_coefficients.len() / 3;
                let center_count = center_coefficients.len() / 3;
                let mut position = [0.0; 3];
                for (component, value) in position.iter_mut().enumerate() {
                    *value = polynomial_difference(
                        &target_coefficients
                            [component * target_count..(component + 1) * target_count],
                        target_tau,
                        &center_coefficients
                            [component * center_count..(component + 1) * center_count],
                        center_tau,
                    )?;
                }
                Ok(Some(position))
            },
        )
    }

    /// Joint derivative of independently selected Cartesian polynomials.
    ///
    /// Both epochs are two-part TT Julian days; the result is ICRS AU/TT-day.
    /// The paired recurrence retains coefficient and normalization residuals.
    /// Channel admission, target-before-center decoding and ordinary errors
    /// follow the relative-position route; no record is extrapolated.
    pub(crate) fn eval_relative_velocity_two_part_cached_source(
        &self,
        bytes: LebSource<'_>,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
        cache: &LebChunkCache,
    ) -> Result<Option<[f64; 3]>> {
        let target_entry = self.find_body(target)?;
        let center_entry = self.find_body(center)?;
        if [target_entry, center_entry]
            .iter()
            .any(|entry| entry.coord_type() != COORD_ICRS_BARY || entry.components() != 3)
        {
            return Ok(None);
        }
        let (target_segment, target_tau) =
            Self::relative_record(target, target_entry, target_epoch)?;
        let mut inline = [0.0; INLINE_TARGET_COEFFICIENTS];
        let target_row = self.relative_target_row(
            bytes,
            target,
            target_entry,
            target_segment,
            target_epoch.value(),
            cache,
            &mut inline,
        )?;
        let target_coefficients = target_row.as_slice();
        let (center_segment, center_tau) =
            Self::relative_record(center, center_entry, center_epoch)?;
        self.with_relative_segment_coefficients(
            bytes,
            center,
            center_entry,
            center_segment,
            center_epoch.value(),
            cache,
            |center_coefficients| {
                let target_count = target_coefficients.len() / 3;
                let center_count = center_coefficients.len() / 3;
                let mut position = [0.0; 3];
                for (component, value) in position.iter_mut().enumerate() {
                    *value = polynomial_velocity_difference(
                        &target_coefficients
                            [component * target_count..(component + 1) * target_count],
                        target_tau,
                        2.0 / target_entry.interval_days(),
                        &center_coefficients
                            [component * center_count..(component + 1) * center_count],
                        center_tau,
                        2.0 / center_entry.interval_days(),
                    )?;
                }
                Ok(Some(position))
            },
        )
    }
}

impl CompositeIndex<'_> {
    /// ICRS AU target-minus-center position at two independently owned TT epochs.
    /// Both channels must use this same selected image; otherwise returns `None`.
    /// Errors and recurrence follow [`FileIndex::eval_relative_position_two_part_cached`].
    pub fn eval_relative_position_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        let a = self
            .serving_index(target)
            .ok_or(Error::LebBodyNotFound { body_id: target })?;
        let b = self
            .serving_index(center)
            .ok_or(Error::LebBodyNotFound { body_id: center })?;
        if a != b {
            return Ok(None);
        }
        self.evaluation(a)
            .relative_position(target, target_epoch, center, center_epoch)
    }

    /// Joint ICRS AU/TT-day velocity at two independently selected TT epochs.
    /// Both channels must share their selected image/tier; otherwise returns
    /// `None`. Missing bodies and selected reader/domain errors propagate.
    pub fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        let a = self
            .serving_index(target)
            .ok_or(Error::LebBodyNotFound { body_id: target })?;
        let b = self
            .serving_index(center)
            .ok_or(Error::LebBodyNotFound { body_id: center })?;
        if a != b {
            return Ok(None);
        }
        self.evaluation(a)
            .relative_velocity(target, target_epoch, center, center_epoch)
    }

    /// TT relative state from aligned ICRS channels in the same selected image.
    ///
    /// Image precedence remains that of `eval_body_any`. Different selected
    /// images return `None`; absent bodies and malformed selected data fail.
    /// Units and polynomial alignment follow [`FileIndex::eval_relative_icrs_cached`].
    pub fn eval_relative_icrs(
        &self,
        target: i32,
        center: i32,
        jd_tt: f64,
    ) -> Result<Option<CartesianState>> {
        let target_image = self
            .serving_index(target)
            .ok_or(Error::LebBodyNotFound { body_id: target })?;
        let center_image = self
            .serving_index(center)
            .ok_or(Error::LebBodyNotFound { body_id: center })?;
        if target_image != center_image {
            return Ok(None);
        }
        self.evaluation(target_image)
            .relative_state(target, center, jd_tt)
    }
}

impl TieredIndex<'_> {
    /// Two-epoch ICRS AU difference retaining each body's actual tier choice.
    /// Different tiers/images return `None`; reader/domain errors propagate.
    pub fn eval_relative_position_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        match (
            self.selected(target, target_epoch.value()),
            self.selected(center, center_epoch.value()),
        ) {
            (Some(a), Some(b)) if a == b => self.tiers[a].1.eval_relative_position_two_part(
                target,
                target_epoch,
                center,
                center_epoch,
            ),
            _ => Ok(None),
        }
    }

    /// Joint ICRS AU/TT-day velocity at two independently selected TT epochs.
    /// Both channels must share their selected image/tier; otherwise returns
    /// `None`. Missing bodies and selected reader/domain errors propagate.
    pub fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        match (
            self.selected(target, target_epoch.value()),
            self.selected(center, center_epoch.value()),
        ) {
            (Some(a), Some(b)) if a == b => self.tiers[a].1.eval_relative_velocity_two_part(
                target,
                target_epoch,
                center,
                center_epoch,
            ),
            _ => Ok(None),
        }
    }

    /// TT relative state using the independently selected tier of each body.
    ///
    /// Only two channels in the same selected tier/image may share polynomial
    /// subtraction. Other pairs return `None` without changing selection.
    /// Units and errors follow [`CompositeIndex::eval_relative_icrs`].
    pub fn eval_relative_icrs(
        &self,
        target: i32,
        center: i32,
        jd_tt: f64,
    ) -> Result<Option<CartesianState>> {
        let target_tier = self.selected(target, jd_tt);
        let center_tier = self.selected(center, jd_tt);
        match (target_tier, center_tier) {
            (Some(target_tier), Some(center_tier)) if target_tier == center_tier => self.tiers
                [target_tier]
                .1
                .eval_relative_icrs(target, center, jd_tt),
            _ => Ok(None),
        }
    }
}

/// Two binary64 words used only by the paired velocity recurrence. Error-free
/// sums/products retain the residual until the final relative rate is rounded.
#[derive(Clone, Copy)]
struct VelocityExpansion {
    hi: f64,
    lo: f64,
}

impl VelocityExpansion {
    const ZERO: Self = Self { hi: 0.0, lo: 0.0 };

    fn sum(a: f64, b: f64) -> Self {
        let hi = a + b;
        let b_part = hi - a;
        let lo = (a - (hi - b_part)) + (b - b_part);
        Self { hi, lo }
    }

    fn product(a: f64, b: f64) -> Self {
        let hi = a * b;
        let lo = a.mul_add(b, -hi);
        Self { hi, lo }
    }

    fn add(self, other: Self) -> Self {
        let high = Self::sum(self.hi, other.hi);
        let low = Self::sum(self.lo, other.lo);
        let middle = Self::sum(high.lo, low.hi);
        let combined = Self::sum(high.hi, middle.hi);
        Self::sum(combined.hi, (combined.lo + middle.lo) + low.lo)
    }

    fn sub(self, other: Self) -> Self {
        self.add(Self {
            hi: -other.hi,
            lo: -other.lo,
        })
    }

    fn scaled(self, factor: f64) -> Self {
        Self::product(self.hi, factor).add(Self::product(self.lo, factor))
    }

    fn multiply(self, other: Self) -> Self {
        Self::product(self.hi, other.hi)
            .add(Self::product(self.hi, other.lo))
            .add(Self::product(self.lo, other.hi))
            .add(Self::product(self.lo, other.lo))
    }

    fn rounded(self) -> f64 {
        self.hi + self.lo
    }
}

/// Paired Clenshaw on the derivative series, retaining the common origin.
/// Coefficient scaling and recurrence registers preserve their residuals until
/// the returned velocity. Scalar positions, normalized arguments and selected
/// records remain independent of this derivative-only arithmetic.
fn polynomial_velocity_difference(
    target: &[f64],
    target_tau: f64,
    target_scale: f64,
    center: &[f64],
    center_tau: f64,
    center_scale: f64,
) -> Result<f64> {
    let count = target.len().max(center.len());
    if count < 2 {
        return Ok(0.0);
    }
    let coefficient = |values: &[f64], index| values.get(index).copied().unwrap_or(0.0);
    let difference_coefficient = |index| {
        VelocityExpansion::product(coefficient(target, index), target_scale).sub(
            VelocityExpansion::product(coefficient(center, index), center_scale),
        )
    };
    // The arguments themselves are unchanged. TwoSum retains the subtraction
    // residual if these already selected binary64 arguments differ in scale.
    let tau_difference = VelocityExpansion::sum(target_tau, -center_tau);
    let (mut center_next, mut center_next2) = (VelocityExpansion::ZERO, VelocityExpansion::ZERO);
    let (mut difference_next, mut difference_next2) =
        (VelocityExpansion::ZERO, VelocityExpansion::ZERO);
    let (mut c1, mut c2, mut d1, mut d2) = (
        VelocityExpansion::ZERO,
        VelocityExpansion::ZERO,
        VelocityExpansion::ZERO,
        VelocityExpansion::ZERO,
    );
    for index in (1..count - 1).rev() {
        let factor = 2.0 * (index + 1) as f64;
        let center_derivative =
            VelocityExpansion::product(coefficient(center, index + 1), center_scale)
                .scaled(factor)
                .add(center_next2);
        let difference_derivative = difference_coefficient(index + 1)
            .scaled(factor)
            .add(difference_next2);
        (center_next2, center_next) = (center_next, center_derivative);
        (difference_next2, difference_next) = (difference_next, difference_derivative);
        let difference = d1
            .scaled(2.0 * target_tau)
            .add(tau_difference.scaled(2.0).multiply(c1))
            .add(difference_derivative)
            .sub(d2);
        let value = c1.scaled(2.0 * center_tau).add(center_derivative).sub(c2);
        (c2, c1, d2, d1) = (c1, value, d1, difference);
    }
    let constant = difference_next2.scaled(0.5).add(difference_coefficient(1));
    let difference = d1
        .scaled(target_tau)
        .add(tau_difference.multiply(c1))
        .add(constant)
        .sub(d2)
        .rounded();
    if difference.is_finite() {
        Ok(difference)
    } else {
        Err(Error::NonFiniteArgument { value: difference })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn state_bits(value: CartesianState) -> ([u64; 3], [u64; 3]) {
        (value.0.map(f64::to_bits), value.1.map(f64::to_bits))
    }

    fn word_bits(value: crate::positions::CartesianWordParts) -> [u64; 12] {
        std::array::from_fn(|i| {
            if i < 6 {
                value.0[i / 2][i % 2].to_bits()
            } else {
                value.1[(i - 6) / 2][i % 2].to_bits()
            }
        })
    }

    #[test]
    fn borrowed_rows_keep_both_recurrences_and_survive_cache_release() {
        for format in [LebFormat::Leb2V1, LebFormat::Leb2V2] {
            for degree in [2, 24] {
                let (index, bytes, _) = scoped_relative_fixture_degree(format, degree);
                let source = LebSource::Slice(&bytes);
                let reference_cache = LebChunkCache::new();
                for hi in [0.0, 4.0_f64.next_down(), 4.0, 4.0_f64.next_up(), 15.5, 16.0] {
                    for lo in [-0.0, 0.0, 1e-18, -1e-18] {
                        let epoch = TwoPartEpoch { hi, lo };
                        // Stay inside the physical stored domain, including its
                        // endpoints; adjacent representable joins remain covered.
                        if epoch.value() < 0.0 || epoch.value() > 16.0 {
                            continue;
                        }
                        for body in [1, 14] {
                            let expected_word = index
                                .eval_raw_word_parts_two_part_source(
                                    source,
                                    body,
                                    epoch,
                                    &reference_cache,
                                )
                                .unwrap()
                                .unwrap();
                            let expected_scalar = index
                                .eval_state_inline_two_part_source(
                                    source,
                                    body,
                                    epoch,
                                    &reference_cache,
                                )
                                .unwrap()
                                .unwrap();
                            let cache = LebChunkCache::with_image_byte_budget(usize::MAX);
                            let actual = index
                                .with_native_word_parts_two_part_source(
                                    source,
                                    body,
                                    epoch,
                                    &cache,
                                    |raw, selected, row| {
                                        assert!(std::ptr::eq(selected, &index));
                                        let coefficients = row.expect("native contiguous LEB2 row");
                                        // This reacquires the coefficient lock. The
                                        // callback keeps its row owner, not that guard.
                                        cache.clear();
                                        assert_eq!(cache.memory_usage().0, 0);
                                        let scalar = selected
                                            .eval_borrowed_state_inline_two_part(
                                                body,
                                                epoch,
                                                coefficients,
                                            )?
                                            .unwrap();
                                        Ok((raw.unwrap(), scalar))
                                    },
                                )
                                .unwrap();
                            assert_eq!(word_bits(actual.0), word_bits(expected_word));
                            assert_eq!(state_bits(actual.1), state_bits(expected_scalar));
                            assert_eq!(cache.stats().1, 1);
                            assert_eq!(cache.memory_usage().0, 0);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn borrowed_rows_exclude_positional_public_and_unsupported_sources() {
        use std::cell::Cell;
        for format in [LebFormat::Leb1, LebFormat::Leb2V1, LebFormat::Leb2V2] {
            let (index, bytes, _) = scoped_relative_fixture(format);
            let epoch = TwoPartEpoch { hi: 4.0, lo: 1e-18 };
            let temporary = ScopedTemporary::new(&bytes);
            let file = access::PositionalFile::open(
                &temporary.0,
                resources::IndexedResources::new(IndexedLimits {
                    cache_bytes: 0,
                    ..Default::default()
                }),
            )
            .unwrap();
            let public = LebChunkCache::with_byte_budget(0);
            for (source, cache) in [
                (LebSource::Slice(&bytes), &public),
                (LebSource::Positional(&file), &file.cache),
            ] {
                let expected = index
                    .eval_raw_word_parts_two_part_source(source, 1, epoch, cache)
                    .unwrap();
                let before = cache.stats().1;
                let actual = index
                    .with_native_word_parts_two_part_source(
                        source,
                        1,
                        epoch,
                        cache,
                        |raw, selected, row| {
                            assert!(row.is_none());
                            // Positional Word ownership has ended before the scalar
                            // acquisition: no extra live working reservation survives.
                            assert_eq!(file.resources.working.used(), 0);
                            let _ = selected
                                .eval_state_inline_two_part_source(source, 1, epoch, cache)?;
                            Ok(raw)
                        },
                    )
                    .unwrap();
                assert_eq!(actual.map(word_bits), expected.map(word_bits));
                if format != LebFormat::Leb1 {
                    assert_eq!(cache.stats().1 - before, 2);
                }
                assert_eq!(file.resources.working.used(), 0);
            }
            let cache = LebChunkCache::with_image_byte_budget(0);
            let called = Cell::new(false);
            for invalid in [
                f64::from_bits(0xfff8_0000_0000_1234),
                f64::INFINITY,
                -1.0,
                17.0,
            ] {
                let epoch = TwoPartEpoch {
                    hi: invalid,
                    lo: 0.0,
                };
                let expected = index
                    .eval_raw_word_parts_two_part_source(LebSource::Slice(&bytes), 1, epoch, &cache)
                    .unwrap_err();
                let actual = index
                    .with_native_word_parts_two_part_source(
                        LebSource::Slice(&bytes),
                        1,
                        epoch,
                        &cache,
                        |_, _, _| {
                            called.set(true);
                            Ok(())
                        },
                    )
                    .unwrap_err();
                match (actual, expected) {
                    (
                        Error::NonFiniteArgument { value: actual },
                        Error::NonFiniteArgument { value: expected },
                    ) => {
                        if expected.is_nan() {
                            assert!(actual.is_nan());
                        } else {
                            assert_eq!(actual.to_bits(), expected.to_bits());
                        }
                    }
                    (actual, expected) => {
                        assert_eq!(format!("{actual:?}"), format!("{expected:?}"))
                    }
                }
            }
            let truncated = LebSource::Slice(&bytes[..bytes.len() - 1]);
            let end = TwoPartEpoch { hi: 16.0, lo: 0.0 };
            let expected = index
                .eval_raw_word_parts_two_part_source(truncated, 1, end, &cache)
                .unwrap_err();
            let actual = index
                .with_native_word_parts_two_part_source(truncated, 1, end, &cache, |_, _, _| {
                    called.set(true);
                    Ok(())
                })
                .unwrap_err();
            assert_eq!(format!("{actual:?}"), format!("{expected:?}"));
            assert!(!called.get());
            let mut unsupported = index.clone();
            unsupported
                .bodies
                .iter_mut()
                .for_each(|entry| entry.coord_type = COORD_ECLIPTIC);
            unsupported
                .with_native_word_parts_two_part_source(
                    LebSource::Slice(&bytes),
                    1,
                    epoch,
                    &cache,
                    |raw, _, row| {
                        assert!(raw.is_none() && row.is_none());
                        Ok(())
                    },
                )
                .unwrap();
        }
    }

    #[test]
    fn native_word_producers_share_one_row_from_the_actual_image_and_tier() {
        use crate::positions::{SingleLeb, StateProvider};
        // The existing synthetic fixture already supplies a decoded index;
        // construct its immutable image owner without adding a second serializer.
        let image = |index, bytes| LebFile {
            path: "synthetic-borrowed-row".into(),
            image: Some(Arc::new(LebImage {
                bytes: storage::ImageStorage::Owned(bytes),
                index,
                cache: Arc::new(LebChunkCache::with_image_byte_budget(0)),
                sha256: std::sync::OnceLock::new(),
            })),
        };
        for format in [LebFormat::Leb2V1, LebFormat::Leb2V2] {
            let (mut short_index, short_bytes, _) = scoped_relative_fixture_degree(format, 2);
            // The base image has body 14 but not body 1; medium supplies body 1.
            // Shared body 14 also distinguishes first-wins from tier priority.
            short_index
                .bodies
                .iter_mut()
                .find(|entry| entry.body_id == 1)
                .unwrap()
                .body_id = 999;
            let short = image(short_index, short_bytes);
            let (long_index, long_bytes, _) = scoped_relative_fixture_degree(format, 24);
            let long = image(long_index, long_bytes);
            let short_image = short.image.as_ref().unwrap();
            let long_image = long.image.as_ref().unwrap();
            let single = SingleLeb::from_file(&long).unwrap();
            let composite = CompositeIndex::from_files([&long, &short]).unwrap();
            let tiered = TieredIndex::new(vec![
                (
                    LebTier::Medium,
                    CompositeIndex::from_files([&long]).unwrap(),
                ),
                (LebTier::Base, CompositeIndex::from_files([&short]).unwrap()),
            ])
            .unwrap();
            for epoch in [
                TwoPartEpoch { hi: 4.0, lo: 1e-18 },
                TwoPartEpoch {
                    hi: 12.0,
                    lo: -1e-18,
                },
            ] {
                for body in [1, 14] {
                    for route in 0..3 {
                        let selected = if route == 2 && body == 14 {
                            short_image
                        } else {
                            long_image
                        };
                        let reference =
                            SingleLeb::new(&selected.index, selected.bytes.bytes().unwrap());
                        let expected = reference
                            .geometric_state_word_parts_two_part(body, epoch)
                            .unwrap();
                        let before = selected.cache.stats().1;
                        let owners = Arc::strong_count(selected);
                        let actual = match route {
                            0 => single.geometric_state_word_parts_two_part(body, epoch),
                            1 => composite.geometric_state_word_parts_two_part(body, epoch),
                            _ => tiered.geometric_state_word_parts_two_part(body, epoch),
                        }
                        .unwrap();
                        assert_eq!(word_bits(actual), word_bits(expected));
                        assert_eq!(selected.cache.stats().1 - before, 1);
                        assert_eq!(selected.cache.memory_usage().0, 0);
                        assert_eq!(Arc::strong_count(selected), owners);
                    }
                }
            }
        }
    }

    #[test]
    fn compact_word_memo_keeps_actual_row_reads_before_reuse_at_zero_chunk_budget() {
        let (index, bytes) = polynomial_pair();
        let epoch = TwoPartEpoch { hi: 2.0, lo: 0.0 };
        let uncached = LebChunkCache::new();
        let cache = LebChunkCache::with_image_byte_budget(0);
        let word_bits = |value: Option<crate::positions::CartesianWordParts>| {
            let (position, velocity) = value.unwrap();
            (
                position.map(|pair| pair.map(f64::to_bits)),
                velocity.map(|pair| pair.map(f64::to_bits)),
            )
        };
        let expected = index
            .eval_raw_word_parts_two_part_source(LebSource::Slice(&bytes), 1, epoch, &uncached)
            .unwrap();
        for _ in 0..2 {
            let actual = index
                .eval_raw_word_parts_two_part_source(LebSource::Slice(&bytes), 1, epoch, &cache)
                .unwrap();
            assert_eq!(word_bits(actual), word_bits(expected));
        }
        // Even after a successful result was memoized, the actual row must be
        // acquired first. This controlled pre-query fault must still be seen.
        let entry = index.find_body(1).unwrap();
        let short_length = entry.data_offset as usize
            + 3 * (entry.degree as usize + 1) * std::mem::size_of::<f64>()
            - 1;
        assert!(matches!(
            index.eval_raw_word_parts_two_part_source(
                LebSource::Slice(&bytes[..short_length]),
                1,
                epoch,
                &cache,
            ),
            Err(Error::LebCorrupt { .. })
        ));
        assert_eq!(cache.memory_usage().0, 0);
    }

    #[test]
    fn compact_word_memo_keeps_nonfinite_error_policy_and_never_retains_errors() {
        let (index, mut bytes) = polynomial_pair();
        let offset = index.find_body(1).unwrap().data_offset as usize;
        bytes[offset..offset + 8].copy_from_slice(&f64::NAN.to_le_bytes());
        let epoch = TwoPartEpoch { hi: 2.0, lo: 0.0 };
        let uncached = LebChunkCache::new();
        let cache = LebChunkCache::with_image_byte_budget(0);
        let expected = index
            .eval_raw_word_parts_two_part_source(LebSource::Slice(&bytes), 1, epoch, &uncached)
            .unwrap_err();
        let key = crate::pure_cache::RawWordStateKey::new(
            cache.word_image_identity().unwrap(),
            1,
            0,
            [epoch.hi.to_bits(), epoch.lo.to_bits()],
        );
        for _ in 0..2 {
            let actual = index
                .eval_raw_word_parts_two_part_source(LebSource::Slice(&bytes), 1, epoch, &cache)
                .unwrap_err();
            assert_eq!(actual.to_string(), expected.to_string());
            assert!(matches!(actual, Error::MathDomain { value } if !value.is_finite()));
            assert!(crate::pure_cache::raw_word_state(key).is_none());
        }
    }

    #[test]
    fn paired_velocity_retains_scaled_coefficient_residual_through_cancellation() {
        let epsilon = f64::EPSILON;
        let scale = 1.0 - epsilon;
        let target = [0.0, -4.0, 1.0 + epsilon];
        let center = [0.0, 0.0, 1.0 + epsilon];
        // Direct differentiation at target tau=1, center tau=0 gives
        // scale * (-4 + 4 * (1 + epsilon)) = 4*epsilon - 4*epsilon^2.
        // All terms and the expected result are exactly representable. Rounding
        // the scaled quadratic coefficient first loses the epsilon^2 term.
        let expected = 4.0 * epsilon - 4.0 * epsilon * epsilon;
        let actual =
            polynomial_velocity_difference(&target, 1.0, scale, &center, 0.0, scale).unwrap();
        assert_eq!(actual.to_bits(), expected.to_bits());
    }

    #[test]
    fn paired_velocity_preserves_common_origin_and_separate_epoch_derivatives() {
        // Values follow direct differentiation of T0..T4 at exact binary64
        // arguments, evaluated independently with 100-decimal arithmetic.
        type VelocityCase<'a> = (&'a [f64], &'a [f64], f64, f64, f64, f64, f64);
        let cases: &[VelocityCase<'_>] = &[
            (
                &[1e16, 1e12, 1.0, 0.125, 0.03125],
                &[1e16, 1e12, 0.5, -0.25, -0.0625],
                0.25,
                0.25,
                0.5,
                0.5,
                -0.3359375,
            ),
            (
                &[1e16, 1e12, 1e8, 0.125, 0.03125],
                &[1e16, 1e12, 1e8, -0.25, -0.0625],
                0.25,
                0.25 + 2.0_f64.powi(-40),
                0.5,
                0.5,
                -0.5861193989399567,
            ),
            (
                &[1e16, 2e12, 2e8, 0.125, 0.03125],
                &[1e16, 1e12, 1e8, -0.25],
                -0.5,
                0.25,
                0.25,
                0.5,
                -150000000.25,
            ),
            (&[1e16, 1e12], &[1e16, 1e12], -0.5, 0.25, 0.5, 0.5, 0.0),
        ];
        for &(target, center, tx, cx, ts, cs, expected) in cases {
            let actual = polynomial_velocity_difference(target, tx, ts, center, cx, cs).unwrap();
            assert!(
                (actual - expected).abs() <= 2.0 * f64::EPSILON * expected.abs(),
                "{actual} != {expected}"
            );
        }
    }

    fn polynomial_pair() -> (FileIndex, Vec<u8>) {
        let body = |id, offset| BodyEntry {
            body_id: id,
            coord_type: COORD_ICRS_BARY,
            segment_count: 1,
            jd_start: 0.0,
            jd_end: 4.0,
            interval_days: 4.0,
            degree: 2,
            components: 3,
            data_offset: offset,
            compressed_size: None,
            uncompressed_size: None,
        };
        let header = FileHeader {
            magic: MAGIC_LEB1,
            version: VERSION_LEB1,
            section_count: 0,
            body_count: 2,
            jd_start: 0.0,
            jd_end: 4.0,
            generation_epoch: 0.0,
            flags: 0,
            format: LebFormat::Leb1,
        };
        let rows: [[f64; 9]; 2] = [
            [1.0, 0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [
                1.0 + 2.0_f64.powi(-20),
                0.25 + 2.0_f64.powi(-30),
                2.0_f64.powi(-50),
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ],
        ];
        let bytes = rows
            .into_iter()
            .flatten()
            .flat_map(f64::to_le_bytes)
            .collect();
        (
            FileIndex {
                header,
                sections: vec![],
                bodies: vec![body(14, 0), body(1, 72)],
            },
            bytes,
        )
    }

    #[test]
    fn aligned_difference_keeps_small_position_and_its_analytic_derivative() {
        let (index, bytes) = polynomial_pair();
        // tau=1/8. The residual polynomial and this evaluation are exact binary
        // fractions; its final term is lost in the independently rounded state.
        let epoch = 2.25;
        let expected_position = 2.0_f64.powi(-20) + 2.0_f64.powi(-33) - 31.0 * 2.0_f64.powi(-55);
        let expected_velocity = (2.0_f64.powi(-30) + 2.0_f64.powi(-51)) / 2.0;
        let actual = index
            .eval_relative_icrs_cached(&bytes, 1, 14, epoch, &LebChunkCache::new())
            .unwrap()
            .unwrap();
        assert_eq!(actual.0, [expected_position, 0.0, 0.0]);
        assert_eq!(actual.1, [expected_velocity, 0.0, 0.0]);
        let target = index.eval_body(&bytes, 1, epoch).unwrap();
        let center = index.eval_body(&bytes, 14, epoch).unwrap();
        assert!((target.0[0] - center.0[0] - expected_position).abs() > 2.0_f64.powi(-60));
        let reverse = index
            .eval_relative_icrs_cached(&bytes, 14, 1, epoch, &LebChunkCache::new())
            .unwrap()
            .unwrap();
        assert_eq!(reverse.0[0], -expected_position);
        assert_eq!(reverse.1[0], -expected_velocity);
    }

    #[test]
    fn two_epoch_difference_retains_sub_julian_day_ulp_displacements() {
        let (mut index, bytes) = polynomial_pair();
        for entry in &mut index.bodies {
            entry.jd_start += 2451544.0;
            entry.jd_end += 2451544.0;
        }
        let center = TwoPartEpoch::split(2451546.25).unwrap();
        let target = center.add_offset(-2.0_f64.powi(-40)).unwrap();
        assert_eq!(center.value(), target.value());
        let tau = 0.125 - 2.0_f64.powi(-41);
        let expected = 2.0_f64.powi(-20) - 2.0_f64.powi(-43)
            + 2.0_f64.powi(-30) * tau
            + 2.0_f64.powi(-50) * (2.0 * tau * tau - 1.0);
        let result = index
            .eval_relative_position_two_part_cached(
                &bytes,
                1,
                target,
                14,
                center,
                &LebChunkCache::new(),
            )
            .unwrap()
            .unwrap();
        assert!((result[0] - expected).abs() <= f64::EPSILON * expected.abs());
        let reverse = index
            .eval_relative_position_two_part_cached(
                &bytes,
                14,
                center,
                1,
                target,
                &LebChunkCache::new(),
            )
            .unwrap()
            .unwrap();
        assert!((reverse[0] + expected).abs() <= f64::EPSILON * expected.abs());
        assert_eq!(&result[1..], &[0.0, 0.0]);
    }

    #[test]
    fn two_epoch_polynomials_keep_independent_records_across_a_join() {
        let (mut index, _) = polynomial_pair();
        let (v, u, separation) = (2.0_f64.powi(-10), 2.0_f64.powi(-20), 2.0_f64.powi(-25));
        let mut bytes = Vec::new();
        for (body, entry) in index.bodies.iter_mut().enumerate() {
            entry.degree = 1;
            entry.segment_count = 2;
            entry.interval_days = 2.0;
            entry.data_offset = bytes.len() as u64;
            let rate = if body == 0 { v } else { v + u };
            let origin = if body == 0 { 1.0 } else { 1.0 + separation };
            for middle in [1.0, 3.0] {
                for coefficient in [origin + rate * (middle - 2.0), rate, 0.0, 0.0, 0.0, 0.0] {
                    bytes.extend_from_slice(&coefficient.to_le_bytes());
                }
            }
        }
        let center = TwoPartEpoch::split(2.0).unwrap();
        let target = center.add_offset(-2.0_f64.powi(-30)).unwrap();
        assert_eq!(select_segment_two_part(0.0, 2.0, 2, target).unwrap().0, 0);
        assert_eq!(select_segment_two_part(0.0, 2.0, 2, center).unwrap().0, 1);
        let expected = separation - (v + u) * 2.0_f64.powi(-30);
        let result = index
            .eval_relative_position_two_part_cached(
                &bytes,
                1,
                target,
                14,
                center,
                &LebChunkCache::new(),
            )
            .unwrap()
            .unwrap();
        assert_eq!(result, [expected, 0.0, 0.0]);
        assert!(matches!(
            index.eval_relative_position_two_part_cached(
                &bytes,
                1,
                TwoPartEpoch::split(-1.0).unwrap(),
                14,
                center,
                &LebChunkCache::new()
            ),
            Err(Error::LebOutOfRange { body_id: 1, .. })
        ));
    }

    #[test]
    fn paired_clenshaw_differentiates_arguments_without_losing_the_common_origin() {
        // T_2 and T_3 identities give an independent exact binary polynomial.
        let a = [1.0, 0.25, 0.125, 0.0625];
        let b = [1.0, 0.25];
        let (x, y) = (0.125_f64, -0.375_f64);
        let expected =
            0.25 * (x - y) + 0.125 * (2.0 * x * x - 1.0) + 0.0625 * (4.0 * x * x * x - 3.0 * x);
        assert_eq!(polynomial_difference(&a, x, &b, y).unwrap(), expected);
        let tiny = 2.0_f64.powi(-45);
        assert_eq!(
            polynomial_difference(&[1.0, 0.25], tiny, &[1.0, 0.25], 0.0).unwrap(),
            0.25 * tiny
        );
        assert!(polynomial_difference(&[f64::NAN], x, &b, y).is_err());
    }

    #[test]
    fn incompatible_frames_grids_and_uncovered_epochs_are_not_combined() {
        let (index, bytes) = polynomial_pair();
        let cache = LebChunkCache::new();
        let mut changed = index.clone();
        changed.bodies[1].coord_type = COORD_ICRS_BARY_SYSTEM;
        assert!(
            changed
                .eval_relative_icrs_cached(&bytes, 1, 14, 2.25, &cache)
                .unwrap()
                .is_none()
        );
        changed = index.clone();
        changed.bodies[1].interval_days = 2.0;
        assert!(
            changed
                .eval_relative_icrs_cached(&bytes, 1, 14, 2.25, &cache)
                .unwrap()
                .is_none()
        );
        assert!(
            index
                .eval_relative_icrs_cached(&bytes, 1, 14, 5.0, &cache)
                .unwrap()
                .is_none()
        );
        assert!(matches!(
            index.eval_relative_icrs_cached(&bytes, 999, 14, 2.25, &cache),
            Err(Error::LebBodyNotFound { body_id: 999 })
        ));
        assert!(matches!(
            index.eval_relative_icrs_cached(&bytes, 1, 14, f64::NAN, &cache),
            Err(Error::NonFiniteArgument { .. })
        ));
    }
    /// Distinct four-record rows retain actual joins in all three formats.
    fn scoped_relative_fixture(format: LebFormat) -> (FileIndex, Vec<u8>, Vec<u64>) {
        scoped_relative_fixture_degree(format, 2)
    }

    fn scoped_relative_fixture_degree(
        format: LebFormat,
        degree: u32,
    ) -> (FileIndex, Vec<u8>, Vec<u64>) {
        let (mut index, original) = polynomial_pair();
        let mut raw_rows = Vec::new();
        for (body, entry) in index.bodies.iter_mut().enumerate() {
            let original_row = &original[body * 72..(body + 1) * 72];
            let mut row = Vec::new();
            for axis in 0..3 {
                row.extend_from_slice(&original_row[axis * 24..(axis + 1) * 24]);
                row.resize((axis + 1) * (degree as usize + 1) * 8, 0);
            }
            if degree > 2 {
                row[degree as usize * 8..(degree as usize + 1) * 8]
                    .copy_from_slice(&(body as f64 * 2.0_f64.powi(-30)).to_le_bytes());
            }
            entry.degree = degree;
            let mut raw = Vec::new();
            for record in 0..4 {
                let mut row = row.clone();
                let coefficient = f64::from_le_bytes(row[..8].try_into().unwrap());
                row[..8].copy_from_slice(&(coefficient + record as f64 * 0.125).to_le_bytes());
                raw.extend(row);
            }
            entry.segment_count = 4;
            entry.jd_end = 16.0;
            raw_rows.push(raw);
        }
        index.header.format = format;
        index.header.jd_end = 16.0;
        let mut bytes = Vec::new();
        let mut payloads = Vec::new();
        for (entry, raw) in index.bodies.iter_mut().zip(raw_rows) {
            entry.data_offset = bytes.len() as u64;
            if format == LebFormat::Leb1 {
                payloads.push(entry.data_offset);
                bytes.extend(raw);
                continue;
            }
            let shuffled =
                shuffle_bytes(&reorder_to_coeff_major(&raw, 4, degree, 3).unwrap()).unwrap();
            let compressed = zstd::encode_all(std::io::Cursor::new(shuffled), 3).unwrap();
            entry.compressed_size = Some(compressed.len() as u64);
            entry.uncompressed_size = Some(raw.len() as u64);
            if format == LebFormat::Leb2V2 {
                let blob = entry.data_offset + 16 + 48;
                bytes.extend(1_u32.to_le_bytes());
                bytes.extend(0_u32.to_le_bytes());
                bytes.extend(16.0_f32.to_le_bytes());
                bytes.extend(0_u32.to_le_bytes());
                bytes.extend(0.0_f64.to_le_bytes());
                bytes.extend(16.0_f64.to_le_bytes());
                bytes.extend(blob.to_le_bytes());
                bytes.extend((compressed.len() as u64).to_le_bytes());
                bytes.extend((raw.len() as u64).to_le_bytes());
                bytes.extend(0_u32.to_le_bytes());
                bytes.extend(4_u32.to_le_bytes());
                payloads.push(blob);
            } else {
                payloads.push(entry.data_offset);
            }
            bytes.extend(compressed);
        }
        (index, bytes, payloads)
    }

    struct ScopedTemporary(std::path::PathBuf);
    impl ScopedTemporary {
        fn new(bytes: &[u8]) -> Self {
            use std::sync::atomic::{AtomicU64, Ordering};
            static NEXT: AtomicU64 = AtomicU64::new(0);
            let path = std::env::temp_dir().join(format!(
                "relative-scoped-{}-{}",
                std::process::id(),
                NEXT.fetch_add(1, Ordering::Relaxed)
            ));
            std::fs::write(&path, bytes).unwrap();
            Self(path)
        }
    }
    impl Drop for ScopedTemporary {
        fn drop(&mut self) {
            let _ = std::fs::remove_file(&self.0);
        }
    }

    #[test]
    fn scoped_relative_rows_preserve_formats_epochs_and_evicted_owners() {
        for degree in [2, 24] {
            let (raw_index, raw, _) = scoped_relative_fixture_degree(LebFormat::Leb1, degree);
            for format in [LebFormat::Leb1, LebFormat::Leb2V1, LebFormat::Leb2V2] {
                let (index, bytes, _) = scoped_relative_fixture_degree(format, degree);
                let temporary = ScopedTemporary::new(&bytes);
                for retention in [0, 288, 576] {
                    let limits = IndexedLimits {
                        cache_bytes: retention,
                        window_log_max: 23,
                        ..Default::default()
                    };
                    let file = access::PositionalFile::open(
                        &temporary.0,
                        resources::IndexedResources::new(limits),
                    )
                    .unwrap();
                    let owned_cache = LebChunkCache::with_byte_budget(retention);
                    for source in [LebSource::Slice(&bytes), LebSource::Positional(&file)] {
                        let cache = if matches!(source, LebSource::Slice(_)) {
                            &owned_cache
                        } else {
                            &file.cache
                        };
                        for jd in [
                            0.0,
                            0.25,
                            4.0_f64.next_down(),
                            4.0,
                            4.0_f64.next_up(),
                            7.5,
                            8.0,
                            15.5,
                            16.0,
                        ] {
                            let tail = if jd == 0.0 || jd == 16.0 { 0.0 } else { 1e-18 };
                            let target = TwoPartEpoch { hi: jd, lo: tail };
                            let center = TwoPartEpoch { hi: jd, lo: -tail };
                            let expected = raw_index
                                .eval_relative_position_two_part_cached_source(
                                    LebSource::Slice(&raw),
                                    1,
                                    target,
                                    14,
                                    center,
                                    &LebChunkCache::new(),
                                )
                                .unwrap();
                            assert_eq!(
                                index
                                    .eval_relative_position_two_part_cached_source(
                                        source, 1, target, 14, center, cache
                                    )
                                    .unwrap(),
                                expected
                            );
                            let expected_velocity = raw_index
                                .eval_relative_velocity_two_part_cached_source(
                                    LebSource::Slice(&raw),
                                    1,
                                    target,
                                    14,
                                    center,
                                    &LebChunkCache::new(),
                                )
                                .unwrap();
                            let actual_velocity = index
                                .eval_relative_velocity_two_part_cached_source(
                                    source, 1, target, 14, center, cache,
                                )
                                .unwrap();
                            assert_eq!(
                                actual_velocity.map(|row| row.map(f64::to_bits)),
                                expected_velocity.map(|row| row.map(f64::to_bits))
                            );
                            let expected_words = raw_index
                                .eval_raw_word_parts_two_part_source(
                                    LebSource::Slice(&raw),
                                    1,
                                    target,
                                    &LebChunkCache::new(),
                                )
                                .unwrap();
                            let actual_words = index
                                .eval_raw_word_parts_two_part_source(source, 1, target, cache)
                                .unwrap();
                            let bits = |state: crate::positions::CartesianWordParts| {
                                (
                                    state.0.map(|row| row.map(f64::to_bits)),
                                    state.1.map(|row| row.map(f64::to_bits)),
                                )
                            };
                            assert_eq!(actual_words.map(bits), expected_words.map(bits));
                            let expected = raw_index
                                .eval_relative_icrs_cached_source(
                                    LebSource::Slice(&raw),
                                    1,
                                    14,
                                    jd,
                                    &LebChunkCache::new(),
                                )
                                .unwrap();
                            assert_eq!(
                                index
                                    .eval_relative_icrs_cached_source(source, 1, 14, jd, cache)
                                    .unwrap(),
                                expected
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn scoped_relative_rows_keep_payload_error_before_later_epoch_validation() {
        for format in [LebFormat::Leb2V1, LebFormat::Leb2V2] {
            let (index, mut bytes, payloads) = scoped_relative_fixture(format);
            let epoch = TwoPartEpoch { hi: 1.0, lo: 0.0 };
            let outside = TwoPartEpoch { hi: 100.0, lo: 0.0 };
            let target_offset = payloads[1] as usize;
            bytes[target_offset..target_offset + 4].fill(0);
            let cache = LebChunkCache::new();
            let expected = index
                .eval_raw_acceleration_two_part_source(LebSource::Slice(&bytes), 1, epoch, &cache)
                .unwrap_err();
            let actual = index
                .eval_relative_position_two_part_cached_source(
                    LebSource::Slice(&bytes),
                    1,
                    epoch,
                    14,
                    outside,
                    &cache,
                )
                .unwrap_err();
            assert_eq!(format!("{actual:?}"), format!("{expected:?}"));
            let actual_velocity = index
                .eval_relative_velocity_two_part_cached_source(
                    LebSource::Slice(&bytes),
                    1,
                    epoch,
                    14,
                    outside,
                    &cache,
                )
                .unwrap_err();
            assert_eq!(format!("{actual_velocity:?}"), format!("{expected:?}"));
            assert!(matches!(
                index.eval_relative_velocity_two_part_cached_source(
                    LebSource::Slice(&bytes),
                    1,
                    outside,
                    14,
                    epoch,
                    &cache
                ),
                Err(Error::LebOutOfRange { body_id: 1, .. })
            ));
            assert!(matches!(
                index.eval_relative_position_two_part_cached_source(
                    LebSource::Slice(&bytes),
                    1,
                    outside,
                    14,
                    epoch,
                    &cache
                ),
                Err(Error::LebOutOfRange { body_id: 1, .. })
            ));
        }
    }

    #[test]
    fn bounded_target_rows_preserve_tight_working_budget_without_retention() {
        for degree in [2, 24] {
            for format in [LebFormat::Leb2V1, LebFormat::Leb2V2] {
                let (index, bytes, _) = scoped_relative_fixture_degree(format, degree);
                let row_bytes = (degree as usize + 1) * 3 * 8;
                let frame_bytes = 4 * row_bytes;
                let largest_compressed = index
                    .bodies
                    .iter()
                    .map(|body| body.compressed_size.unwrap() as usize)
                    .max()
                    .unwrap();
                let maximum_scratch = 4 * frame_bytes + largest_compressed + (1 << 23) + 1;
                // Only the selected target row may survive the next decode. The
                // bounded stack copy releases its chunk and keeps this same tight
                // heap budget sufficient, including disabled coefficient retention.
                let limits = IndexedLimits {
                    cache_bytes: 0,
                    working_bytes: maximum_scratch + row_bytes,
                    window_log_max: 23,
                    ..Default::default()
                };
                let temporary = ScopedTemporary::new(&bytes);
                let file = access::PositionalFile::open(
                    &temporary.0,
                    resources::IndexedResources::new(limits),
                )
                .unwrap();
                let source = LebSource::Positional(&file);
                let copy_row = |body, entry| {
                    index.with_relative_segment_coefficients(
                        source,
                        body,
                        entry,
                        0,
                        1.0,
                        &file.cache,
                        |words| {
                            let allocation = source.reserve_work(words.len() * 8)?;
                            Ok((words.to_vec(), allocation))
                        },
                    )
                };
                {
                    let _target = copy_row(1, index.find_body(1).unwrap()).unwrap();
                    let _center = copy_row(14, index.find_body(14).unwrap()).unwrap();
                }
                assert_eq!(file.resources.working.used(), 0);
                let epoch = TwoPartEpoch { hi: 1.0, lo: 0.0 };
                assert!(
                    index
                        .eval_relative_position_two_part_cached_source(
                            source,
                            1,
                            epoch,
                            14,
                            epoch,
                            &file.cache,
                        )
                        .unwrap()
                        .is_some()
                );
                assert_eq!(file.resources.working.used(), 0);
                assert!(
                    index
                        .eval_relative_velocity_two_part_cached_source(
                            source,
                            1,
                            epoch,
                            14,
                            epoch,
                            &file.cache,
                        )
                        .unwrap()
                        .is_some()
                );
                assert_eq!(file.resources.working.used(), 0);
                assert!(
                    index
                        .eval_raw_word_parts_two_part_source(source, 1, epoch, &file.cache)
                        .unwrap()
                        .is_some()
                );
                assert_eq!(file.resources.working.used(), 0);
            }
        }
    }
}

impl FileIndex {
    /// Residual-preserving state of the same selected Cartesian polynomial.
    /// Pair components are `[high, residual]`; units are ICRS AU and AU/TT-day.
    /// Record ownership uses the existing selector. Its central midpoint is
    /// unchanged; the two-part epoch is not collapsed during normalization.
    pub(crate) fn eval_raw_word_parts_two_part_source(
        &self,
        source: LebSource<'_>,
        body: i32,
        epoch: TwoPartEpoch,
        cache: &LebChunkCache,
    ) -> Result<Option<crate::positions::CartesianWordParts>> {
        self.with_raw_word_parts_two_part_source(source, body, epoch, cache, |raw, _| Ok(raw))
    }

    /// Lend the already acquired row only for an immutable contiguous LEB2
    /// producer. Positional reads, arbitrary index/source pairs and LEB1 keep
    /// their original completed acquisition before native state construction.
    /// No cache lock is held: `consume` borrows the first decoded owner without
    /// another Arc clone. Its row cannot escape this synchronous call.
    pub(crate) fn with_native_word_parts_two_part_source<R>(
        &self,
        source: LebSource<'_>,
        body: i32,
        epoch: TwoPartEpoch,
        cache: &LebChunkCache,
        consume: impl FnOnce(
            Option<crate::positions::CartesianWordParts>,
            &FileIndex,
            Option<&[f64]>,
        ) -> Result<R>,
    ) -> Result<R> {
        if self.header().format() == LebFormat::Leb1
            || !matches!(source, LebSource::Slice(_))
            || cache.word_image_identity().is_none()
        {
            return consume(
                self.eval_raw_word_parts_two_part_source(source, body, epoch, cache)?,
                self,
                None,
            );
        }
        self.with_raw_word_parts_two_part_source(source, body, epoch, cache, |raw, row| {
            consume(raw, self, row)
        })
    }

    fn with_raw_word_parts_two_part_source<R>(
        &self,
        source: LebSource<'_>,
        body: i32,
        epoch: TwoPartEpoch,
        cache: &LebChunkCache,
        consume: impl FnOnce(Option<crate::positions::CartesianWordParts>, Option<&[f64]>) -> Result<R>,
    ) -> Result<R> {
        use crate::positions::derivative::AngularWord as Word;
        let entry = self.find_body(body)?;
        if !matches!(entry.coord_type(), COORD_ICRS_BARY | COORD_ICRS_BARY_SYSTEM)
            || entry.components() != 3
        {
            return consume(None, None);
        }
        let jd = epoch.value();
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        if jd < entry.jd_start() || jd > entry.jd_end() {
            return Err(Error::LebOutOfRange {
                body_id: body,
                jd,
                start: entry.jd_start(),
                end: entry.jd_end(),
            });
        }
        let (segment, _) = select_segment_two_part(
            entry.jd_start(),
            entry.interval_days(),
            entry.segment_count(),
            epoch,
        )?;
        let midpoint = (entry.jd_start() + segment as f64 * entry.interval_days())
            + 0.5 * entry.interval_days();
        let scale = Word::from(2.0) / Word::from(entry.interval_days());
        let raw_tau = (Word::from(epoch.hi) + Word::from(epoch.lo) - Word::from(midpoint)) * scale;
        // Retain the selector's physical [-1, 1] endpoint clamp. A residual
        // never licenses extrapolating the selected polynomial past its edge.
        let tau = if (raw_tau - Word::from(1.0)).rounded() > 0.0 {
            Word::from(1.0)
        } else if (raw_tau + Word::from(1.0)).rounded() < 0.0 {
            Word::from(-1.0)
        } else {
            raw_tau
        };
        let evaluate = |coefficients: &[f64]| {
            let width = entry.degree() as usize + 1;
            // Acquisition, resource admission and the original selected-row
            // shape precede every lookup. An arbitrary public/raw-index cache
            // has no token and always computes the unchanged recurrence.
            let key = if coefficients.len() == 3 * width {
                cache.word_image_identity().and_then(|image| {
                    u32::try_from(segment).ok().map(|record| {
                        crate::pure_cache::RawWordStateKey::new(
                            image,
                            body,
                            record,
                            [epoch.hi.to_bits(), epoch.lo.to_bits()],
                        )
                    })
                })
            } else {
                None
            };
            if let Some(value) = key.and_then(crate::pure_cache::raw_word_state) {
                // Repeat the original first-failing-axis finite/error policy
                // at its ordinary post-acquisition point, even on a hit.
                for (position, velocity) in value.0.iter().zip(value.1.iter()) {
                    let p = Word {
                        hi: position[0],
                        lo: position[1],
                    };
                    let v = Word {
                        hi: velocity[0],
                        lo: velocity[1],
                    };
                    if !p.rounded().is_finite() || !v.rounded().is_finite() {
                        return Err(Error::MathDomain {
                            value: p.rounded() + v.rounded(),
                        });
                    }
                }
                return Ok(Some(value));
            }
            let value = polynomial_xyz_value_rate_words(coefficients, width, tau, scale)?;
            if let Some(key) = key {
                crate::pure_cache::retain_raw_word_state(key, value);
            }
            Ok(Some(value))
        };
        let consume_row = |coefficients: &[f64]| {
            let raw = evaluate(coefficients)?;
            consume(raw, Some(coefficients))
        };
        if self.header().format() == LebFormat::Leb1 {
            consume_row(&read_segment_coefficients(source, entry, segment)?)
        } else {
            self.with_leb2_segment_coefficients(
                source,
                body,
                entry,
                segment,
                jd,
                Some(cache),
                consume_row,
            )
        }
    }
}

impl CompositeIndex<'_> {
    /// Residual-preserving ICRS AU/TT-day state from the same selected image.
    /// Unsupported channels pass None to the consumer; reader and ownership
    /// errors propagate before it runs. Eligible rows are borrowed synchronously.
    pub(crate) fn with_native_word_parts_two_part<R>(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
        consume: impl FnOnce(
            Option<crate::positions::CartesianWordParts>,
            &FileIndex,
            Option<&[f64]>,
        ) -> Result<R>,
    ) -> Result<R> {
        let Some(image) = self.serving_index(body) else {
            return Err(Error::LebBodyNotFound { body_id: body });
        };
        self.evaluation(image)
            .with_native_words(body, epoch, consume)
    }
}

impl TieredIndex<'_> {
    /// Residual-preserving ICRS AU/TT-day state with the existing tier choice.
    /// Unsupported channels pass None to the consumer; selected reader/domain
    /// errors propagate before it runs. No different tier is tried on failure.
    pub(crate) fn with_native_word_parts_two_part<R>(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
        consume: impl FnOnce(
            Option<crate::positions::CartesianWordParts>,
            &FileIndex,
            Option<&[f64]>,
        ) -> Result<R>,
    ) -> Result<R> {
        let candidates = self
            .candidates
            .get(&body)
            .ok_or(Error::LebBodyNotFound { body_id: body })?;
        let selected = self
            .selected(body, epoch.value())
            .unwrap_or(candidates[0].position);
        self.tiers[selected]
            .1
            .with_native_word_parts_two_part(body, epoch, consume)
    }
}

/// Evaluate three equal-width selected Cartesian rows, retaining both words.
/// Independent axes interleave the same Clenshaw operation order; selection,
/// clocks and endpoint normalization remain with the caller. Width is the
/// admitted degree plus one, with no fixed polynomial-degree limit.
fn polynomial_xyz_value_rate_words(
    coefficients: &[f64],
    width: usize,
    tau: crate::positions::derivative::AngularWord,
    scale: crate::positions::derivative::AngularWord,
) -> Result<crate::positions::CartesianWordParts> {
    use crate::positions::derivative::AngularWord as Word;
    if coefficients.len() != 3 * width {
        return Err(Error::LebCorrupt {
            what: "Cartesian state expansion coefficients".into(),
            offset: 0,
        });
    }
    debug_assert_ne!(width, 0);
    let rows: [&[f64]; 3] =
        std::array::from_fn(|axis| &coefficients[axis * width..(axis + 1) * width]);
    let (mut b1, mut b2, mut d1, mut d2) = (
        [Word::from(0.0); 3],
        [Word::from(0.0); 3],
        [Word::from(0.0); 3],
        [Word::from(0.0); 3],
    );
    let twice_tau = Word::from(2.0) * tau;
    for k in (1..width).rev() {
        for (axis, row) in rows.iter().enumerate() {
            let next = twice_tau * b1[axis] + Word::from(row[k]) - b2[axis];
            let derivative = twice_tau * d1[axis] + Word::from(2.0) * b1[axis] - d2[axis];
            (b2[axis], b1[axis], d2[axis], d1[axis]) = (b1[axis], next, d1[axis], derivative);
        }
    }
    let mut position = [[0.0; 2]; 3];
    let mut velocity = [[0.0; 2]; 3];
    // Preserve the original first-failing-axis check and error expression.
    // All earlier interleaved operations are pure per-axis arithmetic.
    for (axis, row) in rows.iter().enumerate() {
        let p = tau * b1[axis] + Word::from(row[0]) - b2[axis];
        let v = (b1[axis] + tau * d1[axis] - d2[axis]) * scale;
        if !p.rounded().is_finite() || !v.rounded().is_finite() {
            return Err(Error::MathDomain {
                value: p.rounded() + v.rounded(),
            });
        }
        position[axis] = [p.hi, p.lo];
        velocity[axis] = [v.hi, v.lo];
    }
    Ok((position, velocity))
}

/// Original serial row sequence used only by exact scheduling regressions.
#[cfg(test)]
fn polynomial_value_rate_words(
    coefficients: &[f64],
    tau: crate::positions::derivative::AngularWord,
    scale: crate::positions::derivative::AngularWord,
) -> (
    crate::positions::derivative::AngularWord,
    crate::positions::derivative::AngularWord,
) {
    use crate::positions::derivative::AngularWord as Word;
    let (mut b1, mut b2, mut d1, mut d2) = (
        Word::from(0.0),
        Word::from(0.0),
        Word::from(0.0),
        Word::from(0.0),
    );
    for coefficient in coefficients.iter().skip(1).rev() {
        let next = Word::from(2.0) * tau * b1 + Word::from(*coefficient) - b2;
        let derivative = Word::from(2.0) * tau * d1 + Word::from(2.0) * b1 - d2;
        (b2, b1, d2, d1) = (b1, next, d1, derivative);
    }
    (
        tau * b1 + Word::from(coefficients[0]) - b2,
        (b1 + tau * d1 - d2) * scale,
    )
}

#[cfg(test)]
mod expanded_source_tests {
    use super::*;
    use crate::positions::derivative::AngularWord as Word;

    #[test]
    fn selected_polynomial_keeps_a_value_residual_below_the_output_ulp() {
        let epsilon = 2.0_f64.powi(-54);
        let (position, rate) =
            polynomial_value_rate_words(&[1.0, epsilon], Word::from(1.0), Word::from(2.0));
        assert_eq!(position.hi, 1.0);
        assert_eq!(position.lo, epsilon);
        assert_eq!(position.rounded(), 1.0);
        assert_eq!(rate.rounded(), 2.0 * epsilon);
    }

    fn original_xyz(
        coefficients: &[f64],
        width: usize,
        tau: Word,
        scale: Word,
    ) -> Result<crate::positions::CartesianWordParts> {
        if coefficients.len() != 3 * width {
            return Err(Error::LebCorrupt {
                what: "Cartesian state expansion coefficients".into(),
                offset: 0,
            });
        }
        let mut position = [[0.0; 2]; 3];
        let mut velocity = [[0.0; 2]; 3];
        for (axis, row) in coefficients.chunks_exact(width).enumerate() {
            let (p, v) = polynomial_value_rate_words(row, tau, scale);
            if !p.rounded().is_finite() || !v.rounded().is_finite() {
                return Err(Error::MathDomain {
                    value: p.rounded() + v.rounded(),
                });
            }
            position[axis] = [p.hi, p.lo];
            velocity[axis] = [v.hi, v.lo];
        }
        Ok((position, velocity))
    }

    #[test]
    fn joint_xyz_retains_every_serial_high_and_low_word() {
        let arguments = [
            Word::from(-1.0),
            Word::from(1.0),
            Word::from(-0.0),
            Word::from(0.0),
            Word::from(f64::from_bits(1)),
            Word::from(-f64::from_bits(1)),
            Word {
                hi: 0.5,
                lo: 2.0_f64.powi(-56),
            },
            Word {
                hi: -0.5,
                lo: -2.0_f64.powi(-56),
            },
            Word {
                hi: f64::from_bits(1.0_f64.to_bits() - 1),
                lo: 2.0_f64.powi(-55),
            },
        ];
        let scales = [
            Word::from(1.0),
            Word::from(0.0625),
            Word::from(2.0) / Word::from(3.0),
        ];
        for degree in (0..=65).chain([127, 256, 1024]) {
            let width = degree + 1;
            for pattern in 0..3 {
                let mut rows = vec![0.0; 3 * width];
                for axis in 0..3 {
                    for k in 0..width {
                        let sign = if (k + degree + axis) % 2 == 0 {
                            1.0
                        } else {
                            -1.0
                        };
                        rows[axis * width + k] = match pattern {
                            0 => {
                                sign * (1.0 + (k + axis) as f64 / 64.0)
                                    * 2.0_f64.powi(-((k * 7 % 48) as i32))
                            }
                            1 => sign * (1e16 + (axis + k) as f64),
                            _ => {
                                if (axis + k) % 2 == 0 {
                                    -0.0
                                } else {
                                    0.0
                                }
                            }
                        };
                    }
                }
                for tau in arguments {
                    for scale in scales {
                        let actual =
                            polynomial_xyz_value_rate_words(&rows, width, tau, scale).unwrap();
                        let expected = original_xyz(&rows, width, tau, scale).unwrap();
                        for (actual, expected) in [(actual.0, expected.0), (actual.1, expected.1)] {
                            for axis in 0..3 {
                                assert_eq!(
                                    actual[axis].map(f64::to_bits),
                                    expected[axis].map(f64::to_bits),
                                    "degree={degree}, pattern={pattern}, axis={axis}, tau=({}, {}), scale=({}, {})",
                                    tau.hi,
                                    tau.lo,
                                    scale.hi,
                                    scale.lo
                                );
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn joint_xyz_preserves_shape_refusal_and_first_axis_error_words() {
        let tau = Word::from(1.0);
        let scale = Word::from(2.0);
        for length in [0, 1, 2, 4, 5, 7, 8] {
            let rows = vec![0.0; length];
            for actual in [
                polynomial_xyz_value_rate_words(&rows, 2, tau, scale),
                original_xyz(&rows, 2, tau, scale),
            ] {
                assert!(
                    matches!(actual, Err(Error::LebCorrupt { what, offset: 0 }) if what == "Cartesian state expansion coefficients")
                );
            }
        }
        for first_bad_axis in 0..3 {
            for invalid in [
                f64::from_bits(0x7ff8_0000_0000_0123),
                f64::INFINITY,
                f64::MAX,
            ] {
                let mut rows = [0.0; 9];
                for axis in first_bad_axis..3 {
                    rows[axis * 3] = invalid;
                    rows[axis * 3 + 1] = if invalid == f64::MAX {
                        f64::MAX
                    } else {
                        axis as f64
                    };
                    rows[axis * 3 + 2] = if invalid.is_nan() {
                        f64::from_bits(0x7ff8_0000_0000_1000 + axis as u64)
                    } else {
                        0.0
                    };
                }
                let actual = polynomial_xyz_value_rate_words(&rows, 3, tau, scale).unwrap_err();
                let expected = original_xyz(&rows, 3, tau, scale).unwrap_err();
                assert_eq!(actual.to_string(), expected.to_string());
                match (actual, expected) {
                    (
                        Error::MathDomain { value: actual },
                        Error::MathDomain { value: expected },
                    ) => {
                        // Rust does not guarantee arithmetic NaN payload/sign preservation:
                        // https://doc.rust-lang.org/std/primitive.f64.html#nan-bit-patterns
                        // Error kind/display and NaN class are contracts; all other
                        // error words remain exact, including infinities and signed zero.
                        assert_eq!(actual.is_nan(), expected.is_nan());
                        if !actual.is_nan() {
                            assert_eq!(
                                actual.to_bits(),
                                expected.to_bits(),
                                "first bad axis={first_bad_axis}"
                            );
                        }
                    }
                    _ => panic!("joint and original first-axis errors must match"),
                }
            }
        }
    }
}
