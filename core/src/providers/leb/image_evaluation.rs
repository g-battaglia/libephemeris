//! Owner-bound evaluation of one immutable index, backing and coefficient cache.
//!
//! Raw public FileIndex adapters keep their caller-managed affinity contract.
//! Engine-owned routes instead obtain this view from a file, a composite slot,
//! or a binding that created its own cache. No caller can supply an unrelated
//! persistent cache when constructing one of these internal evaluation views.

use super::*;

/// Retained binding used by the single-image provider. These are the same
/// index/source/cache components and ownership as the provider's raw API.
pub(crate) struct ImageBinding<'a> {
    index: &'a FileIndex,
    source: LebSource<'a>,
    cache: Arc<LebChunkCache>,
}

impl<'a> ImageBinding<'a> {
    /// Raw compatibility construction always starts a cache for this pair.
    /// The public caller still owns the index/bytes correspondence contract.
    pub(crate) fn from_raw(index: &'a FileIndex, bytes: &'a [u8]) -> Self {
        Self {
            index,
            source: LebSource::Slice(bytes),
            cache: Arc::new(LebChunkCache::new()),
        }
    }

    /// Borrow one admitted file owner; retain only its existing cache Arc.
    /// A closed handle fails before any numerical input is consumed.
    pub(crate) fn from_file(file: &'a LebFile) -> Result<Self> {
        let image = file.image.as_ref().ok_or(Error::LebReaderClosed)?;
        Ok(Self {
            index: &image.index,
            source: image.bytes.source(),
            cache: Arc::clone(&image.cache),
        })
    }

    /// Metadata borrows retain the original file/index lifetime.
    pub(crate) fn index(&self) -> &'a FileIndex {
        self.index
    }

    /// Stack-only read access; no lookup, Arc clone or cache admission occurs.
    pub(crate) fn evaluation(&self) -> ImageEvaluation<'_> {
        ImageEvaluation {
            index: self.index,
            source: self.source,
            cache: &self.cache,
        }
    }
}

/// Borrowed selection whose fields cannot be assembled outside this module.
/// It owns no data or budget and performs no reads until an operation consumes
/// it. TT epochs, channel units, record selection and typed errors are those
/// of the existing FileIndex evaluators; there is no new numerical route.
pub(crate) struct ImageEvaluation<'a> {
    index: &'a FileIndex,
    source: LebSource<'a>,
    cache: &'a LebChunkCache,
}

impl LebFile {
    /// Select the backing/index/cache from this one immutable owner.
    pub(crate) fn evaluation(&self) -> Result<ImageEvaluation<'_>> {
        let image = self.image.as_ref().ok_or(Error::LebReaderClosed)?;
        Ok(ImageEvaluation::from_image(image))
    }
}

impl CompositeIndex<'_> {
    /// Select one already admitted slot, never an independently supplied cache.
    /// Raw constructor entries have their dedicated caches; file-backed entries
    /// take all components directly from the retained immutable image owner.
    pub(super) fn evaluation(&self, slot: usize) -> ImageEvaluation<'_> {
        match &self.entries {
            CompositeSources::Raw(entries) => ImageEvaluation {
                index: &entries[slot].0,
                source: entries[slot].1,
                cache: &self.caches[slot],
            },
            CompositeSources::Shared(files) => {
                let image = files[slot].image.as_ref().expect("validated LEB owner");
                ImageEvaluation::from_image(image)
            }
        }
    }
}

impl<'a> ImageEvaluation<'a> {
    fn from_image(image: &'a LebImage) -> Self {
        Self {
            index: &image.index,
            source: image.bytes.source(),
            cache: &image.cache,
        }
    }

    /// Preserve explicit LEB1 admission rather than silently dispatching formats.
    pub(crate) fn eval_leb1(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.index.eval_body_source(self.source, body, jd)
    }

    /// Split-TT form of explicit LEB1 evaluation.
    pub(crate) fn eval_leb1_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.index
            .eval_body_two_part_source(self.source, body, epoch)
    }

    /// Preserve explicit LEB2 admission and this image's bounded decode cache.
    pub(crate) fn eval_leb2(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.index
            .eval_leb2_body_cached_source(self.source, body, jd, self.cache)
    }

    /// Split-TT form of explicit LEB2 evaluation.
    pub(crate) fn eval_leb2_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.index
            .eval_leb2_body_two_part_cached_source(self.source, body, epoch, self.cache)
    }

    /// One detected-format dispatch shared by file, single and composite routes.
    pub(crate) fn eval_any(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        if self.index.header().format() == LebFormat::Leb1 {
            self.eval_leb1(body, jd)
        } else {
            self.eval_leb2(body, jd)
        }
    }

    /// Detected-format dispatch retaining the original additive TT parts.
    pub(crate) fn eval_any_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        if self.index.header().format() == LebFormat::Leb1 {
            self.eval_leb1_two_part(body, epoch)
        } else {
            self.eval_leb2_two_part(body, epoch)
        }
    }

    /// Inline projection of the same scalar stored-channel evaluation.
    pub(crate) fn eval_inline(&self, body: i32, jd: f64) -> Result<Option<CartesianState>> {
        self.index
            .eval_state_inline_source(self.source, body, jd, self.cache)
    }

    /// Inline projection retaining split-epoch record ownership.
    pub(crate) fn eval_inline_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<CartesianState>> {
        self.index
            .eval_state_inline_two_part_source(self.source, body, epoch, self.cache)
    }

    /// Stored acceleration in AU/TT-day squared, without center composition.
    pub(crate) fn acceleration(&self, body: i32, epoch: TwoPartEpoch) -> Result<Option<[f64; 3]>> {
        self.index
            .eval_raw_acceleration_two_part_source(self.source, body, epoch, self.cache)
    }

    /// Preserve aligned-polynomial refusal and joint-state errors unchanged.
    pub(crate) fn relative_state(
        &self,
        target: i32,
        center: i32,
        jd: f64,
    ) -> Result<Option<CartesianState>> {
        self.index
            .eval_relative_icrs_cached_source(self.source, target, center, jd, self.cache)
    }

    /// Relative AU position with independently owned target/center TT parts.
    pub(crate) fn relative_position(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        self.index.eval_relative_position_two_part_cached_source(
            self.source,
            target,
            target_epoch,
            center,
            center_epoch,
            self.cache,
        )
    }

    /// Relative AU/TT-day velocity with the same selected channel pair.
    pub(crate) fn relative_velocity(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        self.index.eval_relative_velocity_two_part_cached_source(
            self.source,
            target,
            target_epoch,
            center,
            center_epoch,
            self.cache,
        )
    }

    /// Lend an eligible acquired row synchronously; no lock crosses the callback.
    pub(crate) fn with_native_words<R>(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
        consume: impl FnOnce(
            Option<crate::positions::CartesianWordParts>,
            &FileIndex,
            Option<&[f64]>,
        ) -> Result<R>,
    ) -> Result<R> {
        self.index.with_native_word_parts_two_part_source(
            self.source,
            body,
            epoch,
            self.cache,
            consume,
        )
    }

    /// Stored nutation in radians, using the selected immutable image.
    pub(crate) fn nutation(&self, jd: f64) -> Result<(f64, f64)> {
        self.index.eval_nutation_source(self.source, jd)
    }

    /// Exact stored nutation derivative, radians per TT day.
    pub(crate) fn nutation_rate(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        self.index.nutation_rate_source(self.source, jd)
    }

    /// Selected right-owned nutation polynomial bounds, in TT Julian days.
    pub(crate) fn nutation_interval(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        self.index.nutation_smooth_interval_source(self.source, jd)
    }

    /// Header-only presence probe; backing/resource failures remain fallible.
    pub(crate) fn has_nutation(&self) -> Result<bool> {
        self.index.has_nutation_data_source(self.source)
    }

    /// Stored TT-minus-UT1 interpolation, with no clock conversion.
    pub(crate) fn delta_t(&self, jd: f64) -> Result<f64> {
        self.index.delta_t_source(self.source, jd)
    }

    /// Catalog row from the same selected image; no provider lookup is added.
    pub(crate) fn star(&self, star: i32) -> Result<StarEntry> {
        self.index.get_star_source(self.source, star)
    }
}
