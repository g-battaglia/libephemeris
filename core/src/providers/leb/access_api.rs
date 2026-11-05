//! Compatibility adapters for callers supplying contiguous container bytes.
//! Numerical implementations live in the source-based methods of `FileIndex`.

use super::*;

impl FileIndex {
    /// Advisory file ranges for finite TT Julian-day bounds. Checked offset
    /// overflow fails; malformed compressed directories are skipped because
    /// warming never changes numerical data or selects another provider.
    pub fn warm_ranges(&self, bytes: &[u8], start: f64, end: f64) -> Result<Vec<(u64, u64)>> {
        self.warm_ranges_source(LebSource::Slice(bytes), start, end)
    }

    /// Subtract aligned Cartesian ICRS polynomials before rounding their large
    /// common origin. Returns AU/AU-day at TT Julian day, or `None` for
    /// unaligned records; missing bodies and invalid stored data fail.
    /// `bytes`, this index and `cache` must obey [`LebChunkCache`]'s immutable
    /// image-affinity contract; clearing a cache does not change its image.
    pub fn eval_relative_icrs_cached(
        &self,
        bytes: &[u8],
        target: i32,
        center: i32,
        jd: f64,
        cache: &LebChunkCache,
    ) -> Result<Option<CartesianState>> {
        self.eval_relative_icrs_cached_source(LebSource::Slice(bytes), target, center, jd, cache)
    }

    /// Paired Clenshaw difference in ICRS AU, retaining independently owned
    /// two-part TT epochs. System-barycenter channels return `None`; missing
    /// bodies, invalid epochs, coverage and decode errors propagate.
    /// The raw bytes/index/cache must belong to one immutable image, under
    /// [`LebChunkCache`]'s lifetime-long affinity contract.
    pub fn eval_relative_position_two_part_cached(
        &self,
        bytes: &[u8],
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
        cache: &LebChunkCache,
    ) -> Result<Option<[f64; 3]>> {
        self.eval_relative_position_two_part_cached_source(
            LebSource::Slice(bytes),
            target,
            target_epoch,
            center,
            center_epoch,
            cache,
        )
    }

    /// Evaluate an LEB1 body at TT Julian day; returns stored components and
    /// rates per TT day. The returned tuple is `(components, rates)` with
    /// stored-coordinate units: ICRS AU/AU-day or ecliptic deg/deg/AU and their
    /// daily rates, as identified by [`BodyEntry::coord_type`]. The bytes must
    /// match this decoded index. Missing bodies, invalid epochs, coverage and
    /// malformed coefficients fail with the corresponding typed reader error.
    pub fn eval_body(&self, bytes: &[u8], body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_body_source(LebSource::Slice(bytes), body, jd)
    }

    /// Split-TT form of [`Self::eval_body`], retaining the sub-ulp epoch tail.
    /// Record ownership follows the collapsed day; units and errors are unchanged.
    pub fn eval_body_two_part(
        &self,
        bytes: &[u8],
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_body_two_part_source(LebSource::Slice(bytes), body, epoch)
    }

    /// Evaluate compressed LEB2 coefficients at TT Julian day, with no retained
    /// decode cache. Units and domain errors match [`Self::eval_body`]; invalid
    /// compressed geometry, frames or chunk ownership are corruption errors.
    pub fn eval_leb2_body(&self, bytes: &[u8], body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_leb2_body_source(LebSource::Slice(bytes), body, jd)
    }

    /// [`Self::eval_leb2_body`] using this immutable image's decoded cache.
    /// `cache` must remain paired with these immutable bytes and index for its
    /// entire lifetime; see [`LebChunkCache`]. A hit reuses decoded words, not
    /// an evaluation at a different epoch. Decode/resource errors can still occur
    /// on misses; this adapter does not establish image identity.
    pub fn eval_leb2_body_cached(
        &self,
        bytes: &[u8],
        body: i32,
        jd: f64,
        cache: &LebChunkCache,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_leb2_body_cached_source(LebSource::Slice(bytes), body, jd, cache)
    }

    /// Split-TT form of [`Self::eval_leb2_body`]; ownership follows the collapsed
    /// day and normalized polynomial coordinates retain its fractional tail.
    pub fn eval_leb2_body_two_part(
        &self,
        bytes: &[u8],
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_leb2_body_two_part_source(LebSource::Slice(bytes), body, epoch)
    }

    /// Split-TT compressed evaluation sharing an image-specific decoded cache.
    /// Units and errors match [`Self::eval_leb2_body_two_part`]. The same
    /// immutable-image affinity as [`Self::eval_leb2_body_cached`] is required.
    pub fn eval_leb2_body_two_part_cached(
        &self,
        bytes: &[u8],
        body: i32,
        epoch: TwoPartEpoch,
        cache: &LebChunkCache,
    ) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_leb2_body_two_part_cached_source(LebSource::Slice(bytes), body, epoch, cache)
    }

    /// Stored `(dpsi, deps)` in radians at TT Julian day. Missing or malformed
    /// tables, nonfinite epochs and dates outside their stored coverage fail.
    pub fn eval_nutation(&self, bytes: &[u8], jd: f64) -> Result<(f64, f64)> {
        self.eval_nutation_source(LebSource::Slice(bytes), jd)
    }

    /// Split-TT form of [`Self::eval_nutation`], with the same units and errors.
    pub fn eval_nutation_two_part(&self, bytes: &[u8], epoch: TwoPartEpoch) -> Result<(f64, f64)> {
        self.eval_nutation_two_part_source(LebSource::Slice(bytes), epoch)
    }

    /// Right-owned nutation polynomial interval in TT days. Missing/uncovered
    /// data return `None`; malformed table metadata propagates its error.
    pub fn nutation_smooth_interval(&self, bytes: &[u8], jd: f64) -> Result<Option<(f64, f64)>> {
        self.nutation_smooth_interval_source(LebSource::Slice(bytes), jd)
    }

    /// Analytical TT derivative of stored nutation in radians/day. Missing or
    /// uncovered data return `None`; malformed tables fail. Uses the same
    /// right-owned polynomial as [`Self::eval_nutation`].
    pub fn nutation_rate(&self, bytes: &[u8], jd: f64) -> Result<Option<(f64, f64)>> {
        self.nutation_rate_source(LebSource::Slice(bytes), jd)
    }

    /// Stored TT−UT1 in days: right-biased linear interpolation, clamped at
    /// endpoints. `jd` uses the stored table's epoch coordinate; this raw reader
    /// performs no clock conversion or freshness check. Nonfinite requests and
    /// empty/malformed tables are errors. See [`eval_delta_t`] for interpolation.
    pub fn delta_t(&self, bytes: &[u8], jd: f64) -> Result<f64> {
        self.delta_t_source(LebSource::Slice(bytes), jd)
    }

    /// Catalog record with J2000 position and documented stellar units. Last
    /// complete duplicate wins, short tails are ignored; absent IDs or malformed
    /// complete records return their typed reader errors.
    pub fn get_star(&self, bytes: &[u8], star: i32) -> Result<StarEntry> {
        self.get_star_source(LebSource::Slice(bytes), star)
    }
}
