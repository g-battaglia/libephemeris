//! Routing implementation for the parent module.
//! Public entry points and types are exposed through the parent API.

use super::*;

/// Acquired ICRS AU/AU-per-TDB-day state and separate position-polynomial tangent.
type PctrStateWithTdbTangent = ([f64; 3], [f64; 3], Option<[f64; 3]>);

#[path = "model_serving.rs"]
pub(super) mod model_serving;
pub use model_serving::{serving_model_position, serving_model_position_with_models};

#[path = "provider_policy.rs"]
mod provider_policy;
pub use provider_policy::{
    decide_jpl_file, escalate_sealed_provider_miss, forbid_jpl_in_sealed_mode, is_leb_core_body,
    is_leb_model_served_point, is_leb_stored_state_body, is_recognized_de_kernel,
    leb_miss_fallthrough, leb_miss_outcome, plan_route, require_leb, select_de_kernel,
    serving_compute_tt_escalated, tier_de_kernel,
};

#[path = "native_jpl_entry.rs"]
mod native_jpl_entry;
pub(super) use native_jpl_entry::PredictedSpkReads;
use native_jpl_entry::with_local_jpl_entry_outputs_and_fallback_and_range;
pub use native_jpl_entry::{
    calc_native_jpl_entry_with_acquisitions, calc_native_jpl_entry_with_main_image,
    calc_native_jpl_entry_with_minor_assist, calc_native_jpl_entry_with_models,
};

/// Per-call native dispatch hint containing the physical provider label and
/// optional original body alias. It owns no astronomical resource.
struct NativeTraceHint {
    requested: Option<i32>,
    source: Option<crate::trace::TraceSource>,
    serving_source: Option<&'static str>,
    snapshot_fallback: bool,
}

std::thread_local! {
    /// Native hints are separate from the legacy tracer and invocation event
    /// observer. Nested native calls restore their caller's current hint.
    static NATIVE_TRACE_HINTS: std::cell::RefCell<Vec<NativeTraceHint>> = const {
        std::cell::RefCell::new(Vec::new())
    };
}

/// Retain an invocation-local dispatch hint across all output/speed stages.
/// Dropping or unwinding the call restores the enclosing scope. No borrow is
/// held while provider, model or caller work runs.
struct NativeTraceInvocation {
    depth: usize,
}

impl NativeTraceInvocation {
    /// A real public entry retains its original identifier. A generic native
    /// callback factory creates a scope only when no public scope exists.
    fn enter(requested: Option<i32>) -> Self {
        let depth = NATIVE_TRACE_HINTS.with(|cell| {
            let mut stack = cell.borrow_mut();
            let depth = stack.len();
            stack.push(NativeTraceHint {
                requested,
                source: None,
                serving_source: None,
                snapshot_fallback: true,
            });
            depth
        });
        Self { depth }
    }

    fn enter_factory() -> Option<Self> {
        let absent = NATIVE_TRACE_HINTS.with(|cell| cell.borrow().is_empty());
        absent.then(|| Self::enter(None))
    }
}

impl Drop for NativeTraceInvocation {
    fn drop(&mut self) {
        // Hints contain only owned text and scalar values: dropping them runs
        // no caller callback. Restore the parent in place without allocating
        // a temporary suffix; nested scopes and unwinding retain the same depth.
        NATIVE_TRACE_HINTS.with(|cell| cell.borrow_mut().truncate(self.depth));
    }
}

#[cfg(test)]
#[path = "tests/native_trace_hints.rs"]
mod native_trace_hint_tests;

/// Run a known stored or analytical serving leg through the shared public tail.
///
/// `source` belongs to the selected fast/fallback branch, not the configured
/// provider mode or returned values. Fast LEB success has no fallback snapshots;
/// analytical/native retries retain the public fallback's caller/alias rollback.
/// The shared entry publishes only after representation and SPEED3 complete,
/// and transfers a recursive north-node publication to its original south id.
/// Nested calls and unwinding restore the enclosing plain-Rust hint scope.
pub(crate) fn with_serving_trace_source<R>(
    body: i32,
    source: &'static str,
    fast: bool,
    operation: impl FnOnce() -> Result<R>,
) -> Result<R> {
    let _scope = NativeTraceInvocation::enter(Some(body));
    NATIVE_TRACE_HINTS.with(|cell| {
        let mut stack = cell.borrow_mut();
        let hint = stack.last_mut().expect("the invocation owns its hint");
        hint.serving_source = Some(source);
        hint.snapshot_fallback = !fast;
    });
    operation()
}

/// Record a direct stored serving success unless a public branch owns its tail.
///
/// Explicit stored/model entry scopes defer every inner LEB publication until
/// the full result, speed and frame work succeeds. This also prevents support
/// reads from replacing the actual model/fallback label or publishing a
/// remapped identifier. Direct legacy serving calls outside such a scope retain
/// their existing immediate LEB record. No observer callback runs under a hint
/// borrow, and no numerical result or failure is changed.
pub(super) fn record_stored_trace_source(body: i32) -> Result<()> {
    let deferred = NATIVE_TRACE_HINTS.with(|cell| {
        cell.borrow()
            .last()
            .is_some_and(|hint| hint.serving_source.is_some())
    });
    if !deferred {
        crate::trace::try_record(body, "LEB")?;
    }
    Ok(())
}

/// Retain the provider selected by a consumed native branch. Publish success
/// only after the complete output/frame/rate calculation has succeeded.
pub(super) fn mark_native_trace_source(source: &'static str) {
    mark_native_trace_source_owned(crate::trace::TraceSource::named(source));
}

/// Retain an actual owned callback source until the completed public tail.
/// The acquired helper's text is never replaced by a mode/row label.
pub(super) fn mark_native_trace_source_owned(source: crate::trace::TraceSource) {
    let retired = NATIVE_TRACE_HINTS.with(|cell| {
        let mut stack = cell.borrow_mut();
        match stack.last_mut() {
            Some(hint) => hint.source.replace(source),
            None => Some(source),
        }
    });
    drop(retired);
}

/// Publish a completed native body reduction after its output/speed stages.
/// Outside a native hint scope this is a no-op: injected/stored engines are
/// never assigned a guessed native label. Builtin asteroid remaps retain the
/// original caller id; a recursive north-node result keeps its own north id.
pub(super) fn publish_native_trace_source(body: i32) -> Result<()> {
    if let Some((body, source)) = take_native_trace_source(body) {
        crate::trace::try_record_source(body, &source)?;
    }
    Ok(())
}

/// Transfer the completed branch's physical source to the trace recorder.
/// Release the local hint borrow before publishing the value; an absent hint
/// stays absent.
pub(super) fn take_native_trace_source(body: i32) -> Option<(i32, crate::trace::TraceSource)> {
    NATIVE_TRACE_HINTS.with(|cell| {
        let mut stack = cell.borrow_mut();
        let hint = stack.last_mut()?;
        let source = hint.source.take()?;
        let requested = hint
            .requested
            .filter(|id| remap_ast_offset(*id) == body)
            .unwrap_or(body);
        Some((requested, source))
    })
}

/// Snapshot native fallback entries immediately before their original epoch /
/// body/output calculation. The successful result is unchanged; on error the
/// remapped entry is restored first, then the original caller, exactly as the
/// public `calc`/`calc_ut` fallback catch. Recursive north calls own independent
/// snapshots. No snapshot or restore is emitted for a non-native engine.
pub(super) fn with_native_trace_body<R>(
    requested: i32,
    served: i32,
    operation: impl FnOnce() -> Result<R>,
) -> Result<R> {
    let policy = NATIVE_TRACE_HINTS.with(|cell| {
        let mut stack = cell.borrow_mut();
        let hint = stack.last_mut()?;
        let retired = std::mem::replace(
            &mut hint.source,
            hint.serving_source.map(crate::trace::TraceSource::named),
        );
        Some((hint.snapshot_fallback, retired))
    });
    let Some((snapshot_fallback, retired)) = policy else {
        return operation();
    };
    drop(retired);
    if !snapshot_fallback {
        return operation();
    }
    let requested_trace = crate::trace::try_snapshot_record(requested)?;
    let alias_trace = crate::trace::try_snapshot_record(served)?;
    let result = operation();
    if result.is_err() {
        crate::trace::try_restore_record(served, alias_trace)?;
        crate::trace::try_restore_record(requested, requested_trace)?;
    }
    result
}

/// Set a branch hint only after its real calculation succeeds, preserving the
/// exact owned success or failure rather than classifying its display text.
fn native_source_result<R>(result: Result<R>, source: &'static str) -> Result<R> {
    if result.is_ok() {
        mark_native_trace_source(source);
    }
    result
}

impl StateFrame {
    /// Classify a stored `coord_type`. Any other value is
    /// [`Error::LebCorrupt`] (decode accepts any `u32`, but the
    /// calculation path only serves the five inventoried channels).
    pub fn from_coord_type(coord_type: u32) -> Result<Self> {
        match coord_type {
            COORD_ICRS_BARY => Ok(Self::IcrsBarycentric),
            COORD_ICRS_BARY_SYSTEM => Ok(Self::IcrsSystemBarycenter),
            COORD_ECLIPTIC => Ok(Self::EclipticDirect),
            COORD_GEO_ECLIPTIC => Ok(Self::GeocentricEcliptic),
            COORD_HELIO_ECL_RETIRED => Ok(Self::HeliocentricEclipticRetired),
            other => Err(Error::LebCorrupt {
                what: format!("LEB body channel (unknown coord_type {other})"),
                offset: 0,
            }),
        }
    }

    /// True for the two ICRS xyz channels.
    pub fn is_icrs(self) -> bool {
        matches!(self, Self::IcrsBarycentric | Self::IcrsSystemBarycenter)
    }

    /// Position component units, in stored order.
    pub fn position_units(self) -> [&'static str; 3] {
        match self {
            Self::IcrsBarycentric | Self::IcrsSystemBarycenter => ["AU", "AU", "AU"],
            Self::EclipticDirect | Self::GeocentricEcliptic | Self::HeliocentricEclipticRetired => {
                ["deg", "deg", "AU"]
            }
        }
    }

    /// Velocity component units, in stored order.
    pub fn velocity_units(self) -> [&'static str; 3] {
        match self {
            Self::IcrsBarycentric | Self::IcrsSystemBarycenter => ["AU/day", "AU/day", "AU/day"],
            Self::EclipticDirect | Self::GeocentricEcliptic | Self::HeliocentricEclipticRetired => {
                ["deg/day", "deg/day", "AU/day"]
            }
        }
    }
}

impl GeometricState {
    /// Requested body identifier, retained for downstream coverage errors.
    pub fn body_id(self) -> i32 {
        self.body_id
    }

    /// TT Julian Day the state was evaluated at.
    pub fn jd_tt(self) -> f64 {
        self.jd_tt
    }

    /// Stored position: AU xyz for ICRS channels, or
    /// `[lon_deg, lat_deg, dist_au]` for ecliptic channels.
    pub fn pos(self) -> [f64; 3] {
        self.pos
    }

    /// Stored analytical velocity in the matching per-day units.
    pub fn vel(self) -> [f64; 3] {
        self.vel
    }

    /// Stored channel frame (component meaning and units).
    pub fn frame(self) -> StateFrame {
        self.frame
    }

    /// True only when a runtime planet-center offset was applied to a
    /// system-barycenter channel. `false` with
    /// [`StateFrame::IcrsSystemBarycenter`] means the stored system
    /// barycenter is retained as the explicit data-availability fallback
    /// (the LEB body-composition layer), never a silent substitution.
    pub fn is_planet_center(self) -> bool {
        self.planet_center
    }

    /// Tier that served this state, when the provider selects across
    /// tiers (`TieredIndex`); `None` for single/composite providers.
    pub fn tier(self) -> Option<LebTier> {
        self.tier
    }
}

impl EntrySnapshot {
    fn of(entry: &BodyEntry) -> Self {
        Self {
            coord_type: entry.coord_type(),
        }
    }
}

/// Assemble a geometric state from an evaluated stored channel.
///
/// Exactly three components are required; other counts are [`Error::LebCorrupt`].
/// System-barycenter channels add an available physical center correction,
/// including its centered-difference velocity. Missing corrections retain
/// the stored system-barycenter frame; other frames pass through unchanged.
pub(super) fn build_state(
    body_id: i32,
    jd_tt: f64,
    entry: &EntrySnapshot,
    pos: &[f64],
    vel: &[f64],
    tier: Option<LebTier>,
    offset_at: impl Fn(f64) -> Option<[f64; 3]>,
) -> Result<GeometricState> {
    let frame = StateFrame::from_coord_type(entry.coord_type)?;
    let ([p0, p1, p2], [v0, v1, v2]) = match (pos, vel) {
        ([p0, p1, p2], [v0, v1, v2]) => ([*p0, *p1, *p2], [*v0, *v1, *v2]),
        _ => {
            return Err(Error::LebCorrupt {
                what: format!("LEB stored state (body {body_id})"),
                offset: 0,
            });
        }
    };
    let is_system_bary = is_system_barycenter_channel(entry.coord_type);
    let (pos3, vel3) =
        eval_body_center([p0, p1, p2], [v0, v1, v2], is_system_bary, jd_tt, offset_at);
    // `eval_body_center` returns the stored array untouched on the
    // fallback/fast paths, so exact array comparison distinguishes
    // "offset applied" from "barycenter retained" with no extra flag
    // plumbing and no float tolerance.
    let planet_center = is_system_bary && (pos3 != [p0, p1, p2] || vel3 != [v0, v1, v2]);
    Ok(GeometricState {
        body_id,
        jd_tt,
        pos: pos3,
        vel: vel3,
        frame,
        planet_center,
        tier,
    })
}

impl<'a> SingleLeb<'a> {
    /// Bind an already-decoded index to its image, with an empty
    /// decoded-chunk memo (see the `SingleLeb` type docs).
    pub fn new(index: &'a FileIndex, bytes: &'a [u8]) -> Self {
        Self {
            image: crate::leb::ImageBinding::from_raw(index, bytes),
        }
    }

    /// Bind a file's immutable index/bytes and shared image-specific cache.
    /// Rebuilding the provider or cloning the context preserves decoded words;
    /// a closed file fails. No image copy or independent cache is allocated.
    pub fn from_file(file: &'a LebFile) -> Result<Self> {
        Ok(Self {
            image: crate::leb::ImageBinding::from_file(file)?,
        })
    }

    /// The bound index (inventory consumers).
    pub fn index(&self) -> &'a FileIndex {
        self.image.index()
    }
}

impl StateProvider for SingleLeb<'_> {
    fn has_analytic_kinematics(&self) -> bool {
        true
    }

    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        self.image
            .evaluation()
            .with_native_words(body, epoch, |raw, index, row| match row {
                Some(coefficients) => native_state_word_parts_from(self, body, epoch, raw, || {
                    index.eval_borrowed_state_inline_two_part(body, epoch, coefficients)
                }),
                None => native_state_word_parts(self, body, epoch, raw),
            })
    }

    fn geometric_acceleration_tt(&self, body: i32, jd: f64) -> Result<Option<[f64; 3]>> {
        self.image
            .evaluation()
            .acceleration(body, TwoPartEpoch::split(jd)?)
    }

    fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        self.image
            .evaluation()
            .relative_velocity(target, target_epoch, center, center_epoch)
    }

    fn can_reuse_scalar_state(&self, body: i32, jd: f64, epoch: TwoPartEpoch) -> bool {
        self.image
            .index()
            .scalar_state_matches_split(body, jd, epoch)
    }

    fn eval_relative_position_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        self.image
            .evaluation()
            .relative_position(target, target_epoch, center, center_epoch)
    }

    fn eval_relative_icrs(
        &self,
        target: i32,
        center: i32,
        jd: f64,
    ) -> Result<Option<crate::leb::CartesianState>> {
        self.image.evaluation().relative_state(target, center, jd)
    }

    fn provider_name(&self) -> &'static str {
        "sealed-leb/single"
    }

    fn has_body(&self, body_id: i32) -> bool {
        self.image.index().has_body(body_id)
    }

    fn body_coverage(&self, body_id: i32) -> Option<(f64, f64)> {
        self.image.index().body_coverage(body_id)
    }

    fn smooth_interval_tt(&self, body_id: i32, jd_tt: f64) -> Result<Option<(f64, f64)>> {
        self.image.index().body_smooth_interval(body_id, jd_tt)
    }

    fn nutation_smooth_interval_tt(&self, jd_tt: f64) -> Result<Option<(f64, f64)>> {
        self.image.evaluation().nutation_interval(jd_tt)
    }

    fn serving_coord_type(&self, body_id: i32) -> Option<u32> {
        self.image
            .index()
            .find_body(body_id)
            .ok()
            .map(|e| e.coord_type())
    }

    fn eval_raw_inline(&self, body_id: i32, jd: f64) -> Result<Option<crate::leb::CartesianState>> {
        self.image.evaluation().eval_inline(body_id, jd)
    }

    fn eval_raw_inline_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<crate::leb::CartesianState>> {
        self.image.evaluation().eval_inline_two_part(body_id, epoch)
    }

    fn eval_raw(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.image.evaluation().eval_any(body_id, jd)
    }

    fn eval_raw_two_part(&self, body_id: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.image.evaluation().eval_any_two_part(body_id, epoch)
    }

    fn serving_entry_snapshot(&self, body_id: i32) -> Result<EntrySnapshot> {
        self.image.index().find_body(body_id).map(EntrySnapshot::of)
    }
}

impl StateProvider for CompositeIndex<'_> {
    fn has_analytic_kinematics(&self) -> bool {
        true
    }

    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        self.with_native_word_parts_two_part(body, epoch, |raw, index, row| match row {
            Some(coefficients) => native_state_word_parts_from(self, body, epoch, raw, || {
                index.eval_borrowed_state_inline_two_part(body, epoch, coefficients)
            }),
            None => native_state_word_parts(self, body, epoch, raw),
        })
    }

    fn nutation_rate_tt(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        self.nutation_rate(jd)
    }

    fn geometric_acceleration_tt(&self, body: i32, jd: f64) -> Result<Option<[f64; 3]>> {
        self.eval_raw_acceleration_two_part(body, TwoPartEpoch::split(jd)?)
    }

    fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        CompositeIndex::eval_relative_velocity_two_part(
            self,
            target,
            target_epoch,
            center,
            center_epoch,
        )
    }

    fn can_reuse_scalar_state(&self, body: i32, jd: f64, epoch: TwoPartEpoch) -> bool {
        self.scalar_state_matches_split(body, jd, epoch)
    }

    fn eval_relative_position_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        CompositeIndex::eval_relative_position_two_part(
            self,
            target,
            target_epoch,
            center,
            center_epoch,
        )
    }

    fn eval_relative_icrs(
        &self,
        target: i32,
        center: i32,
        jd: f64,
    ) -> Result<Option<crate::leb::CartesianState>> {
        CompositeIndex::eval_relative_icrs(self, target, center, jd)
    }

    fn provider_name(&self) -> &'static str {
        "sealed-leb/composite"
    }

    fn has_body(&self, body_id: i32) -> bool {
        CompositeIndex::has_body(self, body_id)
    }

    fn body_coverage(&self, body_id: i32) -> Option<(f64, f64)> {
        CompositeIndex::body_coverage(self, body_id)
    }

    fn smooth_interval_tt(&self, body_id: i32, jd_tt: f64) -> Result<Option<(f64, f64)>> {
        CompositeIndex::body_smooth_interval(self, body_id, jd_tt)
    }

    fn nutation_smooth_interval_tt(&self, jd_tt: f64) -> Result<Option<(f64, f64)>> {
        CompositeIndex::nutation_smooth_interval(self, jd_tt)
    }

    fn serving_coord_type(&self, body_id: i32) -> Option<u32> {
        self.serving_entry(body_id).map(|e| e.coord_type())
    }

    fn eval_raw_inline(&self, body_id: i32, jd: f64) -> Result<Option<crate::leb::CartesianState>> {
        self.eval_state_inline(body_id, jd)
    }

    fn eval_raw_inline_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<crate::leb::CartesianState>> {
        self.eval_state_inline_two_part(body_id, epoch)
    }

    fn eval_raw(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_body_any(body_id, jd)
    }

    fn eval_raw_two_part(&self, body_id: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_body_any_two_part(body_id, epoch)
    }

    fn serving_entry_snapshot(&self, body_id: i32) -> Result<EntrySnapshot> {
        self.serving_entry(body_id)
            .map(EntrySnapshot::of)
            .ok_or(Error::LebBodyNotFound { body_id })
    }
}

impl StateProvider for TieredIndex<'_> {
    fn has_analytic_kinematics(&self) -> bool {
        true
    }

    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        self.with_native_word_parts_two_part(body, epoch, |raw, index, row| match row {
            Some(coefficients) => native_state_word_parts_from(self, body, epoch, raw, || {
                index.eval_borrowed_state_inline_two_part(body, epoch, coefficients)
            }),
            None => native_state_word_parts(self, body, epoch, raw),
        })
    }

    fn nutation_rate_tt(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        self.nutation_rate(jd)
    }

    fn geometric_acceleration_tt(&self, body: i32, jd: f64) -> Result<Option<[f64; 3]>> {
        self.eval_raw_acceleration_two_part(body, TwoPartEpoch::split(jd)?)
    }

    fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        TieredIndex::eval_relative_velocity_two_part(
            self,
            target,
            target_epoch,
            center,
            center_epoch,
        )
    }

    fn can_reuse_scalar_state(&self, body: i32, jd: f64, epoch: TwoPartEpoch) -> bool {
        self.scalar_state_matches_split(body, jd, epoch)
    }

    fn eval_relative_position_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        TieredIndex::eval_relative_position_two_part(
            self,
            target,
            target_epoch,
            center,
            center_epoch,
        )
    }

    fn eval_relative_icrs(
        &self,
        target: i32,
        center: i32,
        jd: f64,
    ) -> Result<Option<crate::leb::CartesianState>> {
        TieredIndex::eval_relative_icrs(self, target, center, jd)
    }

    fn provider_name(&self) -> &'static str {
        "sealed-leb/tiered"
    }

    fn has_body(&self, body_id: i32) -> bool {
        TieredIndex::has_body(self, body_id)
    }

    fn body_coverage(&self, body_id: i32) -> Option<(f64, f64)> {
        TieredIndex::body_coverage(self, body_id)
    }

    fn smooth_interval_tt(&self, body_id: i32, jd_tt: f64) -> Result<Option<(f64, f64)>> {
        TieredIndex::body_smooth_interval(self, body_id, jd_tt)
    }

    fn nutation_smooth_interval_tt(&self, jd_tt: f64) -> Result<Option<(f64, f64)>> {
        TieredIndex::nutation_smooth_interval(self, jd_tt)
    }

    fn serving_coord_type(&self, body_id: i32) -> Option<u32> {
        self.serving_entry(body_id).map(|e| e.coord_type())
    }

    fn eval_raw_inline(&self, body_id: i32, jd: f64) -> Result<Option<crate::leb::CartesianState>> {
        self.eval_state_inline(body_id, jd)
    }

    fn eval_raw_inline_two_part(
        &self,
        body_id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<crate::leb::CartesianState>> {
        self.eval_state_inline_two_part(body_id, epoch)
    }

    fn eval_raw(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_body_any(body_id, jd)
    }

    fn eval_raw_two_part(&self, body_id: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        self.eval_body_any_two_part(body_id, epoch)
    }

    fn selected_tier(&self, body_id: i32, jd: f64) -> Option<LebTier> {
        TieredIndex::selected_tier(self, body_id, jd)
    }

    fn serving_entry_snapshot(&self, body_id: i32) -> Result<EntrySnapshot> {
        self.serving_entry(body_id)
            .map(EntrySnapshot::of)
            .ok_or(Error::LebBodyNotFound { body_id })
    }
}

/// Assemble one sealed-serving stack over `ctx` and run `f` over the
/// borrowed view.
///
/// Domain: contexts in [`ProviderMode::SealedLeb`] or
/// [`ProviderMode::Auto`] with at least one bound LEB image; the Delta-T
/// model and IERS selections arrive already resolved (binding-side
/// acquisition, network-policy precedence pattern). `skip_closed` selects the closed-image
/// policy: `false` fails fast with the typed reader error (the native
/// [`with_sealed_calc_engine`] contract); `true` counts closed images as
/// absent and lets the empty remainder refuse via [`CompositeIndex::new`]
/// (the bindings sealed-path contract — the bindings re-point adopts this
/// arm, so both policies stay explicit here rather than silently merged).
/// `NativeJpl`/`HorizonsRemote` modes refuse with typed
/// `non_leb_leg_wiring`; sealed/Auto-LEB never opens JPL kernels
/// (`jpl_fallback_active: false` downstream).
pub fn with_sealed_serving_stack<'c, R>(
    ctx: &'c EphemerisContext,
    delta_t_model: DeltaTModel,
    iers_enabled: bool,
    skip_closed: bool,
    f: impl FnOnce(&SealedServingStack<'_, 'c>) -> Result<R>,
) -> Result<R> {
    let mode = ctx.mode();
    match mode {
        ProviderMode::SealedLeb | ProviderMode::Auto => {}
        ProviderMode::NativeJpl | ProviderMode::HorizonsRemote => {
            return Err(Error::ConfigurationMissing {
                missing_config: "non_leb_leg_wiring",
                suggestion: "only sealed/Auto-with-LEB serving is bound yet (POS-02/POS-09 own the native-JPL/Horizons legs)",
            });
        }
    }
    if !skip_closed {
        require_leb(!ctx.leb_files().is_empty())?;
    }
    let composite = ctx.leb_composite(TierPolicy::AttachOrder, skip_closed)?;
    let nutation_at = |jd: f64| composite.eval_nutation(jd);
    let serving = SealedServing::new_stored(composite, &nutation_at);
    let delta_cfg = ctx.delta_t_config(delta_t_model, iers_enabled);
    let stack = SealedServingStack {
        mode,
        composite,
        serving: &serving,
        nutation_at: &nutation_at,
        delta_cfg: &delta_cfg,
    };
    f(&stack)
}

/// Native-contract [`CalcEngine`] assembly over `ctx`: the
/// [`with_shared_calc_engine`] core with the native [`ServingRequest`]
/// (no explicit observer, [`TierPolicy::AttachOrder`], fail-fast
/// closed-image policy).
///
/// Domain: [`ProviderMode::SealedLeb`] or [`ProviderMode::Auto`] with
/// bound LEB images. Model/IERS selections arrive already resolved. Other
/// modes refuse before callback construction with `non_leb_leg_wiring`;
/// this explicitly sealed factory never consults SPK or substitutes another
/// provider. Use [`with_shared_calc_engine`] for a NativeJpl request.
/// This metadata-free wrapper treats site samples as TT entries; use
/// [`with_shared_calc_engine_at_entry`] with a sealed request and the original
/// UT1 scalar for a topocentric UT calculation.
pub fn with_sealed_calc_engine<R>(
    ctx: &EphemerisContext,
    delta_t_model: DeltaTModel,
    iers_enabled: bool,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    match ctx.mode() {
        ProviderMode::SealedLeb | ProviderMode::Auto => {}
        ProviderMode::NativeJpl | ProviderMode::HorizonsRemote => {
            return Err(Error::ConfigurationMissing {
                missing_config: "non_leb_leg_wiring",
                suggestion: "the sealed factory requires a sealed/Auto LEB context; use with_shared_calc_engine for local SPK",
            });
        }
    }
    with_shared_calc_engine(
        ServingRequest {
            ctx,
            delta_t_model,
            iers_enabled,
            observer: None,
            tier_policy: TierPolicy::AttachOrder,
            skip_closed: false,
        },
        f,
    )
}

/// Use the reader's canonical tier-token grammar, including original LEB1
/// `ephemeris_{tier}[_group]` names. Directory names never supply the tier.
/// An unrecognized basename stays explicit: no attach-order fallback.
pub(super) fn tier_token_of(path: &str) -> Option<LebTier> {
    file_tier(path.rsplit(['/', '\\']).next().unwrap_or(path))
}

/// Map the context tier to its LEB tier for the installed-set check.
pub(super) fn active_leb_tier(tier: Tier) -> LebTier {
    match tier {
        Tier::Base => LebTier::Base,
        Tier::Medium => LebTier::Medium,
        Tier::Extended => LebTier::Extended,
    }
}

/// Context-to-engine assembly for explicitly selected local providers.
///
/// SealedLeb and Auto-with-LEB requests select image priority ([`TierPolicy`]),
/// build the first-wins composite and stored-nutation feed, and bind the
/// context Delta-T chain. Missing/closed images follow the request's policy.
/// Topocentric serving uses the explicit request observer in preference to
/// the context site; an absent required site retains its typed refusal.
///
/// NativeJpl requests instead bind the context's registered local SPK chains,
/// without inspecting or substituting bound LEB images. Construction is lazy;
/// an actual core calculation reports missing target registrations, coverage
/// misses or sidereal anchor failures. NativeJpl UT
/// entries acquire TT from the independent bundled timescale, while explicit
/// TT entries retain the supplied epoch. Horizons remains an explicitly
/// unbound leg with `non_leb_leg_wiring`. No branch performs network I/O.
///
/// This metadata-free factory pairs site samples as TT entries. Use
/// [`with_shared_calc_engine_at_entry`] with the original UT1 scalar for a
/// topocentric UT entry; an inverse cannot reconstruct that word in general.
/// The named [`with_sealed_calc_engine`] wrapper deliberately restricts this
/// general factory to its original sealed/Auto LEB domain.
pub fn with_shared_calc_engine<R>(
    req: ServingRequest<'_>,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    with_shared_calc_engine_at_entry(req, None, f)
}

/// Assemble one calculation engine with immutable original UT1 metadata.
///
/// `original_ut` is `Some(U)` for one UT1 entry and `None` for TT. For LEB
/// topocentric callbacks the central sample retains U exactly and SPEED uses
/// its frozen central TT-U difference. TT samples invert configured Delta-T
/// with one refinement, `U = T - D(T - D(T))`; D is in days. Observer resolution
/// precedes this site-clock lookup in the TT callback. A caller of `calc_ut`
/// still owns any required observer preflight before the engine's forward
/// time callback; factory construction has no body or flags to validate.
///
/// Provider selection and synchronous borrowed callback lifetimes are those
/// of [`with_shared_calc_engine`]. One factory invocation belongs to one
/// public epoch: do not reuse Some(U) for another UT entry. NativeJpl uses
/// the independent bundled timescale and existing original-UT range/model
/// metadata. No LEB clock replaces native, lunar, frame or defining-epoch
/// clocks. Canonical Delta-T failures propagate with their original typed
/// identity. For a TT site inverse, only a NaN calendar argument under
/// Espenak-Meeus with IERS disabled and no override uses the selected reader's
/// Delta-T table; no unclassified failure activates that alternative.
pub fn with_shared_calc_engine_at_entry<R>(
    req: ServingRequest<'_>,
    original_ut: Option<f64>,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    with_shared_calc_outputs_at_entry(req, original_ut, |engine, _| f(engine))
}

/// Assemble spherical and direct Cartesian context callbacks.
///
/// Selection, observer precedence, clocks and borrowed ownership follow
/// [`with_shared_calc_engine`]. Stored ICRS channels can return direct
/// Cartesian AU/AU-day components; other channels return `None` for spherical
/// shaping. Native JPL uses its spherical public reduction. Provider failures
/// remain errors, and assembly performs no network I/O. Use
/// [`with_shared_calc_outputs_at_entry`] to retain an original UT1 entry.
pub fn with_shared_calc_outputs<R>(
    req: ServingRequest<'_>,
    f: impl FnOnce(&CalcEngine<'_>, &CartesianCompute<'_>) -> Result<R>,
) -> Result<R> {
    with_shared_calc_outputs_at_entry(req, None, f)
}

/// Bind spherical and direct Cartesian callbacks for one original TT/UT entry.
///
/// Inputs, immutable entry metadata, time-scale domains, borrows, observer
/// priority and failures follow [`with_shared_calc_engine_at_entry`]. Stored
/// ICRS channels provide direct selected-frame `(x,y,z,vx,vy,vz)` in AU and
/// AU/day; other channels return None for the existing spherical reducer.
/// Both callbacks use the same configured site UT1 partner and unchanged
/// WGS84/ERFA rotation. The topocentric reducer freezes the central TT-UT1
/// difference across its internal SPEED stencil. The classified TT calendar
/// alternative reads the actual composite's Delta-T table at original TT.
/// All other time, coverage, channel, frame and rotation errors retain their
/// typed identity; a consumed table failure propagates.
/// NativeJpl keeps its spherical output channel and independent clocks.
/// Returned six-slot shaping/retflags are owned by the calculation entry.
/// There is no retained mutable clock state or network acquisition here.
pub fn with_shared_calc_outputs_at_entry<R>(
    req: ServingRequest<'_>,
    original_ut: Option<f64>,
    f: impl FnOnce(&CalcEngine<'_>, &CartesianCompute<'_>) -> Result<R>,
) -> Result<R> {
    let ctx = req.ctx;
    let mode = ctx.mode();
    match mode {
        ProviderMode::NativeJpl => return with_native_jpl_entry_outputs(req, original_ut, f),
        ProviderMode::SealedLeb | ProviderMode::Auto => {}
        ProviderMode::HorizonsRemote => {
            return Err(Error::ConfigurationMissing {
                missing_config: "non_leb_leg_wiring",
                suggestion: "bind the configured Horizons vector client to the apparent-position pipeline",
            });
        }
    }
    if !req.skip_closed {
        require_leb(!ctx.leb_files().is_empty())?;
    }
    let composite = ctx.leb_composite(req.tier_policy, req.skip_closed)?;
    let nutation_at = |jd: f64| composite.eval_nutation(jd);
    let serving = SealedServing::new_stored(composite, &nutation_at);
    let delta_cfg = ctx.delta_t_config(req.delta_t_model, req.iers_enabled);
    let observer = req.observer.or_else(|| ctx.observer());
    let site = observer.map(|obs| (obs.lon_deg(), obs.lat_deg(), obs.height_m()));
    let ut1_to_tt_jd = |tjdut: f64| Ok(tjdut + deltat_with_config(tjdut, &delta_cfg)?);
    let c2t_at = |jd_tt: f64, jd_ut1: f64| crate::erfa::c2t06a(jd_tt, jd_ut1, 0.0, 0.0);
    let compute_in_frame = |jd_tt: f64, body_id: i32, calc_flags: i32, xyz: bool| {
        let geocentric = || {
            if xyz {
                serving_compute_xyz_tt(&serving, jd_tt, body_id, calc_flags).map_err(|err| {
                    escalate_sealed_provider_miss(mode, composite, body_id, jd_tt, err)
                })
            } else {
                serving_compute_tt_escalated(mode, &serving, jd_tt, body_id, calc_flags)
            }
        };
        if calc_flags & FLG_TOPOCTR != 0 {
            match resolve_topo_request(calc_flags, body_id, site, None)? {
                TopoPlan::Geocentric | TopoPlan::TopoNoop => geocentric(),
                TopoPlan::Topocentric {
                    lon_deg,
                    lat_deg,
                    height_m,
                } => {
                    let partner = configured_site_ut1_partner_with_table(
                        jd_tt,
                        original_ut,
                        &delta_cfg,
                        &|jd| composite.delta_t(jd),
                    )?;
                    let ut1_at = |tt| {
                        super::site_clock::configured_site_sample_ut1(
                            tt,
                            original_ut,
                            &delta_cfg,
                            &|jd| composite.delta_t(jd),
                        )
                    };
                    super::serving::serving_compute_topocentric_with_clock(
                        &serving,
                        jd_tt,
                        partner,
                        body_id,
                        calc_flags,
                        (lon_deg, lat_deg, height_m),
                        &c2t_at,
                        xyz,
                        &ut1_at,
                        if calc_flags & FLG_SPEED != 0 {
                            super::site_clock::configured_site_clock_support(
                                jd_tt,
                                original_ut,
                                &delta_cfg,
                            )?
                        } else {
                            (f64::NEG_INFINITY, f64::INFINITY)
                        },
                    )
                    .map_err(|err| {
                        escalate_sealed_provider_miss(mode, composite, body_id, jd_tt, err)
                    })
                }
            }
        } else {
            geocentric()
        }
    };
    let compute_tt = |jd, body, flags| compute_in_frame(jd, body, flags, false);
    let cartesian_tt = |jd, body, flags| match serving_channel_gate(composite, body) {
        Err(Error::SealedChannelNotServed { .. }) => Ok(None),
        Err(err) => Err(escalate_sealed_provider_miss(
            mode, composite, body, jd, err,
        )),
        Ok(_) => compute_in_frame(jd, body, flags, true).map(Some),
    };
    let engine = CalcEngine {
        ut1_to_tt_jd: &ut1_to_tt_jd,
        compute_tt: &compute_tt,
        nutation_obliquity_ut: &ecl_nut_result_ut,
        nutation_obliquity_tt: &|jd_tt: f64, nut_flags: i32| {
            ecl_nut_result(jd_tt, nut_flags, &crate::erfa::nut06a)
        },
        jpl_fallback_active: false,
    };
    f(&engine, &cartesian_tt)
}

/// Bind native position callbacks with an explicit original public UT1 epoch.
///
/// `original_ut` is `Some` only for a public UT entry and retains its scalar
/// Julian-day word; a TT entry passes `None`. The engine keeps its existing
/// bundled UT1-to-TT conversion, SPK split-TDB sampling, flag normalization and
/// speed stencils. A classic model's ecliptic sidereal tail lazily compares its
/// current TT word with the forward conversion of the original UT: the central
/// sample uses that original UT unchanged, while a fresh TT sample derives its
/// own UT partner. Tropical and equatorial tails acquire no unused inverse.
/// Entry metadata is immutable and local; no implicit clock state is retained.
///
/// Public scalar-range validation and promoted raw SPK coverage failures use
/// the original UT requested JD, including a derivative support read failure.
/// An already typed analytical-model range keeps its own TT requested JD.
/// Bodies and owning-image envelopes remain attached to their acquisition.
/// TT entries and raw provider clocks/errors keep their TT/TDB contracts.
/// Registry, observer, clock, model and frame
/// failures otherwise propagate as in the shared native factory. Outputs remain
/// six public coordinates/rates and retflags. Classic model routes perform no
/// network acquisition; other body routes retain the configured network policy.
pub fn with_native_jpl_entry_outputs<R>(
    req: ServingRequest<'_>,
    original_ut: Option<f64>,
    f: impl FnOnce(&CalcEngine<'_>, &CartesianCompute<'_>) -> Result<R>,
) -> Result<R> {
    with_local_jpl_entry_outputs(req, original_ut, f)
}

/// Calculate a native position at a TT or UT1 Julian day.
///
/// Returns six coordinates/rates plus the returned flags. Positions are degrees
/// and AU, or AU for Cartesian output; velocities use the corresponding daily
/// units. The context owns all ephemeris images. Invalid flags, non-finite time,
/// missing data, coverage failures and numerical errors are returned unchanged.
pub fn calc_native_jpl_entry(
    req: ServingRequest<'_>,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
) -> Result<(PositionResult, i32)> {
    calc_native_jpl_entry_with_main_image(None, None, None, req, epoch, body, flags, ut_entry, None)
}

/// Select the UT1 consumed by a native classic model's sidereal tail.
///
/// Both inputs are scalar Julian days. Only an exact binary64 match with the
/// bundled forward UT1-to-TT conversion identifies the public central sample;
/// side samples use a fresh TT inverse. Callers invoke this only when UT1 is
/// actually consumed, so tropical output introduces no extra clock failures.
/// Timescale domain and non-finite errors propagate without fallback.
fn native_model_ut_at(tt: f64, original_ut: Option<f64>) -> Result<f64> {
    if let Some(ut) = original_ut {
        if crate::time::default_timescale_tt_from_ut1(ut)?.to_bits() == tt.to_bits() {
            return Ok(ut);
        }
    }
    crate::time::default_timescale_ut1_from_tt(tt)
}

/// Assemble native SPK calculation callbacks over an owned context.
///
/// Registered NAIF chains supply ICRS AU/AU-day states. TT is retained; UT1
/// uses the bundled observation timescale. Configured Delta-T/IERS remains
/// attached to nested model and sidereal defining epochs. The public entry
/// applies output shaping after the scalar place calculation.
///
/// Callbacks are lazy: construction does not read a body or require a kernel.
/// Missing aliases, coverage and frame errors arise when their state is used.
/// Data-independent nutation output needs no DE image. Unregistered asteroids
/// may obtain SBDB elements under the effective network policy.
fn with_local_jpl_calc_outputs<R>(
    req: ServingRequest<'_>,
    f: impl FnOnce(&CalcEngine<'_>, &CartesianCompute<'_>) -> Result<R>,
) -> Result<R> {
    with_local_jpl_entry_outputs(req, None, f)
}

/// Assemble native callbacks while retaining immutable public-entry metadata.
/// The optional UT1 JD supplies consumed classic sidereal clocks and the
/// public range epoch. Raw SPK samples and independent surface observer rotation remain
/// attached to their TT/TDB times; typed model ranges keep their owning TT.
/// Callback lifetimes are bounded by `f`, and owned registry images remain
/// borrowed through the request context. Acquisition errors remain typed.
fn with_local_jpl_entry_outputs<R>(
    req: ServingRequest<'_>,
    original_ut: Option<f64>,
    f: impl FnOnce(&CalcEngine<'_>, &CartesianCompute<'_>) -> Result<R>,
) -> Result<R> {
    with_local_jpl_entry_outputs_and_fallback(req, original_ut, false, f)
}

/// Assemble native callbacks with a resolved main-filename fallback decision.
///
/// `fallback_active` is true only after an absent unrecognized requested file
/// was replaced by its tier main kernel. The acquisition context owns bytes
/// and registrations throughout `f`; a supplied registry uses the legacy false
/// decision. Only public flag normalization consumes this bit. Original UT,
/// scalar/split provider sampling, models, frames and raw error clocks keep
/// their existing contracts. Callback references never escape the call.
pub(super) fn with_local_jpl_entry_outputs_and_fallback<R>(
    req: ServingRequest<'_>,
    original_ut: Option<f64>,
    fallback_active: bool,
    f: impl FnOnce(&CalcEngine<'_>, &CartesianCompute<'_>) -> Result<R>,
) -> Result<R> {
    with_local_jpl_entry_outputs_and_fallback_and_range(req, original_ut, fallback_active, false, f)
}

/// Apply a fixed-epoch or SIDBIT frame to a native sidereal calculation.
///
/// The finite observation epoch is TT JD. Nested tropical output is spherical
/// degrees/AU and daily rates in J2000 mean axes. ECL_T0 uses the selected
/// model's defining plane; SSY_PLANE projects onto the invariable plane.
/// When both bits are set, ECL_T0 takes precedence. Equatorial ECL_T0 requires
/// only its mean equatorial plane. Configured Delta-T supplies user-UT and
/// Valens defining epochs. Model, clock, frame and provider errors propagate;
/// a frame rotation is never replaced by a scalar longitude offset.
fn native_sidereal_frame_request(
    tropical: &dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>,
    ayanamsha: &crate::sidereal::SealedAyanamshaServing<'_>,
    selection: crate::sidereal::SiderealSelection,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    use crate::sidereal::SiderealRoute;
    let mode = selection.mode();
    let bits = selection.bits();
    let route = crate::sidereal::resolve_sidereal_route(flags, mode, bits);
    if !matches!(
        route,
        Some(
            SiderealRoute::FixedEpoch
                | SiderealRoute::SidbitEcliptic
                | SiderealRoute::SidbitEquatorialT0
        )
    ) {
        return tropical(tt, body, flags);
    }
    (|| {
        let sub_flags = strip_output_flags(crate::sidereal::fixed_epoch_request_flags(flags));
        if route == Some(SiderealRoute::FixedEpoch) {
            let (position, retflag) = tropical(tt, body, sub_flags)?;
            return Ok((
                crate::sidereal::transform_fixed_epoch_result(position, flags, mode)?,
                crate::sidereal::fixed_epoch_retflag(retflag, flags),
            ));
        }
        let projection_t0 = || -> Result<f64> {
            if bits & crate::constants::SIDBIT_ECL_T0 == 0 {
                return Ok(crate::constants::J2000);
            }
            let user_t0 = if mode == crate::constants::SIDM_USER {
                ayanamsha.user_pair_tt()?.0
            } else {
                0.0
            };
            let valens_t0 = if mode == crate::constants::SIDM_VALENS_MOON {
                ayanamsha.valens_epoch_tt()?
            } else {
                0.0
            };
            Ok(crate::sidereal::ecl_t0_epoch_jd(mode, user_t0, valens_t0))
        };
        if route == Some(SiderealRoute::SidbitEquatorialT0) {
            let (position, retflag) = tropical(tt, body, sub_flags)?;
            return Ok((
                crate::sidereal::transform_equatorial_epoch_result(
                    position,
                    flags,
                    projection_t0()?,
                )?,
                crate::sidereal::fixed_epoch_retflag(retflag, flags),
            ));
        }
        let t0 = projection_t0()?;
        let zero_tt = crate::time::default_timescale_tt_from_ut1(t0)?;
        let (mean, _) = ayanamsha.mean_tt_with(zero_tt, false, false)?;
        let zero = if bits & crate::constants::SIDBIT_ECL_T0 != 0 {
            mean
        } else {
            crate::sidereal::ssy_plane_zero_point_deg(mean)?
        };
        let matrix = crate::sidereal::sidbit_ecliptic_matrix(bits, t0, zero)?.ok_or(
            Error::ConfigurationMissing {
                missing_config: "sidereal_projection_matrix",
                suggestion: "the selected native route requires a SIDBIT projection matrix",
            },
        )?;
        let (position, retflag) = tropical(tt, body, sub_flags)?;
        Ok((
            crate::sidereal::transform_sidbit_result(position, flags, matrix)?,
            crate::sidereal::sidbit_ecliptic_retflag(retflag, flags),
        ))
    })()
}

/// Calculate a context-bound Horizons position at scalar TT.
///
/// The remote route subtracts configured Delta-T at TT, then recalculates it
/// at the resulting UT before scalar TDB-labelled requests. Output and errors
/// follow the public calc contract. Unsupported remote requests and connection/
/// lookup failures use only the explicitly owned native registry; sealed policy
/// is never relaxed. No LEB/SPK identity is substituted for remote vectors.
pub fn calc_horizons(
    req: ServingRequest<'_>,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    context_horizons_position(req, tt, body, flags, false, None, None)
}

/// Calculate a context-bound Horizons position at the caller's scalar UT.
///
/// This entry retains the exact UT word for frame/model reductions and computes
/// query TT with configured Delta-T once. It must not be implemented by feeding
/// its computed TT through the approximate public TT inverse. Native fallback
/// uses the native provider's independent bundled UT1-to-TT channel instead.
pub fn calc_horizons_ut(
    req: ServingRequest<'_>,
    ut: f64,
    body: i32,
    flags: i32,
) -> Result<(PositionResult, i32)> {
    context_horizons_position(req, ut, body, flags, true, None, None)
}

/// Native position handoff retaining the current recursive scalar operands.
/// Epoch is TT or UT1 according to the final boolean; body and flags are the
/// owning fallback branch's pre-normalization operands. The synchronous
/// callback returns six public coordinates/rates and flags or its actual error.
type NativePositionFallback<'a> = dyn Fn(f64, i32, i32, bool) -> Result<(PositionResult, i32)> + 'a;

/// Complete a remote route with the Session's current native acquisitions.
/// The original remote clocks, cache, policy, aliases and catch classification
/// remain in the shared implementation. Only an actual existing native
/// fallback branch calls `fallback`; no completed LEB attempt is replayed.
/// Recursive south-node requests pass their north-body operands before the
/// existing alias transform. Callback failures propagate unchanged.
pub(crate) fn calc_horizons_with_fallback(
    req: ServingRequest<'_>,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
    fallback: &NativePositionFallback<'_>,
) -> Result<(PositionResult, i32)> {
    context_horizons_position(req, epoch, body, flags, ut_entry, None, Some(fallback))
}

/// Public flags/aliases surround the complete scalar remote reduction. Each
/// calculation owns its entry clock, so nested requests and SPEED3 samples do
/// not depend on thread-local state or an unconsumed mutable clock marker.
/// Direct Horizons reduction at UT1 JD with an explicitly owned client.
///
/// Output follows the selected physical frames and flags; positions are AU
/// and angles degrees unless the flags request another layout. Unsupported
/// combinations return an explicit refusal, never a local-provider fallback.
/// Invalid times, model domains and transport failures remain typed errors.
pub fn calc_horizons_backend_ut(
    req: ServingRequest<'_>,
    client: &crate::horizons::HorizonsClient,
    ut: f64,
    body: i32,
    flags: i32,
    coverage: Option<crate::lunar::LunarCoverage>,
) -> Result<super::horizons_serving::HorizonsOutcome> {
    if flags & FLG_TOPOCTR != 0 {
        return Ok(super::horizons_serving::HorizonsOutcome::Fallback(
            "FLG_TOPOCTR not supported by Horizons backend",
        ));
    }
    if flags & FLG_ICRS != 0 {
        return Ok(super::horizons_serving::HorizonsOutcome::Fallback(
            "FLG_ICRS not supported by Horizons backend",
        ));
    }
    horizons_backend_position(req, client, ut, body, flags, true, None, coverage)
}

#[allow(clippy::too_many_arguments)]
fn horizons_backend_position(
    req: ServingRequest<'_>,
    client: &crate::horizons::HorizonsClient,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
    fetch_override: Option<&super::horizons_serving::HorizonsFetch<'_>>,
    coverage_override: Option<crate::lunar::LunarCoverage>,
) -> Result<super::horizons_serving::HorizonsOutcome> {
    let delta_cfg = req.ctx.delta_t_config(req.delta_t_model, req.iers_enabled);
    let delta = |ut| deltat_with_config(ut, &delta_cfg);
    let states = LocalSpkStateProvider::new(req.ctx.spk_registry(), req.ctx.mode())?;
    let anchors = crate::stars::SealedAnchorProvider::new(&states);
    let (selection, _) = crate::houses::sidereal_selection_for_context(req.ctx);
    let ayanamsha = crate::sidereal::SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
    let aya = |ut, flags| crate::planetary_moons::moon_spk_ayanamsha(ut, flags, &ayanamsha);
    let coverage = if let Some(coverage) = coverage_override {
        coverage
    } else if crate::horizons::is_analytical_body(body) && flags & (FLG_HELCTR | FLG_BARYCTR) == 0 {
        local_jpl_lunar_coverage(
            &JplLeg::new(req.ctx.spk_registry()),
            body,
            if ut_entry {
                epoch + delta(epoch)?
            } else {
                let ut = epoch - delta(epoch)?;
                ut + delta(ut)?
            },
        )?
    } else {
        crate::lunar::LunarCoverage::Unenforced
    };
    let numerics =
        super::horizons_output::HorizonsNumerics::new(&delta, &aya, selection.full().0, coverage);
    let output = |p, v, tt, ut, flags| numerics.output(p, v, tt, ut, flags).map(|value| value.0);
    let local = |ut, tt, body, flags| numerics.local(ut, tt, body, flags).map(|value| value.0);
    let serving = super::horizons_serving::HorizonsServing {
        context: req.ctx,
        client,
        fetch_override,
        delta_days: &delta,
        output: &output,
        local: &local,
    };
    if ut_entry {
        serving.calc_ut(epoch, body, flags)
    } else {
        serving.calc_tt(epoch, body, flags)
    }
}

fn context_horizons_position(
    req: ServingRequest<'_>,
    epoch: f64,
    body: i32,
    flags: i32,
    ut_entry: bool,
    fetch_override: Option<&super::horizons_serving::HorizonsFetch<'_>>,
    native_fallback: Option<&NativePositionFallback<'_>>,
) -> Result<(PositionResult, i32)> {
    forbid_jpl_in_sealed_mode(req.ctx.mode())?;
    if !epoch.is_finite() {
        return Err(Error::NonFiniteArgument { value: epoch });
    }
    validate_barycentric_moseph(flags, if ut_entry { "calc_ut" } else { "calc" })?;
    let fallback = || {
        let _trace_scope = NativeTraceInvocation::enter(Some(body));
        if let Some(native) = native_fallback {
            return native(epoch, body, flags, ut_entry);
        }
        with_local_jpl_calc_outputs(req, |engine, _| {
            if ut_entry {
                calc_ut(engine, epoch, body, flags)
            } else {
                calc(engine, epoch, body, flags)
            }
        })
    };
    if body == ECL_NUT {
        return fallback();
    }
    let raw_flags = consume_non_echoed_flags(flags, body);
    let flags = normalize_calc_flags(raw_flags, false);
    let requested_body = body;
    let body = remap_ast_offset(body);
    if body == -MEAN_NODE || body == -TRUE_NODE {
        let previous = crate::trace::try_snapshot_record(-body)?;
        let (north, retflag) = context_horizons_position(
            req,
            epoch,
            -body,
            flags,
            ut_entry,
            fetch_override,
            native_fallback,
        )?;
        crate::trace::try_record_alias(body, -body, previous)?;
        let output = south_node_from_north(north, flags);
        // TT output echoes only requested ephemeris bits; UT1 output includes its default.
        let retflag = echo_request_bits(retflag, raw_flags);
        return Ok((
            output,
            if ut_entry {
                retflag
            } else {
                calc_tt_epheflag_echo(retflag, raw_flags)
            },
        ));
    }
    let Some(client) = req.ctx.horizons_client() else {
        return fallback();
    };
    let outcome = horizons_backend_position(
        req,
        client,
        epoch,
        body,
        flags,
        ut_entry,
        fetch_override,
        None,
    );
    match outcome {
        Ok(super::horizons_serving::HorizonsOutcome::Served(position, retflag)) => {
            let source = if crate::horizons::is_analytical_body(body)
                || crate::horizons::is_hypothetical_body(body)
                || body == SUN && flags & FLG_HELCTR != 0
            {
                "Analytical"
            } else {
                "Horizons"
            };
            crate::trace::try_record(requested_body, source)?;
            let retflag = echo_request_bits(retflag | implied_retflag_bits(flags), raw_flags);
            Ok((
                position,
                if ut_entry {
                    retflag
                } else {
                    calc_tt_epheflag_echo(retflag, raw_flags)
                },
            ))
        }
        Ok(super::horizons_serving::HorizonsOutcome::Fallback(_)) => fallback(),
        Err(error @ Error::NetworkSealed { .. }) => Err(error),
        Err(error) if horizons_fallback_error(&error) => fallback(),
        Err(error) => Err(error),
    }
}

/// Exception-family classification for the public remote-to-native handoff.
/// Lookup/transport failures and numerical ValueError counterparts may retry
/// through the declared native pipeline. Configuration, sealed policy and
/// custom calculation/range errors do not become a permissive fallback.
fn horizons_fallback_error(error: &Error) -> bool {
    matches!(
        error,
        Error::HorizonsApiError { .. }
            | Error::HorizonsMalformed { .. }
            | Error::HorizonsRequestFailed { .. }
            | Error::HttpTransport { .. }
            | Error::HttpTimeout { .. }
            | Error::HttpStatus { .. }
            | Error::UnknownBody { .. }
            | Error::SpkFileNotFound { .. }
            | Error::LebFileNotFound { .. }
            | Error::NonFiniteArgument { .. }
            | Error::ZeroVectorNormalize
            | Error::ZeroVectorToSpherical
            | Error::SingularMatrix
            | Error::NormalizedTimeOutOfRange { .. }
            | Error::InvalidSegment { .. }
            | Error::InvalidDerivativeOrder { .. }
    )
}

/// Read the configured main DE image's metadata without sampling a state.
///
/// The NativeJpl context identifies its main planetary image through the Sun
/// registration. A separately registered Moon image does not replace this
/// selection. All descriptors contribute to the file envelope, in scalar JD;
/// neither a target directory lookup nor SPK evaluation is performed. Missing
/// registration/image errors propagate, with no download or arbitrary image
/// selection. The returned window does not itself validate a model epoch.
fn local_jpl_main_window(leg: &JplLeg<'_>) -> Result<(f64, f64)> {
    let key = local_jpl_target(leg, SUN)?;
    let info = leg
        .registry()
        .body_info(key)
        .ok_or(Error::SpkUnregisteredBody { body: SUN })?;
    let kernel = leg
        .registry()
        .kernel(&info.kernel)
        .ok_or_else(|| Error::SpkFileNotFound {
            path: info.kernel.clone(),
        })?;
    let (start, end) =
        kernel
            .index()
            .all_segments()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |span, segment| {
                (
                    span.0.min(segment.descriptor.start_et),
                    span.1.max(segment.descriptor.end_et),
                )
            });
    let (start_jd, end_jd) = (FRAME_J2000 + start / 86400.0, FRAME_J2000 + end / 86400.0);
    Ok((start_jd, end_jd))
}

/// Validate a mean lunar model's TT against the main DE file envelope.
/// Both scalar descriptor endpoints are included. The already typed model
/// range retains this TT JD even when its public caller entered at UT1.
fn local_jpl_lunar_coverage(
    leg: &JplLeg<'_>,
    body: i32,
    tt: f64,
) -> Result<crate::lunar::LunarCoverage> {
    let (start_jd, end_jd) = local_jpl_main_window(leg)?;
    if tt < start_jd || tt > end_jd {
        return Err(Error::EphemerisRange {
            requested_jd: tt,
            start_jd,
            end_jd,
            body: Some(body),
        });
    }
    Ok(crate::lunar::LunarCoverage::Window { start_jd, end_jd })
}

/// Promote provider coverage failures to the public calculation domain.
///
/// The requested epoch remains the caller's central TT or UT1 even if an
/// emitted or derivative sample failed. The selected image's descriptor
/// envelope supplies bounds; lunar models use Moon support and hypothetical
/// models use Sun support. Gaps remain enforced by the underlying reader.
/// Convergence and non-coverage errors retain their original type.
fn local_jpl_public_range_error(leg: &JplLeg<'_>, tt: f64, body: i32, error: Error) -> Error {
    let provider_span = match &error {
        Error::SpkOutOfRange { start, end, .. } => Some((*start, *end)),
        _ => return error,
    };
    let acquisition_body = if crate::lunar::is_lunar_point(body) {
        MOON
    } else if (CUPIDO..=WALDEMATH).contains(&body) {
        SUN
    } else {
        body
    };
    let kernel = local_jpl_target(leg, acquisition_body)
        .ok()
        .and_then(|key| leg.registry().body_info(key))
        .and_then(|info| leg.registry().kernel(&info.kernel));
    let main_span = kernel.map(|kernel| {
        kernel
            .index()
            .all_segments()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |span, segment| {
                (
                    span.0.min(segment.descriptor.start_et),
                    span.1.max(segment.descriptor.end_et),
                )
            })
    });
    let Some((start, end)) = main_span.or(provider_span) else {
        return error;
    };
    Error::EphemerisRange {
        requested_jd: tt,
        start_jd: FRAME_J2000 + start / 86400.0,
        end_jd: FRAME_J2000 + end / 86400.0,
        body: Some(body),
    }
}

/// Geometric Moon-minus-Earth state in the true ecliptic of date.
///
/// A split TT epoch is converted to TDB for all SPK constituents. Position
/// and physical velocity, in AU and AU/day, are rotated using IAU 2000A
/// nutation, P03 precession and true obliquity; no frame-motion derivative is
/// added to this osculating state. Missing aliases, corrupt states, invalid
/// frames and actual record-coverage errors propagate.
fn local_jpl_true_moon_state(leg: &JplLeg<'_>, jd: f64) -> Result<crate::lunar::EclipticMoonState> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    // Split the lunar-model TT observation into its integer day and fraction.
    let whole = jd.floor();
    let fraction = jd - whole;
    local_jpl_true_moon_state_parts(leg, whole, fraction)
}

/// True lunar state retaining an acquired scalar epoch's TT split.
/// Provider reads retain both input words; scalar TT/TDB observables alone
/// supply nutation/precession. No second normalization changes the fraction.
fn local_jpl_true_moon_state_parts(
    leg: &JplLeg<'_>,
    whole: f64,
    fraction: f64,
) -> Result<crate::lunar::EclipticMoonState> {
    for value in [whole, fraction] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let jd = whole + fraction;
    let earth_key = local_jpl_target(leg, EARTH)?;
    let moon_key = local_jpl_target(leg, MOON)?;
    let tdb_fraction =
        fraction + reference_frames::tdb_minus_tt_seconds(whole, fraction)? / 86400.0;
    let (relative, relative_vel) = leg.registry().eval_geometric_difference_au_two_part_jd(
        moon_key,
        earth_key,
        whole,
        tdb_fraction,
    )?;
    let tdb =
        whole + (fraction + reference_frames::tdb_minus_tt_seconds(whole, fraction)? / 86400.0);
    let eps = (reference_frames::mean_obliquity_arcsec(tdb)? / 3600.0).to_radians();
    let (psi, deps) = reference_frames::iau2000a_tenth_microarcsec(jd)?;
    let scale = reference_frames::FRAME_TENTH_USEC_2_RAD;
    let nut = reference_frames::nutation_matrix(eps, eps + deps * scale, psi * scale)?;
    let equinox = reference_frames::equinox_matrix(nut, reference_frames::precession_matrix(tdb)?)?;
    let frame = reference_frames::ecliptic_date_matrix(equinox, eps, deps * scale)?;
    let r = reference_frames::matvec3(frame, relative)?;
    let v = reference_frames::matvec3(frame, relative_vel)?;
    Ok(((r[0], r[1], r[2]), (v[0], v[1], v[2])))
}

/// Compute a geocentric lunar point from its native analytical or osculating
/// model, without looking up the point identifier as an SPK target.
///
/// Mean and interpolated models use TT and native nutation. Osculating models
/// acquire Moon/Earth geometry over their actual split-TDB record domains.
/// Mean-model validity uses the configured DE envelope; interpolated apsides
/// own their tapered formula's domain. Output is degrees, AU and per-day rates;
/// sidereal model clocks remain separate from the frame epoch. Model,
/// registration, nutation and coverage errors propagate.
fn local_jpl_lunar_position(
    leg: &JplLeg<'_>,
    tt: f64,
    original_ut: Option<f64>,
    body: i32,
    flags: i32,
    ayanamsha_at_ut: &dyn Fn(f64) -> Result<f64>,
) -> Result<(PositionResult, i32)> {
    if matches!(body, MEAN_NODE | MEAN_APOG) {
        local_jpl_lunar_coverage(leg, body, tt)?;
    } else if matches!(body, INTP_APOG | INTP_PERG) {
        local_jpl_main_window(leg)?;
    }
    let dpsi = |jd| crate::erfa::nut06a(jd).map(|(psi, _)| psi);
    let rate =
        |jd: f64, half: f64| Ok((dpsi(jd + half)? - dpsi(jd - half)?).to_degrees() / (2.0 * half));
    let mut output = match body {
        MEAN_NODE | MEAN_APOG => {
            crate::lunar::mean_lunar_state_with_providers(body, tt, flags, &dpsi, &rate)?
        }
        INTP_APOG | INTP_PERG => {
            let (a, b, c, d, e, f) =
                crate::lunar::interp_lunar_point_state(body, tt, flags, &dpsi, &rate)?;
            [a, b, c, d, e, f]
        }
        TRUE_NODE | OSCU_APOG => crate::lunar::true_lunar_state_with_providers(
            body,
            tt,
            flags,
            &|jd| local_jpl_true_moon_state(leg, jd),
            &dpsi,
            &rate,
        )?,
        _ => return Err(Error::UnknownBody { body_id: body }),
    };
    // SIDEREAL+J2000 subtracts ayanamsha only after the epoch transform.
    if flags & FLG_J2000 == 0 {
        local_jpl_sidereal_word(&mut output, tt, original_ut, flags, ayanamsha_at_ut)?;
    }
    local_jpl_lunar_frame(&mut output, tt, flags)?;
    if flags & FLG_J2000 != 0 {
        local_jpl_sidereal_word(&mut output, tt, original_ut, flags, ayanamsha_at_ut)?;
    }
    Ok((output, flags))
}

/// Project an ecliptic model state into the requested epoch/equator.
/// The six slots are longitude/latitude degrees, radius AU, angular rates
/// degrees/day and radial rate AU/day. `tt` is the entry-owned finite scalar
/// TT Julian date; this helper acquires no model or independent clock.
/// Undefined angular coordinates pass through unchanged. J2000 angular rates
/// differentiate the complete varying-epoch precession map at a half-day step;
/// radial components remain unchanged. Precession rotates the unit direction
/// with bias-free Vondrák matrices, retaining the right-transposed operand;
/// latitude is the arcsine of the clamped rotated z component, without a
/// further vector-norm division. The shared [`crate::hypothetical::precess_ecliptic`]
/// contract specifies the degree/radian bridge and platform matrix backend.
/// Equatorial rotation uses Vondrák-2011 of-date mean obliquity plus IAU
/// nutation except for NONUT or sidereal requests, in degrees. The angular
/// derivative folds a longitude difference outside ±180 degrees before
/// dividing by the one-day full stencil. This fixed stencil has no iterative
/// stopping rule. Required epoch/matrix/angle and equatorial conversion errors
/// propagate; no extra extrapolation or accuracy bound is introduced.
fn local_jpl_lunar_frame(output: &mut PositionResult, tt: f64, flags: i32) -> Result<()> {
    if !output[0].is_finite() || !output[1].is_finite() {
        return Ok(());
    }
    if flags & FLG_J2000 != 0 {
        let precess = |lon, lat, jd| {
            crate::hypothetical::precess_ecliptic(lon, lat, jd, crate::constants::J2000)
        };
        let (lon, lat) = precess(output[0], output[1], tt)?;
        if flags & FLG_SPEED != 0 {
            let (lo_m, la_m) = precess(
                output[0] - output[3] * 0.5,
                output[1] - output[4] * 0.5,
                tt - 0.5,
            )?;
            let (lo_p, la_p) = precess(
                output[0] + output[3] * 0.5,
                output[1] + output[4] * 0.5,
                tt + 0.5,
            )?;
            let mut delta = lo_p - lo_m;
            if delta > 180.0 {
                delta -= 360.0;
            } else if delta < -180.0 {
                delta += 360.0;
            }
            output[3] = delta;
            output[4] = la_p - la_m;
        }
        output[0] = lon;
        output[1] = lat;
    }
    if flags & FLG_EQUATORIAL != 0 {
        let eps = if flags & FLG_J2000 != 0 {
            84381.406 / 3600.0
        } else {
            vondrak_mean_obliquity_deg(tt)?
                + if flags & (FLG_NONUT | FLG_SIDEREAL) == 0 {
                    crate::erfa::nut06a(tt)?.1.to_degrees()
                } else {
                    0.0
                }
        };
        let (position, speed) = crate::utils::cotrans_sp(
            output[0], output[1], output[2], output[3], output[4], output[5], -eps,
        )?;
        *output = [
            position[0],
            position[1],
            position[2],
            speed[0],
            speed[1],
            speed[2],
        ];
    }
    Ok(())
}

/// Resolve a core body's centre first, then its registered barycentre.
/// Body keys are not assumed to be NAIF ids: both ordinary library-body
/// registration and `register_all_naif_targets` are accepted. Missing targets
/// remain typed registration failures. A centre present but out of coverage is
/// not silently exchanged for a different target.
pub(super) fn local_jpl_target(leg: &JplLeg<'_>, body: i32) -> Result<i32> {
    let candidates: &[i32] = match body {
        SUN => &[10],
        MOON => &[301],
        MERCURY => &[199, 1],
        VENUS => &[299, 2],
        MARS => &[499, 4],
        JUPITER => &[599, 5],
        SATURN => &[699, 6],
        URANUS => &[799, 7],
        NEPTUNE => &[899, 8],
        PLUTO => &[999, 9],
        EARTH => &[399],
        _ => {
            // Caller-registered target identifiers enter the shared native vector
            // calculation. Unregistered identifiers remain unknown-body errors.
            if leg.registry().body_info(body).is_some() {
                return Ok(body);
            }
            return Err(Error::UnknownBody { body_id: body });
        }
    };
    let registrations = leg.registry().registered_bodies();
    for target in candidates {
        if let Some(entry) = registrations.iter().find(|entry| entry.naif_id == *target) {
            return Ok(entry.body);
        }
    }
    Err(Error::SpkUnregisteredBody { body })
}

/// Evaluate a registered SSB state at split TT, returning AU and AU/day.
/// Convert TT to TDB at fraction magnitude and retain both epoch parts through
/// actual record selection, light-time sampling and polynomial evaluation.
pub(super) fn local_jpl_state(
    leg: &JplLeg<'_>,
    body: i32,
    whole: f64,
    fraction: f64,
) -> Result<([f64; 3], [f64; 3])> {
    let target = local_jpl_target(leg, body)?;
    let tdb_seconds = reference_frames::tdb_minus_tt_seconds(whole, fraction)?;
    let tdb_fraction = fraction + tdb_seconds / 86400.0;
    leg.registry()
        .eval_geometric_ssb_au_two_part_jd(target, whole, tdb_fraction)
}

/// Acquire an SSB state at an already validated split TT epoch.
/// The retained fraction receives the TDB correction before strict record
/// selection and polynomial evaluation. Units are ICRS AU and AU/day;
/// missing aliases, unsupported states and coverage errors propagate.
fn local_jpl_source_time_state(
    leg: &JplLeg<'_>,
    body: i32,
    whole: f64,
    fraction: f64,
) -> Result<([f64; 3], [f64; 3])> {
    let target = local_jpl_target(leg, body)?;
    let correction = crate::time::source_tdb_minus_tt_seconds(whole, fraction);
    leg.registry()
        .eval_geometric_ssb_au_two_part_jd(target, whole, fraction + correction / 86400.0)
}

/// Call-owned SPK segment metadata for analytical-model support reads.
/// The immutable image remains owned when its range failure is returned.
struct PublicModelDirectory<'a, 'registry> {
    leg: &'a JplLeg<'registry>,
    raw: std::cell::RefCell<Option<crate::spk::OwnedSpkRange>>,
}
impl PublicModelDirectory<'_, '_> {
    /// Evaluate an SSB target at split TDB, retaining any failed segment.
    fn tdb(&self, target: i32, whole: f64, fraction: f64) -> Result<([f64; 3], [f64; 3])> {
        #[cfg(test)]
        let _observed_query = super::context_serving::hypothetical_support_observations::spk_query(
            target, whole, fraction,
        );
        *self.raw.borrow_mut() = None;
        let mut captured = None;
        let result = self
            .leg
            .registry()
            .eval_geometric_ssb_au_two_part_jd_with_range(target, whole, fraction, &mut captured);
        *self.raw.borrow_mut() = captured;
        result
    }
    /// Evaluate an SSB body at split TT, applying its TDB correction in days.
    fn state(&self, body: i32, whole: f64, fraction: f64) -> Result<([f64; 3], [f64; 3])> {
        let target = local_jpl_target(self.leg, body)?;
        let correction = crate::time::source_tdb_minus_tt_seconds(whole, fraction);
        let tdb_fraction = fraction + correction / 86400.0;
        #[cfg(test)]
        super::context_serving::hypothetical_support_observations::tt_to_tdb(
            body,
            whole,
            fraction,
            tdb_fraction,
        );
        self.tdb(target, whole, tdb_fraction)
    }
    /// Retarded Earth-minus-Sun direction in AU, with the native light-time solver.
    fn earth(&self, tt: f64) -> Result<[f64; 3]> {
        retarded_earth_minus_sun_source_time_via(self.leg, tt, &|target, day, fraction| {
            self.tdb(target, day, fraction)
        })
    }
    /// Retain the same reception Sun for the predicted model's PPN origin.
    fn earth_observation(&self, tt: f64) -> Result<model_serving::RetardedEarthObservation> {
        retarded_earth_observation_source_time_via(self.leg, tt, &|target, day, fraction| {
            self.tdb(target, day, fraction)
        })
    }
    /// Report a real SPK coverage failure with its native cause and segment owner.
    /// Convergence and other numerical failures retain their own typed category.
    fn promote(&self, epoch: f64, body: i32, error: Error) -> Error {
        let raw = self.raw.borrow_mut().take();
        if !matches!(error.primary_error(), Error::SpkOutOfRange { .. }) {
            return error;
        }
        let cause = error.clone();
        let promoted = local_jpl_public_range_error(self.leg, epoch, body, error);
        let promoted = if matches!(promoted.primary_error(), Error::EphemerisRange { .. }) {
            promoted.with_cause(cause)
        } else {
            promoted
        };
        match raw {
            Some(raw) => promoted.with_backend_range(raw),
            None => promoted,
        }
    }
}

/// IERS 2010 geodetic site offset in GCRS AU and nominal AU/day velocity.
///
/// Longitude/latitude are east-positive degrees and altitude is metres.
/// IAU 2000A nutation, P03 precession and bundled TT-minus-UT1 define the
/// independent Earth-rotation clock. TT fractions are retained in UT1 rotation;
/// no configured model Delta-T or polar-motion offset is substituted.
pub(crate) fn local_jpl_topo(
    whole: f64,
    fraction: f64,
    site: (f64, f64, f64),
) -> Result<([f64; 3], [f64; 3])> {
    use reference_frames::*;
    let itrs = local_jpl_site_itrs(site)?;
    let (equinox, gast) = local_jpl_site_orientation(whole, fraction)?;
    let rotation = itrs_rotation_matrix(equinox, gast)?;
    let state = itrs_to_gcrs(rotation, itrs)?;
    #[cfg(test)]
    super::context_serving::hypothetical_support_observations::topocentric_site(
        whole, fraction, site, itrs, rotation, state,
    );
    Ok(state)
}

/// Fixed IERS-2010 geodetic site in terrestrial AU.
fn local_jpl_site_itrs(site: (f64, f64, f64)) -> Result<[f64; 3]> {
    use reference_frames::*;
    let (lon, lat, height) = site;
    check_finite_vec3([lon, lat, height])?;
    let phi = lat / 360.0 * FRAME_TAU;
    let lambda = lon / 360.0 * FRAME_TAU;
    let (s, c) = (scalar_sin(phi), scalar_cos(phi));
    let omf = (298.25642 - 1.0) / 298.25642;
    let omf2 = omf * omf;
    let n = 1.0 / (c * c + s * s * omf2).sqrt();
    let a = 6378136.6 / FRAME_AU_M;
    let h = height / FRAME_AU_M;
    let xy = (a * n + h) * c;
    Ok([
        xy * scalar_cos(lambda),
        xy * scalar_sin(lambda),
        (a * omf2 * n + h) * s,
    ])
}

type SiteOrientation = ([[f64; 3]; 3], f64);

/// The actual slow equinox matrix and GAST in hours at retained TT parts.
fn local_jpl_site_orientation(whole: f64, fraction: f64) -> Result<SiteOrientation> {
    use reference_frames::*;
    let tt = whole + fraction;
    let tdb = whole + (fraction + tdb_minus_tt_seconds(whole, fraction)? / 86400.0);
    let (psi, deps) = iau2000a_tenth_microarcsec(tt)?;
    let (psi, deps) = (psi * FRAME_TENTH_USEC_2_RAD, deps * FRAME_TENTH_USEC_2_RAD);
    let eps = mean_obliquity_arcsec(tdb)? * FRAME_ASEC2RAD;
    let nut = nutation_matrix(eps, eps + deps, psi)?;
    let equinox = equinox_matrix(nut, precession_matrix(tdb)?)?;
    let delta_t = crate::smh2016::bundled_deltat_seconds(tt)?;
    let gast = gast_hours_two_part(whole, fraction, delta_t, tdb)?;
    #[cfg(test)]
    super::context_serving::hypothetical_support_observations::topocentric_orientation(
        (whole, fraction),
        tdb,
        delta_t,
        (psi, deps),
        eps,
        equinox,
        gast,
    );
    Ok((equinox, gast))
}

/// One successful central SPK read and the offsets owned by its actual chain.
struct LocalJplReadDomain {
    body: i32,
    epoch: TwoPartEpoch,
    offsets: (f64, f64),
}

#[derive(Default)]
struct LocalJplReadDomains {
    ordinary: std::cell::RefCell<Vec<LocalJplReadDomain>>,
    retarded: std::cell::RefCell<Vec<LocalJplReadDomain>>,
}

impl LocalJplReadDomains {
    fn offsets(&self) -> (f64, f64) {
        self.ordinary
            .borrow()
            .iter()
            .chain(self.retarded.borrow().iter())
            .fold((f64::NEG_INFINITY, f64::INFINITY), |bounds, read| {
                (bounds.0.max(read.offsets.0), bounds.1.min(read.offsets.1))
            })
    }
}

/// Admission of one complete physical place on its central SPK records.
/// Light-time iteration counts can differ by one; additional converged target
/// reads must remain on the final target record. Ordinary observer/deflector
/// reads retain their separate acquisition roles and split TDB epochs.
struct LocalJplReadAdmission<'a> {
    domains: &'a LocalJplReadDomains,
    ordinary: std::cell::Cell<usize>,
    retarded: std::cell::Cell<usize>,
    last_retarded: std::cell::Cell<Option<(i32, TwoPartEpoch)>>,
    refused: std::cell::Cell<bool>,
}

impl<'a> LocalJplReadAdmission<'a> {
    fn new(domains: &'a LocalJplReadDomains) -> Self {
        Self {
            domains,
            ordinary: std::cell::Cell::new(0),
            retarded: std::cell::Cell::new(0),
            last_retarded: std::cell::Cell::new(None),
            refused: std::cell::Cell::new(false),
        }
    }

    fn complete(&self) -> bool {
        if self.refused.get() || self.ordinary.get() != self.domains.ordinary.borrow().len() {
            return false;
        }
        let domains = self.domains.retarded.borrow();
        match (domains.last(), self.last_retarded.get()) {
            (None, None) => true,
            (Some(domain), Some((body, epoch))) => {
                let offset = TwoPartEpoch {
                    hi: epoch.hi - domain.epoch.hi,
                    lo: epoch.lo - domain.epoch.lo,
                }
                .value();
                body == domain.body && domain.offsets.0 <= offset && offset <= domain.offsets.1
            }
            _ => false,
        }
    }

    fn check(&self, body: i32, epoch: TwoPartEpoch, retarded: bool) -> Result<()> {
        let (reads, cursor) = if retarded {
            (&self.domains.retarded, &self.retarded)
        } else {
            (&self.domains.ordinary, &self.ordinary)
        };
        let index = cursor.get();
        let reads = reads.borrow();
        let domain = reads
            .get(index)
            .or_else(|| retarded.then(|| reads.last()).flatten());
        if let Some(domain) = domain {
            let offset = TwoPartEpoch {
                hi: epoch.hi - domain.epoch.hi,
                lo: epoch.lo - domain.epoch.lo,
            }
            .value();
            if body == domain.body && domain.offsets.0 <= offset && offset <= domain.offsets.1 {
                cursor.set(index + 1);
                if retarded {
                    self.last_retarded.set(Some((body, epoch)));
                }
                return Ok(());
            }
            self.refused.set(true);
            // This is an internal owned-record admission refusal, caught only
            // by the derivative sampler. No substitute state is evaluated.
            return Err(Error::SpkOutOfRange {
                target: body,
                et: offset * 86400.0,
                start: domain.offsets.0 * 86400.0,
                end: domain.offsets.1 * 86400.0,
            });
        }
        self.refused.set(true);
        Err(Error::MathDomain {
            value: epoch.value(),
        })
    }

    /// A stable unavailable optional deflector consumes no state/role. If a
    /// centrally consumed read becomes unavailable *outside its owned record*,
    /// the whole side place is inadmissible even when the existing optional
    /// deflector policy skips that range error. In-domain failures are untouched.
    fn failed_read(&self, body: i32, epoch: TwoPartEpoch, retarded: bool) {
        let (reads, index) = if retarded {
            (&self.domains.retarded, self.retarded.get())
        } else {
            (&self.domains.ordinary, self.ordinary.get())
        };
        let reads = reads.borrow();
        if let Some(domain) = reads
            .get(index)
            .or_else(|| retarded.then(|| reads.last()).flatten())
        {
            let offset = TwoPartEpoch {
                hi: epoch.hi - domain.epoch.hi,
                lo: epoch.lo - domain.epoch.lo,
            }
            .value();
            if body == domain.body && !(domain.offsets.0 <= offset && offset <= domain.offsets.1) {
                self.refused.set(true);
            }
        }
    }
}

/// One physical SPK place at split TT, with optional owned-record admission.
/// Target identity precedes observer reads. Every state callback retains the
/// supplied TT/TDB fractions through the physical observer/light-path/frame
/// model; registration, domain and numerical errors propagate unchanged.
/// `site_state_override` supplies the actual site position/nominal velocity
/// on its independent TT clock while retaining live ephemeris light paths.
#[allow(clippy::too_many_arguments)]
fn local_jpl_place(
    leg: &JplLeg<'_>,
    whole: f64,
    fraction: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    capture: Option<&LocalJplReadDomains>,
    admission: Option<&LocalJplReadAdmission<'_>>,
    site_state_override: Option<([f64; 3], [f64; 3])>,
) -> Result<NativeVectorPlace> {
    local_jpl_target(leg, body)?;
    let state_at_tdb = |id, w, f, retarded| {
        let target = local_jpl_target(leg, id)?;
        let value = match leg
            .registry()
            .eval_geometric_ssb_au_two_part_jd(target, w, f)
        {
            Ok(value) => value,
            Err(error) => {
                if let Some(admission) = admission {
                    admission.failed_read(target, TwoPartEpoch { hi: w, lo: f }, retarded);
                }
                return Err(error);
            }
        };
        if let Some(admission) = admission {
            admission.check(target, TwoPartEpoch { hi: w, lo: f }, retarded)?;
        }
        if let Some(capture) = capture {
            let read = LocalJplReadDomain {
                body: target,
                epoch: TwoPartEpoch { hi: w, lo: f },
                offsets: leg.registry().geometric_ssb_record_offsets(target, w, f)?,
            };
            (if retarded {
                &capture.retarded
            } else {
                &capture.ordinary
            })
            .borrow_mut()
            .push(read);
        }
        Ok(value)
    };
    // Retain both acquired TT parts, even when the fraction is zero.
    let retarded_at = |id, w, f, delay| {
        // Subtract delay directly from the fractional TDB reception epoch.
        let tdb_fraction = f + reference_frames::tdb_minus_tt_seconds(w, f)? / 86400.0;
        state_at_tdb(id, w, tdb_fraction - delay, true)
    };
    let site_velocity_tt =
        if capture.is_some() && flags & FLG_TOPOCTR != 0 && flags & (FLG_HELCTR | FLG_BARYCTR) == 0
        {
            site.map(|site| local_jpl_topo_velocity_tt(whole, fraction, site))
                .transpose()?
        } else {
            None
        };
    vector_place_from_state_routes_with_policy(
        &|id, w, f| {
            let correction = reference_frames::tdb_minus_tt_seconds(w, f)? / 86400.0;
            state_at_tdb(id, w, f + correction, false)
        },
        Some(&retarded_at),
        whole,
        fraction,
        body,
        flags,
        site,
        false,
        Some(crate::time::tdb_rate(whole + fraction)?),
        site_velocity_tt,
        site_state_override,
    )
}

/// An output direction with the owned physical astrometric range derivative.
struct NativeVectorPlace {
    /// Longitude/latitude in degrees and distance in AU.
    coordinates: [f64; 3],
    /// Actual apparent ray in the output frame, AU, before spherical shaping.
    output_ray: [f64; 3],
    /// Final apparent ray in ICRS, AU, before the output-frame rotation.
    apparent_ray: [f64; 3],
    /// Astrometric light-path velocity in ICRS AU per observation TT day.
    ray_velocity_tt: [f64; 3],
    /// Small actual PPN/Lorentz increment for the native optional-support map.
    /// Required-support legacy callers do not consume this derivative channel.
    correction_increment: Option<[f64; 3]>,
    /// Physical AU/AU-per-TT-day state, before direction corrections/frames.
    range: Option<super::planet_range::RangeState>,
}

/// Source early sentinel for geocentric Earth, including topocentric Earth.
/// Sun heliocentric is deliberately excluded: it must acquire its central
/// Sun state before the resulting zero geometry can bypass derivative reads.
fn vector_self_observation(body: i32, flags: i32) -> bool {
    body == EARTH && flags & (FLG_BARYCTR | FLG_HELCTR) == 0
}

/// Immutable SSB-state callback: body id and split TT days to AU/AU-day.
/// Implementations enforce coverage and finite states at their own boundary.
pub(super) type VectorStateAt<'a> = dyn Fn(i32, f64, f64) -> Result<([f64; 3], [f64; 3])> + 'a;

/// Stored SSB ICRS state at split TT, with an explicit retarded-target role.
pub(super) type RoleVectorStateAt<'a> =
    dyn Fn(i32, f64, f64, bool) -> Result<([f64; 3], [f64; 3])> + 'a;

/// Retarded state: body id, observation TT parts and delay in days.
/// A registered-SPK implementation subtracts the delay from observation TDB
/// for every center; the ordinary callback retains its provider's TT route.
type RetardedVectorStateAt<'a> = dyn Fn(i32, f64, f64, f64) -> Result<([f64; 3], [f64; 3])> + 'a;

/// Native physical-observation vector, distinct from public calc shaping.
/// The central/closest deflector Times preserve their TT parts, while target
/// retardation subtracts delay from observation TDB through `retarded_at`.
/// Strict Sun/Jupiter/Saturn aliases use the public center/barycenter selection.
/// Missing supports and raw provider failures propagate; no scalar envelope
/// or optional-deflector policy is imposed. Every supplied state is finite.
/// Geometric input returns before any retarded/deflector acquisition.
///
/// Observer position is SSB ICRS AU and velocity is AU/day at the observation
/// Time. The optional surface vector is geocenter-to-observer ICRS AU. The
/// retarded callback receives that original split TT Time and delay in days,
/// and subtracts the delay from observation TDB. Light time is iterated for at
/// most ten steps with a 1e-12-day change threshold; exhaustion returns
/// `Error::ConvergenceFailed`. Other acquisition errors propagate in read order.
/// `apparent=false` returns the astrometric observer-relative AU ray; otherwise
/// strict deflection and aberration follow. Raw states and caller parts are
/// finite-checked, while derived IEEE NaNs/infinities remain output values.
#[allow(clippy::too_many_arguments)]
fn native_pheno_observe(
    state_at: &dyn Fn(i32, TwoPartEpoch) -> Result<[f64; 3]>,
    retarded_at: &dyn Fn(TwoPartEpoch, f64) -> Result<[f64; 3]>,
    epoch: TwoPartEpoch,
    body: i32,
    observer: [f64; 3],
    velocity: [f64; 3],
    surface: Option<[f64; 3]>,
    geometric: bool,
    apparent: bool,
) -> Result<[f64; 3]> {
    check_finite_vec3([epoch.hi, epoch.lo, epoch.value()])?;
    check_finite_vec3(observer)?;
    check_finite_vec3(velocity)?;
    if let Some(surface) = surface {
        check_finite_vec3(surface)?;
    }
    let target = state_at(body, epoch)?;
    check_finite_vec3(target)?;
    let vector = std::array::from_fn(|i| target[i] - observer[i]);
    if geometric {
        return Ok(vector);
    }
    let (vector, light_time) = native_astrometric_ray(retarded_at, epoch, vector, observer)?;
    if !apparent {
        return Ok(vector);
    }
    native_apparent_ray(
        state_at, epoch, vector, light_time, observer, velocity, surface,
    )
}

/// Retard an already acquired observer-relative ray with the accepted native
/// ten-step light-time loop. The observer is held at observation time; every
/// target callback receives the original split TT and a delay in days to be
/// subtracted from its TDB fraction. Acquired target vectors must be finite.
/// Returns the ICRS AU ray and its final astrometric light time in days.
/// This helper performs no flag/provider selection or apparent reduction.
/// Exhaustion and target acquisition failures preserve the original errors.
pub(super) fn native_astrometric_ray(
    retarded_at: &dyn Fn(TwoPartEpoch, f64) -> Result<[f64; 3]>,
    epoch: TwoPartEpoch,
    mut vector: [f64; 3],
    observer: [f64; 3],
) -> Result<([f64; 3], f64)> {
    let mut previous = 0.0;
    let mut light_time = 0.0;
    let mut converged = false;
    for _ in 0..10 {
        light_time = vec3_dist(vector) / FRAME_LIGHT_SPEED_AU_DAY;
        if (light_time - previous).abs() < 1e-12 {
            converged = true;
            break;
        }
        let target = retarded_at(epoch, light_time)?;
        check_finite_vec3(target)?;
        vector = std::array::from_fn(|i| target[i] - observer[i]);
        previous = light_time;
    }
    if !converged {
        return Err(Error::ConvergenceFailed {
            algorithm: "spk_light_time",
            iterations: 10,
            max_iterations: 10,
            tolerance: 1e-12,
            last_value: light_time,
        });
    }
    Ok((vector, light_time))
}

/// Apply native apparent-place corrections to an acquired astrometric ray.
///
/// Ray and observer position are ICRS AU, observer velocity is AU/day, and
/// epoch is split TT. The one-way delay is days. Each planetary deflector is
/// sampled at the projected closest-approach delay, bounded by reception and
/// emission. Sun/Jupiter/Saturn PPN terms precede optional surface-Earth
/// deflection and Lorentz aberration. Missing support and numerical errors
/// propagate; the caller validates acquired raw states and epoch parts.
#[allow(clippy::too_many_arguments)]
pub(super) fn native_apparent_ray(
    state_at: &dyn Fn(i32, TwoPartEpoch) -> Result<[f64; 3]>,
    epoch: TwoPartEpoch,
    vector: [f64; 3],
    light_time: f64,
    observer: [f64; 3],
    velocity: [f64; 3],
    surface: Option<[f64; 3]>,
) -> Result<[f64; 3]> {
    native_apparent_ray_with_deflectors(
        state_at, epoch, vector, light_time, observer, velocity, surface, true,
    )
}

/// Native apparent ray with an explicit planetary deflector-list selection.
///
/// Inputs and outputs follow [`native_apparent_ray`]. `include_deflectors`
/// selects the default Sun/Jupiter/Saturn list or the empty list; the optional
/// Earth term and aberration still run for an empty list. The arithmetic tail
/// and read order are identical for the default list. This switch supplies
/// the public star NOGDEFL and differential NOABERR branches; it does not
/// relax finite acquired-state checks or silently skip absent/ranged supports.
#[allow(clippy::too_many_arguments)]
pub(super) fn native_apparent_ray_with_deflectors(
    state_at: &dyn Fn(i32, TwoPartEpoch) -> Result<[f64; 3]>,
    epoch: TwoPartEpoch,
    mut vector: [f64; 3],
    light_time: f64,
    observer: [f64; 3],
    velocity: [f64; 3],
    surface: Option<[f64; 3]>,
    include_deflectors: bool,
) -> Result<[f64; 3]> {
    // Astrometric.apparent freezes its initial ray length for all three
    // closest-approach clamps; each later ray includes earlier deflection.
    let tlt = vec3_dist(vector) / FRAME_LIGHT_SPEED_AU_DAY;
    if include_deflectors {
        for (id, reciprocal_mass) in DEFLECTORS {
            let position = state_at(id, epoch)?;
            check_finite_vec3(position)?;
            let gpv: [f64; 3] = std::array::from_fn(|i| position[i] - observer[i]);
            let length = vec3_dist(vector);
            let unit: [f64; 3] = std::array::from_fn(|i| vector[i] / (length + f64::MIN_POSITIVE));
            let dlt =
                (unit[0] * gpv[0] + unit[1] * gpv[1] + unit[2] * gpv[2]) / FRAME_LIGHT_SPEED_AU_DAY;
            let delay = dlt.clamp(0.0, tlt);
            // Subtract the projected delay from the retained observation TT parts.
            let whole_delay = delay.floor();
            let fraction_delay = delay - whole_delay;
            let closest = TwoPartEpoch {
                hi: epoch.hi - whole_delay,
                lo: epoch.lo - fraction_delay,
            };
            let position = state_at(id, closest)?;
            check_finite_vec3(position)?;
            let pe = std::array::from_fn(|i| observer[i] - position[i]);
            let deflection = native_pheno_deflection(vector, pe, reciprocal_mass);
            for i in 0..3 {
                vector[i] += deflection[i];
            }
        }
    }
    if let Some(surface) = surface {
        let deflection = native_pheno_deflection(vector, surface, EARTH_DEFLECTION_RMASS);
        let aprad = ((IERS_EARTH_RADIUS_M / AU_IN_M) / vec3_dist(surface))
            .min(1.0)
            .asin();
        let cosine = ((vector[0] * surface[0] + vector[1] * surface[1] + vector[2] * surface[2])
            / (vec3_dist(vector) * vec3_dist(surface)))
        .clamp(-1.0, 1.0);
        let nadir = (std::f64::consts::PI - cosine.acos()) / aprad;
        let include = if nadir >= EARTH_DEFLECTION_LIMB_GATE {
            1.0
        } else {
            0.0
        };
        for i in 0..3 {
            vector[i] += deflection[i] * include;
        }
    }
    Ok(native_pheno_aberration(vector, velocity, light_time))
}

/// Apply the Lorentz aberration to an already acquired apparent ray.
/// Position is in AU, velocity in AU/day and light time in days. The original
/// astrometric delay fixes the ray scale; it is not recomputed after deflection.
/// This raw primitive preserves IEEE results for singular/non-finite inputs.
/// Checked native callers retain their acquired-state validation before entry.
pub fn native_pheno_aberration(
    mut vector: [f64; 3],
    velocity: [f64; 3],
    light_time: f64,
) -> [f64; 3] {
    let p1mag = light_time * FRAME_LIGHT_SPEED_AU_DAY;
    let vemag = vec3_dist(velocity);
    let beta = vemag / FRAME_LIGHT_SPEED_AU_DAY;
    let dot = vector[0] * velocity[0] + vector[1] * velocity[1] + vector[2] * velocity[2];
    let cosd = dot / (p1mag * vemag + f64::MIN_POSITIVE);
    let gammai = (1.0 - beta * beta).sqrt();
    let p = beta * cosd;
    let q = (1.0 + p / (1.0 + gammai)) * light_time;
    let r = 1.0 + p;
    for i in 0..3 {
        vector[i] *= gammai;
        vector[i] += q * velocity[i];
        vector[i] /= r;
    }
    vector
}

/// Finite-source PPN deflection of an ICRS ray by one body.
///
/// Ray and observer-minus-deflector vector are AU; `reciprocal_mass` is the
/// Sun-to-deflector mass ratio. The gamma-one scale is
/// `2*GM_sun/(c²*|pe|*AU)`. The alignment gate suppresses collinear geometry.
/// This low-level kernel performs no acquisition or finite validation;
/// singular derived geometry can produce non-finite values. Its caller owns
/// state validation, body support and the surface-Earth applicability gate.
pub fn native_pheno_deflection(ray: [f64; 3], pe: [f64; 3], reciprocal_mass: f64) -> [f64; 3] {
    let pq: [f64; 3] = std::array::from_fn(|i| ray[i] + pe[i]);
    let pmag = vec3_dist(ray);
    let qmag = vec3_dist(pq);
    let emag = vec3_dist(pe);
    let phat: [f64; 3] = std::array::from_fn(|i| ray[i] / if pmag != 0.0 { pmag } else { 1.0 });
    let qhat: [f64; 3] = std::array::from_fn(|i| pq[i] / if qmag != 0.0 { qmag } else { 1.0 });
    let ehat: [f64; 3] = std::array::from_fn(|i| pe[i] / if emag != 0.0 { emag } else { 1.0 });
    let pdotq = phat[0] * qhat[0] + phat[1] * qhat[1] + phat[2] * qhat[2];
    let qdote = qhat[0] * ehat[0] + qhat[1] * ehat[1] + qhat[2] * ehat[2];
    let edotp = ehat[0] * phat[0] + ehat[1] * phat[1] + ehat[2] * phat[2];
    let flag = if edotp.abs() <= DEFLECTION_EDOTP_LIMIT {
        1.0
    } else {
        0.0
    };
    let fac1 = 2.0 * GRAV_GS
        / (SPEED_OF_LIGHT_M_S * SPEED_OF_LIGHT_M_S * emag * AU_IN_M * reciprocal_mass);
    let fac2 = 1.0 + qdote;
    std::array::from_fn(|i| flag * fac1 * (pdotq * ehat[i] - edotp * qhat[i]) / fac2 * pmag)
}

/// Reduce finite SSB ICRS AU/AU-day states at split TT to a scalar place.
///
/// The callback retains epoch parts and enforces actual body/record domains.
/// Only the target is retarded; geocentric light time uses a ten-iteration,
/// 1e-12-day update criterion, and centered work uses three iterations.
/// Planetary deflection precedes aberration. Surface observers use IERS 2010
/// geometry and the bundled rotation clock. IAU nutation and long-term
/// precession supply the frame. Output is degrees/degrees/AU; sidereal
/// longitude correction remains the caller's separate model-clock operation.
/// Origin sentinels require no side-epoch reads. Errors propagate unchanged.
#[cfg(test)]
pub(super) fn vector_place_from_states(
    state_at: &VectorStateAt<'_>,
    whole: f64,
    fraction: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
) -> Result<[f64; 3]> {
    vector_place_from_state_routes(state_at, None, whole, fraction, body, flags, site)
}

/// Scalar vector place with explicit central and geocentric retarded channels.
/// Inputs/units/coverage obligations follow [`vector_place_from_states`].
/// Optional retardation serves every apparent registered-SPK target read.
/// TRUEPOS never invokes it; providers without that callback retain TT sampling.
#[allow(clippy::too_many_arguments)]
#[cfg(test)]
fn vector_place_from_state_routes(
    state_at: &VectorStateAt<'_>,
    retarded_at: Option<&RetardedVectorStateAt<'_>>,
    whole: f64,
    fraction: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
) -> Result<[f64; 3]> {
    vector_place_from_state_routes_with_policy(
        state_at,
        retarded_at,
        whole,
        fraction,
        body,
        flags,
        site,
        false,
        None,
        None,
        None,
    )
    .map(|place| place.coordinates)
}

/// Stored scalar place with explicit ordinary and retarded acquisition roles.
///
/// The callback returns SSB ICRS AU/AU-day at split TT; `true` identifies only
/// a target light-time read, after the caller's selected retardation clock.
/// Observer and closest-approach deflectors use `false`. Separating these roles
/// permits exact record admission even when light-time iteration counts vary.
/// `strict_support` retains required/optional apparent-support policy and its
/// TDB-delay/TT or direct TT-retardation route, respectively. Provider errors
/// retain that policy without an alternate numerical model.
/// Output is longitude/latitude in degrees and distance in AU. Site input is
/// east longitude/north latitude in degrees and height in metres; all epochs
/// are TT days. No derivative or polynomial extrapolation is performed here.
#[allow(clippy::too_many_arguments)]
pub(super) fn stored_vector_place_with_read_roles(
    state_at: &RoleVectorStateAt<'_>,
    whole: f64,
    fraction: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    strict_support: bool,
) -> Result<[f64; 3]> {
    let retarded = |id, w, f, delay| {
        if !strict_support {
            return state_at(id, w, f - delay, true);
        }
        let tdb_fraction = f + reference_frames::tdb_minus_tt_seconds(w, f)? / 86400.0;
        let retarded_fraction = tdb_fraction - delay;
        let tt_fraction = retarded_fraction
            - reference_frames::tdb_minus_tt_seconds(w, retarded_fraction)? / 86400.0;
        state_at(id, w, tt_fraction, true)
    };
    vector_place_from_state_routes_with_policy(
        &|id, w, f| state_at(id, w, f, false),
        Some(&retarded),
        whole,
        fraction,
        body,
        flags,
        site,
        strict_support,
        None,
        None,
        None,
    )
    .map(|place| place.coordinates)
}

/// Shared geometry; the policy changes only consumed apparent support roles.
#[allow(clippy::too_many_arguments)]
fn vector_place_from_state_routes_with_policy(
    state_at: &VectorStateAt<'_>,
    retarded_at: Option<&RetardedVectorStateAt<'_>>,
    whole: f64,
    fraction: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    strict_support: bool,
    range_clock_rate: Option<f64>,
    site_velocity_tt: Option<[f64; 3]>,
    site_state_override: Option<([f64; 3], [f64; 3])>,
) -> Result<NativeVectorPlace> {
    let tt = whole + fraction;
    check_finite_vec3([whole, fraction, tt])?;
    if vector_self_observation(body, flags) {
        return Ok(NativeVectorPlace {
            coordinates: [0.0; 3],
            output_ray: [0.0; 3],
            apparent_ray: [0.0; 3],
            ray_velocity_tt: [0.0; 3],
            correction_increment: Some([0.0; 3]),
            range: range_clock_rate
                .map(|_| super::planet_range::RangeState::from_cartesian([0.0; 3], [0.0; 3])),
        });
    }
    let bary = flags & FLG_BARYCTR != 0;
    let helio = !bary && flags & FLG_HELCTR != 0;
    let target_first = strict_support && (flags & FLG_TRUEPOS != 0 || bary || helio);
    let early_target = if target_first {
        Some(state_at(body, whole, fraction)?)
    } else {
        None
    };
    let (mut observer, mut velocity) = if bary {
        ([0.0; 3], [0.0; 3])
    } else {
        state_at(if helio { SUN } else { EARTH }, whole, fraction)?
    };
    check_finite_vec3(observer)?;
    check_finite_vec3(velocity)?;
    let observation_rate = range_clock_rate.unwrap_or(1.0);
    let mut observer_velocity_tt = velocity.map(|v| v * observation_rate);
    let mut angular_site_adjustment = [0.0; 3];
    let mut surface = None;
    if !bary && !helio && flags & FLG_TOPOCTR != 0 {
        let site = site.ok_or(Error::ConfigurationMissing {
            missing_config: "observer_location",
            suggestion: "Call set_topo(lon, lat, alt) first",
        })?;
        let (pos, vel) = match site_state_override {
            Some(state) => state,
            None => local_jpl_topo(whole, fraction, site)?,
        };
        check_finite_vec3(pos)?;
        check_finite_vec3(vel)?;
        surface = Some(pos);
        for i in 0..3 {
            observer[i] += pos[i];
            velocity[i] += vel[i];
            // The SPK velocity is per TDB day, while the rotating site's
            // nominal GCRS velocity is already per TT day.
            observer_velocity_tt[i] += vel[i];
            angular_site_adjustment[i] = site_velocity_tt.unwrap_or(vel)[i] - vel[i];
        }
    }
    let (target, target_velocity) = match early_target {
        Some(state) => state,
        None => state_at(body, whole, fraction)?,
    };
    check_finite_vec3(target)?;
    let mut geo = std::array::from_fn(|i| target[i] - observer[i]);
    let ray_velocity_tt = std::cell::Cell::new([0.0; 3]);
    let angular_range = std::cell::Cell::new(None::<super::planet_range::RangeState>);
    let range_at = |ray: [f64; 3],
                    target_velocity: [f64; 3],
                    emission_rate: f64,
                    angular_emission_rate: f64,
                    observer_velocity_tt: [f64; 3]| {
        let derivative =
            std::array::from_fn(|i| target_velocity[i] * emission_rate - observer_velocity_tt[i]);
        let angular_derivative = std::array::from_fn(|i| {
            target_velocity[i] * angular_emission_rate
                - (observer_velocity_tt[i] + angular_site_adjustment[i])
        });
        ray_velocity_tt.set(angular_derivative);
        angular_range.set(Some(super::planet_range::RangeState::from_cartesian(
            ray,
            angular_derivative,
        )));
        range_clock_rate.map(|_| super::planet_range::RangeState::from_cartesian(ray, derivative))
    };
    let mut range = range_at(
        geo,
        target_velocity,
        observation_rate,
        observation_rate,
        observer_velocity_tt,
    );
    if vec3_dist(geo) == 0.0 {
        return Ok(NativeVectorPlace {
            coordinates: [0.0; 3],
            output_ray: [0.0; 3],
            apparent_ray: [0.0; 3],
            ray_velocity_tt: ray_velocity_tt.get(),
            correction_increment: Some([0.0; 3]),
            range,
        });
    }
    let mut lt = 0.0;
    let mut correction_increment = Some([0.0; 3]);
    if flags & FLG_TRUEPOS == 0 {
        if strict_support && helio {
            (observer, velocity) = state_at(SUN, whole, fraction)?;
            observer_velocity_tt = velocity.map(|v| v * observation_rate);
        }
        let iterations = if bary || helio { 3 } else { 10 };
        let mut converged = bary || helio;
        for _ in 0..iterations {
            let c = if bary || helio {
                C_LIGHT_AU_DAY
            } else {
                FRAME_LIGHT_SPEED_AU_DAY
            };
            let next_lt = vec3_dist(geo) / c;
            if !bary && !helio && (next_lt - lt).abs() < 1e-12 {
                converged = true;
                break;
            }
            lt = next_lt;
            let emission_rate = range.map_or(observation_rate, |state| {
                state.next_emission_rate(observation_rate, c)
            });
            let angular_emission_rate = angular_range.get().map_or(observation_rate, |state| {
                state.next_emission_rate(observation_rate, c)
            });
            let (retarded, retarded_velocity) = match retarded_at {
                Some(at) => at(body, whole, fraction, lt)?,
                None if strict_support => {
                    let tdb_fraction = fraction
                        + reference_frames::tdb_minus_tt_seconds(whole, fraction)? / 86400.0;
                    let retarded_fraction = tdb_fraction - lt;
                    let tt_fraction = retarded_fraction
                        - reference_frames::tdb_minus_tt_seconds(whole, retarded_fraction)?
                            / 86400.0;
                    state_at(body, whole, tt_fraction)?
                }
                None => state_at(body, whole, fraction - lt)?,
            };
            check_finite_vec3(retarded)?;
            geo = std::array::from_fn(|i| retarded[i] - observer[i]);
            range = range_at(
                geo,
                retarded_velocity,
                emission_rate,
                angular_emission_rate,
                observer_velocity_tt,
            );
        }
        if !converged {
            return Err(Error::ConvergenceFailed {
                algorithm: "spk_light_time",
                iterations: 10,
                max_iterations: 10,
                tolerance: 1e-12,
                last_value: lt,
            });
        }
        if !bary && !helio && strict_support {
            correction_increment = None;
            // An empty planetary-deflector list still permits Earth surface deflection and aberration.
            if flags & FLG_NOGDEFL != 0 && flags & FLG_NOABERR == 0 {
                geo = native_apparent_ray_with_deflectors(
                    &|id, time| state_at(id, time.hi, time.lo).map(|state| state.0),
                    TwoPartEpoch {
                        hi: whole,
                        lo: fraction,
                    },
                    geo,
                    lt,
                    observer,
                    velocity,
                    surface,
                    false,
                )?;
            } else {
                if flags & FLG_NOGDEFL == 0 {
                    geo = apply_gravitational_deflection_with_lookup(
                        geo,
                        observer,
                        tt,
                        lt,
                        |id, epoch| match state_at(id, epoch, 0.0) {
                            Ok(state) => Ok(Some(state.0)),
                            Err(error) => {
                                let skipped = matches!(
                                    error.primary_error(),
                                    Error::LebBodyNotFound { .. } | Error::LebOutOfRange { .. }
                                );
                                if skipped { Ok(None) } else { Err(error) }
                            }
                        },
                    )?;
                }
                if flags & FLG_NOABERR == 0 {
                    geo = apply_aberration_with_c(
                        geo,
                        velocity,
                        vec3_dist(geo) / FRAME_LIGHT_SPEED_AU_DAY,
                        FRAME_LIGHT_SPEED_AU_DAY,
                    )?;
                }
            }
        } else if !bary && !helio {
            if flags & FLG_NOGDEFL == 0 {
                let corrected = super::frames::apply_gravitational_deflection_with_correction(
                    geo,
                    observer,
                    tt,
                    lt,
                    |id, epoch| state_at(id, epoch, 0.0).map(|state| state.0),
                )?;
                geo = corrected.position;
                correction_increment = Some(corrected.increment);
            }
            if flags & FLG_NOABERR == 0 {
                let corrected = super::frames::apply_aberration_with_correction(
                    geo,
                    velocity,
                    vec3_dist(geo) / FRAME_LIGHT_SPEED_AU_DAY,
                    FRAME_LIGHT_SPEED_AU_DAY,
                )?;
                geo = corrected.position;
                let increment = correction_increment
                    .as_mut()
                    .expect("optional-support increment is retained");
                for (total, correction) in increment.iter_mut().zip(corrected.increment) {
                    *total += correction;
                }
            }
        }
    }
    let (psi, deps) = reference_frames::iau2000a_tenth_microarcsec(tt)?;
    let scale = reference_frames::FRAME_TENTH_USEC_2_RAD;
    let frames = frame_matrices_for_date(tt, psi * scale, deps * scale)?;
    let output_ray = if flags & FLG_NONUT != 0 && flags & (FLG_EQUATORIAL | FLG_J2000) == 0 {
        let matrix = vondrak_precession_matrix(tt, flags & FLG_ICRS == 0)?;
        let eq = mat3_vec3(matrix, geo)?;
        rotate_equatorial_to_ecliptic(eq, vondrak_mean_obliquity_rad(tt)?)?
    } else {
        frame_transform(geo, tt, flags, true, frames)?
    };
    let mut coordinates = cartesian_to_spherical_deg(output_ray)?;
    if let Some(range) = range {
        coordinates[2] = range.distance;
    }
    Ok(NativeVectorPlace {
        coordinates,
        output_ray,
        apparent_ray: geo,
        ray_velocity_tt: ray_velocity_tt.get(),
        correction_increment,
        range,
    })
}

/// Validate a finite scalar public request against its selected main-image span.
/// Descriptor endpoints form the public closed TT/UT1 envelope before any TDB
/// conversion. Body/registration and range errors retain their request identity;
/// derived split TDB reads enforce their separate directory domains later.
fn validate_local_jpl_public_epoch(
    leg: &JplLeg<'_>,
    epoch: f64,
    body: i32,
    public_entry: bool,
) -> Result<()> {
    if !epoch.is_finite() {
        return Err(Error::NonFiniteArgument { value: epoch });
    }
    let target = local_jpl_target(leg, body)?;
    let target = leg.registry().primary_body_for_coverage(target)?;
    if let Some(kernel) = leg
        .registry()
        .body_info(target)
        .and_then(|info| leg.registry().kernel(&info.kernel))
    {
        let (start, end) = kernel.index().all_segments().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |span, segment| {
                (
                    span.0.min(segment.descriptor.start_et),
                    span.1.max(segment.descriptor.end_et),
                )
            },
        );
        let start_jd = FRAME_J2000 + start / 86400.0;
        let end_jd = FRAME_J2000 + end / 86400.0;
        if public_entry && start_jd == 0.0 && end_jd == 0.0 {
            return Ok(());
        }
        if epoch < start_jd || epoch > end_jd {
            return Err(Error::EphemerisRange {
                requested_jd: epoch,
                start_jd,
                end_jd,
                body: Some(body),
            });
        }
    }
    Ok(())
}

/// Validate a scalar TT entry without public UT1 metadata.
/// This adapter retains the direct TT and legacy provider request contract.
fn validate_local_jpl_scalar_tt(leg: &JplLeg<'_>, tt: f64, body: i32) -> Result<()> {
    validate_local_jpl_public_epoch(leg, tt, body, false)
}

/// Select the central split TT epoch for a native public calculation.
///
/// A matching UT1 request retains its original UT1 JD plus the twice-refined
/// bundled TT-minus-UT1 fraction in days. A TT request uses its integer day
/// and remainder. Configured model Delta-T is separate. SPEED displacements
/// act on this retained fractional part, not a rounded full JD. Invalid clock
/// results propagate; this helper acquires no provider, site or frame.
pub(super) fn native_position_central_time(
    tt: f64,
    original_ut: Option<f64>,
) -> Result<(f64, f64)> {
    if let Some(ut) = original_ut {
        let pair = crate::smh2016::timescale_parts_from_ut1(ut)?;
        if (pair.0 + pair.1).to_bits() == tt.to_bits() {
            return Ok(pair);
        }
    }
    let whole = tt.floor();
    Ok((whole, tt - whole))
}

/// Public six-component SPK core with a domain-limited TT place derivative.
/// Angular rates differentiate the same retarded ray, observer and output frame
/// as the central position; the central light-path Jacobian supplies range rate.
/// Non-SPEED tails are exact zeros. Ecliptic sidereal reduction is applied by
/// the context-bound caller after this tropical derivative.
/// Equatorial SIDEREAL selects the mean-equator convention without ayanamsha
/// subtraction. Public scalar TT and derived split TDB domains stay separate.
/// `public_epoch` is the original scalar TT or UT1 JD; only the closed main
/// file guard consumes it. Split-TDB geometry, retardation and speed samples
/// consume the selected observation epoch parts. `original_ut` retains the central
/// public UT1 metadata: a matching central sample uses that epoch's full TT
/// pair; fresh TT entries use an integer-day/fraction TT split. SPEED samples retain that
/// same pair and displace its small TT fraction. Registration, scalar-range
/// and provider failures propagate.
#[allow(clippy::too_many_arguments)]
fn local_jpl_position(
    leg: &JplLeg<'_>,
    tt: f64,
    public_epoch: f64,
    original_ut: Option<f64>,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
) -> Result<(PositionResult, i32)> {
    local_jpl_position_with_public_entry(
        leg,
        tt,
        public_epoch,
        original_ut,
        body,
        flags,
        site,
        false,
    )
}

/// Preserve the original geometry while retaining genuine public range metadata.
/// The optional constructor is consumed solely at the existing scalar guard.
#[allow(clippy::too_many_arguments)]
fn local_jpl_position_with_public_entry(
    leg: &JplLeg<'_>,
    tt: f64,
    public_epoch: f64,
    original_ut: Option<f64>,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    public_entry: bool,
) -> Result<(PositionResult, i32)> {
    validate_local_jpl_public_epoch(leg, public_epoch, body, public_entry)?;
    if vector_self_observation(body, flags) {
        return native_source_result(Ok(([0.0; 6], flags)), "Analytical");
    }
    local_jpl_target(leg, body)?;
    let (central_whole, central_fraction) = native_position_central_time(tt, original_ut)?;
    let domains = LocalJplReadDomains::default();
    let central = local_jpl_place(
        leg,
        central_whole,
        central_fraction,
        body,
        flags & !FLG_SPEED,
        site,
        (flags & FLG_SPEED != 0).then_some(&domains),
        None,
        None,
    )?;
    if body == SUN && flags & FLG_HELCTR != 0 && flags & FLG_BARYCTR == 0 {
        // Unlike Earth, Sun self-observation has no origin sentinel.
        // Read central geometry even when its zero result needs no ayanamsha.
        if central.coordinates[2] == 0.0 {
            return Ok(([0.0; 6], flags));
        }
    }
    let mut output = [
        central.coordinates[0],
        central.coordinates[1],
        central.coordinates[2],
        0.0,
        0.0,
        0.0,
    ];
    if flags & FLG_SPEED != 0 && output[2] != 0.0 {
        let tdb_rate = crate::time::tdb_rate(central_whole + central_fraction)?;
        let bounds = domains.offsets();
        // This linear mapping selects an initial local side/window only.
        // Every full sampled light path independently verifies actual record
        // ownership, so nonlinear emission/closest-approach clocks cannot
        // permit a sample beyond the consumed record's real domain.
        let support = (bounds.0 / tdb_rate, bounds.1 / tdb_rate);
        let site_clock = if flags & FLG_TOPOCTR != 0 && flags & (FLG_HELCTR | FLG_BARYCTR) == 0 {
            Some((
                crate::smh2016::bundled_deltat_piece(central_whole + central_fraction)?,
                local_jpl_site_clock_support(central_whole, central_fraction),
            ))
        } else {
            None
        };
        let cache = std::cell::RefCell::new(std::collections::HashMap::new());
        let sample = |ephemeris_offset: f64, site_offset: f64| -> Result<Option<[f64; 3]>> {
            if ephemeris_offset == 0.0 && site_offset == 0.0 {
                return Ok(central.correction_increment);
            }
            let key = (ephemeris_offset.to_bits(), site_offset.to_bits());
            if let Some(cached) = cache.borrow().get(&key) {
                return Ok(*cached);
            }
            let actual_site_state = if let Some((central_clock, _)) = site_clock {
                let sampled_clock = crate::smh2016::bundled_deltat_piece(
                    central_whole + (central_fraction + site_offset),
                )?;
                if sampled_clock != central_clock {
                    cache.borrow_mut().insert(key, None);
                    return Ok(None);
                }
                Some(local_jpl_topo(
                    central_whole,
                    central_fraction + site_offset,
                    site.expect("the central topocentric place admitted its site"),
                )?)
            } else {
                None
            };
            let admission = LocalJplReadAdmission::new(&domains);
            let sampled = local_jpl_place(
                leg,
                central_whole,
                central_fraction + ephemeris_offset,
                body,
                flags & !FLG_SPEED,
                site,
                None,
                Some(&admission),
                actual_site_state,
            );
            let value = if admission.refused.get() {
                None
            } else {
                let place = sampled?;
                if admission.complete() {
                    place.correction_increment
                } else {
                    None
                }
            };
            cache.borrow_mut().insert(key, value);
            Ok(value)
        };
        let correction_rate = if flags & (FLG_TRUEPOS | FLG_HELCTR | FLG_BARYCTR) != 0
            || flags & (FLG_NOGDEFL | FLG_NOABERR) == (FLG_NOGDEFL | FLG_NOABERR)
        {
            [0.0; 3]
        } else {
            correction_rates_from_independent_clocks(
                &sample,
                support,
                site_clock.map(|(_, bounds)| bounds),
            )?
        };
        let apparent_velocity =
            std::array::from_fn(|i| central.ray_velocity_tt[i] + correction_rate[i]);
        let (matrix, matrix_rate) = super::frame_velocity::native_output_frame_matrix_rate(
            central_whole + central_fraction,
            flags,
        )?;
        let vector_rate = mat3_vec3(matrix, apparent_velocity)?;
        let frame_rate = mat3_vec3(matrix_rate, central.apparent_ray)?;
        let output_velocity = std::array::from_fn(|i| vector_rate[i] + frame_rate[i]);
        [output[3], output[4]] =
            super::derivative::spherical_angular_velocity(central.output_ray, output_velocity)?;
        if let Some(range) = central.range {
            output[5] = range.rate;
        }
    }
    Ok((output, flags))
}

/// Total correction derivative when the ephemeris and rotating-site clocks
/// own independent continuous domains. Each partial keeps the other physical
/// input fixed, while live retarded/closest reads still enforce their actual
/// records through the sampler. No derivative crosses either model boundary;
/// an empty common window therefore does not invent an ill-conditioned step.
fn correction_rates_from_independent_clocks(
    sample: &dyn Fn(f64, f64) -> Result<Option<[f64; 3]>>,
    ephemeris_support: (f64, f64),
    site_support: Option<(f64, f64)>,
) -> Result<[f64; 3]> {
    let ephemeris = super::derivative::cartesian_rates_on_owned_domain(
        &|offset| sample(offset, 0.0),
        1.0 / 64.0,
        ephemeris_support,
    )?;
    let Some(site_support) = site_support else {
        return Ok(ephemeris);
    };
    let site = super::derivative::cartesian_rates_on_owned_domain(
        &|offset| sample(0.0, offset),
        1.0 / 64.0,
        site_support,
    )?;
    Ok(std::array::from_fn(|i| ephemeris[i] + site[i]))
}

/// TT offsets of the bundled Delta-T piece actually used by the local site.
/// The site reads this table once at observation TT. Its interpolation knots,
/// not Gregorian month boundaries or a caller's configured Delta-T, bound the
/// local rotating-observer derivative. IAU frame series are smooth here.
pub(super) fn local_jpl_site_clock_support(whole: f64, fraction: f64) -> (f64, f64) {
    use crate::smh2016_tables::{SMH2016_CURVE_ROWS, SMH2016_TABLE_TT};
    let tt = whole + fraction;
    let knots = &SMH2016_TABLE_TT;
    let relative = |boundary: f64| (boundary - whole) - fraction;
    if knots[0] <= tt && tt <= knots[knots.len() - 1] {
        let index = knots
            .partition_point(|&knot| knot <= tt)
            .saturating_sub(1)
            .min(knots.len() - 2);
        let upper = relative(knots[index + 1]);
        return (
            relative(knots[index]).min(0.0),
            if index + 2 == knots.len() {
                upper.max(0.0)
            } else {
                upper.next_down().max(0.0)
            },
        );
    }
    let mut support = if tt < knots[0] {
        (f64::NEG_INFINITY, relative(knots[0]).next_down())
    } else {
        (relative(knots[knots.len() - 1]).next_up(), f64::INFINITY)
    };
    for row in SMH2016_CURVE_ROWS {
        let boundary = row[0] * 365.25 + 1721045.0;
        if boundary <= tt {
            support.0 = support.0.max(relative(boundary));
        } else {
            support.1 = support.1.min(relative(boundary).next_down());
        }
    }
    (support.0.min(0.0), support.1.max(0.0))
}

/// TT derivative of the actual rotating-site position, including its bundled
/// UT1 and precession/nutation clocks. This derivative channel does not replace
/// the site's nominal velocity consumed by the central aberration formula.
fn local_jpl_topo_velocity_tt(
    whole: f64,
    fraction: f64,
    site: (f64, f64, f64),
) -> Result<[f64; 3]> {
    let piece = crate::smh2016::bundled_deltat_piece(whole + fraction)?;
    rotating_site_velocity_tt(
        local_jpl_site_itrs(site)?,
        &|offset| {
            if crate::smh2016::bundled_deltat_piece(whole + (fraction + offset))? != piece {
                return Ok(None);
            }
            local_jpl_site_orientation(whole, fraction + offset).map(Some)
        },
        local_jpl_site_clock_support(whole, fraction),
    )
}

/// Chain-rule site velocity for `r = Eᵀ Rz(GAST) r_ITRS`, in AU/TT-day.
///
/// The rapid diurnal rotation is differentiated analytically. Only the actual
/// slow equinox matrix and unwrapped sidereal angle/clock are numerically
/// differentiated on their owned support. This avoids fourth-order truncation
/// of the complete sinusoidal site trajectory at a one-sided clock boundary.
/// `None` is an actual clock-domain refusal; calculation errors propagate.
fn rotating_site_velocity_tt(
    itrs: [f64; 3],
    orientation_at: &dyn Fn(f64) -> Result<Option<SiteOrientation>>,
    support: (f64, f64),
) -> Result<[f64; 3]> {
    use reference_frames::*;
    let cache = std::cell::RefCell::new(std::collections::HashMap::new());
    let at = |offset: f64| -> Result<Option<SiteOrientation>> {
        if let Some(value) = cache.borrow().get(&offset.to_bits()).copied() {
            return Ok(value);
        }
        let value = orientation_at(offset)?;
        cache.borrow_mut().insert(offset.to_bits(), value);
        Ok(value)
    };
    let (equinox, gast) = at(0.0)?.ok_or(Error::MathDomain { value: 0.0 })?;
    let rotated = matvec3(rot_z(gast * FRAME_TAU / 24.0), itrs)?;
    let angular_rate = super::derivative::angular_rates_on_owned_domain(
        &|offset| at(offset).map(|value| value.map(|(_, angle)| [angle * 15.0, 0.0])),
        1.0 / 64.0,
        support,
    )?[0]
        .to_radians();
    let mut equinox_rate = [[0.0; 3]; 3];
    for column in 0..3 {
        let rate = super::derivative::cartesian_rates_on_owned_domain(
            &|offset| {
                at(offset).map(|value| {
                    value.map(|(matrix, _)| std::array::from_fn(|row| matrix[row][column]))
                })
            },
            1.0 / 64.0,
            support,
        )?;
        for row in 0..3 {
            equinox_rate[row][column] = rate[row];
        }
    }
    let moving_equinox = matvec3(transpose3(equinox_rate), rotated)?;
    let spin = matvec3(
        transpose3(equinox),
        [-angular_rate * rotated[1], angular_rate * rotated[0], 0.0],
    )?;
    Ok(std::array::from_fn(|i| moving_equinox[i] + spin[i]))
}

/// Admit a retained minor-body SPK before acquiring a main observer image.
/// Missing registration, image or target returns false. Type 2/3 files use
/// their closed whole-file descriptor envelope compared with scalar TT; type
/// 21 uses its half-open usable TDB span with the maximum speed stencil. A range
/// miss is an error, while inside-coverage provider failures remain distinct.
pub fn registered_spk_heliocentric_probe(
    registry: &SpkRegistry,
    jd_tt: f64,
    body: i32,
) -> Result<bool> {
    registered_spk_probe_inner(registry, jd_tt, body, None)
}

/// Admit a direct request at TT with its actual geometric/apparent rate domain.
/// Geometric type-21 requests require stored speed samples but no light time.
pub(crate) fn registered_spk_position_probe(
    registry: &SpkRegistry,
    jd_tt: f64,
    body: i32,
    flags: i32,
) -> Result<bool> {
    registered_spk_probe_inner(registry, jd_tt, body, Some(flags))
}

fn registered_spk_probe_inner(
    registry: &SpkRegistry,
    jd_tt: f64,
    body: i32,
    flags: Option<i32>,
) -> Result<bool> {
    let Some(registration) = registry.body_info(body) else {
        return Ok(false);
    };
    let Some(kernel) = registry.kernel(&registration.kernel) else {
        return Ok(false);
    };
    let type21 = kernel
        .index()
        .all_segments()
        .any(|s| s.descriptor.seg_type == crate::spk::TYPE_EMDA);
    if let Some(flags) = flags.filter(|value| type21 && value & FLG_TRUEPOS != 0) {
        crate::spk::declared_center(kernel.index(), registration.naif_id)?;
        let stencil = crate::spk::type21_request_stencil_days(flags);
        let low = crate::time::tt_to_tdb_jd(jd_tt - stencil)?;
        let high = crate::time::tt_to_tdb_jd(jd_tt + stencil)?;
        let spans = registry.coverage_jd_tdb(body)?;
        if spans.iter().any(|&(start, end)| start <= low && high < end) {
            return Ok(true);
        }
        let (start, end) = spans
            .iter()
            .find(|&&(start, end)| start <= high && low < end)
            .or_else(|| spans.first())
            .copied()
            .expect("registered coverage is nonempty");
        return Err(Error::EphemerisRange {
            requested_jd: crate::time::tt_to_tdb_jd(jd_tt)?,
            start_jd: start,
            end_jd: end,
            body: Some(body),
        });
    }
    let (epoch, span, half_open) = if type21 {
        crate::spk::declared_center(kernel.index(), registration.naif_id)?;
        let whole = jd_tt.floor();
        let fraction = jd_tt - whole;
        let tdb = whole
            + (fraction + crate::time::source_tdb_minus_tt_seconds(whole, fraction) / 86400.0);
        (
            tdb,
            registry.usable_coverage_jd_tdb(
                body,
                flags.map_or_else(
                    crate::spk::type21_default_stencil_days,
                    crate::spk::type21_request_stencil_days,
                ),
            ),
            true,
        )
    } else {
        (
            jd_tt,
            registry.file_coverage_jd_tdb(&registration.kernel),
            false,
        )
    };
    if let Some((start, end)) = span {
        if epoch < start || epoch > end || half_open && epoch == end {
            return Err(Error::EphemerisRange {
                requested_jd: epoch,
                start_jd: start,
                end_jd: end,
                body: Some(body),
            });
        }
    }
    if type21 {
        return Ok(kernel
            .index()
            .segments_for(registration.naif_id)
            .next()
            .is_some());
    }
    if registration.naif_id == crate::spk::CENTER_SSB {
        return Ok(true);
    }
    let mut target = registration.naif_id;
    for _ in 0..32 {
        if kernel.index().segments_for(target).next().is_none() {
            return Ok(false);
        }
        target = crate::spk::resolve_chain(kernel.index(), target)?.center;
        if target == crate::spk::CENTER_SSB {
            return Ok(true);
        }
    }
    Err(Error::SpkCorrupt {
        what: "SPK center cycle".to_string(),
        offset: 0,
    })
}

/// Registered SPK heliocentric place using one already acquired observer store.
/// The target and main stores may differ; the retained target mapping wins.
/// Type 2/3 performs exactly three target reads at freshly constructed scalar
/// TT emission Times and one Sun read at observation TT. True ecliptic-of-date
/// shaping returns degrees/degrees/AU. Type 21 instead retains three split-TDB
/// retardations plus a final read, fixed J2000 ecliptic, Vondrak precession
/// at observation TT and IAU 2000B longitude nutation.
/// No model fallback or implicit main-image acquisition occurs here.
pub fn registered_spk_heliocentric_with_main(
    registry: &SpkRegistry,
    main_registry: &SpkRegistry,
    jd_tt: f64,
    body: i32,
) -> Result<Option<(f64, f64, f64)>> {
    if !registered_spk_heliocentric_probe(registry, jd_tt, body)? {
        return Ok(None);
    }
    let registration = registry
        .body_info(body)
        .expect("probe retained registration");
    let kernel = registry
        .kernel(&registration.kernel)
        .expect("probe retained image");
    if kernel
        .index()
        .all_segments()
        .any(|s| s.descriptor.seg_type == crate::spk::TYPE_EMDA)
    {
        return direct_type21_heliocentric(registry, main_registry, jd_tt, body).map(Some);
    }
    let main = JplLeg::new(main_registry);
    let whole = jd_tt.floor();
    let fraction = jd_tt - whole;
    let observer = local_jpl_source_time_state(&main, SUN, whole, fraction)?.0;
    let relative = minor_spk_relative_at(jd_tt, observer, &|tt| {
        let whole = tt.floor();
        let fraction = tt - whole;
        registry
            .eval_geometric_ssb_au_two_part_jd_with_range(
                body,
                whole,
                fraction + crate::time::source_tdb_minus_tt_seconds(whole, fraction) / 86400.0,
                &mut None,
            )
            .map(|state| state.0)
    })?;
    let tdb =
        whole + (fraction + crate::time::source_tdb_minus_tt_seconds(whole, fraction) / 86400.0);
    let [x, y, z] = true_ecliptic_vector_at(jd_tt, tdb, relative)?;
    let distance = (x * x + y * y + z * z).sqrt();
    let latitude = z.atan2(x.hypot(y)).to_degrees();
    let longitude = y.atan2(x).to_degrees().rem_euclid(360.0);
    Ok(Some((longitude, latitude, distance)))
}

/// Direct type-21 heliocentric output, separate from planetary vector shaping.
/// The reader evaluates center-relative km at split TDB; the independently
/// retained main supplies the declared center and the observation Sun in AU.
/// Three Euclidean-norm light-time iterations precede a fourth target read. Fixed
/// ICRS-to-J2000-ecliptic rotation, spherical reduction, TT Vondrak precession
/// and IAU 2000B nutation produce degrees/degrees/AU. Provider/frame failures
/// propagate; no checked generic planet route is substituted.
fn direct_type21_heliocentric(
    registry: &SpkRegistry,
    main_registry: &SpkRegistry,
    tt: f64,
    body: i32,
) -> Result<(f64, f64, f64)> {
    let registration = registry
        .body_info(body)
        .ok_or(Error::SpkUnregisteredBody { body })?;
    let kernel = registry
        .kernel(&registration.kernel)
        .ok_or(Error::SpkFileNotFound {
            path: registration.kernel.clone(),
        })?;
    let center = crate::spk::declared_center(kernel.index(), registration.naif_id)?;
    let main = JplLeg::new(main_registry);
    let whole = tt.floor();
    let fraction = tt - whole;
    let tdb_fraction =
        fraction + crate::time::source_tdb_minus_tt_seconds(whole, fraction) / 86400.0;
    let observer = local_jpl_source_time_state(&main, SUN, whole, fraction)?.0;
    let center_body = if center == crate::spk::CENTER_SSB {
        None
    } else {
        Some(
            main_registry
                .registered_bodies()
                .into_iter()
                .find(|entry| entry.naif_id == center)
                .ok_or(Error::SpkBodyNotFound { target: center })?
                .body,
        )
    };
    let target = |emission_fraction: f64| -> Result<[f64; 3]> {
        let relative = registry
            .eval_body_relative_two_part_jd(body, whole, emission_fraction)?
            .0
            .position_au();
        let center = match center_body {
            None => return Ok(relative),
            Some(key) => {
                main_registry
                    .eval_geometric_ssb_au_two_part_jd_with_range(
                        key,
                        whole,
                        emission_fraction,
                        &mut None,
                    )?
                    .0
            }
        };
        Ok([
            center[0] + relative[0],
            center[1] + relative[1],
            center[2] + relative[2],
        ])
    };
    let relative = direct_type21_relative_at(tdb_fraction, observer, &target)?;
    let epsilon = crate::spk::OBLIQUITY_J2000_DEG.to_radians();
    let [x, y, z] = relative;
    let ecliptic = [
        x,
        y * epsilon.cos() + z * epsilon.sin(),
        -y * epsilon.sin() + z * epsilon.cos(),
    ];
    let [x, y, z] = ecliptic;
    let radius = (x * x + y * y + z * z).sqrt();
    let longitude = y.atan2(x).to_degrees().rem_euclid(360.0);
    let latitude = if radius > 0.0 {
        (z / radius).asin().to_degrees()
    } else {
        0.0
    };
    let (longitude, latitude) =
        crate::hypothetical::precess_ecliptic_source_angles(longitude, latitude, FRAME_J2000, tt)?;
    let dpsi = crate::erfa::nut00b(tt)?.0.to_degrees();
    Ok(((longitude + dpsi).rem_euclid(360.0), latitude, radius))
}

/// Three split-TDB type-21 retardations followed by the final target read.
/// Original observer and TDB fraction remain fixed. The Euclidean vector norm
/// supplies distance; the type-21 speed-of-light constant has its own identity.
/// Callback failures stop in acquisition order, without a degraded result.
fn direct_type21_relative_at(
    fraction: f64,
    observer: [f64; 3],
    target: &dyn Fn(f64) -> Result<[f64; 3]>,
) -> Result<[f64; 3]> {
    let relative_at = |sample: f64| -> Result<[f64; 3]> {
        let position = target(sample)?;
        Ok([
            position[0] - observer[0],
            position[1] - observer[1],
            position[2] - observer[2],
        ])
    };
    let mut emission = fraction;
    for _ in 0..3 {
        let relative = relative_at(emission)?;
        let distance = reference_frames::series_dot(&relative, &relative)?.sqrt();
        emission = fraction - distance / crate::spk::LIGHT_SPEED_AU_PER_DAY;
    }
    relative_at(emission)
}

/// Three scalar-TT light-time reads for the direct minor-SPK route.
/// The observation Sun vector remains fixed. Each iteration constructs a new
/// scalar TT from the original observation TT minus distance/c, including the
/// final unused constructor; the third relative vector is returned. Target
/// failures propagate in read order; this route adds no aberration/deflection.
fn minor_spk_relative_at(
    tt: f64,
    observer: [f64; 3],
    target: &dyn Fn(f64) -> Result<[f64; 3]>,
) -> Result<[f64; 3]> {
    let mut emission = tt;
    let mut relative = [0.0; 3];
    for _ in 0..3 {
        let position = target(emission)?;
        relative = [
            position[0] - observer[0],
            position[1] - observer[1],
            position[2] - observer[2],
        ];
        let distance =
            (relative[0] * relative[0] + relative[1] * relative[1] + relative[2] * relative[2])
                .sqrt();
        emission = tt - distance / 173.144632674240;
    }
    Ok(relative)
}

/// Registered SPK heliocentric place over one retained registry containing both
/// target and Sun supports. Outputs and failure semantics are those of
/// [`registered_spk_heliocentric_with_main`]; no main is opened implicitly.
pub fn registered_spk_heliocentric(
    registry: &SpkRegistry,
    jd_tt: f64,
    body: i32,
) -> Result<Option<(f64, f64, f64)>> {
    registered_spk_heliocentric_with_main(registry, registry, jd_tt, body)
}

/// Whether a body may use the minor-body SPK/element branch.
///
/// Curated IDs and numbered asteroids strictly between `AST_OFFSET` and
/// `FIXSTAR_OFFSET` participate. Offset sentinels and fixed-star IDs do not.
fn local_jpl_minor_body_candidate(body: i32) -> bool {
    crate::minor_bodies::minor_body_elements(body).is_some()
        || AST_OFFSET < body && body < FIXSTAR_OFFSET
}

/// Whether the registered minor-body SPK covers the requested TT epoch
/// and required rate samples after conversion to TDB. Type-21 usable spans
/// exclude their final endpoint; unsupported span analysis uses the reader's
/// ordinary coverage check. Non-finite epochs do not serve.
fn local_jpl_minor_body_servable(leg: &JplLeg<'_>, body: i32, tt: f64, flags: i32) -> bool {
    if !tt.is_finite() {
        return false;
    }
    let stencil = crate::spk::type21_request_stencil_days(flags);
    if let Some((start, end)) = leg.registry().usable_coverage_jd_tdb(body, stencil) {
        if flags & FLG_TRUEPOS != 0 {
            // Geometric positions never read a retarded target epoch. Admit
            // the actual TT speed stencil after converting both ends to TDB,
            // without the unrelated apparent-place light-time exclusion.
            return registered_spk_position_probe(leg.registry(), tt, body, flags).unwrap_or(false);
        }
        if let Ok(tdb) = crate::time::tt_to_tdb_jd(tt) {
            return start <= tdb && tdb < end;
        }
        return false;
    }
    leg.covers_tdb(body, tt)
}

// Minor-body numerical fallbacks in heliocentric coordinates need no main
// image. A served target instead uses the barycentric SPK pipeline: heliocentric
// output needs the Sun, and a non-SSB center needs its planetary chain even for
// barycentric output. Geocentric projection also needs Earth and Sun.
fn minor_requires_main_image(leg: &JplLeg<'_>, body: i32, tt: f64, flags: i32) -> Result<bool> {
    let registered = local_jpl_minor_body_servable(leg, body, tt, flags);
    if flags & FLG_HELCTR != 0 {
        return Ok(registered);
    }
    if flags & FLG_BARYCTR == 0 || !registered {
        return Ok(true);
    }
    let registration = leg
        .registry()
        .body_info(body)
        .ok_or(Error::SpkUnregisteredBody { body })?;
    let kernel =
        leg.registry()
            .kernel(&registration.kernel)
            .ok_or_else(|| Error::SpkFileNotFound {
                path: registration.kernel,
            })?;
    Ok(crate::spk::declared_center(kernel.index(), registration.naif_id)? != 0)
}

/// J2000 equatorial-to-ecliptic rotation for the IAU 1976 obliquity.
/// The SPICE ECLIPJ2000 orientation uses 84381.448 arcseconds without frame
/// bias, distinct from the IAU 2006 value 84381.406 arcseconds. Entries are
/// dimensionless and row-major.
pub(super) fn iau1976_ecliptic_j2000_matrix() -> [[f64; 3]; 3] {
    [
        [1.0, 0.0, 0.0],
        [0.0, 0.9174820620691818, 0.3977771559319137],
        [0.0, -0.3977771559319137, 0.9174820620691818],
    ]
}

/// Apply the IAU 1976 J2000 ecliptic rotation to an AU or AU/day vector.
fn rotate_icrs_to_iau1976_ecliptic_j2000(v: [f64; 3]) -> [f64; 3] {
    let m = iau1976_ecliptic_j2000_matrix();
    [
        (m[0][0] * v[0] + m[0][1] * v[1]) + m[0][2] * v[2],
        (m[1][0] * v[0] + m[1][1] * v[1]) + m[1][2] * v[2],
        (m[2][0] * v[0] + m[2][1] * v[1]) + m[2][2] * v[2],
    ]
}
/// Select the site for hypothetical-body parallax.
///
/// Geocentric topocentric requests require a configured longitude/latitude in
/// degrees and altitude in metres. Other center selections return no site.
/// Parallax is applied in the frame of date before sidereal/J2000 output.
fn hypothetical_topo_site(
    flags: i32,
    site: Option<(f64, f64, f64)>,
) -> Result<Option<(f64, f64, f64)>> {
    if flags & FLG_TOPOCTR == 0 || flags & (FLG_HELCTR | FLG_BARYCTR) != 0 {
        return Ok(None);
    }
    Ok(Some(site.ok_or(Error::ConfigurationMissing {
        missing_config: "observer_location",
        suggestion: "FLG_TOPOCTR requires a configured site",
    })?))
}

/// IERS 2010 site displacement in J2000 ecliptic axes, in AU and AU/day.
///
/// `tt` is a finite TT JD and site coordinates are degrees/metres. The native
/// IAU ecliptic-J2000 matrix rotates the geocentric offset; no DE states are
/// needed. Non-finite inputs and frame errors propagate.
fn hypothetical_topo_offset_j2000(tt: f64, site: (f64, f64, f64)) -> Result<([f64; 3], [f64; 3])> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    let day = tt.floor();
    let (pos, vel) = local_jpl_topo(day, tt - day, site)?;
    Ok((
        rotate_icrs_to_iau1976_ecliptic_j2000(pos),
        rotate_icrs_to_iau1976_ecliptic_j2000(vel),
    ))
}

/// Recenter a symbolic Selena/Waldemath state on the requested origin.
///
/// The finite epoch is TT JD. Angles are degrees, radius AU and rates per day.
/// Registered Earth support is geometric for barycentric output and retarded
/// relative to the Sun for heliocentric output. Geocentric work can apply
/// first-order annual aberration. The caller applies nutation, sidereal and
/// representation shaping. Waldemath's defining UT1 anchor uses configured
/// Delta-T independently of the observation frame clock. Missing support,
/// used anchor and model errors propagate.
#[cfg(test)]
fn local_jpl_symbolic_center_state(
    leg: &JplLeg<'_>,
    tt: f64,
    body: i32,
    flags: i32,
    model_deltat: &dyn Fn(f64) -> Result<f64>,
) -> Result<PositionResult> {
    local_jpl_symbolic_center_state_via(
        tt,
        body,
        flags,
        model_deltat,
        &|id, whole, fraction| local_jpl_source_time_state(leg, id, whole, fraction),
        &|epoch| retarded_earth_minus_sun_source_time(leg, epoch),
    )
}

/// Select the same symbolic/anchor model and reduce it through owned support
/// callbacks. Epochs, model Delta-T, nutation, flags and center math are shared.
fn local_jpl_symbolic_center_state_via(
    tt: f64,
    body: i32,
    flags: i32,
    model_deltat: &dyn Fn(f64) -> Result<f64>,
    state_at: &VectorStateAt<'_>,
    earth_at: &dyn Fn(f64) -> Result<[f64; 3]>,
) -> Result<PositionResult> {
    let symbolic = if body == WALDEMATH {
        let anchor_ut = crate::hypothetical::WALDEMATH_ANCHOR_JD_UT;
        let anchor_tt = anchor_ut + model_deltat(anchor_ut)?;
        let mean_lon = crate::hypothetical::WALDEMATH_ANCHOR_APPARENT_LON_DEG
            - crate::erfa::nut06a(anchor_tt)?.0.to_degrees();
        crate::hypothetical::waldemath_source_time(tt, anchor_tt, mean_lon)
    } else {
        crate::hypothetical::white_moon_source_time(tt)
    };
    model_serving::symbolic_center_state(
        &|id, epoch| state_at(id, epoch.floor(), epoch - epoch.floor()),
        earth_at,
        tt,
        flags,
        symbolic,
    )
}

/// Reduce a heliocentric predicted-planet model at a TT Julian Day.
///
/// Vulcan/Proserpina use mean-ecliptic-of-date elements; classical predictions
/// use J2000 ecliptic axes. Output is degrees, AU and per-day rates. Centered
/// rates differentiate that reported model. Barycentric output adds the Sun;
/// geocentric output subtracts retarded Earth and may apply PPN deflection and
/// Bradley aberration. Geometric flags suppress these apparent terms.
/// Missing DE support and used frame/model errors propagate.
#[cfg(test)]
fn local_jpl_predicted_center_state(
    leg: &JplLeg<'_>,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<PositionResult> {
    model_serving::predicted_center_state_source_time(
        &|id, epoch| {
            if id != SUN && id != EARTH && local_jpl_target(leg, id).is_err() {
                return Err(Error::UnknownBody { body_id: id });
            }
            local_jpl_source_time_state(leg, id, epoch, 0.0)
        },
        &|epoch| {
            retarded_earth_observation_source_time_via(leg, epoch, &|target, day, fraction| {
                leg.registry()
                    .eval_geometric_ssb_au_two_part_jd(target, day, fraction)
            })
        },
        false,
        tt,
        body,
        flags,
    )
}

/// Subtract light time in days from a two-part TDB epoch's fraction.
///
/// Retain the whole-day part unchanged, including when the resulting fraction
/// crosses a day boundary. Every input and the resulting fraction must be
/// finite. Keeping the delay at fraction magnitude avoids full-JD cancellation.
fn retarded_tdb_split(whole: f64, tdb_fraction: f64, light_time: f64) -> Result<(f64, f64)> {
    for value in [whole, tdb_fraction, light_time] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let retarded = tdb_fraction - light_time;
    if !retarded.is_finite() {
        return Err(Error::NonFiniteArgument { value: retarded });
    }
    Ok((whole, retarded))
}

/// Observe Earth from a Sun position fixed at the reception epoch.
///
/// `tt` is a finite TT Julian Day. Convert its day/fraction split to TDB, then
/// subtract each light-time delay directly from that fractional component.
/// The solve uses the exact-SI AU/day light speed and stops when consecutive
/// delays differ by less than 1e-12 days, for at most ten iterations.
/// The returned ray is ICRS AU. Missing support and provider errors propagate.
fn retarded_earth_minus_sun(leg: &JplLeg<'_>, tt: f64) -> Result<[f64; 3]> {
    use reference_frames::FRAME_DAY_S;
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    // Retain a whole-day and fractional TT epoch before converting to TDB.
    let whole = tt.floor();
    let tt_fraction = tt - whole;
    let tdb_fraction =
        tt_fraction + reference_frames::tdb_minus_tt_seconds(whole, tt_fraction)? / FRAME_DAY_S;
    let sun_key = local_jpl_target(leg, SUN)?;
    let (sun, _) =
        leg.registry()
            .eval_geometric_ssb_au_two_part_jd(sun_key, whole, tdb_fraction)?;
    let earth_key = local_jpl_target(leg, EARTH)?;
    let earth = leg
        .registry()
        .eval_geometric_ssb_au_two_part_jd(earth_key, whole, tdb_fraction)?
        .0;
    retarded_earth_minus_sun_samples(whole, tdb_fraction, sun, earth, &|day, fraction| {
        leg.registry()
            .eval_geometric_ssb_au_two_part_jd(earth_key, day, fraction)
            .map(|state| state.0)
    })
}

/// Test adapter observing Earth from the Sun at an already admitted TT epoch.
/// The same split-TDB geometric acquisition and light-time convergence
/// criterion apply as in [`retarded_earth_minus_sun`]. Provider errors retain
/// their actual record metadata for the outer calculation to classify.
#[cfg(test)]
fn retarded_earth_minus_sun_source_time(leg: &JplLeg<'_>, tt: f64) -> Result<[f64; 3]> {
    retarded_earth_minus_sun_source_time_via(leg, tt, &|target, day, fraction| {
        leg.registry()
            .eval_geometric_ssb_au_two_part_jd(target, day, fraction)
    })
}

/// Inject immutable ICRS AU/AU-day support for the shared split-epoch propagation solve.
fn retarded_earth_minus_sun_source_time_via(
    leg: &JplLeg<'_>,
    tt: f64,
    state_at: &VectorStateAt<'_>,
) -> Result<[f64; 3]> {
    retarded_earth_observation_source_time_via(leg, tt, state_at)
        .map(|observation| observation.relative)
}

/// Keep the already acquired reception Sun beside the retarded Earth ray.
fn retarded_earth_observation_source_time_via(
    leg: &JplLeg<'_>,
    tt: f64,
    state_at: &VectorStateAt<'_>,
) -> Result<model_serving::RetardedEarthObservation> {
    let whole = tt.floor();
    let fraction = tt - whole;
    let tdb_fraction =
        fraction + crate::time::source_tdb_minus_tt_seconds(whole, fraction) / 86400.0;
    let sun_key = local_jpl_target(leg, SUN)?;
    let sun = state_at(sun_key, whole, tdb_fraction)?.0;
    let earth_key = local_jpl_target(leg, EARTH)?;
    let earth = state_at(earth_key, whole, tdb_fraction)?.0;
    let relative = retarded_earth_minus_sun_samples_with_split(
        whole,
        tdb_fraction,
        sun,
        earth,
        &|day, fraction| state_at(earth_key, day, fraction).map(|s| s.0),
        // The emission callback already owns the retarded split TDB epoch.
        &|day, fraction, light_time| Ok((day, fraction - light_time)),
    )?;
    Ok(model_serving::RetardedEarthObservation {
        relative,
        sun_reception: sun,
    })
}

/// Retard an already acquired heliocentric Earth ray at an original TDB pair.
/// `sun` and `earth` are SSB ICRS positions in AU. The fixed Sun observer is
/// subtracted from each Earth sample, with the existing scalar squared norm.
/// Ten successive delay updates must settle within 1e-12 day; exhaustion is
/// a convergence error, never an unconverged position. The injected reader
/// receives retained whole/fraction TDB parts after each delay subtraction;
/// its errors propagate immediately, and no output frame is applied here.
fn retarded_earth_minus_sun_samples(
    whole: f64,
    tdb_fraction: f64,
    sun: [f64; 3],
    earth: [f64; 3],
    state_at: &dyn Fn(f64, f64) -> Result<[f64; 3]>,
) -> Result<[f64; 3]> {
    retarded_earth_minus_sun_samples_with_split(
        whole,
        tdb_fraction,
        sun,
        earth,
        state_at,
        &retarded_tdb_split,
    )
}

/// Solve Earth-from-Sun propagation with an explicit retarded-epoch callback.
/// States are ICRS AU and epoch parts are TDB days. Return only a converged
/// ray; acquisition or convergence errors propagate unchanged.
fn retarded_earth_minus_sun_samples_with_split(
    whole: f64,
    tdb_fraction: f64,
    sun: [f64; 3],
    mut earth: [f64; 3],
    state_at: &dyn Fn(f64, f64) -> Result<[f64; 3]>,
    split_at: &dyn Fn(f64, f64, f64) -> Result<(f64, f64)>,
) -> Result<[f64; 3]> {
    use reference_frames::FRAME_LIGHT_SPEED_AU_DAY_LOCAL as C_FRAME_AU_DAY;
    let mut light_time0 = 0.0;
    let mut last_delta = f64::NAN;
    for _ in 0..10 {
        let distance = ((earth[0] - sun[0]).powi(2)
            + (earth[1] - sun[1]).powi(2)
            + (earth[2] - sun[2]).powi(2))
        .sqrt();
        let light_time = distance / C_FRAME_AU_DAY;
        last_delta = light_time - light_time0;
        if last_delta.abs() < 1e-12 {
            return Ok([earth[0] - sun[0], earth[1] - sun[1], earth[2] - sun[2]]);
        }
        // Evaluate the emitted Earth state at its retained TDB day and fraction.
        let (whole_ret, tdb_retarded) = split_at(whole, tdb_fraction, light_time)?;
        earth = state_at(whole_ret, tdb_retarded)?;
        light_time0 = light_time;
    }
    Err(Error::ConvergenceFailed {
        algorithm: "spk_light_time",
        iterations: 10,
        max_iterations: 10,
        tolerance: 1e-12,
        last_value: last_delta,
    })
}

/// True-ecliptic-of-date heliocentric Earth vector from registered DE
/// states: the light-time-retarded Earth-minus-Sun vector rotated into date axes with
/// the shared ICRS-to-date construction. Scalar `jd` is TT JD; AU
/// components. Non-finite epochs and unregistered support propagate.
fn earth_helio_true_of_date(leg: &JplLeg<'_>, jd: f64) -> Result<[f64; 3]> {
    if !jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd });
    }
    // This scalar minor-body route resolves Earth before Sun. Other users of
    // the shared observation primitive retain their Sun-first acquisition.
    local_jpl_target(leg, EARTH)?;
    let relative = retarded_earth_minus_sun(leg, jd)?;
    let tdb = jd + reference_frames::tdb_minus_tt_seconds(jd, 0.0)? / 86400.0;
    true_ecliptic_vector_at(jd, tdb, relative)
}

/// Rotate an ICRS AU ray to the true ecliptic of date.
/// TT supplies IAU 2000A nutation; TDB supplies P03 precession and mean
/// obliquity. The caller provides both observation clocks. Non-finite input
/// and frame errors propagate.
fn true_ecliptic_vector_at(tt: f64, tdb: f64, relative: [f64; 3]) -> Result<[f64; 3]> {
    let eps = (reference_frames::mean_obliquity_arcsec(tdb)? / 3600.0).to_radians();
    let (psi, deps) = reference_frames::iau2000a_tenth_microarcsec(tt)?;
    let scale = reference_frames::FRAME_TENTH_USEC_2_RAD;
    let nut = reference_frames::nutation_matrix(eps, eps + deps * scale, psi * scale)?;
    let equinox = reference_frames::equinox_matrix(nut, reference_frames::precession_matrix(tdb)?)?;
    let frame = reference_frames::ecliptic_date_matrix(equinox, eps, deps * scale)?;
    reference_frames::matvec3(frame, relative)
}

/// One minor-body Keplerian sample at scalar TT: curated heliocentric
/// true-of-date place, or its geocentric composition with the date-frame
/// Earth vector. Barycentric requests take the geocentric arm.
/// Nutation stripping follows the selected nutation-free flags
/// (NONUT/J2000/SIDEREAL+EQUATORIAL). Angles are degrees, distance AU.
fn keplerian_position_at(
    leg: &JplLeg<'_>,
    jd_tt: f64,
    body: i32,
    flags: i32,
) -> Result<(f64, f64, f64)> {
    keplerian_position_with_earth(jd_tt, body, flags, &|tt| earth_helio_true_of_date(leg, tt))
}

/// One curated Keplerian sample over a caller-owned Earth observation.
/// `earth_at` returns the light-time-retarded heliocentric Earth in the
/// true ecliptic of date, in AU at scalar TT. Heliocentric requests do
/// not consume it; every other center uses the geocentric reduction.
/// Model evaluation, nutation stripping and arithmetic are shared with
/// the registered-SPK wrapper, without changing the provider boundary.
/// Returns longitude in [0,360) degrees, latitude in degrees and radius
/// in AU. NONUT, J2000, or the combination SIDEREAL|EQUATORIAL removes
/// the longitude nutation term; other requests retain true-of-date axes.
/// Unknown model IDs, non-finite TT and required observer errors propagate.
fn keplerian_position_with_earth(
    jd_tt: f64,
    body: i32,
    flags: i32,
    earth_at: &dyn Fn(f64) -> Result<[f64; 3]>,
) -> Result<(f64, f64, f64)> {
    let dpsi = crate::erfa::nut06a(jd_tt)?.0;
    let (lon_hel, lat_hel, r_hel) =
        crate::minor_bodies::calc_minor_body_heliocentric(body, jd_tt, false, None, dpsi)?;
    keplerian_sample_with_earth(jd_tt, flags, dpsi, (lon_hel, lat_hel, r_hel), earth_at)
}

/// Shape a heliocentric true-of-date minor-body place over scalar TT.
///
/// `heliocentric` is longitude/latitude in degrees and radius AU; `dpsi` is
/// longitude nutation in radians at the same TT. Heliocentric flags bypass
/// Earth. Every other center subtracts the caller's true-of-date retarded
/// Earth vector before spherical reduction. Nutation-free flags remove dpsi
/// after center composition. The returned units are degrees, degrees and AU;
/// required Earth errors propagate without another provider or frame choice.
fn keplerian_sample_with_earth(
    jd_tt: f64,
    flags: i32,
    dpsi: f64,
    heliocentric: (f64, f64, f64),
    earth_at: &dyn Fn(f64) -> Result<[f64; 3]>,
) -> Result<(f64, f64, f64)> {
    let (lon_hel, lat_hel, r_hel) = heliocentric;
    let nut_free = flags & (FLG_NONUT | FLG_J2000) != 0
        || flags & (FLG_SIDEREAL | FLG_EQUATORIAL) == (FLG_SIDEREAL | FLG_EQUATORIAL);
    if flags & FLG_HELCTR != 0 {
        let lon = if nut_free {
            (lon_hel - dpsi.to_degrees()).rem_euclid(360.0)
        } else {
            lon_hel
        };
        return Ok((lon, lat_hel, r_hel));
    }
    let lon_r = lon_hel.to_radians();
    let lat_r = lat_hel.to_radians();
    let cl = lat_r.cos();
    let helio = [
        r_hel * cl * lon_r.cos(),
        r_hel * cl * lon_r.sin(),
        r_hel * lat_r.sin(),
    ];
    let earth = earth_at(jd_tt)?;
    let geo = [
        helio[0] - earth[0],
        helio[1] - earth[1],
        helio[2] - earth[2],
    ];
    let r_geo = (geo[0] * geo[0] + geo[1] * geo[1] + geo[2] * geo[2]).sqrt();
    let lon = geo[1].atan2(geo[0]).to_degrees().rem_euclid(360.0);
    let lat = if r_geo > 0.0 {
        (geo[2] / r_geo).asin().to_degrees()
    } else {
        0.0
    };
    let lon = if nut_free {
        (lon - dpsi.to_degrees()).rem_euclid(360.0)
    } else {
        lon
    };
    Ok((lon, lat, r_geo))
}

/// Evaluate one numbered asteroid using curated or SBDB orbital elements.
///
/// `body` is `AST_OFFSET + catalog_number`, below `FIXSTAR_OFFSET`.
/// Element acquisition precedes required frame and Earth geometry. Missing
/// or transport-failed elements become [`Error::UnknownBody`]; sealed-network
/// refusal propagates. The element model does not perform a second SPK lookup.
/// TT, center selection, degrees/AU output and rates follow
/// [`keplerian_sample_with_earth`].
fn numbered_asteroid_position_with_earth(
    jd_tt: f64,
    body: i32,
    flags: i32,
    sbdb_elements: crate::minor_bodies::SbdbElementsProvider<'_>,
    earth_at: &dyn Fn(f64) -> Result<[f64; 3]>,
) -> Result<(f64, f64, f64)> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let number = i64::from(body) - i64::from(AST_OFFSET);
    let fetched = if crate::minor_bodies::resolve_asteroid_number(number)?.is_none() {
        Some(sbdb_elements(number))
    } else {
        None
    };
    let preloaded = |number| match &fetched {
        Some(result) => result.clone(),
        None => sbdb_elements(number),
    };
    if matches!(&fetched, Some(result) if !matches!(result, Ok(Some(_)))) {
        // The existing numbered-asteroid kernel owns the exception mapping.
        // A failed acquisition cannot reach its shaping or nutation legs.
        return crate::minor_bodies::calc_asteroid_by_number(
            number,
            jd_tt,
            false,
            None,
            Some(&preloaded),
            0.0,
        );
    }
    let dpsi = crate::erfa::nut06a(jd_tt)?.0;
    let helio = crate::minor_bodies::calc_asteroid_by_number(
        number,
        jd_tt,
        false,
        None,
        Some(&preloaded),
        dpsi,
    )?;
    keplerian_sample_with_earth(jd_tt, flags, dpsi, helio, earth_at)
}

/// Numbered-asteroid fallback over the same one-second SPEED stencil and
/// registered Earth support as curated models. Element lookup uses the
/// injected per-entry SBDB provider; frames and output shaping stay with the
/// public caller. Required source/model/observer errors propagate unchanged.
fn local_jpl_numbered_asteroid_fallback(
    leg: &JplLeg<'_>,
    tt: f64,
    body: i32,
    flags: i32,
    sbdb_elements: crate::minor_bodies::SbdbElementsProvider<'_>,
) -> Result<PositionResult> {
    keplerian_fallback_with_samples(tt, flags, &|jd| {
        numbered_asteroid_position_with_earth(jd, body, flags, sbdb_elements, &|epoch| {
            earth_helio_true_of_date(leg, epoch)
        })
    })
}

/// Minor-body Keplerian fallback place: the composed sample with a
/// centered one-second difference stencil for SPEED.
/// The stencil differentiates the reported (nutated) longitude, so no
/// separate nutation-rate term is added. Only the longitude rate folds.
/// Sidereal/equatorial/topocentric reduction belongs to the caller; the
/// fallback itself carries no parallax.
fn local_jpl_keplerian_fallback(
    leg: &JplLeg<'_>,
    tt: f64,
    body: i32,
    flags: i32,
) -> Result<PositionResult> {
    keplerian_fallback_with_samples(tt, flags, &|jd| keplerian_position_at(leg, jd, body, flags))
}

/// Raw six-component curated Keplerian fallback over scalar TT samples.
/// The injected observer follows [`keplerian_position_with_earth`]. SPEED
/// consumes the central and ±one-second places, with wrapped longitude
/// folding after division; absent SPEED leaves exact zero rates. Errors
/// propagate from every required sample. Frames, sidereal and output
/// representation remain the caller's tail.
fn keplerian_fallback_with_earth(
    tt: f64,
    body: i32,
    flags: i32,
    earth_at: &dyn Fn(f64) -> Result<[f64; 3]>,
) -> Result<PositionResult> {
    keplerian_fallback_with_samples(tt, flags, &|jd| {
        keplerian_position_with_earth(jd, body, flags, earth_at)
    })
}

/// Common one-second finite difference over complete Keplerian places.
/// Sample acquisition owns the provider and flags; derivative arithmetic
/// and error propagation are identical for scalar LEB and split SPK support.
/// Outputs are longitude/latitude degrees and radius AU followed by rates
/// in degrees/day and AU/day. For TT K and h=1/86400 day, SPEED evaluates
/// q'=[q(K+h)-q(K-h)]/(2h); longitude folds by 360/(2h) when its rate
/// crosses ±180/(2h). Without SPEED all rates are positive zero.
/// Any error in the central or required side samples propagates; no accuracy
/// or convergence bound is inferred from this finite-difference operation.
fn keplerian_fallback_with_samples(
    tt: f64,
    flags: i32,
    sample: &dyn Fn(f64) -> Result<(f64, f64, f64)>,
) -> Result<PositionResult> {
    const DT_DAYS: f64 = 1.0 / 86400.0;
    let (lon, lat, dist) = sample(tt)?;
    let (dlon, dlat, ddist) = if flags & FLG_SPEED != 0 {
        let prev = sample(tt - DT_DAYS)?;
        let next = sample(tt + DT_DAYS)?;
        let state = crate::hypothetical::keplerian_state_from_samples((lon, lat, dist), prev, next);
        (state.dlon, state.dlat, state.ddist)
    } else {
        (0.0, 0.0, 0.0)
    };
    Ok([lon, lat, dist, dlon, dlat, ddist])
}

/// Subtract an ecliptic sidereal offset from a tropical SPK result.
///
/// The observation is TT JD. The supplied ayanamsha is degrees at UT1; an
/// original UT1 entry is retained centrally, while other samples use the
/// bundled timescale. Configured model Delta-T applies to defining epochs.
/// A one-second scalar-UT1 stencil supplies the daily drift. Equatorial or
/// zero-origin output does not acquire ayanamsha. Clock, anchor and invalid
/// numerical errors propagate; other position components are unchanged.
fn local_jpl_sidereal_word(
    output: &mut PositionResult,
    tt: f64,
    original_ut: Option<f64>,
    flags: i32,
    ayanamsha_at: &dyn Fn(f64) -> Result<f64>,
) -> Result<()> {
    if flags & FLG_SIDEREAL == 0 || flags & FLG_EQUATORIAL != 0 || output[2] == 0.0 {
        return Ok(());
    }
    let ut = native_model_ut_at(tt, original_ut)?;
    let (lon, rate) =
        crate::planetary_moons::moon_sidereal_word(output[0], output[3], ut, flags, ayanamsha_at)?;
    output[0] = lon;
    output[3] = rate;
    Ok(())
}

/// Six-component vector position differentiated on the supplied place model.
///
/// Inputs and provider obligations follow [`vector_place_from_states`]. SPEED
/// samples retain `whole` and adjust only the small TT `fraction`. Angular
/// rates use the adaptive fourth-order physical-place derivative; the legacy
/// state-only callback supplies no additional record metadata. Registered SPK
/// public calls add their actual record admission in `local_jpl_position` and
/// use the analytic light-path range derivative. This low-level helper retains
/// its finite-width range derivative. Output is degrees/degrees/AU followed by
/// per-TT-day rates, with exact zero non-SPEED tails. A zero-distance central
/// place requires no side read. Callback/provider errors propagate unchanged.
#[cfg(test)]
pub(super) fn vector_position_from_states(
    state_at: &VectorStateAt<'_>,
    whole: f64,
    fraction: f64,
    body: i32,
    flags: i32,
    site: Option<(f64, f64, f64)>,
) -> Result<PositionResult> {
    vector_position_from_places(
        &|w, f| vector_place_from_states(state_at, w, f, body, flags & !FLG_SPEED, site),
        whole,
        fraction,
        body,
        flags,
    )
}

/// Differentiate the actual scalar place callback at split TT epochs.
/// The callback must enforce its own physical domain; public SPK calls use
/// explicit consumed-record admission as well. No side read serves an origin.
#[cfg(test)]
pub(super) fn vector_position_from_places(
    place_at: &dyn Fn(f64, f64) -> Result<[f64; 3]>,
    whole: f64,
    fraction: f64,
    body: i32,
    flags: i32,
) -> Result<PositionResult> {
    let pos = place_at(whole, fraction)?;
    let mut output = [pos[0], pos[1], pos[2], 0.0, 0.0, 0.0];
    if flags & FLG_SPEED != 0 && pos[2] != 0.0 {
        [output[3], output[4]] = super::derivative::angular_rates(
            &|offset| place_at(whole, fraction + offset).map(|value| [value[0], value[1]]),
            1.0 / 64.0,
            (f64::NEG_INFINITY, f64::INFINITY),
        )?;
        let h = if body == MOON {
            6.0 / 86400.0
        } else {
            1.0 / 86400.0
        };
        let minus = place_at(whole, fraction - h)?;
        let plus = place_at(whole, fraction + h)?;
        output[5] = (plus[2] - minus[2]) / (2.0 * h);
    }
    Ok(output)
}

/// Adapter from an immutable local SPK registry to geometric states.
///
/// The center-first alias policy is shared with native JPL calculations.
/// Scalar and split TT inputs are converted to split TDB, retained through
/// record evaluation. Each relative segment converts km/km/s to AU/AU-day
/// before center-chain composition. Reported coverage is a TT envelope of
/// stored TDB spans; actual gaps remain checked by the reader.
/// Construction refuses sealed LEB mode and never mutates the registry.
pub(crate) struct LocalSpkStateProvider<'a> {
    registry: &'a SpkRegistry,
}

/// Public trace identity of a successfully acquired local SPK target.
/// Standard planets retain their historical public label; explicitly
/// registered minor bodies and planetary satellites identify the SPK route.
fn local_spk_trace_source(body: i32) -> &'static str {
    if local_jpl_minor_body_candidate(body) || crate::planetary_moons::is_planetary_moon(body) {
        "SPK"
    } else {
        "Skyfield"
    }
}

impl<'a> LocalSpkStateProvider<'a> {
    /// Apparent ray for an arbitrary zero-motion J2000 star.
    ///
    /// RA and declination are degrees; the scalar angles are converted through
    /// hours and turns. Declination is not range-clamped: finite angles beyond
    /// the poles retain their trigonometric meaning. The model has zero proper
    /// motion, radial velocity and parallax; the stellar parallax floor is
    /// 1e-6 milliarcseconds. `epoch` is the original global Time split in TT,
    /// `observer` its SSB ICRS AU position, `velocity` its AU/day velocity and
    /// `surface` its GCRS site vector. The star propagates to scalar TDB plus
    /// projected light-time difference before observer subtraction. Deflection
    /// and aberration retain the original astrometric light time. Output is an
    /// observer-relative ICRS AU ray. Finite input/state checks and native raw
    /// provider acquisition errors propagate without public scalar range guards.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn arbitrary_star_apparent_vector(
        &self,
        epoch: TwoPartEpoch,
        ra_deg: f64,
        dec_deg: f64,
        observer: [f64; 3],
        velocity: [f64; 3],
        surface: [f64; 3],
    ) -> Result<[f64; 3]> {
        check_finite_vec3([epoch.hi, epoch.lo, epoch.value()])?;
        check_finite_vec3(observer)?;
        check_finite_vec3(velocity)?;
        check_finite_vec3(surface)?;
        check_finite_vec3([ra_deg, dec_deg, 0.0])?;
        let data = crate::stars::StarData::from_row(ra_deg, dec_deg, 0.0, 0.0, 0.0, 0.0);
        let (position, stellar_velocity) =
            crate::stars::star_icrs_initial_state_catalog_angles(&data);
        let distance = vec3_dist(position);
        let unit: [f64; 3] = std::array::from_fn(|i| position[i] / distance);
        let projection = (unit[0] * observer[0] + unit[1] * observer[1] + unit[2] * observer[2])
            / FRAME_LIGHT_SPEED_AU_DAY;
        let tdb = epoch.hi
            + (epoch.lo + reference_frames::tdb_minus_tt_seconds(epoch.hi, epoch.lo)? / 86400.0);
        let span = tdb + projection - crate::constants::J2000;
        let vector = std::array::from_fn(|i| {
            let propagated = position[i] + stellar_velocity[i] * span;
            propagated - observer[i]
        });
        let light_time = vec3_dist(vector) / FRAME_LIGHT_SPEED_AU_DAY;
        let state_at = |id, time: TwoPartEpoch| {
            self.geometric_state_two_part(id, time)
                .map(|state| state.pos())
        };
        native_apparent_ray(
            &state_at,
            epoch,
            vector,
            light_time,
            observer,
            velocity,
            Some(surface),
        )
    }

    /// Bind a local registry under an explicit non-sealed calculation mode.
    /// Empty or partial registries are allowed; body evaluation retains its
    /// typed missing-registration/coverage errors, rather than inventing data.
    pub(crate) fn new(registry: &'a SpkRegistry, mode: ProviderMode) -> Result<Self> {
        forbid_jpl_in_sealed_mode(mode)?;
        Ok(Self { registry })
    }

    /// Raw retained-Time SSB state for a catalog-star observer or apparent support.
    ///
    /// `epoch` holds TT whole/fraction parts. The seven-term TDB correction feeds
    /// the split directory-index/polynomial evaluator; per-segment AU conversion
    /// precedes center-chain addition. Output is ICRS position AU and velocity
    /// AU/day with center-first target alias selection. The registry borrow is
    /// immutable and call-scoped; the returned vectors own their scalar words.
    ///
    /// Whole-chain, alias, frame/type/directory errors and derived IEEE arithmetic
    /// retain the raw acquisition contract. This helper does not apply the public
    /// scalar planet range envelope, caller model Delta-T or a LEB substitution.
    /// Original caller epoch admission belongs to the enclosing star entry, and
    /// public range promotion belongs to its explicit outer callback.
    pub(super) fn star_source_time_state(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<([f64; 3], [f64; 3])> {
        local_jpl_source_time_state(&JplLeg::new(self.registry), body, epoch.hi, epoch.lo)
    }

    /// Node/apse geometric state with the actual failed segment retained.
    ///
    /// Inputs are finite split TT Julian days; the checked TDB correction and
    /// stored directory sampling retain AU/AU-day conversion and center-first
    /// composition. The geometric entry keeps its original frame, body and
    /// scalar epoch. Only a directory refusal attaches owned image/descriptor
    /// facts; aliases, unsupported types and other errors stay unchanged.
    /// Strict public SPK and scalar calculation envelopes are separate routes.
    pub(crate) fn nodaps_state_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<GeometricState> {
        self.source_state_two_part_captured(body, epoch)
    }

    /// Physical observation state at its original split TT, in AU/AU-day.
    /// The actual selected descriptor owns a raw directory refusal; successful
    /// checked state construction and component/frame arithmetic are shared.
    pub(crate) fn pheno_state_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
        source_time: bool,
    ) -> Result<GeometricState> {
        if source_time {
            self.source_state_two_part_captured(body, epoch)
        } else {
            self.geometric_state_two_part(body, epoch)
        }
    }

    /// Shared checked state builder for callers that retain failed physical
    /// source descriptors. The public scalar-TT envelope remains separate.
    fn source_state_two_part_captured(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<GeometricState> {
        let entry = self.serving_entry_snapshot(body)?;
        let target = local_jpl_target(&JplLeg::new(self.registry), body)?;
        let correction = reference_frames::tdb_minus_tt_seconds(epoch.hi, epoch.lo)?;
        let fraction = epoch.lo + correction / 86400.0;
        check_finite_vec3([epoch.hi, fraction, epoch.hi + fraction])?;
        let mut captured = None;
        let result = self.registry.eval_geometric_ssb_au_two_part_jd_with_range(
            target,
            epoch.hi,
            fraction,
            &mut captured,
        );
        let (pos, vel) = result.map_err(|error| match captured {
            Some(raw) => error.with_backend_range(raw),
            None => error,
        })?;
        build_state(
            body,
            epoch.value(),
            &entry,
            &pos,
            &vel,
            self.selected_tier(body, epoch.value()),
            |_| None,
        )
    }

    /// Promote a consumed raw fixed-star range error using the main Sun image.
    /// `tt` is the original scalar TT of the enclosing scalar star sample,
    /// shared-observer batch or dynamic-anchor request. The main planetary
    /// descriptor envelope is independent of the failed Earth/deflector image
    /// and observing center; `_center` is retained for the existing callback
    /// signature. Only raw SpkOutOfRange is promoted, with absent star-body id.
    /// Missing aliases and all non-range errors remain unchanged. The numeric
    /// error contains no filename, calendar, message or cause fields.
    pub(super) fn star_range_error(&self, _center: i32, tt: f64, error: Error) -> Error {
        if !matches!(&error, Error::SpkOutOfRange { .. }) {
            return error;
        }
        match self.public_range_error(SUN, tt, error) {
            Error::EphemerisRange {
                requested_jd,
                start_jd,
                end_jd,
                ..
            } => Error::EphemerisRange {
                requested_jd,
                start_jd,
                end_jd,
                body: None,
            },
            other => other,
        }
    }

    /// Resolve one native vector alias without evaluating an epoch. Planet
    /// centers retain the existing center-first/barycenter selection; Earth
    /// here means the geocenter (NAIF 399), independently of orbital EMB.
    /// Missing aliases propagate before state/type/frame/directory reads.
    pub(crate) fn validate_source_alias(&self, body: i32) -> Result<()> {
        local_jpl_target(&JplLeg::new(self.registry), body).map(|_| ())
    }

    /// Ordered target-minus-center vector at the original split TT Time.
    ///
    /// The optional geocenter-to-observer AU segment is negated first, then
    /// center constituents in reverse order and target constituents forward.
    /// Raw SPK directory domains, complete-chain guards and per-segment AU
    /// conversion are retained; no completed barycentric sums are subtracted.
    /// Only position is returned. Registration, finite and raw provider errors
    /// propagate before any output frame or scalar-TT envelope is applied.
    pub(crate) fn pheno_geometric_difference(
        &self,
        body: i32,
        center: i32,
        epoch: TwoPartEpoch,
        surface: Option<[f64; 3]>,
        source_time: bool,
    ) -> Result<[f64; 3]> {
        let leg = JplLeg::new(self.registry);
        let center_key = local_jpl_target(&leg, center)?;
        let target_key = local_jpl_target(&leg, body)?;
        let fraction =
            epoch.lo + reference_frames::tdb_minus_tt_seconds(epoch.hi, epoch.lo)? / 86400.0;
        let mut captured = None;
        let result = if source_time {
            self.registry
                .eval_geometric_difference_with_observer_and_range(
                    target_key,
                    center_key,
                    epoch.hi,
                    fraction,
                    surface.map(|position| (position, [0.0; 3])),
                    &mut captured,
                )
        } else {
            self.registry
                .eval_geometric_difference_with_observer_au_two_part_jd(
                    target_key,
                    center_key,
                    epoch.hi,
                    fraction,
                    surface.map(|position| (position, [0.0; 3])),
                )
        };
        result
            .map_err(|error| match captured {
                Some(raw) => error.with_backend_range(raw),
                None => error,
            })
            .map(|state| state.0)
    }

    /// Acquire a geometric, astrometric or apparent native ICRS target ray.
    ///
    /// The split TT epoch and finite SSB observer state use AU/AU-day. Geometric
    /// requests skip retardation and corrections; otherwise target emission is
    /// solved at split TDB with ten iterations and a 1e-12-day update criterion.
    /// `apparent` additionally applies strict planetary deflection, optional
    /// surface-Earth deflection and Lorentz aberration. No frame projection or
    /// public scalar-JD coverage promotion occurs here. Actual chain, segment,
    /// state and convergence errors propagate without substitution.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn pheno_observed_vector(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
        observer: [f64; 3],
        velocity: [f64; 3],
        surface: Option<[f64; 3]>,
        geometric: bool,
        apparent: bool,
        source_time: bool,
    ) -> Result<[f64; 3]> {
        self.validate_source_alias(body)?;
        let state_at = |id, time: TwoPartEpoch| {
            self.pheno_state_two_part(id, time, source_time)
                .map(|state| state.pos())
        };
        let retarded_at = |time: TwoPartEpoch, delay| {
            let fraction =
                time.lo + reference_frames::tdb_minus_tt_seconds(time.hi, time.lo)? / 86400.0;
            let epoch = TwoPartEpoch {
                hi: time.hi,
                lo: fraction - delay,
            };
            if source_time {
                self.pctr_source_tdb_state(body, epoch)
            } else {
                self.state_tdb_two_part(body, epoch)
            }
            .map(|state| state.0)
        };
        native_pheno_observe(
            &state_at,
            &retarded_at,
            epoch,
            body,
            observer,
            velocity,
            surface,
            geometric,
            apparent,
        )
    }

    /// Apparent target ray for an observer on the bundled UT1/TT timescale.
    ///
    /// `epoch` retains split TT. Observer position/velocity are finite SSB ICRS
    /// AU and AU/day. Optional `surface` is the geocentric site vector in AU,
    /// used for Earth's surface deflection. Target retardation uses split TDB
    /// and a ten-iteration, 1e-12-day light-time update criterion. Closest-approach
    /// deflectors use their actual TT epochs. Strict planetary deflection, optional
    /// Earth deflection and aberration precede any frame or horizon projection.
    ///
    /// The returned ray is observer-relative ICRS AU. Missing aliases, coverage,
    /// unsupported segments, corrupt data and convergence errors propagate.
    /// Acquisition neither substitutes LEB nor changes the caller's kernel.
    pub(crate) fn eclipse_source_observed_vector(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
        observer: [f64; 3],
        velocity: [f64; 3],
        surface: Option<[f64; 3]>,
    ) -> Result<[f64; 3]> {
        self.validate_source_alias(body)?;
        let state_at = |id, time| self.star_source_time_state(id, time).map(|state| state.0);
        let retarded_at = |time: TwoPartEpoch, delay| {
            let key = local_jpl_target(&JplLeg::new(self.registry), body)?;
            let fraction = time.lo
                + crate::time::source_tdb_minus_tt_seconds(time.hi, time.lo)
                    / reference_frames::FRAME_DAY_S;
            self.registry
                .eval_geometric_ssb_au_two_part_jd(key, time.hi, fraction - delay)
                .map(|state| state.0)
        };
        let ray = native_pheno_observe(
            &state_at,
            &retarded_at,
            epoch,
            body,
            observer,
            velocity,
            surface,
            false,
            true,
        )?;
        // These direct event/visibility rays bypass the public calc tail.
        // Publish only the successfully consumed SPK target, after light-time
        // and apparent-place corrections, never the configured provider mode.
        crate::trace::try_record(body, local_spk_trace_source(body))?;
        Ok(ray)
    }

    /// Resolve an orbital target without evaluating an epoch.
    /// Earth selects NAIF 3 EMB; ordinary targets prefer physical centers over
    /// barycenters. The selected reader enforces actual record coverage later.
    pub(crate) fn validate_orbital_target_alias(&self, body: i32) -> Result<()> {
        if body == EARTH {
            self.earth_barycenter_key().map(|_| ())
        } else {
            local_jpl_target(&JplLeg::new(self.registry), body).map(|_| ())
        }
    }

    /// Resolve a node/apse raw target without an epoch read. Earth is the
    /// geocenter NAIF 399 here, unlike the orbital-element EMB alias; ordinary
    /// planets retain the native center-first/barycenter lookup. Missing
    /// registrations and unsupported body identities remain typed.
    pub(crate) fn validate_nodaps_target_alias(&self, body: i32) -> Result<()> {
        local_jpl_target(&JplLeg::new(self.registry), body).map(|_| ())
    }

    /// Fresh scalar-TT Moon-minus-Earth state for native conic reduction.
    /// Returns true-ecliptic-of-date AU/AU-day using ordered raw constituents,
    /// independent TDB/frame clocks and whole-target guards. No public calc
    /// envelope or caller Delta-T is substituted; acquisition errors propagate.
    pub(crate) fn nodaps_true_moon_state(
        &self,
        tt: f64,
    ) -> Result<crate::lunar::EclipticMoonState> {
        local_jpl_true_moon_state(&JplLeg::new(self.registry), tt)
    }

    /// Scalar Time lunar state from an already acquired TT whole/fraction.
    /// Both words reach the raw SPK directory; frame observables are scalar.
    pub(crate) fn public_true_moon_state_parts(
        &self,
        whole: f64,
        fraction: f64,
    ) -> Result<crate::lunar::EclipticMoonState> {
        local_jpl_true_moon_state_parts(&JplLeg::new(self.registry), whole, fraction)
    }

    /// Acquire a node/apse PPN deflector at a finite TT Julian Day.
    ///
    /// Sun/Jupiter/Saturn select NAIF 10/5/6 independently of ordinary planet-center
    /// aliases. Return finite SSB ICRS AU. Missing aliases or uncovered optional
    /// segments return `None`; corrupt data and other typed errors propagate.
    /// No public calculation-envelope guard or alternate kernel is applied.
    pub(crate) fn nodaps_deflector_position(&self, body: i32, tt: f64) -> Result<Option<[f64; 3]>> {
        if !tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: tt });
        }
        let target = match body {
            SUN => 10,
            JUPITER => 5,
            SATURN => 6,
            _ => return Err(Error::UnknownBody { body_id: body }),
        };
        let Some(key) = self
            .registry
            .registered_bodies()
            .into_iter()
            .find(|entry| entry.naif_id == target)
            .map(|entry| entry.body)
        else {
            return Ok(None);
        };
        let whole = tt.floor();
        let fraction = tt - whole;
        let tdb_fraction =
            fraction + reference_frames::tdb_minus_tt_seconds(whole, fraction)? / 86400.0;
        let (position, velocity) =
            match self
                .registry
                .eval_geometric_ssb_au_two_part_jd(key, whole, tdb_fraction)
            {
                Ok(state) => state,
                Err(Error::SpkOutOfRange { .. }) => return Ok(None),
                Err(other) => return Err(other),
            };
        check_finite_vec3(position)?;
        check_finite_vec3(velocity)?;
        Ok(Some(position))
    }

    /// Resolve the source `earth barycenter` alias to an immutable registry
    /// entry carrying NAIF 3. The library body's numeric key is not assumed to
    /// equal its NAIF target. A missing target is a typed SPK target failure;
    /// no Earth/Moon reconstruction or center substitute is admitted.
    fn earth_barycenter_key(&self) -> Result<i32> {
        self.registry
            .registered_bodies()
            .into_iter()
            .find(|entry| entry.naif_id == 3)
            .map(|entry| entry.body)
            .ok_or(Error::SpkBodyNotFound { target: 3 })
    }

    /// Earth-Moon barycenter SSB state from a finite split TT epoch.
    ///
    /// Returns ICRS AU and AU/day with `IcrsSystemBarycenter` identity. The
    /// TT-to-TDB correction retains both epoch parts; each segment converts units
    /// before chain composition. Missing registration, invalid segment type/frame,
    /// corruption and actual directory-coverage errors propagate. No apparent
    /// correction or public calculation-envelope promotion is applied.
    pub(crate) fn earth_barycenter_state_two_part(
        &self,
        epoch: TwoPartEpoch,
    ) -> Result<GeometricState> {
        self.earth_barycenter_state_with_range(epoch, false)
    }

    /// Direct EMB state for node/apse geometry, retaining an actual raw refusal.
    /// Units, split clocks and barycenter identity match the ordinary EMB feed;
    /// only this native vector channel attaches its owned failed descriptor.
    pub(crate) fn nodaps_earth_barycenter_state_two_part(
        &self,
        epoch: TwoPartEpoch,
    ) -> Result<GeometricState> {
        self.earth_barycenter_state_with_range(epoch, true)
    }

    /// Share checked TT/TDB and state construction while preserving the
    /// ordinary noncapturing EMB acquisition for all existing other callers.
    fn earth_barycenter_state_with_range(
        &self,
        epoch: TwoPartEpoch,
        retain_range: bool,
    ) -> Result<GeometricState> {
        check_finite_vec3([epoch.hi, epoch.lo, epoch.value()])?;
        let key = self.earth_barycenter_key()?;
        let tdb_fraction =
            epoch.lo + reference_frames::tdb_minus_tt_seconds(epoch.hi, epoch.lo)? / 86400.0;
        let (pos, vel) = if retain_range {
            check_finite_vec3([epoch.hi, tdb_fraction, epoch.hi + tdb_fraction])?;
            let mut captured = None;
            self.registry
                .eval_geometric_ssb_au_two_part_jd_with_range(
                    key,
                    epoch.hi,
                    tdb_fraction,
                    &mut captured,
                )
                .map_err(|error| match captured {
                    Some(raw) => error.with_backend_range(raw),
                    None => error,
                })?
        } else {
            self.registry
                .eval_geometric_ssb_au_two_part_jd(key, epoch.hi, tdb_fraction)?
        };
        check_finite_vec3(pos)?;
        check_finite_vec3(vel)?;
        Ok(GeometricState {
            body_id: EARTH,
            jd_tt: epoch.value(),
            pos,
            vel,
            frame: StateFrame::IcrsSystemBarycenter,
            planet_center: false,
            tier: None,
        })
    }

    /// Check the public scalar-TT envelope independently of split TDB reads.
    /// The body's selected main kernel provides the descriptor envelope;
    /// finite/registration/range failures retain the normal calc identity.
    pub(crate) fn validate_scalar_tt(&self, body: i32, tt: f64) -> Result<()> {
        validate_local_jpl_scalar_tt(&JplLeg::new(self.registry), tt, body)
    }

    /// Translate a derived SPK range miss to the public central-TT contract.
    /// The requested target's main file supplies the complete descriptor
    /// envelope; registration and other provider errors retain their identity.
    pub(crate) fn public_range_error(&self, body: i32, tt: f64, error: Error) -> Error {
        local_jpl_public_range_error(&JplLeg::new(self.registry), tt, body, error)
    }

    /// Acquire a planet-centric constituent at split TT, in ICRS AU/AU-day.
    /// The reader retains the physical alias, per-segment unit conversion and
    /// actual record-domain errors without a public scalar-JD envelope.
    /// The final optional vector is the position-polynomial dP/dTDB of the
    /// complete successfully evaluated chain. Original velocity remains separate,
    /// including Type 3's independent stored channels, in AU per TDB day.
    pub(crate) fn pctr_source_time_state(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<PctrStateWithTdbTangent> {
        let target = local_jpl_target(&JplLeg::new(self.registry), body)?;
        let correction = crate::time::source_tdb_minus_tt_seconds(epoch.hi, epoch.lo);
        let mut captured = None;
        self.registry
            .eval_geometric_ssb_au_two_part_jd_with_range_and_derivative(
                target,
                epoch.hi,
                epoch.lo + correction / 86400.0,
                &mut captured,
            )
            .map_err(|error| match captured {
                Some(raw) => error.with_backend_range(raw),
                None => error,
            })
    }

    /// Open offsets, in TDB days, owned by the actual PCTR SSB source chain.
    /// This metadata-only query evaluates no state or coefficient. Callers
    /// retain the original physical read/error before collecting a central
    /// domain and check supplementary epochs against it before another read.
    pub(crate) fn pctr_record_offsets_tdb(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(f64, f64)> {
        let target = local_jpl_target(&JplLeg::new(self.registry), body)?;
        self.registry
            .geometric_ssb_record_offsets(target, epoch.hi, epoch.lo)
    }

    /// Module-TT PCTR state and actual full-chain position-derivative status.
    ///
    /// Preserve the ordinary entry snapshot, checked split clock and native
    /// state construction. The separate optional position tangent certifies
    /// every evaluated link; it does not change the acquired AU or AU/TDB-day words. This route
    /// retains the module's original noncapturing range-error contract.
    pub(crate) fn pctr_state_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<(GeometricState, Option<[f64; 3]>)> {
        let entry = self.serving_entry_snapshot(body)?;
        let target = local_jpl_target(&JplLeg::new(self.registry), body)?;
        let correction = reference_frames::tdb_minus_tt_seconds(epoch.hi, epoch.lo)?;
        let (pos, vel, derivative) = self
            .registry
            .eval_geometric_ssb_au_two_part_jd_with_derivative(
                target,
                epoch.hi,
                epoch.lo + correction / 86400.0,
            )?;
        let state = build_state(body, epoch.value(), &entry, &pos, &vel, None, |_| None)?;
        Ok((state, derivative))
    }

    /// PCTR retarded target at its already acquired observation-TDB parts.
    /// The unchanged raw native directory evaluator consumes this explicit
    /// fraction directly; TT conversion is not applied twice. Returned vectors
    /// and any refused image/descriptor are independently owned scalar data.
    pub(crate) fn pctr_source_tdb_state(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<([f64; 3], [f64; 3])> {
        let target = local_jpl_target(&JplLeg::new(self.registry), body)?;
        let mut captured = None;
        self.registry
            .eval_geometric_ssb_au_two_part_jd_with_range(target, epoch.hi, epoch.lo, &mut captured)
            .map_err(|error| match captured {
                Some(raw) => error.with_backend_range(raw),
                None => error,
            })
    }

    /// Report a planetocentric SPK coverage failure with its actual segment metadata.
    /// The original UT1 JD identifies the request; its cause retains backend units.
    pub(crate) fn pctr_range_error(&self, body: i32, ut1: f64, error: Error) -> Error {
        let Some(raw) = error.backend_range().cloned() else {
            return error;
        };
        if !matches!(error.primary_error(), Error::SpkOutOfRange { .. }) {
            return error;
        }
        let Ok(key) = local_jpl_target(&JplLeg::new(self.registry), SUN) else {
            return error;
        };
        let Some(main) = self
            .registry
            .body_info(key)
            .and_then(|info| self.registry.kernel(&info.kernel))
        else {
            return error;
        };
        let (start, end) = main.index().all_segments().fold(
            (f64::INFINITY, f64::NEG_INFINITY),
            |span, segment| {
                (
                    span.0.min(segment.descriptor.start_et),
                    span.1.max(segment.descriptor.end_et),
                )
            },
        );
        Error::EphemerisRange {
            requested_jd: ut1,
            start_jd: FRAME_J2000 + start / 86400.0,
            end_jd: FRAME_J2000 + end / 86400.0,
            body: Some(body),
        }
        .with_cause(error)
        .with_backend_range(raw)
    }

    /// Read an explicit split TDB epoch as an SSB ICRS AU/AU-day state.
    /// Actual aliases, segment domains and per-segment unit conversion are
    /// retained. No second TT conversion or public scalar-JD guard is applied;
    /// non-finite, missing-registration and coverage errors propagate.
    pub(crate) fn state_tdb_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<([f64; 3], [f64; 3])> {
        let target = local_jpl_target(&JplLeg::new(self.registry), body)?;
        self.registry
            .eval_geometric_ssb_au_two_part_jd(target, epoch.hi, epoch.lo)
    }

    /// Classify the selected NAIF target without treating a retained system
    /// barycentre as a reconstructed planet centre.
    fn coord_type(&self, body: i32) -> Option<u32> {
        let key = local_jpl_target(&JplLeg::new(self.registry), body).ok()?;
        let target = self.registry.body_info(key)?.naif_id;
        Some(if (1..=9).contains(&target) {
            COORD_ICRS_BARY_SYSTEM
        } else {
            COORD_ICRS_BARY
        })
    }
}

impl StateProvider for LocalSpkStateProvider<'_> {
    fn provider_name(&self) -> &'static str {
        "native-jpl/local-spk"
    }

    fn has_body(&self, body_id: i32) -> bool {
        local_jpl_target(&JplLeg::new(self.registry), body_id).is_ok()
    }

    fn body_coverage(&self, body_id: i32) -> Option<(f64, f64)> {
        let key = local_jpl_target(&JplLeg::new(self.registry), body_id).ok()?;
        let key = self.registry.primary_body_for_coverage(key).ok()?;
        let spans = self.registry.coverage_jd_tdb(key).ok()?;
        let start = spans.first()?.0;
        let end = spans.last()?.1;
        Some((
            crate::time::tdb_to_tt_jd(start).ok()?,
            crate::time::tdb_to_tt_jd(end).ok()?,
        ))
    }

    fn serving_coord_type(&self, body_id: i32) -> Option<u32> {
        self.coord_type(body_id)
    }

    fn eval_raw(&self, body_id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        let day = jd.floor();
        let (pos, vel) = local_jpl_state(&JplLeg::new(self.registry), body_id, day, jd - day)?;
        Ok((pos.to_vec(), vel.to_vec()))
    }

    fn eval_raw_two_part(&self, body_id: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        let (pos, vel) = local_jpl_state(&JplLeg::new(self.registry), body_id, epoch.hi, epoch.lo)?;
        Ok((pos.to_vec(), vel.to_vec()))
    }

    fn serving_entry_snapshot(&self, body_id: i32) -> Result<EntrySnapshot> {
        // Preserve UnknownBody for unsupported ids, distinct from a missing
        // registration for an admitted core target.
        local_jpl_target(&JplLeg::new(self.registry), body_id)?;
        let coord_type = self
            .coord_type(body_id)
            .ok_or(Error::SpkUnregisteredBody { body: body_id })?;
        Ok(EntrySnapshot { coord_type })
    }
}

#[cfg(test)]
#[path = "tests/local_jpl_contract.rs"]
mod local_jpl_contract_tests;

impl<'a> JplLeg<'a> {
    /// Bind the context-local kernel store (never shared across
    /// contexts; SPK registry-ownership layer).
    pub fn new(registry: &'a SpkRegistry) -> Self {
        Self { registry }
    }

    /// The bound kernel store.
    pub fn registry(&self) -> &'a SpkRegistry {
        self.registry
    }

    /// Whether the registry contains any local SPK image.
    ///
    /// This is an availability hint, not a per-body coverage test. Actual target
    /// support and stored segment domains are checked during evaluation.
    pub fn has_local_source(&self) -> bool {
        !self.registry.is_empty()
    }

    /// Whether the registry's merged stored window for `body_id` covers
    /// `jd_tdb` (closed interval, TDB Julian Days). Unknown bodies and
    /// non-finite epochs simply do not cover — this is a routing
    /// predicate, not the evaluator's typed error.
    pub fn covers_tdb(&self, body_id: i32, jd_tdb: f64) -> bool {
        if !jd_tdb.is_finite() {
            return false;
        }
        match self.registry.coverage_jd_tdb(body_id) {
            Ok(spans) => spans
                .iter()
                .any(|&(start, end)| start <= jd_tdb && jd_tdb <= end),
            Err(_) => false,
        }
    }

    /// Geometric SSB-centered state of library body `body_id` at
    /// ephemeris seconds past J2000 (`et`; TDB-based TT/TDB convention), in km
    /// and km/s (chain-composition [`eval_ssb`](crate::spk::eval_ssb) over the
    /// registry body map, centers resolved recursively from the same
    /// registry). Unregistered bodies, coverage misses, frame/type
    /// refusals, and corrupt records propagate unchanged — never
    /// invented data, never a silent fallback.
    pub fn eval_geometric_ssb(&self, body_id: i32, et_seconds: f64) -> Result<ChebState> {
        let provider = NaifBaseProvider::new(self.registry);
        self.registry.eval_body_ssb(body_id, &provider, et_seconds)
    }
}

#[cfg(test)]
#[path = "tests/native_pheno_vectors.rs"]
mod native_pheno_vector_contract_tests;

#[cfg(test)]
#[path = "tests/native_pheno_ieee.rs"]
mod native_pheno_review_observed_ieee_tests;

#[cfg(test)]
#[path = "tests/source_time_domains.rs"]
mod source_time_domain_tests;

#[cfg(test)]
#[path = "tests/minor_model_sources.rs"]
mod minor_model_j2000_source_fixture_tests;

#[cfg(test)]
#[path = "tests/astropy_time.rs"]
mod astropy_core_time_tests;

/// Return the main DE image's owned registration without evaluating a state.
///
/// Selection reuses the NativeJpl Sun target resolver: library body keys are
/// independent of NAIF ids, and the first NAIF10 registration in body-key
/// order wins. No kernel image is copied, acquired or replaced. Missing Sun
/// registration is a typed refusal; another same-span image is never selected.
pub fn native_main_image_registration(
    registry: &SpkRegistry,
) -> Result<crate::spk::BodyRegistration> {
    let key = local_jpl_target(&JplLeg::new(registry), SUN)?;
    registry
        .body_info(key)
        .ok_or(Error::SpkUnregisteredBody { body: SUN })
}

#[cfg(test)]
#[path = "tests/main_image_registration.rs"]
mod main_image_registration_tests;

#[cfg(test)]
#[path = "tests/public_model_ranges.rs"]
mod public_model_backend_range_tests;

#[cfg(test)]
#[path = "tests/native_minor_provider_priority.rs"]
mod native_minor_provider_priority_tests;

#[cfg(test)]
#[path = "routing_derivative_tests.rs"]
mod routing_derivative_tests;

#[cfg(test)]
mod ppn_reception_origin_tests {
    use super::*;
    use std::cell::RefCell;

    #[test]
    fn native_observation_keeps_actual_sun_split_and_retarded_earth() {
        // Index-only synthetic registrations; the injected state reader below
        // supplies all positions. No stored astronomical coefficients are read.
        let mut image = crate::spk::moon_test_kernel_bytes();
        let first_target = crate::spk::RECORD_BYTES + 24 + 16;
        image[first_target..first_target + 4].copy_from_slice(&10_i32.to_le_bytes());
        image[first_target + 40..first_target + 44].copy_from_slice(&399_i32.to_le_bytes());
        let mut registry = SpkRegistry::new();
        registry.register_image("synthetic-origin", image).unwrap();
        registry.register_body(SUN, "synthetic-origin", 10).unwrap();
        registry
            .register_body(EARTH, "synthetic-origin", 399)
            .unwrap();
        let leg = JplLeg::new(&registry);
        let tt = 2451545.125;
        let reads = RefCell::new(Vec::new());
        let state = |body, day, fraction| {
            let mut reads = reads.borrow_mut();
            let position = match body {
                SUN => {
                    assert!(
                        reads.is_empty(),
                        "the origin is read only once, before Earth"
                    );
                    [8.0, -4.0, 2.0]
                }
                EARTH if reads.len() == 1 => [9.0, -4.0, 2.0],
                EARTH => [10.0, -4.0, 2.0],
                _ => panic!("unexpected support"),
            };
            reads.push((body, day, fraction));
            Ok((position, [0.0; 3]))
        };
        let observation = retarded_earth_observation_source_time_via(&leg, tt, &state).unwrap();
        assert_eq!(observation.sun_reception, [8.0, -4.0, 2.0]);
        assert_eq!(observation.relative, [2.0, 0.0, 0.0]);
        let complete = reads.take();
        assert_eq!(
            complete.iter().map(|read| read.0).collect::<Vec<_>>(),
            [SUN, EARTH, EARTH, EARTH]
        );
        assert!(complete.iter().all(|read| read.1 == 2451545.0));
        assert_eq!(complete[0].2, complete[1].2);
        assert_ne!(complete[0].2, 0.125); // The origin is the actual converted TDB sample.
        assert!(complete[2].2 < complete[1].2 && complete[3].2 < complete[2].2);
        assert_eq!(
            retarded_earth_minus_sun_source_time_via(&leg, tt, &state).unwrap(),
            observation.relative
        );
        assert_eq!(*reads.borrow(), complete); // Other model wrappers retain their exact reads.
    }
}
