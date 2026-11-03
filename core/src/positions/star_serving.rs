//! Provider-aware public fixed-star serving over immutable context images.
//!
//! The stored LEB pipeline and the native SPK vector pipeline keep distinct
//! acquisition channels. Catalog resolution, public retflags, coordinate and
//! sidereal shaping are shared. Native scalar calls create fresh scalar-TT
//! Times; native geocentric batches retain the original bundled UT1 Time for
//! their shared observer. Empty NativeJpl registries acquire their configured main image at consumed reads;
//! explicit LEB never substitutes JPL.

use super::reference_frames;
#[cfg(test)]
use super::routing::LocalSpkStateProvider;
use super::routing::native_apparent_ray_with_deflectors;
use super::star_main::{NativeStarMainCache, NativeStarProvider};
use super::{EntrySnapshot, GeometricState, StateProvider};
use crate::constants::{
    EARTH, FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NOABERR, FLG_NOGDEFL, FLG_SIDEREAL, FLG_SPEED,
    FLG_TOPOCTR, FLG_TRUEPOS, SIDBIT_ECL_T0, SIDBIT_SSY_PLANE, SUN,
};
use crate::context::{EphemerisContext, ProviderMode};
use crate::errors::{Error, Result};
use crate::stars::{
    AcquiredStarRow, SealedStarServing, StarData, StarMatch, StarResolveError, StarRowAcquisition,
    StarRowBackend, StarSample, StarTopoSite, WrappedStarEntry, fixed_star_data,
    fixstar_observation_center, fixstar_ret_flags, preprocess_star_flags, resolve_star_ref,
    star_geometric_range_state, star_velocity_samples_via,
};
use crate::time::{DeltaTModel, deltat_with_config};
use crate::types::TwoPartEpoch;
use std::sync::Arc;

/// Shared borrowed state interface for the public star flag/anchor layer.
///
/// `inner` is the actual date-aware LEB provider or native SPK provider. Every
/// raw/geometric operation delegates unchanged, including stored metadata,
/// coverage, split epochs, units and typed errors. This adapter contains no
/// state arithmetic and cannot outlive the selected
/// immutable provider. Its NativeJpl view borrows a lazy owned main selector.
/// Native catalog places and dynamic anchors use separate
/// split-epoch ray callbacks; the stored LEB channel stays unchanged.
pub struct ContextStarProvider<'a> {
    inner: ContextStarBackend<'a>,
}

/// Selected immutable provider borrow. No provider-wide generic method is
/// erased to a trait object; the existing StateProvider contract stays intact.
enum ContextStarBackend<'a> {
    Leb(&'a super::context_serving::BoundProvider<'a>),
    #[cfg(test)]
    Native(&'a LocalSpkStateProvider<'a>),
    AcquiringNative(&'a NativeStarProvider<'a>),
}

/// Forward an operation to the actual selected provider without arithmetic.
macro_rules! selected_provider {
    ($this:expr, $provider:ident, $operation:expr) => {
        match &$this.inner {
            ContextStarBackend::Leb($provider) => $operation,
            #[cfg(test)]
            ContextStarBackend::Native($provider) => $operation,
            ContextStarBackend::AcquiringNative($provider) => $operation,
        }
    };
}

impl StateProvider for ContextStarProvider<'_> {
    fn has_analytic_kinematics(&self) -> bool {
        selected_provider!(self, provider, provider.has_analytic_kinematics())
    }
    fn geometric_acceleration_tt(&self, id: i32, jd: f64) -> Result<Option<[f64; 3]>> {
        selected_provider!(self, provider, provider.geometric_acceleration_tt(id, jd))
    }
    fn geometric_state_word_parts_two_part(
        &self,
        id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<super::CartesianWordParts> {
        selected_provider!(
            self,
            provider,
            provider.geometric_state_word_parts_two_part(id, epoch)
        )
    }

    fn provider_name(&self) -> &'static str {
        selected_provider!(self, provider, provider.provider_name())
    }
    fn has_body(&self, id: i32) -> bool {
        selected_provider!(self, provider, provider.has_body(id))
    }
    fn body_coverage(&self, id: i32) -> Option<(f64, f64)> {
        selected_provider!(self, provider, provider.body_coverage(id))
    }
    fn smooth_interval_tt(&self, id: i32, jd: f64) -> Result<Option<(f64, f64)>> {
        selected_provider!(self, provider, provider.smooth_interval_tt(id, jd))
    }
    fn nutation_smooth_interval_tt(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        selected_provider!(self, provider, provider.nutation_smooth_interval_tt(jd))
    }
    fn serving_coord_type(&self, id: i32) -> Option<u32> {
        selected_provider!(self, provider, provider.serving_coord_type(id))
    }
    fn eval_raw(&self, id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        selected_provider!(self, provider, provider.eval_raw(id, jd))
    }
    fn eval_raw_two_part(&self, id: i32, time: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        selected_provider!(self, provider, provider.eval_raw_two_part(id, time))
    }
    fn selected_tier(&self, id: i32, jd: f64) -> Option<crate::leb::LebTier> {
        selected_provider!(self, provider, provider.selected_tier(id, jd))
    }
    fn geometric_state(&self, id: i32, jd: f64) -> Result<GeometricState> {
        selected_provider!(self, provider, provider.geometric_state(id, jd))
    }
    fn geometric_state_two_part(&self, id: i32, time: TwoPartEpoch) -> Result<GeometricState> {
        selected_provider!(self, provider, provider.geometric_state_two_part(id, time))
    }
    fn serving_entry_snapshot(&self, id: i32) -> Result<EntrySnapshot> {
        selected_provider!(self, provider, provider.serving_entry_snapshot(id))
    }
}

impl ContextStarProvider<'_> {
    /// Preserve the selected stored reader's existing vector-range graph.
    /// This mapper is bound only to live stored-backed vector reductions.
    fn stored_vector_error(&self, error: Error) -> Error {
        match &self.inner {
            ContextStarBackend::Leb(index) => index.vector_error(error),
            #[cfg(test)]
            ContextStarBackend::Native(_) => error,
            ContextStarBackend::AcquiringNative(_) => error,
        }
    }
}

/// Cache-first live stored anchors with acquisition only at a consumed miss.
/// This factory retains configuration and row owners without assembling a
/// reader, checking file presence, evaluating Delta-T or reading a catalog.
/// A miss has already acquired six fields and constructed initial vectors;
/// its helper role determines whether scalar TT Time precedes or follows the
/// stored vector-stack assembly. Every vector query consumes scalar epoch.value
/// and keeps the actual reader error graph. Static/fixed-star routes are
/// unchanged. Auto cache hits need no image; Auto cold JPL selection is a
/// separate caller policy and is not manufactured by this stored-only factory.
pub(crate) fn with_context_live_stored_anchor_provider<R>(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    rows: &dyn crate::stars::AnchorRowAcquisition,
    run: impl FnOnce(&dyn crate::sidereal::DynamicAnchorProvider) -> Result<R>,
) -> Result<R> {
    struct Unconsumed;
    impl StateProvider for Unconsumed {
        fn provider_name(&self) -> &'static str {
            "unconsumed-live-anchor"
        }
        fn has_body(&self, _: i32) -> bool {
            false
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            None
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            None
        }
        fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            Err(unconsumed())
        }
        fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            Err(unconsumed())
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<super::EntrySnapshot> {
            Err(unconsumed())
        }
    }
    fn unconsumed() -> Error {
        Error::ConfigurationMissing {
            missing_config: "live_anchor_acquired_ray",
            suggestion: "live stored anchors must use their consumed ray callback",
        }
    }
    let acquired = |initial, tt, noaberr, nogdefl, role| {
        let time = if role == crate::stars::AnchorRayRole::GalacticRightAscension {
            Some(scalar_tt(tt)?)
        } else {
            None
        };
        super::context_serving::with_sealed_serving_stack(ctx, model, iers, true, |stack| {
            let provider = ContextStarProvider {
                inner: ContextStarBackend::Leb(stack.serving.provider()),
            };
            let map = |error| provider.stored_vector_error(error);
            let state = |body, time| live_stored_anchor_state_via(&provider, body, time, &map);
            let time = match time {
                Some(time) => time,
                None => scalar_tt(tt)?,
            };
            live_anchor_ray_at_time_via(&state, initial, time, noaberr, nogdefl)
        })
    };
    let range = |_, error| error;
    let unused = Unconsumed;
    let anchors = crate::stars::SealedAnchorProvider::new(&unused)
        .with_live_rows(rows)
        .with_acquired_ray(&acquired, &range);
    run(&anchors)
}

/// Borrow the actual independent computation provider for live nested anchors.
///
/// Forced LEB uses its lazy scalar-Time vector adapter. Auto selects the
/// Native main computation view even when the fixed-star body was served by
/// LEB; the configuration clone opens no image. Native main acquisition stays
/// inside a consumed uncached ray after its row fields and constructor. A
/// completed anchor cache hit therefore needs no main or reader image. The
/// synchronous callback cannot retain providers, and genuine terminal Native
/// errors propagate from the provider actually used for the anchor.
#[allow(clippy::too_many_arguments)]
pub(crate) fn with_context_live_anchor_provider_captured<R>(
    ctx: &EphemerisContext,
    cache: &NativeStarMainCache,
    resolve: Option<&StarPositionContextResolver<'_>>,
    model: DeltaTModel,
    iers: bool,
    rows: &dyn crate::stars::AnchorRowAcquisition,
    run: impl FnOnce(&dyn crate::sidereal::DynamicAnchorProvider) -> Result<R>,
) -> CapturedStarResult<R> {
    if ctx.mode() == ProviderMode::SealedLeb {
        return CapturedStarResult {
            result: with_context_live_stored_anchor_provider(ctx, model, iers, rows, run),
        };
    }
    let native_context =
        (ctx.mode() == ProviderMode::Auto).then(|| ctx.clone().with_mode(ProviderMode::NativeJpl));
    with_context_star_serving_captured_with_main(
        native_context.as_ref().unwrap_or(ctx),
        cache,
        resolve,
        model,
        iers,
        0,
        None,
        |serving, _| {
            let anchors = serving.anchor_provider().with_live_rows(rows);
            run(&anchors)
        },
    )
}

/// Run a synchronous stellar calculation over the context's actual provider.
///
/// LEB and Auto select covered stored states and nutation; closed images count
/// as absent. Native JPL retains the immutable registry or a lazily acquired
/// main image, with strict planetary deflectors and the bundled WGS84 observer
/// clock. Remote Horizons stellar serving is unsupported. `site` is east
/// longitude/latitude in degrees and height in metres, used only for TOPOCTR.
/// A required missing site is an error. Output includes six position/rate
/// components, catalog name and returned flags.
///
/// Configured Delta-T/IERS and user/tidal settings remain available for
/// sidereal defining epochs. Provider owners and borrows stay within this
/// call; image acquisition respects the effective network gate. Missing aliases,
/// lookup, coverage, corrupt-state, clock and frame errors propagate unchanged.
pub fn with_context_star_serving<R>(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    run: impl FnOnce(&SealedStarServing<'_, ContextStarProvider<'_>>, Option<StarTopoSite>) -> Result<R>,
) -> Result<R> {
    with_context_star_serving_captured(ctx, model, iers, flags, site, run).result
}

/// Owned result of a synchronous stellar calculation.
/// The native error owns its data and no provider/context borrow escapes.
pub struct CapturedStarResult<R> {
    /// Completed native calculation or its typed failure.
    pub result: Result<R>,
}

/// Run one synchronous stellar calculation with context-owned resources.
/// Stored and SPK providers retain their normal data, clock and range checks.
pub fn with_context_star_serving_captured<R>(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    run: impl FnOnce(&SealedStarServing<'_, ContextStarProvider<'_>>, Option<StarTopoSite>) -> Result<R>,
) -> CapturedStarResult<R> {
    let cache = NativeStarMainCache::default();
    with_context_star_serving_captured_with_main(ctx, &cache, None, model, iers, flags, site, run)
}

/// Borrowed configuration-only main-image resolver at a consumed stellar read.
/// `body` is Sun, identifying the main image. Return an owned NativeJpl context
/// without opening data or replacing supplied registry images. The callback is
/// synchronous, runs without the image-owner lock held, and may fail normally.
pub type StarPositionContextResolver<'a> =
    dyn Fn(&EphemerisContext, i32) -> Result<EphemerisContext> + 'a;

/// The actual public observer-clock stage. Scalar samples construct TT before
/// Main. Shared batches retain their original UT1 and prepare/select side Times
/// separately from resource acquisition and the autonomous numeric clocks.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StarObservationClock {
    /// Scalar TT Julian day, including each independent speed-side sample.
    Tt(f64),
    /// Original shared UT1 Julian day, before observer/row evaluation.
    Ut1(f64),
    /// Construct both shared speed-side Times, previous then next, after the
    /// central observer and before either side observer. No new Main is read.
    SharedSpeedTimes,
    /// Select an already constructed shared Time for an observer or row read.
    /// Selection performs no clock construction or Main acquisition.
    SharedSample(StarObservationSample),
}

/// A retained shared-batch observer and its actual Time owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum StarObservationSample {
    /// Original UT1 Time used by the central shared observer.
    Central,
    /// Scalar TT Time half a day before the original epoch's TT observable.
    Previous,
    /// Scalar TT Time half a day after the original epoch's TT observable.
    Next,
}

/// Acquire the actual Time and optional Main at their consumed stages. A
/// returned context owns the selected immutable registry. None retains the
/// independent receiver's existing native Main/cache. Neither choice replaces
/// model/frame clocks or caller flags. The synchronous callback runs without a
/// provider lock, including warm-cache calls; shared selection reuses its owner.
pub type StarObservationAcquisition<'a> =
    dyn Fn(StarObservationClock) -> Result<Option<Arc<EphemerisContext>>> + Send + Sync + 'a;

/// Public star serving with lazy binding configuration and an owned main image.
/// Inputs, clocks, units and outputs follow [`with_context_star_serving`]. The
/// resolver is used only for an empty NativeJpl registry, at the actual scalar,
/// shared batch or dynamic-anchor acquisition stage. LEB/Auto do not consume it.
/// It supplies configuration; the ordinary core main loader owns file/network
/// acquisition. The call retains its success snapshot through every state read.
#[allow(clippy::too_many_arguments)]
pub fn with_context_star_serving_with_position_context<R>(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    resolve: &StarPositionContextResolver<'_>,
    run: impl FnOnce(&SealedStarServing<'_, ContextStarProvider<'_>>, Option<StarTopoSite>) -> Result<R>,
) -> Result<R> {
    let cache = NativeStarMainCache::default();
    with_context_star_serving_captured_with_main(
        ctx,
        &cache,
        Some(resolve),
        model,
        iers,
        flags,
        site,
        run,
    )
    .result
}

/// Retain actual failure ownership over a Session's success-only main snapshot.
/// The optional configuration callback is borrowed only through this call. It
/// executes at consumed acquisition, after the scalar clock or before the anchor
/// clock; completed values and errors remain native owned results.
#[allow(clippy::too_many_arguments)]
pub(crate) fn with_context_star_serving_captured_with_main<R>(
    ctx: &EphemerisContext,
    cache: &NativeStarMainCache,
    resolve: Option<&StarPositionContextResolver<'_>>,
    model: DeltaTModel,
    iers: bool,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    run: impl FnOnce(&SealedStarServing<'_, ContextStarProvider<'_>>, Option<StarTopoSite>) -> Result<R>,
) -> CapturedStarResult<R> {
    with_context_star_serving_captured_with_rows(
        ctx, cache, resolve, model, iers, flags, site, None, run,
    )
}

/// Add live row leases to the existing context-owned star acquisition.
/// No row is acquired at construction; scalar and shared observer callbacks
/// consume rows at the provider's physical sample epochs.
#[allow(clippy::too_many_arguments)]
pub(crate) fn with_context_star_serving_captured_with_rows<R>(
    ctx: &EphemerisContext,
    cache: &NativeStarMainCache,
    resolve: Option<&StarPositionContextResolver<'_>>,
    model: DeltaTModel,
    iers: bool,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    rows: Option<&dyn StarRowAcquisition>,
    run: impl FnOnce(&SealedStarServing<'_, ContextStarProvider<'_>>, Option<StarTopoSite>) -> Result<R>,
) -> CapturedStarResult<R> {
    with_context_star_serving_captured_with_rows_and_acquisition(
        ctx, cache, resolve, model, iers, flags, site, rows, None, run,
    )
}

/// Live acquired-row variant retaining the exact observer-stage Main producer.
/// Absent acquisition preserves the original provider/cache behavior; stored
/// channels never call the native producer.
#[allow(clippy::too_many_arguments)]
pub(crate) fn with_context_star_serving_captured_with_rows_and_acquisition<R>(
    ctx: &EphemerisContext,
    cache: &NativeStarMainCache,
    resolve: Option<&StarPositionContextResolver<'_>>,
    model: DeltaTModel,
    iers: bool,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    rows: Option<&dyn StarRowAcquisition>,
    acquisition: Option<&StarObservationAcquisition<'_>>,
    run: impl FnOnce(&SealedStarServing<'_, ContextStarProvider<'_>>, Option<StarTopoSite>) -> Result<R>,
) -> CapturedStarResult<R> {
    CapturedStarResult {
        result: with_context_star_serving_impl(
            ctx,
            cache,
            resolve,
            model,
            iers,
            flags,
            site,
            rows,
            acquisition,
            run,
        ),
    }
}

/// Resolve Auto stellar resources from actual support at the consumed TT epoch.
/// Returns a NativeJpl context only when the stored observer/frame samples
/// are unavailable, or when Auto requests the vector topocentric provider.
/// Forced LEB is never changed. Metadata corruption propagates; no numerical
/// failure is caught or converted into provider availability. The immutable
/// images/registrations are shared with the original context.
pub(crate) fn auto_star_context(
    context: &EphemerisContext,
    epoch: f64,
    is_ut: bool,
    flags: i32,
) -> Result<Option<EphemerisContext>> {
    if context.mode() != ProviderMode::Auto {
        return Ok(None);
    }
    let tt = if is_ut {
        let (whole, fraction) = crate::smh2016::timescale_parts_from_ut1(epoch)?;
        whole + fraction
    } else {
        scalar_tt(epoch)?.value()
    };
    let native = || Some(context.clone().with_mode(ProviderMode::NativeJpl));
    let flags = preprocess_star_flags(flags);
    if flags & FLG_TOPOCTR != 0 {
        return Ok(native());
    }
    let files: Vec<_> = context
        .leb_files()
        .iter()
        .filter(|file| !file.is_closed())
        .collect();
    if files.is_empty() {
        return Ok(native());
    }
    let provider = super::context_serving::provider_from_images(
        files.iter().map(|file| (file.path(), *file)),
        context.tier(),
        true,
    )?;
    let center = fixstar_observation_center(flags);
    let observer = match center {
        FLG_BARYCTR => None,
        FLG_HELCTR => Some(SUN),
        _ => Some(EARTH),
    };
    let mut epochs = vec![tt];
    if flags & FLG_SPEED != 0 {
        // The raw angular samples and the flag-frame samples both consume
        // half-day epochs. Check each actual polynomial, including seams.
        epochs.extend([tt - 0.5, tt + 0.5]);
    }
    for sample in epochs {
        if let Some(body) = observer {
            if provider.smooth_interval_tt(body, sample)?.is_none() {
                return Ok(native());
            }
        }
        if flags & FLG_J2000 == 0 && provider.nutation_smooth_interval_tt(sample)?.is_none() {
            return Ok(native());
        }
    }
    Ok(None)
}

/// Borrow the selected provider synchronously for a native stellar calculation.
#[allow(clippy::too_many_arguments)]
fn with_context_star_serving_impl<R>(
    ctx: &EphemerisContext,
    cache: &NativeStarMainCache,
    resolve: Option<&StarPositionContextResolver<'_>>,
    model: DeltaTModel,
    iers: bool,
    flags: i32,
    site: Option<(f64, f64, f64)>,
    rows: Option<&dyn StarRowAcquisition>,
    acquisition: Option<&StarObservationAcquisition<'_>>,
    run: impl FnOnce(&SealedStarServing<'_, ContextStarProvider<'_>>, Option<StarTopoSite>) -> Result<R>,
) -> Result<R> {
    if ctx.mode() != ProviderMode::NativeJpl {
        return super::context_serving::with_sealed_serving_stack(
            ctx,
            model,
            iers,
            true,
            |stack| {
                let provider = ContextStarProvider {
                    inner: ContextStarBackend::Leb(stack.serving.provider()),
                };
                let delta = |ut| deltat_with_config(ut, stack.delta_cfg);
                let (selection, _) = crate::houses::sidereal_selection_for_context(ctx);
                let vector_error = |error| provider.stored_vector_error(error);
                let anchor_state =
                    |body, time| live_stored_anchor_state_via(&provider, body, time, &vector_error);
                let acquired_anchor_at = |initial, tt, noaberr, nogdefl, role| {
                    // The immutable LEB vector adapter is already retained by this
                    // stack; no JPL main acquisition belongs to either helper.
                    live_anchor_ray_via(&anchor_state, initial, tt, noaberr, nogdefl, role, || {
                        Ok(())
                    })
                };
                let anchor_range = |_, error| error;
                let mut serving =
                    SealedStarServing::new(&provider, stack.nutation_at, &delta, selection)
                        .with_default_timescale()
                        .with_stored_vector_error(&vector_error)
                        .with_acquired_anchors(&acquired_anchor_at, &anchor_range);
                if let Some(rows) = rows {
                    serving = serving.with_live_rows(rows);
                }
                let topo = public_site(flags, site)?;
                run(&serving, topo)
            },
        );
    }
    let native = NativeStarProvider::new(ctx, cache, resolve);
    let acquire_at = |clock| {
        if let Some(selected) = acquisition
            .map(|acquire| acquire(clock))
            .transpose()?
            .flatten()
        {
            native.select_acquired(selected)?;
        }
        native.ensure_main()
    };
    let provider = ContextStarProvider {
        inner: ContextStarBackend::AcquiringNative(&native),
    };
    let config = ctx.delta_t_config(model, iers);
    let delta = |ut| deltat_with_config(ut, &config);
    let nutation = |tt| crate::erfa::nut06a(tt);
    let (selection, _) = crate::houses::sidereal_selection_for_context(ctx);
    let state_at = |id, time| native.star_source_time_state(id, time);
    let range_at = |center, tt, error| native.star_range_error(center, tt, error);
    let anchor_at = |data: &StarData, tt, noaberr, nogdefl| {
        native.ensure_main()?;
        native_anchor_ray_via(&state_at, data, tt, noaberr, nogdefl)
    };
    let acquired_anchor_at = |initial, tt, noaberr, nogdefl, role| {
        live_anchor_ray_via(&state_at, initial, tt, noaberr, nogdefl, role, || {
            native.ensure_main()
        })
    };
    let anchor_range = |tt, error| range_at(0, tt, error);
    let raw_at = |id, tt, noaberr, nogdefl, j2000, center, topo, speed| {
        if let Some(rows) = rows {
            if speed {
                rows.velocity_admission(id)?;
            }
        } else if fixed_star_data(id).is_some() {
            scalar_tt(tt)?;
            native.ensure_main()?;
        }
        if rows.is_some() {
            native_raw_place_with_rows_and_main_via(
                &state_at,
                &range_at,
                id,
                tt,
                noaberr,
                nogdefl,
                j2000,
                center,
                topo,
                speed,
                rows,
                Some(&|tt| acquire_at(StarObservationClock::Tt(tt))),
            )
        } else {
            native_raw_place_with_rows_via(
                &state_at, &range_at, id, tt, noaberr, nogdefl, j2000, center, topo, speed, rows,
            )
        }
    };
    // The delegate has the scalar callback but no batch callback, preventing
    // recursion when the source delegates topo/fixed-epoch/SIDBIT batches.
    let mut delegate = SealedStarServing::new(&provider, &nutation, &delta, selection)
        .with_default_timescale()
        .with_native_places(&raw_at, None)
        .with_native_anchors(&anchor_at, &anchor_range)
        .with_acquired_anchors(&acquired_anchor_at, &anchor_range);
    if let Some(rows) = rows {
        delegate = delegate.with_live_rows(rows);
    }
    let batch_at = |names: &[&str], ut, flags, topo, skip| {
        native_batch_with_main_via(
            &|| native.ensure_main(),
            &state_at,
            &range_at,
            &delegate,
            selection,
            names,
            ut,
            flags,
            topo,
            skip,
        )
    };
    let resolved_batch_at =
        |resolved: &[Option<StarMatch>],
         row_owners: &[Option<std::sync::Arc<dyn StarRowAcquisition>>],
         ut,
         flags,
         _topo,
         skip,
         anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>| {
            let live = anchors.map(|anchors| delegate.borrowing_dynamic_anchors(anchors));
            let delegate = live.as_ref().unwrap_or(&delegate);
            native_resolved_batch_with_main_via(
                &|| acquire_at(StarObservationClock::Ut1(ut)),
                acquisition.map(|_| &acquire_at as &dyn Fn(StarObservationClock) -> Result<()>),
                &state_at,
                &range_at,
                delegate,
                resolved,
                row_owners,
                ut,
                flags,
                skip,
            )
        };
    let mut serving = SealedStarServing::new(&provider, &nutation, &delta, selection)
        .with_default_timescale()
        .with_native_places(&raw_at, Some(&batch_at))
        .with_native_anchors(&anchor_at, &anchor_range)
        .with_acquired_anchors(&acquired_anchor_at, &anchor_range);
    serving = serving.with_resolved_native_anchor_batch(&resolved_batch_at);
    if let Some(rows) = rows {
        serving = serving.with_live_rows(rows);
    }
    run(&serving, public_site(flags, site)?)
}

/// Admit a per-call WGS84 site only when the raw request selects TOPOCTR.
/// Other calls leave even unused site range values unconsumed. The binding
/// separately checks finite scalar inputs on every call; StarTopoSite checks
/// finite consumed values here. Coordinates keep their trigonometric meaning.
fn public_site(flags: i32, site: Option<(f64, f64, f64)>) -> Result<Option<StarTopoSite>> {
    if flags & FLG_TOPOCTR == 0 {
        return Ok(None);
    }
    let site = site.ok_or(Error::ConfigurationMissing {
        missing_config: "observer_location",
        suggestion: "topocentric star places need an explicit site",
    })?;
    StarTopoSite::new(site.0, site.1, site.2).map(Some)
}

/// Star scalar/shared-observer range promotion callback: observing center,
/// failed sample scalar TT and original typed error to the public error.
/// Only raw SPK range failures are promoted; other errors stay unchanged.
type StarRangeAt<'a> = dyn Fn(i32, f64, Error) -> Error + 'a;

/// Native split-epoch SSB state callback in ICRS AU and AU/day.
/// The callback owns alias/whole-chain/segment/frame/type/domain admission;
/// it preserves both TT parts and per-segment conversion before addition.
type StarStateAt<'a> = dyn Fn(i32, TwoPartEpoch) -> Result<([f64; 3], [f64; 3])> + 'a;

/// A source observer frozen at one epoch: SSB ICRS position AU, velocity
/// AU/day and optional geocenter-to-observer GCRS AU surface vector. It is
/// owned scalar state, not a borrowed provider or cached mutable observer.
#[derive(Clone, Copy)]
pub(super) struct NativeStarObserver {
    pub(super) time: TwoPartEpoch,
    pub(super) position: [f64; 3],
    pub(super) velocity: [f64; 3],
    /// Observer-position derivative in AU per TT day.
    pub(super) velocity_per_tt: [f64; 3],
    pub(super) surface: Option<[f64; 3]>,
}

/// Split a finite TT Julian Day into integer day and fraction.
/// Configured model Delta-T does not enter this observation conversion;
/// non-finite arguments return [`Error::NonFiniteArgument`].
fn scalar_tt(tt: f64) -> Result<TwoPartEpoch> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    let hi = tt.floor();
    Ok(TwoPartEpoch { hi, lo: tt - hi })
}

/// Select the physical stellar observer at a split TT epoch.
///
/// Heliocentric work reads the Sun, barycentric work uses the SSB origin,
/// and geocentric/topocentric work reads Earth with an optional WGS84 site.
/// A required missing site is refused before state acquisition. Site rotation
/// uses bundled TT-minus-UT1, native precession and nutation without polar
/// motion. Missing states and frame/clock errors propagate without substitution.
fn observer_via(
    state_at: &StarStateAt<'_>,
    time: TwoPartEpoch,
    center: i32,
    topo: Option<StarTopoSite>,
) -> Result<NativeStarObserver> {
    if center == FLG_BARYCTR {
        return Ok(NativeStarObserver {
            time,
            position: [0.0; 3],
            velocity: [0.0; 3],
            velocity_per_tt: [0.0; 3],
            surface: None,
        });
    }
    let site = if center == FLG_TOPOCTR {
        Some(topo.ok_or(Error::ConfigurationMissing {
            missing_config: "observer_location",
            suggestion: "topocentric star places need an explicit site",
        })?)
    } else {
        None
    };
    let (position, velocity) = state_at(if center == FLG_HELCTR { SUN } else { EARTH }, time)?;
    let barycentric_velocity = velocity;
    let tdb_rate = crate::time::tdb_rate(time.value())?;
    if let Some(site) = site {
        let tt = time.value();
        let tdb = time.hi
            + (time.lo + reference_frames::tdb_minus_tt_seconds(time.hi, time.lo)? / 86400.0);
        let delta = crate::smh2016::bundled_deltat_seconds(tt)?;
        let (position, velocity, surface) = reference_frames::barycentric_site_state(
            position,
            velocity,
            site.lon_deg(),
            site.lat_deg(),
            site.height_m(),
            tt,
            delta,
            tdb,
        )?;
        Ok(NativeStarObserver {
            time,
            position,
            velocity,
            velocity_per_tt: std::array::from_fn(|i| {
                velocity[i] + barycentric_velocity[i] * (tdb_rate - 1.0)
            }),
            surface: Some(surface),
        })
    } else {
        Ok(NativeStarObserver {
            time,
            position,
            velocity,
            velocity_per_tt: velocity.map(|v| v * tdb_rate),
            surface: None,
        })
    }
}

/// Astrometric stellar ray and light time for a precomputed observer.
///
/// Catalog angles, parallax and Doppler-corrected space motion define the
/// initial ICRS AU/AU-day state. Propagate over
/// `TDB + dot(unit(initial), observer)/c − J2000` days, then subtract the
/// observer. This linear stellar motion does not use the planetary
/// retarded-target iteration. Distance divided by the exact-SI AU/day light
/// speed supplies the astrometric delay. The enclosing entry validates
/// catalog acquisition, observation epoch and provider states.
pub(super) fn astrometric_ray(data: &StarData, observer: NativeStarObserver) -> ([f64; 3], f64) {
    let (initial, velocity) = crate::stars::star_icrs_initial_state_catalog_angles(data);
    astrometric_ray_from_initial((initial, velocity), observer)
}

/// Propagate an owned constructor at the retained observer's original epoch.
/// Initial vectors are ICRS AU/AU-day; no catalog or clock is reacquired here.
fn astrometric_ray_from_initial(
    (initial, velocity): ([f64; 3], [f64; 3]),
    observer: NativeStarObserver,
) -> ([f64; 3], f64) {
    let norm = super::vec3_dist(initial);
    let unit: [f64; 3] = std::array::from_fn(|i| initial[i] / (norm + f64::MIN_POSITIVE));
    let projection = (unit[0] * observer.position[0]
        + unit[1] * observer.position[1]
        + unit[2] * observer.position[2])
        / super::FRAME_LIGHT_SPEED_AU_DAY;
    let tdb = observer.time.hi
        + (observer.time.lo
            + crate::time::source_tdb_minus_tt_seconds(observer.time.hi, observer.time.lo)
                / 86400.0);
    let span = tdb + projection - crate::constants::J2000;
    let ray = std::array::from_fn(|i| {
        let propagated = initial[i] + velocity[i] * span;
        propagated - observer.position[i]
    });
    (ray, super::vec3_dist(ray) / super::FRAME_LIGHT_SPEED_AU_DAY)
}

/// Propagate the catalog state with its TDB/projected-light-time Jacobian.
fn astrometric_state_from_initial(
    (initial, velocity): ([f64; 3], [f64; 3]),
    observer: NativeStarObserver,
) -> Result<([f64; 3], f64, f64)> {
    let (ray, light_time) = astrometric_ray_from_initial((initial, velocity), observer);
    let norm = super::vec3_dist(initial);
    let epoch_rate = crate::time::tdb_rate(observer.time.value())?
        + (0..3)
            .map(|i| (initial[i] / norm) * observer.velocity_per_tt[i])
            .sum::<f64>()
            / super::FRAME_LIGHT_SPEED_AU_DAY;
    let relative_velocity =
        std::array::from_fn(|i| velocity[i] * epoch_rate - observer.velocity_per_tt[i]);
    let (_, radial_rate) = star_geometric_range_state(ray, relative_velocity);
    Ok((ray, light_time, radial_rate))
}

/// Anchor vector correction in ICRS AU before mean-frame projection.
/// Both toggles select the unchanged astrometric ray; `nogdefl` selects only
/// aberration. `noaberr` computes full apparent minus (aberrated minus astro)
/// component by component, with those calls in that order. This is a vector
/// differential, distinct from public fixed-star spherical NOABERR. The
/// callback receives default or empty deflector-list selection; consumed
/// errors propagate and unused branches remain lazy. Derived IEEE values are
/// retained rather than reclassified as invalid original caller epochs.
fn native_anchor_correction_via(
    astro: [f64; 3],
    noaberr: bool,
    nogdefl: bool,
    apparent: impl Fn(bool) -> Result<[f64; 3]>,
) -> Result<[f64; 3]> {
    if noaberr && nogdefl {
        return Ok(astro);
    }
    if nogdefl {
        return apparent(false);
    }
    let full = apparent(true)?;
    if noaberr {
        let aberrated = apparent(false)?;
        return Ok(std::array::from_fn(|i| full[i] - (aberrated[i] - astro[i])));
    }
    Ok(full)
}

/// Strict geocentric native catalog-anchor ICRS ray at one scalar TT sample.
/// Fresh scalar TT parts feed actual Earth split-epoch acquisition and catalog
/// TDB/projection propagation. The acquired observer velocity and astrometric
/// distance/c drive default or empty-list apparent reduction. Sun/Jupiter/Saturn
/// original and closest-approach reads preserve split epoch and strict failures.
/// Output is a reduced observer-relative AU vector without Vondrák projection,
/// interpolation or public error promotion; those belong to the outer anchor.
/// Finite caller admission and raw state errors keep their existing contracts;
/// no configured model clock, optional LEB supports or spherical differential
/// is substituted. There is no planetary target light-time iteration here.
fn native_anchor_ray_via(
    state_at: &StarStateAt<'_>,
    data: &StarData,
    tt: f64,
    noaberr: bool,
    nogdefl: bool,
) -> Result<[f64; 3]> {
    let observer = observer_via(state_at, scalar_tt(tt)?, 0, None)?;
    let (astro, light_time) = astrometric_ray(data, observer);
    let positions = |id, time| state_at(id, time).map(|state| state.0);
    native_anchor_correction_via(astro, noaberr, nogdefl, |include_deflectors| {
        native_apparent_ray_with_deflectors(
            &positions,
            observer.time,
            astro,
            light_time,
            observer.position,
            observer.velocity,
            None,
            include_deflectors,
        )
    })
}

/// Read a sidereal anchor's stored state at the scalar TT Julian day.
/// Missing aliases return the native UnknownBody error; coverage and malformed
/// data failures are mapped by the selected strict reader.
fn live_stored_anchor_state_via<P: StateProvider>(
    provider: &P,
    body: i32,
    time: TwoPartEpoch,
    vector_error: &impl Fn(Error) -> Error,
) -> Result<([f64; 3], [f64; 3])> {
    if !provider.has_body(body) {
        return Err(Error::UnknownBody { body_id: body });
    }
    provider
        .geometric_state(body, time.value())
        .map(|state| (state.pos(), state.vel()))
        .map_err(vector_error)
}

/// Reduce a live sidereal anchor constructor at its actual helper boundary.
/// `initial` owns ICRS AU/AU-day vectors computed from original native keyword
/// units before this call. The ecliptic helper acquires its main ephemeris
/// before scalar TT Time; Galactic Center RA acquires Time before main. Earth
/// then consumes the original split epoch, followed by selected strict native
/// deflectors. State callbacks for a stored vector adapter must consume its
/// scalar epoch value instead, preserving the independent frame clock.
/// No row or cache is accessed here, and consumed failures are not replaced.
pub(crate) fn live_anchor_ray_via(
    state_at: &StarStateAt<'_>,
    initial: ([f64; 3], [f64; 3]),
    tt: f64,
    noaberr: bool,
    nogdefl: bool,
    role: crate::stars::AnchorRayRole,
    acquire_main: impl FnOnce() -> Result<()>,
) -> Result<[f64; 3]> {
    let time = match role {
        crate::stars::AnchorRayRole::EclipticLongitude => {
            acquire_main()?;
            scalar_tt(tt)?
        }
        crate::stars::AnchorRayRole::GalacticRightAscension => {
            let time = scalar_tt(tt)?;
            acquire_main()?;
            time
        }
    };
    live_anchor_ray_at_time_via(state_at, initial, time, noaberr, nogdefl)
}

/// Reduce a catalog initial ICRS state through the retained observer epoch,
/// apparent corrections and selected output frame. Positions are AU and rates
/// are per TT day; provider and frame errors propagate.
fn live_anchor_ray_at_time_via(
    state_at: &StarStateAt<'_>,
    initial: ([f64; 3], [f64; 3]),
    time: TwoPartEpoch,
    noaberr: bool,
    nogdefl: bool,
) -> Result<[f64; 3]> {
    let observer = observer_via(state_at, time, 0, None)?;
    let (astro, light_time) = astrometric_ray_from_initial(initial, observer);
    let positions = |id, time| state_at(id, time).map(|state| state.0);
    native_anchor_correction_via(astro, noaberr, nogdefl, |include_deflectors| {
        native_apparent_ray_with_deflectors(
            &positions,
            observer.time,
            astro,
            light_time,
            observer.position,
            observer.velocity,
            None,
            include_deflectors,
        )
    })
}

/// Ecliptic longitude/latitude degrees and AU distance for a retained Time/ray.
/// J2000 uses the IAU2006 mean-obliquity matrix without frame bias, followed
/// by the shared spherical coordinate conversion. Of-date uses TT nutation
/// and TDB precession/obliquity with the same native frame transformation.
/// Original Time parts remain intact for TDB; scalar TT drives nutation.
/// Frame/conversion guards keep their typed errors and derived IEEE behavior;
/// this function does not apply parallax floors, output flags or rate stencils.
fn ecliptic_place(ray: [f64; 3], time: TwoPartEpoch, j2000: bool) -> Result<[f64; 3]> {
    // Frame rotations change direction, not range. Keeping the ICRS norm
    // prevents rotation roundoff from becoming a radial speed when distant
    // stellar ranges are differentiated at neighboring epochs.
    let distance = super::vec3_dist(ray);
    if j2000 {
        let mut position =
            super::cartesian_to_spherical_deg(super::rotate_icrs_to_ecliptic_j2000(ray)?)?;
        position[2] = distance;
        return Ok(position);
    }
    let tdb =
        time.hi + (time.lo + reference_frames::tdb_minus_tt_seconds(time.hi, time.lo)? / 86400.0);
    let (psi, deps) = reference_frames::iau2000a_tenth_microarcsec(time.value())?;
    let psi = psi * reference_frames::FRAME_TENTH_USEC_2_RAD;
    let deps = deps * reference_frames::FRAME_TENTH_USEC_2_RAD;
    let eps = reference_frames::mean_obliquity_arcsec(tdb)? * reference_frames::FRAME_ASEC2RAD;
    let n = reference_frames::nutation_matrix(eps, eps + deps, psi)?;
    let m = reference_frames::equinox_matrix(n, reference_frames::precession_matrix(tdb)?)?;
    let frame = reference_frames::ecliptic_date_matrix(m, eps, deps)?;
    let mut position = reference_frames::frame_latlon_deg(frame, ray)?;
    position[2] = distance;
    Ok(position)
}

/// Select the public star correction branch using its spherical counterpart.
/// Both-disabled returns geometric; NOGDEFL returns the empty-list apparent
/// place (retaining Earth deflection for a surface observer). NOABERR alone
/// returns `g_lon+wrap180(a_lon-b_lon), g_lat+(a_lat-b_lat), g_distance`.
/// It is a reported-coordinate differential, not a deflection-only ray.
/// `apparent` receives true/default or false/empty deflector-list selection.
/// Errors propagate in the actual branch call order; unused branches are lazy.
pub(crate) fn corrected_place_via(
    noaberr: bool,
    nogdefl: bool,
    geometric: impl Fn() -> Result<[f64; 3]>,
    apparent: impl Fn(bool) -> Result<[f64; 3]>,
) -> Result<[f64; 3]> {
    if noaberr && nogdefl {
        return geometric();
    }
    if nogdefl {
        return apparent(false);
    }
    if noaberr {
        let g = geometric()?;
        let a = apparent(true)?;
        let b = apparent(false)?;
        let dlon = reference_frames::floored_mod(a[0] - b[0] + 180.0, 360.0) - 180.0;
        return Ok([
            reference_frames::floored_mod(g[0] + dlon, 360.0),
            g[1] + (a[1] - b[1]),
            g[2],
        ]);
    }
    apparent(true)
}

/// Catalog place from a frozen native observer through actual source branches.
/// `data` is the existing pinned catalog row. Observer epoch/state and optional
/// surface are reused unchanged across the differential calls. Default PPN
/// resolves Sun then Jupiter-center/barycenter then Saturn-center/barycenter,
/// with original and closest-approach split-epoch reads; missing/ranged supports
/// propagate. The existing shared native apparent tail handles Earth gating
/// and aberration without a planetary light-time solve. Output is ecliptic
/// degrees/AU; no sidereal/output shaping or synthetic catalog data is added.
fn position_from_observer_via(
    state_at: &StarStateAt<'_>,
    data: &StarData,
    observer: NativeStarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
) -> Result<StarSample> {
    let initial = crate::stars::star_icrs_initial_state_catalog_angles(data);
    let (ray, light_time, radial_rate) = astrometric_state_from_initial(initial, observer)?;
    position_from_astrometric_via(
        state_at,
        (ray, light_time, radial_rate),
        observer,
        noaberr,
        nogdefl,
        j2000,
    )
}

/// Apparent/frame tail over the actual ray and retained observer.
/// This factors unchanged correction calls without selecting a row or Time.
fn position_from_astrometric_via(
    state_at: &StarStateAt<'_>,
    astrometric: ([f64; 3], f64, f64),
    observer: NativeStarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
) -> Result<StarSample> {
    let (ray, light_time, radial_rate) = astrometric;
    let positions = |id, time| state_at(id, time).map(|state| state.0);
    let mut position = corrected_place_via(
        noaberr,
        nogdefl,
        || ecliptic_place(ray, observer.time, j2000),
        |default_deflectors| {
            let apparent = native_apparent_ray_with_deflectors(
                &positions,
                observer.time,
                ray,
                light_time,
                observer.position,
                observer.velocity,
                observer.surface,
                default_deflectors,
            )?;
            ecliptic_place(apparent, observer.time, j2000)
        },
    )?;
    position[2] = super::vec3_dist(ray);
    Ok(StarSample {
        position,
        radial_rate,
    })
}

/// Acquire the complete native astrometric tuple before computing its
/// initial catalog state. The owned row remains retained for the calculation.
#[allow(clippy::too_many_arguments)]
fn position_from_live_observer_via(
    state_at: &StarStateAt<'_>,
    rows: &dyn StarRowAcquisition,
    id: i32,
    observer: NativeStarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
) -> Result<StarSample> {
    let row = rows.acquire_row(id, StarRowBackend::Native)?;
    position_from_row_observer_via(state_at, row.as_ref(), observer, noaberr, nogdefl, j2000)
}

/// Compute native catalog geometry and apparent reduction over one retained row.
fn position_from_row_observer_via(
    state_at: &StarStateAt<'_>,
    row: &dyn AcquiredStarRow,
    observer: NativeStarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
) -> Result<StarSample> {
    let initial = crate::stars::star_initial_state_from_native(row.native_inputs()?);
    let (ray, lt, radial_rate) = astrometric_state_from_initial(initial, observer)?;
    position_from_astrometric_via(
        state_at,
        (ray, lt, radial_rate),
        observer,
        noaberr,
        nogdefl,
        j2000,
    )
}

/// Live row through the stored-backed topocentric vector channel.
/// The caller already owns its surface observer. Construct its scalar epoch
/// before selecting the native row. The frame/observer retains its split epoch;
/// the LEB vector adapter consumes its scalar TT sum for each state query.
#[allow(clippy::too_many_arguments)]
pub(crate) fn live_topocentric_place_via(
    provider: &impl StateProvider,
    vector_error: Option<&dyn Fn(Error) -> Error>,
    rows: &dyn StarRowAcquisition,
    id: i32,
    tt: f64,
    observer: &crate::stars::StarObserver,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
) -> Result<StarSample> {
    let observer = NativeStarObserver {
        time: scalar_tt(tt)?,
        position: observer.pos_au(),
        velocity: observer.vel_au_day(),
        velocity_per_tt: observer.vel_au_day(),
        surface: observer.site_gcrs_au(),
    };
    let state_at = |id, time: TwoPartEpoch| {
        let state = provider
            .geometric_state(id, time.value())
            .map_err(|error| {
                if let Some(mapper) = vector_error {
                    mapper(error)
                } else {
                    error
                }
            })?;
        Ok((state.pos(), state.vel()))
    };
    position_from_live_observer_via(&state_at, rows, id, observer, noaberr, nogdefl, j2000)
}

/// Scalar vector retry over the exact already selected stored provider.
/// It constructs a fresh scalar-TT Time and observer, then reacquires the row
/// and original native keyword units. Explicit LEB has no JPL/Main substitution.
#[allow(clippy::too_many_arguments)]
pub(crate) fn live_stored_vector_place_via(
    provider: &impl StateProvider,
    vector_error: Option<&dyn Fn(Error) -> Error>,
    rows: &dyn StarRowAcquisition,
    id: i32,
    tt: f64,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
    center: i32,
) -> Result<StarSample> {
    let state_at = |id, time: TwoPartEpoch| {
        let state = provider
            .geometric_state(id, time.value())
            .map_err(|error| {
                if let Some(mapper) = vector_error {
                    mapper(error)
                } else {
                    error
                }
            })?;
        Ok((state.pos(), state.vel()))
    };
    let mut observer = observer_via(&state_at, scalar_tt(tt)?, center, None)?;
    // The stored adapter is already parameterized by TT, unlike SPK TDB.
    observer.velocity_per_tt = observer.velocity;
    position_from_live_observer_via(&state_at, rows, id, observer, noaberr, nogdefl, j2000)
}

/// Whole stored-batch vector retry after the owning caught failure/log stage.
/// The original shared UT1 Time is retained; all vector observers precede fresh
/// per-occurrence rows. Successful vector retries overwrite stored slots;
/// skipped failures preserve their earlier assigned stored value. The
/// retained provider supplies scalar-TT LEB states and cannot acquire JPL.
#[allow(clippy::too_many_arguments)]
pub(crate) fn live_stored_vector_batch_via(
    provider: &impl StateProvider,
    vector_error: Option<&dyn Fn(Error) -> Error>,
    serving: &SealedStarServing<'_, impl StateProvider>,
    resolved: &[Option<StarMatch>],
    row_owners: &[Option<std::sync::Arc<dyn StarRowAcquisition>>],
    partial: Vec<Option<WrappedStarEntry>>,
    time: TwoPartEpoch,
    flags: i32,
    skip: bool,
) -> Result<Vec<Option<WrappedStarEntry>>> {
    let state_at = |id, time: TwoPartEpoch| {
        let state = provider
            .geometric_state(id, time.value())
            .map_err(|error| {
                if let Some(mapper) = vector_error {
                    mapper(error)
                } else {
                    error
                }
            })?;
        Ok((state.pos(), state.vel()))
    };
    let selected = resolved
        .iter()
        .enumerate()
        .filter_map(|(index, star)| star.clone().map(|star| (index, star)))
        .collect();
    native_batch_samples_at_time_via(
        &|| Ok(()),
        None,
        &state_at,
        &|_, _, error| error,
        serving,
        selected,
        partial,
        time,
        flags,
        skip,
        Some(row_owners),
    )
}

/// Raw public scalar star place/rates at fresh scalar TT central/side Times.
///
/// Catalog miss precedes observer reads. Each position uses a fresh epoch and
/// observer; central is evaluated before previous and next. SPEED differentiates
/// reported positions at h=0.005 days for surface observer or 0.5 otherwise, with a single
/// longitude ±360 fold, divided by 2h. Other calls return positive-zero rates.
/// All source callback/frame/clock errors propagate; output shaping is separate.
#[allow(clippy::too_many_arguments)]
#[cfg(test)]
fn native_raw_place_via(
    state_at: &StarStateAt<'_>,
    range_at: &StarRangeAt<'_>,
    id: i32,
    tt: f64,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
    center: i32,
    topo: Option<StarTopoSite>,
    speed: bool,
) -> Result<[f64; 6]> {
    native_raw_place_with_rows_via(
        state_at, range_at, id, tt, noaberr, nogdefl, j2000, center, topo, speed, None,
    )
}

/// Scalar samples with optional live row acquisition after each observer.
/// The enclosing scalar callback owns velocity's extra admission and Main
/// selection. This helper does not perform those stages a second time.
#[allow(clippy::too_many_arguments)]
fn native_raw_place_with_rows_via(
    state_at: &StarStateAt<'_>,
    range_at: &StarRangeAt<'_>,
    id: i32,
    tt: f64,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
    center: i32,
    topo: Option<StarTopoSite>,
    speed: bool,
    rows: Option<&dyn StarRowAcquisition>,
) -> Result<[f64; 6]> {
    native_raw_place_with_rows_and_main_via(
        state_at, range_at, id, tt, noaberr, nogdefl, j2000, center, topo, speed, rows, None,
    )
}

/// Evaluate catalog places with an optional shared main-image owner.
/// Each consumed observation obtains its actual retained image before state
/// acquisition. Static catalog rows and live native rows share the same
/// physical observer and apparent-place reduction.
#[allow(clippy::too_many_arguments)]
fn native_raw_place_with_rows_and_main_via(
    state_at: &StarStateAt<'_>,
    range_at: &StarRangeAt<'_>,
    id: i32,
    tt: f64,
    noaberr: bool,
    nogdefl: bool,
    j2000: bool,
    center: i32,
    topo: Option<StarTopoSite>,
    speed: bool,
    rows: Option<&dyn StarRowAcquisition>,
    main: Option<&dyn Fn(f64) -> Result<()>>,
) -> Result<[f64; 6]> {
    let data = if rows.is_none() {
        let data = fixed_star_data(id).ok_or_else(|| Error::StarNotFound {
            star_id: id.to_string(),
            search_type: Some("sequential"),
        })?;
        Some(data)
    } else {
        None
    };
    let at = |tt| {
        (|| {
            let time = scalar_tt(tt)?;
            if let Some(main) = main {
                main(tt)?;
            }
            let observer = observer_via(state_at, time, center, topo)?;
            if let Some(rows) = rows {
                position_from_live_observer_via(
                    state_at, rows, id, observer, noaberr, nogdefl, j2000,
                )
            } else {
                position_from_observer_via(
                    state_at,
                    data.as_ref().expect("static row acquired"),
                    observer,
                    noaberr,
                    nogdefl,
                    j2000,
                )
            }
        })()
        .map_err(|error| range_at(center, tt, error))
    };
    if speed {
        let h = if center == FLG_TOPOCTR { 0.005 } else { 0.5 };
        star_velocity_samples_via(at, tt, h)
    } else {
        let p = at(tt)?.position;
        Ok([p[0], p[1], p[2], 0.0, 0.0, 0.0])
    }
}

/// Whether batch calculation needs the scalar stellar reduction.
/// Topocentric, fixed-epoch sidereal and active SIDBIT planes use that route;
/// ordinary geocentric batches share their observation state.
pub(crate) fn scalar_batch_required(
    flags: i32,
    selection: crate::sidereal::SiderealSelection,
) -> bool {
    flags & FLG_TOPOCTR != 0
        || (flags & FLG_SIDEREAL != 0
            && (crate::sidereal::is_fixed_epoch_request(flags, selection.mode())
                || (selection.bits() & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE) != 0
                    && !crate::sidereal::SIDBIT_PROJECTION_SUPPRESS_MODES
                        .contains(&selection.mode()))))
}

/// Geocentric stellar batch with shared UT1/TT observation and observer.
///
/// Resolve names in input order before acquiring a clock or image. With no
/// resolved stars, return immediately. The central observer retains split TT;
/// SPEED uses fresh TT observers at ±0.5 day. Each star then applies its own
/// catalog motion and corrections through the common flag layer. Requested
/// skip-errors behavior applies inside the stellar loop; shared observation
/// or provider errors remain failures for the whole batch.
#[allow(clippy::too_many_arguments)]
fn native_batch_with_main_via(
    ensure_main: &dyn Fn() -> Result<()>,
    state_at: &StarStateAt<'_>,
    range_at: &StarRangeAt<'_>,
    serving: &SealedStarServing<'_, ContextStarProvider<'_>>,
    selection: crate::sidereal::SiderealSelection,
    names: &[&str],
    ut: f64,
    flags_in: i32,
    topo: Option<StarTopoSite>,
    skip: bool,
) -> Result<Vec<Option<WrappedStarEntry>>> {
    if scalar_batch_required(flags_in, selection) {
        return serving.batch_fixstars_ut(names, ut, flags_in, topo, skip);
    }
    let mut resolved = Vec::new();
    let output = vec![None; names.len()];
    for (index, name) in names.iter().enumerate() {
        match resolve_star_ref(name) {
            Ok(star) => resolved.push((index, star)),
            Err(StarResolveError::Miss(failure)) if !skip => return Err(failure.to_error(name)),
            Err(StarResolveError::Miss(_)) => {}
            Err(StarResolveError::Refused(error)) => return Err(error),
        }
    }
    if resolved.is_empty() {
        return Ok(output);
    }
    native_batch_samples_via(
        ensure_main,
        None,
        state_at,
        range_at,
        serving,
        resolved,
        output,
        ut,
        flags_in,
        skip,
        None,
    )
}

/// Shared observer acquisition for externally resolved, live catalog entries.
/// Empty slots never acquire a provider; duplicate numerical ids still select
/// separate row leases for each consumed central/minus/plus sample.
#[allow(clippy::too_many_arguments)]
fn native_resolved_batch_with_main_via(
    ensure_main: &dyn Fn() -> Result<()>,
    observation: Option<&dyn Fn(StarObservationClock) -> Result<()>>,
    state_at: &StarStateAt<'_>,
    range_at: &StarRangeAt<'_>,
    serving: &SealedStarServing<'_, ContextStarProvider<'_>>,
    resolved: &[Option<StarMatch>],
    row_owners: &[Option<std::sync::Arc<dyn StarRowAcquisition>>],
    ut: f64,
    flags: i32,
    skip: bool,
) -> Result<Vec<Option<WrappedStarEntry>>> {
    let selected = resolved
        .iter()
        .enumerate()
        .filter_map(|(index, star)| star.clone().map(|star| (index, star)))
        .collect::<Vec<_>>();
    let output = vec![None; resolved.len()];
    if selected.is_empty() {
        return Ok(output);
    }
    native_batch_samples_via(
        ensure_main,
        observation,
        state_at,
        range_at,
        serving,
        selected,
        output,
        ut,
        flags,
        skip,
        Some(row_owners),
    )
}

/// Unchanged shared original-epoch observer and public flag tail.
/// Live selection occurs only inside the per-star loop, after all shared
/// observers. An adapter's actual exception classification controls skipping;
/// the independent legacy route retains its prior typed catch policy.
#[allow(clippy::too_many_arguments)]
fn native_batch_samples_via(
    ensure_main: &dyn Fn() -> Result<()>,
    observation: Option<&dyn Fn(StarObservationClock) -> Result<()>>,
    state_at: &StarStateAt<'_>,
    range_at: &StarRangeAt<'_>,
    serving: &SealedStarServing<'_, ContextStarProvider<'_>>,
    resolved: Vec<(usize, StarMatch)>,
    output: Vec<Option<WrappedStarEntry>>,
    ut: f64,
    flags_in: i32,
    skip: bool,
    row_owners: Option<&[Option<std::sync::Arc<dyn StarRowAcquisition>>]>,
) -> Result<Vec<Option<WrappedStarEntry>>> {
    let (hi, lo) = crate::smh2016::timescale_parts_from_ut1(ut)?;
    native_batch_samples_at_time_via(
        ensure_main,
        observation,
        state_at,
        range_at,
        serving,
        resolved,
        output,
        TwoPartEpoch { hi, lo },
        flags_in,
        skip,
        row_owners,
    )
}

/// Shared-vector tail retaining an already constructed original epoch.
/// Stored fallback passes its existing UT1 Time here and uses only its retained
/// LEB provider. Main acquisition is supplied only by the native caller.
#[allow(clippy::too_many_arguments)]
fn native_batch_samples_at_time_via(
    ensure_main: &dyn Fn() -> Result<()>,
    observation: Option<&dyn Fn(StarObservationClock) -> Result<()>>,
    state_at: &StarStateAt<'_>,
    range_at: &StarRangeAt<'_>,
    serving: &SealedStarServing<'_, impl StateProvider>,
    resolved: Vec<(usize, StarMatch)>,
    mut output: Vec<Option<WrappedStarEntry>>,
    time: TwoPartEpoch,
    flags_in: i32,
    skip: bool,
    row_owners: Option<&[Option<std::sync::Arc<dyn StarRowAcquisition>>]>,
) -> Result<Vec<Option<WrappedStarEntry>>> {
    let ret = fixstar_ret_flags(flags_in, true);
    let flags = preprocess_star_flags(flags_in);
    let center = fixstar_observation_center(flags);
    let noaberr = flags & (FLG_NOABERR | FLG_TRUEPOS) != 0 || center != 0;
    let nogdefl = flags & (FLG_NOGDEFL | FLG_TRUEPOS) != 0 || center != 0;
    let j2000 = flags & FLG_J2000 != 0;
    ensure_main()?;
    let observed = |sample: TwoPartEpoch| {
        observer_via(state_at, sample, center, None)
            .map_err(|error| range_at(center, time.value(), error))
    };
    let central = observed(time)?;
    let select = |sample| match observation {
        Some(observation) => observation(StarObservationClock::SharedSample(sample)),
        None => Ok(()),
    };
    let sides = if flags & FLG_SPEED != 0 {
        if let Some(observation) = observation {
            observation(StarObservationClock::SharedSpeedTimes)?;
        }
        select(StarObservationSample::Previous)?;
        let previous = observed(scalar_tt(time.value() - 0.5)?)?;
        select(StarObservationSample::Next)?;
        let next = observed(scalar_tt(time.value() + 0.5)?)?;
        Some((previous, next))
    } else {
        None
    };
    for (index, star) in resolved {
        let rows = row_owners.map(|owners| {
            owners[index]
                .as_ref()
                .expect("validated occurrence owner")
                .as_ref()
        });
        let result = (|| {
            let data = if rows.is_none() {
                let data = fixed_star_data(star.star_id).ok_or_else(|| Error::StarNotFound {
                    star_id: star.star_id.to_string(),
                    search_type: Some("sequential"),
                })?;
                Some(data)
            } else {
                None
            };
            let place_at = |observer| {
                if let Some(rows) = rows {
                    position_from_live_observer_via(
                        state_at,
                        rows,
                        star.star_id,
                        observer,
                        noaberr,
                        nogdefl,
                        j2000,
                    )
                } else {
                    position_from_observer_via(
                        state_at,
                        data.as_ref().expect("static row acquired"),
                        observer,
                        noaberr,
                        nogdefl,
                        j2000,
                    )
                }
            };
            select(StarObservationSample::Central)?;
            let central_place = place_at(central)?;
            let p = central_place.position;
            let rates = if let Some((previous, next)) = sides {
                select(StarObservationSample::Previous)?;
                let previous = place_at(previous)?;
                select(StarObservationSample::Next)?;
                let next = place_at(next)?;
                let mut longitude = next.position[0] - previous.position[0];
                if longitude > 180.0 {
                    longitude -= 360.0;
                } else if longitude < -180.0 {
                    longitude += 360.0;
                }
                [
                    longitude,
                    next.position[1] - previous.position[1],
                    central_place.radial_rate,
                ]
            } else {
                [0.0; 3]
            };
            let raw = [p[0], p[1], p[2], rates[0], rates[1], rates[2]];
            let shaped = serving.shape_native_place(raw, time.value(), flags, j2000)?;
            Ok((shaped, star.canonical, ret))
        })();
        match result {
            Ok(entry) => output[index] = Some(entry),
            Err(error) => {
                let skippable = rows
                    .map(|rows| rows.is_skippable(&error))
                    .transpose()?
                    .unwrap_or(true);
                if !skip || !skippable {
                    return Err(error);
                }
            }
        }
    }
    Ok(output)
}

/// Test-only delegation over an injected state channel, without file acquisition.
/// The production callback supplies its consumed main selector to the sibling.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn native_batch_via(
    state_at: &StarStateAt<'_>,
    range_at: &StarRangeAt<'_>,
    serving: &SealedStarServing<'_, ContextStarProvider<'_>>,
    selection: crate::sidereal::SiderealSelection,
    names: &[&str],
    ut: f64,
    flags_in: i32,
    topo: Option<StarTopoSite>,
    skip: bool,
) -> Result<Vec<Option<WrappedStarEntry>>> {
    native_batch_with_main_via(
        &|| Ok(()),
        state_at,
        range_at,
        serving,
        selection,
        names,
        ut,
        flags_in,
        topo,
        skip,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::stars::star_velocity_via;
    use std::cell::RefCell;

    #[test]
    fn auto_stars_select_actual_observer_and_frame_coverage_without_acquisition() {
        let path = crate::test_data::core_leb_path();
        let mut context = EphemerisContext::new().with_mode(ProviderMode::Auto);
        context.attach_leb_file(crate::leb::LebFile::open(&path).unwrap());
        let composite = crate::leb::CompositeIndex::from_files(context.leb_files()).unwrap();
        let (start, end) = composite.body_coverage(EARTH).unwrap();
        assert!(
            auto_star_context(&context, 2451545.0, false, FLG_SPEED)
                .unwrap()
                .is_none()
        );
        for tt in [start - 4.0, end + 4.0] {
            let selected = auto_star_context(&context, tt, false, FLG_SPEED)
                .unwrap()
                .unwrap();
            assert_eq!(selected.mode(), ProviderMode::NativeJpl);
            assert_eq!(selected.leb_files().len(), context.leb_files().len());
            assert!(
                selected.spk_registry().is_empty(),
                "selection opens no main image"
            );
        }
        let ut = crate::time::default_timescale_ut1_from_tt(start - 4.0).unwrap();
        assert_eq!(
            auto_star_context(&context, ut, true, 0)
                .unwrap()
                .unwrap()
                .mode(),
            ProviderMode::NativeJpl
        );
        let forced = context.clone().with_mode(ProviderMode::SealedLeb);
        assert!(
            auto_star_context(&forced, start - 4.0, false, FLG_SPEED)
                .unwrap()
                .is_none()
        );
        assert_eq!(context.mode(), ProviderMode::Auto);
    }

    #[test]
    fn auto_stellar_topocentric_and_empty_inputs_select_owned_native_resources() {
        let context = EphemerisContext::new().with_mode(ProviderMode::Auto);
        let selected = auto_star_context(&context, 2451545.0, false, FLG_TOPOCTR)
            .unwrap()
            .unwrap();
        assert_eq!(selected.mode(), ProviderMode::NativeJpl);
        assert!(selected.spk_registry().is_empty());
        assert!(matches!(
            auto_star_context(&context, f64::NAN, false, 0),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn geometric_stellar_range_rate_matches_eighty_digit_dot_reference() {
        let ray = [1e7, 2e7, 3e7];
        let velocity = [0.001, -0.002, 0.003];
        // Dot(r,v)/norm(r), with exact binary64 inputs, at 80 decimal digits.
        let reference = 0.001_603_567_451_474_546_3;
        let (range, rate) = star_geometric_range_state(ray, velocity);
        assert!((range - 37_416_573.867_739_42).abs() <= f64::EPSILON * range);
        assert!((rate - reference).abs() <= 8.0 * f64::EPSILON * reference);
        assert_eq!(star_geometric_range_state(ray, [0.0; 3]).1, 0.0);
        assert_eq!(star_geometric_range_state([0.0; 3], velocity), (0.0, 0.0));
    }

    #[test]
    fn changing_reference_frames_cannot_create_stellar_radial_motion() {
        // The exact Euclidean norm is sqrt(14) × 10^7 AU; this decimal was
        // evaluated independently at 80 decimal digits, then rounded once.
        let ray = [1e7, 2e7, 3e7];
        let reference = 37_416_573.867_739_42;
        let distance = super::super::vec3_dist(ray);
        assert!((distance - reference).abs() <= f64::EPSILON * reference);
        for epoch in [2415020.5, 2433282.5, 2451545.0, 2460310.5, 2488069.5] {
            for j2000 in [false, true] {
                let at = |jd| ecliptic_place(ray, scalar_tt(jd)?, j2000);
                let result = star_velocity_via(at, epoch, 0.5).unwrap();
                assert_eq!(result[2], distance);
                assert_eq!(result[5], 0.0, "fixed ICRS ray at {epoch}, J2000={j2000}");
            }
        }
    }
    use std::sync::{
        Arc, Mutex,
        atomic::{AtomicUsize, Ordering},
    };

    struct LiveRows {
        trace: Arc<Mutex<Vec<String>>>,
        label: &'static str,
        next: AtomicUsize,
        refuse: bool,
        refuse_classification: bool,
    }
    struct LiveRow {
        trace: Arc<Mutex<Vec<String>>>,
        label: &'static str,
        serial: usize,
    }
    impl Drop for LiveRow {
        fn drop(&mut self) {
            self.trace
                .lock()
                .unwrap()
                .push(format!("{}:drop:{}", self.label, self.serial));
        }
    }
    impl AcquiredStarRow for LiveRow {
        fn native_inputs(&self) -> Result<crate::stars::NativeStarInputs> {
            self.trace
                .lock()
                .unwrap()
                .push(format!("{}:keywords:{}", self.label, self.serial));
            Ok(crate::stars::NativeStarInputs {
                ra_hours: self.serial as f64 + if self.label == "a" { 0.0 } else { 8.0 },
                dec_degrees: 2.0,
                ra_mas_per_year: 0.0,
                dec_mas_per_year: 0.0,
                parallax_mas: 1000.0,
                radial_km_per_s: 0.0,
            })
        }
    }
    impl StarRowAcquisition for LiveRows {
        fn velocity_admission(&self, _: i32) -> Result<()> {
            panic!("shared native batch has no scalar velocity admission")
        }
        fn acquire_row(&self, _: i32, backend: StarRowBackend) -> Result<Arc<dyn AcquiredStarRow>> {
            assert_eq!(backend, StarRowBackend::Native);
            let serial = self.next.fetch_add(1, Ordering::SeqCst) + 1;
            self.trace
                .lock()
                .unwrap()
                .push(format!("{}:select:{serial}", self.label));
            if self.refuse {
                return Err(Error::UnknownBody { body_id: 77 });
            }
            Ok(Arc::new(LiveRow {
                trace: Arc::clone(&self.trace),
                label: self.label,
                serial,
            }))
        }
        fn is_skippable(&self, _: &Error) -> Result<bool> {
            self.trace
                .lock()
                .unwrap()
                .push(format!("{}:catch", self.label));
            if self.refuse_classification {
                return Err(Error::ConfigurationMissing {
                    missing_config: "actual_catch_class",
                    suggestion: "same owned classifier refusal",
                });
            }
            Ok(self.refuse)
        }
    }

    #[test]
    fn live_native_scalar_reacquires_after_each_observer_and_releases_each_lease() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let rows = LiveRows {
            trace: Arc::clone(&trace),
            label: "a",
            next: AtomicUsize::new(0),
            refuse: false,
            refuse_classification: false,
        };
        let state_at = |id, _| {
            assert_eq!(id, EARTH);
            trace.lock().unwrap().push("observer".into());
            Ok(([0.0; 3], [0.0; 3]))
        };
        let raw = native_raw_place_with_rows_via(
            &state_at,
            &|_, _, error| error,
            90000001,
            2451545.0,
            true,
            true,
            true,
            0,
            None,
            true,
            Some(&rows),
        )
        .unwrap();
        assert!(raw[3].abs() > 1.0);
        assert_eq!(
            *trace.lock().unwrap(),
            [
                "observer",
                "a:select:1",
                "a:keywords:1",
                "a:drop:1",
                "observer",
                "a:select:2",
                "a:keywords:2",
                "a:drop:2",
                "observer",
                "a:select:3",
                "a:keywords:3",
                "a:drop:3",
            ]
        );
    }

    #[test]
    fn shared_batch_keeps_occurrence_owners_and_observers_before_all_leases() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let registry = crate::spk::SpkRegistry::new();
        let missing = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
        let provider = ContextStarProvider {
            inner: ContextStarBackend::Native(&missing),
        };
        let nutation = |_| Ok((0.0, 0.0));
        let delta = |_| Ok(0.0);
        let serving = SealedStarServing::new(
            &provider,
            &nutation,
            &delta,
            crate::sidereal::SiderealSelection::unset(),
        );
        let owners = ["a", "b"].map(|label| {
            Some(Arc::new(LiveRows {
                trace: Arc::clone(&trace),
                label,
                next: AtomicUsize::new(0),
                refuse: false,
                refuse_classification: false,
            }) as Arc<dyn StarRowAcquisition>)
        });
        let resolved = ["first", "second"].map(|name| {
            Some(StarMatch {
                star_id: 90000001,
                canonical: name.into(),
            })
        });
        let state_at = |id, _| {
            assert_eq!(id, EARTH);
            trace.lock().unwrap().push("observer".into());
            Ok(([0.0; 3], [0.0; 3]))
        };
        let output = native_resolved_batch_with_main_via(
            &|| {
                trace.lock().unwrap().push("main".into());
                Ok(())
            },
            Some(&|clock| {
                trace.lock().unwrap().push(format!("clock:{clock:?}"));
                Ok(())
            }),
            &state_at,
            &|_, _, error| error,
            &serving,
            &resolved,
            &owners,
            2451545.0,
            FLG_TRUEPOS | FLG_J2000 | FLG_SPEED,
            false,
        )
        .unwrap();
        assert_ne!(
            output[0].as_ref().unwrap().0[0],
            output[1].as_ref().unwrap().0[0]
        );
        assert_eq!(
            *trace.lock().unwrap(),
            [
                "main",
                "observer",
                "clock:SharedSpeedTimes",
                "clock:SharedSample(Previous)",
                "observer",
                "clock:SharedSample(Next)",
                "observer",
                "clock:SharedSample(Central)",
                "a:select:1",
                "a:keywords:1",
                "a:drop:1",
                "clock:SharedSample(Previous)",
                "a:select:2",
                "a:keywords:2",
                "a:drop:2",
                "clock:SharedSample(Next)",
                "a:select:3",
                "a:keywords:3",
                "a:drop:3",
                "clock:SharedSample(Central)",
                "b:select:1",
                "b:keywords:1",
                "b:drop:1",
                "clock:SharedSample(Previous)",
                "b:select:2",
                "b:keywords:2",
                "b:drop:2",
                "clock:SharedSample(Next)",
                "b:select:3",
                "b:keywords:3",
                "b:drop:3",
            ]
        );
    }

    #[test]
    fn shared_side_time_refusal_precedes_side_observers_and_keeps_selected_error_owner() {
        use std::cell::Cell;
        for refuse_preparation in [true, false] {
            let registry = crate::spk::SpkRegistry::new();
            let missing = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
            let provider = ContextStarProvider {
                inner: ContextStarBackend::Native(&missing),
            };
            let nutation = |_| panic!("observer refusal precedes frame");
            let delta = |_| panic!("retained original Time");
            let serving = SealedStarServing::new(
                &provider,
                &nutation,
                &delta,
                crate::sidereal::SiderealSelection::unset(),
            );
            let selected = Cell::new(StarObservationSample::Central);
            let states = Cell::new(0);
            let events = RefCell::new(Vec::new());
            let result = native_batch_samples_at_time_via(
                &|| {
                    events.borrow_mut().push("main");
                    Ok(())
                },
                Some(&|clock| match clock {
                    StarObservationClock::SharedSpeedTimes => {
                        events.borrow_mut().push("prepare");
                        if refuse_preparation {
                            Err(Error::UnknownBody { body_id: 71 })
                        } else {
                            Ok(())
                        }
                    }
                    StarObservationClock::SharedSample(sample) => {
                        selected.set(sample);
                        Ok(())
                    }
                    _ => panic!("only already acquired shared times are consumed"),
                }),
                &|_, _| {
                    states.set(states.get() + 1);
                    if selected.get() == StarObservationSample::Central {
                        events.borrow_mut().push("central");
                        Ok(([0.; 3], [0.; 3]))
                    } else {
                        assert_eq!(selected.get(), StarObservationSample::Previous);
                        events.borrow_mut().push("previous");
                        Err(Error::UnknownBody { body_id: 72 })
                    }
                },
                &|_, _, error| {
                    assert_eq!(selected.get(), StarObservationSample::Previous);
                    error
                },
                &serving,
                Vec::new(),
                Vec::new(),
                TwoPartEpoch {
                    hi: 2451545.,
                    lo: 0.125,
                },
                FLG_TRUEPOS | FLG_J2000 | FLG_SPEED,
                false,
                None,
            );
            assert!(
                matches!(result,Err(Error::UnknownBody { body_id }) if body_id==if refuse_preparation {71} else {72})
            );
            assert_eq!(states.get(), if refuse_preparation { 1 } else { 2 });
            assert_eq!(
                *events.borrow(),
                if refuse_preparation {
                    vec!["main", "central", "prepare"]
                } else {
                    vec!["main", "central", "prepare", "previous"]
                }
            );
        }
    }

    #[test]
    fn batch_classifier_refusal_is_terminal_even_when_skipping_is_disabled() {
        for skip in [false, true] {
            let trace = Arc::new(Mutex::new(Vec::new()));
            let registry = crate::spk::SpkRegistry::new();
            let missing = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
            let provider = ContextStarProvider {
                inner: ContextStarBackend::Native(&missing),
            };
            let nutation = |_| Ok((0.0, 0.0));
            let delta = |_| Ok(0.0);
            let serving = SealedStarServing::new(
                &provider,
                &nutation,
                &delta,
                crate::sidereal::SiderealSelection::unset(),
            );
            let rows = [Some(Arc::new(LiveRows {
                trace: Arc::clone(&trace),
                label: "a",
                next: AtomicUsize::new(0),
                refuse: true,
                refuse_classification: true,
            }) as Arc<dyn StarRowAcquisition>)];
            let resolved = [Some(StarMatch {
                star_id: 90000001,
                canonical: "actual".into(),
            })];
            let result = native_resolved_batch_with_main_via(
                &|| Ok(()),
                None,
                &|_, _| Ok(([0.0; 3], [0.0; 3])),
                &|_, _, error| error,
                &serving,
                &resolved,
                &rows,
                2451545.0,
                FLG_TRUEPOS | FLG_J2000,
                skip,
            );
            assert!(matches!(
                result,
                Err(Error::ConfigurationMissing {
                    missing_config: "actual_catch_class",
                    ..
                })
            ));
            assert_eq!(*trace.lock().unwrap(), ["a:select:1", "a:catch"]);
        }
    }

    #[test]
    fn skipped_vector_retry_preserves_an_actual_previously_stored_slot() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let registry = crate::spk::SpkRegistry::new();
        let missing = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
        let provider = ContextStarProvider {
            inner: ContextStarBackend::Native(&missing),
        };
        let nutation = |_| panic!("failed row must not reach frame");
        let delta = |_| panic!("already acquired original Time");
        let serving = SealedStarServing::new(
            &provider,
            &nutation,
            &delta,
            crate::sidereal::SiderealSelection::unset(),
        );
        let rows = [Some(Arc::new(LiveRows {
            trace: Arc::clone(&trace),
            label: "a",
            next: AtomicUsize::new(0),
            refuse: true,
            refuse_classification: false,
        }) as Arc<dyn StarRowAcquisition>)];
        let stored = (
            [11.0, 12.0, 13.0, 14.0, 15.0, 16.0],
            "actual earlier stored slot".into(),
            7,
        );
        let output = native_batch_samples_at_time_via(
            &|| Ok(()),
            None,
            &|_, _| panic!("BARYCTR has no provider state"),
            &|_, _, error| error,
            &serving,
            vec![(
                0,
                StarMatch {
                    star_id: 90000001,
                    canonical: "resolver name".into(),
                },
            )],
            vec![Some(stored.clone())],
            TwoPartEpoch {
                hi: 2451545.0,
                lo: 0.25,
            },
            FLG_BARYCTR | FLG_J2000,
            true,
            Some(&rows),
        )
        .unwrap();
        assert_eq!(output, vec![Some(stored)]);
        assert_eq!(*trace.lock().unwrap(), ["a:select:1", "a:catch"]);
    }

    #[test]
    fn noaberr_is_a_spherical_differential_with_geometric_distance() {
        let order = RefCell::new(Vec::new());
        let result = corrected_place_via(
            true,
            false,
            || {
                order.borrow_mut().push("geometric");
                Ok([359.0, 2.0, 7.0])
            },
            |default| {
                order
                    .borrow_mut()
                    .push(if default { "full" } else { "empty" });
                Ok(if default {
                    [1.0, 6.0, 90.0]
                } else {
                    [359.0, 5.0, 80.0]
                })
            },
        )
        .unwrap();
        assert_eq!(*order.borrow(), ["geometric", "full", "empty"]);
        assert_eq!(result, [1.0, 3.0, 7.0]);
    }

    #[test]
    fn correction_branches_are_lazy_and_empty_is_still_apparent() {
        assert_eq!(
            corrected_place_via(true, true, || Ok([1.0; 3]), |_| panic!("unused")),
            Ok([1.0; 3])
        );
        assert_eq!(
            corrected_place_via(
                false,
                true,
                || panic!("unused"),
                |default| {
                    assert!(!default);
                    Ok([2.0; 3])
                }
            ),
            Ok([2.0; 3])
        );
        let error = corrected_place_via(
            true,
            false,
            || Ok([1.0; 3]),
            |_| Err(Error::SpkBodyNotFound { target: 10 }),
        );
        assert!(matches!(error, Err(Error::SpkBodyNotFound { target: 10 })));
    }

    #[test]
    fn native_raw_unknown_star_and_missing_site_precede_image_reads() {
        let absent = |_, _| -> Result<([f64; 3], [f64; 3])> { panic!("no read") };
        assert!(matches!(
            native_raw_place_via(
                &absent,
                &|_, _, error| error,
                -999,
                2451545.0,
                true,
                true,
                true,
                0,
                None,
                false
            ),
            Err(Error::StarNotFound { .. })
        ));
        assert!(matches!(
            observer_via(&absent, scalar_tt(2451545.0).unwrap(), FLG_TOPOCTR, None),
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            })
        ));
        let origin =
            observer_via(&absent, scalar_tt(2451545.0).unwrap(), FLG_BARYCTR, None).unwrap();
        assert_eq!(origin.position, [0.0; 3]);
    }

    #[test]
    fn scalar_clock_parts_and_unused_site_guards_are_distinct() {
        assert!(matches!(
            scalar_tt(f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
        let scalar = scalar_tt(2451545.125).unwrap();
        assert_eq!((scalar.hi, scalar.lo), (2451545.0, 0.125));
        assert_eq!(public_site(0, Some((1.0, 91.0, -1e9))), Ok(None));
        assert!(matches!(
            public_site(FLG_TOPOCTR, None),
            Err(Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            })
        ));
        assert!(public_site(FLG_TOPOCTR, Some((1.0, f64::NAN, 0.0))).is_err());
    }

    #[test]
    fn batch_resolves_all_names_before_shared_observer_and_keeps_empty_misses() {
        let registry = crate::spk::SpkRegistry::new();
        let missing = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
        let provider = ContextStarProvider {
            inner: ContextStarBackend::Native(&missing),
        };
        let nutation = |_| Ok((0.0, 0.0));
        let delta = |_| Ok(0.0);
        let selection = crate::sidereal::SiderealSelection::default();
        let serving = SealedStarServing::new(&provider, &nutation, &delta, selection);
        let absent = |_, _| -> Result<([f64; 3], [f64; 3])> { panic!("must resolve first") };
        assert!(matches!(
            native_batch_via(
                &absent,
                &|_, _, error| error,
                &serving,
                selection,
                &["Regulus", "no-such-star-xyz"],
                2451545.125,
                2,
                None,
                false
            ),
            Err(Error::StarNotFound { .. })
        ));
        assert_eq!(
            native_batch_via(
                &absent,
                &|_, _, error| error,
                &serving,
                selection,
                &["no-such-star-xyz"],
                2451545.125,
                2,
                None,
                true
            )
            .unwrap(),
            vec![None]
        );
        assert_eq!(
            native_batch_via(
                &absent,
                &|_, _, error| error,
                &serving,
                selection,
                &["²"],
                2451545.125,
                2,
                None,
                true
            )
            .unwrap(),
            vec![None]
        );
        assert!(
            native_batch_via(
                &absent,
                &|_, _, error| error,
                &serving,
                selection,
                &[],
                2451545.125,
                2,
                None,
                true
            )
            .unwrap()
            .is_empty()
        );
    }

    #[test]
    fn batch_shared_observer_failure_is_not_a_per_star_skipped_miss() {
        let registry = crate::spk::SpkRegistry::new();
        let missing = LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
        let provider = ContextStarProvider {
            inner: ContextStarBackend::Native(&missing),
        };
        let nutation = |_| Ok((0.0, 0.0));
        let delta = |_| Ok(0.0);
        let selection = crate::sidereal::SiderealSelection::default();
        let serving = SealedStarServing::new(&provider, &nutation, &delta, selection);
        let failed = |id, _: TwoPartEpoch| -> Result<([f64; 3], [f64; 3])> {
            assert_eq!(id, EARTH);
            Err(Error::SpkBodyNotFound { target: 399 })
        };
        assert!(matches!(
            native_batch_via(
                &failed,
                &|_, _, error| error,
                &serving,
                selection,
                &["Regulus"],
                2451545.125,
                2,
                None,
                true
            ),
            Err(Error::SpkBodyNotFound { target: 399 })
        ));
    }

    #[test]
    fn source_batch_delegation_preserves_topo_and_projection_rules() {
        let base = crate::sidereal::SiderealSelection::default();
        assert!(scalar_batch_required(FLG_TOPOCTR, base));
        assert!(!scalar_batch_required(FLG_SIDEREAL, base));
        let mut fixed = base;
        fixed.set_sidereal_mode(18, 0.0, 0.0);
        assert!(scalar_batch_required(FLG_SIDEREAL, fixed));
    }
}

#[cfg(test)]
#[path = "star_native_anchor_tests.rs"]
mod native_anchor_tests;

#[cfg(test)]
mod live_anchor_ray_tests {
    use super::*;
    use crate::stars::AnchorRayRole;
    use std::cell::RefCell;

    struct Stored {
        absent: Option<i32>,
        fail: Option<i32>,
        reads: RefCell<Vec<(i32, f64)>>,
    }
    impl StateProvider for Stored {
        fn provider_name(&self) -> &'static str {
            "synthetic-stored-vector-anchor"
        }
        fn has_body(&self, body: i32) -> bool {
            self.absent != Some(body)
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            Some((2451544.0, 2451546.0))
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            Some(crate::leb::COORD_ICRS_BARY)
        }
        fn serving_entry_snapshot(&self, _: i32) -> Result<super::super::EntrySnapshot> {
            Ok(super::super::EntrySnapshot {
                coord_type: crate::leb::COORD_ICRS_BARY,
            })
        }
        fn eval_raw(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            self.reads.borrow_mut().push((body, jd));
            if self.fail == Some(body) {
                return Err(Error::LebOutOfRange {
                    body_id: body,
                    jd,
                    start: 2451544.0,
                    end: 2451546.0,
                });
            }
            let position = if body == EARTH {
                vec![0.0; 3]
            } else {
                vec![2.0 + f64::from(body), 1.0, 0.5]
            };
            Ok((position, vec![0.0; 3]))
        }
        fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("live LEB vector must consume scalar Time.value")
        }
    }

    #[test]
    fn actual_live_anchor_main_and_time_refusals_follow_helper_role() {
        let state = |_, _| panic!("observer must remain unused");
        let marker = Error::SpkBodyNotFound { target: 10 };
        for role in [
            AnchorRayRole::EclipticLongitude,
            AnchorRayRole::GalacticRightAscension,
        ] {
            let entered = RefCell::new(false);
            let main = || {
                *entered.borrow_mut() = true;
                Err(marker.clone())
            };
            let result = live_anchor_ray_via(
                &state,
                ([1.0; 3], [0.0; 3]),
                f64::NAN,
                true,
                true,
                role,
                main,
            );
            if role == AnchorRayRole::EclipticLongitude {
                assert_eq!(result, Err(marker.clone()));
                assert!(*entered.borrow());
            } else {
                assert!(matches!(result, Err(Error::NonFiniteArgument { .. })));
                assert!(!*entered.borrow());
            }
        }
    }

    #[test]
    fn actual_live_anchor_retains_constructor_and_original_earth_time() {
        let log = RefCell::new(Vec::new());
        let time = scalar_tt(2451545.125).unwrap();
        let state = |id, epoch| {
            assert_eq!(id, EARTH);
            assert_eq!(epoch, time);
            log.borrow_mut().push("earth");
            Ok(([0.0; 3], [0.0; 3]))
        };
        let initial = ([1.0, 2.0, 3.0], [0.01, 0.02, 0.03]);
        let main = || {
            log.borrow_mut().push("main");
            Ok(())
        };
        let ray = live_anchor_ray_via(
            &state,
            initial,
            time.value(),
            true,
            true,
            AnchorRayRole::EclipticLongitude,
            main,
        )
        .unwrap();
        let observer = NativeStarObserver {
            time,
            position: [0.0; 3],
            velocity: [0.0; 3],
            velocity_per_tt: [0.0; 3],
            surface: None,
        };
        assert_eq!(
            ray.map(f64::to_bits),
            astrometric_ray_from_initial(initial, observer)
                .0
                .map(f64::to_bits)
        );
        assert_eq!(*log.borrow(), ["main", "earth"]);
    }

    #[test]
    fn actual_stored_anchor_required_deflector_is_not_optional_raw_lookup() {
        for body in [SUN, crate::constants::JUPITER, crate::constants::SATURN] {
            let provider = Stored {
                absent: Some(body),
                fail: None,
                reads: RefCell::new(Vec::new()),
            };
            let state = |id, time| live_stored_anchor_state_via(&provider, id, time, &|e| e);
            let result = live_anchor_ray_via(
                &state,
                ([1.0, 2.0, 3.0], [0.0; 3]),
                2451545.125,
                false,
                false,
                AnchorRayRole::EclipticLongitude,
                || Ok(()),
            );
            let error = result.unwrap_err();
            assert!(
                matches!(error.primary_error(), Error::UnknownBody { body_id } if *body_id == body)
            );
            assert_eq!(provider.reads.borrow()[0].0, EARTH);
        }
    }
}
