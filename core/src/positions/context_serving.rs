//! Context-owned provider selection and synchronous serving bindings.
//!
//! Retained LEB dispatch shares immutable images and bounded caches. Call-scoped
//! factories bind that dispatch, clocks and explicit callbacks to the existing
//! position, orbital and phenomenon calculations; support-state and error
//! adapters retain each family's distinct physical and acquisition conventions.
//! Single-tier and unnamed images use attach order. Recognized multiple tiers
//! use `TieredIndex` so a narrow image cannot shadow wider date coverage.

use crate::context::{EphemerisContext, ProviderMode, Tier};
use crate::errors::{Error, Result};
use crate::leb::{CompositeIndex, LebFile, LebTier, TieredIndex, check_tier_subset, file_tier};
use crate::positions::{
    self, CalcEngine, CartesianWordParts, EntrySnapshot, PositionResult, SealedServing,
    ServingRequest, StateProvider, TierPolicy, TopoPlan, shift_state_word_parts,
};
use crate::time::{DeltaTConfig, DeltaTModel, deltat_with_config};
use crate::types::TwoPartEpoch;
use std::sync::{Arc, OnceLock};

/// Borrowed scalar center-offset query in split TDB, returning km and km/s
/// plus the actual selected descriptor. The callback never outlives its call.
type OffsetRead<'a> = dyn Fn(
        &crate::spk::SpkRegistry,
        i32,
        f64,
        f64,
    ) -> Result<(crate::spk::ChebState, crate::spk::Descriptor)>
    + 'a;

#[path = "context_serving/stack.rs"]
mod stack;
pub use stack::with_sealed_serving_stack;

#[path = "context_serving/predicted_rates.rs"]
mod predicted_rates;
pub(crate) use predicted_rates::model_position_with_owned_support;
pub(in crate::positions) use predicted_rates::pctr_state_with_owned_support;

#[cfg(test)]
#[path = "context_serving/hypothetical_support_observations.rs"]
pub(crate) mod hypothetical_support_observations;

/// Derived immutable dispatch shared by contexts with identical image settings.
/// It owns only image handles and metadata; coefficient storage stays shared.
#[derive(Default)]
pub(crate) struct LebProviderCache(OnceLock<Arc<RetainedLebProviders>>);

#[derive(Default)]
struct RetainedLebProviders {
    attach_order: OnceLock<CompositeIndex<'static>>,
    tier_priority: OnceLock<CompositeIndex<'static>>,
    date_aware: OnceLock<BoundProvider<'static>>,
}

impl Clone for LebProviderCache {
    fn clone(&self) -> Self {
        // Forks must share even an empty set of lazy slots. Copying an empty
        // OnceLock would rebuild the same image indexes independently later.
        // Allocating this owner reads no image and preserves lazy validation.
        Self(OnceLock::from(self.0.get_or_init(Default::default).clone()))
    }
}

impl std::fmt::Debug for LebProviderCache {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("LebProviderCache")
            .field("initialized", &self.0.get().is_some())
            .finish()
    }
}

impl PartialEq for LebProviderCache {
    fn eq(&self, _other: &Self) -> bool {
        // Derived memoization is not part of a context's physical settings.
        true
    }
}

impl LebProviderCache {
    fn providers(&self) -> &RetainedLebProviders {
        self.0.get_or_init(Default::default)
    }

    /// Retain first-wins dispatch with the requested image ordering.
    /// File, tier-label and subset errors keep their initial validation order;
    /// a cached view still enforces the caller's closed-handle policy.
    pub(crate) fn composite(
        &self,
        ctx: &EphemerisContext,
        policy: TierPolicy,
        skip_closed: bool,
    ) -> Result<&CompositeIndex<'static>> {
        let slot = match policy {
            TierPolicy::AttachOrder => &self.providers().attach_order,
            TierPolicy::TierPriority => &self.providers().tier_priority,
        };
        if let Some(composite) = slot.get() {
            if !skip_closed && ctx.leb_files().iter().any(LebFile::is_closed) {
                return Err(Error::LebReaderClosed);
            }
            return Ok(composite);
        }
        // Initial construction retains file/tier validation order. A missing
        // tier label can precede a later closed handle's refusal.
        {
            let mut files: Vec<&LebFile> = Vec::with_capacity(ctx.leb_files().len());
            let mut installed: Vec<LebTier> = Vec::with_capacity(ctx.leb_files().len());
            for file in ctx.leb_files() {
                if file.is_closed() {
                    if skip_closed {
                        continue;
                    }
                    return Err(Error::LebReaderClosed);
                }
                if policy == TierPolicy::TierPriority {
                    match super::routing::tier_token_of(file.path()) {
                        Some(tier) => installed.push(tier),
                        None => {
                            return Err(Error::ConfigurationMissing {
                                missing_config: "tier_priority_dispatch",
                                suggestion: "tier priority needs a recognized tier token on every bound image",
                            });
                        }
                    }
                }
                files.push(file);
            }
            if policy == TierPolicy::TierPriority {
                check_tier_subset(&installed, super::routing::active_leb_tier(ctx.tier()))?;
                let mut order: Vec<usize> = (0..files.len()).collect();
                order.sort_by_key(|&i| {
                    LebTier::PRIORITY
                        .iter()
                        .position(|t| *t == installed[i])
                        .unwrap_or(usize::MAX)
                });
                files = order.into_iter().map(|i| files[i]).collect();
            }
            let composite = CompositeIndex::from_shared_files(files)?;
            let _ = slot.set(composite);
        }
        Ok(slot.get().expect("initialized LEB composite"))
    }

    /// Retain date-aware dispatch and the context's optional center images.
    /// Closed/invalid images and ineligible tiers fail without caching a failure.
    pub(crate) fn provider(
        &self,
        ctx: &EphemerisContext,
        skip_closed: bool,
    ) -> Result<&BoundProvider<'static>> {
        // A view initialized with skip_closed must not hide a later strict
        // request's refusal. Attached handles are otherwise immutable.
        if !skip_closed && ctx.leb_files().iter().any(LebFile::is_closed) {
            return Err(Error::LebReaderClosed);
        }
        let slot = &self.providers().date_aware;
        if slot.get().is_none() {
            let provider = provider_from_images(
                ctx.leb_files().iter().map(|file| (file.path(), file)),
                ctx.tier(),
                skip_closed,
            )?;
            let provider = if ctx.mode() == ProviderMode::Auto && !ctx.spk_registry().is_empty() {
                BoundProvider::WithCenters(Box::new(provider), ctx.spk_registry().clone())
            } else {
                provider
            };
            // Racing initializers build the same immutable dispatch. Only one
            // is retained; failures are never cached or substituted.
            let _ = slot.set(provider);
        }
        Ok(slot.get().expect("initialized LEB provider"))
    }
}

#[cfg(test)]
mod lazy_provider_tests {
    use super::*;

    #[test]
    fn cold_clones_share_lazy_immutable_indexes() {
        let original = LebProviderCache::default();
        let fork = original.clone();
        let original = original.providers();
        let fork = fork.providers();
        assert!(std::ptr::eq(original, fork));
        assert!(!std::ptr::eq(
            original,
            LebProviderCache::default().providers()
        ));
        // Sharing the preparation owner must not load or validate any image.
        assert!(original.attach_order.get().is_none());
        assert!(original.tier_priority.get().is_none());
        assert!(original.date_aware.get().is_none());
    }
}

/// Apparent unrefracted altitude and azimuth of an arbitrary ICRS star.
///
/// `jd_ut1` is the original finite UT1 Julian day, independent of caller
/// Delta-T/model settings. RA/declination and east longitude/latitude are
/// finite degrees, with a WGS84 site at zero metres elevation. Finite angles
/// are used directly without a pole clamp. The star is fixed at J2000 with
/// zero proper motion/radial velocity/parallax and the 1e-6 mas stellar
/// floor. One owned NativeJpl context supplies Earth and strict apparent
/// deflector states at original split TT/TDB. No scalar public TT envelope
/// replaces raw provider directory checks. The bundled default global Time
/// supplies TT and original UT1 fractions, with no polar-motion table.
///
/// Return `(altitude_deg, azimuth_deg)`, altitude above the horizon and
/// azimuth north through east in [0,360). Construct the WGS84 site vector,
/// observer position/velocity and apparent ICRS ray, then rotate with
/// `R_latlon @ R_ITRS`, where `R_latlon = reverse_rows(R_y(latitude)) @
/// R_z(-longitude)`. No refraction or configured observer is consumed.
/// Nonfinite inputs, missing NativeJpl context/aliases, frames and actual
/// provider range/type failures retain typed errors in acquisition order.
pub fn native_jpl_arbitrary_star_altaz_ut1(
    ctx: &EphemerisContext,
    jd_ut1: f64,
    ra_deg: f64,
    dec_deg: f64,
    lon_deg: f64,
    lat_deg: f64,
) -> Result<(f64, f64)> {
    use positions::reference_frames::*;
    let (epochs, tt) = VectorRotationEpochs::from_ut1(jd_ut1)?;
    positions::forbid_jpl_in_sealed_mode(ctx.mode())?;
    if ctx.mode() != ProviderMode::NativeJpl {
        return Err(Error::ConfigurationMissing {
            missing_config: "native_star_context",
            suggestion: "select NativeJpl and register the intended local SPK image",
        });
    }
    let provider = positions::LocalSpkStateProvider::new(ctx.spk_registry(), ctx.mode())?;
    provider.validate_source_alias(crate::constants::EARTH)?;
    for value in [ra_deg, dec_deg, lon_deg, lat_deg] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let site =
        native_pheno_itrs_position((lon_deg, lat_deg, 0.0), positions::SiteEllipsoid::Wgs84)?;
    let phi = lat_deg / 360.0 * FRAME_TAU;
    let zero = phi * 0.0;
    let one = zero + 1.0;
    let lat = [
        [-scalar_sin(phi), zero, scalar_cos(phi)],
        [zero, one, zero],
        [scalar_cos(phi), zero, scalar_sin(phi)],
    ];
    let latlon = matmul3(lat, rot_z(-(lon_deg / 360.0 * FRAME_TAU)))?;
    let time = TwoPartEpoch {
        hi: epochs.whole,
        lo: epochs.tt_fraction,
    };
    let earth = provider.geometric_state_two_part(crate::constants::EARTH, time)?;
    let rotation =
        itrs_rotation_at_parts(epochs.whole, epochs.tt_fraction, epochs.ut1_fraction, tt)?;
    let (surface, site_velocity) = itrs_to_gcrs(rotation, site)?;
    let mut observer = earth.pos();
    let mut velocity = earth.vel();
    for i in 0..3 {
        observer[i] += surface[i];
        velocity[i] += site_velocity[i];
    }
    let ray = provider
        .arbitrary_star_apparent_vector(time, ra_deg, dec_deg, observer, velocity, surface)?;
    let horizon = matmul3(latlon, rotation)?;
    let (azimuth, altitude) = native_pheno_frame_angles(horizon, ray)?;
    Ok((altitude, azimuth))
}

/// Acquire and own the actual native eclipse main image at its consumed
/// boundary. Supplied registries retain their mappings without a file probe;
/// otherwise the existing factory registers all selected targets. The caller
/// retains this configuration through its scalar/vector acquisitions. No clock,
/// observer, astronomical sample, process setting or callback is consumed here.
pub(crate) fn native_eclipse_main_context(ctx: &EphemerisContext) -> Result<EphemerisContext> {
    if ctx.mode() != ProviderMode::NativeJpl {
        return Err(Error::ConfigurationMissing {
            missing_config: "native_eclipse_context",
            suggestion: "select NativeJpl with the intended main image",
        });
    }
    if ctx.spk_registry().is_empty() {
        Ok(positions::main_image::acquire_main_image(ctx)?.context)
    } else {
        Ok(ctx.clone())
    }
}

/// Acquire the configured main image for an analytic lunar range query.
/// Forced LEB retains its sealed-mode refusal before loader/policy acquisition,
/// even after a reader bootstrap miss. Other modes retain their actual policy
/// and selected images. No coefficient is evaluated or caller mode mutated.
pub(crate) fn mean_lunar_main_context(ctx: &EphemerisContext) -> Result<EphemerisContext> {
    positions::forbid_jpl_in_sealed_mode(ctx.mode())?;
    Ok(positions::main_image::acquire_missing_main_image(ctx)?.context)
}

/// Apparent Sun then Moon ecliptic-of-date places for lunar circumstances.
///
/// `jd_ut1` is the finite original UT1 scalar, `site` is east longitude and
/// latitude in degrees plus WGS84 elevation in metres. Main-image selection
/// precedes bundled global UT1 Time; Earth alias precedes WGS84 construction,
/// Earth state precedes independent ITRS rotation. TT and original UT1
/// fractions are retained separately; no inverse Delta-T or configured model
/// clock replaces them. Add the site AU/AU-day segments to the SSB Earth state,
/// then observe Sun and Moon in that order with strict apparent corrections.
/// Frame nutation consumes scalar TT; mean obliquity/precession consumes TDB.
///
/// Return `(sun, moon)`, each `[longitude_deg, latitude_deg, distance_AU]`.
/// Longitude is in [0,360); distance is the left-associated norm of the actual
/// frame-rotated ray. Finite sites retain ellipsoid arithmetic without a new
/// latitude bound. Missing aliases, range/type/directory/frame and ten-step
/// light-time failures retain raw acquisition order and identity. This channel
/// applies no public scalar-TT envelope or caller correction/output flags.
/// The owning context/snapshot and all provider borrows end with this call.
/// Explicit LEB/Auto does not enter this native helper.
pub(crate) fn native_lunar_eclipse_topocentric_pair(
    ctx: &EphemerisContext,
    jd_ut1: f64,
    site: [f64; 3],
) -> Result<([f64; 3], [f64; 3])> {
    use positions::reference_frames::*;
    if ctx.mode() != ProviderMode::NativeJpl {
        return Err(Error::ConfigurationMissing {
            missing_config: "native_lunar_eclipse_context",
            suggestion: "select NativeJpl with the intended main-image configuration",
        });
    }
    // The main image is acquired before constructing the independent clock.
    let acquired = if ctx.spk_registry().is_empty() {
        Some(positions::main_image::acquire_main_image(ctx)?)
    } else {
        None
    };
    let context = acquired.as_ref().map_or(ctx, |owned| &owned.context);
    let (epochs, tt) = VectorRotationEpochs::from_ut1(jd_ut1)?;
    let provider = positions::LocalSpkStateProvider::new(context.spk_registry(), context.mode())?;
    provider.validate_source_alias(crate::constants::EARTH)?;
    let itrs =
        native_pheno_itrs_position((site[0], site[1], site[2]), positions::SiteEllipsoid::Wgs84)?;
    let time = TwoPartEpoch {
        hi: epochs.whole,
        lo: epochs.tt_fraction,
    };
    let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
    let rotation =
        itrs_rotation_at_parts(epochs.whole, epochs.tt_fraction, epochs.ut1_fraction, tt)?;
    let (surface, site_velocity) = itrs_to_gcrs(rotation, itrs)?;
    let observer = std::array::from_fn(|i| earth.0[i] + surface[i]);
    let velocity = std::array::from_fn(|i| earth.1[i] + site_velocity[i]);
    let place = |body| -> Result<[f64; 3]> {
        let ray = provider.eclipse_source_observed_vector(
            body,
            time,
            observer,
            velocity,
            Some(surface),
        )?;
        let frame = native_pheno_ecliptic_frame(time)?;
        let (lon, lat) = native_pheno_frame_angles(frame, ray)?;
        let distance = native_pheno_vector_length(matvec3(frame, ray)?);
        Ok([lon, lat, distance])
    };
    let sun = place(crate::constants::SUN)?;
    let moon = place(crate::constants::MOON)?;
    Ok((sun, moon))
}

/// Apparent geocentric Sun then Moon ICRS rays in AU at original UT1.
/// Main selection precedes independent bundled Time; Earth state precedes both
/// light-time observations. No surface vector is installed, so Earth
/// deflection is omitted. Actual provider directory/range/aberration and
/// per-segment units are retained. No observer, scalar public envelope, output
/// flag, equatorial/ecliptic reconstruction or configured model clock is used.
pub(crate) fn native_solar_gamma_pair(
    ctx: &EphemerisContext,
    jd_ut1: f64,
) -> Result<([f64; 3], [f64; 3])> {
    use positions::reference_frames::VectorRotationEpochs;
    if ctx.mode() != ProviderMode::NativeJpl {
        return Err(Error::ConfigurationMissing {
            missing_config: "native_solar_gamma_context",
            suggestion: "select NativeJpl with the intended main-image configuration",
        });
    }
    let acquired = if ctx.spk_registry().is_empty() {
        Some(positions::main_image::acquire_main_image(ctx)?)
    } else {
        None
    };
    let context = acquired.as_ref().map_or(ctx, |owned| &owned.context);
    let (epochs, _) = VectorRotationEpochs::from_ut1(jd_ut1)?;
    let provider = positions::LocalSpkStateProvider::new(context.spk_registry(), context.mode())?;
    provider.validate_source_alias(crate::constants::EARTH)?;
    let time = TwoPartEpoch {
        hi: epochs.whole,
        lo: epochs.tt_fraction,
    };
    let earth = provider.star_source_time_state(crate::constants::EARTH, time)?;
    let sun = provider.eclipse_source_observed_vector(
        crate::constants::SUN,
        time,
        earth.0,
        earth.1,
        None,
    )?;
    let moon = provider.eclipse_source_observed_vector(
        crate::constants::MOON,
        time,
        earth.0,
        earth.1,
        None,
    )?;
    Ok((sun, moon))
}

/// Real provider view over immutable images, retaining legacy single-tier
/// attach order or the producer's cumulative date-aware selection.
pub enum BoundProvider<'a> {
    /// Existing first-wins dispatch within one tier/unnamed explicit input.
    Composite(CompositeIndex<'a>),
    /// Image priority and per-evaluation TT date selection across tiers.
    Tiered(TieredIndex<'a>),
    /// The same stored index with owned optional planet-center SPK images.
    WithCenters(Box<BoundProvider<'a>>, crate::spk::SpkRegistry),
}

impl BoundProvider<'_> {
    /// Whether any owned coefficient group contains a nutation series.
    pub fn has_nutation(&self) -> bool {
        match self {
            Self::Composite(index) => index.has_nutation(),
            Self::Tiered(index) => index.has_nutation(),
            Self::WithCenters(inner, _) => inner.has_nutation(),
        }
    }
    /// Report a real stored-body coverage failure with its original reader cause.
    /// An in-range support failure is not relabeled as missing requested coverage.
    pub(crate) fn public_calc_range_error(
        &self,
        error: &Error,
        requested_body: i32,
        epoch: f64,
    ) -> Option<Error> {
        if !matches!(error.primary_error(), Error::LebOutOfRange { .. })
            || !crate::positions::is_leb_core_body(requested_body)
        {
            return None;
        }
        let (start_jd, end_jd) = self.body_coverage(requested_body)?;
        if !(epoch < start_jd || epoch > end_jd) {
            return None;
        }
        Some(
            Error::EphemerisRange {
                requested_jd: epoch,
                start_jd,
                end_jd,
                body: Some(requested_body),
            }
            .with_cause(error.clone()),
        )
    }

    /// Report an exterior stored-vector miss with this provider's union coverage.
    /// An interior tier gap retains the reader's precise missing interval.
    pub(crate) fn vector_error(&self, error: Error) -> Error {
        vector_reader_error(error, |body| self.body_coverage(body))
    }

    /// Stored Delta-T in days at scalar TT, over this acquired reader view.
    /// Composite readers retain their first table carrier; tiered readers
    /// retain date-aware table priority and fall-through. Table absence,
    /// malformed data and exhausted coverage follow each view's existing
    /// lookup policy; a final failure is returned unchanged.
    /// This query does not use a configured model or another image source.
    pub(crate) fn delta_t(&self, jd_tt: f64) -> Result<f64> {
        match self {
            Self::Composite(index) => index.delta_t(jd_tt),
            Self::Tiered(index) => index.delta_t(jd_tt),
            Self::WithCenters(index, _) => index.delta_t(jd_tt),
        }
    }

    /// Stored nutation in radians at TT JD, selected by the same provider.
    fn nutation_at(&self, jd: f64) -> Result<(f64, f64)> {
        match self {
            Self::Composite(index) => index.eval_nutation(jd),
            Self::Tiered(index) => index.eval_nutation(jd),
            Self::WithCenters(index, _) => index.nutation_at(jd),
        }
    }
}

impl BoundProvider<'_> {
    /// Shared scalar composition. The borrowed read changes only call-scoped
    /// ownership observation; clock conversion, optional absence and errors are
    /// identical to the ordinary provider path.
    fn center_offset_using(
        &self,
        body: i32,
        jd_tt: f64,
        read: &OffsetRead<'_>,
    ) -> Result<Option<[f64; 3]>> {
        if !self
            .serving_coord_type(body)
            .is_some_and(crate::leb::is_system_barycenter_channel)
        {
            return Ok(None);
        }
        let Self::WithCenters(_, registry) = self else {
            return Ok(None);
        };
        let (target, parent) = match body {
            crate::constants::JUPITER => (599, 5),
            crate::constants::SATURN => (699, 6),
            crate::constants::URANUS => (799, 7),
            crate::constants::NEPTUNE => (899, 8),
            crate::constants::PLUTO => (999, 9),
            _ => return Ok(None),
        };
        if registry.body_info(target).is_none() {
            return Ok(None);
        }
        let whole = jd_tt.floor();
        let fraction = crate::time::tt_to_tdb_fraction(whole, jd_tt - whole)?;
        let result = read(registry, target, whole, fraction);
        #[cfg(test)]
        hypothetical_support_observations::center_read(
            self, registry, body, target, jd_tt, whole, fraction, &result,
        );
        match result {
            Ok((_, descriptor)) if descriptor.frame != 1 => Err(Error::SpkUnsupportedFrame {
                frame: descriptor.frame,
            }),
            Ok((state, descriptor)) if descriptor.center == parent => Ok(Some(state.position_au())),
            Ok(_) => Err(Error::SpkCorrupt {
                what: "planet-center offset has the wrong system barycenter".to_owned(),
                offset: 0,
            }),
            Err(Error::SpkOutOfRange { .. }) | Err(Error::SpkBodyNotFound { .. }) => Ok(None),
            Err(error) => Err(error),
        }
    }
}

impl StateProvider for BoundProvider<'_> {
    fn has_analytic_kinematics(&self) -> bool {
        match self {
            Self::Composite(index) => index.has_analytic_kinematics(),
            Self::Tiered(index) => index.has_analytic_kinematics(),
            // A valid Earth acceleration does not qualify the SPK offsets of
            // a target or deflector. Their composed derivatives remain sampled.
            Self::WithCenters(_, _) => false,
        }
    }

    fn geometric_state_word_parts_two_part(
        &self,
        body: i32,
        epoch: TwoPartEpoch,
    ) -> Result<CartesianWordParts> {
        match self {
            Self::Composite(index) => index.geometric_state_word_parts_two_part(body, epoch),
            Self::Tiered(index) => index.geometric_state_word_parts_two_part(body, epoch),
            Self::WithCenters(index, _) => {
                let words = index.geometric_state_word_parts_two_part(body, epoch)?;
                let previous = index.geometric_state_two_part(body, epoch)?;
                let current = self.geometric_state_two_part(body, epoch)?;
                Ok(shift_state_word_parts(
                    words,
                    (previous.pos(), previous.vel()),
                    (current.pos(), current.vel()),
                ))
            }
        }
    }

    fn geometric_acceleration_tt(&self, body: i32, jd: f64) -> Result<Option<[f64; 3]>> {
        match self {
            Self::Composite(index) => index.geometric_acceleration_tt(body, jd),
            Self::Tiered(index) => index.geometric_acceleration_tt(body, jd),
            Self::WithCenters(index, _) => {
                // A stored barycenter is not the physical center once an
                // offset is applied. Until that source supplies acceleration,
                // retain the established complete calculation for this case.
                if self.body_center_offset_at(body, jd)?.is_some() {
                    return Ok(None);
                }
                index.geometric_acceleration_tt(body, jd)
            }
        }
    }

    fn can_reuse_scalar_state(&self, body: i32, jd: f64, epoch: TwoPartEpoch) -> bool {
        match self {
            Self::Composite(index) => index.can_reuse_scalar_state(body, jd, epoch),
            Self::Tiered(index) => index.can_reuse_scalar_state(body, jd, epoch),
            // A raw LEB argument proof does not establish equivalence of the
            // additional SPK center composition or its separate error route.
            Self::WithCenters(_, _) => false,
        }
    }

    fn provider_name(&self) -> &'static str {
        match self {
            Self::Composite(index) => index.provider_name(),
            Self::Tiered(index) => index.provider_name(),
            Self::WithCenters(index, _) => index.provider_name(),
        }
    }
    fn has_body(&self, id: i32) -> bool {
        match self {
            Self::Composite(index) => index.has_body(id),
            Self::Tiered(index) => index.has_body(id),
            Self::WithCenters(index, _) => index.has_body(id),
        }
    }
    fn body_coverage(&self, id: i32) -> Option<(f64, f64)> {
        match self {
            Self::Composite(index) => index.body_coverage(id),
            Self::Tiered(index) => index.body_coverage(id),
            Self::WithCenters(index, _) => index.body_coverage(id),
        }
    }
    fn smooth_interval_tt(&self, id: i32, jd: f64) -> Result<Option<(f64, f64)>> {
        match self {
            Self::Composite(index) => index.smooth_interval_tt(id, jd),
            Self::Tiered(index) => index.smooth_interval_tt(id, jd),
            Self::WithCenters(index, _) => index.smooth_interval_tt(id, jd),
        }
    }
    fn nutation_smooth_interval_tt(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        match self {
            Self::Composite(index) => index.nutation_smooth_interval_tt(jd),
            Self::Tiered(index) => index.nutation_smooth_interval_tt(jd),
            Self::WithCenters(index, _) => index.nutation_smooth_interval_tt(jd),
        }
    }
    fn nutation_rate_tt(&self, jd: f64) -> Result<Option<(f64, f64)>> {
        match self {
            Self::Composite(index) => index.nutation_rate(jd),
            Self::Tiered(index) => index.nutation_rate(jd),
            Self::WithCenters(index, _) => index.nutation_rate_tt(jd),
        }
    }
    fn serving_coord_type(&self, id: i32) -> Option<u32> {
        match self {
            Self::Composite(index) => index.serving_coord_type(id),
            Self::Tiered(index) => index.serving_coord_type(id),
            Self::WithCenters(index, _) => index.serving_coord_type(id),
        }
    }
    fn eval_raw_inline(&self, id: i32, jd: f64) -> Result<Option<crate::leb::CartesianState>> {
        let result = match self {
            Self::Composite(index) => index.eval_raw_inline(id, jd),
            Self::Tiered(index) => index.eval_raw_inline(id, jd),
            Self::WithCenters(index, _) => index.eval_raw_inline(id, jd),
        };
        #[cfg(test)]
        hypothetical_support_observations::leb_read(self, id, jd, &result);
        result
    }

    fn eval_raw_inline_two_part(
        &self,
        id: i32,
        epoch: TwoPartEpoch,
    ) -> Result<Option<crate::leb::CartesianState>> {
        match self {
            Self::Composite(index) => index.eval_raw_inline_two_part(id, epoch),
            Self::Tiered(index) => index.eval_raw_inline_two_part(id, epoch),
            Self::WithCenters(index, _) => index.eval_raw_inline_two_part(id, epoch),
        }
    }

    fn eval_raw(&self, id: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
        match self {
            Self::Composite(index) => index.eval_raw(id, jd),
            Self::Tiered(index) => index.eval_raw(id, jd),
            Self::WithCenters(index, _) => index.eval_raw(id, jd),
        }
    }
    fn eval_relative_position_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        match self {
            Self::Composite(index) => {
                index.eval_relative_position_two_part(target, target_epoch, center, center_epoch)
            }
            Self::Tiered(index) => {
                index.eval_relative_position_two_part(target, target_epoch, center, center_epoch)
            }
            Self::WithCenters(index, _) => {
                index.eval_relative_position_two_part(target, target_epoch, center, center_epoch)
            }
        }
    }
    fn eval_relative_velocity_two_part(
        &self,
        target: i32,
        target_epoch: TwoPartEpoch,
        center: i32,
        center_epoch: TwoPartEpoch,
    ) -> Result<Option<[f64; 3]>> {
        match self {
            Self::Composite(index) => {
                index.eval_relative_velocity_two_part(target, target_epoch, center, center_epoch)
            }
            Self::Tiered(index) => {
                index.eval_relative_velocity_two_part(target, target_epoch, center, center_epoch)
            }
            Self::WithCenters(index, _) => {
                index.eval_relative_velocity_two_part(target, target_epoch, center, center_epoch)
            }
        }
    }

    fn eval_relative_icrs(
        &self,
        target: i32,
        center: i32,
        jd_tt: f64,
    ) -> Result<Option<([f64; 3], [f64; 3])>> {
        match self {
            Self::Composite(index) => index.eval_relative_icrs(target, center, jd_tt),
            Self::Tiered(index) => index.eval_relative_icrs(target, center, jd_tt),
            Self::WithCenters(index, _) => index.eval_relative_icrs(target, center, jd_tt),
        }
    }
    fn eval_raw_two_part(&self, id: i32, epoch: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
        match self {
            Self::Composite(index) => index.eval_raw_two_part(id, epoch),
            Self::Tiered(index) => index.eval_raw_two_part(id, epoch),
            Self::WithCenters(index, _) => index.eval_raw_two_part(id, epoch),
        }
    }
    fn selected_tier(&self, id: i32, jd: f64) -> Option<LebTier> {
        match self {
            Self::Composite(_) => None,
            Self::Tiered(index) => index.selected_tier(id, jd),
            Self::WithCenters(index, _) => index.selected_tier(id, jd),
        }
    }
    fn body_center_offset_at(&self, body: i32, jd_tt: f64) -> Result<Option<[f64; 3]>> {
        self.center_offset_using(body, jd_tt, &|registry, target, whole, fraction| {
            registry.eval_body_relative_two_part_jd(target, whole, fraction)
        })
    }

    fn serving_entry_snapshot(&self, id: i32) -> Result<EntrySnapshot> {
        match self {
            Self::Composite(index) => index.serving_entry_snapshot(id),
            Self::Tiered(index) => index.serving_entry_snapshot(id),
            Self::WithCenters(index, _) => index.serving_entry_snapshot(id),
        }
    }
}

/// Convert actual context tier to the producer's reader-tier vocabulary.
fn active_tier(tier: Tier) -> LebTier {
    match tier {
        Tier::Base => LebTier::Base,
        Tier::Medium => LebTier::Medium,
        Tier::Extended => LebTier::Extended,
    }
}

/// Build an index over shared immutable image owners. Multi-tier input requires every
/// image to name its tier and the installed set to be eligible under `active`.
/// Single-tier explicit images retain attach order/unknown-label behavior.
/// Required closed/invalid containers fail; `skip_closed` is the explicit SDK
/// policy, not automatic data substitution. No image bytes are copied here.
pub fn provider_from_images<'a>(
    images: impl Iterator<Item = (&'a str, &'a LebFile)>,
    active: Tier,
    skip_closed: bool,
) -> Result<BoundProvider<'static>> {
    let mut files = Vec::new();
    let mut tiers = Vec::new();
    for (path, image) in images {
        if image.is_closed() {
            if skip_closed {
                continue;
            }
            return Err(Error::LebReaderClosed);
        }
        files.push(image);
        tiers.push(file_tier(path.rsplit(['/', '\\']).next().unwrap_or(path)));
    }
    let distinct: Vec<_> = LebTier::PRIORITY
        .into_iter()
        .filter(|tier| tiers.contains(&Some(*tier)))
        .collect();
    if distinct.len() <= 1 {
        return CompositeIndex::from_shared_files(files).map(BoundProvider::Composite);
    }
    let installed: Vec<_> =
        tiers
            .iter()
            .copied()
            .collect::<Option<Vec<_>>>()
            .ok_or(Error::ConfigurationMissing {
                missing_config: "tier_priority_dispatch",
                suggestion: "every image in a cumulative dataset must have a tier label",
            })?;
    check_tier_subset(&installed, active_tier(active))?;
    let mut composites = Vec::new();
    for tier in distinct {
        let group = files
            .iter()
            .zip(&installed)
            .filter(|(_, t)| **t == tier)
            .map(|(file, _)| *file);
        composites.push((tier, CompositeIndex::from_shared_files(group)?));
    }
    TieredIndex::new(composites).map(BoundProvider::Tiered)
}

/// Borrowed serving view. All numerical processing stays in the published
/// producer's generic SealedServing, with this real provider and its nutation.
pub struct ServingStack<'s, 'c> {
    /// Actual explicitly configured provider mode.
    mode: ProviderMode,
    /// Public producer serving, borrowing the selected index and nutation feed.
    pub serving: &'s SealedServing<'s, BoundProvider<'c>>,
    /// TT JD → nutation longitude/obliquity radians, selected stored carrier.
    pub nutation_at: &'s dyn Fn(f64) -> Result<(f64, f64)>,
    /// Whether the angle callback is the retained provider's stored series.
    stored_nutation_feed: bool,
    /// Actual context time-scale policy; UT input is paired with TT by this SDK.
    pub delta_cfg: &'s DeltaTConfig<'c>,
}

impl<'s, 'c> ServingStack<'s, 'c> {
    /// Alternate derivative serving over the same acquired source/feed/clock.
    /// This borrows stack storage and performs no reads or acquisitions.
    pub(crate) fn with_derivative_serving<'a>(
        &'a self,
        serving: &'a SealedServing<'a, BoundProvider<'c>>,
    ) -> ServingStack<'a, 'c> {
        ServingStack {
            mode: self.mode,
            serving,
            nutation_at: self.nutation_at,
            stored_nutation_feed: self.stored_nutation_feed,
            delta_cfg: self.delta_cfg,
        }
    }

    /// Borrow the same provider and time policy with a selected TT frame feed.
    /// The callback returns nutation longitude and obliquity in radians and
    /// propagates its frame-domain errors. The view retains the checked provider
    /// mode and cannot outlive either the original stack or the replacement feed.
    pub(crate) fn with_nutation<'a>(
        &'a self,
        nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    ) -> ServingStack<'a, 'c> {
        ServingStack {
            mode: self.mode,
            serving: self.serving,
            nutation_at,
            stored_nutation_feed: false,
            delta_cfg: self.delta_cfg,
        }
    }

    /// Longitude-nutation derivative of the actual consumed feed, degrees/day.
    /// Stored angles use their same owned polynomial's analytical TT derivative;
    /// an explicit replacement differentiates that callback, without consulting
    /// the stored rate or substituting a different astronomical model. Missing
    /// data, invalid metadata, sample/domain and nonfinite errors propagate.
    fn nutation_longitude_rate_deg_day(&self, jd_tt: f64) -> Result<f64> {
        if self.stored_nutation_feed {
            if let Some((rate, _)) = self.serving.provider.nutation_rate_tt(jd_tt)? {
                if !rate.is_finite() {
                    return Err(Error::NonFiniteArgument { value: rate });
                }
                return Ok(rate.to_degrees());
            }
        }
        // This path also preserves the actual value error when a stored
        // provider has no derivative capability at the requested epoch.
        (self.nutation_at)(jd_tt)?;
        let span = if self.stored_nutation_feed {
            self.serving.provider.nutation_smooth_interval_tt(jd_tt)?
        } else {
            None
        };
        let support = span.map_or((f64::NEG_INFINITY, f64::INFINITY), |(start, end)| {
            (start - jd_tt, end - jd_tt)
        });
        let sample = |offset| {
            (self.nutation_at)(jd_tt + offset).map(|(longitude, _)| [longitude.to_degrees(), 0.0])
        };
        super::derivative::angular_rates(&sample, 1.0 / 64.0, support).map(|rate| rate[0])
    }
}

/// Osculating orbital elements over an owned native JPL/SPK context.
///
/// `epoch` is a finite Julian Day in UT1 when `entry_is_ut1`, otherwise TT.
/// The bundled timescale maps UT1 observations to split TT/TDB; configured
/// Delta-T and IERS settings remain attached to nested analytical-model calls.
/// Raw ICRS states are AU and AU/day. Earth's osculating reduction uses NAIF 3
/// (EMB); each SPK segment is converted to AU before chain composition.
/// The 50 output slots and units are documented by
/// [`crate::orbital::orbital_elements_from_ecliptic_state`].
///
/// The registry snapshot is retained throughout the synchronous calculation.
/// Missing aliases, invalid states, coverage, configuration and model errors
/// propagate. This entry requires native JPL mode and does not acquire a new
/// kernel or substitute LEB data.
pub fn native_jpl_orbital_elements(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    epoch: f64,
    entry_is_ut1: bool,
    body: i32,
    flags: i32,
) -> Result<[f64; 50]> {
    let entry = if entry_is_ut1 {
        "get_orbital_elements_ut"
    } else {
        "get_orbital_elements"
    };
    crate::orbital::validate_orbital_object(entry, body, false)?;
    let original_time = orbital_source_time(epoch, entry_is_ut1)?;
    with_native_orbital_serving(ctx, model, iers, original_time, |orbital| {
        if entry_is_ut1 {
            orbital.orbital_elements_ut(body, epoch, flags)
        } else {
            orbital.orbital_elements_tt(body, original_time.value(), flags)
        }
    })
}

/// Maximum, minimum and current orbital-distance components at finite TT.
///
/// The extrema use osculating ellipse geometry over native ICRS AU/AU-day
/// states, with NAIF 3 for Earth's EMB. The current component is slot 2 of
/// geometric public output: radius for spherical flags, signed z for XYZ.
/// The shared solver is
/// [`crate::orbital::SealedOrbitalServing::orbit_max_min_true_distance_tt`].
///
/// The observation stays in TT. Configured time policies apply to any nested
/// model-defining epochs. The owned native registry is not replaced or
/// extended; missing support, coverage, configuration and flag errors propagate.
pub fn native_jpl_orbit_max_min_true_distance_tt(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    jd_tt: f64,
    body: i32,
    flags: i32,
) -> Result<(f64, f64, f64)> {
    crate::orbital::validate_orbital_object("orbit_max_min_true_distance", body, true)?;
    let original_time = orbital_source_time(jd_tt, false)?;
    with_native_orbital_serving(ctx, model, iers, original_time, |orbital| {
        orbital.orbit_max_min_true_distance_tt(body, original_time.value(), flags)
    })
}

/// Nodes and apsides over an owned native SPK context.
///
/// Returns ascending node, descending node, perihelion and aphelion (or second
/// focus with FOPOINT), each in the six-component public layout. Spherical
/// angles are degrees, radius is AU and rates are per TT day; representation,
/// equatorial and sidereal flags apply through the shared output reducer.
/// Method bits select mean/osculating geometry, barycentric elements and focus.
/// Undefined interpolated apsides retain NaN components.
///
/// `epoch` is TT or UT1 according to `entry_is_ut1`. The bundled observation
/// timescale supplies raw SPK TT/TDB; configured Delta-T/IERS settings apply
/// to nested lunar/minor-body models and sidereal defining epochs. Observer
/// parallax uses the snapshotted site and native IAU/Vondrák frames. Each read
/// retains its actual segment coverage. No new kernel or LEB source is acquired.
/// Missing aliases, invalid states, uncovered epochs and callback errors propagate;
/// explicit sealed LEB mode refuses this JPL-only entry.
#[allow(clippy::too_many_arguments)]
pub fn native_jpl_nod_aps(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    epoch: f64,
    entry_is_ut1: bool,
    body: i32,
    method: i32,
    flags: i32,
) -> Result<[PositionResult; 4]> {
    let original_time = orbital_source_time(epoch, entry_is_ut1)?;
    positions::forbid_jpl_in_sealed_mode(ctx.mode())?;
    if ctx.mode() != ProviderMode::NativeJpl {
        return Err(Error::ConfigurationMissing {
            missing_config: "native_nodaps_context",
            suggestion: "select NativeJpl and register the intended local SPK image",
        });
    }
    // Own mappings and bytes before borrowing the geometry/public feeds.
    let selected = ctx
        .clone()
        .with_delta_t_model(Some(model))
        .with_iers_delta_t(Some(iers));
    let provider = positions::LocalSpkStateProvider::new(selected.spk_registry(), selected.mode())?;
    let delta_cfg = selected.delta_t_config(model, iers);
    let delta = |ut| deltat_with_config(ut, &delta_cfg);
    let nutation = |tt| crate::erfa::nut06a(tt);
    let serving = SealedServing::new(&provider, &nutation);
    let session = crate::engine_api::Session::from_context(selected.clone());
    let calc_tt = |tt, body, flags| {
        session
            .calc(tt, body, flags)
            .map(|position| (position.values, position.returned_flags))
    };
    let calc_ut = |ut, body, flags| {
        session
            .calc_ut(ut, body, flags)
            .map(|position| (position.values, position.returned_flags))
    };
    let time_at = |tt| nodaps_source_time_at(original_time, tt);
    let icrs_at = |body, tt| provider.nodaps_state_two_part(body, time_at(tt)?);
    let emb_at = |tt| provider.nodaps_earth_barycenter_state_two_part(time_at(tt)?);
    let target_lookup = |body| provider.validate_nodaps_target_alias(body);
    let ut1_at = |tt| nodaps_source_ut1_at(original_time, epoch, entry_is_ut1, tt);
    let deflector_at = |body, tt| provider.nodaps_deflector_position(body, tt);
    let moon_state_at = |tt| provider.nodaps_true_moon_state(tt);
    let site_state_at = |tt, site| {
        // Each point's formatter creates a fresh scalar-TT Time, not the
        // original UT-created central Time. Earth alias priority precedes
        // topocenter then geocenter evaluation on that fresh clock.
        provider.validate_nodaps_target_alias(crate::constants::EARTH)?;
        nodaps_site_state_via(
            tt,
            site,
            |time| provider.nodaps_state_two_part(crate::constants::EARTH, time),
            |time, site| super::routing::local_jpl_topo(time.hi, time.lo, site),
        )
    };
    let (selection, _) = crate::houses::sidereal_selection_for_context(&selected);
    let anchors = crate::stars::SealedAnchorProvider::new(&provider);
    let aya = crate::sidereal::SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
    let nodaps = crate::orbital::SealedOrbitalServing::new(&serving, &delta)
        .with_nodaps_acquisition(
            &icrs_at,
            &emb_at,
            &target_lookup,
            &ut1_at,
            &calc_ut,
            &deflector_at,
        )
        .with_calc_tt(&calc_tt)
        .with_nodaps_site_state(&site_state_at)
        .with_nodaps_moon_state(&moon_state_at)
        .with_ayanamsha(aya);
    let tt = original_time.value();
    if flags & crate::constants::FLG_TOPOCTR != 0
        && flags & (crate::constants::FLG_HELCTR | crate::constants::FLG_BARYCTR) == 0
    {
        if let Some(observer) = selected.observer() {
            // The site callback owns its UT1 acquisition; this unused partner needs no conversion.
            return nodaps.nod_aps_topo_tt(
                body,
                tt,
                tt,
                method,
                flags,
                (observer.lon_deg(), observer.lat_deg(), observer.height_m()),
            );
        }
    }
    nodaps.nod_aps_tt(body, tt, method, flags)
}

/// Geocentric site position and velocity from native Earth support.
///
/// `tt` is a finite TT JD. Acquire Earth plus its surface offset and Earth
/// alone, then subtract in ICRS AU/AU-day. The site callback supplies the
/// bundled rotation clock and validates its geometry. Actual split-TDB
/// coverage and frame errors propagate without public range substitution.
pub(crate) fn nodaps_site_state_via(
    tt: f64,
    site: (f64, f64, f64),
    earth_at: impl Fn(TwoPartEpoch) -> Result<positions::GeometricState>,
    topo_at: impl Fn(TwoPartEpoch, (f64, f64, f64)) -> Result<([f64; 3], [f64; 3])>,
) -> Result<([f64; 3], [f64; 3])> {
    let time = orbital_source_time(tt, false)?;
    let earth_topo = earth_at(time)?;
    let (offset_pos, offset_vel) = topo_at(time, site)?;
    let topo_pos = std::array::from_fn::<_, 3, _>(|i| earth_topo.pos()[i] + offset_pos[i]);
    let topo_vel = std::array::from_fn::<_, 3, _>(|i| earth_topo.vel()[i] + offset_vel[i]);
    let geocenter = earth_at(time)?;
    Ok((
        std::array::from_fn(|i| topo_pos[i] - geocenter.pos()[i]),
        std::array::from_fn(|i| topo_vel[i] - geocenter.vel()[i]),
    ))
}

/// Reuse the central split observation epoch; other samples split their finite scalar TT.
fn nodaps_source_time_at(original: TwoPartEpoch, tt: f64) -> Result<TwoPartEpoch> {
    if tt.to_bits() == original.value().to_bits() {
        Ok(original)
    } else {
        orbital_source_time(tt, false)
    }
}

/// Read scalar epoch.ut1 for a native body/model/sidereal request. Original
/// UT1 construction stores an exact zero UT fraction, so central UT1 is the
/// caller word. All scalar TT constructors use the independent global inverse.
/// This channel is separate from observer formatting's fresh scalar epoch.
fn nodaps_source_ut1_at(original: TwoPartEpoch, epoch: f64, is_ut: bool, tt: f64) -> Result<f64> {
    if !tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tt });
    }
    if is_ut && tt.to_bits() == original.value().to_bits() {
        Ok(epoch)
    } else {
        crate::time::default_timescale_ut1_from_tt(tt)
    }
}

/// Construct the native orbital observation's finite split TT epoch.
/// TT inputs use integer day plus fraction. UT1 inputs retain the original
/// UT1 day and the refined bundled TT-minus-UT1 fraction. Configured
/// model Delta-T does not alter this observation clock.
fn orbital_source_time(epoch: f64, entry_is_ut1: bool) -> Result<TwoPartEpoch> {
    if !epoch.is_finite() {
        return Err(Error::NonFiniteArgument { value: epoch });
    }
    let (hi, lo) = if entry_is_ut1 {
        crate::smh2016::timescale_parts_from_ut1(epoch)?
    } else {
        let whole = epoch.floor();
        (whole, epoch - whole)
    };
    let original = TwoPartEpoch { hi, lo };
    for value in [hi, lo, original.value()] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    Ok(original)
}

/// Retain the original split observation when its scalar TT matches the
/// central epoch; split other scalar TT samples independently.
fn with_native_orbital_serving<R>(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    original_time: TwoPartEpoch,
    f: impl FnOnce(
        &crate::orbital::SealedOrbitalServing<'_, positions::LocalSpkStateProvider<'_>>,
    ) -> Result<R>,
) -> Result<R> {
    positions::forbid_jpl_in_sealed_mode(ctx.mode())?;
    if ctx.mode() != ProviderMode::NativeJpl {
        return Err(Error::ConfigurationMissing {
            missing_config: "native_orbital_context",
            suggestion: "select NativeJpl and register the intended local SPK image",
        });
    }
    // Clone retains owned image bytes and an immutable mapping snapshot before
    // any provider or nested public callback is borrowed. No I/O occurs here.
    let selected = ctx
        .clone()
        .with_delta_t_model(Some(model))
        .with_iers_delta_t(Some(iers));
    let provider = positions::LocalSpkStateProvider::new(selected.spk_registry(), selected.mode())?;
    let delta_cfg = selected.delta_t_config(model, iers);
    let delta = |ut| deltat_with_config(ut, &delta_cfg);
    let nutation = |tt| crate::erfa::nut06a(tt);
    let serving = SealedServing::new(&provider, &nutation);
    let session = crate::engine_api::Session::from_context(selected.clone());
    let calc_tt = |tt, body, flags| {
        session
            .calc(tt, body, flags)
            .map(|position| (position.values, position.returned_flags))
    };
    let time_at = |tt: f64| {
        if tt.to_bits() == original_time.value().to_bits() {
            Ok(original_time)
        } else {
            orbital_source_time(tt, false)
        }
    };
    let icrs_at = |body, tt| provider.geometric_state_two_part(body, time_at(tt)?);
    let emb_at = |tt| provider.earth_barycenter_state_two_part(time_at(tt)?);
    let target_lookup = |body| provider.validate_orbital_target_alias(body);
    let orbital = crate::orbital::SealedOrbitalServing::new(&serving, &delta)
        .with_orbital_geometry(&icrs_at, &emb_at, &target_lookup)
        .with_calc_tt(&calc_tt)
        .with_current_calc_tt(&calc_tt)
        .with_orbital_ut1_to_tt(&crate::time::default_timescale_tt_from_ut1);
    f(&orbital)
}

/// Physical phenomena from an owned native JPL/SPK context.
///
/// The finite epoch is TT or original UT1 according to `entry_is_ut1`.
/// The bundled timescale supplies the observation's TT/TDB and UT1; configured
/// Delta-T/IERS settings apply to nested lunar/model/asteroid calculations.
/// Surface observers use the IERS 2010 ellipsoid and native celestial frames.
///
/// Output contains phase angle (degrees), illuminated fraction, elongation
/// (degrees), diameter (degrees), magnitude and lunar parallax (degrees), then
/// reserved zeros. Physical-body geometry uses ICRS vectors, independent of
/// representation and sidereal flags. TRUEPOS selects geometric target/Sun
/// geometry except that solar observation remains apparent. Heliocentric target
/// distance uses the retarded state. Analytical families enter configured
/// [`crate::native::Session::calc_ut`] callbacks.
///
/// The owned registry snapshot serves every chain, with per-segment AU
/// conversion and strict support aliases. Missing registration, coverage,
/// configuration and corrupt-state errors propagate. No alternate kernel or
/// LEB image is substituted; sealed LEB mode refuses this native entry.
#[allow(clippy::too_many_arguments)]
pub fn native_jpl_pheno(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    epoch: f64,
    entry_is_ut1: bool,
    body: i32,
    flags: i32,
) -> Result<[f64; 20]> {
    native_jpl_pheno_with_calc_ut(ctx, model, iers, epoch, entry_is_ut1, body, flags, None)
}

/// Native phenomena with an optional call-borrowed public position callback.
///
/// Only analytical/asteroid branches consume the six-component degrees/AU
/// result, and only its position components enter photometry. The callback,
/// registry snapshot and errors are confined to this synchronous call.
/// Without a callback, acquisition follows [`native_jpl_pheno`].
#[allow(clippy::too_many_arguments)]
pub fn native_jpl_pheno_with_calc_ut(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    epoch: f64,
    entry_is_ut1: bool,
    body: i32,
    flags: i32,
    acquire: Option<&crate::engine_api::PhenoPosition<'_>>,
) -> Result<[f64; 20]> {
    let callbacks = acquire.map(crate::engine_api::PhenoCallbacks::new);
    native_jpl_pheno_with_callbacks(
        ctx,
        model,
        iers,
        epoch,
        entry_is_ut1,
        body,
        flags,
        callbacks.as_ref(),
    )
}

/// Native phenomena with borrowed positions and optional raw scalar math.
/// Only the completed positional callback branch consumes raw arithmetic;
/// physical vectors and the ordinary None route retain their exact kernels.
/// Original UT1, flags, provider/observer admission and callback errors remain
/// owned by this synchronous call. Native never invokes the direct-LEB retry.
#[allow(clippy::too_many_arguments)]
pub fn native_jpl_pheno_with_callbacks(
    ctx: &EphemerisContext,
    model: DeltaTModel,
    iers: bool,
    epoch: f64,
    entry_is_ut1: bool,
    body: i32,
    flags: i32,
    callbacks: Option<&crate::engine_api::PhenoCallbacks<'_>>,
) -> Result<[f64; 20]> {
    let acquire = callbacks.map(crate::engine_api::PhenoCallbacks::position);
    let source_time = callbacks
        .and_then(crate::engine_api::PhenoCallbacks::observation_ut1)
        .is_some();
    positions::forbid_jpl_in_sealed_mode(ctx.mode())?;
    if ctx.mode() != ProviderMode::NativeJpl {
        return Err(Error::ConfigurationMissing {
            missing_config: "native_pheno_context",
            suggestion: "select NativeJpl and register the intended local SPK image",
        });
    }
    if !epoch.is_finite() {
        return Err(Error::NonFiniteArgument { value: epoch });
    }
    if let Some(values) = crate::photometry::degenerate_pheno(body) {
        return Ok(values);
    }
    let observation = NativePhenoObservation::from_entry(epoch, entry_is_ut1)?;
    let time = observation.split_tt;
    let selected = ctx
        .clone()
        .with_delta_t_model(Some(model))
        .with_iers_delta_t(Some(iers));
    let provider = positions::LocalSpkStateProvider::new(selected.spk_registry(), selected.mode())?;
    let body = positions::remap_ast_offset(body);
    // Alias resolution owns source priority; no state is sampled here.
    provider.validate_source_alias(crate::constants::EARTH)?;
    provider.validate_source_alias(crate::constants::SUN)?;
    let site = if flags & crate::constants::FLG_TOPOCTR != 0 {
        Some(selected.observer().ok_or(Error::ConfigurationMissing {
            missing_config: "observer_location",
            suggestion: "FLG_TOPOCTR requires a configured observer site",
        })?)
    } else {
        None
    };
    let session = crate::engine_api::Session::from_context(selected.clone());
    let calc_ut = |ut1, target, calc_flags| {
        let values = match acquire {
            Some(calc_ut) => calc_ut(ut1, target, calc_flags)?,
            None => session.calc_ut(ut1, target, calc_flags)?.values,
        };
        crate::photometry::PhenoState::new(values[0], values[1], values[2])
    };
    if crate::photometry::asteroid_photometry(body).is_some() {
        let ut1 = match callbacks.and_then(crate::engine_api::PhenoCallbacks::observation_ut1) {
            Some(observation) => observation()?,
            None => observation.ut1()?,
        };
        let base = flags
            & (crate::constants::FLG_TRUEPOS
                | crate::constants::FLG_NONUT
                | crate::constants::FLG_TOPOCTR);
        let helio_flags = (base & !crate::constants::FLG_TOPOCTR) | crate::constants::FLG_HELCTR;
        let (target, sun, helio) = if let Some(acquire) = acquire {
            // Acquire all required public states before reducing their physical geometry.
            let (target, sun, helio) =
                native_pheno_public_states(ut1, body, base, Some(helio_flags), acquire)?;
            (
                target,
                sun,
                helio.expect("asteroid stage supplies heliocentric flags"),
            )
        } else {
            (
                calc_ut(ut1, body, base)?,
                calc_ut(ut1, crate::constants::SUN, base)?,
                calc_ut(ut1, body, helio_flags)?,
            )
        };
        return crate::photometry::asteroid_pheno_from_public(body, target, sun, helio);
    }
    if body != crate::constants::SUN
        && body != crate::constants::MOON
        && !crate::photometry::is_disc_body(body)
    {
        let source_body = crate::lunar::south_node_north_source(body).unwrap_or(body);
        if !source_time
            && positions::classify_body(source_body) == positions::BodyKind::Unknown
            && selected.spk_registry().body_info(body).is_none()
        {
            return Err(Error::UnknownBody { body_id: body });
        }
        let ut1 = match callbacks.and_then(crate::engine_api::PhenoCallbacks::observation_ut1) {
            Some(observation) => observation()?,
            None => observation.ut1()?,
        };
        let base = flags & (crate::constants::FLG_TRUEPOS | crate::constants::FLG_TOPOCTR);
        if let Some(callbacks) = callbacks {
            let target = (callbacks.position())(ut1, body, base)?;
            let sun = (callbacks.position())(ut1, crate::constants::SUN, base)?;
            return crate::photometry::positional_pheno_from_raw(body, target, sun);
        }
        let (target, sun) = (
            calc_ut(ut1, body, base)?,
            calc_ut(ut1, crate::constants::SUN, base)?,
        );
        return crate::photometry::positional_pheno(body, target, sun);
    }
    provider.validate_source_alias(body)?;
    let geometric = flags & crate::constants::FLG_TRUEPOS != 0 && body != crate::constants::SUN;
    let earth = provider
        .pheno_state_two_part(crate::constants::EARTH, time, source_time)
        .map_err(|error| native_pheno_observation_error(time, geometric, error))?;
    let (mut observer, mut velocity) = (earth.pos(), earth.vel());
    let surface = if let Some(site) = site {
        let ut1_fraction = observation.ut1_fraction()?;
        let rotation = positions::reference_frames::itrs_rotation_at_parts(
            time.hi,
            time.lo,
            ut1_fraction,
            time.value(),
        )?;
        let itrs = native_pheno_itrs_position(
            (site.lon_deg(), site.lat_deg(), site.height_m()),
            positions::SiteEllipsoid::Iers2010,
        )?;
        let (position, speed) = positions::reference_frames::itrs_to_gcrs(rotation, itrs)?;
        for i in 0..3 {
            observer[i] += position[i];
            velocity[i] += speed[i];
        }
        Some(position)
    } else {
        None
    };
    let target = if geometric {
        provider.pheno_geometric_difference(
            body,
            crate::constants::EARTH,
            time,
            surface,
            source_time,
        )?
    } else {
        provider.pheno_observed_vector(
            body,
            time,
            observer,
            velocity,
            surface,
            false,
            true,
            source_time,
        )?
    };
    if body == crate::constants::SUN {
        return crate::photometry::pheno_from_native_vectors(
            body,
            target,
            [0.0; 3],
            [0.0; 3],
            time.value(),
            (0.0, 0.0),
            (0.0, 0.0),
            0.0,
        );
    }
    let sun = if geometric {
        provider.pheno_geometric_difference(
            crate::constants::SUN,
            crate::constants::EARTH,
            time,
            surface,
            source_time,
        )?
    } else {
        provider.pheno_observed_vector(
            crate::constants::SUN,
            time,
            observer,
            velocity,
            surface,
            false,
            true,
            source_time,
        )?
    };
    let sun_state = provider.pheno_state_two_part(crate::constants::SUN, time, source_time)?;
    let helio = provider.pheno_observed_vector(
        body,
        time,
        sun_state.pos(),
        sun_state.vel(),
        None,
        false,
        false,
        source_time,
    )?;
    let parallax = if body == crate::constants::MOON {
        if site.is_some() {
            let geo = provider.pheno_observed_vector(
                body,
                time,
                earth.pos(),
                earth.vel(),
                None,
                false,
                true,
                source_time,
            )?;
            let topo = provider.pheno_observed_vector(
                body,
                time,
                observer,
                velocity,
                surface,
                false,
                true,
                source_time,
            )?;
            native_pheno_parallax(geo, topo)
        } else {
            let geom = provider.pheno_geometric_difference(
                body,
                crate::constants::EARTH,
                time,
                None,
                source_time,
            )?;
            crate::photometry::moon_horizontal_parallax_deg(native_pheno_vector_length(geom))
        }
    } else {
        0.0
    };
    // Saturn/Uranus ring geometry requires a valid ecliptic frame and direction.
    let (geo_ecliptic, helio_ecliptic) =
        if matches!(body, crate::constants::SATURN | crate::constants::URANUS) {
            let frame = native_pheno_ecliptic_frame(time);
            let coords = |vector| match &frame {
                Ok(frame) => native_pheno_frame_angles(*frame, vector),
                Err(error) => Err(error.clone()),
            };
            (coords(target)?, coords(helio)?)
        } else {
            ((0.0, 0.0), (0.0, 0.0))
        };
    crate::photometry::pheno_from_native_vectors(
        body,
        target,
        sun,
        helio,
        time.value(),
        geo_ecliptic,
        helio_ecliptic,
        parallax,
    )
}

/// Preserve the observation algorithm's error when a finite Time cannot be
/// represented as ephemeris seconds. The split TT-to-TDB series is evaluated
/// at the original epoch; each TDB day part is scaled by 86400 seconds. If a
/// rounded product or their sum overflows, an apparent light-time solution
/// cannot converge. This is a convergence error with the ten-step,
/// 1e-12-day iteration contract, rather than a non-finite caller argument.
/// Geometric requests and all other provider failures retain their identity;
/// the normal provider is invoked first, retaining its acquisition and guard
/// precedence. A seconds overflow can prevent later directory or state reads;
/// only that derived-seconds failure is translated. No state is synthesized.
fn native_pheno_observation_error(time: TwoPartEpoch, geometric: bool, error: Error) -> Error {
    if !geometric
        && matches!(error, Error::NonFiniteArgument { .. })
        && time.hi.is_finite()
        && time.lo.is_finite()
    {
        if let Ok(delta) = positions::reference_frames::tdb_minus_tt_seconds(time.hi, time.lo) {
            let fraction = time.lo + delta / positions::reference_frames::FRAME_DAY_S;
            let a = (time.hi - crate::spk::J2000_JD_TDB) * positions::reference_frames::FRAME_DAY_S;
            let b = fraction * positions::reference_frames::FRAME_DAY_S;
            if !a.is_finite() || !b.is_finite() || !(a + b).is_finite() {
                return Error::ConvergenceFailed {
                    algorithm: "spk_light_time",
                    iterations: 10,
                    max_iterations: 10,
                    tolerance: 1e-12,
                    last_value: f64::NAN,
                };
            }
        }
    }
    error
}

/// Complete the public target/Sun and optional heliocentric target acquisitions
/// before checking any output state. All six raw components are retained until
/// the calls complete; later acquisition errors precede derived geometry errors.
/// The returned third state is Some exactly when heliocentric flags are Some.
/// Checked nonfinite/distance geometry domains remain unchanged after staging.
fn native_pheno_public_states(
    ut1: f64,
    body: i32,
    base_flags: i32,
    heliocentric_flags: Option<i32>,
    acquire: &crate::engine_api::PhenoPosition<'_>,
) -> Result<(
    crate::photometry::PhenoState,
    crate::photometry::PhenoState,
    Option<crate::photometry::PhenoState>,
)> {
    let target = acquire(ut1, body, base_flags)?;
    let sun = acquire(ut1, crate::constants::SUN, base_flags)?;
    let helio = heliocentric_flags
        .map(|flags| acquire(ut1, body, flags))
        .transpose()?;
    Ok((
        crate::photometry::PhenoState::new(target[0], target[1], target[2])?,
        crate::photometry::PhenoState::new(sun[0], sun[1], sun[2])?,
        helio
            .map(|p| crate::photometry::PhenoState::new(p[0], p[1], p[2]))
            .transpose()?,
    ))
}

/// Matching observation clocks for one native phenomenon request.
///
/// `split_tt` retains the exact original TT words, which need not be a
/// canonical day/fraction split. `original_ut1_fraction`, when present, is a
/// fraction relative to that split's leading word, never a full UT1 JD.
/// Keeping them together prevents mixing this fraction with other entries'
/// optional scalar UT1 dates. Inverse clock work remains lazy.
#[derive(Clone, Copy)]
struct NativePhenoObservation {
    split_tt: TwoPartEpoch,
    original_ut1_fraction: Option<f64>,
}

impl NativePhenoObservation {
    /// Preserve the entry's existing finite check and physical clock choice.
    /// TT retains floor/remainder only; UT1 acquires the independent refined
    /// global clock and its known fraction. No inverse Delta-T domain is added.
    fn from_entry(epoch: f64, entry_is_ut1: bool) -> Result<Self> {
        if !epoch.is_finite() {
            return Err(Error::NonFiniteArgument { value: epoch });
        }
        if entry_is_ut1 {
            let (epochs, _) = positions::reference_frames::VectorRotationEpochs::from_ut1(epoch)?;
            Ok(Self {
                split_tt: TwoPartEpoch {
                    hi: epochs.whole,
                    lo: epochs.tt_fraction,
                },
                original_ut1_fraction: Some(epochs.ut1_fraction),
            })
        } else {
            let hi = epoch.floor();
            Ok(Self {
                split_tt: TwoPartEpoch { hi, lo: epoch - hi },
                original_ut1_fraction: None,
            })
        }
    }

    /// Acquire the UT1 fraction only for a consuming site or public callback.
    /// An original UT1 entry reuses its known fraction. TT observations consult
    /// the bundled timescale without replacing or canonicalizing either word.
    fn ut1_fraction(self) -> Result<f64> {
        match self.original_ut1_fraction {
            Some(fraction) => Ok(fraction),
            None => Ok(self.split_tt.lo
                - crate::smh2016::bundled_deltat_seconds(self.split_tt.value())?
                    / positions::reference_frames::FRAME_DAY_S),
        }
    }

    /// Scalar UT1 for a genuine public callback, not geocentric TT composition.
    fn ut1(self) -> Result<f64> {
        Ok(self.split_tt.hi + self.ut1_fraction()?)
    }
}

/// Geodetic site ITRS vector for the chosen reference ellipsoid.
/// East longitude/latitude are degrees, finite elevation is metres, output
/// is AU. IERS2010 and WGS84 each retain their independently specified
/// model. Latitude bounds belong to Observer, while every consumed argument
/// is finite-checked here. With `k=(f_inv-1)/f_inv`, curvature factors are
/// `c=1/sqrt(cos²(phi)+k²*sin²(phi))` and `s=k²*c`; height is added to each
/// radius component before the latitude projection. Finite negative heights,
/// including a zero equatorial site vector, use the same ellipsoid formula.
pub(super) fn native_pheno_itrs_position(
    site: (f64, f64, f64),
    ellipsoid: positions::SiteEllipsoid,
) -> Result<[f64; 3]> {
    use positions::reference_frames::{FRAME_AU_M, FRAME_TAU, scalar_cos, scalar_sin};
    let (lon, lat, height) = site;
    for value in [lon, lat, height] {
        if !value.is_finite() {
            return Err(Error::NonFiniteArgument { value });
        }
    }
    let phi = lat / 360.0 * FRAME_TAU;
    let lambda = lon / 360.0 * FRAME_TAU;
    let radius = ellipsoid.semi_major_m() / FRAME_AU_M;
    let elevation = height / FRAME_AU_M;
    let sinphi = scalar_sin(phi);
    let cosphi = scalar_cos(phi);
    let inverse = ellipsoid.inv_flattening();
    let omf = (inverse - 1.0) / inverse;
    let omf2 = omf * omf;
    let c = 1.0 / (cosphi * cosphi + sinphi * sinphi * omf2).sqrt();
    let s = omf2 * c;
    let radius_xy = radius * c;
    let xy = (radius_xy + elevation) * cosphi;
    let radius_z = radius * s;
    Ok([
        xy * scalar_cos(lambda),
        xy * scalar_sin(lambda),
        (radius_z + elevation) * sinphi,
    ])
}

/// Transform a finite direction and return longitude/latitude in degrees.
/// Non-finite matrix/vector inputs and arithmetic overflow are typed failures.
pub(super) fn native_pheno_frame_angles(
    frame: [[f64; 3]; 3],
    vector: [f64; 3],
) -> Result<(f64, f64)> {
    use positions::reference_frames::{FRAME_TAU, floored_mod, matvec3};
    let [x, y, z] = matvec3(frame, vector)?;
    let latitude = z.atan2(x.hypot(y));
    let longitude = floored_mod(y.atan2(x), FRAME_TAU);
    Ok((longitude * 360.0 / FRAME_TAU, latitude * 360.0 / FRAME_TAU))
}

/// Left-associated three-vector length (AU), shared by native parallax reads.
fn native_pheno_vector_length(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Angular geocenter-to-topocenter lunar displacement in degrees.
/// Vectors are observer-relative ICRS AU at the same epoch. Positive norms
/// use the clamped direction-cosine acos; a zero vector returns zero.
fn native_pheno_parallax(geo: [f64; 3], topo: [f64; 3]) -> f64 {
    let g = native_pheno_vector_length(geo);
    let t = native_pheno_vector_length(topo);
    if g > 0.0 && t > 0.0 {
        crate::photometry::native_pheno_cosine_clamp(
            (geo[0] * topo[0] + geo[1] * topo[1] + geo[2] * topo[2]) / (g * t),
        )
        .acos()
        .to_degrees()
    } else {
        0.0
    }
}

/// True-ecliptic-of-date rotation from native IAU 2000A/P03 frames.
/// Nutation uses TT; precession and mean obliquity use TDB obtained from the
/// retained observation split. Non-finite arguments and frame errors propagate.
pub(super) fn native_pheno_ecliptic_frame(time: TwoPartEpoch) -> Result<[[f64; 3]; 3]> {
    use positions::reference_frames::*;
    let tt = time.value();
    let tdb = time.hi + (time.lo + tdb_minus_tt_seconds(time.hi, time.lo)? / FRAME_DAY_S);
    let (dpsi, deps) = iau2000a_tenth_microarcsec(tt)?;
    let dpsi = dpsi * FRAME_TENTH_USEC_2_RAD;
    let deps = deps * FRAME_TENTH_USEC_2_RAD;
    let eps = mean_obliquity_arcsec(tdb)? * FRAME_ASEC2RAD;
    let nut = nutation_matrix(eps, eps + deps, dpsi)?;
    let equinox = equinox_matrix(nut, precession_matrix(tdb)?)?;
    ecliptic_date_matrix(equinox, eps, deps)
}

/// Borrow a fixed-epoch engine for one complete UT or TT entry over the
/// caller's tropical core. Shared sidereal kernels own epoch selection,
/// rotation, nutation and speed stencils; the callback cannot retain providers.
pub fn with_fixed_epoch_engine<R>(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tropical: &dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    with_fixed_epoch_engine_with_anchors(stack, ctx, tropical, None, f)
}

/// The fixed-frame engine with an optional live dynamic anchor provider.
/// Only consumed ayanamsha callbacks read this provider; tropical geometry,
/// stored frame samples and configured defining-epoch clocks remain unchanged.
pub fn with_fixed_epoch_engine_with_anchors<R>(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tropical: &dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>,
    anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    with_fixed_epoch_engine_with_anchor_clocks(stack, ctx, tropical, anchors, anchors, f)
}

/// Fixed-frame serving with distinct fast-base and plain-UT projection owners.
/// The plain projection helper owns independent global Time; its base helper
/// may own fast scalar-TT delegation. None retains all original typed clocks.
#[allow(clippy::too_many_arguments)]
pub fn with_fixed_epoch_engine_with_anchor_clocks<R>(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tropical: &dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>,
    anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    projection_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    with_fixed_epoch_engine_impl(stack, ctx, tropical, anchors, projection_anchors, false, f)
}

/// Private native binding: the supplied closure is the same stored calculation.
/// Public opaque closures never acquire the local word-motion capability.
pub(crate) fn with_fixed_epoch_engine_with_native_motion<R>(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tropical: &dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>,
    anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    projection_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    with_fixed_epoch_engine_impl(stack, ctx, tropical, anchors, projection_anchors, true, f)
}

#[allow(clippy::too_many_arguments)]
fn with_fixed_epoch_engine_impl<R>(
    stack: &ServingStack<'_, '_>,
    ctx: &EphemerisContext,
    tropical: &dyn Fn(f64, i32, i32) -> Result<(PositionResult, i32)>,
    anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    projection_anchors: Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    native_motion: bool,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    let to_tt = |jd_ut: f64| Ok(jd_ut + deltat_with_config(jd_ut, stack.delta_cfg)?);
    let (selection, _) = crate::houses::sidereal_selection_for_context(ctx);
    let deltat_days = |jd: f64| deltat_with_config(jd, stack.delta_cfg);
    let ayanamsha = crate::sidereal::SealedAyanamshaServing::new(selection, &deltat_days, anchors);
    let mean_tt = |jd_tt: f64, noaberr: bool, nogdefl: bool| {
        ayanamsha
            .mean_tt_with(jd_tt, noaberr, nogdefl)
            .map(|(value, _)| value)
    };
    let projection =
        crate::sidereal::SealedAyanamshaServing::new(selection, &deltat_days, projection_anchors);
    let projection_tt = std::cell::Cell::new(None);
    let mean_ut_plain = |jd_ut: f64| {
        if projection_anchors.is_some() {
            let tt = crate::time::default_timescale_tt_from_ut1(jd_ut)?;
            projection_tt.set(Some(tt));
            projection.mean_tt(tt).map(|(value, _)| value)
        } else {
            ayanamsha.mean_ut(jd_ut).map(|(value, _)| value)
        }
    };
    let dpsi_deg = |jd_tt: f64| (stack.nutation_at)(jd_tt).map(|(dpsi, _)| dpsi.to_degrees());
    let dpsi_rate = |jd_tt: f64| stack.nutation_longitude_rate_deg_day(jd_tt);
    let (_, t0, _) = selection.full();
    let user_t0_tt = if selection.bits() & crate::constants::SIDBIT_USER_UT != 0 {
        t0 + deltat_days(t0)?
    } else {
        t0
    };
    let valens_epoch_tt = ayanamsha.valens_epoch_tt()?;
    let projection_words = |route, tt, equatorial| {
        use crate::positions::sidereal_frame_words::{self as words, EpochPairZero};
        use crate::sidereal::SiderealRoute;
        let mode = crate::sidereal::strip_to_base_mode(selection.mode());
        match route {
            SiderealRoute::FixedEpoch if mode == crate::constants::SIDM_J2000 => Ok(None),
            SiderealRoute::FixedEpoch => words::fixed_epoch(mode, equatorial).map(Some),
            SiderealRoute::SidbitEquatorialT0 => words::equatorial_epoch(tt).map(Some),
            SiderealRoute::SidbitEcliptic => {
                // Only the actual bound native zero clock may supply words.
                // Opaque/no-anchor projection callbacks keep their old route.
                let Some(zero_tt) = projection_tt.get() else {
                    return Ok(None);
                };
                let pair = if mode == crate::constants::SIDM_USER {
                    let (_, _, angle) = selection.full();
                    Some((angle, user_t0_tt))
                } else if mode == crate::constants::SIDM_VALENS_MOON {
                    Some((crate::sidereal::VALENS_MOON_AYAN_T0_DEG, valens_epoch_tt))
                } else {
                    crate::sidereal::AYANAMSHA_DEFINING
                        .iter()
                        .find(|(id, _, _)| *id == mode)
                        .map(|(_, angle, epoch)| (*angle, *epoch))
                };
                let Some((angle, defining_tt)) = pair else {
                    return Ok(None);
                };
                words::ecliptic_projection(
                    selection.bits(),
                    tt,
                    EpochPairZero {
                        tt: zero_tt,
                        defining_tt,
                        defining_degrees: angle,
                        ecliptic_date: selection.bits() & crate::constants::SIDBIT_ECL_DATE != 0,
                    },
                )
                .map(Some)
            }
            _ => Ok(None),
        }
    };
    let sidereal = crate::sidereal::SiderealFastServing::new(
        selection,
        &tropical,
        &mean_tt,
        &mean_ut_plain,
        &dpsi_deg,
        &dpsi_rate,
        user_t0_tt,
        valens_epoch_tt,
    );
    let epoch_pair_rate = &super::corrections::epoch_pair_ayanamsha_rate_deg_day;
    let sidereal = if stack.serving.provider.has_analytic_kinematics() {
        sidereal
            .with_epoch_pair_rate(epoch_pair_rate)
            .with_projection_words(&projection_words)
    } else {
        sidereal
    };
    let tropical_words = |tt, body, flags| {
        stack
            .serving
            .with_output_motion(tt, body, flags, || tropical(tt, body, flags))
            .map(|(result, motion)| crate::sidereal::TropicalWordResult { result, motion })
    };
    let sidereal = if native_motion && stack.serving.provider.has_analytic_kinematics() {
        sidereal.with_tropical_words(&tropical_words)
    } else {
        sidereal
    };
    let live_rate = |tt, noaberr, nogdefl| ayanamsha.live_mean_rate_deg_day(tt, noaberr, nogdefl);
    let sidereal = sidereal.with_live_anchor_rate(&live_rate);
    let sid_compute = |tt: f64, id: i32, calc_flags: i32| sidereal.compute_tt(tt, id, calc_flags);
    let engine = CalcEngine {
        ut1_to_tt_jd: &to_tt,
        compute_tt: &sid_compute,
        nutation_obliquity_ut: &positions::ecl_nut_result_ut,
        nutation_obliquity_tt: &|jd_tt: f64, nut_flags: i32| {
            positions::ecl_nut_result(jd_tt, nut_flags, &crate::erfa::nut06a)
        },
        jpl_fallback_active: false,
    };
    f(&engine)
}

/// Assemble a borrowed calculation engine over the selected LEB tier set.
///
/// Multiple recognized LEB tiers use date-aware selection; other requests
/// delegate to the SDK factory. Explicit sites precede the context observer.
/// Stored channels, nutation, WGS84 rotation, output stencils and typed errors
/// retain their provider contracts. This metadata-free factory treats site
/// samples as TT entries, using configured Delta-T's two-evaluation inverse
/// `U = T - D(T - D(T))` and the selected reader's table at the original `T`
/// only for the classified calendar-NaN fallback. Times and Delta-T use days.
/// Use [`with_shared_calc_engine_at_entry`] to retain an original UT1 scalar
/// for one topocentric UT entry. NativeJpl ignores LEB images and retains its
/// separately bundled clocks; there is no provider or network substitution.
pub fn with_shared_calc_engine<R>(
    request: ServingRequest<'_>,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    with_shared_calc_engine_at_entry(request, None, f)
}

/// Bind date-aware tier callbacks with immutable original UT1 entry metadata.
///
/// `original_ut` is Some(U) for one UT1 entry and None for TT. LEB site clocks
/// use configured Delta-T in days: TT gives `U = T - D(T - D(T))`; a UT central
/// sample retains U exactly and other TT samples freeze `(U + D(U)) - U`.
/// The topocentric reducer reuses its central difference for its internal
/// SPEED samples. Clock lookup is lazy after required observer resolution in
/// the TT callback. A UT calculation caller still owns observer preflight
/// before forward time conversion; the factory has no body/flags at entry.
/// Only a TT inverse's NaN calendar argument under Espenak-Meeus, disabled
/// IERS and no override uses this reader's Delta-T table at original TT.
/// Table failures and all unclassified clock failures retain typed errors.
///
/// Provider priority, coverage, six spherical coordinates/rates, retflags and
/// synchronous borrowed callback lifetime follow [`with_shared_calc_engine`].
/// NativeJpl and single-tier requests use the entry-aware SDK factory without
/// changing native scalar/split, frame, lunar or model acquisition clocks.
/// Entry metadata is neither mutated nor retained after the callback returns.
pub fn with_shared_calc_engine_at_entry<R>(
    request: ServingRequest<'_>,
    original_ut: Option<f64>,
    f: impl FnOnce(&CalcEngine<'_>) -> Result<R>,
) -> Result<R> {
    if !matches!(
        request.ctx.mode(),
        ProviderMode::SealedLeb | ProviderMode::Auto
    ) {
        return positions::with_shared_calc_engine_at_entry(request, original_ut, f);
    }
    let distinct = LebTier::PRIORITY
        .into_iter()
        .filter(|tier| {
            request.ctx.leb_files().iter().any(|file| {
                file_tier(
                    file.path()
                        .rsplit(['/', '\\'])
                        .next()
                        .unwrap_or(file.path()),
                ) == Some(*tier)
            })
        })
        .count();
    if distinct <= 1 {
        return positions::with_shared_calc_engine_at_entry(request, original_ut, f);
    }
    with_sealed_serving_stack(
        request.ctx,
        request.delta_t_model,
        request.iers_enabled,
        request.skip_closed,
        |stack| {
            let observer = request.observer.or_else(|| request.ctx.observer());
            let site = observer
                .map(|observer| (observer.lon_deg(), observer.lat_deg(), observer.height_m()));
            let to_tt = |jd| Ok(jd + deltat_with_config(jd, stack.delta_cfg)?);
            let rotation = |tt, ut| crate::erfa::c2t06a(tt, ut, 0.0, 0.0);
            let compute = |tt: f64, id: i32, flags: i32| {
                if flags & crate::constants::FLG_TOPOCTR != 0 {
                    match positions::resolve_topo_request(flags, id, site, None)? {
                        TopoPlan::Geocentric | TopoPlan::TopoNoop => {
                            positions::serving_compute_tt_escalated(
                                stack.mode,
                                stack.serving,
                                tt,
                                id,
                                flags,
                            )
                        }
                        TopoPlan::Topocentric {
                            lon_deg,
                            lat_deg,
                            height_m,
                        } => {
                            let ut = super::site_clock::configured_site_ut1_partner_with_table(
                                tt,
                                original_ut,
                                stack.delta_cfg,
                                &|jd| stack.serving.provider().delta_t(jd),
                            )?;
                            let ut1_at = |sample_tt| {
                                super::site_clock::configured_site_sample_ut1(
                                    sample_tt,
                                    original_ut,
                                    stack.delta_cfg,
                                    &|jd| stack.serving.provider().delta_t(jd),
                                )
                            };
                            super::serving::serving_compute_topocentric_with_clock(
                                stack.serving,
                                tt,
                                ut,
                                id,
                                flags,
                                (lon_deg, lat_deg, height_m),
                                &rotation,
                                false,
                                &ut1_at,
                                if flags & crate::constants::FLG_SPEED != 0 {
                                    super::site_clock::configured_site_clock_support(
                                        tt,
                                        original_ut,
                                        stack.delta_cfg,
                                    )?
                                } else {
                                    (f64::NEG_INFINITY, f64::INFINITY)
                                },
                            )
                        }
                    }
                } else {
                    positions::serving_compute_tt_escalated(
                        stack.mode,
                        stack.serving,
                        tt,
                        id,
                        flags,
                    )
                }
            };
            let nut_tt = |tt, flags| positions::ecl_nut_result(tt, flags, &crate::erfa::nut06a);
            let engine = CalcEngine {
                ut1_to_tt_jd: &to_tt,
                compute_tt: &compute,
                nutation_obliquity_ut: &positions::ecl_nut_result_ut,
                nutation_obliquity_tt: &nut_tt,
                jpl_fallback_active: false,
            };
            f(&engine)
        },
    )
}

#[cfg(test)]
mod native_orbital_context_tests {
    use super::*;

    #[test]
    fn native_pheno_observation_error_separates_derived_overflow_from_raw_failures() {
        let failed = || Error::NonFiniteArgument {
            value: f64::INFINITY,
        };
        let extreme = TwoPartEpoch {
            hi: f64::MAX,
            lo: 0.0,
        };
        assert!(matches!(
            native_pheno_observation_error(extreme, false, failed()),
            Error::ConvergenceFailed { iterations: 10, .. }
        ));
        assert!(matches!(
            native_pheno_observation_error(extreme, true, failed()),
            Error::NonFiniteArgument { .. }
        ));
        assert!(matches!(
            native_pheno_observation_error(
                TwoPartEpoch {
                    hi: 2451545.0,
                    lo: 0.125
                },
                false,
                failed()
            ),
            Error::NonFiniteArgument { .. }
        ));
        assert!(matches!(
            native_pheno_observation_error(extreme, false, Error::SpkBodyNotFound { target: 399 }),
            Error::SpkBodyNotFound { target: 399 }
        ));
    }

    #[test]
    fn original_ut_time_retains_the_global_pair() {
        for ut in [2451545.125, 2460400.25] {
            let original = orbital_source_time(ut, true).unwrap();
            let (whole, fraction) = crate::smh2016::timescale_parts_from_ut1(ut).unwrap();
            assert_eq!(original.hi.to_bits(), ut.to_bits());
            assert_eq!(original.hi.to_bits(), whole.to_bits());
            assert_eq!(original.lo.to_bits(), fraction.to_bits());
            assert_eq!(
                original.value().to_bits(),
                crate::time::default_timescale_tt_from_ut1(ut)
                    .unwrap()
                    .to_bits()
            );
            // The scalar TT sum alone does not retain the original split observation parts.
            let scalar = orbital_source_time(original.value(), false).unwrap();
            assert_ne!(original.lo.to_bits(), scalar.lo.to_bits());
        }
    }

    #[test]
    fn tt_time_preserves_scalar_construction_and_rejects_nonfinite() {
        for tt in [-1.125_f64, 0.0, 2451545.125, 2460400.25] {
            let original = orbital_source_time(tt, false).unwrap();
            assert_eq!(original.hi.to_bits(), tt.floor().to_bits());
            assert_eq!(original.lo.to_bits(), (tt - tt.floor()).to_bits());
            assert_eq!(original.value().to_bits(), tt.to_bits());
        }
        for epoch in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            for is_ut in [false, true] {
                assert!(matches!(
                    orbital_source_time(epoch, is_ut),
                    Err(Error::NonFiniteArgument { .. })
                ));
            }
        }
    }

    #[test]
    fn signed_zero_tt_reduces_at_the_time_observable() {
        let input = -0.0_f64;
        let original = orbital_source_time(input, false).unwrap();
        assert_eq!(original.hi.to_bits(), input.to_bits());
        assert_eq!(original.lo.to_bits(), 0.0_f64.to_bits());
        assert_eq!(original.value().to_bits(), 0.0_f64.to_bits());
        assert_ne!(original.value().to_bits(), input.to_bits());
        // Raw acquisition and reducer use the same Time observable, rather
        // than treating its original negative-zero caller word as its epoch.
        let fresh = orbital_source_time(original.value(), false).unwrap();
        assert_eq!(fresh.value().to_bits(), original.value().to_bits());
    }

    #[test]
    fn public_object_refusal_precedes_time_or_registry() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        for is_ut in [false, true] {
            assert!(matches!(
                native_jpl_orbital_elements(
                    &ctx,
                    DeltaTModel::EspenakMeeus,
                    false,
                    f64::NAN,
                    is_ut,
                    crate::constants::SUN,
                    0
                ),
                Err(Error::InvalidOrbitalObject { .. })
            ));
        }
    }

    #[test]
    fn native_alias_refusal_names_target_before_missing_center() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        assert!(matches!(
            native_jpl_orbital_elements(
                &ctx,
                DeltaTModel::EspenakMeeus,
                false,
                2451545.125,
                false,
                crate::constants::MARS,
                0
            ),
            Err(Error::SpkUnregisteredBody {
                body: crate::constants::MARS
            })
        ));
        assert!(matches!(
            native_jpl_orbital_elements(
                &ctx,
                DeltaTModel::EspenakMeeus,
                false,
                2451545.125,
                false,
                crate::constants::EARTH,
                0
            ),
            Err(Error::SpkBodyNotFound { target: 3 })
        ));
    }

    #[test]
    fn direct_emb_feed_does_not_reconstruct_from_earth_and_moon() {
        struct NoStoredReads;
        impl StateProvider for NoStoredReads {
            fn provider_name(&self) -> &'static str {
                "test/direct-emb"
            }
            fn has_body(&self, _: i32) -> bool {
                panic!("direct feed owns aliases")
            }
            fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
                panic!("direct feed owns domains")
            }
            fn serving_coord_type(&self, _: i32) -> Option<u32> {
                panic!("direct feed owns frame")
            }
            fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
                panic!("direct feed owns entry identity")
            }
            fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>)> {
                panic!("no stored Earth/Moon reads")
            }
            fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
                panic!("no stored Earth/Moon reads")
            }
        }
        let provider = NoStoredReads;
        let nutation = |_: f64| -> Result<(f64, f64)> {
            panic!("fixed J2000 geometry does not acquire nutation")
        };
        let serving = SealedServing::new(&provider, &nutation);
        let delta = |_: f64| -> Result<f64> { panic!("TT geometry does not acquire Delta-T") };
        let tt = 2451545.125;
        let radius = 1.25_f64;
        let gm = crate::orbital::planet_system_gm(crate::constants::EARTH).unwrap();
        let calls = std::cell::RefCell::new(Vec::new());
        let raw = |body, epoch| {
            calls.borrow_mut().push(body);
            assert_eq!(body, crate::constants::SUN);
            Ok(crate::positions::GeometricState {
                body_id: body,
                jd_tt: epoch,
                pos: [0.0; 3],
                vel: [0.0; 3],
                frame: crate::positions::StateFrame::IcrsBarycentric,
                planet_center: false,
                tier: None,
            })
        };
        let emb = |epoch| {
            calls.borrow_mut().push(3);
            Ok(crate::positions::GeometricState {
                body_id: crate::constants::EARTH,
                jd_tt: epoch,
                pos: [radius, 0.0, 0.0],
                vel: [0.0, (gm / radius).sqrt(), 0.0],
                frame: crate::positions::StateFrame::IcrsSystemBarycenter,
                planet_center: false,
                tier: None,
            })
        };
        let target = |body| {
            assert_eq!(body, crate::constants::EARTH);
            calls.borrow_mut().push(-3);
            Ok(())
        };
        let orbital = crate::orbital::SealedOrbitalServing::new(&serving, &delta)
            .with_orbital_geometry(&raw, &emb, &target);
        let elements = orbital
            .orbital_elements_tt(crate::constants::EARTH, tt, 0)
            .unwrap();
        assert!((elements[0] - radius).abs() < 1e-12);
        assert_eq!(*calls.borrow(), vec![-3, crate::constants::SUN, 3]);
    }

    #[test]
    fn explicit_leb_never_crosses_into_native_orbital_acquisition() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert!(matches!(
            native_jpl_orbital_elements(
                &ctx,
                DeltaTModel::EspenakMeeus,
                false,
                2451545.125,
                false,
                crate::constants::MARS,
                0
            ),
            Err(Error::ConfigurationMissing {
                missing_config: "jpl_source_access",
                ..
            })
        ));
    }
}

#[cfg(test)]
mod native_nodaps_context_tests {
    use super::*;
    use crate::constants::*;
    use crate::positions::{GeometricState, StateFrame};
    use std::cell::RefCell;

    /// Native feeds must not accidentally acquire any stored point channel.
    struct NoStoredNodeReads;
    impl StateProvider for NoStoredNodeReads {
        fn serving_entry_snapshot(&self, _: i32) -> Result<EntrySnapshot> {
            panic!("native geometry owns entry identity")
        }
        fn provider_name(&self) -> &'static str {
            "test/native-nodaps-feeds"
        }
        fn has_body(&self, _: i32) -> bool {
            panic!("native alias feed owns lookup")
        }
        fn body_coverage(&self, _: i32) -> Option<(f64, f64)> {
            panic!("native geometry owns coverage")
        }
        fn serving_coord_type(&self, _: i32) -> Option<u32> {
            panic!("native geometry owns frame")
        }
        fn eval_raw(&self, _: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("no stored geometry")
        }
        fn eval_raw_two_part(&self, _: i32, _: TwoPartEpoch) -> Result<(Vec<f64>, Vec<f64>)> {
            panic!("no stored geometry")
        }
    }

    /// Finite test ICRS state with explicit body and scalar epoch identity.
    fn state(body: i32, tt: f64, p: [f64; 3], v: [f64; 3]) -> GeometricState {
        GeometricState {
            body_id: body,
            jd_tt: tt,
            pos: p,
            vel: v,
            frame: StateFrame::IcrsBarycentric,
            planet_center: false,
            tier: None,
        }
    }

    #[test]
    fn failed_node_support_owns_its_actual_segment_after_registry_clear() {
        let mut image = crate::spk::moon_test_kernel_bytes();
        for (index, target) in [10_i32, 399, 3].into_iter().enumerate() {
            let offset = 1048 + index * 40;
            image[offset + 16..offset + 20].copy_from_slice(&target.to_le_bytes());
            image[offset + 20..offset + 24].copy_from_slice(&0_i32.to_le_bytes());
        }
        // Deliberately wider descriptor than the actual coefficient directory:
        // the query reaches a selected segment, whose real directory refuses.
        image[1048..1056].copy_from_slice(&(-50_000_000.0_f64).to_le_bytes());
        // The index fixture has zeroed coefficient storage. Supply one valid
        // type-2 record and its directory before testing a real range refusal.
        image[1080..1084].copy_from_slice(&385_i32.to_le_bytes());
        image[1084..1088].copy_from_slice(&393_i32.to_le_bytes());
        for (index, word) in [
            5_000_000.0_f64,
            5_000_000.0,
            1.0,
            0.0,
            0.0,
            0.0,
            10_000_000.0,
            5.0,
            1.0,
        ]
        .into_iter()
        .enumerate()
        {
            let offset = 3072 + index * 8;
            image[offset..offset + 8].copy_from_slice(&word.to_le_bytes());
        }
        let mut registry = crate::spk::SpkRegistry::new();
        registry
            .register_image("node-support.bsp", image.clone())
            .unwrap();
        registry.register_body(SUN, "node-support.bsp", 10).unwrap();
        let error = {
            let provider =
                positions::LocalSpkStateProvider::new(&registry, ProviderMode::NativeJpl).unwrap();
            provider
                .nodaps_state_two_part(SUN, orbital_source_time(2451000.0, false).unwrap())
                .unwrap_err()
        };
        assert!(matches!(
            error.primary_error(),
            Error::SpkOutOfRange {
                target: 10,
                start: 0.0,
                end: 10_000_000.0,
                ..
            }
        ));
        let raw = error
            .backend_range()
            .unwrap_or_else(|| panic!("actual directory refusal retained: {error:?}"));
        assert_eq!(raw.descriptor.target, 10);
        assert_eq!(raw.segment_ordinal, 0);
        registry.clear();
        assert_eq!(raw.image.bytes(), image);
        assert_eq!(raw.image.path, "node-support.bsp");
    }

    #[test]
    fn original_ut_and_scalar_speed_sides_keep_distinct_time_parts() {
        for ut in [2451545.125, 2460400.25] {
            let original = orbital_source_time(ut, true).unwrap();
            let tt = original.value();
            let central = nodaps_source_time_at(original, tt).unwrap();
            assert_eq!(central.hi.to_bits(), ut.to_bits());
            assert_eq!(central.lo.to_bits(), original.lo.to_bits());
            assert_eq!(
                nodaps_source_ut1_at(original, ut, true, tt)
                    .unwrap()
                    .to_bits(),
                ut.to_bits()
            );
            for side in [tt - 0.001, tt + 0.001] {
                let sample = nodaps_source_time_at(original, side).unwrap();
                assert_eq!(sample.hi.to_bits(), side.floor().to_bits());
                assert_eq!(sample.lo.to_bits(), (side - side.floor()).to_bits());
                assert_eq!(
                    nodaps_source_ut1_at(original, ut, true, side)
                        .unwrap()
                        .to_bits(),
                    crate::time::default_timescale_ut1_from_tt(side)
                        .unwrap()
                        .to_bits()
                );
            }
            assert_eq!(
                nodaps_source_ut1_at(original, tt, false, tt)
                    .unwrap()
                    .to_bits(),
                crate::time::default_timescale_ut1_from_tt(tt)
                    .unwrap()
                    .to_bits()
            );
        }
    }

    #[test]
    fn missing_support_alias_precedes_the_mean_target() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        assert!(matches!(
            native_jpl_nod_aps(
                &ctx,
                DeltaTModel::EspenakMeeus,
                false,
                2451545.125,
                false,
                MARS,
                NODBIT_MEAN,
                FLG_TRUEPOS
            ),
            Err(Error::SpkUnregisteredBody { body: SUN })
        ));
    }

    #[test]
    fn undefined_and_unsupported_bodies_precede_registry_acquisition() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        for body in [INTP_APOG, INTP_PERG] {
            for flags in [0, FLG_SPEED] {
                let output = native_jpl_nod_aps(
                    &ctx,
                    DeltaTModel::EspenakMeeus,
                    false,
                    2451545.125,
                    true,
                    body,
                    NODBIT_OSCU,
                    flags,
                )
                .unwrap();
                for point in output {
                    assert!(point[..3].iter().all(|value| value.is_nan()));
                    for rate in &point[3..] {
                        if flags & FLG_SPEED != 0 {
                            assert!(rate.is_nan());
                        } else {
                            assert_eq!(rate.to_bits(), 0.0_f64.to_bits());
                        }
                    }
                }
            }
            assert!(matches!(
                native_jpl_nod_aps(
                    &ctx,
                    DeltaTModel::EspenakMeeus,
                    false,
                    2451545.125,
                    false,
                    body,
                    1,
                    FLG_TOPOCTR
                ),
                Err(Error::ConfigurationMissing {
                    missing_config: "observer_location",
                    ..
                })
            ));
        }
        assert!(matches!(
            native_jpl_nod_aps(
                &ctx,
                DeltaTModel::EspenakMeeus,
                false,
                2451545.125,
                false,
                CUPIDO,
                1,
                0
            ),
            Err(Error::InvalidOrbitalObject { .. })
        ));
    }

    #[test]
    fn explicit_leb_refuses_native_node_acquisition() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert!(matches!(
            native_jpl_nod_aps(
                &ctx,
                DeltaTModel::EspenakMeeus,
                false,
                2451545.125,
                false,
                MARS,
                1,
                FLG_TRUEPOS
            ),
            Err(Error::ConfigurationMissing {
                missing_config: "jpl_source_access",
                ..
            })
        ));
    }

    #[test]
    fn native_raw_deflectors_are_optional_and_have_explicit_body_identity() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        let provider =
            positions::LocalSpkStateProvider::new(ctx.spk_registry(), ctx.mode()).unwrap();
        for body in [SUN, JUPITER, SATURN] {
            assert!(
                provider
                    .nodaps_deflector_position(body, 2451545.125)
                    .unwrap()
                    .is_none()
            );
        }
        assert!(matches!(
            provider.nodaps_deflector_position(MARS, 2451545.125),
            Err(Error::UnknownBody { body_id: MARS })
        ));
        assert!(matches!(
            provider.nodaps_deflector_position(SUN, f64::NAN),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn osculating_earth_uses_direct_emb_after_geocenter_observation() {
        let provider = NoStoredNodeReads;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("native NONUT has no stored frame") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("ordinary geometric TT has no model clock") };
        let calls = RefCell::new(Vec::new());
        let raw = |body, tt| {
            calls.borrow_mut().push(("read", body));
            match body {
                SUN => Ok(state(body, tt, [0.0; 3], [0.0; 3])),
                EARTH => Ok(state(body, tt, [1.0, 0.0, 0.0], [0.0, 0.017, 0.0])),
                _ => panic!("Earth OSCU must not reconstruct EMB from Moon"),
            }
        };
        let radius = 1.25_f64;
        let gm = crate::orbital::planet_system_gm(EARTH).unwrap();
        let emb = |tt| {
            calls.borrow_mut().push(("emb", 3));
            Ok(state(
                EARTH,
                tt,
                [radius, 0.0, 0.0],
                [0.0, (gm / radius).sqrt(), 0.0],
            ))
        };
        let lookup = |body| {
            calls.borrow_mut().push(("alias", body));
            Ok(())
        };
        let ut = |_: f64| -> Result<f64> { panic!("non-sidereal planet never acquires UT") };
        let calc = |_: f64, _: i32, _: i32| -> Result<(PositionResult, i32)> {
            panic!("no public model request")
        };
        let defl = |_: i32, _: f64| -> Result<Option<[f64; 3]>> {
            panic!("TRUEPOS suppresses deflection")
        };
        let engine = crate::orbital::SealedOrbitalServing::new(&serving, &delta)
            .with_nodaps_acquisition(&raw, &emb, &lookup, &ut, &calc, &defl);
        let output = engine
            .nod_aps_tt(
                EARTH,
                2451545.125,
                NODBIT_OSCU,
                FLG_TRUEPOS | FLG_NONUT | FLG_HELCTR,
            )
            .unwrap();
        assert_eq!(
            *calls.borrow(),
            vec![
                ("alias", SUN),
                ("alias", EARTH),
                ("read", SUN),
                ("read", EARTH),
                ("emb", 3)
            ]
        );
        assert!(
            output[0]
                .iter()
                .all(|word| word.to_bits() == 0.0_f64.to_bits())
        );
        assert!(
            output[1]
                .iter()
                .all(|word| word.to_bits() == 0.0_f64.to_bits())
        );
        assert!((output[2][2] - radius).abs() < 1e-12);
        assert!((output[3][2] - radius).abs() < 1e-12);
    }

    #[test]
    fn observer_aliases_resolve_before_a_failing_sun_epoch_read() {
        let provider = NoStoredNodeReads;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("NONUT has no stored frame") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("no model clock") };
        let raw = |_: i32, _: f64| -> Result<GeometricState> {
            panic!("Earth alias precedes Sun epoch read")
        };
        let emb = |_: f64| -> Result<GeometricState> { panic!("no EMB read") };
        let calls = RefCell::new(Vec::new());
        let lookup = |body| {
            calls.borrow_mut().push(body);
            if body == EARTH {
                Err(Error::SpkUnregisteredBody { body })
            } else {
                Ok(())
            }
        };
        let ut = |_: f64| -> Result<f64> { panic!("no UT read") };
        let calc =
            |_: f64, _: i32, _: i32| -> Result<(PositionResult, i32)> { panic!("no model read") };
        let defl = |_: i32, _: f64| -> Result<Option<[f64; 3]>> { panic!("no deflection read") };
        let engine = crate::orbital::SealedOrbitalServing::new(&serving, &delta)
            .with_nodaps_acquisition(&raw, &emb, &lookup, &ut, &calc, &defl);
        assert!(matches!(
            engine.nod_aps_tt(MARS, 2451545.125, NODBIT_MEAN, FLG_TRUEPOS | FLG_NONUT),
            Err(Error::SpkUnregisteredBody { body: EARTH })
        ));
        assert_eq!(*calls.borrow(), vec![SUN, EARTH]);
    }

    #[test]
    fn missing_moon_alias_preserves_source_zero_tuple_before_point_callbacks() {
        let provider = NoStoredNodeReads;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("NONUT has no stored frame") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("no model clock") };
        let reads = RefCell::new(Vec::new());
        let raw = |body, tt| {
            reads.borrow_mut().push(body);
            Ok(state(body, tt, [1.0, 0.0, 0.0], [0.0; 3]))
        };
        let emb = |_: f64| -> Result<GeometricState> { panic!("no EMB read") };
        let lookup = |body| {
            if body == MOON {
                Err(Error::SpkUnregisteredBody { body })
            } else {
                Ok(())
            }
        };
        let ut = |_: f64| -> Result<f64> { panic!("no UT read after missing Moon") };
        let calc = |_: f64, _: i32, _: i32| -> Result<(PositionResult, i32)> {
            panic!("missing Moon precedes point callbacks")
        };
        let defl = |_: i32, _: f64| -> Result<Option<[f64; 3]>> { panic!("no deflection") };
        let engine = crate::orbital::SealedOrbitalServing::new(&serving, &delta)
            .with_nodaps_acquisition(&raw, &emb, &lookup, &ut, &calc, &defl);
        let output = engine
            .nod_aps_tt(MOON, 2451545.125, NODBIT_MEAN, FLG_NONUT | FLG_TOPOCTR)
            .unwrap();
        assert_eq!(*reads.borrow(), vec![SUN, EARTH]);
        assert!(
            output
                .iter()
                .flatten()
                .all(|word| word.to_bits() == 0.0_f64.to_bits())
        );
    }

    #[test]
    fn moon_queries_genuine_public_ut_after_supports_with_source_flags() {
        let provider = NoStoredNodeReads;
        let nut = |_: f64| -> Result<(f64, f64)> { panic!("native NONUT has no stored frame") };
        let serving = SealedServing::new(&provider, &nut);
        let delta = |_: f64| -> Result<f64> { panic!("public lunar callback owns model clocks") };
        let calls = RefCell::new(Vec::new());
        let raw = |body, tt| {
            calls.borrow_mut().push(("read", body));
            Ok(state(body, tt, [1.0, 0.0, 0.0], [0.0; 3]))
        };
        let emb = |_: f64| -> Result<GeometricState> { panic!("mean Moon does not use EMB") };
        let lookup = |body| {
            calls.borrow_mut().push(("alias", body));
            Ok(())
        };
        let ut = |_: f64| Ok(2451545.125);
        let calc = |jd: f64, body, flags| {
            assert_eq!(jd.to_bits(), 2451545.125_f64.to_bits());
            assert_eq!(flags, FLG_NONUT | FLG_TRUEPOS);
            calls.borrow_mut().push(("calc_ut", body));
            Ok((
                if body == MEAN_NODE {
                    [10.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                } else {
                    assert_eq!(body, MEAN_APOG);
                    [30.0, 2.0, 0.0027, 0.0, 0.0, 0.0]
                },
                flags,
            ))
        };
        let defl = |_: i32, _: f64| -> Result<Option<[f64; 3]>> {
            panic!("Moon point path does not deflect")
        };
        let engine = crate::orbital::SealedOrbitalServing::new(&serving, &delta)
            .with_nodaps_acquisition(&raw, &emb, &lookup, &ut, &calc, &defl);
        let flags = FLG_NONUT
            | FLG_TRUEPOS
            | FLG_J2000
            | FLG_HELCTR
            | FLG_BARYCTR
            | FLG_XYZ
            | FLG_EQUATORIAL
            | FLG_RADIANS;
        let output = engine
            .nod_aps_tt(MOON, 2451545.125, NODBIT_MEAN, flags)
            .unwrap();
        assert_eq!(
            *calls.borrow(),
            vec![
                ("alias", SUN),
                ("alias", EARTH),
                ("read", SUN),
                ("read", EARTH),
                ("alias", MOON),
                ("read", MOON),
                ("calc_ut", MEAN_NODE),
                ("calc_ut", MEAN_APOG)
            ]
        );
        assert!(output.iter().flatten().all(|value| value.is_finite()));
    }
}

#[cfg(test)]
mod native_nodaps_site_review_tests {
    use super::*;
    use crate::constants::*;
    use crate::positions::{GeometricState, StateFrame};
    use std::cell::RefCell;

    fn earth(tt: f64) -> GeometricState {
        GeometricState {
            body_id: EARTH,
            jd_tt: tt,
            pos: [1.0, 2.0, 3.0],
            vel: [0.01, 0.02, 0.03],
            frame: StateFrame::IcrsBarycentric,
            planet_center: false,
            tier: None,
        }
    }

    #[test]
    fn fresh_site_time_reads_topocenter_before_geocenter_with_roundtrip_rounding() {
        let original = orbital_source_time(2460400.25, true).unwrap();
        let central = original.value();
        for tt in [central, central - 0.0001, central + 0.0001] {
            let calls = RefCell::new(Vec::new());
            let offset = [0.00003, 0.00002, 0.00001];
            let spin = [0.0001, 0.0002, 0.0003];
            let raw = |time: TwoPartEpoch| {
                calls.borrow_mut().push("earth");
                assert_eq!(time.hi.to_bits(), tt.floor().to_bits());
                assert_eq!(time.lo.to_bits(), (tt - tt.floor()).to_bits());
                Ok(earth(time.value()))
            };
            let topocenter = |time: TwoPartEpoch, site: (f64, f64, f64)| {
                calls.borrow_mut().push("topos");
                assert_eq!(time.hi.to_bits(), tt.floor().to_bits());
                assert_eq!(time.lo.to_bits(), (tt - tt.floor()).to_bits());
                assert_eq!(site, (13.405, 52.52, 34.0));
                Ok((offset, spin))
            };
            let (position, velocity) =
                nodaps_site_state_via(tt, (13.405, 52.52, 34.0), raw, topocenter).unwrap();
            assert_eq!(*calls.borrow(), vec!["earth", "topos", "earth"]);
            for i in 0..3 {
                assert_eq!(
                    position[i].to_bits(),
                    ((earth(tt).pos()[i] + offset[i]) - earth(tt).pos()[i]).to_bits()
                );
                assert_eq!(
                    velocity[i].to_bits(),
                    ((earth(tt).vel()[i] + spin[i]) - earth(tt).vel()[i]).to_bits()
                );
            }
        }
    }

    #[test]
    fn raw_site_support_error_precedes_offset_and_has_original_domain_identity() {
        let raw = |_: TwoPartEpoch| -> Result<GeometricState> {
            Err(Error::SpkOutOfRange {
                target: 399,
                et: 4.0,
                start: 0.0,
                end: 3.0,
            })
        };
        let topocenter = |_: TwoPartEpoch, _: (f64, f64, f64)| -> Result<([f64; 3], [f64; 3])> {
            panic!("failing Earth precedes site-frame work")
        };
        assert!(matches!(
            nodaps_site_state_via(2451545.0, (0.0, 0.0, 0.0), raw, topocenter),
            Err(Error::SpkOutOfRange {
                target: 399,
                et: 4.0,
                start: 0.0,
                end: 3.0
            })
        ));
    }

    #[test]
    fn configured_native_site_defers_global_inverse_until_source_support_reads() {
        let ctx = EphemerisContext::new()
            .with_mode(ProviderMode::NativeJpl)
            .with_observer(crate::context::Observer::new(13.405, 52.52, 34.0).unwrap());
        assert!(crate::time::default_timescale_ut1_from_tt(f64::MAX).is_err());
        for body in [CUPIDO, INTP_APOG] {
            let result = native_jpl_nod_aps(
                &ctx,
                DeltaTModel::EspenakMeeus,
                false,
                f64::MAX,
                false,
                body,
                NODBIT_MEAN,
                FLG_TOPOCTR,
            );
            if body == CUPIDO {
                assert!(matches!(result, Err(Error::InvalidOrbitalObject { .. })));
            } else {
                assert!(matches!(
                    result,
                    Err(Error::SpkUnregisteredBody { body: EARTH })
                ));
            }
        }
    }

    #[test]
    fn configured_native_undefined_apsides_require_earth_but_missing_site_has_priority() {
        let empty = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        let configured = empty
            .clone()
            .with_observer(crate::context::Observer::new(13.405, 52.52, 34.0).unwrap());
        for body in [INTP_APOG, INTP_PERG] {
            for flags in [
                FLG_TOPOCTR,
                FLG_TOPOCTR | FLG_TRUEPOS,
                FLG_TOPOCTR | FLG_NOABERR,
            ] {
                assert!(matches!(
                    native_jpl_nod_aps(
                        &empty,
                        DeltaTModel::EspenakMeeus,
                        false,
                        2460400.25,
                        true,
                        body,
                        1,
                        flags
                    ),
                    Err(Error::ConfigurationMissing {
                        missing_config: "observer_location",
                        ..
                    })
                ));
                assert!(matches!(
                    native_jpl_nod_aps(
                        &configured,
                        DeltaTModel::EspenakMeeus,
                        false,
                        2460400.25,
                        true,
                        body,
                        1,
                        flags
                    ),
                    Err(Error::SpkUnregisteredBody { body: EARTH })
                ));
            }
        }
    }

    #[test]
    fn all_thirteen_live_anchor_modes_read_actual_native_support() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        let provider =
            positions::LocalSpkStateProvider::new(ctx.spk_registry(), ctx.mode()).unwrap();
        let anchors = crate::stars::SealedAnchorProvider::new(&provider);
        let delta = |_: f64| Ok(0.25);
        let mut count = 0;
        for mode in crate::sidereal::DYNAMIC_AYANAMSHA_MODES {
            if [SIDM_GALALIGN_MARDYKS, SIDM_VALENS_MOON].contains(&mode) {
                continue;
            }
            let mut selection = crate::sidereal::SiderealSelection::unset();
            assert!(selection.set_sidereal_mode(mode, 0.0, 0.0).is_empty());
            let aya =
                crate::sidereal::SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
            assert!(
                matches!(
                    aya.mean_tt_with(2460400.25, false, false),
                    Err(Error::SpkUnregisteredBody { body: EARTH })
                ),
                "mode {mode}"
            );
            count += 1;
        }
        assert_eq!(count, 13);
    }
}

#[cfg(test)]
mod native_pheno_context_contract_tests {
    use super::*;

    #[test]
    fn public_pheno_stages_complete_before_derived_target_validation() {
        use std::cell::RefCell;
        for helio in [None, Some(crate::constants::FLG_HELCTR)] {
            let calls = RefCell::new(Vec::new());
            let acquire = |ut: f64, body, flags| {
                assert_eq!(ut.to_bits(), 2460400.25f64.to_bits());
                calls.borrow_mut().push((body, flags));
                if body == crate::constants::SUN {
                    Err(Error::UnknownBody { body_id: 991 })
                } else {
                    Ok([f64::NAN, 0.0, 1.0, 0.0, 0.0, 0.0])
                }
            };
            assert!(matches!(
                native_pheno_public_states(
                    2460400.25,
                    crate::constants::CUPIDO,
                    0,
                    helio,
                    &acquire
                ),
                Err(Error::UnknownBody { body_id: 991 })
            ));
            assert_eq!(
                *calls.borrow(),
                [(crate::constants::CUPIDO, 0), (crate::constants::SUN, 0)]
            );
        }
    }

    #[test]
    fn public_asteroid_third_acquisition_precedes_derived_target_validation() {
        use std::cell::RefCell;
        let calls = RefCell::new(Vec::new());
        let acquire = |_: f64, body, flags| {
            calls.borrow_mut().push((body, flags));
            if flags == crate::constants::FLG_HELCTR {
                Err(Error::UnknownBody { body_id: 992 })
            } else {
                Ok([f64::NAN, 0.0, 1.0, 0.0, 0.0, 0.0])
            }
        };
        assert!(matches!(
            native_pheno_public_states(
                2460400.25,
                crate::constants::CERES,
                0,
                Some(crate::constants::FLG_HELCTR),
                &acquire
            ),
            Err(Error::UnknownBody { body_id: 992 })
        ));
        assert_eq!(
            *calls.borrow(),
            [
                (crate::constants::CERES, 0),
                (crate::constants::SUN, 0),
                (crate::constants::CERES, crate::constants::FLG_HELCTR)
            ]
        );
    }

    #[test]
    fn native_pheno_keeps_explicit_leb_gate_and_finite_entry_guard() {
        let sealed = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert!(
            native_jpl_pheno(
                &sealed,
                DeltaTModel::EspenakMeeus,
                false,
                2451545.125,
                false,
                crate::constants::EARTH,
                0
            )
            .is_err()
        );
        let native = sealed.with_mode(ProviderMode::NativeJpl);
        for is_ut in [false, true] {
            assert!(matches!(
                native_jpl_pheno(
                    &native,
                    DeltaTModel::EspenakMeeus,
                    false,
                    f64::NAN,
                    is_ut,
                    crate::constants::EARTH,
                    0
                ),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
    }

    #[test]
    fn native_pheno_degenerate_sentinels_need_no_registry_model_or_site() {
        let native = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        for is_ut in [false, true] {
            for flags in [
                0,
                crate::constants::FLG_TOPOCTR
                    | crate::constants::FLG_TRUEPOS
                    | crate::constants::FLG_SIDEREAL
                    | crate::constants::FLG_XYZ,
            ] {
                let values = native_jpl_pheno(
                    &native,
                    DeltaTModel::Smh2016,
                    true,
                    2451545.125,
                    is_ut,
                    crate::constants::EARTH,
                    flags,
                )
                .unwrap();
                for (i, value) in values.iter().enumerate() {
                    assert_eq!(
                        value.to_bits(),
                        if i == 3 { 180.0_f64.to_bits() } else { 0 }
                    );
                }
                let values = native_jpl_pheno(
                    &native,
                    DeltaTModel::Smh2016,
                    true,
                    2451545.125,
                    is_ut,
                    crate::constants::ECL_NUT,
                    flags,
                )
                .unwrap();
                assert!(values[..3].iter().all(|value| value.is_nan()));
                assert!(values[3..].iter().all(|value| value.to_bits() == 0));
            }
        }
    }

    #[test]
    fn native_pheno_earth_alias_precedes_observer_and_target() {
        let native = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        assert!(matches!(
            native_jpl_pheno(
                &native,
                DeltaTModel::EspenakMeeus,
                false,
                2451545.125,
                true,
                crate::constants::MARS,
                crate::constants::FLG_TOPOCTR
            ),
            Err(Error::SpkUnregisteredBody {
                body: crate::constants::EARTH
            })
        ));
    }

    #[test]
    fn native_pheno_unknown_family_does_not_consume_extreme_tt_clock() {
        use crate::constants::{EARTH, FLG_TOPOCTR, SUN};
        let mut native = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        for (body, naif) in [(EARTH, 399_i32), (SUN, 10_i32)] {
            let mut image = crate::spk::moon_test_kernel_bytes();
            let target_offset =
                crate::spk::RECORD_BYTES + (3 + 2 * crate::spk::SPK_SUMMARY_WORDS) * 8 + 16;
            image[target_offset..target_offset + 4].copy_from_slice(&naif.to_le_bytes());
            let name = format!("pheno-alias-{body}");
            native
                .spk_registry_mut()
                .register_image(&name, image)
                .unwrap();
            native
                .spk_registry_mut()
                .register_body(body, &name, naif)
                .unwrap();
        }
        // These images deliberately have no numerical directory: positional
        // refusal requires alias resolution, but neither support states nor a
        // Delta-T evaluation whose extrapolation overflows at this epoch.
        for epoch in [2451545.125, 1e308] {
            for body in [-999, 59, 999] {
                assert_eq!(
                    native_jpl_pheno(
                        &native,
                        DeltaTModel::EspenakMeeus,
                        false,
                        epoch,
                        false,
                        body,
                        0,
                    ),
                    Err(Error::UnknownBody { body_id: body }),
                );
                assert!(matches!(
                    native_jpl_pheno(
                        &native,
                        DeltaTModel::EspenakMeeus,
                        false,
                        epoch,
                        false,
                        body,
                        FLG_TOPOCTR,
                    ),
                    Err(Error::ConfigurationMissing {
                        missing_config: "observer_location",
                        ..
                    }),
                ));
            }
        }
        let absent = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        assert_eq!(
            native_jpl_pheno(
                &absent,
                DeltaTModel::EspenakMeeus,
                false,
                1e308,
                false,
                999,
                0
            ),
            Err(Error::SpkUnregisteredBody { body: EARTH }),
        );
    }

    #[test]
    fn native_pheno_actual_observation_precedes_public_callback_without_core_ut1_guard() {
        use crate::constants::{EARTH, FLG_TOPOCTR, SUN};
        use crate::engine_api::PhenoCallbacks;
        use std::cell::RefCell;
        let mut native = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        for (body, naif) in [(EARTH, 399_i32), (SUN, 10_i32)] {
            let mut image = crate::spk::moon_test_kernel_bytes();
            let target_offset =
                crate::spk::RECORD_BYTES + (3 + 2 * crate::spk::SPK_SUMMARY_WORDS) * 8 + 16;
            image[target_offset..target_offset + 4].copy_from_slice(&naif.to_le_bytes());
            let name = format!("owned-pheno-alias-{body}");
            native
                .spk_registry_mut()
                .register_image(&name, image)
                .unwrap();
            native
                .spk_registry_mut()
                .register_body(body, &name, naif)
                .unwrap();
        }
        let events = RefCell::new(Vec::new());
        let ut1 = || {
            events.borrow_mut().push("observation.ut1");
            Ok(f64::NEG_INFINITY)
        };
        let position = |epoch: f64, body: i32, flags: i32| {
            assert_eq!(epoch, f64::NEG_INFINITY);
            assert_eq!(body, 999);
            assert_eq!(flags, 0);
            events.borrow_mut().push("public.calc_ut");
            Err(Error::ConfigurationMissing {
                missing_config: "actual_public_callback",
                suggestion: "retain this source-owned error",
            })
        };
        let callbacks = PhenoCallbacks::new(&position).with_observation_ut1(&ut1);
        let error = native_jpl_pheno_with_callbacks(
            &native,
            DeltaTModel::EspenakMeeus,
            false,
            1e308,
            false,
            999,
            0,
            Some(&callbacks),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            Error::ConfigurationMissing {
                missing_config: "actual_public_callback",
                ..
            }
        ));
        assert_eq!(*events.borrow(), ["observation.ut1", "public.calc_ut"]);
        events.borrow_mut().clear();
        let error = native_jpl_pheno_with_callbacks(
            &native,
            DeltaTModel::EspenakMeeus,
            false,
            1e308,
            false,
            999,
            FLG_TOPOCTR,
            Some(&callbacks),
        )
        .unwrap_err();
        assert!(matches!(
            error,
            Error::ConfigurationMissing {
                missing_config: "observer_location",
                ..
            }
        ));
        assert!(events.borrow().is_empty());
    }
}

#[cfg(test)]
mod native_pheno_review_domain_tests {
    use super::*;

    #[test]
    fn finite_extreme_original_tt_preserves_earth_alias_priority() {
        let ctx = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        let epoch = 1e308;
        let observation = NativePhenoObservation::from_entry(epoch, false).unwrap();
        assert_eq!(observation.split_tt.value().to_bits(), epoch.to_bits());
        assert!(observation.original_ut1_fraction.is_none());
        // Missing Earth must precede an unused long-term Delta-T overflow.
        // This is the original finite TT input, never a normalized surrogate.
        assert!(matches!(
            native_jpl_pheno(
                &ctx,
                DeltaTModel::EspenakMeeus,
                false,
                epoch,
                false,
                crate::constants::MARS,
                crate::constants::FLG_TOPOCTR
            ),
            Err(Error::SpkUnregisteredBody {
                body: crate::constants::EARTH
            })
        ));
    }

    #[test]
    fn original_ut_time_keeps_refined_tt_and_exact_zero_ut1_fraction() {
        for epoch in [2451545.125, 2460400.25] {
            let observation = NativePhenoObservation::from_entry(epoch, true).unwrap();
            let (whole, refined) = crate::smh2016::timescale_parts_from_ut1(epoch).unwrap();
            assert_eq!(observation.split_tt.hi.to_bits(), whole.to_bits());
            assert_eq!(observation.split_tt.lo.to_bits(), refined.to_bits());
            assert_eq!(
                observation.original_ut1_fraction.unwrap().to_bits(),
                0.0_f64.to_bits()
            );
            assert_eq!(observation.ut1().unwrap().to_bits(), epoch.to_bits());
        }
    }

    #[test]
    fn iers2010_equatorial_zero_site_preserves_accepted_finite_height() {
        let height = -positions::IERS2010_SEMI_MAJOR_M;
        let _observer = crate::context::Observer::new(0.0, 0.0, height).unwrap();
        let vector =
            native_pheno_itrs_position((0.0, 0.0, height), positions::SiteEllipsoid::Iers2010)
                .unwrap();
        assert!(vector.into_iter().all(|value| value == 0.0));
        let wgs = native_pheno_itrs_position((0.0, 0.0, height), positions::SiteEllipsoid::Wgs84)
            .unwrap();
        assert!(wgs[0] > 0.0);
    }

    #[test]
    fn iers2010_midlatitude_obeys_ellipsoid_and_wgs84_stays_independent() {
        let site = (23.5, 45.0, 0.0);
        let iers = native_pheno_itrs_position(site, positions::SiteEllipsoid::Iers2010).unwrap();
        let wgs = native_pheno_itrs_position(site, positions::SiteEllipsoid::Wgs84).unwrap();
        let original_wgs =
            positions::reference_frames::itrs_position_au(site.0, site.1, site.2).unwrap();
        assert_eq!(wgs.map(f64::to_bits), original_wgs.map(f64::to_bits));
        assert_ne!(iers.map(f64::to_bits), wgs.map(f64::to_bits));
        let a = positions::IERS2010_SEMI_MAJOR_M / positions::reference_frames::FRAME_AU_M;
        let b = a * (1.0 - 1.0 / positions::IERS2010_INV_FLATTENING);
        let surface =
            (iers[0] * iers[0] + iers[1] * iers[1]) / (a * a) + iers[2] * iers[2] / (b * b);
        assert!((surface - 1.0).abs() < 16.0 * f64::EPSILON);
    }

    #[test]
    fn frame_projection_rejects_nonfinite_vectors() {
        let frame = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert!(matches!(
            native_pheno_frame_angles(frame, [f64::NAN; 3]),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn parallax_preserves_undefined_norm_overflow() {
        assert_eq!(native_pheno_parallax([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]), 0.0);
        assert_eq!(
            native_pheno_parallax([1.0, 0.0, 0.0], [0.0, 1.0, 0.0]),
            90.0
        );
        let finite = [f64::MAX, 0.0, 0.0];
        assert!(native_pheno_parallax(finite, finite).is_nan());
    }
}

#[cfg(test)]
mod astropy_core_time_tests {
    use super::*;

    #[test]
    fn astropy_core_time_owned_star_route_keeps_clock_mode_and_alias_priority() {
        let absent = EphemerisContext::new().with_mode(ProviderMode::NativeJpl);
        assert!(matches!(
            native_jpl_arbitrary_star_altaz_ut1(&absent, f64::NAN, 0.0, 0.0, 0.0, 0.0),
            Err(Error::NonFiniteArgument { .. })
        ));
        // Source obtains Earth before constructing the Star and WGS84 site.
        assert_eq!(
            native_jpl_arbitrary_star_altaz_ut1(&absent, 2451545.0, f64::NAN, 0.0, 0.0, 0.0),
            Err(Error::SpkUnregisteredBody {
                body: crate::constants::EARTH
            })
        );
        let sealed = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert!(matches!(
            native_jpl_arbitrary_star_altaz_ut1(&sealed, 2451545.0, 0.0, 0.0, 0.0, 0.0),
            Err(Error::ConfigurationMissing { .. })
        ));
    }
}
/// Bundled UT1 Time and its collapsed TT, retained for the width observer.
/// Main-image selection precedes this constructor; body aliases follow it.
pub(crate) struct NativeSolarPathTime {
    epochs: positions::reference_frames::VectorRotationEpochs,
    tt: f64,
}

/// Construct the width Time then look up Earth, Sun and Moon in that order.
/// These failures lie outside the observer/observation zero-width catch.
pub(crate) fn native_solar_path_time(
    ctx: &EphemerisContext,
    jd: f64,
) -> Result<NativeSolarPathTime> {
    let (epochs, tt) = positions::reference_frames::VectorRotationEpochs::from_ut1(jd)?;
    let provider = positions::LocalSpkStateProvider::new(ctx.spk_registry(), ctx.mode())?;
    provider.validate_source_alias(crate::constants::EARTH)?;
    provider.validate_source_alias(crate::constants::SUN)?;
    provider.validate_source_alias(crate::constants::MOON)?;
    Ok(NativeSolarPathTime { epochs, tt })
}

/// Retained apparent ICRS rays in AU and the local unrefracted horizon frame.
/// All data are owned; no provider, caller path or callback is retained.
pub(crate) struct NativeSolarPathRays {
    sun: [f64; 3],
    moon: [f64; 3],
    rotation: [[f64; 3]; 3],
    lon: f64,
    lat: f64,
}

/// Build a true WGS84 height0 observer, then observe Sun and Moon at one
/// bundled UT1 Time. The Earth surface is Some(surface) for both rays.
/// Caller position/output flags and configured model clocks are not consumed.
/// Alias/Time admission belongs to native_solar_path_time outside this arm.
pub(crate) fn native_solar_path_rays(
    ctx: &EphemerisContext,
    time: &NativeSolarPathTime,
    lon: f64,
    lat: f64,
) -> Result<NativeSolarPathRays> {
    use positions::reference_frames::*;
    let provider = positions::LocalSpkStateProvider::new(ctx.spk_registry(), ctx.mode())?;
    let itrs = native_pheno_itrs_position((lon, lat, 0.0), positions::SiteEllipsoid::Wgs84)?;
    let epochs = &time.epochs;
    let epoch = TwoPartEpoch {
        hi: epochs.whole,
        lo: epochs.tt_fraction,
    };
    let earth = provider.star_source_time_state(crate::constants::EARTH, epoch)?;
    let rotation = itrs_rotation_at_parts(
        epochs.whole,
        epochs.tt_fraction,
        epochs.ut1_fraction,
        time.tt,
    )?;
    let (surface, velocity) = itrs_to_gcrs(rotation, itrs)?;
    let observer = std::array::from_fn(|i| earth.0[i] + surface[i]);
    let observer_velocity = std::array::from_fn(|i| earth.1[i] + velocity[i]);
    let sun = provider.eclipse_source_observed_vector(
        crate::constants::SUN,
        epoch,
        observer,
        observer_velocity,
        Some(surface),
    )?;
    let moon = provider.eclipse_source_observed_vector(
        crate::constants::MOON,
        epoch,
        observer,
        observer_velocity,
        Some(surface),
    )?;
    Ok(NativeSolarPathRays {
        sun,
        moon,
        rotation,
        lon,
        lat,
    })
}

impl NativeSolarPathRays {
    /// True Sun altitude, degrees, through R_latlon @ R_ITRS. No refraction,
    /// atmospheric parameters, ecliptic frame or configured sidereal clock.
    pub(crate) fn altitude(&self) -> Result<f64> {
        use positions::reference_frames::*;
        let phi = self.lat / 360.0 * FRAME_TAU;
        let zero = phi * 0.0;
        let one = zero + 1.0;
        let latitude = [
            [-scalar_sin(phi), zero, scalar_cos(phi)],
            [zero, one, zero],
            [scalar_cos(phi), zero, scalar_sin(phi)],
        ];
        let latlon = matmul3(latitude, rot_z(-(self.lon / 360.0 * FRAME_TAU)))?;
        let horizon = matmul3(latlon, self.rotation)?;
        Ok(native_pheno_frame_angles(horizon, self.sun)?.1)
    }
    /// Left-associated observer-relative Sun and Moon vector norms, AU.
    pub(crate) fn distances(&self) -> (f64, f64) {
        (
            native_pheno_vector_length(self.sun),
            native_pheno_vector_length(self.moon),
        )
    }
    /// Stable vector separation, degrees; normalize by each opposing norm,
    /// then use twice atan2(norm(a-b),norm(a+b)) without an acos clamp.
    pub(crate) fn separation(&self) -> f64 {
        let a = self.sun.map(|v| v * native_pheno_vector_length(self.moon));
        let b = self.moon.map(|v| v * native_pheno_vector_length(self.sun));
        let sub = std::array::from_fn(|i| a[i] - b[i]);
        let add = std::array::from_fn(|i| a[i] + b[i]);
        (2.0 * native_pheno_vector_length(sub).atan2(native_pheno_vector_length(add))) * 360.0
            / positions::reference_frames::FRAME_TAU
    }
}

#[path = "context_serving/native_occultation.rs"]
mod native_occultation;
pub(crate) use native_occultation::native_occultation_topocentric_place;
pub(crate) use native_occultation::{NativePlanetObserver, PlanetObservedTarget};

#[path = "context_serving/eclipse_utilities.rs"]
mod eclipse_utilities;
pub(crate) use eclipse_utilities::NativeFixedSolarObserver;
pub(crate) use eclipse_utilities::fixed_solar_vector_separation;

/// Relabel an exterior stored-vector miss using that body's union coverage.
/// Interior tier gaps retain the original reader failure and precise stored bounds.
fn vector_reader_error(error: Error, coverage: impl FnOnce(i32) -> Option<(f64, f64)>) -> Error {
    let Error::LebOutOfRange { body_id, jd, .. } = error.primary_error() else {
        return error;
    };
    let Some((start_jd, end_jd)) = coverage(*body_id) else {
        return error;
    };
    if *jd < start_jd || *jd > end_jd {
        Error::EphemerisRange {
            requested_jd: *jd,
            start_jd,
            end_jd,
            body: Some(*body_id),
        }
        .with_cause(error)
    } else {
        error
    }
}

#[cfg(test)]
mod leb_center_offset_tests {
    use super::*;

    /// The established sampler sees the same physical states and domains, but
    /// has no analytical capability. This is independent of the bound-provider
    /// opt-in so the regression catches an accidental analytical route.
    struct SampledProvider<'a, 'b>(&'a BoundProvider<'b>);

    impl StateProvider for SampledProvider<'_, '_> {
        fn provider_name(&self) -> &'static str {
            self.0.provider_name()
        }
        fn has_body(&self, body: i32) -> bool {
            self.0.has_body(body)
        }
        fn body_coverage(&self, body: i32) -> Option<(f64, f64)> {
            self.0.body_coverage(body)
        }
        fn smooth_interval_tt(&self, body: i32, jd: f64) -> Result<Option<(f64, f64)>> {
            self.0.smooth_interval_tt(body, jd)
        }
        fn nutation_smooth_interval_tt(&self, jd: f64) -> Result<Option<(f64, f64)>> {
            self.0.nutation_smooth_interval_tt(jd)
        }
        fn serving_coord_type(&self, body: i32) -> Option<u32> {
            self.0.serving_coord_type(body)
        }
        fn serving_entry_snapshot(&self, body: i32) -> Result<EntrySnapshot> {
            self.0.serving_entry_snapshot(body)
        }
        fn eval_raw(&self, body: i32, jd: f64) -> Result<(Vec<f64>, Vec<f64>)> {
            self.0.eval_raw(body, jd)
        }
        fn eval_raw_two_part(
            &self,
            body: i32,
            epoch: TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>)> {
            self.0.eval_raw_two_part(body, epoch)
        }
        fn geometric_state(&self, body: i32, jd: f64) -> Result<positions::GeometricState> {
            self.0.geometric_state(body, jd)
        }
        fn geometric_state_two_part(
            &self,
            body: i32,
            epoch: TwoPartEpoch,
        ) -> Result<positions::GeometricState> {
            self.0.geometric_state_two_part(body, epoch)
        }
        fn eval_relative_position_two_part(
            &self,
            target: i32,
            target_epoch: TwoPartEpoch,
            center: i32,
            center_epoch: TwoPartEpoch,
        ) -> Result<Option<[f64; 3]>> {
            self.0
                .eval_relative_position_two_part(target, target_epoch, center, center_epoch)
        }
        fn geometric_acceleration_tt(&self, _: i32, _: f64) -> Result<Option<[f64; 3]>> {
            panic!("an unqualified provider must not acquire analytical acceleration")
        }
    }

    fn linear_center_registry() -> crate::spk::SpkRegistry {
        let mut image = crate::spk::moon_test_kernel_bytes();
        image[1064..1068].copy_from_slice(&599_i32.to_le_bytes());
        image[1080..1084].copy_from_slice(&385_i32.to_le_bytes());
        image[1084..1088].copy_from_slice(&396_i32.to_le_bytes());
        // One linear type-2 record, km and seconds past J2000 TDB.
        for (index, word) in [
            5_000_000.0_f64,
            5_000_000.0,
            10.0,
            2.0,
            -3.0,
            4.0,
            5.0,
            -1.0,
            0.0,
            10_000_000.0,
            8.0,
            1.0,
        ]
        .into_iter()
        .enumerate()
        {
            let offset = 3072 + index * 8;
            image[offset..offset + 8].copy_from_slice(&word.to_le_bytes());
        }
        let mut registry = crate::spk::SpkRegistry::new();
        registry.register_image("linear-centers", image).unwrap();
        registry
            .register_automatic_planet_centers("linear-centers")
            .unwrap();
        registry
    }

    #[test]
    fn stored_system_barycenter_gains_its_real_offset_and_rate_without_de() {
        let file = LebFile::open(&crate::test_data::core_leb_path()).unwrap();
        let base =
            provider_from_images(std::iter::once((file.path(), &file)), Tier::Base, false).unwrap();
        let mut registry = linear_center_registry();
        let corrected = BoundProvider::WithCenters(
            Box::new(
                provider_from_images(std::iter::once((file.path(), &file)), Tier::Base, false)
                    .unwrap(),
            ),
            registry.clone(),
        );
        // Retiring the source registry cannot invalidate an active provider.
        registry.clear();
        let jd = 2451546.0;
        // The underlying stored argument may agree, but the SPK-composed
        // provider must retain its independent second-read/error route.
        assert!(!corrected.can_reuse_scalar_state(
            crate::constants::JUPITER,
            jd,
            TwoPartEpoch::split(jd).unwrap(),
        ));
        let raw = base.geometric_state(crate::constants::JUPITER, jd).unwrap();
        let state = corrected
            .geometric_state(crate::constants::JUPITER, jd)
            .unwrap();
        assert!(state.is_planet_center());
        let tdb = crate::time::tt_to_tdb_jd(jd).unwrap();
        let tau = ((tdb - 2451545.0) * 86400.0 - 5_000_000.0) / 5_000_000.0;
        let au_km = 149_597_870.7;
        for (axis, (constant, slope)) in [(10.0, 2.0), (-3.0, 4.0), (5.0, -1.0)]
            .into_iter()
            .enumerate()
        {
            let expected_offset = (constant + slope * tau) / au_km;
            assert!((state.pos()[axis] - raw.pos()[axis] - expected_offset).abs() < 3e-15);
            let expected_rate = slope / 5_000_000.0 * 86400.0 / au_km;
            assert!((state.vel()[axis] - raw.vel()[axis] - expected_rate).abs() < 1e-16);
        }
        let split = corrected
            .geometric_state_two_part(crate::constants::JUPITER, TwoPartEpoch::split(jd).unwrap())
            .unwrap();
        assert_eq!(split.pos(), state.pos());
        assert_eq!(split.vel(), state.vel());
        // A real exterior TDB epoch retains the stored barycenter, never an
        // extrapolated or zero-filled offset state.
        let before = 2451544.9;
        assert_eq!(
            corrected
                .geometric_state(crate::constants::JUPITER, before)
                .unwrap(),
            base.geometric_state(crate::constants::JUPITER, before)
                .unwrap()
        );
    }

    #[test]
    fn composed_centers_keep_sampled_speed_despite_valid_earth_acceleration() {
        use crate::constants::{
            EARTH, FLG_EQUATORIAL, FLG_ICRS, FLG_J2000, FLG_NONUT, FLG_SPEED, JUPITER, MARS,
        };
        let file = LebFile::open(&crate::test_data::core_leb_path()).unwrap();
        assert!(
            positions::SingleLeb::from_file(&file)
                .unwrap()
                .has_analytic_kinematics()
        );
        // Declining the composed-source capability is structural: it must not
        // depend on whether this request happens to find a center correction.
        for registry in [crate::spk::SpkRegistry::new(), linear_center_registry()] {
            let base =
                provider_from_images(std::iter::once((file.path(), &file)), Tier::Base, false)
                    .unwrap();
            assert!(base.has_analytic_kinematics());
            let corrected = BoundProvider::WithCenters(Box::new(base), registry);
            assert!(!corrected.has_analytic_kinematics());
            let sampled = SampledProvider(&corrected);
            assert!(!sampled.has_analytic_kinematics());
            for jd in [2451546.0, 2451555.0] {
                assert_eq!(corrected.body_center_offset_at(EARTH, jd).unwrap(), None);
                assert!(
                    corrected
                        .geometric_acceleration_tt(EARTH, jd)
                        .unwrap()
                        .is_some()
                );
                for body in [JUPITER, MARS] {
                    for frame in [
                        0,
                        FLG_EQUATORIAL,
                        FLG_NONUT,
                        FLG_EQUATORIAL | FLG_NONUT,
                        FLG_J2000,
                        FLG_J2000 | FLG_EQUATORIAL | FLG_ICRS,
                    ] {
                        for layout in 0..3 {
                            let actual_reads = std::cell::RefCell::new(Vec::new());
                            let expected_reads = std::cell::RefCell::new(Vec::new());
                            let actual_nutation = |epoch: f64| {
                                actual_reads.borrow_mut().push(epoch.to_bits());
                                corrected.nutation_at(epoch)
                            };
                            let expected_nutation = |epoch: f64| {
                                expected_reads.borrow_mut().push(epoch.to_bits());
                                corrected.nutation_at(epoch)
                            };
                            let actual = SealedServing::new_stored(&corrected, &actual_nutation);
                            let expected = SealedServing::new_stored(&sampled, &expected_nutation);
                            let flags = FLG_SPEED | frame;
                            let actual = match layout {
                                0 => positions::serving_compute_tt(&actual, jd, body, flags),
                                1 => positions::serving_compute_xyz_tt(&actual, jd, body, flags),
                                _ => positions::serving_compute_tt_split(
                                    &actual,
                                    jd,
                                    body,
                                    flags,
                                    Some((jd, 1e-12)),
                                ),
                            }
                            .unwrap();
                            let expected = match layout {
                                0 => positions::serving_compute_tt(&expected, jd, body, flags),
                                1 => positions::serving_compute_xyz_tt(&expected, jd, body, flags),
                                _ => positions::serving_compute_tt_split(
                                    &expected,
                                    jd,
                                    body,
                                    flags,
                                    Some((jd, 1e-12)),
                                ),
                            }
                            .unwrap();
                            assert_eq!(actual.0.map(f64::to_bits), expected.0.map(f64::to_bits));
                            assert_eq!(actual.1, expected.1);
                            assert_eq!(*actual_reads.borrow(), *expected_reads.borrow());
                        }
                    }
                }
            }
        }
    }
}

#[cfg(test)]
mod actual_nutation_feed_rate_tests {
    use super::*;

    fn piecewise_nutation_image() -> LebFile {
        use crate::leb::{
            HEADER_SIZE, LEB1_BODY_ENTRY_SIZE, MAGIC_LEB1, SECTION_DIR_SIZE, SECTION_NUTATION,
        };
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&MAGIC_LEB1);
        for value in [1_u32, 2, 1] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        for value in [0.0_f64, 2.0, 0.0] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        bytes.extend_from_slice(&0_u32.to_le_bytes());
        bytes.extend_from_slice(&[0; 20]);
        let body_index = HEADER_SIZE + 2 * SECTION_DIR_SIZE;
        let body_data = body_index + LEB1_BODY_ENTRY_SIZE;
        let nutation_data = body_data + 24;
        for (kind, offset, size) in [
            (0, body_index, LEB1_BODY_ENTRY_SIZE),
            (SECTION_NUTATION, nutation_data, 104),
        ] {
            bytes.extend_from_slice(&kind.to_le_bytes());
            bytes.extend_from_slice(&0_u32.to_le_bytes());
            bytes.extend_from_slice(&(offset as u64).to_le_bytes());
            bytes.extend_from_slice(&(size as u64).to_le_bytes());
        }
        bytes.extend_from_slice(&0_i32.to_le_bytes());
        bytes.extend_from_slice(&0_u32.to_le_bytes());
        bytes.extend_from_slice(&1_u32.to_le_bytes());
        for value in [0.0_f64, 2.0, 2.0] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        for value in [0_u32, 3] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        bytes.extend_from_slice(&(body_data as u64).to_le_bytes());
        bytes.extend_from_slice(&[0; 24]);
        for value in [0.0_f64, 2.0, 1.0] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        for value in [1_u32, 2, 2, 0] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        // Independently fitted records have unequal slopes and a small value
        // mismatch at their join. A cross-record central difference is not
        // the derivative of either owned polynomial.
        for value in [
            0.0001_f64,
            0.00001,
            0.0,
            0.0,
            0.000130001,
            0.00002,
            0.0,
            0.0,
        ] {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        LebFile::from_bytes("synthetic-nutation-seam.leb".into(), bytes).unwrap()
    }

    #[test]
    fn stored_rate_uses_the_owned_polynomial_at_seams_and_interiors() {
        let mut context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        context.attach_leb_file(piecewise_nutation_image());
        with_sealed_serving_stack(&context, DeltaTModel::EspenakMeeus, false, false, |stack| {
            for (epoch, expected) in [
                (0.25, 0.00002_f64),
                (0.75, 0.00002),
                (1.0, 0.00004),
                (1.5, 0.00004),
            ] {
                let actual = stack.nutation_longitude_rate_deg_day(epoch)?;
                let expected = expected.to_degrees();
                assert!((actual - expected).abs() <= 16.0 * f64::EPSILON * expected.abs());
            }
            let old_cross_record = positions::nutation_rate_deg_day(1.0, |t| {
                (stack.nutation_at)(t).map(|(angle, _)| angle)
            })?;
            let owned = stack.nutation_longitude_rate_deg_day(1.0)?;
            assert!((old_cross_record - owned).abs() > owned.abs() * 0.1);
            assert!(matches!(
                stack.nutation_longitude_rate_deg_day(2.5),
                Err(Error::LebOutOfRange { body_id: -1, .. })
            ));
            Ok(())
        })
        .unwrap();
    }

    #[test]
    fn replacement_feed_differentiates_its_actual_cubic_not_the_stored_series() {
        let mut context = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        context.attach_leb_file(piecewise_nutation_image());
        with_sealed_serving_stack(&context, DeltaTModel::EspenakMeeus, false, false, |stack| {
            let replacement = |epoch: f64| {
                let t = epoch - 1.0;
                Ok((0.0001 + 0.00002 * t + 0.000001 * t * t * t, 0.0))
            };
            let actual = stack
                .with_nutation(&replacement)
                .nutation_longitude_rate_deg_day(1.0)?;
            let expected = 0.00002_f64.to_degrees();
            // A cubic has no fourth-order truncation error. This bound covers
            // binary64 coefficient evaluation and weighted sample subtraction.
            let roundoff = (512.0 * f64::EPSILON * 0.0002 / (1.0 / 128.0)).to_degrees();
            assert!((actual - expected).abs() <= roundoff);
            let unavailable = |_: f64| {
                Err(Error::NonFiniteArgument {
                    value: f64::INFINITY,
                })
            };
            assert!(matches!(
                stack
                    .with_nutation(&unavailable)
                    .nutation_longitude_rate_deg_day(1.0),
                Err(Error::NonFiniteArgument { .. })
            ));
            Ok(())
        })
        .unwrap();
    }
}
