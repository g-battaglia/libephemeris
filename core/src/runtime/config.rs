//! Thread-safe configuration and shared, lazily opened calculation resources.
//!
//! A runtime owns explicit settings. Each calculation retains an immutable
//! context snapshot, so subsequent configuration changes cannot invalidate
//! its files. File I/O and numerical errors use the ordinary core error type.
//!
//! Explicit overrides take priority over captured environment, dotenv, TOML
//! and defaults. [`Runtime::default`] does not discover those sources; use
//! [`Runtime::with_sources`] with [`ConfigurationSources::discover`] when that
//! is intended. An optional setter's `None` usually removes its explicit
//! override, allowing captured values to apply again; it need not mean `auto`,
//! `false` or a built-in default. Each setter documents exceptions.
//!
//! # Captured configuration example
//!
//! This example supplies one environment-style value directly, without reading
//! the process environment or any configuration file. No body file is needed.
//!
//! ```
//! use libephemeris_core::{config::{ConfigurationSources, Runtime}, errors::Error};
//! use std::collections::BTreeMap;
//!
//! fn main() -> Result<(), Error> {
//!     let sources = ConfigurationSources::load(
//!         BTreeMap::from([("LIBEPHEMERIS_MODE".into(), "leb".into())]),
//!         None,
//!         None,
//!     );
//!     let runtime = Runtime::with_sources(sources);
//!     runtime.set_calc_mode(Some("auto"))?; // Explicit override.
//!     runtime.set_calc_mode(None)?;         // Captured "leb" applies again.
//!     assert_eq!(runtime.effective_context()?.mode().as_calc_mode_name(), "leb");
//!     Ok(())
//! }
//! ```

use crate::context::{EphemerisContext, Observer, ProviderMode};
use crate::engine_api::{CurrentFileMetadata, MainImageLifecycle, Position, Session};
use crate::errors::{Error, Result, Warning};
use crate::leb::LebFile;
use std::collections::HashMap;
#[path = "config/data.rs"]
pub mod data;
#[path = "config/download_surface.rs"]
pub mod download_surface;
#[path = "config/hypothetical.rs"]
mod hypothetical;
#[path = "config/minor_bodies.rs"]
mod minor_bodies;
#[path = "config/moons.rs"]
mod moons;
#[path = "config/orbital.rs"]
mod orbital;
#[path = "config/reader_acquisition.rs"]
mod reader_acquisition;
#[path = "config/snapshots.rs"]
mod snapshots;
#[path = "config/sources.rs"]
mod sources;
#[path = "config/spk_auto.rs"]
pub mod spk_auto;
pub use sources::ConfigurationSources;
use std::path::PathBuf;
use std::sync::{Arc, Mutex, MutexGuard, OnceLock, Weak};

#[cfg(test)]
use reader_acquisition::partition_is_unused;

#[cfg(test)]
#[path = "config/observer_model_session_tests.rs"]
mod observer_model_session_tests;

/// Independently configured, thread-safe calculation runtime.
///
/// Numeric units and flags follow [`Session`]. Explicit LEB files are opened on
/// the first consuming calculation and their immutable images/caches are shared across
/// calls. The default context uses the core defaults; this constructor performs
/// no environment, TOML, network or filesystem access. Use [`Self::from_context`]
/// to supply settings resolved by the application.
#[derive(Debug, Default)]
pub struct Runtime {
    state: Mutex<RuntimeState>,
    /// Serialize LEB acquisition, never configuration edits or calculations.
    reader_acquisition: Mutex<()>,
    sources: ConfigurationSources,
    main_image: MainImageLifecycle,
}

#[derive(Debug, Default)]
struct RuntimeState {
    context: EphemerisContext,
    // Exactly one current snapshot; old generations survive only in active calls.
    cached_session: Option<CachedSession>,
    horizons_client: Option<Arc<crate::horizons::HorizonsClient>>,
    mode: Option<ProviderMode>,
    tier: Option<crate::context::Tier>,
    reader: Vec<(PathBuf, LebFile)>,
    reader_policy: Option<crate::leb::LebReadPolicy>,
    reader_generation: u64,
    reader_paths: Option<(u64, Vec<PathBuf>)>,
    // Only consumed admission identities are retained here. Full provenance
    // queries identify images through their own shared digest cells.
    reader_manifest_matches: HashMap<PathBuf, bool>,
    center_images: Option<(
        Option<PathBuf>,
        crate::context::Tier,
        crate::spk::SpkRegistry,
    )>,
    reloaded_sources: Option<ConfigurationSources>,
    spk_definitions: HashMap<i32, (PathBuf, i32)>,
    moons: crate::planetary_moons::MoonSpkRegistry,
    sbdb_elements: crate::minor_bodies::SbdbElementsCache,
    asteroid_names: crate::minor_bodies::AsteroidNameCache,
    auto_spk: HashMap<i32, Arc<Mutex<spk_auto::AutoSpkSettings>>>,
}

/// One configuration snapshot, promoted when a calculation actually needs data.
#[derive(Debug)]
struct CachedSession {
    session: Arc<Session>,
    resources: SessionResources,
}

/// Increasing acquisition scope; promotion retains existing image owners.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
enum SessionResources {
    Models,
    PlanetaryCore,
    // Adds the stored true-node channel; ordinary states include Moon/Earth
    // and the observer/deflector dependencies of custom Cartesian channels.
    LunarCore,
    All,
}

// Keep this emission in the configuration module so logger targets and
// module-path filters retain their existing meaning.
fn warn_optional_leb_open(path: &std::path::Path, error: &Error) {
    log::warn!("Failed to open optional LEB {}: {error}", path.display());
}

impl Runtime {
    /// Construct an independent context from resolved provider/time settings.
    ///
    /// Observer, sidereal state, angle results and explicit SPK registrations
    /// start fresh. Matching immutable LEB images share their existing storage.
    /// Kernel acquisition remains lazy.
    ///
    /// # Arguments
    ///
    /// * `base` - Runtime whose effective provider/time settings are inherited.
    /// * `ephe_path` - Kernel directory override; `None` keeps the inherited directory.
    /// * `ephe_file` - Kernel filename override. The default `"de440.bsp"` keeps
    ///   the inherited filename rather than forcing a DE440 selection.
    ///
    /// # Returns
    ///
    /// A shared handle to a runtime with independent observing state and any
    /// matching immutable LEB storage retained from `base`.
    ///
    /// # Errors
    ///
    /// Propagates errors from resolving the base runtime's effective configuration.
    pub fn independent_context(
        base: &Self,
        ephe_path: Option<PathBuf>,
        ephe_file: &str,
    ) -> Result<Arc<Self>> {
        let mut context = base.effective_context()?;
        context.clear_observer();
        context.set_sidereal_mode(0, 2451545.0, 0.0);
        context.clear_angles_cache();
        context.spk_registry_mut().clear();
        if ephe_path.is_some() {
            context.set_jpl_dir(ephe_path);
        }
        if ephe_file != "de440.bsp" {
            context.set_ephemeris_file(Some(ephe_file));
        }
        let state = base.lock();
        let sources = state.reloaded_sources.as_ref().unwrap_or(&base.sources);
        let selected = sources.leb_paths(&context);
        if context.leb_files().is_empty() {
            for (path, image) in &state.reader {
                if selected.contains(path) && state.reader_policy == Some(context.leb_read_policy())
                {
                    context.attach_leb_file(image.clone());
                }
            }
        }
        Ok(Self::shared(context))
    }

    /// Create an independently configured shared runtime and register a weak
    /// resource-lifecycle handle. The registry never extends a runtime's life.
    pub fn shared(context: EphemerisContext) -> Arc<Self> {
        let runtime = Arc::new(Self::from_context(context));
        let mut entries = shared_runtimes()
            .lock()
            .unwrap_or_else(|error| error.into_inner());
        entries.retain(|entry| entry.strong_count() > 0);
        entries.push(Arc::downgrade(&runtime));
        runtime
    }

    /// Release cached resources of still-live shared contexts. Their observing
    /// settings and active call-owned snapshots remain independent and valid.
    pub fn release_shared_resources() {
        let runtimes: Vec<_> = shared_runtimes()
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .iter()
            .filter_map(Weak::upgrade)
            .collect();
        for runtime in runtimes {
            runtime.release_resources();
        }
    }

    /// Acquire the selected native planetary SPK image through its shared cache.
    /// Explicit LEB mode refuses JPL access before any filesystem/network work.
    pub fn main_image(&self) -> Result<crate::positions::MainImageOwner> {
        self.main_image.acquire(&self.effective_context()?)
    }

    /// Independent Earth-orientation clock at UT1, returned as two-part TT JD.
    /// Model Delta-T overrides do not replace this modern frame clock.
    pub fn frame_tt_parts(&self, jd_ut1: f64) -> Result<(f64, f64)> {
        crate::smh2016::timescale_parts_from_ut1(jd_ut1)
    }

    /// Split a finite TT Julian day without converting its physical timescale.
    pub fn tt_parts(jd_tt: f64) -> Result<(f64, f64)> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let whole = jd_tt.floor();
        Ok((whole, jd_tt - whole))
    }

    /// Open a named planetary SPK using this runtime's data roots and policy,
    /// without replacing its current main-image selection.
    pub fn open_planetary_image(&self, filename: &str) -> Result<crate::positions::MainImageOwner> {
        let mut context = self.effective_context()?;
        context.set_ephemeris_file(Some(filename));
        crate::positions::acquire_current_main_image(&context)
    }

    /// Model Delta-T in days at UT1, using this runtime's resolved settings.
    pub fn delta_t(&self, jd_ut1: f64) -> Result<f64> {
        self.model_session()?.deltat(jd_ut1)
    }

    /// Set independent sidereal projection bits without changing mode/epoch.
    pub fn set_sidereal_bits(&self, bits: i32) {
        self.edit().context.set_sidereal_bits(bits);
    }

    /// Own a context without copying its kernel data or performing I/O.
    pub fn from_context(context: EphemerisContext) -> Self {
        Self {
            state: Mutex::new(RuntimeState {
                mode: Some(context.mode()),
                tier: Some(context.tier()),
                context: context.clone(),
                cached_session: None,
                horizons_client: context.shared_horizons_client(),
                reader: Vec::new(),
                reader_policy: None,
                reader_generation: 0,
                reader_paths: None,
                reader_manifest_matches: HashMap::new(),
                center_images: None,
                reloaded_sources: None,
                moons: crate::planetary_moons::MoonSpkRegistry::new(),
                sbdb_elements: crate::minor_bodies::SbdbElementsCache::default(),
                asteroid_names: crate::minor_bodies::AsteroidNameCache::default(),
                auto_spk: HashMap::new(),
                spk_definitions: context
                    .spk_registry()
                    .registered_bodies()
                    .into_iter()
                    .map(|record| (record.body, (PathBuf::from(record.kernel), record.naif_id)))
                    .collect(),
            }),
            reader_acquisition: Mutex::new(()),
            sources: ConfigurationSources::default(),
            main_image: MainImageLifecycle::new(),
        }
    }

    /// Retain configuration and any opened images independently of future edits.
    /// No filesystem access or numerical evaluation occurs.
    pub fn context(&self) -> EphemerisContext {
        let state = self.lock();
        let mut context = state.context.clone();
        for (_, reader) in &state.reader {
            context.attach_leb_file(reader.clone());
        }
        context
    }

    /// Construct a runtime whose unset fields resolve from acquired settings.
    /// No kernel is opened and invalid enum settings are checked when consumed.
    pub fn with_sources(sources: ConfigurationSources) -> Self {
        Self {
            state: Mutex::new(RuntimeState::default()),
            reader_acquisition: Mutex::new(()),
            sources,
            main_image: MainImageLifecycle::new(),
        }
    }

    /// Borrow the original captured input sources and their warnings.
    ///
    /// This is the capture supplied at construction, not a view of effective
    /// settings. [`Self::reload_configuration`] and [`Self::load_dotenv`] may
    /// install a newer internal capture without replacing this returned value.
    /// Use [`Self::effective_context`] to inspect currently resolved settings;
    /// warnings in a later internal capture are not exposed by this accessor.
    pub fn sources(&self) -> &ConfigurationSources {
        &self.sources
    }

    /// Resolve effective configuration without opening a file or calculating.
    /// Invalid configured mode/network-policy values remain ordinary errors.
    pub fn effective_context(&self) -> Result<EphemerisContext> {
        let state = self.lock();
        self.resolve(&state)
    }

    /// Set or remove the explicit provider-mode override, without opening data.
    ///
    /// Names are `auto`, `leb`, `skyfield` (native JPL/SPK) and `horizons`,
    /// normalized by Unicode lowercase and surrounding-whitespace trimming.
    /// `Some("auto")` explicitly selects automatic routing. `None` instead
    /// removes the override: captured environment/dotenv/TOML may select a
    /// different effective mode, with `auto` used only when none does.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidCalcMode`] for an unknown name without changing
    /// configuration. Resource/policy/coverage checks occur when consumed.
    pub fn set_calc_mode(&self, mode: Option<&str>) -> Result<()> {
        let selected =
            ProviderMode::from_calc_mode_name(mode.unwrap_or("auto")).ok_or_else(|| {
                Error::InvalidCalcMode {
                    value: mode.unwrap_or("auto").to_owned(),
                    source: crate::context::CALC_MODE_SOURCE_SETTER,
                }
            })?;
        let explicit = mode.map(|_| selected);
        let mut state = self.lock();
        if state.context.mode() == selected && state.mode == explicit {
            return Ok(());
        }
        state.cached_session = None;
        state.context.set_mode_name(mode.unwrap_or("auto"))?;
        state.mode = explicit;
        Ok(())
    }

    /// Set the explicit data tier: exactly `base`, `medium` or `extended`.
    ///
    /// A successful call clears the context's explicit JPL filename, even
    /// when the tier is unchanged. Captured filename configuration may then
    /// apply; otherwise the tier chooses its default kernel. An explicit LEB
    /// file remains selected. A changed tier retires runtime discovery/main
    /// selections; active snapshots retain their old images. No file is opened.
    /// A tier does not extend the actual coverage of a body or model.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidTier`] for an unknown spelling, preserving the
    /// selected tier and filename. The cached calculation snapshot may still
    /// be retired by this call.
    pub fn set_precision_tier(&self, tier: &str) -> Result<()> {
        let mut state = self.edit();
        state.context.set_tier_name(tier)?;
        if state.tier == Some(state.context.tier()) {
            return Ok(());
        }
        state.tier = Some(state.context.tier());
        state.center_images = None;
        state.context = self.main_image.invalidate(state.context.clone());
        state.invalidate_readers();
        Ok(())
    }

    /// Select an explicit local LEB path, or remove that override with `None`.
    ///
    /// `None` resumes captured file selection and then ordinary discovery; it
    /// does not disable LEB use. Opening and validation are lazy. An unchanged
    /// path keeps current readers; a changed selection retires them for future
    /// queries. Previously acquired snapshots keep their old immutable images.
    pub fn set_leb_file(&self, path: Option<PathBuf>) {
        let mut state = self.lock();
        if !state.context.set_leb_file(path) {
            return;
        }
        state.invalidate_readers();
    }

    /// Override observed-IERS Delta-T use; `None` resumes captured/default policy.
    /// Observational rows are acquired only when their ordinary consumer needs
    /// them. This setting does not replace independent bundled frame clocks.
    pub fn set_iers_delta_t_enabled(&self, enabled: Option<bool>) {
        self.edit().context.set_iers_delta_t(enabled);
    }

    /// Set model TT−UT1 in days, not seconds.
    ///
    /// `None` and `Some(DELTAT_AUTOMATIC)` resume computed Delta-T. Other
    /// values are stored without finite/range validation; consuming operations
    /// retain their own domains. Independent frame/observer clocks are not
    /// replaced by this override.
    pub fn set_delta_t_userdef(&self, days: Option<f64>) {
        self.edit().context.set_user_delta_t(days);
    }

    /// Select a named model, or remove its explicit override with `None`.
    ///
    /// Names follow [`crate::time::DeltaTModel::from_name`]. `None` resumes
    /// captured environment/dotenv/TOML selection before the core default.
    /// Model Delta-T values are days; a user value override has its own role.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidDeltaTModel`] for an unknown name, preserving
    /// the selected model. The cached snapshot may nevertheless be retired.
    pub fn set_delta_t_model(&self, model: Option<&str>) -> Result<()> {
        self.edit().context.set_delta_t_model_name(model)
    }

    /// Select a `SIDM_*` mode and supported `SIDBIT_*` projection/epoch flags.
    ///
    /// USER mode reads `t0` as a TT Julian day, or UT1 when `SIDBIT_USER_UT`
    /// is retained, and `ayan_t0` as degrees. Values are stored without a
    /// finite-date check here. Unsupported projection bits produce returned
    /// advisories; the selection itself is stored. Identical stored IEEE words
    /// keep the current snapshot. Return warnings are produced on every call;
    /// callers decide how to present them.
    pub fn set_sid_mode(&self, mode: i32, t0: f64, ayan_t0: f64) -> Vec<Warning> {
        let mut state = self.lock();
        let before = state.context.sidereal();
        let warnings = state.context.set_sidereal_mode(mode, t0, ayan_t0);
        let unchanged = match (before, state.context.sidereal()) {
            (
                crate::context::SiderealConfig::Sidereal(before),
                crate::context::SiderealConfig::Sidereal(after),
            ) => {
                let (old_mode, old_t0, old_ayan) = before.full();
                let (new_mode, new_t0, new_ayan) = after.full();
                old_mode == new_mode
                    && before.bits() == after.bits()
                    && old_t0.to_bits() == new_t0.to_bits()
                    && old_ayan.to_bits() == new_ayan.to_bits()
            }
            // A never-set tropical configuration is not an explicit mode-zero
            // selection, even when their public scalar readbacks coincide.
            _ => false,
        };
        if !unchanged {
            state.cached_session = None;
        }
        warnings
    }

    /// Select a geographic observer in degrees east/north and metres.
    /// Non-finite or out-of-range coordinates are rejected without changing state.
    /// Bitwise-identical observer values retain the current snapshot.
    pub fn set_topo(&self, lon: f64, lat: f64, alt: f64) -> Result<()> {
        let observer = Observer::new(lon, lat, alt)?;
        let mut state = self.lock();
        if state.context.observer().is_some_and(|current| {
            current.lon_deg().to_bits() == observer.lon_deg().to_bits()
                && current.lat_deg().to_bits() == observer.lat_deg().to_bits()
                && current.height_m().to_bits() == observer.height_m().to_bits()
        }) {
            return Ok(());
        }
        state.cached_session = None;
        state.context.set_observer(observer);
        Ok(())
    }

    /// Set the explicit JPL search directory, or clear it with `None`.
    ///
    /// An unchanged path is a no-op. A changed path retires runtime LEB
    /// discovery, optional centers and the current main-image association,
    /// while existing snapshots keep their images and unrelated explicit SPK
    /// registrations remain. `None` resumes other resolved data-root choices.
    ///
    /// Unlike a store-only context setter, `Some(path)` also **eagerly scans
    /// local SPKs** and may open/register images. Discovery errors are discarded
    /// by this infallible setter; successful registrations are retained. Later
    /// consumers still report missing or invalid required data. No download is
    /// requested by this setter.
    pub fn set_ephe_path(&self, path: Option<PathBuf>) {
        let mut state = self.lock();
        if state.context.paths().jpl_dir == path {
            return;
        }
        state.invalidate_discovery();
        let registry = state.context.spk_registry().clone();
        state.context.set_jpl_dir(path.clone());
        *state.context.spk_registry_mut() = registry;
        state.context = self.main_image.invalidate(state.context.clone());
        if let Some(path) = path {
            let _ = crate::spk::discover_local_spks(state.context.spk_registry_mut(), &path);
            for record in state.context.spk_registry().registered_bodies() {
                state
                    .spk_definitions
                    .insert(record.body, (PathBuf::from(record.kernel), record.naif_id));
            }
        }
    }

    /// Select future LEB storage without reopening already acquired snapshots.
    /// `None` restores captured configuration. Invalid settings preserve state;
    /// changing policy retires discovery/image caches, not explicit attachments.
    pub fn set_leb_read_policy(&self, policy: Option<crate::leb::LebReadPolicy>) -> Result<()> {
        let mut state = self.edit();
        if state.context.set_leb_read_policy(policy)? {
            state.invalidate_readers();
        }
        Ok(())
    }

    /// Select an explicit JPL filename. Opening occurs only when it is consumed.
    pub fn set_ephemeris_file(&self, filename: &str) {
        let mut state = self.edit();
        if state.context.ephemeris_file_override() == Some(filename) {
            return;
        }
        state.context.set_ephemeris_file(Some(filename));
        if self
            .resolve(&state)
            .is_ok_and(|resolved| self.main_image.has_current(&resolved))
        {
            return;
        }
        state.context = self.main_image.invalidate(state.context.clone());
    }

    /// Set lunar tidal acceleration in arcseconds per Julian century squared.
    /// The automatic sentinel restores the model's ordinary acceleration policy.
    pub fn set_tid_acc(&self, acceleration: f64) {
        self.edit()
            .context
            .set_tidal_acceleration(Some(acceleration));
    }

    /// Set atmospheric lapse rate in K/m; `None` restores the physical default.
    pub fn set_lapse_rate(&self, lapse_rate: Option<f64>) {
        self.edit().context.set_lapse_rate(lapse_rate);
    }

    /// Set automatic SPK provisioning; `None` restores acquired configuration.
    pub fn set_auto_spk_download(&self, enabled: Option<bool>) {
        self.edit().context.set_auto_spk_download(enabled);
    }

    /// Set an explicit SPK cache directory, without opening or creating files.
    /// `None` removes the override and resumes captured/default directory
    /// resolution. Existing explicit registrations are not reinterpreted.
    pub fn set_spk_cache_dir(&self, path: Option<PathBuf>) {
        self.edit().context.set_spk_cache_dir(path);
    }

    /// Set non-negative padding in days around requested SPK download coverage.
    /// Negative padding is rejected without changing the previous value.
    pub fn set_spk_date_padding(&self, days: i64) -> Result<()> {
        let days = u64::try_from(days).map_err(|_| Error::InvalidPadding {
            entry: "set_spk_date_padding",
            value: days as f64,
        })?;
        self.edit().context.set_spk_date_padding_days(days);
        Ok(())
    }

    /// Set strict minor-body precision; `None` restores acquired configuration.
    pub fn set_strict_precision(&self, enabled: Option<bool>) {
        self.edit().context.set_strict_precision(enabled);
    }

    /// Clear observer, sidereal selection, angle results and the runtime's
    /// retained Horizons client, while retaining data images and other settings.
    ///
    /// Active callers retain their snapshots and client owners. With no such
    /// transient state, the current snapshot is preserved; an actual reset
    /// retires it for future callers. This is not cancellation or full resource
    /// teardown; see [`Self::release_resources`] and [`Self::close`].
    pub fn reset_session(&self) {
        let mut state = self.lock();
        if state.context.observer().is_none()
            && matches!(
                state.context.sidereal(),
                crate::context::SiderealConfig::Tropical
            )
            && state.context.angles_cache().is_empty()
            && state.horizons_client.is_none()
        {
            return;
        }
        state.cached_session = None;
        state.horizons_client = None;
        state.context.reset_session();
    }

    /// Reset transient settings and retire runtime-owned data selections.
    ///
    /// Applies [`EphemerisContext::close`], clears explicit runtime SPK body
    /// definitions and registered moon kernels, and drops the current main,
    /// reader and Horizons owners. The explicit tier override is removed, so
    /// captured tier settings may apply again. Provider mode, explicit LEB
    /// selection, user Delta-T, model and strictness persist as documented by
    /// the context lifecycle. Existing call-owned snapshots remain usable;
    /// this does not interrupt their calculations.
    pub fn close(&self) {
        let mut state = self.edit();
        state.center_images = None;
        state.horizons_client = None;
        state.context.close();
        state.context = self.main_image.close(state.context.clone());
        state.spk_definitions.clear();
        state.moons.close_moon_kernels();
        state.tier = None;
        state.invalidate_readers();
    }

    /// Retire acquired resources without resetting observing/configuration settings.
    ///
    /// Drops runtime LEB readers, optional centers, the current main-image
    /// association, context SPK resources and acquired IERS rows. Explicit SPK
    /// definitions remain for lazy reopening. It also releases the **process-
    /// default ASSIST session**, not just this runtime's readers. The separate
    /// runtime Horizons client and moon registry are not cleared here.
    ///
    /// Live computations retain their own immutable owners; releasing a
    /// runtime reference neither guarantees immediate memory reclamation nor
    /// cancels an active calculation. Future consuming operations reacquire
    /// released data and can fail under their ordinary I/O/policy checks.
    pub fn release_resources(&self) {
        crate::nbody::release_default_assist_session();
        let mut state = self.edit();
        state.invalidate_discovery();
        let selection = state.context.leb_file().cloned();
        state.context.set_leb_file(Some(PathBuf::new()));
        state.context.set_leb_file(selection);
        state.context = self.main_image.close(state.context.clone());
        state.context.spk_registry_mut().clear();
        state.context.clear_iers_cache();
    }

    /// Capture current environment/dotenv/TOML settings for future resolution.
    ///
    /// Explicit overrides remain authoritative. Discovery reads configuration
    /// files but opens no ephemeris; unreadable files are ignored and malformed
    /// TOML is retained as a warning in the internal capture. Reader discovery,
    /// optional-center state and the cached session are retired. No process
    /// environment variable is mutated, and active snapshots remain valid.
    /// [`Self::sources`] still returns the original construction-time capture.
    pub fn reload_configuration(&self) {
        let sources = ConfigurationSources::discover();
        let mut state = self.edit();
        state.center_images = None;
        state.reloaded_sources = Some(sources);
        state.invalidate_readers();
    }

    /// Open and retain the active native LEB reader only for automatic or LEB mode.
    /// LEB mode propagates file errors; auto returns `None` for an unavailable or
    /// invalid optional container. Resource-policy failures always propagate.
    pub fn leb_reader(&self) -> Result<Option<LebFile>> {
        let mut state = self.lock();
        let mut context = self.resolve(&state)?;
        self.attach_readers(&mut state, &mut context)?;
        Ok(context.leb_files().first().cloned())
    }

    /// Inspect selected Sun/Moon file metadata, acquiring configured LEB readers.
    ///
    /// # Arguments
    ///
    /// `index` is 0 for the Sun or 1 for the Moon. Other values return the
    /// empty record without configuration or resource acquisition.
    ///
    /// # Returns
    ///
    /// `(path, start_jd, end_jd, de_number)`, or `("", 0, 0, 0)` when there is
    /// no applicable current image. Selection follows
    /// [`Session::current_file_metadata`]:
    ///
    /// - Stored LEB bounds are the selected body's **TT Julian days**. Source
    ///   number is 441 when coverage exceeds the DE440 span, otherwise zero
    ///   (unknown); this heuristic does not authenticate the coefficients.
    /// - An already retained main SPK supplies its complete **TDB Julian-day
    ///   descriptor envelope**, not guaranteed gap-free body/center coverage.
    ///   Its DE number is parsed from the actual effective opened name when
    ///   possible, otherwise zero; it is not inferred from coefficient data.
    ///
    /// This query may open/hash configured LEB images, unlike the Session
    /// metadata-only method. It does not open a not-yet-acquired main SPK.
    ///
    /// # Errors
    ///
    /// Returns configuration, required LEB acquisition/resource-limit, closed
    /// reader or invalid metadata errors. Auto mode retains its ordinary
    /// optional-reader policy. No body coefficients or timescale are evaluated.
    pub fn current_file_data(&self, index: i32) -> Result<(String, f64, f64, i32)> {
        if !(0..=1).contains(&index) {
            return Ok((String::new(), 0.0, 0.0, 0));
        }
        let mut state = self.lock();
        let mut context = self.resolve(&state)?;
        self.attach_readers(&mut state, &mut context)?;
        let owner = self.main_image.current_owner();
        let main = owner
            .as_ref()
            .filter(|_| self.main_image.has_current(&context));
        if let Some(owner) = main {
            context = owner.with_current_context(context)?.context().clone();
        }
        let session = Session::from_context(context);
        match session.current_file_metadata(index, main.map(|owner| owner.registry_key()))? {
            Some(CurrentFileMetadata::Main {
                path,
                start_jd,
                end_jd,
            }) => {
                let name = main
                    .map(|owner| owner.effective_name())
                    .unwrap_or("")
                    .to_ascii_lowercase();
                let de_number = name
                    .strip_prefix("de")
                    .and_then(|name| name.split('.').next())
                    .and_then(|digits| digits.parse::<i32>().ok())
                    .unwrap_or(0);
                Ok((path, start_jd, end_jd, de_number))
            }
            Some(CurrentFileMetadata::Stored {
                path,
                start_jd,
                end_jd,
                source_de_number,
                ..
            }) => Ok((path, start_jd, end_jd, source_de_number)),
            None => Ok((String::new(), 0.0, 0.0, 0)),
        }
    }

    /// Resolve the absolute provider data directory without creating it.
    pub fn library_path(&self) -> Result<PathBuf> {
        let context = self.effective_context()?;
        let path = context
            .paths()
            .jpl_dir
            .as_ref()
            .or(context.jpl_data_dir())
            .cloned()
            .unwrap_or_else(|| PathBuf::from("."));
        std::path::absolute(&path).map_err(|error| Error::EphemerisLoaderIo {
            path: path.display().to_string(),
            kind: error.kind(),
            errno: error.raw_os_error(),
            message: error.to_string(),
        })
    }

    /// Advise the operating system that library data pages may be reclaimed.
    /// This is best effort and never invalidates immutable image owners.
    pub fn release_data_cache(&self) -> Result<()> {
        for (_, reader) in &self.lock().reader {
            reader.release_file_pages();
        }
        crate::leb::release_data_cache(&self.library_path()?);
        Ok(())
    }

    /// Replace the context-local angle cache with owned scalar values.
    pub fn set_angles_cache(&self, values: &HashMap<String, f64>) {
        self.edit().context.set_angles_cache(values);
    }

    /// Clear the context-local angle cache without changing observing settings.
    pub fn clear_angles_cache(&self) {
        self.edit().context.clear_angles_cache();
    }

    /// Register an immutable SPK image and a library-body to NAIF target mapping.
    /// Missing files, corrupt directories and absent targets propagate unchanged.
    pub fn register_spk_body(&self, body: i32, path: &std::path::Path, target: i32) -> Result<()> {
        let path = if path.is_absolute() || path.exists() {
            path.to_path_buf()
        } else {
            self.library_path()?.join(path)
        };
        let mut state = self.edit();
        let registry = state.context.spk_registry_mut();
        let label = registry.register_file(&path)?;
        registry.register_body(body, &label, target)?;
        state.spk_definitions.insert(body, (path, target));
        Ok(())
    }

    /// Remove a context-local body registration, preserving images in snapshots.
    pub fn unregister_spk_body(&self, body: i32) {
        let mut state = self.edit();
        state.context.spk_registry_mut().unregister_body(body);
        state.spk_definitions.remove(&body);
    }

    /// Retained explicit SPK registrations, including images released for lazy
    /// reopening. Returned paths and target IDs are owned scalar configuration.
    pub fn registered_spks(&self) -> HashMap<i32, (String, i32)> {
        self.lock()
            .spk_definitions
            .iter()
            .map(|(body, (path, target))| (*body, (path.to_string_lossy().into_owned(), *target)))
            .collect()
    }

    fn resolve(&self, state: &RuntimeState) -> Result<EphemerisContext> {
        state
            .reloaded_sources
            .as_ref()
            .unwrap_or(&self.sources)
            .resolve(&state.context, state.mode, state.tier)
    }

    /// Calculate a position at UT1 Julian day, with the units/flags of [`Session::calc_ut`].
    /// Configuration and file failures propagate before numerical evaluation.
    pub fn calc_ut(&self, jd_ut: f64, body: i32, flags: i32) -> Result<Position> {
        if let Some(position) = self.registered_moon_position(jd_ut, body, flags, true)? {
            return Ok(position);
        }
        self.position_with_native_minor(jd_ut, body, flags, true)
    }

    /// Calculate a position at TT Julian day, with the units/flags of [`Session::calc`].
    /// Configuration and file failures propagate before numerical evaluation.
    pub fn calc(&self, jd_tt: f64, body: i32, flags: i32) -> Result<Position> {
        if let Some(position) = self.registered_moon_position(jd_tt, body, flags, false)? {
            return Ok(position);
        }
        self.position_with_native_minor(jd_tt, body, flags, false)
    }

    /// Physical phenomena at UT1 Julian day, with the twenty-component layout
    /// and units of [`Session::pheno_ut`]. Native calculations share this
    /// runtime's immutable main image; configuration and data errors propagate.
    pub fn pheno_ut(
        &self,
        jd_ut: f64,
        body: i32,
        flags: i32,
    ) -> Result<crate::engine_api::Phenomena> {
        self.shared_session()?
            .pheno_ut_with_main_image(jd_ut, body, flags, &self.main_image)
    }

    /// Physical phenomena at TT Julian day. Angles, magnitudes and fractions
    /// follow [`Session::pheno`]; native support images remain owned throughout
    /// the calculation and are reused until their selection is invalidated.
    pub fn pheno(&self, jd_tt: f64, body: i32, flags: i32) -> Result<crate::engine_api::Phenomena> {
        self.shared_session()?
            .pheno_with_main_image(jd_tt, body, flags, &self.main_image)
    }

    fn position_with_native_minor(
        &self,
        epoch: f64,
        body: i32,
        flags: i32,
        is_ut: bool,
    ) -> Result<Position> {
        let mut supplier = crate::engine_api::NativeMinorAssist::for_runtime(self);
        let resources = if flags & crate::constants::FLG_SIDEREAL != 0 {
            SessionResources::All
        } else if body == crate::constants::TRUE_NODE || body == -crate::constants::TRUE_NODE {
            SessionResources::LunarCore
        } else if crate::leb::is_stored_state_body(body) {
            SessionResources::PlanetaryCore
        } else {
            SessionResources::All
        };
        // ECL_NUT uses independent clocks and frame formulas. It consumes no
        // planetary image, even when the surrounding runtime selects LEB/JPL.
        let session = if body == crate::constants::ECL_NUT {
            self.model_session()?
        } else {
            self.shared_session_for(resources)?
        };
        let mut position = session.position_with_native_minor_and_main_image(
            &supplier,
            &self.main_image,
            epoch,
            body,
            flags,
            is_ut,
        )?;
        let diagnostics = supplier.drain_diagnostics_exclusive();
        position.warnings.extend(
            diagnostics
                .warnings
                .into_iter()
                .map(crate::engine_api::Advisory::from_engine),
        );
        position.messages.extend(diagnostics.messages);
        Ok(position)
    }

    /// Planet-centric position at UT1 Julian day. The configured physical
    /// Delta-T model converts UT1 to TT before the shared planet-centric solver.
    pub fn calc_pctr_ut(&self, jd_ut: f64, body: i32, center: i32, flags: i32) -> Result<Position> {
        let session = self.shared_session()?;
        let jd_tt = jd_ut + session.deltat(jd_ut)?;
        session.calc_pctr(jd_tt, body, center, flags)
    }

    /// House cusps and principal angles at UT1, in degrees.
    /// Latitude/longitude are degrees north/east; the byte selects a house system.
    pub fn houses(
        &self,
        jd_ut: f64,
        lat: f64,
        lon: f64,
        system: u8,
    ) -> Result<crate::engine_api::Houses> {
        self.house_session(char::from(system), 0)?
            .houses(jd_ut, lat, lon, system)
    }

    /// Calculate and cache eight geometric angles and four planetary longitudes.
    ///
    /// `jd_ut` is UT1 Julian day, latitude/longitude are degrees (east positive).
    /// Selects a zero-height observer; Placidus angles use the Porphyry fallback
    /// at polar latitudes. Planetary values are Sun, Moon, Mercury and Venus.
    /// Coordinate, time, geometry and data failures propagate unchanged.
    pub fn calc_angles(&self, jd_ut: f64, lat: f64, lon: f64) -> Result<HashMap<String, f64>> {
        self.edit()
            .context
            .set_observer(Observer::new(lon, lat, 0.0)?);
        let values = self.shared_session()?.public_calc_angles(jd_ut, lat, lon)?;
        let result: HashMap<_, _> = crate::utils::CALC_ANGLES_KEYS
            .iter()
            .zip(values)
            .map(|(name, value)| ((*name).to_string(), value))
            .collect();
        self.edit().context.set_angles_cache(&result);
        Ok(result)
    }

    /// Convert true/apparent altitude in degrees using pressure in hPa and °C.
    /// `flag == 0` selects true to apparent; other values select the reverse.
    /// Uses the closed-form atmospheric model, including its horizon limits.
    pub fn refrac(&self, alt: f64, atpress: f64, attemp: f64, flag: i32) -> f64 {
        if flag == 0 {
            crate::refraction::refrac_true_to_app(alt, atpress, attemp)
        } else {
            crate::refraction::refrac_app_to_true(alt, atpress, attemp)
        }
    }

    /// Extended atmospheric altitude conversion, with the four physical details.
    ///
    /// Altitude is degrees, observer height metres, pressure hPa, temperature
    /// °C and lapse rate K/m. `None` uses the context's lapse rate. The output
    /// is `(converted, (true, apparent, refraction, horizon_dip))`, in degrees.
    /// A zero direction flag selects true to apparent, other flags the reverse.
    #[allow(
        clippy::too_many_arguments,
        reason = "Independent public atmospheric quantities"
    )]
    pub fn refrac_extended(
        &self,
        alt: f64,
        geoalt: f64,
        atpress: f64,
        attemp: f64,
        lapserate: Option<f64>,
        flag: i32,
    ) -> (f64, (f64, f64, f64, f64)) {
        let lapse = lapserate.unwrap_or_else(|| self.lock().context.effective_lapse_rate());
        if flag == 0 {
            crate::refraction::refrac_extended_true_to_app(alt, geoalt, atpress, attemp, lapse)
        } else {
            crate::refraction::refrac_extended_app_to_true(alt, geoalt, atpress, attemp, lapse)
        }
    }

    /// Retire the cached snapshot before changing any configuration or resource.
    /// Existing callers retain their own Arc; no historical snapshot is stored.
    fn edit(&self) -> MutexGuard<'_, RuntimeState> {
        let mut state = self.lock();
        state.cached_session = None;
        state
    }

    /// Recover the guarded state after poisoning. Partial edits are not rolled
    /// back; this preserves the existing policy rather than promising a transaction.
    fn lock(&self) -> MutexGuard<'_, RuntimeState> {
        self.state.lock().unwrap_or_else(|error| error.into_inner())
    }
}

fn shared_runtimes() -> &'static Mutex<Vec<Weak<Runtime>>> {
    static RUNTIMES: OnceLock<Mutex<Vec<Weak<Runtime>>>> = OnceLock::new();
    RUNTIMES.get_or_init(|| Mutex::new(Vec::new()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn runtime_position_entry_retains_native_minor_ownership_but_public_acquisitions_are_opaque() {
        let runtime = Runtime::default();
        let session = runtime.model_session().unwrap();
        let supplier = crate::engine_api::NativeMinorAssist::for_runtime(&runtime);
        let acquisitions = crate::engine_api::PositionAcquisitions {
            minor_assist: Some(&supplier),
            ..Default::default()
        };
        // ECL_NUT exercises both real Runtime entries without opening orbit
        // data. This probes ownership wiring, not emission geometry or accuracy.
        let ((), observations) = Session::observe_minor_entry_ownership(|| {
            runtime
                .calc(2451545.0, crate::constants::ECL_NUT, 0)
                .unwrap();
            runtime
                .calc_ut(2451545.0, crate::constants::ECL_NUT, 0)
                .unwrap();
            session
                .calc_with_acquisitions_and_main_image(
                    2451545.0,
                    crate::constants::ECL_NUT,
                    0,
                    &acquisitions,
                    &runtime.main_image,
                )
                .unwrap();
            session
                .calc_ut_with_acquisitions_and_main_image(
                    2451545.0,
                    crate::constants::ECL_NUT,
                    0,
                    &acquisitions,
                    &runtime.main_image,
                )
                .unwrap();
        });
        assert_eq!(observations, [true, true, false, false]);
    }

    /// One constant ICRS channel, optionally with an unconsumed auxiliary byte.
    fn partition_fixture(body: i32, auxiliary: Option<u32>) -> Vec<u8> {
        use crate::leb::{HEADER_SIZE, LEB1_BODY_ENTRY_SIZE, SECTION_DIR_SIZE};
        let sections = if auxiliary.is_some() { 2_u32 } else { 1 };
        let body_at = HEADER_SIZE + sections as usize * SECTION_DIR_SIZE;
        let data_at = body_at + LEB1_BODY_ENTRY_SIZE;
        let mut bytes = vec![0; data_at + 3 * 8 + usize::from(auxiliary.is_some())];
        bytes[..4].copy_from_slice(&crate::leb::MAGIC_LEB1);
        bytes[4..8].copy_from_slice(&1_u32.to_le_bytes());
        bytes[8..12].copy_from_slice(&sections.to_le_bytes());
        bytes[12..16].copy_from_slice(&1_u32.to_le_bytes());
        bytes[16..24].copy_from_slice(&2451545.0_f64.to_le_bytes());
        bytes[24..32].copy_from_slice(&2451546.0_f64.to_le_bytes());
        bytes[72..80].copy_from_slice(&(body_at as u64).to_le_bytes());
        bytes[80..88].copy_from_slice(&(LEB1_BODY_ENTRY_SIZE as u64).to_le_bytes());
        if let Some(id) = auxiliary {
            let payload = bytes.len() as u64 - 1;
            bytes[88..92].copy_from_slice(&id.to_le_bytes());
            bytes[96..104].copy_from_slice(&payload.to_le_bytes());
            bytes[104..112].copy_from_slice(&1_u64.to_le_bytes());
        }
        bytes[body_at..body_at + 4].copy_from_slice(&body.to_le_bytes());
        bytes[body_at + 8..body_at + 12].copy_from_slice(&1_u32.to_le_bytes());
        bytes[body_at + 12..body_at + 20].copy_from_slice(&2451545.0_f64.to_le_bytes());
        bytes[body_at + 20..body_at + 28].copy_from_slice(&2451546.0_f64.to_le_bytes());
        bytes[body_at + 28..body_at + 36].copy_from_slice(&1.0_f64.to_le_bytes());
        bytes[body_at + 40..body_at + 44].copy_from_slice(&3_u32.to_le_bytes());
        bytes[body_at + 44..body_at + 52].copy_from_slice(&(data_at as u64).to_le_bytes());
        bytes
    }

    #[test]
    fn planetary_acquisition_inspects_contents_and_promotes_unused_partitions() {
        let directory = discovery_directory("planetary-partitions");
        let core = directory.join("base_core.leb2");
        let mut unreviewed = std::fs::read(crate::test_data::core_leb_path()).unwrap();
        // Same valid coefficients, a different complete-file identity: the
        // optimization must work without a manifest/name-content assumption.
        unreviewed.push(0);
        std::fs::write(&core, unreviewed).unwrap();
        let companion = directory.join("base_exotics.leb2");
        for body in (SUN..=crate::constants::PLUTO).chain(std::iter::once(crate::constants::EARTH))
        {
            std::fs::write(&companion, partition_fixture(body, None)).unwrap();
            assert!(!partition_is_unused(
                &companion,
                None,
                SessionResources::PlanetaryCore
            ));
        }
        std::fs::write(
            &companion,
            partition_fixture(crate::constants::TRUE_NODE, None),
        )
        .unwrap();
        assert!(partition_is_unused(
            &companion,
            None,
            SessionResources::PlanetaryCore
        ));
        assert!(!partition_is_unused(
            &companion,
            None,
            SessionResources::LunarCore
        ));
        for auxiliary in [
            crate::leb::SECTION_NUTATION,
            crate::leb::SECTION_DELTA_T,
            crate::leb::SECTION_STARS,
            99,
        ] {
            std::fs::write(&companion, partition_fixture(13, Some(auxiliary))).unwrap();
            assert!(!partition_is_unused(
                &companion,
                None,
                SessionResources::PlanetaryCore
            ));
        }
        std::fs::write(&companion, b"corrupt").unwrap();
        assert!(!partition_is_unused(
            &companion,
            None,
            SessionResources::PlanetaryCore
        ));
        std::fs::write(&companion, partition_fixture(13, None)).unwrap();
        assert!(partition_is_unused(
            &companion,
            None,
            SessionResources::PlanetaryCore
        ));
        let runtime = discovery_runtime(&directory);
        let planetary = runtime
            .shared_session_for(SessionResources::PlanetaryCore)
            .unwrap();
        assert_eq!(planetary.context().leb_files().len(), 1);
        assert_eq!(runtime.lock().reader.len(), 1);
        let original = planetary
            .calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
            .unwrap();
        let complete = runtime.shared_session().unwrap();
        assert!(
            runtime
                .lock()
                .reader
                .iter()
                .any(|(path, _)| path == &companion)
        );
        let promoted = complete
            .calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
            .unwrap();
        assert_eq!(
            original.values.map(f64::to_bits),
            promoted.values.map(f64::to_bits)
        );
        assert_eq!(original.returned_flags, promoted.returned_flags);
        assert_eq!(planetary.context().leb_files().len(), 1);
        // Inspect the retained image after atomic replacement, not the new
        // pathname's index. Both directions preserve snapshot ownership.
        let retained_minor = complete
            .context()
            .leb_files()
            .iter()
            .find(|image| image.path() == companion.to_string_lossy())
            .unwrap()
            .clone();
        let replacement = directory.join("replacement");
        std::fs::write(&replacement, partition_fixture(SUN, None)).unwrap();
        std::fs::rename(&replacement, &companion).unwrap();
        assert!(partition_is_unused(
            &companion,
            Some(&retained_minor),
            SessionResources::PlanetaryCore
        ));
        assert!(!partition_is_unused(
            &companion,
            None,
            SessionResources::PlanetaryCore
        ));
        let retained_planet = LebFile::open(&companion).unwrap();
        std::fs::write(&replacement, partition_fixture(13, None)).unwrap();
        std::fs::rename(&replacement, &companion).unwrap();
        assert!(!partition_is_unused(
            &companion,
            Some(&retained_planet),
            SessionResources::PlanetaryCore
        ));
        assert!(partition_is_unused(
            &companion,
            None,
            SessionResources::PlanetaryCore
        ));
        runtime.set_leb_file(Some(companion.clone()));
        let explicit = runtime
            .shared_session_for(SessionResources::PlanetaryCore)
            .unwrap();
        assert!(
            explicit
                .context()
                .leb_files()
                .iter()
                .any(|file| file.path() == companion.to_string_lossy())
        );
        runtime.set_leb_file(None);
        std::fs::write(&companion, b"corrupt").unwrap();
        assert!(
            runtime
                .shared_session_for(SessionResources::PlanetaryCore)
                .is_err()
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn true_node_acquisition_preserves_values_and_promotes_dependencies() {
        use crate::constants::{FLG_EQUATORIAL, FLG_J2000, FLG_TOPOCTR, TRUE_NODE};
        let directory = discovery_directory("true-node-partitions");
        let mut core = std::fs::read(crate::test_data::core_leb_path()).unwrap();
        // Disable the known-file manifest guard without changing coefficients.
        core.push(0);
        std::fs::write(directory.join("base_core.leb2"), core).unwrap();
        let companion = directory.join("base_exotics.leb2");
        std::fs::write(&companion, partition_fixture(13, None)).unwrap();
        for topocentric in [false, true] {
            let runtime = discovery_runtime(&directory);
            runtime.set_iers_delta_t_enabled(Some(false));
            if topocentric {
                runtime.set_topo(12.5, 41.9, 100.0).unwrap();
            }
            let mut expected = Vec::new();
            for body in [TRUE_NODE, -TRUE_NODE] {
                for layout in [0, FLG_EQUATORIAL, FLG_J2000] {
                    for speed in [0, FLG_SPEED] {
                        let flags =
                            FLG_SWIEPH | layout | speed | if topocentric { FLG_TOPOCTR } else { 0 };
                        for ut in [false, true] {
                            let result = if ut {
                                runtime.calc_ut(2451545.0, body, flags)
                            } else {
                                runtime.calc(2451545.0, body, flags)
                            }
                            .unwrap();
                            expected.push((body, flags, ut, result));
                        }
                    }
                }
            }
            let lunar = runtime
                .shared_session_for(SessionResources::LunarCore)
                .unwrap();
            assert_eq!(
                runtime.lock().cached_session.as_ref().unwrap().resources,
                SessionResources::LunarCore
            );
            assert_eq!(lunar.context().leb_files().len(), 1);
            let all = runtime.shared_session().unwrap();
            assert_eq!(all.context().leb_files().len(), 2);
            for (body, flags, ut, expected) in expected {
                let actual = if ut {
                    runtime.calc_ut(2451545.0, body, flags)
                } else {
                    runtime.calc(2451545.0, body, flags)
                }
                .unwrap();
                assert_eq!(
                    actual.values.map(f64::to_bits),
                    expected.values.map(f64::to_bits)
                );
                assert_eq!(actual.returned_flags, expected.returned_flags);
            }
            // Promotion preserves the immutable snapshot already held by a caller.
            assert_eq!(lunar.context().leb_files().len(), 1);
        }
        // A true-node channel in a companion is an actual dependency, regardless
        // of filename. Inspecting metadata never substitutes a calculated state.
        std::fs::write(&companion, partition_fixture(TRUE_NODE, None)).unwrap();
        let runtime = discovery_runtime(&directory);
        let lunar = runtime
            .shared_session_for(SessionResources::LunarCore)
            .unwrap();
        assert_eq!(lunar.context().leb_files().len(), 2);
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn cold_houses_acquire_only_consumed_resources_and_promote_on_body_use() {
        let directory = discovery_directory("cold-houses");
        let path = directory.join("explicit.leb2");
        let mut context = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_iers_delta_t(Some(false));
        context.set_leb_file(Some(path.clone()));
        let runtime = Runtime::from_context(context);
        let expected_delta = runtime.delta_t(2451545.0).unwrap();
        let expected_nutation = runtime
            .calc_ut(2451545.0, crate::constants::ECL_NUT, 0)
            .unwrap();
        let first = runtime.house_session('P', 0).unwrap();
        let expected = first
            .public_houses_ex(2451545.0, 41.9, 12.5, 'P', 0)
            .unwrap();
        assert!(Arc::ptr_eq(&first, &runtime.house_session('P', 0).unwrap()));
        assert!(runtime.lock().reader.is_empty());
        assert!(first.context().leb_files().is_empty());

        // A geometry snapshot must never satisfy a later data acquisition.
        // The absent file remains strict, including Sunshine and live anchors.
        for (system, flags) in [('I', 0), ('i', 0), ('P', crate::constants::FLG_SIDEREAL)] {
            assert!(matches!(
                runtime.house_session(system, flags),
                Err(Error::LebFileNotFound { .. })
            ));
        }
        assert!(matches!(
            runtime.calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED),
            Err(Error::LebFileNotFound { .. })
        ));
        std::fs::copy(crate::test_data::core_leb_path(), &path).unwrap();
        let acquired = runtime.shared_session().unwrap();
        assert_eq!(acquired.context().leb_files().len(), 1);
        assert!(Arc::ptr_eq(
            &acquired,
            &runtime.house_session('P', 0).unwrap()
        ));
        assert_eq!(
            acquired
                .public_houses_ex(2451545.0, 41.9, 12.5, 'P', 0)
                .unwrap(),
            expected
        );
        assert_eq!(runtime.delta_t(2451545.0).unwrap(), expected_delta);
        assert_eq!(
            runtime
                .calc_ut(2451545.0, crate::constants::ECL_NUT, 0)
                .unwrap()
                .values
                .map(f64::to_bits),
            expected_nutation.values.map(f64::to_bits),
        );
        runtime
            .calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
            .unwrap();
        // Promotion cannot change an active caller's image-free snapshot.
        assert!(first.context().leb_files().is_empty());
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn metadata_identifies_the_retained_image_after_lazy_planetary_acquisition() {
        let directory = discovery_directory("lazy-provenance");
        let path = directory.join("base_core.leb2");
        let bytes = std::fs::read(crate::test_data::published_core_path()).unwrap();
        std::fs::write(&path, &bytes).unwrap();
        let runtime = discovery_runtime(&directory);
        let initial = runtime
            .calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
            .unwrap();
        assert!(!runtime.lock().reader_manifest_matches.contains_key(&path));

        let mut changed = bytes;
        changed.push(0);
        let replacement = directory.join("replacement");
        std::fs::write(&replacement, changed).unwrap();
        std::fs::rename(&replacement, &path).unwrap();
        assert!(
            runtime
                .body_coverage(SUN, Some(2451545.0))
                .unwrap()
                .unwrap()
                .reviewed
        );
        assert_eq!(
            runtime
                .calc_ut(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
                .unwrap()
                .values
                .map(f64::to_bits),
            initial.values.map(f64::to_bits),
        );
        runtime.reload_configuration();
        assert!(
            !runtime
                .body_coverage(SUN, Some(2451545.0))
                .unwrap()
                .unwrap()
                .reviewed
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn companion_admission_consumes_core_identity_before_ignoring_mismatches() {
        let directory = discovery_directory("lazy-admission");
        let core = directory.join("base_core.leb2");
        let bytes = std::fs::read(crate::test_data::published_core_path()).unwrap();
        std::fs::write(&core, &bytes).unwrap();
        let companion = directory.join("base_exotics.leb2");
        std::fs::write(&companion, partition_fixture(SUN, None)).unwrap();
        let runtime = discovery_runtime(&directory);
        let snapshot = runtime
            .shared_session_for(SessionResources::PlanetaryCore)
            .unwrap();
        assert_eq!(snapshot.context().leb_files().len(), 1);
        assert_eq!(
            runtime.lock().reader_manifest_matches.get(&core),
            Some(&true)
        );
        std::fs::write(&companion, b"corrupt").unwrap();
        runtime.reload_configuration();
        assert!(runtime.shared_session().is_ok());

        // An unreviewed but structurally valid core grants no manifest guard.
        let mut unreviewed = bytes;
        unreviewed.push(0);
        std::fs::write(&core, unreviewed).unwrap();
        runtime.reload_configuration();
        assert!(runtime.shared_session().is_err());
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn cold_house_snapshot_retires_on_clock_edits_and_is_shared_across_threads() {
        let runtime = Runtime::from_context(EphemerisContext::new().with_iers_delta_t(Some(false)));
        let first = runtime.house_session('P', 0).unwrap();
        std::thread::scope(|scope| {
            for _ in 0..8 {
                let runtime = &runtime;
                let first = &first;
                scope.spawn(move || {
                    assert!(Arc::ptr_eq(first, &runtime.house_session('K', 0).unwrap()));
                });
            }
        });
        let old_delta = first.deltat(2451545.0).unwrap();
        runtime.set_delta_t_userdef(Some(0.125));
        let changed = runtime.house_session('P', 0).unwrap();
        assert!(!Arc::ptr_eq(&first, &changed));
        assert_eq!(changed.deltat(2451545.0).unwrap(), 0.125);
        assert_eq!(first.deltat(2451545.0).unwrap(), old_delta);
        let retired = Arc::downgrade(&first);
        drop(first);
        assert!(retired.upgrade().is_none());
    }
    use crate::constants::{FLG_SPEED, FLG_SWIEPH, SUN};

    #[test]
    fn noop_reset_keeps_every_cached_resource_scope() {
        let path = crate::test_data::core_leb_path();
        let image = LebFile::open(&path).unwrap();
        for resources in [
            SessionResources::Models,
            SessionResources::PlanetaryCore,
            SessionResources::LunarCore,
            SessionResources::All,
        ] {
            let mut context = EphemerisContext::new()
                .with_mode(ProviderMode::SealedLeb)
                .with_iers_delta_t(Some(false));
            context.attach_leb_file(image.clone());
            let runtime = Runtime::from_context(context);
            let before = if resources == SessionResources::Models {
                runtime.model_session().unwrap()
            } else {
                runtime.shared_session_for(resources).unwrap()
            };
            let generation = runtime.lock().reader_generation;
            runtime.reset_session();
            let state = runtime.lock();
            let retained = state.cached_session.as_ref().unwrap();
            assert!(Arc::ptr_eq(&before, &retained.session));
            assert_eq!(retained.resources, resources);
            assert_eq!(state.reader_generation, generation);
            assert!(state.context.observer().is_none());
            assert!(state.context.angles_cache().is_empty());
        }
    }

    #[test]
    fn repeated_topo_preserves_only_identical_valid_observer_words() {
        let runtime =
            Runtime::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        runtime.set_topo(0.0, 0.0, -0.0).unwrap();
        let mut before = runtime.model_session().unwrap();
        runtime.set_topo(0.0, 0.0, -0.0).unwrap();
        assert!(Arc::ptr_eq(&before, &runtime.model_session().unwrap()));
        for (lon, lat, alt) in [
            (f64::NAN, 0.0, -0.0),
            (0.0, 91.0, -0.0),
            (0.0, 0.0, f64::INFINITY),
        ] {
            assert!(runtime.set_topo(lon, lat, alt).is_err());
            assert!(Arc::ptr_eq(&before, &runtime.model_session().unwrap()));
        }
        for (lon, lat, alt) in [(-0.0, 0.0, -0.0), (-0.0, -0.0, -0.0), (-0.0, -0.0, 0.0)] {
            let old_observer = before.context().observer().unwrap();
            runtime.set_topo(lon, lat, alt).unwrap();
            let after = runtime.model_session().unwrap();
            assert!(!Arc::ptr_eq(&before, &after));
            let observer = after.context().observer().unwrap();
            assert_eq!(observer.lon_deg().to_bits(), lon.to_bits());
            assert_eq!(observer.lat_deg().to_bits(), lat.to_bits());
            assert_eq!(observer.height_m().to_bits(), alt.to_bits());
            assert_eq!(before.context().observer(), Some(old_observer));
            runtime.set_topo(lon, lat, alt).unwrap();
            assert!(Arc::ptr_eq(&after, &runtime.model_session().unwrap()));
            before = after;
        }
    }

    #[test]
    fn repeated_sidereal_preserves_normalized_words_and_repeats_advisories() {
        use crate::constants::{SIDBIT_ECL_T0, SIDM_USER};
        let runtime =
            Runtime::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let tropical = runtime.model_session().unwrap();
        assert!(runtime.set_sid_mode(0, 0.0, 0.0).is_empty());
        let explicit = runtime.model_session().unwrap();
        assert!(!Arc::ptr_eq(&tropical, &explicit));
        assert!(matches!(
            tropical.context().sidereal(),
            crate::context::SiderealConfig::Tropical
        ));
        assert!(runtime.set_sid_mode(0, 0.0, 0.0).is_empty());
        assert!(Arc::ptr_eq(&explicit, &runtime.model_session().unwrap()));

        let nan = f64::from_bits(0x7ff8_0000_0000_0042);
        let mode = SIDM_USER | SIDBIT_ECL_T0;
        let warning = vec![Warning::UnsupportedSiderealProjection { bits: 0x4000 }];
        assert_eq!(runtime.set_sid_mode(mode | 0x4000, nan, -0.0), warning);
        let mut before = runtime.model_session().unwrap();
        for _ in 0..2 {
            assert_eq!(runtime.set_sid_mode(mode | 0x4000, nan, -0.0), warning);
            assert!(Arc::ptr_eq(&before, &runtime.model_session().unwrap()));
        }
        assert!(runtime.set_sid_mode(mode, nan, -0.0).is_empty());
        assert!(Arc::ptr_eq(&before, &runtime.model_session().unwrap()));
        for (mode, t0, ayan) in [
            (mode, f64::from_bits(nan.to_bits() + 1), -0.0),
            (mode, -0.0, -0.0),
            (mode, 0.0, -0.0),
            (mode, 0.0, 0.0),
            (SIDM_USER, 0.0, 0.0),
            (-1, 0.0, 0.0),
        ] {
            let old = before.context().sidereal_mode_full();
            assert!(runtime.set_sid_mode(mode, t0, ayan).is_empty());
            let after = runtime.model_session().unwrap();
            assert!(!Arc::ptr_eq(&before, &after));
            assert_eq!(before.context().sidereal_mode_full().0, old.0);
            assert_eq!(
                before.context().sidereal_mode_full().1.to_bits(),
                old.1.to_bits()
            );
            assert_eq!(
                before.context().sidereal_mode_full().2.to_bits(),
                old.2.to_bits()
            );
            before = after;
        }
    }

    #[test]
    fn actual_reset_clears_each_transient_role_and_retires_only_its_snapshot() {
        for role in 0..4 {
            let runtime =
                Runtime::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
            match role {
                0 => runtime.set_topo(12.5, 41.9, 50.0).unwrap(),
                1 => {
                    runtime.set_sid_mode(0, 0.0, 0.0);
                }
                2 => runtime.set_angles_cache(&HashMap::from([("mc".to_owned(), -0.0)])),
                _ => {
                    // Client construction performs no I/O. Keep only the runtime's
                    // strong reference so the reset must actually release it.
                    runtime.lock().horizons_client =
                        Some(Arc::new(crate::horizons::HorizonsClient::new().unwrap()));
                }
            }
            let before = runtime.model_session().unwrap();
            let retired = Arc::downgrade(&before);
            let client = runtime.lock().horizons_client.as_ref().map(Arc::downgrade);
            runtime.reset_session();
            let after = runtime.model_session().unwrap();
            assert!(!Arc::ptr_eq(&before, &after));
            assert!(after.context().observer().is_none());
            assert!(matches!(
                after.context().sidereal(),
                crate::context::SiderealConfig::Tropical
            ));
            assert!(after.context().angles_cache().is_empty());
            assert!(runtime.lock().horizons_client.is_none());
            if let Some(client) = client {
                assert!(client.upgrade().is_none());
            }
            match role {
                0 => assert!(before.context().observer().is_some()),
                1 => assert!(matches!(
                    before.context().sidereal(),
                    crate::context::SiderealConfig::Sidereal(_)
                )),
                2 => assert!(!before.context().angles_cache().is_empty()),
                _ => (),
            }
            runtime.reset_session();
            assert!(Arc::ptr_eq(&after, &runtime.model_session().unwrap()));
            drop(before);
            assert!(retired.upgrade().is_none());
        }
    }

    #[test]
    fn metadata_and_unchanged_selections_retain_the_same_snapshot() {
        let path = crate::test_data::core_leb_path();
        let runtime = Runtime::default();
        runtime.set_calc_mode(Some("leb")).unwrap();
        runtime.set_leb_file(Some(path.clone()));
        runtime.set_ephe_path(Some(path.parent().unwrap().to_owned()));
        runtime.set_iers_delta_t_enabled(Some(false));
        let snapshot = runtime.shared_session().unwrap();
        runtime.leb_reader().unwrap();
        assert!(Arc::ptr_eq(&snapshot, &runtime.shared_session().unwrap()));
        runtime.current_file_data(0).unwrap();
        runtime.current_file_data(1).unwrap();
        runtime.current_file_data(9).unwrap();
        assert!(Arc::ptr_eq(&snapshot, &runtime.shared_session().unwrap()));
        runtime.set_calc_mode(Some(" LEB ")).unwrap();
        runtime.set_leb_file(Some(path.clone()));
        runtime.set_ephe_path(Some(path.parent().unwrap().to_owned()));
        assert!(Arc::ptr_eq(&snapshot, &runtime.shared_session().unwrap()));
        assert!(runtime.set_calc_mode(Some("invalid")).is_err());
        assert!(Arc::ptr_eq(&snapshot, &runtime.shared_session().unwrap()));
        runtime.set_ephe_path(None);
        assert!(!Arc::ptr_eq(&snapshot, &runtime.shared_session().unwrap()));
        // A retained caller still owns its selected image after real edits.
        assert_eq!(
            PathBuf::from(snapshot.context().leb_files()[0].path()),
            path
        );
    }

    #[test]
    fn explicit_and_deferred_auto_mode_remain_distinct_settings() {
        let runtime =
            Runtime::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let first = runtime.shared_session().unwrap();
        runtime.set_calc_mode(Some("skyfield")).unwrap();
        let explicit = runtime.shared_session().unwrap();
        assert!(Arc::ptr_eq(&first, &explicit));
        runtime.set_calc_mode(Some("SKYFIELD")).unwrap();
        assert!(Arc::ptr_eq(&explicit, &runtime.shared_session().unwrap()));
        runtime.set_calc_mode(None).unwrap();
        assert!(runtime.lock().mode.is_none());
        runtime.set_calc_mode(Some("auto")).unwrap();
        assert_eq!(runtime.lock().mode, Some(ProviderMode::Auto));
    }

    #[test]
    fn shared_snapshot_is_reused_and_old_generations_are_released() {
        let runtime =
            Runtime::from_context(EphemerisContext::new().with_mode(ProviderMode::NativeJpl));
        let first = runtime.shared_session().unwrap();
        let again = runtime.shared_session().unwrap();
        assert!(Arc::ptr_eq(&first, &again));
        runtime.effective_context().unwrap();
        assert!(Arc::ptr_eq(&first, &runtime.shared_session().unwrap()));
        let retired = Arc::downgrade(&first);
        runtime.set_topo(12.5, 41.9, 50.0).unwrap();
        let changed = runtime.shared_session().unwrap();
        assert!(!Arc::ptr_eq(&first, &changed));
        assert!(first.context().observer().is_none());
        assert_eq!(changed.context().observer().unwrap().height_m(), 50.0);
        drop(first);
        drop(again);
        assert!(retired.upgrade().is_none());
    }

    #[test]
    fn cached_and_fresh_snapshots_match_after_configuration_sequences() {
        let image = LebFile::open(&crate::test_data::core_leb_path()).unwrap();
        let mut context = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_iers_delta_t(Some(false));
        context.attach_leb_file(image);
        let runtime = Runtime::from_context(context);
        for index in 0..32 {
            match index % 6 {
                0 => runtime
                    .set_topo(12.5 + index as f64 / 100.0, 41.9, index as f64)
                    .unwrap(),
                1 => {
                    runtime.set_sid_mode(index % 5, 2451545.0, 0.0);
                }
                2 => runtime.set_delta_t_userdef(Some((60.0 + index as f64) / 86400.0)),
                3 => runtime.set_tid_acc(-25.8 + index as f64 / 100.0),
                4 => runtime.set_calc_mode(Some("leb")).unwrap(),
                _ => runtime.set_delta_t_userdef(None),
            }
            let fresh = Runtime::from_context(runtime.effective_context().unwrap());
            for flags in [
                FLG_SWIEPH,
                FLG_SWIEPH | FLG_SPEED,
                FLG_SWIEPH | FLG_SPEED | crate::constants::FLG_EQUATORIAL,
                FLG_SWIEPH | FLG_SPEED | crate::constants::FLG_SIDEREAL,
            ] {
                let epoch = 2451545.0 + index as f64 / 64.0;
                let cached = runtime
                    .calc_ut(epoch, crate::constants::MARS, flags)
                    .unwrap();
                let uncached = fresh.calc_ut(epoch, crate::constants::MARS, flags).unwrap();
                assert_eq!(
                    cached.values.map(f64::to_bits),
                    uncached.values.map(f64::to_bits)
                );
                assert_eq!(cached.returned_flags, uncached.returned_flags);
            }
        }
    }

    #[test]
    fn configuration_edit_during_acquisition_is_nonblocking_and_discards_stale_failure() {
        let directory = discovery_directory("concurrent-acquisition");
        let good = directory.join("base_core.leb2");
        std::fs::copy(crate::test_data::core_leb_path(), &good).unwrap();
        let runtime = Arc::new(discovery_runtime(&directory));
        runtime.set_leb_file(Some(directory.join("missing.leb2")));
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (resume_tx, resume_rx) = std::sync::mpsc::channel();
        let worker_runtime = Arc::clone(&runtime);
        let worker = std::thread::spawn(move || {
            let calls = std::cell::Cell::new(0);
            let (_state, context) = worker_runtime
                .reader_snapshot_with(|state, context| {
                    if calls.get() == 0 {
                        started_tx.send(()).unwrap();
                        resume_rx.recv().unwrap();
                    }
                    calls.set(calls.get() + 1);
                    worker_runtime.attach_readers(state, context)
                })
                .unwrap();
            (context.leb_file().cloned(), calls.get())
        });
        started_rx
            .recv_timeout(std::time::Duration::from_secs(5))
            .unwrap();
        let (edited_tx, edited_rx) = std::sync::mpsc::channel();
        let setter_runtime = Arc::clone(&runtime);
        let selected = good.clone();
        let setter = std::thread::spawn(move || {
            setter_runtime.set_leb_file(Some(selected));
            edited_tx.send(()).unwrap();
        });
        let nonblocking = edited_rx
            .recv_timeout(std::time::Duration::from_secs(2))
            .is_ok();
        // Always release/join the workers, even if a regression held the lock.
        resume_tx.send(()).unwrap();
        setter.join().unwrap();
        let (selected, calls) = worker.join().unwrap();
        assert!(
            nonblocking,
            "configuration lock was held during file acquisition"
        );
        assert_eq!(selected, Some(good));
        assert_eq!(calls, 2, "the stale missing-file error must be retried");
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn indexed_policy_changes_future_acquisitions_not_existing_snapshots() {
        use crate::leb::{IndexedLimits, LebReadPolicy};
        let directory = discovery_directory("indexed-policy");
        std::fs::copy(
            crate::test_data::core_leb_path(),
            directory.join("base_core.leb2"),
        )
        .unwrap();
        let runtime = discovery_runtime(&directory);
        let owned = runtime.session().unwrap();
        let expected = owned.calc(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED).unwrap();
        assert!(owned.context().leb_files()[0].allocated_image_bytes() > 0);
        let expected_inventory = runtime.leb_inventory().unwrap();
        runtime
            .set_leb_read_policy(Some(LebReadPolicy::Indexed(IndexedLimits::default())))
            .unwrap();
        let indexed = runtime.session().unwrap();
        assert_eq!(indexed.context().leb_files()[0].allocated_image_bytes(), 0);
        assert!(matches!(
            indexed.context().leb_files()[0].bytes(),
            Err(Error::LebBytesUnavailable)
        ));
        assert_eq!(runtime.leb_inventory().unwrap(), expected_inventory);
        assert_eq!(
            indexed
                .calc(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
                .unwrap(),
            expected
        );
        assert_eq!(
            owned.calc(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED).unwrap(),
            expected
        );
        runtime
            .set_leb_read_policy(Some(LebReadPolicy::Owned))
            .unwrap();
        assert!(runtime.session().unwrap().context().leb_files()[0].allocated_image_bytes() > 0);
        assert_eq!(
            indexed
                .calc(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
                .unwrap(),
            expected
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn auto_mode_does_not_turn_an_indexed_resource_refusal_into_provider_fallback() {
        use crate::leb::{IndexedLimits, LebReadPolicy};
        let directory = discovery_directory("indexed-refusal");
        let path = directory.join("base_core.leb2");
        std::fs::copy(crate::test_data::core_leb_path(), &path).unwrap();
        let mut context = EphemerisContext::new().with_mode(ProviderMode::Auto);
        context.set_leb_file(Some(path));
        context.set_network_policy(Some(crate::context::NetworkPolicy::Sealed));
        context
            .set_leb_read_policy(Some(LebReadPolicy::Indexed(IndexedLimits {
                frame_bytes: 0,
                ..IndexedLimits::default()
            })))
            .unwrap();
        let error = Runtime::from_context(context)
            .calc(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
            .unwrap_err();
        assert!(
            matches!(error.primary_error(), Error::ResourceLimit { .. }),
            "{error:?}"
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn indexed_configuration_survives_resolution_without_changing_scalar_results() {
        let directory = discovery_directory("indexed-sources");
        std::fs::create_dir(directory.join("leb")).unwrap();
        std::fs::copy(
            crate::test_data::published_core_path(),
            directory.join("leb/base_core.leb2"),
        )
        .unwrap();
        let sources = ConfigurationSources::load(
            std::collections::BTreeMap::from([
                (
                    "LIBEPHEMERIS_DATA_DIR".into(),
                    directory.display().to_string(),
                ),
                ("LIBEPHEMERIS_LEB_READER".into(), "indexed".into()),
                ("LIBEPHEMERIS_MODE".into(), "leb".into()),
                ("LIBEPHEMERIS_PRECISION".into(), "base".into()),
                ("LIBEPHEMERIS_IERS_DELTA_T".into(), "0".into()),
            ]),
            None,
            None,
        );
        let runtime = Runtime::with_sources(sources);
        let session = runtime.session().unwrap();
        assert_eq!(session.context().leb_files()[0].allocated_image_bytes(), 0);
        assert!(matches!(
            session.context().leb_files()[0].bytes(),
            Err(Error::LebBytesUnavailable)
        ));
        assert!(
            runtime
                .body_coverage(SUN, Some(2451545.0))
                .unwrap()
                .unwrap()
                .reviewed
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn remote_mode_does_not_open_an_explicit_leb_and_reset_retains_clock_settings() {
        let runtime = Runtime::default();
        runtime.set_calc_mode(Some("horizons")).unwrap();
        runtime.set_leb_file(Some(PathBuf::from("missing.leb2")));
        assert!(runtime.session().is_ok());
        assert!(runtime.leb_reader().unwrap().is_none());
        runtime.set_topo(12.0, 42.0, 100.0).unwrap();
        runtime.set_delta_t_userdef(Some(64.0 / 86400.0));
        runtime.reset_session();
        assert!(runtime.context().observer().is_none());
        assert_eq!(runtime.context().user_delta_t_days(), Some(64.0 / 86400.0));
        assert!(runtime.set_spk_date_padding(-1).is_err());
        assert_eq!(runtime.context().spk_date_padding_days(), 0);
    }

    #[test]
    fn class_resource_release_preserves_independent_configuration() {
        let runtime = Runtime::shared(EphemerisContext::new());
        runtime.set_topo(12.0, 42.0, 100.0).unwrap();
        runtime.set_sid_mode(255, 2451545.0, 24.0);
        let before = runtime.context();
        Runtime::release_shared_resources();
        assert_eq!(runtime.context().observer(), before.observer());
        assert_eq!(
            runtime.context().sidereal_mode_full(),
            before.sidereal_mode_full()
        );
    }

    #[test]
    fn invalid_settings_leave_independent_contexts_unchanged() {
        let left = Runtime::default();
        let right = Runtime::default();
        left.set_delta_t_userdef(Some(64.0 / 86400.0));
        assert_eq!(right.context().user_delta_t_days(), None);
        let before = left.context();
        assert!(left.set_calc_mode(Some("invalid")).is_err());
        assert!(left.set_precision_tier("invalid").is_err());
        assert_eq!(left.context(), before);
    }

    #[test]
    fn restoring_the_same_resolved_filename_keeps_the_owned_main_image() {
        let directory =
            std::env::temp_dir().join(format!("libephemeris-main-setting-{}", std::process::id()));
        std::fs::create_dir_all(&directory).unwrap();
        let path = directory.join("de440.bsp");
        let mut bytes = crate::spk::moon_test_kernel_bytes();
        bytes[1064..1068].copy_from_slice(&10_i32.to_le_bytes());
        bytes[1068..1072].copy_from_slice(&0_i32.to_le_bytes());
        std::fs::write(&path, bytes).unwrap();
        let sources = ConfigurationSources::load(
            std::collections::BTreeMap::from([(
                "LIBEPHEMERIS_DATA_DIR".into(),
                directory.display().to_string(),
            )]),
            None,
            None,
        );
        let runtime = Runtime::with_sources(sources);
        runtime.set_calc_mode(Some("skyfield")).unwrap();
        runtime.set_precision_tier("base").unwrap();
        runtime.set_ephemeris_file("de440.bsp");
        let first = runtime.main_image().unwrap();
        std::fs::remove_file(path).unwrap();
        runtime.set_precision_tier("base").unwrap();
        assert!(runtime.lock().context.ephemeris_file_override().is_none());
        runtime.set_ephemeris_file("de440.bsp");
        let second = runtime.main_image().unwrap();
        assert_eq!(first.registry_key(), second.registry_key());
        assert!(runtime.main_image.current_owner().is_some());
        std::fs::remove_dir(directory).unwrap();
    }

    fn discovery_runtime(directory: &std::path::Path) -> Runtime {
        let mut context = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_tier(crate::context::Tier::Base);
        context.set_leb_dir(Some(directory.to_path_buf()));
        Runtime::from_context(context)
    }

    fn discovery_directory(label: &str) -> PathBuf {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static NEXT: AtomicUsize = AtomicUsize::new(0);
        let path = std::env::temp_dir().join(format!(
            "libephemeris-discovery-{label}-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).unwrap();
        path
    }

    #[test]
    fn discovered_images_and_digests_are_shared_until_explicit_reload() {
        let directory = discovery_directory("owners");
        let path = directory.join("base_core.leb2");
        let original = std::fs::read(crate::test_data::published_core_path()).unwrap();
        std::fs::write(&path, &original).unwrap();
        let runtime = discovery_runtime(&directory);
        let first = runtime.session().unwrap();
        let file = &first.context().leb_files()[0];
        let expected_digest = file.sha256().unwrap().to_owned();
        let expected = first.calc(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED).unwrap();
        assert!(
            runtime
                .body_coverage(SUN, Some(2451545.0))
                .unwrap()
                .unwrap()
                .reviewed
        );
        let pointer = file.bytes().unwrap().as_ptr() as usize;

        // The disappeared directory must not be probed again by warm calls.
        // All calls use the already-owned immutable image, including its digest.
        let hidden = directory.with_extension("retired");
        std::fs::rename(&directory, &hidden).unwrap();
        std::thread::scope(|scope| {
            let calls: Vec<_> = (0..8)
                .map(|_| {
                    scope.spawn(|| {
                        let snapshot = runtime.session().unwrap();
                        assert_eq!(
                            snapshot.context().leb_files()[0].bytes().unwrap().as_ptr() as usize,
                            pointer
                        );
                        assert!(
                            runtime
                                .body_coverage(SUN, Some(2451545.0))
                                .unwrap()
                                .unwrap()
                                .reviewed
                        );
                        snapshot
                            .calc(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED)
                            .unwrap()
                    })
                })
                .collect();
            for call in calls {
                assert_eq!(call.join().unwrap(), expected);
            }
        });
        std::fs::rename(&hidden, &directory).unwrap();

        // Atomic pathname replacement cannot change an active image's identity.
        let replacement = directory.join("replacement.leb2");
        let mut changed = original;
        changed.push(0); // Valid container, different complete-file identity.
        std::fs::write(&replacement, changed).unwrap();
        std::fs::rename(replacement, &path).unwrap();
        assert!(
            runtime
                .body_coverage(SUN, Some(2451545.0))
                .unwrap()
                .unwrap()
                .reviewed
        );
        runtime.reload_configuration();
        let replaced = runtime.session().unwrap();
        assert_ne!(
            replaced.context().leb_files()[0].bytes().unwrap().as_ptr() as usize,
            pointer
        );
        assert_ne!(
            replaced.context().leb_files()[0].sha256().unwrap(),
            expected_digest
        );
        assert!(
            !runtime
                .body_coverage(SUN, Some(2451545.0))
                .unwrap()
                .unwrap()
                .reviewed
        );
        assert_eq!(file.sha256().unwrap(), expected_digest);
        assert_eq!(
            first.calc(2451545.0, SUN, FLG_SWIEPH | FLG_SPEED).unwrap(),
            expected
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn discovery_invalidates_for_tier_path_and_failed_explicit_files_are_retried() {
        let directory = discovery_directory("invalidation");
        let runtime = discovery_runtime(&directory);
        runtime.session().unwrap();
        let generation = runtime.lock().reader_generation;
        runtime.set_precision_tier("medium").unwrap();
        assert_ne!(runtime.lock().reader_generation, generation);
        assert!(runtime.lock().reader_paths.is_none());
        runtime.session().unwrap();
        runtime.set_ephe_path(Some(directory.join("kernels")));
        assert!(runtime.lock().reader_paths.is_none());

        let missing = directory.join("missing.leb2");
        runtime.set_leb_file(Some(missing.clone()));
        assert!(matches!(
            runtime.session(),
            Err(Error::LebFileNotFound { .. })
        ));
        std::fs::write(&missing, b"corrupt").unwrap();
        assert!(runtime.session().is_err());
        std::fs::copy(crate::test_data::core_leb_path(), &missing).unwrap();
        assert_eq!(
            runtime.session().unwrap().context().leb_files()[0].path(),
            missing.to_string_lossy()
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn auto_skips_unavailable_explicit_leb_and_retries_the_same_path_after_repair() {
        let directory = discovery_directory("auto-explicit-repair");
        let runtime = discovery_runtime(&directory);
        runtime.set_calc_mode(Some("auto")).unwrap();
        let path = directory.join("selected.leb2");
        runtime.set_leb_file(Some(path.clone()));
        assert!(runtime.session().unwrap().context().leb_files().is_empty());
        assert!(runtime.leb_reader().unwrap().is_none());
        std::fs::write(&path, b"corrupt").unwrap();
        assert!(runtime.session().unwrap().context().leb_files().is_empty());
        assert!(runtime.leb_reader().unwrap().is_none());
        std::fs::copy(crate::test_data::core_leb_path(), &path).unwrap();
        assert_eq!(
            runtime.leb_reader().unwrap().unwrap().path(),
            path.to_string_lossy()
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[cfg(feature = "bundled-data")]
    #[test]
    fn installing_a_core_invalidates_an_empty_discovery_without_activation() {
        let directory = discovery_directory("install");
        let mut context = EphemerisContext::new()
            .with_mode(ProviderMode::SealedLeb)
            .with_tier(crate::context::Tier::Base);
        context.set_jpl_data_dir(Some(directory.clone()));
        let runtime = Runtime::from_context(context);
        runtime.session().unwrap();
        assert!(runtime.lock().reader_paths.as_ref().unwrap().1.is_empty());
        let installed = runtime.download_leb_for_tier("base", false, false).unwrap();
        assert!(runtime.lock().reader_paths.is_none());
        assert_eq!(
            runtime.session().unwrap().context().leb_files()[0].path(),
            installed.to_string_lossy()
        );
        std::fs::remove_dir_all(directory).unwrap();
    }

    #[test]
    fn snapshots_keep_images_after_selection_changes_and_calls_are_parallel() {
        let path = crate::test_data::core_leb_path();
        let runtime = Runtime::default();
        runtime.set_calc_mode(Some("leb")).unwrap();
        runtime.set_leb_file(Some(path.clone()));
        runtime.set_iers_delta_t_enabled(Some(false));
        runtime.set_delta_t_userdef(Some(64.0 / 86400.0));
        let snapshot = runtime.session().unwrap();
        runtime.set_precision_tier("base").unwrap();
        let selected = runtime.session().unwrap();
        runtime.set_precision_tier("base").unwrap();
        runtime.set_leb_file(Some(path));
        let reused = runtime.session().unwrap();
        assert!(std::ptr::eq(
            selected.context().leb_files()[0].bytes().unwrap().as_ptr(),
            reused.context().leb_files()[0].bytes().unwrap().as_ptr(),
        ));

        let flags = FLG_SWIEPH | FLG_SPEED;
        let expected = snapshot.calc_ut(2451545.0, SUN, flags).unwrap();
        std::thread::scope(|scope| {
            let calls: Vec<_> = (0..8)
                .map(|_| scope.spawn(|| runtime.calc_ut(2451545.0, SUN, flags).unwrap()))
                .collect();
            for call in calls {
                assert_eq!(call.join().unwrap(), expected);
            }
        });
        runtime.set_leb_file(None);
        assert!(runtime.context().leb_files().is_empty());
        assert_eq!(snapshot.calc_ut(2451545.0, SUN, flags).unwrap(), expected);
    }
}

#[path = "config/module_positions.rs"]
mod module_positions;
#[path = "config/module_state.rs"]
pub mod module_state;
