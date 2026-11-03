//! Independent ephemeris configuration and retained data resources.
//!
//! [`EphemerisContext`] owns observer, sidereal selection, precision tier,
//! provider mode, data paths, Delta-T/IERS configuration and refraction
//! settings. Cloned contexts retain immutable resources while mutable
//! configuration is independently owned. [`crate::config`] resolves explicit
//! settings, environment, dotenv, TOML and defaults in Rust.
//!
//! Observer longitude and latitude are degrees, east and north positive;
//! altitude is metres. Delta-T overrides are days, tidal acceleration is
//! arcseconds per century squared, and lapse rate is K/m. The context itself
//! does not convert epochs or calculate positions.
//!
//! Validated setter/builder pairs preserve prior state on error. Repeating
//! an unchanged path retains its resource owners; changing a path invalidates
//! resources selected from that path. Session reset clears per-calculation
//! state; close also retires retained images without invalidating resources
//! already owned by active calculations.
//!
//! The public provider spelling `skyfield` selects [`ProviderMode::NativeJpl`],
//! the native JPL/SPK route. Sidereal modes, warnings, automatic sentinels
//! and numeric domains are documented on their corresponding entries.

use crate::constants::{DELTAT_AUTOMATIC, TIDAL_AUTOMATIC, TIDAL_DEFAULT};
use crate::errors::{Error, Result, Warning};
use crate::horizons::{HorizonsClient, HorizonsConfig};
use crate::leb::LebFile;
use crate::sidereal::SiderealSelection;
use crate::spk::SpkRegistry;
use crate::time::{DeltaTConfig, DeltaTModel};
use std::collections::HashMap;
use std::path::PathBuf;

/// Default atmospheric lapse rate, 0.0065 K/m (standard atmosphere).
/// Used when no explicit lapse rate is selected.
pub const LAPSE_RATE_DEFAULT: f64 = f64::from_bits(0x3f7a_9fbe_76c8_b439);

/// Environment variable selecting the calculation mode.
/// TOML key: [`CALC_MODE_TOML_KEY`].
pub const CALC_MODE_ENV_VAR: &str = "LIBEPHEMERIS_MODE";
/// TOML configuration key selecting the calculation mode.
pub const CALC_MODE_TOML_KEY: &str = "mode";
/// Environment variable selecting the precision tier.
/// TOML key: [`TIER_TOML_KEY`].
pub const TIER_ENV_VAR: &str = "LIBEPHEMERIS_PRECISION";
/// TOML configuration key selecting the precision tier (`precision`).
pub const TIER_TOML_KEY: &str = "precision";
/// Environment variable selecting the Delta-T model.
/// TOML key: [`DELTA_T_MODEL_TOML_KEY`].
pub const DELTA_T_MODEL_ENV_VAR: &str = "LIBEPHEMERIS_DELTAT_MODEL";
/// TOML configuration key selecting the Delta-T model.
pub const DELTA_T_MODEL_TOML_KEY: &str = "deltat_model";
/// Environment variable enabling observed-IERS Delta T.
/// TOML key: [`IERS_DELTA_T_TOML_KEY`].
pub const IERS_DELTA_T_ENV_VAR: &str = "LIBEPHEMERIS_IERS_DELTA_T";
/// TOML configuration key enabling observed-IERS Delta T.
pub const IERS_DELTA_T_TOML_KEY: &str = "iers_delta_t";
/// Environment variable selecting the network policy
/// (`LIBEPHEMERIS_NETWORK_POLICY`).
/// TOML key: [`NETWORK_POLICY_TOML_KEY`].
pub const NETWORK_POLICY_ENV_VAR: &str = "LIBEPHEMERIS_NETWORK_POLICY";
/// Alias of [`NETWORK_POLICY_ENV_VAR`] for the same configuration variable.
pub const NETWORK_POLICY_ENV: &str = NETWORK_POLICY_ENV_VAR;
/// TOML configuration key selecting the network policy
/// (`network_policy`, one of `"auto"` / `"allow"` / `"sealed"`).
pub const NETWORK_POLICY_TOML_KEY: &str = "network_policy";
/// Environment variable enabling automatic SPK download for minor bodies.
/// TOML key: [`AUTO_SPK_TOML_KEY`].
pub const AUTO_SPK_ENV_VAR: &str = "LIBEPHEMERIS_AUTO_SPK";
/// TOML configuration key enabling automatic SPK download.
pub const AUTO_SPK_TOML_KEY: &str = "auto_spk";
/// Environment variable enabling strict precision mode for major
/// asteroids. TOML key: [`STRICT_PRECISION_TOML_KEY`].
pub const STRICT_PRECISION_ENV_VAR: &str = "LIBEPHEMERIS_STRICT_PRECISION";
/// TOML configuration key enabling strict precision mode.
pub const STRICT_PRECISION_TOML_KEY: &str = "strict_precision";
/// Environment variable enabling automatic IERS data download.
/// TOML key: [`IERS_AUTO_DOWNLOAD_TOML_KEY`].
pub const IERS_AUTO_DOWNLOAD_ENV_VAR: &str = "LIBEPHEMERIS_IERS_AUTO_DOWNLOAD";
/// TOML configuration key enabling automatic IERS data download.
pub const IERS_AUTO_DOWNLOAD_TOML_KEY: &str = "iers_auto_download";

/// Source name for a programmatic (setter) configuration value in
/// [`Error::InvalidCalcMode`] messages.
pub const CALC_MODE_SOURCE_SETTER: &str = "set_calc_mode";
/// Source name for a programmatic (setter) configuration value in
/// [`Error::InvalidNetworkPolicy`] messages.
pub const NETWORK_POLICY_SOURCE_SETTER: &str = "set_network_policy";

/// Data-coverage tier (source tiers are `base`/`medium`/`extended`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Tier {
    /// Smallest coverage tier.
    Base,
    /// Intermediate coverage tier and default selection.
    #[default]
    Medium,
    /// Widest coverage tier.
    Extended,
}

impl Tier {
    /// Canonical tier name (`"base"` / `"medium"` / `"extended"`).
    pub fn as_name(self) -> &'static str {
        match self {
            Self::Base => "base",
            Self::Medium => "medium",
            Self::Extended => "extended",
        }
    }

    /// Parse an exact tier name without case or whitespace normalization.
    /// Unknown spellings, including `"Medium"` or `" medium"`, return `None`.
    /// Configuration-file resolution has its own normalization policy; see
    /// [`resolve_tier`].
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "base" => Some(Self::Base),
            "medium" => Some(Self::Medium),
            "extended" => Some(Self::Extended),
            _ => None,
        }
    }
}

/// Provider routing mode.
///
/// LEB reads only locally bound data. Native JPL uses owned SPK chains, and
/// Horizons uses the network policy configured for the request. The public
/// mode spelling `skyfield` selects the native JPL implementation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ProviderMode {
    /// Sealed local LEB data only.
    SealedLeb,
    /// Native JPL-backed path (SPK registry chains).
    NativeJpl,
    /// Explicit remote Horizons provider (NET transport).
    HorizonsRemote,
    /// Routed selection per the accepted `auto` contract (route-planning layer):
    /// sealed LEB when bound, remote Horizons exactly when no local JPL
    /// source exists, native JPL otherwise. This is the release default.
    #[default]
    Auto,
}

/// Normalize a mode name with native Unicode whitespace and lowercase rules.
fn normalized_calc_mode_name(value: &str) -> String {
    value.trim().to_lowercase()
}

impl ProviderMode {
    /// Public calculation-mode name.
    /// `skyfield` is the retained spelling for the native JPL/SPK route.
    pub fn as_calc_mode_name(self) -> &'static str {
        match self {
            Self::SealedLeb => "leb",
            Self::NativeJpl => "skyfield",
            Self::HorizonsRemote => "horizons",
            Self::Auto => "auto",
        }
    }

    /// Parse a release calc-mode name with the `set_calc_mode`
    /// normalization (Unicode lowercase/Unicode-whitespace strip). Returns `None`
    /// for unknown spellings; the caller maps that to
    /// [`Error::InvalidCalcMode`], mirroring the release `ValueError`
    /// (a typo must fail loudly, never widen a sealed deployment back to
    /// `auto`).
    pub fn from_calc_mode_name(name: &str) -> Option<Self> {
        match normalized_calc_mode_name(name).as_str() {
            "auto" => Some(Self::Auto),
            "skyfield" => Some(Self::NativeJpl),
            "leb" => Some(Self::SealedLeb),
            "horizons" => Some(Self::HorizonsRemote),
            _ => None,
        }
    }
}

/// Process network policy.
///
/// `"auto"` seals the sealed-LEB mode and allows
/// every other mode; `"sealed"` and `"allow"` are explicit overrides.
/// The `download` provisioning command selects `"allow"` because a
/// download is an explicit provisioning action, not calculation-time
/// egress (the download-provisioning layer owns that call site).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum NetworkPolicy {
    /// Mode-dependent default: sealed under sealed LEB, allowed otherwise.
    #[default]
    Auto,
    /// Network access is permitted.
    Allow,
    /// No network access: every gated operation fails before socket work.
    Sealed,
}

impl NetworkPolicy {
    /// Canonical release policy spelling (`"auto"` / `"allow"` / `"sealed"`).
    pub fn as_policy_name(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Allow => "allow",
            Self::Sealed => "sealed",
        }
    }

    /// Parse a release policy spelling with the `set_network_policy`
    /// normalization (ASCII case/whitespace-insensitive). Returns `None`
    /// for unknown spellings; the caller maps that to
    /// [`Error::InvalidNetworkPolicy`], mirroring the release `ValueError`
    /// (a typo in a security-relevant boundary must never silently widen
    /// it to `"auto"`/`"allow"`).
    pub fn from_policy_name(name: &str) -> Option<Self> {
        match name.trim().to_ascii_lowercase().as_str() {
            "auto" => Some(Self::Auto),
            "allow" => Some(Self::Allow),
            "sealed" => Some(Self::Sealed),
            _ => None,
        }
    }
}

/// Effective network gate for one calculation mode.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EffectiveNetworkPolicy {
    /// Network access is permitted.
    Allow,
    /// Network I/O must fail before any socket work.
    Sealed,
}

/// Validated geodetic site: east longitude/north latitude in degrees, height
/// in metres above the reference ellipsoid used by the consuming operation.
///
/// Construction goes through [`Observer::new`]. This value stores a site; it
/// neither selects topocentric output nor computes an Earth-rotation clock.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Observer {
    lon_deg: f64,
    lat_deg: f64,
    height_m: f64,
}

impl Observer {
    /// Validate and retain a geographic site without normalizing its values.
    ///
    /// # Arguments
    ///
    /// - `lon_deg`: east-positive longitude in the inclusive range
    ///   `[-180, 360]` degrees. Signed west longitudes and east-only longitudes
    ///   are both accepted; the consuming geometry is periodic in 360 degrees.
    /// - `lat_deg`: north-positive latitude in `[-90, 90]` degrees.
    /// - `height_m`: finite ellipsoidal height in metres; negative is allowed.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidObserver`] for nonfinite or out-of-range input,
    /// checking latitude, longitude and height in that order. No data are read
    /// and no coordinate conversion or physical site survey is performed.
    pub fn new(lon_deg: f64, lat_deg: f64, height_m: f64) -> Result<Self> {
        if !lat_deg.is_finite() || lat_deg < -90.0 || lat_deg > 90.0 {
            return Err(Error::InvalidObserver {
                field: "lat_deg",
                value: lat_deg,
            });
        }
        if !lon_deg.is_finite() || lon_deg < -180.0 || lon_deg > 360.0 {
            return Err(Error::InvalidObserver {
                field: "lon_deg",
                value: lon_deg,
            });
        }
        if !height_m.is_finite() {
            return Err(Error::InvalidObserver {
                field: "height_m",
                value: height_m,
            });
        }
        Ok(Self {
            lon_deg,
            lat_deg,
            height_m,
        })
    }

    /// East longitude in degrees, `[-180, 360]`.
    pub fn lon_deg(self) -> f64 {
        self.lon_deg
    }

    /// Latitude in degrees, `[-90, 90]`.
    pub fn lat_deg(self) -> f64 {
        self.lat_deg
    }

    /// Ellipsoidal height in metres (finite; negative values are allowed).
    pub fn height_m(self) -> f64 {
        self.height_m
    }
}

/// Sidereal configuration: tropical or a selected ayanamsha mode.
///
/// The [`Sidereal`](SiderealConfig::Sidereal) variant carries the full
/// sidereal-selection registry selection (canonical `SIDM_*` base ID, `SIDM_USER`
/// epoch pair, retained `SIDBIT` flags); see [`crate::sidereal`]. The
/// default [`Tropical`](SiderealConfig::Tropical) reads back exactly like
/// the release never-set state (Fagan/Bradley ID, zero epoch pair, no
/// flags), and selecting it again clears any previous selection.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub enum SiderealConfig {
    /// Tropical zodiac (no ayanamsha shift).
    #[default]
    Tropical,
    /// Selected sidereal mode (ayanamsha system for `FLG_SIDEREAL` work).
    Sidereal(SiderealSelection),
}

/// Stored filesystem locations for provider assets.
///
/// The value itself probes no paths. Context path setters retire resources
/// selected from a changed path; actual acquisition checks file existence,
/// structure, identity and coverage according to the selected operation.
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct DataPaths {
    /// Directory holding LEB/LEB2 assets, if configured.
    pub leb_dir: Option<PathBuf>,
    /// Directory holding JPL/SPK kernels, if configured.
    pub jpl_dir: Option<PathBuf>,
}

/// Owned calculation settings with retained, shareable data resources.
///
/// Construct with [`Self::new`]. Consuming `with_*` builders return an updated
/// value; `set_*` methods update one in place. Clones have independent settings
/// and retain shared immutable image/cache owners. Resetting or closing one
/// context does not invalidate images held by an existing snapshot.
///
/// Optional configuration fields distinguish an explicit value from an unset
/// override. A [`crate::config::Runtime`] resolves unset fields against captured
/// environment, dotenv, TOML and defaults. A bare context does not discover
/// those inputs automatically: callers may instead supply values through the
/// `resolve_*` functions. Concrete mode/tier fields have their own defaults.
/// `None` therefore removes an override; it is not uniformly equivalent to
/// `false`, zero or a particular provider.
///
/// Observer angles are degrees, heights metres, model Delta-T days, tidal
/// acceleration arcseconds per Julian century squared and lapse rate K/m.
/// Numeric store-only setters may accept nonfinite values: a context is not a
/// proof that every setting or requested calculation is valid. Checked methods
/// and consuming operations document their own errors. Context construction
/// performs no file/network I/O or astronomical calculation.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct EphemerisContext {
    observer: Option<Observer>,
    sidereal: SiderealConfig,
    tier: Tier,
    mode: ProviderMode,
    paths: DataPaths,
    angles_cache: HashMap<String, f64>,
    user_delta_t_days: Option<f64>,
    tidal_acceleration: Option<f64>,
    delta_t_model: Option<DeltaTModel>,
    iers_delta_t: Option<bool>,
    lapse_rate: Option<f64>,
    network_policy: Option<NetworkPolicy>,
    /// Explicit LEB file selected for binary-ephemeris work: `None` (the default) disables the explicit-file
    /// leg and leaves source selection to the mode/tier/discovery chain.
    leb_file: Option<PathBuf>,
    /// Optional override for future file acquisition; attached images remain owned
    /// by their original storage policy and active snapshots stay valid.
    leb_read_policy: Option<crate::leb::LebReadPolicy>,
    /// Explicit JPL ephemeris filename: `Some`
    /// takes priority over the tier kernel; `None` (the default)
    /// resumes tier-based selection. Cleared by a tier change and by
    /// [`EphemerisContext::close`], exactly like the release.
    ephemeris_file: Option<String>,
    /// Already resolved main-image loader directory. No environment, home or
    /// TOML lookup occurs in the core; callers supply the filesystem root.
    /// Preserved across reset/close like the source process data-dir setting.
    jpl_data_dir: Option<PathBuf>,
    /// Programmatic SPK cache-directory override: `None` (the default) defers to the
    /// environment/TOML/default chain ([`crate::spk::resolve_spk_cache_dir`]).
    /// A plain store with no invalidation, exactly like the release.
    spk_cache_dir: Option<PathBuf>,
    /// Programmatic automatic-SPK-download toggle: `None` (the default)
    /// defers to the environment/TOML/default chain
    /// ([`resolve_auto_spk_download`]). Reset to `None` by
    /// [`EphemerisContext::close`], exactly like the release.
    auto_spk_download: Option<bool>,
    /// SPK download date padding in days on each side of the requested
    /// range: `0` (the default)
    /// means no padding. Stored as `u64`; negativity is refused at the
    /// consuming boundary, mirroring the release `ValueError`. Reset to
    /// `0` by [`EphemerisContext::close`], exactly like the release.
    spk_date_padding_days: u64,
    /// Programmatic strict-precision toggle: `None` (the default) defers to the
    /// environment/TOML/default chain ([`resolve_strict_precision`]).
    /// Preserved by [`EphemerisContext::close`], exactly like the release
    strict_precision: Option<bool>,
    /// Programmatic IERS cache-directory override: `None` (the default)
    /// uses the release default cache location. A plain store with no
    /// validation and no directory creation, exactly like the release
    /// (creation happens only when an enabled lookup consumes the lazy
    /// acquisition). Preserved by [`EphemerisContext::close`],
    /// exactly like the release.
    iers_cache_dir: Option<PathBuf>,
    /// Lazily acquired immutable IERS rows; configuration stores do not drop rows.
    iers_source: crate::iers::IersContextSource,
    /// Programmatic automatic-IERS-download toggle: `None` (the default)
    /// defers to the environment/TOML/default chain
    /// ([`resolve_iers_auto_download`]). Preserved by
    /// [`EphemerisContext::close`], exactly like the release (the
    /// release `_close_inner` never resets `_IERS_AUTO_DOWNLOAD`).
    iers_auto_download: Option<bool>,
    /// Owned LEB images bound to this context (path-configuration layer).
    ///
    /// Loaded from the configured [`DataPaths::leb_dir`], from an
    /// explicit file, or from memory via
    /// [`LebFile::from_bytes`](crate::leb::LebFile::from_bytes);
    /// composite/tiered indexes borrow these images. A changed LEB
    /// path drops every entry ([`EphemerisContext::set_leb_dir`]);
    /// [`EphemerisContext::reset_session`] preserves them and
    /// [`EphemerisContext::close`] releases them.
    leb_files: Vec<LebFile>,
    /// Reusable dispatch; clones share it until images, tier, mode or centers change.
    leb_provider: crate::positions::context_serving::LebProviderCache,
    /// Context-local SPK kernel store and body map (path-configuration layer).
    ///
    /// The registry is context-independent by design (registry-lifetime layer): each
    /// context owns its kernels, so no kernel is ever shared across
    /// contexts through a global. A changed JPL path clears it
    /// ([`EphemerisContext::set_jpl_dir`]);
    /// [`EphemerisContext::reset_session`] preserves it and
    /// [`EphemerisContext::close`] releases it.
    spk_registry: SpkRegistry,
    /// Optional shared Horizons transport/cache owner. Context clones retain
    /// this same Arc; changing a clone's settings does not duplicate the cache.
    /// Acquiring a client and deciding whether a cache miss may use the network
    /// belong to the consuming route. Attaching it performs no network I/O.
    /// Context reset preserves it; close installs an empty independent client
    /// with the same configuration while retained owners remain valid.
    horizons: Option<std::sync::Arc<HorizonsClient>>,
}

impl EphemerisContext {
    /// Create a primary API builder from this explicitly configured native context.
    ///
    /// # Returns
    ///
    /// A builder sharing this context's immutable data owners. No environment
    /// configuration is resolved and no resource is acquired by this adapter.
    pub fn primary_session_builder(&self) -> crate::api::SessionBuilder {
        crate::api::SessionBuilder::from_context(self.clone())
    }

    /// Explicit storage-policy override, or `None` for captured configuration.
    /// This accessor performs no I/O and does not resolve environment values.
    pub fn leb_read_policy_override(&self) -> Option<crate::leb::LebReadPolicy> {
        self.leb_read_policy
    }

    /// Effective native file policy; unresolved contexts default to owned.
    /// Runtime configuration resolution may select an indexed override.
    pub fn leb_read_policy(&self) -> crate::leb::LebReadPolicy {
        self.leb_read_policy.unwrap_or_default()
    }

    /// Select the storage of future LEB opens. Invalid limits leave all state
    /// unchanged. Explicitly attached files and existing snapshots are retained;
    /// no inode is reopened and no environment lookup occurs here.
    pub fn set_leb_read_policy(
        &mut self,
        policy: Option<crate::leb::LebReadPolicy>,
    ) -> Result<bool> {
        if let Some(policy) = policy {
            policy.validate()?;
        }
        let changed = self.leb_read_policy != policy;
        self.leb_read_policy = policy;
        Ok(changed)
    }

    /// Open under the resolved policy, sharing the first compatible attached
    /// indexed image's resource domain. No process-global store is involved.
    pub(crate) fn open_leb_file(&self, path: &std::path::Path) -> Result<LebFile> {
        if let crate::leb::LebReadPolicy::Indexed(limits) = self.leb_read_policy() {
            for file in &self.leb_files {
                if let crate::leb::LebSource::Positional(source) = file.source()? {
                    if source.limits == limits {
                        return LebFile::open_indexed_with_resources(
                            path,
                            std::sync::Arc::clone(&source.resources),
                        );
                    }
                }
            }
        }
        self.leb_read_policy().open(path)
    }

    /// Default context: no observer (geocentric), tropical, medium tier
    /// routed `auto` mode (the route-planning layer plans LEB → Horizons → native-JPL from the
    /// bound stores), no data paths, and every Delta-T/refraction
    /// override unset (automatic).
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the observer (replaces any previous one).
    pub fn with_observer(mut self, observer: Observer) -> Self {
        self.observer = Some(observer);
        self
    }

    /// Clear the observer (geocentric calculations).
    pub fn without_observer(mut self) -> Self {
        self.observer = None;
        self
    }

    /// Set the sidereal configuration.
    pub fn with_sidereal(mut self, sidereal: SiderealConfig) -> Self {
        self.sidereal = sidereal;
        self
    }

    /// Set the data tier and clear the explicit JPL filename, even if unchanged.
    ///
    /// Tier-based kernel selection resumes until a caller/configuration resolver
    /// supplies another filename. An explicit LEB file is preserved. The tier
    /// does not extend a body's actual stored/model coverage; no file is opened.
    pub fn with_tier(mut self, tier: Tier) -> Self {
        self.leb_provider = Default::default();
        self.tier = tier;
        self.ephemeris_file = None;
        self
    }

    /// Set a tier by exact name: `base`, `medium` or `extended`.
    ///
    /// Like [`Self::with_tier`], every successful call clears the explicit JPL
    /// filename, even for the current tier, while preserving the LEB selection.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidTier`] for an unknown spelling. The in-place
    /// counterpart [`Self::set_tier_name`] preserves state on this error.
    pub fn with_tier_name(mut self, name: &str) -> Result<Self> {
        self.tier = Tier::from_name(name).ok_or_else(|| Error::InvalidTier {
            value: name.to_string(),
        })?;
        self.leb_provider = Default::default();
        self.ephemeris_file = None;
        Ok(self)
    }

    /// Set the provider mode.
    pub fn with_mode(mut self, mode: ProviderMode) -> Self {
        self.leb_provider = Default::default();
        self.mode = mode;
        self.refresh_iers_source();
        self
    }

    /// Set the provider mode by release calc-mode name (normalized
    /// Unicode lowercase/Unicode-whitespace-insensitive: `"auto"` / `"skyfield"` / `"leb"`
    /// / `"horizons"`). Normalization uses Unicode lowercase and source strip; unknown spellings are [`Error::InvalidCalcMode`]
    /// naming `set_calc_mode` as the source, mirroring the release
    /// `ValueError`.
    pub fn with_mode_name(mut self, name: &str) -> Result<Self> {
        self.mode =
            ProviderMode::from_calc_mode_name(name).ok_or_else(|| Error::InvalidCalcMode {
                value: name.to_string(),
                source: CALC_MODE_SOURCE_SETTER,
            })?;
        self.leb_provider = Default::default();
        self.refresh_iers_source();
        Ok(self)
    }

    /// Attach a Horizons network configuration (consuming builder for
    /// construction).
    ///
    /// Materializes the context-owned [`HorizonsClient`] eagerly from
    /// `config`; the retry budget stays fixed and construction
    /// performs no I/O. A refused timeout ([`Error::InvalidTimeout`])
    /// aborts the builder chain with no partial attachment. Replaces any
    /// previously attached client (and its cached vectors).
    pub fn with_horizons_config(mut self, config: HorizonsConfig) -> Result<Self> {
        self.horizons = Some(std::sync::Arc::new(HorizonsClient::from_config(config)?));
        Ok(self)
    }

    /// Retain one explicitly shared Horizons vector cache and transport owner.
    /// Cloned calculation snapshots retain this same owner. Each cache miss
    /// still obeys the consuming context's network policy; no I/O occurs here.
    pub fn with_shared_horizons_client(mut self, client: std::sync::Arc<HorizonsClient>) -> Self {
        self.horizons = Some(client);
        self
    }

    /// Remove this context's retained Horizons client without network I/O.
    ///
    /// Other clones and active users retain their existing Arc owners. This
    /// does not select a provider mode or clear another owner's vector cache;
    /// a consuming route decides whether and how to acquire a later client.
    pub fn without_horizons_client(mut self) -> Self {
        self.horizons = None;
        self
    }

    /// Borrow the retained Horizons transport/cache, if attached.
    ///
    /// `None` means only that this context owns no client; the consuming route
    /// decides any later acquisition. Borrowing performs no I/O. Context clones
    /// share an attached owner, while network misses retain the consuming
    /// route's effective-policy checks.
    pub fn horizons_client(&self) -> Option<&HorizonsClient> {
        self.horizons.as_deref()
    }

    /// Retain the current transport/cache owner without copying its entries.
    pub(crate) fn shared_horizons_client(&self) -> Option<std::sync::Arc<HorizonsClient>> {
        self.horizons.as_ref().map(std::sync::Arc::clone)
    }

    /// Set the data paths (consuming builder for construction).
    ///
    /// When the context already holds bound stores, this runs the
    /// path-configuration invalidation rule: an unchanged path keeps its stores
    /// warm, a changed path drops exactly the stores resolved from it
    /// (see [`EphemerisContext::set_paths`]).
    pub fn with_paths(mut self, paths: DataPaths) -> Self {
        self.set_paths(paths);
        self
    }

    /// Select the explicit LEB file (consuming builder for construction).
    ///
    /// Same contract as [`EphemerisContext::set_leb_file`].
    pub fn with_leb_file(mut self, file: Option<PathBuf>) -> Self {
        self.set_leb_file(file);
        self
    }

    /// Select the explicit JPL ephemeris filename (consuming builder for
    /// construction).
    ///
    /// Same contract as [`EphemerisContext::set_ephemeris_file`].
    pub fn with_ephemeris_file(mut self, name: Option<&str>) -> Self {
        self.set_ephemeris_file(name);
        self
    }

    /// Bind an already resolved directory for main JPL image acquisition.
    /// Same store-only and lifecycle contract as [`Self::set_jpl_data_dir`].
    pub fn with_jpl_data_dir(mut self, directory: Option<PathBuf>) -> Self {
        self.set_jpl_data_dir(directory);
        self
    }

    /// Set the programmatic SPK cache-directory override (consuming
    /// builder for construction).
    ///
    /// Same contract as [`EphemerisContext::set_spk_cache_dir`].
    pub fn with_spk_cache_dir(mut self, dir: Option<PathBuf>) -> Self {
        self.set_spk_cache_dir(dir);
        self
    }

    /// Store a model TT−UT1 override in **days**, not seconds.
    ///
    /// `None` and `Some(DELTAT_AUTOMATIC)` resume computed values. Other
    /// values, including nonfinite ones, are stored without validation and
    /// returned by the model Delta-T path without computing. Consuming
    /// calculations retain their own domain checks. Independent frame and
    /// observer clocks are not replaced by this override.
    pub fn with_user_delta_t(mut self, days: Option<f64>) -> Self {
        self.user_delta_t_days = match days {
            None | Some(DELTAT_AUTOMATIC) => None,
            Some(v) => Some(v),
        };
        self
    }

    /// Store lunar tidal acceleration in arcseconds per Julian century squared.
    /// `None` and `TIDAL_AUTOMATIC` (999999) resume automatic selection.
    /// Other values are stored without validation; consuming operations retain
    /// their own numerical domains.
    pub fn with_tidal_acceleration(mut self, acc: Option<f64>) -> Self {
        self.tidal_acceleration = match acc {
            None => None,
            Some(v) if v == TIDAL_AUTOMATIC as f64 => None,
            Some(v) => Some(v),
        };
        self
    }

    /// Select the Delta-T model programmatically: `None` defers to the environment/TOML/default chain
    /// ([`resolve_delta_t_model`]).
    pub fn with_delta_t_model(mut self, model: Option<DeltaTModel>) -> Self {
        self.delta_t_model = model;
        self
    }

    /// Select the Delta-T model by release name (`None` defers, like the
    /// release `None`). Names normalize ASCII case/whitespace-insensitively;
    /// unknown names are [`Error::InvalidDeltaTModel`], mirroring the
    /// release `ValueError`.
    pub fn with_delta_t_model_name(mut self, name: Option<&str>) -> Result<Self> {
        self.delta_t_model = match name {
            None => None,
            Some(n) => {
                Some(
                    DeltaTModel::from_name(n).ok_or_else(|| Error::InvalidDeltaTModel {
                        value: n.to_string(),
                    })?,
                )
            }
        };
        Ok(self)
    }

    /// Enable or disable observed-IERS Delta T: `None` defers to the
    /// environment/TOML/default chain ([`resolve_iers_enabled`]).
    pub fn with_iers_delta_t(mut self, enabled: Option<bool>) -> Self {
        self.iers_delta_t = enabled;
        self
    }

    /// Enable or disable automatic SPK download for minor bodies: `None` defers to the
    /// environment/TOML/default chain ([`resolve_auto_spk_download`]).
    pub fn with_auto_spk_download(mut self, enabled: Option<bool>) -> Self {
        self.auto_spk_download = enabled;
        self
    }

    /// Set the SPK download date padding in days on each side of the
    /// requested range.
    pub fn with_spk_date_padding_days(mut self, days: u64) -> Self {
        self.spk_date_padding_days = days;
        self
    }

    /// Enable or disable strict precision mode for major asteroids
    /// `None` defers to the
    /// environment/TOML/default chain ([`resolve_strict_precision`]).
    pub fn with_strict_precision(mut self, enabled: Option<bool>) -> Self {
        self.strict_precision = enabled;
        self
    }

    /// Set the programmatic IERS cache-directory override (consuming
    /// builder for construction).
    ///
    /// Same contract as [`EphemerisContext::set_iers_cache_dir`].
    pub fn with_iers_cache_dir(mut self, dir: Option<PathBuf>) -> Self {
        self.set_iers_cache_dir(dir);
        self
    }

    /// Enable or disable automatic IERS data download: `None` defers to the
    /// environment/TOML/default chain ([`resolve_iers_auto_download`]).
    pub fn with_iers_auto_download(mut self, enabled: Option<bool>) -> Self {
        self.iers_auto_download = enabled;
        self.refresh_iers_source();
        self
    }

    /// Store the atmospheric lapse rate in K/m.
    /// `None` resumes [`LAPSE_RATE_DEFAULT`]. Values are stored without
    /// validation; consuming atmospheric operations retain their own domains.
    pub fn with_lapse_rate(mut self, rate: Option<f64>) -> Self {
        self.lapse_rate = rate;
        self
    }

    /// Set the programmatic network-policy override: `Some` policy wins over every configured source in
    /// [`resolve_network_policy`]; `None` (the default) defers to the
    /// environment/TOML/default chain. Preserved by both
    /// [`EphemerisContext::reset_session`] and
    /// [`EphemerisContext::close`]: the release process global survives
    /// both calls, so the override is process-level configuration, not
    /// session state.
    pub fn with_network_policy(mut self, policy: Option<NetworkPolicy>) -> Self {
        self.network_policy = policy;
        self.refresh_iers_source();
        self
    }

    /// Set the programmatic network-policy override by release policy name
    /// (`"auto"` / `"allow"` / `"sealed"`, ASCII case/whitespace
    /// normalized). Unknown spellings are [`Error::InvalidNetworkPolicy`]
    /// naming `set_network_policy` as the source, mirroring the release
    /// `ValueError`. `None` clears the override back to the
    /// environment/TOML/default chain.
    pub fn with_network_policy_name(mut self, name: Option<&str>) -> Result<Self> {
        self.network_policy = match name {
            None => None,
            Some(raw) => Some(NetworkPolicy::from_policy_name(raw).ok_or_else(|| {
                Error::InvalidNetworkPolicy {
                    value: raw.to_string(),
                    source: NETWORK_POLICY_SOURCE_SETTER,
                }
            })?),
        };
        self.refresh_iers_source();
        Ok(self)
    }

    /// Programmatic network-policy override, if set.
    pub fn network_policy_override(&self) -> Option<NetworkPolicy> {
        self.network_policy
    }

    // ------------------------------------------------------------------
    // In-place setters for handle-backed contexts (C-ABI serving).
    //
    // The builder (`with_*`) methods above consume `self`; the C-ABI
    // context handle owns a single boxed `EphemerisContext` whose
    // identity must survive configuration, so these route the same
    // validated logic through `&mut self`. Each one clones, applies
    // the corresponding builder, and assigns back only on success — a
    // refused name leaves the stored value untouched (a take-then-`?`
    // formulation would clobber state with the default on the error
    // path). Contracts, normalization and errors are the builder's.
    // ------------------------------------------------------------------

    /// Set the observer in place (replaces any previous one).
    /// Same contract as [`EphemerisContext::with_observer`].
    pub fn set_observer(&mut self, observer: Observer) {
        *self = self.clone().with_observer(observer);
    }

    /// Clear the observer in place (geocentric calculations).
    /// Same contract as [`EphemerisContext::without_observer`].
    pub fn clear_observer(&mut self) {
        *self = self.clone().without_observer();
    }

    /// Set the data tier by release name in place.
    /// Same contract as [`EphemerisContext::with_tier_name`].
    pub fn set_tier_name(&mut self, name: &str) -> Result<()> {
        *self = self.clone().with_tier_name(name)?;
        Ok(())
    }

    /// Set the provider mode by release calc-mode name in place.
    /// Same contract as [`EphemerisContext::with_mode_name`].
    pub fn set_mode_name(&mut self, name: &str) -> Result<()> {
        *self = self.clone().with_mode_name(name)?;
        Ok(())
    }

    /// Select the Delta-T model by release name in place (`None`
    /// defers, like the release `None`).
    /// Same contract as [`EphemerisContext::with_delta_t_model_name`].
    pub fn set_delta_t_model_name(&mut self, name: Option<&str>) -> Result<()> {
        *self = self.clone().with_delta_t_model_name(name)?;
        Ok(())
    }

    /// Set the user Delta-T override in days in place (`None` clears;
    /// the `DELTAT_AUTOMATIC` sentinel clears like the release).
    /// Same contract as [`EphemerisContext::with_user_delta_t`].
    pub fn set_user_delta_t(&mut self, days: Option<f64>) {
        *self = self.clone().with_user_delta_t(days);
    }

    /// Enable or disable observed-IERS Delta T in place (`None`
    /// defers).
    /// Same contract as [`EphemerisContext::with_iers_delta_t`].
    pub fn set_iers_delta_t(&mut self, enabled: Option<bool>) {
        *self = self.clone().with_iers_delta_t(enabled);
    }

    /// Enable or disable automatic SPK download in place (`None`
    /// defers).
    /// Same contract as [`EphemerisContext::with_auto_spk_download`].
    pub fn set_auto_spk_download(&mut self, enabled: Option<bool>) {
        *self = self.clone().with_auto_spk_download(enabled);
    }

    /// Set the SPK download date padding in days in place.
    /// Same contract as
    /// [`EphemerisContext::with_spk_date_padding_days`].
    pub fn set_spk_date_padding_days(&mut self, days: u64) {
        *self = self.clone().with_spk_date_padding_days(days);
    }

    /// Enable or disable strict precision mode in place (`None`
    /// defers).
    /// Same contract as [`EphemerisContext::with_strict_precision`].
    pub fn set_strict_precision(&mut self, enabled: Option<bool>) {
        *self = self.clone().with_strict_precision(enabled);
    }

    /// Set the programmatic IERS cache-directory override in place
    /// (`None` reverts to the release default cache location).
    ///
    /// Release `set_iers_cache_dir` semantics (context-local instead of
    /// the process-wide `_IERS_CACHE_DIR` global): a plain store with no
    /// validation and no directory creation — creation happens at first
    /// cache use in the filesystem leg, not here (the store-not-probe
    /// rule). Infallible on a live context.
    pub fn set_iers_cache_dir(&mut self, dir: Option<PathBuf>) {
        self.iers_cache_dir = dir;
        self.refresh_iers_source();
    }

    /// Enable or disable automatic IERS data download in place (`None`
    /// defers).
    /// Same contract as [`EphemerisContext::with_iers_auto_download`].
    pub fn set_iers_auto_download(&mut self, enabled: Option<bool>) {
        *self = self.clone().with_iers_auto_download(enabled);
    }

    /// Set the tidal-acceleration override in place (`None` resumes
    /// automatic mode, like the `TIDAL_AUTOMATIC` sentinel).
    /// Same contract as [`EphemerisContext::with_tidal_acceleration`].
    pub fn set_tidal_acceleration(&mut self, acc: Option<f64>) {
        *self = self.clone().with_tidal_acceleration(acc);
    }

    /// Set the atmospheric lapse-rate override in place (`None` resets
    /// to the default).
    /// Same contract as [`EphemerisContext::with_lapse_rate`].
    pub fn set_lapse_rate(&mut self, rate: Option<f64>) {
        *self = self.clone().with_lapse_rate(rate);
    }

    /// Set the programmatic network-policy override in place (`None`
    /// defers to the environment/TOML/default chain).
    /// Same contract as [`EphemerisContext::with_network_policy`].
    pub fn set_network_policy(&mut self, policy: Option<NetworkPolicy>) {
        *self = self.clone().with_network_policy(policy);
    }

    /// Set the programmatic network-policy override by release name in
    /// place (`None` clears the override; a refused spelling leaves the
    /// stored value untouched).
    /// Same contract as [`EphemerisContext::with_network_policy_name`].
    pub fn set_network_policy_name(&mut self, name: Option<&str>) -> Result<()> {
        *self = self.clone().with_network_policy_name(name)?;
        Ok(())
    }

    /// Select the explicit LEB file in place (`None` disables the
    /// explicit-file leg).
    ///
    /// Release `set_leb_file` file-leg semantics (context-local instead
    /// of the process-wide `_LEB_FILE` global): a changed selection
    /// stores the new value and drops every bound LEB image, forcing
    /// re-resolution on the next calculation — the engine analogue of
    /// closing the old reader and resetting the active-reader/frame
    /// caches. Follows the path-configuration idempotency rule like
    /// [`EphemerisContext::set_leb_dir`]: an unchanged selection is a
    /// no-op returning `false` and keeps every bound image warm (the
    /// release always closes and re-creates its reader, even unchanged;
    /// the re-created reader serves identical data, so the no-op is
    /// semantically behavior-preserving). Returns `true` when the stored
    /// selection changed. Comparison is exact on the stored [`PathBuf`]
    /// (no canonicalization, no filesystem probing — the store-not-probe
    /// rule); directory scanning for companion assets stays a caller-side
    /// operation and never runs here.
    pub fn set_leb_file(&mut self, file: Option<PathBuf>) -> bool {
        if self.leb_file == file {
            return false;
        }
        self.leb_file = file;
        self.release_leb_images();
        true
    }

    /// Select the explicit JPL ephemeris filename in place (`None`
    /// clears the override and resumes tier-based selection).
    ///
    /// Release `set_ephemeris_file` / `set_jpl_file` semantics
    /// (context-local instead of the process-wide `_EPHEMERIS_FILE` +
    /// `_EPHEMERIS_FILE_EXPLICIT` globals): any name, including an
    /// unrecognized one, is stored as given with no validation —
    /// recognition, local presence, and download decisions belong to
    /// filename resolution ([`crate::positions::select_de_kernel`],
    /// [`crate::positions::decide_jpl_file`]), not to configuration.
    /// Like the release (which closes only the DE kernel and clears
    /// computation caches, leaving SPK body registrations intact), this
    /// drops no registry state: the next load re-resolves the kernel
    /// from the stored override. Infallible on a live context.
    pub fn set_ephemeris_file(&mut self, name: Option<&str>) {
        self.ephemeris_file = name.map(str::to_string);
    }

    /// Set the programmatic SPK cache-directory override in place
    /// (`None` reverts to the environment/TOML/default chain).
    ///
    /// Release `set_spk_cache_dir` semantics (context-local instead of
    /// the process-wide `_SPK_CACHE_DIR` global): a plain store with no
    /// validation and no invalidation — already-registered SPK bodies
    /// are unaffected, exactly like the release. Directory creation
    /// happens at first cache use, not here (the store-not-probe rule).
    /// Infallible on a live context.
    pub fn set_spk_cache_dir(&mut self, dir: Option<PathBuf>) {
        self.spk_cache_dir = dir;
    }

    /// Effective network gate for this context's provider mode: the
    /// override when set, else the release `"auto"` default, resolved
    /// against [`EphemerisContext::mode`] (see
    /// [`effective_network_policy`]). Environment/TOML inputs arrive via
    /// [`resolve_network_policy`] on the binding side.
    pub fn effective_network_policy(&self) -> EffectiveNetworkPolicy {
        effective_network_policy(
            self.network_policy.unwrap_or(NetworkPolicy::Auto),
            self.mode,
        )
    }

    /// Refuse network I/O before any socket work when this context is
    /// sealed: the NET provider call sites
    /// pass their operation label as `purpose`.
    pub fn require_network(&self, purpose: &str) -> Result<()> {
        require_network(self.effective_network_policy(), purpose)
    }

    /// Current observer, if any.
    pub fn observer(&self) -> Option<Observer> {
        self.observer
    }

    /// Current sidereal configuration.
    pub fn sidereal(&self) -> SiderealConfig {
        self.sidereal
    }

    /// Select a sidereal (ayanamsha) mode and retain its normalization warnings.
    ///
    /// `mode` combines a `SIDM_*` identifier with supported `SIDBIT_*` flags;
    /// see [`SiderealSelection::set_sidereal_mode`]. Implemented flags are
    /// retained. Unknown high bits are stripped with
    /// [`Warning::UnsupportedSiderealProjection`]; negative IDs remain stored
    /// for the calculation's invalid-mode fallback.
    ///
    /// USER mode reads `t0` as a TT Julian day, or UT1 with `SIDBIT_USER_UT`,
    /// and `ayan_t0` as degrees. Both are stored without finite validation.
    /// Returns warnings rather than an error; the selection is always stored.
    /// This does not itself apply a sidereal correction to a position.
    pub fn set_sidereal_mode(&mut self, mode: i32, t0: f64, ayan_t0: f64) -> Vec<Warning> {
        let mut sel = match self.sidereal {
            SiderealConfig::Sidereal(sel) => sel,
            SiderealConfig::Tropical => SiderealSelection::unset(),
        };
        let warnings = sel.set_sidereal_mode(mode, t0, ayan_t0);
        self.sidereal = SiderealConfig::Sidereal(sel);
        warnings
    }

    /// Store implemented sidereal flags without changing the base or epoch pair.
    ///
    /// `bits` selects projection and defining-epoch timescale flags through
    /// [`crate::sidereal::IMPLEMENTED_SIDBITS`]. Unknown bits are discarded
    /// silently; this operation emits no warnings and validates no mode or
    /// date. The stored mode may be negative. The defining epoch is a Julian
    /// Day and its offset is in degrees; both retain their original IEEE words.
    /// A tropical context supplies the default base and zero epoch pair.
    ///
    /// This field update supports warning delivery between the flags and base
    /// assignments. Immutable Sessions built earlier retain their own settings.
    pub fn set_sidereal_bits(&mut self, bits: i32) {
        let mut selection = match self.sidereal {
            SiderealConfig::Sidereal(selection) => selection,
            SiderealConfig::Tropical => SiderealSelection::unset(),
        };
        selection.set_retained_bits(bits);
        self.sidereal = SiderealConfig::Sidereal(selection);
    }

    /// Configure independently stored sidereal flags on an owned context.
    ///
    /// Preserves the base and epoch pair and applies the same mask and
    /// warning-free behavior as [`EphemerisContext::set_sidereal_bits`].
    pub fn with_sidereal_bits(mut self, bits: i32) -> Self {
        self.set_sidereal_bits(bits);
        self
    }

    /// Stored sidereal base mode ID: the
    /// Fagan/Bradley default when never set or tropical.
    pub fn sidereal_mode(&self) -> i32 {
        match self.sidereal {
            SiderealConfig::Sidereal(sel) => sel.mode(),
            SiderealConfig::Tropical => SiderealSelection::unset().mode(),
        }
    }

    /// Full `(mode, t0, ayan_t0)` triple: USER epoch is TT JD unless
    /// `SIDBIT_USER_UT` selects UT1, and the defining offset is degrees.
    pub fn sidereal_mode_full(&self) -> (i32, f64, f64) {
        match self.sidereal {
            SiderealConfig::Sidereal(sel) => sel.full(),
            SiderealConfig::Tropical => SiderealSelection::unset().full(),
        }
    }

    /// Retained implemented `SIDBIT` flags; zero when never set or tropical.
    pub fn sidereal_bits(&self) -> i32 {
        match self.sidereal {
            SiderealConfig::Sidereal(sel) => sel.bits(),
            SiderealConfig::Tropical => 0,
        }
    }

    /// Current data tier.
    pub fn tier(&self) -> Tier {
        self.tier
    }

    /// Current provider mode.
    pub fn mode(&self) -> ProviderMode {
        self.mode
    }

    /// Current data paths.
    pub fn paths(&self) -> &DataPaths {
        &self.paths
    }

    /// Explicit LEB file selection, if set.
    pub fn leb_file(&self) -> Option<&PathBuf> {
        self.leb_file.as_ref()
    }

    /// Explicit JPL ephemeris filename override, if set. Feeds the
    /// `explicit` leg of [`crate::positions::select_de_kernel`]; `None`
    /// while tier-based selection applies.
    pub fn ephemeris_file_override(&self) -> Option<&str> {
        self.ephemeris_file.as_deref()
    }

    /// Explicitly resolved main-image loader directory, if supplied.
    /// Custom [`DataPaths::jpl_dir`] candidates take priority. Absolute
    /// filenames override both prefixes. An empty native registry can read a
    /// present local file, or acquire a genuinely missing selected kernel
    /// under this context's network policy, before a public clock is formed.
    pub fn jpl_data_dir(&self) -> Option<&PathBuf> {
        self.jpl_data_dir.as_ref()
    }

    /// Store the caller-resolved main-image loader directory without I/O.
    ///
    /// The boundary caller resolves environment/TOML/home inputs before this setter;
    /// the SDK receives an explicit path. No canonicalization, environment
    /// read, directory creation or resource invalidation occurs here. Existing
    /// owned registry images retain their identity; only a future empty-registry
    /// acquisition consumes the directory. Reset and close preserve it, as the
    /// source preserves process-level data-directory configuration.
    /// An unset directory permits absolute local filenames; a permitted
    /// download needing a relative destination refuses with missing config.
    pub fn set_jpl_data_dir(&mut self, directory: Option<PathBuf>) {
        self.jpl_data_dir = directory;
        self.refresh_iers_source();
    }

    /// Programmatic SPK cache-directory override, if set. Feeds the
    /// `explicit` leg of [`crate::spk::resolve_spk_cache_dir`]; `None`
    /// while the environment/TOML/default chain applies.
    pub fn spk_cache_dir(&self) -> Option<&PathBuf> {
        self.spk_cache_dir.as_ref()
    }

    /// Programmatic IERS cache-directory override, if set. `None` while
    /// the release default cache location applies (the consumed acquisition layer resolves the default directory).
    pub fn iers_cache_dir(&self) -> Option<&PathBuf> {
        self.iers_cache_dir.as_ref()
    }

    /// Set the LEB asset directory without probing or creating it.
    ///
    /// Repeating the same directory retains attached images; changing the
    /// directory retires cached readers so subsequent
    /// calculations resolve the new selection. Active owners remain valid.
    pub fn set_leb_dir(&mut self, dir: Option<PathBuf>) -> bool {
        if self.paths.leb_dir == dir {
            return false;
        }
        self.paths.leb_dir = dir;
        self.release_leb_images();
        true
    }

    /// Set the JPL/SPK kernel directory (path-configuration layer).
    ///
    /// Same contract as [`EphemerisContext::set_leb_dir`]: an unchanged
    /// path is a no-op returning `false`; a changed path stores the new
    /// value, clears the context-local [`SpkRegistry`]
    /// (kernels and body registrations resolved from the old directory),
    /// and returns `true`.
    pub fn set_jpl_dir(&mut self, dir: Option<PathBuf>) -> bool {
        if self.paths.jpl_dir == dir {
            return false;
        }
        self.paths.jpl_dir = dir;
        self.leb_provider = Default::default();
        self.spk_registry.clear();
        true
    }

    /// Set both data paths through the path-configuration invalidation rule.
    ///
    /// Returns `true` when either path changed (and exactly the stores
    /// resolved from a changed path were dropped); an unchanged pair is
    /// a no-op returning `false`.
    pub fn set_paths(&mut self, paths: DataPaths) -> bool {
        let leb_changed = self.set_leb_dir(paths.leb_dir);
        let jpl_changed = self.set_jpl_dir(paths.jpl_dir);
        leb_changed || jpl_changed
    }

    /// Bind an owned LEB image to this context.
    ///
    /// The image typically comes from the configured LEB directory, an
    /// explicit file ([`LebFile::open`](crate::leb::LebFile::open)), or
    /// memory ([`LebFile::from_bytes`](crate::leb::LebFile::from_bytes)).
    /// Composite/tiered indexes share immutable image owners. Bound images survive
    /// [`EphemerisContext::reset_session`]; a changed LEB path or
    /// [`EphemerisContext::close`] releases them.
    pub fn attach_leb_file(&mut self, file: LebFile) {
        self.leb_provider = Default::default();
        self.leb_files.push(file);
    }

    /// Owned LEB images bound to this context, in attach order.
    pub fn leb_files(&self) -> &[LebFile] {
        &self.leb_files
    }

    /// Release this context's loaded LEB images while preserving selection.
    ///
    /// Configured file/directory paths, provider mode, observer, models and
    /// other caches remain unchanged. The derived provider is retired. Previously
    /// cloned contexts and active Sessions retain their immutable image owners;
    /// this operation does not close those shared handles or revoke a query.
    /// A later acquisition can reopen the same configured paths.
    pub fn release_leb_images(&mut self) {
        self.leb_provider = Default::default();
        self.leb_files.clear();
    }

    /// Borrow the image/tier dispatch shared by this immutable context snapshot.
    /// Failed construction is retried; no image or coefficient data are copied.
    pub(crate) fn leb_provider(
        &self,
        skip_closed: bool,
    ) -> Result<&crate::positions::context_serving::BoundProvider<'static>> {
        self.leb_provider.provider(self, skip_closed)
    }

    /// Borrow the retained attach-order or tier-priority composite dispatch.
    /// Its image owners, coefficients and error policies match `leb_provider`.
    pub(crate) fn leb_composite(
        &self,
        policy: crate::positions::TierPolicy,
        skip_closed: bool,
    ) -> Result<&crate::leb::CompositeIndex<'static>> {
        self.leb_provider.composite(self, policy, skip_closed)
    }

    /// Context-local SPK kernel store and body map.
    pub fn spk_registry(&self) -> &SpkRegistry {
        &self.spk_registry
    }

    /// Mutable access to the context-local SPK kernel store (register
    /// images and bodies; see [`SpkRegistry`]).
    /// Registrations survive [`EphemerisContext::reset_session`]; a
    /// changed JPL path or [`EphemerisContext::close`] releases them.
    pub fn spk_registry_mut(&mut self) -> &mut SpkRegistry {
        self.leb_provider = Default::default();
        &mut self.spk_registry
    }

    /// Cached derived angles for the current calculation session.
    ///
    /// Release `EphemerisContext.get_angles_cache` semantics: a borrowed
    /// view of the per-context `{name: longitude}` map.
    pub fn angles_cache(&self) -> &HashMap<String, f64> {
        &self.angles_cache
    }

    /// Replace the cached derived angles.
    ///
    /// Release `EphemerisContext.set_angles_cache` semantics: the map is
    /// copied in, so later external mutation cannot leak into the context.
    pub fn set_angles_cache(&mut self, angles: &HashMap<String, f64>) {
        self.angles_cache = angles.clone();
    }

    /// Clear the cached derived angles.
    ///
    /// Release `EphemerisContext.clear_angles_cache` semantics.
    pub fn clear_angles_cache(&mut self) {
        self.angles_cache.clear();
    }

    /// User Delta-T override in days, if set.
    pub fn user_delta_t_days(&self) -> Option<f64> {
        self.user_delta_t_days
    }

    /// Explicit tidal-acceleration override, if set.
    pub fn tidal_acceleration_override(&self) -> Option<f64> {
        self.tidal_acceleration
    }

    /// Whether a tidal acceleration was explicitly set.
    pub fn tidal_acceleration_is_set(&self) -> bool {
        self.tidal_acceleration.is_some()
    }

    /// Effective tidal acceleration:
    /// the override when set, else [`TIDAL_DEFAULT`].
    pub fn effective_tidal_acceleration(&self) -> f64 {
        self.tidal_acceleration.unwrap_or(TIDAL_DEFAULT)
    }

    /// Programmatic Delta-T model override, if set (`None` defers to the
    /// environment/TOML/default chain).
    pub fn delta_t_model_override(&self) -> Option<DeltaTModel> {
        self.delta_t_model
    }

    /// Programmatic IERS Delta-T toggle, if set (`None` defers to the
    /// environment/TOML/default chain).
    pub fn iers_delta_t_override(&self) -> Option<bool> {
        self.iers_delta_t
    }

    /// Programmatic automatic-SPK-download toggle, if set (`None` defers
    /// to the environment/TOML/default chain).
    pub fn auto_spk_download_override(&self) -> Option<bool> {
        self.auto_spk_download
    }

    /// Programmatic automatic-IERS-download toggle, if set (`None`
    /// defers to the environment/TOML/default chain).
    pub fn iers_auto_download_override(&self) -> Option<bool> {
        self.iers_auto_download
    }

    /// Effective automatic-IERS-download setting: the programmatic value when set, else the `false`
    /// default (the environment/TOML legs resolve binding-side via
    /// [`resolve_iers_auto_download`]).
    pub fn effective_iers_auto_download(&self) -> bool {
        self.iers_auto_download.unwrap_or(false)
    }

    /// Effective automatic-SPK-download setting:
    /// the programmatic value when set, else the `true` default (the
    /// environment/TOML legs resolve binding-side via
    /// [`resolve_auto_spk_download`]).
    pub fn effective_auto_spk_download(&self) -> bool {
        self.auto_spk_download.unwrap_or(true)
    }

    /// SPK download date padding in days on each side (`0` = none).
    pub fn spk_date_padding_days(&self) -> u64 {
        self.spk_date_padding_days
    }

    /// Programmatic strict-precision toggle, if set (`None` defers to
    /// the environment/TOML/default chain).
    pub fn strict_precision_override(&self) -> Option<bool> {
        self.strict_precision
    }

    /// Effective strict-precision setting: the programmatic value
    /// when set, else the `true` default (the environment/TOML legs
    /// resolve binding-side via [`resolve_strict_precision`]).
    pub fn effective_strict_precision(&self) -> bool {
        self.strict_precision.unwrap_or(true)
    }

    /// Effective observed-IERS Delta-T setting: the programmatic value when set, else the `false` default
    /// (the environment/TOML legs resolve binding-side via
    /// [`resolve_iers_enabled`]).
    pub fn effective_iers_delta_t(&self) -> bool {
        self.iers_delta_t.unwrap_or(false)
    }

    /// Explicit lapse-rate override, if set.
    pub fn lapse_rate_override(&self) -> Option<f64> {
        self.lapse_rate
    }

    /// Effective lapse rate in K/m:
    /// the override when set, else [`LAPSE_RATE_DEFAULT`].
    pub fn effective_lapse_rate(&self) -> f64 {
        self.lapse_rate.unwrap_or(LAPSE_RATE_DEFAULT)
    }

    /// Install a lazy boundary resolver for IERS acquisition configuration.
    ///
    /// The resolver is consumed only by an enabled IERS lookup after the user
    /// Delta-T override and finite-epoch gate. It may resolve caller process
    /// settings at an outer binding boundary; SDK callers can instead supply
    /// an explicit cache directory or resolved data directory. This setter
    /// does not invoke the resolver or discard already acquired rows.
    pub fn set_iers_boundary_resolver(
        &mut self,
        resolver: Option<crate::iers::IersBoundaryResolver>,
    ) {
        self.iers_source = self.iers_source.with_boundary(resolver);
    }

    /// Attach immutable parsed IERS rows without filesystem or network work.
    ///
    /// The snapshot owns its finals, leap-second and observed Delta-T inputs.
    /// Context clones retain it; model Delta-T uses it only while IERS is enabled.
    /// Independent frame/observer timescale tables do not consume these rows.
    pub fn attach_iers_rows(&mut self, rows: crate::iers::IersDataSnapshot) {
        self.iers_source = self.iers_source.with_rows(rows);
    }

    /// Clear this context's acquired IERS rows while preserving configuration.
    ///
    /// A later enabled lookup reloads from the configured source. Already
    /// borrowed/cloned serving contexts retain their immutable prior snapshot;
    /// no process-global data store is involved. Reset preserves the cache;
    /// close calls this method to release this context's rows.
    pub fn clear_iers_cache(&mut self) {
        self.iers_source = self.iers_source.cleared();
    }

    /// Borrow this context's immutable-row acquisition owner, not process state.
    /// Its original-UT1 seconds callback retains the context's row lifetime;
    /// no filesystem or environment read is performed by this accessor.
    pub fn iers_source(&self) -> &crate::iers::IersContextSource {
        &self.iers_source
    }

    /// Install a retained row generation while preserving this context's policy.
    /// The passed owner supplies its boundary/cache; stored context acquisition
    /// settings remain authoritative. No acquisition runs during installation.
    /// Already cloned contexts retain their preceding cache allocation.
    pub fn attach_iers_source(&mut self, source: crate::iers::IersContextSource) {
        self.iers_source = source;
        self.refresh_iers_source();
    }

    /// Explicitly load/reload owned IERS rows and return the nonempty verdict.
    /// Force requests refresh only under enabled automatic-download policy.
    /// Successful publication replaces this context's row generation; a typed
    /// acquisition failure leaves it unchanged. Model/frame settings are intact.
    pub fn load_iers_data(&mut self, force_download: bool) -> crate::errors::Result<bool> {
        self.iers_source = self.iers_source.loaded(force_download)?;
        Ok(self.iers_source.loaded_counts().verdict())
    }

    /// Inspect the selected filesystem cache and owned row counts/ranges.
    /// Counts precede cold range acquisition. The explicit clock is Unix seconds;
    /// reported ranges are JD UT1. This can acquire rows without enabling their
    /// use by the model, and never changes independent frame/observer clocks.
    pub fn iers_cache_info(
        &self,
        now_unix_seconds: f64,
    ) -> crate::errors::Result<crate::iers::IersCacheInfo> {
        self.iers_source.cache_info(now_unix_seconds)
    }

    /// Explicitly delete the three configured disk files and clear this owner.
    /// Individual file errors are swallowed; resolution/directory failure is
    /// returned before replacement. Other contexts and in-flight row Arcs stay
    /// alive. Use only a caller-selected writable cache for this operation.
    pub fn delete_iers_cache_files(&mut self) -> crate::errors::Result<usize> {
        let (source, count) = self.iers_source.deleted_files()?;
        self.iers_source = source;
        Ok(count)
    }

    /// Synchronize stored acquisition inputs without resolving or reading them.
    fn refresh_iers_source(&mut self) {
        self.iers_source = self
            .iers_source
            .with_request(crate::iers::IersAcquisitionRequest {
                cache_dir: self.iers_cache_dir.clone(),
                data_dir: self.jpl_data_dir.clone(),
                auto_download: self.iers_auto_download,
                network_policy: self.network_policy,
                mode: self.mode,
            });
    }

    /// Build model Delta-T configuration for a UT1 query.
    ///
    /// User overrides are returned in days without any row acquisition.
    /// When enabled, IERS first supplies observed monthly seconds, then
    /// finals-derived TAI-UTC + 32.184 - UT1-UTC. Missing coverage or an ordinary
    /// acquisition/query error falls through to the selected model. Tidal
    /// adjustment is added to computed seconds before division by 86400.
    /// The callback borrows this context, retains immutable owned rows and
    /// performs acquisition only if consumed. Panics are not caught.
    /// Independent frame/observer TT and bundled timescale inputs are unchanged.
    pub fn delta_t_config(&self, model: DeltaTModel, iers_enabled: bool) -> DeltaTConfig<'_> {
        DeltaTConfig {
            model,
            user_override_days: self.user_delta_t_days,
            tidal_acceleration: self.effective_tidal_acceleration(),
            iers_enabled,
            iers_lookup_seconds: iers_enabled
                .then(|| self.iers_source.provider() as &dyn Fn(f64) -> Option<f64>),
            smh2016_lookup_seconds: Some(&crate::smh2016::smh2016_lookup_seconds),
        }
    }

    /// One native configured Delta-T value and derivative, days and days/day.
    /// A user constant wins before input validation or row acquisition. IERS
    /// value/rate share one actual immutable owner and one acquisition attempt;
    /// a missing source selects the model once. Espenak-Meeus is month-constant
    /// apart from the actual tidal term; SMH differentiates its three lookups.
    pub(crate) fn delta_t_value_and_rate_days(
        &self,
        jd: f64,
        model: DeltaTModel,
        iers: bool,
    ) -> Result<(f64, f64)> {
        let cfg = self.delta_t_config(model, iers);
        if let Some(value) = cfg.user_override_days {
            return Ok((value, 0.0));
        }
        if !jd.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd });
        }
        let adjust = crate::time::tid_acc_adjustment_seconds(jd, cfg.tidal_acceleration);
        let tidal = if cfg.tidal_acceleration == crate::constants::TIDAL_DEFAULT
            || jd >= crate::time::TID_ACC_REFERENCE_JD
        {
            0.0
        } else {
            2.0 * crate::time::tid_acc_adjustment_coefficient()
                * (cfg.tidal_acceleration - crate::constants::TIDAL_DEFAULT)
                * (jd - crate::time::TID_ACC_REFERENCE_JD)
                / (36525.0 * 36525.0)
        };
        if iers {
            if let Some((value, rate)) = self.iers_source.delta_t_seconds_and_rate(jd) {
                return Ok(((value + adjust) / 86400.0, (rate + tidal) / 86400.0));
            }
        }
        // The model callback is bundled pure mathematics, never a new IERS read.
        let model_cfg = DeltaTConfig {
            iers_enabled: false,
            ..cfg
        };
        let value = crate::time::deltat_with_config(jd, &model_cfg)?;
        let model_rate = match model {
            DeltaTModel::EspenakMeeus => 0.0,
            DeltaTModel::Smh2016 => crate::smh2016::smh2016_lookup_rate_seconds_per_day(jd).ok_or(
                Error::DeltaTModelUnavailable {
                    model: model.as_name(),
                },
            )?,
        };
        Ok((value, (model_rate + tidal) / 86400.0))
    }

    /// Reset per-calculation state without releasing owned resources or
    /// clearing configuration.
    ///
    /// Clears the observer (topocentric location),
    /// the sidereal selection, and the angles cache. Preserves the data
    /// tier, the provider mode, the data paths, every Delta-T /
    /// refraction override, and the network-policy override — and every
    /// bound store: attached [`LebFile`]
    /// images and the context-local [`SpkRegistry`]
    /// stay warm across the reset.
    ///
    /// Use this between independent calculations that may use different
    /// observer/sidereal settings. Use [`EphemerisContext::close`] only
    /// for full teardown.
    pub fn reset_session(&mut self) {
        self.observer = None;
        self.sidereal = SiderealConfig::default();
        self.angles_cache.clear();
    }

    /// Reset session state and release resources while retaining process settings.
    ///
    /// First apply [`Self::reset_session`] to the observer, sidereal selection
    /// and angles cache. Then restore the default tier, clear the configured
    /// JPL directory/filename and SPK cache-directory/download overrides, and
    /// reset SPK padding to zero. Tidal acceleration, atmospheric lapse rate
    /// and observed-IERS Delta-T selection return to automatic settings.
    ///
    /// Preserve the provider mode, LEB directory/file selections, resolved main
    /// data directory, Delta-T model/value, strict-precision selection, IERS
    /// cache-directory/automatic-download settings and network policy.
    /// Release acquired IERS rows, bound [`LebFile`] images and [`SpkRegistry`]
    /// kernels/registrations; the next calculation acquires resources again
    /// from the retained settings.
    ///
    /// The context receives an empty Horizons cache with the same settings.
    /// Retained owners and in-flight calculations keep the previous cache alive.
    pub fn close(&mut self) {
        self.reset_session();
        self.tier = Tier::default();
        self.paths.jpl_dir = None;
        self.ephemeris_file = None;
        self.spk_cache_dir = None;
        self.auto_spk_download = None;
        self.spk_date_padding_days = 0;
        self.tidal_acceleration = None;
        self.iers_delta_t = None;
        self.clear_iers_cache();
        self.refresh_iers_source();
        self.lapse_rate = None;
        self.release_leb_images();
        self.spk_registry.clear();
        if let Some(client) = self.horizons.take() {
            let fresh = client.as_ref().clone();
            fresh.clear_cache();
            self.horizons = Some(std::sync::Arc::new(fresh));
        }
    }
}

/// Resolve the active precision tier: the programmatic override wins, then the
/// `LIBEPHEMERIS_PRECISION` environment value (ASCII case/whitespace
/// normalized, silently skipped when unknown or blank), then the TOML
/// `precision` key (same normalization, same silent skip), then the
/// `medium` default. Unknown configured spellings never fail here — the
/// loud rejection happens only in `set_precision_tier`
/// ([`Tier::from_name`] / [`EphemerisContext::with_tier_name`]).
pub fn resolve_tier(programmatic: Option<Tier>, env: Option<&str>, toml: Option<&str>) -> Tier {
    if let Some(tier) = programmatic {
        return tier;
    }
    for candidate in [env, toml].into_iter().flatten() {
        let normalized = candidate.trim().to_ascii_lowercase();
        if let Some(tier) = Tier::from_name(normalized.as_str()) {
            return tier;
        }
    }
    Tier::Medium
}

/// One precision-tier identity row: the canonical
/// tier name, the tier kernel filename, the minor-body SPK download window
/// (`"YYYY-MM-DD"`), and the release description. All strings use the
/// release wording in the release entry order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct TierRow {
    /// Canonical tier name (`"base"` / `"medium"` / `"extended"`).
    pub name: &'static str,
    /// Tier kernel filename (`"de440s.bsp"` / `"de440.bsp"` /
    /// `"de441.bsp"`).
    pub ephemeris_file: &'static str,
    /// SPK download window start (`"YYYY-MM-DD"`).
    pub spk_start: &'static str,
    /// SPK download window end (`"YYYY-MM-DD"`).
    pub spk_end: &'static str,
    /// Release description (e.g. `"Modern usage (1850-2150), ~31 MB"`).
    pub description: &'static str,
}

/// Static precision-tier identity table, in canonical tier order.
pub const TIER_TABLE: [TierRow; 3] = [
    TierRow {
        name: "base",
        ephemeris_file: "de440s.bsp",
        spk_start: "1850-01-01",
        spk_end: "2150-01-01",
        description: "Modern usage (1850-2150), ~31 MB",
    },
    TierRow {
        name: "medium",
        ephemeris_file: "de440.bsp",
        spk_start: "1900-01-01",
        spk_end: "2100-01-01",
        description: "General purpose (1550-2650), ~114 MB",
    },
    TierRow {
        name: "extended",
        ephemeris_file: "de441.bsp",
        spk_start: "1600-01-01",
        spk_end: "2500-01-01",
        description: "Extended range (-13200 to +17191), ~3.1 GB",
    },
];

/// Priority order of the cumulative precision tiers: a higher tier covers at least every lower tier's
/// minor-body SPK window.
pub const TIER_PRIORITY: [Tier; 3] = [Tier::Base, Tier::Medium, Tier::Extended];

/// Identity row for a tier.
pub fn tier_row(tier: Tier) -> &'static TierRow {
    &TIER_TABLE[tier as usize]
}

/// One tier identity string by field index: `0` name, `1` ephemeris file,
/// `2` SPK window start, `3` SPK window end, `4` description; every other
/// value reads as `""` (never an error — the two-call caller fetches only
/// after a positive length).
pub fn tier_field(tier: Tier, field: i32) -> &'static str {
    let row = tier_row(tier);
    match field {
        0 => row.name,
        1 => row.ephemeris_file,
        2 => row.spk_start,
        3 => row.spk_end,
        4 => row.description,
        _ => "",
    }
}

/// SPK date range for a tier: the union of the requested tier's own window
/// with every lower-priority tier's window, so the download window is
/// monotonic in tier order. ISO `"YYYY-MM-DD"` strings order
/// lexicographically, so `min`/`max` are correct. A tier outside the
/// priority order (impossible for [`Tier`]) would read its own window.
pub fn spk_date_range_for_tier(tier: Tier) -> (&'static str, &'static str) {
    let position = TIER_PRIORITY
        .iter()
        .position(|entry| *entry == tier)
        .map_or(1, |index| index + 1);
    let eligible = &TIER_PRIORITY[..position];
    let start = eligible
        .iter()
        .map(|entry| tier_row(*entry).spk_start)
        .min()
        .unwrap_or(tier_row(tier).spk_start);
    let end = eligible
        .iter()
        .map(|entry| tier_row(*entry).spk_end)
        .max()
        .unwrap_or(tier_row(tier).spk_end);
    (start, end)
}

/// Resolve the active calculation mode: the programmatic override wins, then the
/// `LIBEPHEMERIS_MODE` environment value (blank means unset), then the TOML
/// `mode` key (blank means unset), then `"auto"`. A non-blank but unknown
/// configured spelling fails loudly with [`Error::InvalidCalcMode`]
/// naming its source — a typo must never widen a sealed deployment back to
/// `auto` — mirroring the release `_validated_calc_mode` guard.
pub fn resolve_calc_mode(
    programmatic: Option<ProviderMode>,
    env: Option<&str>,
    toml: Option<&str>,
) -> Result<ProviderMode> {
    if let Some(mode) = programmatic {
        return Ok(mode);
    }
    for (candidate, source) in [(env, CALC_MODE_ENV_VAR), (toml, "TOML key 'mode'")] {
        if let Some(raw) = candidate {
            if raw.trim().is_empty() {
                continue;
            }
            return ProviderMode::from_calc_mode_name(raw).ok_or_else(|| Error::InvalidCalcMode {
                value: raw.to_string(),
                source,
            });
        }
    }
    Ok(ProviderMode::Auto)
}

/// Resolve the active Delta-T model: the programmatic override wins, then the
/// `LIBEPHEMERIS_DELTAT_MODEL` environment value, then the TOML
/// `deltat_model` key, then the `smh2016` default. Unknown configured
/// spellings are silently skipped — the loud rejection lives in
/// [`EphemerisContext::with_delta_t_model_name`].
pub fn resolve_delta_t_model(
    programmatic: Option<DeltaTModel>,
    env: Option<&str>,
    toml: Option<&str>,
) -> DeltaTModel {
    if let Some(model) = programmatic {
        return model;
    }
    for candidate in [env, toml].into_iter().flatten() {
        if let Some(model) = DeltaTModel::from_name(candidate) {
            return model;
        }
    }
    DeltaTModel::Smh2016
}

/// Resolve whether observed-IERS Delta T outranks the selected model.
/// The programmatic value wins, then the `LIBEPHEMERIS_IERS_DELTA_T` environment value
/// (`1/true/yes/on/enabled` enable; `0/false/no/off/disabled` disable;
/// anything else, including blank, falls through), then the TOML
/// `iers_delta_t` key, then the `false` default.
pub fn resolve_iers_enabled(
    programmatic: Option<bool>,
    env: Option<&str>,
    toml: Option<bool>,
) -> bool {
    resolve_boolean_setting(programmatic, env, toml, false)
}

/// Resolve the automatic-SPK-download setting: the programmatic value wins, then
/// the `LIBEPHEMERIS_AUTO_SPK` environment value (`0/false/no/off/disabled`
/// disable; `1/true/yes/on/enabled` enable; anything else, including blank,
/// falls through), then the TOML `auto_spk` key, then the `true` default.
///
/// Inputs are already captured. Resolution performs no environment access,
/// file I/O or resource acquisition. Environment values use Unicode whitespace
/// trimming and ASCII case normalization; unknown values do not produce errors.
pub fn resolve_auto_spk_download(
    programmatic: Option<bool>,
    env: Option<&str>,
    toml: Option<bool>,
) -> bool {
    resolve_boolean_setting(programmatic, env, toml, true)
}

/// Resolve the strict-precision setting: the programmatic value wins, then the
/// `LIBEPHEMERIS_STRICT_PRECISION` environment value (same word lists as
/// [`resolve_auto_spk_download`]; anything else, including blank, falls
/// through), then the TOML `strict_precision` key, then the `true`
/// default. Same already-read-inputs contract as
/// [`resolve_auto_spk_download`].
pub fn resolve_strict_precision(
    programmatic: Option<bool>,
    env: Option<&str>,
    toml: Option<bool>,
) -> bool {
    resolve_boolean_setting(programmatic, env, toml, true)
}

/// Resolve the automatic-IERS-download setting: the programmatic value wins,
/// then the `LIBEPHEMERIS_IERS_AUTO_DOWNLOAD` environment value (same
/// word lists as [`resolve_auto_spk_download`]; anything else,
/// including blank, falls through), then the TOML
/// `iers_auto_download` key, then the `false` default. Same
/// already-read-inputs contract as [`resolve_auto_spk_download`].
pub fn resolve_iers_auto_download(
    programmatic: Option<bool>,
    env: Option<&str>,
    toml: Option<bool>,
) -> bool {
    resolve_boolean_setting(programmatic, env, toml, false)
}

/// Apply the common boolean precedence while keeping each setting's default explicit.
/// An explicit false is just as authoritative as true; neither consumes `env`.
fn resolve_boolean_setting(
    programmatic: Option<bool>,
    env: Option<&str>,
    toml: Option<bool>,
    default: bool,
) -> bool {
    if let Some(enabled) = programmatic {
        return enabled;
    }
    if let Some(raw) = env {
        match raw.trim().to_ascii_lowercase().as_str() {
            "1" | "true" | "yes" | "on" | "enabled" => return true,
            "0" | "false" | "no" | "off" | "disabled" => return false,
            _ => {}
        }
    }
    toml.unwrap_or(default)
}

/// Resolve the configured network policy: the programmatic override
/// wins, then the `LIBEPHEMERIS_NETWORK_POLICY` environment value (blank
/// means unset), then the TOML `network_policy` key (blank means unset),
/// then `"auto"`. A non-blank but unknown spelling fails loudly with
/// [`Error::InvalidNetworkPolicy`] naming its source — a typo in a
/// security-relevant boundary must never silently widen it to
/// `"auto"`/`"allow"` — mirroring the release `_validated_policy` guard.
///
/// The inputs are already-read values (the binding layer owns reading the
/// process environment, `.env` files, and TOML documents), so resolution
/// stays deterministic and testable: tests pass explicit values and never
/// touch the process environment.
pub fn resolve_network_policy(
    programmatic: Option<NetworkPolicy>,
    env: Option<&str>,
    toml: Option<&str>,
) -> Result<NetworkPolicy> {
    if let Some(policy) = programmatic {
        return Ok(policy);
    }
    for (candidate, source) in [
        (env, NETWORK_POLICY_ENV_VAR),
        (toml, "TOML key 'network_policy'"),
    ] {
        if let Some(raw) = candidate {
            if raw.trim().is_empty() {
                continue;
            }
            return NetworkPolicy::from_policy_name(raw).ok_or_else(|| {
                Error::InvalidNetworkPolicy {
                    value: raw.to_string(),
                    source,
                }
            });
        }
    }
    Ok(NetworkPolicy::Auto)
}

/// Resolve the effective network gate:
/// an explicit `"allow"`/`"sealed"` passes through; `"auto"` seals the
/// sealed-LEB mode and allows every other mode. The `download`
/// provisioning command selects [`NetworkPolicy::Allow`] explicitly
/// (download-provisioning layer), since a download is provisioning, not calculation egress.
pub fn effective_network_policy(
    configured: NetworkPolicy,
    mode: ProviderMode,
) -> EffectiveNetworkPolicy {
    match configured {
        NetworkPolicy::Allow => EffectiveNetworkPolicy::Allow,
        NetworkPolicy::Sealed => EffectiveNetworkPolicy::Sealed,
        NetworkPolicy::Auto => match mode {
            ProviderMode::SealedLeb => EffectiveNetworkPolicy::Sealed,
            _ => EffectiveNetworkPolicy::Allow,
        },
    }
}

/// Refuse network I/O before any socket work when the effective policy is
/// sealed: the NET provider call sites pass
/// their operation label as `purpose`. Mirrors the release
/// `NetworkSealedError`, whose message names the blocked operation and
/// points at explicit provisioning.
///
/// Kernel selection may provision a genuinely missing configured image under
/// an allowed policy. Such acquisition, explicit asset provisioning and the
/// remote provider must pass this gate before socket work; attached LEB/SPK
/// snapshots and existing local images do not require network permission.
/// The check authorizes the labelled operation only. It does not download,
/// select a kernel, mutate a registry or replace a configured provider.
pub fn require_network(effective: EffectiveNetworkPolicy, purpose: &str) -> Result<()> {
    match effective {
        EffectiveNetworkPolicy::Allow => Ok(()),
        EffectiveNetworkPolicy::Sealed => Err(Error::NetworkSealed {
            purpose: purpose.to_string(),
        }),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn independent_sidereal_bits_preserve_negative_base_and_raw_epoch_words() {
        use crate::constants::SIDBIT_ECL_T0;
        use crate::sidereal::IMPLEMENTED_SIDBITS;
        let epoch = f64::from_bits(0x7ff8_0000_0000_0042);
        let mut context = EphemerisContext::new();
        assert!(context.set_sidereal_mode(-1, epoch, -0.0).is_empty());
        context.set_sidereal_bits(SIDBIT_ECL_T0);
        let (base, t0, offset) = context.sidereal_mode_full();
        assert_eq!(base, -1);
        assert_eq!(t0.to_bits(), epoch.to_bits());
        assert_eq!(offset.to_bits(), (-0.0_f64).to_bits());
        assert_eq!(context.sidereal_bits(), SIDBIT_ECL_T0);
        context.set_sidereal_bits(i32::MAX);
        assert_eq!(context.sidereal_bits(), IMPLEMENTED_SIDBITS);
        context.set_sidereal_bits(0);
        assert_eq!(context.sidereal_bits(), 0);
        assert_eq!(context.sidereal_mode_full().0, -1);
    }

    #[test]
    fn independent_sidereal_bits_builder_and_session_snapshots_are_owned() {
        use crate::constants::{SIDBIT_ECL_T0, SIDBIT_SSY_PLANE};
        let mut context = EphemerisContext::new();
        context.set_sidereal_mode(-5, f64::INFINITY, -0.0);
        let mut context = context.with_sidereal_bits(SIDBIT_ECL_T0);
        let session = crate::engine_api::Session::from_context(context.clone());
        context.set_sidereal_bits(SIDBIT_SSY_PLANE);
        assert_eq!(session.context().sidereal_mode(), -5);
        assert_eq!(session.context().sidereal_bits(), SIDBIT_ECL_T0);
        assert_eq!(context.sidereal_bits(), SIDBIT_SSY_PLANE);
        assert_eq!(session.context().sidereal_mode_full().1, f64::INFINITY);
        assert_eq!(
            session.context().sidereal_mode_full().2.to_bits(),
            (-0.0_f64).to_bits()
        );
        context.reset_session();
        assert_eq!(context.sidereal_bits(), 0);
        assert_eq!(session.context().sidereal_bits(), SIDBIT_ECL_T0);
    }

    fn assert_send_sync<T: Send + Sync>() {}

    #[test]
    fn defaults_are_documented() {
        let ctx = EphemerisContext::new();
        assert_eq!(ctx.observer(), None);
        assert_eq!(ctx.sidereal(), SiderealConfig::Tropical);
        // Release `_DEFAULT_TIER` is "medium" ( reconciliation; the
        // The default precision tier is medium. The provider mode is
        // the release `"auto"` global default ( routing).
        assert_eq!(ctx.tier(), Tier::Medium);
        assert_eq!(ctx.mode(), ProviderMode::Auto);
        assert_eq!(ctx.paths().leb_dir, None);
        assert_eq!(ctx.paths().jpl_dir, None);
        // Delta-T / refraction overrides all unset (automatic).
        assert_eq!(ctx.user_delta_t_days(), None);
        assert_eq!(ctx.tidal_acceleration_override(), None);
        assert!(!ctx.tidal_acceleration_is_set());
        assert_eq!(ctx.effective_tidal_acceleration(), TIDAL_DEFAULT);
        assert_eq!(ctx.delta_t_model_override(), None);
        assert_eq!(ctx.iers_delta_t_override(), None);
        assert_eq!(ctx.lapse_rate_override(), None);
        assert_eq!(ctx.effective_lapse_rate(), LAPSE_RATE_DEFAULT);
        assert_send_sync::<EphemerisContext>();
    }

    #[test]
    fn sidereal_mode_setting_is_context_local() {
        use crate::errors::Warning;
        let mut ctx = EphemerisContext::new();
        // Never-set readback matches the release default (mode 0, zero
        // epoch pair, no flags).
        assert_eq!(ctx.sidereal(), SiderealConfig::Tropical);
        assert_eq!(ctx.sidereal_mode(), 0);
        assert_eq!(ctx.sidereal_mode_full(), (0, 0.0, 0.0));
        assert_eq!(ctx.sidereal_bits(), 0);
        // Plain mode + SIDM_USER epoch pair round-trip silently.
        assert!(ctx.set_sidereal_mode(1, 0.0, 0.0).is_empty());
        assert_eq!(ctx.sidereal_mode_full(), (1, 0.0, 0.0));
        assert!(ctx.set_sidereal_mode(255, 2451545.0, 24.0).is_empty());
        assert_eq!(ctx.sidereal_mode_full(), (255, 2451545.0, 24.0));
        // SIDBIT composite: implemented flags retained, base stored.
        assert!(ctx.set_sidereal_mode(256 | 1, 0.0, 0.0).is_empty());
        assert_eq!(ctx.sidereal_mode(), 1);
        assert_eq!(ctx.sidereal_bits(), 256);
        // Unknown high bits: stripped with a diagnostic, selection stored.
        assert_eq!(
            ctx.set_sidereal_mode(0x4000 | 2, 0.0, 0.0),
            vec![Warning::UnsupportedSiderealProjection { bits: 0x4000 }]
        );
        assert_eq!(ctx.sidereal_mode(), 2);
        // Selecting Tropical clears back to the never-set readback.
        ctx.sidereal = SiderealConfig::Tropical;
        assert_eq!(ctx.sidereal_mode_full(), (0, 0.0, 0.0));
        assert_eq!(ctx.sidereal_bits(), 0);
        // Independent contexts never share the selection.
        let mut other = EphemerisContext::new();
        other.set_sidereal_mode(3, 0.0, 0.0);
        assert_eq!(ctx.sidereal_mode(), 0);
        assert_eq!(other.sidereal_mode(), 3);
    }

    #[test]
    fn reset_and_close_clear_the_sidereal_selection() {
        let mut ctx = EphemerisContext::new();
        ctx.set_sidereal_mode(256 | 5, 2415020.0, 21.0);
        assert_eq!(ctx.sidereal_mode(), 5);
        ctx.reset_session();
        assert_eq!(ctx.sidereal(), SiderealConfig::Tropical);
        assert_eq!(ctx.sidereal_mode_full(), (0, 0.0, 0.0));
        assert_eq!(ctx.sidereal_bits(), 0);
        ctx.set_sidereal_mode(1, 0.0, 0.0);
        ctx.close();
        assert_eq!(ctx.sidereal(), SiderealConfig::Tropical);
        assert_eq!(ctx.sidereal_mode_full(), (0, 0.0, 0.0));
    }

    #[test]
    fn sidereal_triples_and_bits_are_isolated_between_contexts() {
        // SID-05: two live contexts hold different full selections —
        // mode, epoch pair and retained SIDBITs — without interference,
        // and clearing one never touches the other.
        let mut lahiri = EphemerisContext::new();
        let mut user = EphemerisContext::new();
        assert!(lahiri.set_sidereal_mode(256 | 1, 0.0, 0.0).is_empty());
        assert!(user.set_sidereal_mode(255, 2451545.0, 24.0).is_empty());
        assert_eq!(lahiri.sidereal_mode_full(), (1, 0.0, 0.0));
        assert_eq!(lahiri.sidereal_bits(), 256);
        assert_eq!(user.sidereal_mode_full(), (255, 2451545.0, 24.0));
        assert_eq!(user.sidereal_bits(), 0);
        // Re-selecting on one context leaves the other intact.
        assert!(lahiri.set_sidereal_mode(2, 0.0, 0.0).is_empty());
        assert_eq!(lahiri.sidereal_mode_full(), (2, 0.0, 0.0));
        assert_eq!(lahiri.sidereal_bits(), 0);
        assert_eq!(user.sidereal_mode_full(), (255, 2451545.0, 24.0));
        // Resetting one context clears only that context.
        lahiri.reset_session();
        assert_eq!(lahiri.sidereal_mode_full(), (0, 0.0, 0.0));
        assert_eq!(lahiri.sidereal_bits(), 0);
        assert_eq!(user.sidereal_mode_full(), (255, 2451545.0, 24.0));
        user.close();
        assert_eq!(user.sidereal(), SiderealConfig::Tropical);
        assert_eq!(lahiri.sidereal(), SiderealConfig::Tropical);
    }

    #[test]
    fn builders_set_each_option() {
        let obs = Observer::new(13.4, 52.5, 34.0).unwrap();
        let ctx = EphemerisContext::new()
            .with_observer(obs)
            .with_tier(Tier::Extended)
            .with_mode(ProviderMode::Auto)
            .with_paths(DataPaths {
                leb_dir: Some(PathBuf::from("/data/leb")),
                jpl_dir: None,
            });
        assert_eq!(ctx.observer(), Some(obs));
        assert_eq!(ctx.tier(), Tier::Extended);
        assert_eq!(ctx.mode(), ProviderMode::Auto);
        assert_eq!(ctx.paths().leb_dir, Some(PathBuf::from("/data/leb")));
        // Clearing restores geocentric.
        assert_eq!(ctx.without_observer().observer(), None);
    }

    #[test]
    fn contexts_are_independent_values() {
        let a = EphemerisContext::new().with_tier(Tier::Base);
        let b = EphemerisContext::new();
        // Mutating (rebuilding) one context never affects the other: there
        // is no shared mutable state to leak through.
        let a2 = a
            .clone()
            .with_observer(Observer::new(0.0, 51.5, 11.0).unwrap())
            .with_mode(ProviderMode::HorizonsRemote);
        assert_eq!(a.tier(), Tier::Base);
        assert_eq!(a.observer(), None);
        assert_eq!(b.tier(), Tier::Medium);
        assert_eq!(b.mode(), ProviderMode::Auto);
        assert_eq!(a2.mode(), ProviderMode::HorizonsRemote);
    }

    #[test]
    fn observer_validation_retains_latitude_first_error() {
        let cases: [(f64, f64, f64, &str, f64); 4] = [
            (361.0, 91.0, 0.0, "lat_deg", 91.0),
            (-181.0, -91.0, 0.0, "lat_deg", -91.0),
            (361.0, 0.0, f64::NAN, "lon_deg", 361.0),
            (0.0, 91.0, f64::NAN, "lat_deg", 91.0),
        ];
        for (lon, lat, height, field, value) in cases {
            assert!(matches!(Observer::new(lon, lat, height),
                Err(Error::InvalidObserver { field: actual, value: actual_value })
                    if actual == field && actual_value.to_bits() == value.to_bits()));
        }
    }

    #[test]
    fn observer_validation() {
        assert!(Observer::new(0.0, 0.0, 0.0).is_ok());
        assert!(Observer::new(-180.0, -90.0, -430.0).is_ok());
        assert!(Observer::new(180.0, 90.0, 8848.0).is_ok());
        // Release `set_topo` contract : the east-positive
        // convention passes through as-is (periodic geometry downstream).
        assert!(Observer::new(200.0, 0.0, 0.0).is_ok());
        assert!(Observer::new(360.0, 0.0, 0.0).is_ok());
        for (lon, lat, h) in [
            (361.0, 0.0, 0.0),
            (-181.0, 0.0, 0.0),
            (400.0, 0.0, 0.0),
            (0.0, 91.0, 0.0),
            (0.0, -91.0, 0.0),
            (0.0, 0.0, f64::NAN),
            (f64::INFINITY, 0.0, 0.0),
        ] {
            assert!(
                matches!(
                    Observer::new(lon, lat, h),
                    Err(Error::InvalidObserver { .. })
                ),
                "accepted invalid observer ({lon}, {lat}, {h})"
            );
        }
    }

    #[test]
    fn tier_names_are_exact_and_round_trip() {
        // Release `set_precision_tier` matches exactly: no normalization.
        for (tier, name) in [
            (Tier::Base, "base"),
            (Tier::Medium, "medium"),
            (Tier::Extended, "extended"),
        ] {
            assert_eq!(tier.as_name(), name);
            assert_eq!(Tier::from_name(name), Some(tier));
            assert_eq!(
                EphemerisContext::new().with_tier_name(name).unwrap().tier(),
                tier
            );
        }
        // Near-misses the release rejects with `ValueError`.
        for bad in ["", "Medium", " medium", "medium ", "ultra", "bas"] {
            assert_eq!(Tier::from_name(bad), None);
            assert!(
                matches!(
                    EphemerisContext::new().with_tier_name(bad),
                    Err(Error::InvalidTier { .. })
                ),
                "accepted tier {bad:?}"
            );
        }
        // Release default pin.
        assert_eq!(Tier::default(), Tier::Medium);
    }

    #[test]
    fn calc_mode_names_normalize_and_round_trip() {
        for (mode, name) in [
            (ProviderMode::Auto, "auto"),
            (ProviderMode::NativeJpl, "skyfield"),
            (ProviderMode::SealedLeb, "leb"),
            (ProviderMode::HorizonsRemote, "horizons"),
        ] {
            assert_eq!(mode.as_calc_mode_name(), name);
            assert_eq!(ProviderMode::from_calc_mode_name(name), Some(mode));
        }
        // Release `set_calc_mode` lower/strips before matching.
        assert_eq!(
            ProviderMode::from_calc_mode_name(" LEB "),
            Some(ProviderMode::SealedLeb)
        );
        assert_eq!(
            ProviderMode::from_calc_mode_name("SkyField"),
            Some(ProviderMode::NativeJpl)
        );
        assert_eq!(
            EphemerisContext::new()
                .with_mode_name("horizons")
                .unwrap()
                .mode(),
            ProviderMode::HorizonsRemote
        );
        // Typos fail loudly naming the setter source.
        match EphemerisContext::new().with_mode_name("lebb") {
            Err(Error::InvalidCalcMode { value, source }) => {
                assert_eq!(value, "lebb");
                assert_eq!(source, CALC_MODE_SOURCE_SETTER);
            }
            other => panic!("expected InvalidCalcMode, got {other:?}"),
        }
    }

    #[test]
    fn resolve_calc_mode_precedence_and_loud_typos() {
        // Programmatic override wins over everything, even a bad env value.
        assert_eq!(
            resolve_calc_mode(Some(ProviderMode::SealedLeb), Some("bogus"), Some("bogus")),
            Ok(ProviderMode::SealedLeb)
        );
        // Env wins over TOML; TOML wins over the `auto` default.
        assert_eq!(
            resolve_calc_mode(None, Some("leb"), Some("horizons")),
            Ok(ProviderMode::SealedLeb)
        );
        assert_eq!(
            resolve_calc_mode(None, None, Some("horizons")),
            Ok(ProviderMode::HorizonsRemote)
        );
        assert_eq!(resolve_calc_mode(None, None, None), Ok(ProviderMode::Auto));
        // Blank configured values mean "unset", never a typo.
        assert_eq!(
            resolve_calc_mode(None, Some("  "), Some("leb")),
            Ok(ProviderMode::SealedLeb)
        );
        assert_eq!(
            resolve_calc_mode(None, Some(""), None),
            Ok(ProviderMode::Auto)
        );
        // Non-blank typos fail loudly naming their source.
        match resolve_calc_mode(None, Some("lebb"), None) {
            Err(Error::InvalidCalcMode { value, source }) => {
                assert_eq!(value, "lebb");
                assert_eq!(source, CALC_MODE_ENV_VAR);
            }
            other => panic!("expected InvalidCalcMode, got {other:?}"),
        }
        match resolve_calc_mode(None, None, Some("auto ")) {
            // Trailing space normalizes away: valid.
            Ok(ProviderMode::Auto) => {}
            other => panic!("expected normalized auto, got {other:?}"),
        }
        match resolve_calc_mode(None, None, Some("skifield")) {
            Err(Error::InvalidCalcMode { source, .. }) => {
                assert_eq!(source, "TOML key 'mode'");
            }
            other => panic!("expected InvalidCalcMode, got {other:?}"),
        }
    }

    #[test]
    fn network_policy_names_normalize_and_round_trip() {
        assert_eq!(NetworkPolicy::Auto.as_policy_name(), "auto");
        assert_eq!(NetworkPolicy::Allow.as_policy_name(), "allow");
        assert_eq!(NetworkPolicy::Sealed.as_policy_name(), "sealed");
        assert_eq!(
            NetworkPolicy::from_policy_name(" SEALED "),
            Some(NetworkPolicy::Sealed)
        );
        assert_eq!(
            NetworkPolicy::from_policy_name("Allow"),
            Some(NetworkPolicy::Allow)
        );
        assert_eq!(
            NetworkPolicy::from_policy_name("auto"),
            Some(NetworkPolicy::Auto)
        );
        // Unknown spellings are `None` here; the loud rejection lives in
        // the resolver and the validated name setter.
        assert_eq!(NetworkPolicy::from_policy_name("opne"), None);
        assert_eq!(NetworkPolicy::from_policy_name(""), None);
    }

    #[test]
    fn resolve_network_policy_precedence_and_loud_typos() {
        // Programmatic override wins over everything, even a bad env value.
        assert_eq!(
            resolve_network_policy(Some(NetworkPolicy::Allow), Some("bogus"), Some("bogus")),
            Ok(NetworkPolicy::Allow)
        );
        // Env wins over TOML; TOML wins over the `auto` default.
        assert_eq!(
            resolve_network_policy(None, Some("sealed"), Some("allow")),
            Ok(NetworkPolicy::Sealed)
        );
        assert_eq!(
            resolve_network_policy(None, None, Some("allow")),
            Ok(NetworkPolicy::Allow)
        );
        assert_eq!(
            resolve_network_policy(None, None, None),
            Ok(NetworkPolicy::Auto)
        );
        // Blank configured values mean "unset", never a typo.
        assert_eq!(
            resolve_network_policy(None, Some("  "), Some("sealed")),
            Ok(NetworkPolicy::Sealed)
        );
        assert_eq!(
            resolve_network_policy(None, Some(""), None),
            Ok(NetworkPolicy::Auto)
        );
        // Non-blank typos fail loudly naming their source: a typo must
        // never silently widen a sealed deployment.
        match resolve_network_policy(None, Some("opne"), None) {
            Err(Error::InvalidNetworkPolicy { value, source }) => {
                assert_eq!(value, "opne");
                assert_eq!(source, NETWORK_POLICY_ENV_VAR);
            }
            other => panic!("expected InvalidNetworkPolicy, got {other:?}"),
        }
        match resolve_network_policy(None, None, Some("AUTO ")) {
            // Trailing space normalizes away: valid.
            Ok(NetworkPolicy::Auto) => {}
            other => panic!("expected normalized auto, got {other:?}"),
        }
        match resolve_network_policy(None, None, Some("selad")) {
            Err(Error::InvalidNetworkPolicy { source, .. }) => {
                assert_eq!(source, "TOML key 'network_policy'");
            }
            other => panic!("expected InvalidNetworkPolicy, got {other:?}"),
        }
    }

    #[test]
    fn effective_network_policy_seals_only_leb_under_auto() {
        // `auto` seals the sealed-LEB mode and allows every other mode.
        assert_eq!(
            effective_network_policy(NetworkPolicy::Auto, ProviderMode::SealedLeb),
            EffectiveNetworkPolicy::Sealed
        );
        for mode in [
            ProviderMode::Auto,
            ProviderMode::NativeJpl,
            ProviderMode::HorizonsRemote,
        ] {
            assert_eq!(
                effective_network_policy(NetworkPolicy::Auto, mode),
                EffectiveNetworkPolicy::Allow,
                "mode {mode:?} must be allowed under auto"
            );
        }
        // Explicit overrides pass through on every mode, including the
        // `download` provisioning selection (`allow`) under sealed LEB.
        for mode in [
            ProviderMode::SealedLeb,
            ProviderMode::Auto,
            ProviderMode::NativeJpl,
            ProviderMode::HorizonsRemote,
        ] {
            assert_eq!(
                effective_network_policy(NetworkPolicy::Allow, mode),
                EffectiveNetworkPolicy::Allow
            );
            assert_eq!(
                effective_network_policy(NetworkPolicy::Sealed, mode),
                EffectiveNetworkPolicy::Sealed
            );
        }
    }

    #[test]
    fn require_network_gates_sealed_before_socket_work() {
        assert!(require_network(EffectiveNetworkPolicy::Allow, "Horizons query").is_ok());
        match require_network(EffectiveNetworkPolicy::Sealed, "Horizons query") {
            Err(Error::NetworkSealed { purpose }) => {
                assert_eq!(purpose, "Horizons query");
            }
            other => panic!("expected NetworkSealed, got {other:?}"),
        }
        // The message names the blocked operation and points at explicit
        // provisioning, mirroring the release `NetworkSealedError`.
        let message = Error::NetworkSealed {
            purpose: "IERS refresh".to_string(),
        }
        .to_string();
        assert!(message.contains("sealed"), "got {message}");
        assert!(message.contains("IERS refresh"), "got {message}");
    }

    #[test]
    fn network_policy_override_survives_reset_and_close() {
        // Default: no override; the default
        // context allows because its mode is the release `auto` default
        // `auto` seals sealed-LEB only.
        let ctx = EphemerisContext::new();
        assert_eq!(ctx.network_policy_override(), None);
        assert_eq!(
            ctx.effective_network_policy(),
            EffectiveNetworkPolicy::Allow
        );
        assert!(ctx.require_network("probe").is_ok());
        // An explicit sealed-LEB mode still seals under the default chain.
        let sealed = EphemerisContext::new().with_mode(ProviderMode::SealedLeb);
        assert_eq!(
            sealed.effective_network_policy(),
            EffectiveNetworkPolicy::Sealed
        );
        assert!(sealed.require_network("probe").is_err());
        // An explicit override is authoritative for the gate.
        let ctx = EphemerisContext::new().with_network_policy(Some(NetworkPolicy::Allow));
        assert_eq!(ctx.network_policy_override(), Some(NetworkPolicy::Allow));
        assert_eq!(
            ctx.effective_network_policy(),
            EffectiveNetworkPolicy::Allow
        );
        assert!(ctx.require_network("probe").is_ok());
        // Validated name setter: normalization + loud typo + clearing.
        let ctx = EphemerisContext::new()
            .with_network_policy_name(Some(" Sealed "))
            .unwrap();
        assert_eq!(ctx.network_policy_override(), Some(NetworkPolicy::Sealed));
        match EphemerisContext::new().with_network_policy_name(Some("opne")) {
            Err(Error::InvalidNetworkPolicy { value, source }) => {
                assert_eq!(value, "opne");
                assert_eq!(source, NETWORK_POLICY_SOURCE_SETTER);
            }
            other => panic!("expected InvalidNetworkPolicy, got {other:?}"),
        }
        let cleared = ctx.with_network_policy_name(None).unwrap();
        assert_eq!(cleared.network_policy_override(), None);
        // The override is process-level configuration, not session state:
        // both `reset_session` and `close` preserve it.
        let mut ctx = EphemerisContext::new().with_network_policy(Some(NetworkPolicy::Sealed));
        ctx.reset_session();
        assert_eq!(ctx.network_policy_override(), Some(NetworkPolicy::Sealed));
        ctx.close();
        assert_eq!(ctx.network_policy_override(), Some(NetworkPolicy::Sealed));
        // Context-local: a sibling never sees the override.
        let sibling = EphemerisContext::new();
        assert_eq!(sibling.network_policy_override(), None);
    }

    #[test]
    fn resolve_tier_silently_skips_unknown_spellings() {
        // Programmatic override wins unconditionally.
        assert_eq!(
            resolve_tier(Some(Tier::Base), Some("extended"), Some("extended")),
            Tier::Base
        );
        // Env normalizes case/whitespace and wins over TOML.
        assert_eq!(
            resolve_tier(None, Some(" Extended "), Some("base")),
            Tier::Extended
        );
        assert_eq!(resolve_tier(None, None, Some("base")), Tier::Base);
        assert_eq!(resolve_tier(None, None, None), Tier::Medium);
        // Unknown configured spellings fall through silently.
        assert_eq!(resolve_tier(None, Some("ultra"), Some("base")), Tier::Base);
        assert_eq!(
            resolve_tier(None, Some("ultra"), Some("ultra")),
            Tier::Medium
        );
        assert_eq!(resolve_tier(None, Some(""), None), Tier::Medium);
    }

    #[test]
    fn user_delta_t_sentinel_clears_and_verbatim_stores() {
        // Default: no override (computed values).
        assert_eq!(EphemerisContext::new().user_delta_t_days(), None);
        // None clears; the DELTAT_AUTOMATIC sentinel clears too (it must
        // never be stored literally).
        let cleared = EphemerisContext::new()
            .with_user_delta_t(Some(0.001))
            .with_user_delta_t(None);
        assert_eq!(cleared.user_delta_t_days(), None);
        let sentinel = EphemerisContext::new()
            .with_user_delta_t(Some(0.001))
            .with_user_delta_t(Some(DELTAT_AUTOMATIC));
        assert_eq!(sentinel.user_delta_t_days(), None);
        // Any other value stores.
        let stored = EphemerisContext::new().with_user_delta_t(Some(0.001));
        assert_eq!(stored.user_delta_t_days(), Some(0.001));
        let non_finite = EphemerisContext::new().with_user_delta_t(Some(f64::NAN));
        assert_eq!(non_finite.user_delta_t_days().map(f64::is_nan), Some(true));
    }

    #[test]
    fn tidal_acceleration_automatic_and_effective_value() {
        // Default: automatic mode, effective value is the release default.
        let ctx = EphemerisContext::new();
        assert!(!ctx.tidal_acceleration_is_set());
        assert_eq!(ctx.effective_tidal_acceleration(), TIDAL_DEFAULT);
        // None resumes automatic; the TIDAL_AUTOMATIC sentinel does too.
        let cleared = EphemerisContext::new()
            .with_tidal_acceleration(Some(-25.85))
            .with_tidal_acceleration(None);
        assert!(!cleared.tidal_acceleration_is_set());
        let sentinel = EphemerisContext::new()
            .with_tidal_acceleration(Some(-25.85))
            .with_tidal_acceleration(Some(TIDAL_AUTOMATIC as f64));
        assert!(!sentinel.tidal_acceleration_is_set());
        assert_eq!(sentinel.effective_tidal_acceleration(), TIDAL_DEFAULT);
        // An explicit value is authoritative.
        let explicit = EphemerisContext::new().with_tidal_acceleration(Some(-25.85));
        assert!(explicit.tidal_acceleration_is_set());
        assert_eq!(explicit.tidal_acceleration_override(), Some(-25.85));
        assert_eq!(explicit.effective_tidal_acceleration(), -25.85);
    }

    #[test]
    fn delta_t_model_names_and_resolution() {
        use crate::time::DeltaTModel;
        // None defers; names normalize; unknown names fail loudly.
        let deferred = EphemerisContext::new()
            .with_delta_t_model(Some(DeltaTModel::EspenakMeeus))
            .with_delta_t_model_name(None)
            .unwrap();
        assert_eq!(deferred.delta_t_model_override(), None);
        let named = EphemerisContext::new()
            .with_delta_t_model_name(Some(" ESPENAK_MEEUS "))
            .unwrap();
        assert_eq!(
            named.delta_t_model_override(),
            Some(DeltaTModel::EspenakMeeus)
        );
        assert!(
            matches!(
                EphemerisContext::new().with_delta_t_model_name(Some("smh2017")),
                Err(Error::InvalidDeltaTModel { .. })
            ),
            "accepted unknown Delta-T model"
        );
        // Resolution: programmatic, env, TOML, `smh2016` default; unknown
        // configured spellings fall through silently.
        assert_eq!(
            resolve_delta_t_model(Some(DeltaTModel::EspenakMeeus), Some("smh2016"), None),
            DeltaTModel::EspenakMeeus
        );
        assert_eq!(
            resolve_delta_t_model(None, Some("espenak_meeus"), Some("smh2016")),
            DeltaTModel::EspenakMeeus
        );
        assert_eq!(
            resolve_delta_t_model(None, None, Some(" Espenak_Meeus ")),
            DeltaTModel::EspenakMeeus
        );
        assert_eq!(
            resolve_delta_t_model(None, Some("bogus"), Some("also-bogus")),
            DeltaTModel::Smh2016
        );
        assert_eq!(
            resolve_delta_t_model(None, None, None),
            DeltaTModel::Smh2016
        );
    }

    #[test]
    fn in_place_setters_match_builders() {
        use crate::time::DeltaTModel;
        // Observer set/clear round-trips through the same validation.
        let mut ctx = EphemerisContext::new();
        assert_eq!(ctx.observer(), None);
        let obs = Observer::new(13.4, 52.5, 34.0).unwrap();
        ctx.set_observer(obs);
        assert_eq!(ctx.observer(), Some(obs));
        ctx.clear_observer();
        assert_eq!(ctx.observer(), None);
        // Tier/mode names accept public spellings, refuse unknown ones
        // without mutating the stored value.
        let mut ctx = EphemerisContext::new();
        ctx.set_tier_name("extended").unwrap();
        assert_eq!(ctx.tier(), Tier::Extended);
        assert!(ctx.set_tier_name("ultra").is_err());
        assert_eq!(ctx.tier(), Tier::Extended);
        ctx.set_mode_name("leb").unwrap();
        assert_eq!(ctx.mode(), ProviderMode::SealedLeb);
        assert!(ctx.set_mode_name("nope").is_err());
        assert_eq!(ctx.mode(), ProviderMode::SealedLeb);
        // Delta-T legs: programmatic model, user override, IERS toggle.
        let mut ctx = EphemerisContext::new();
        ctx.set_delta_t_model_name(Some("espenak_meeus")).unwrap();
        assert_eq!(
            ctx.delta_t_model_override(),
            Some(DeltaTModel::EspenakMeeus)
        );
        assert!(ctx.set_delta_t_model_name(Some("smh2017")).is_err());
        assert_eq!(
            ctx.delta_t_model_override(),
            Some(DeltaTModel::EspenakMeeus)
        );
        ctx.set_delta_t_model_name(None).unwrap();
        assert_eq!(ctx.delta_t_model_override(), None);
        ctx.set_user_delta_t(Some(0.001));
        assert_eq!(ctx.user_delta_t_days(), Some(0.001));
        ctx.set_user_delta_t(None);
        assert_eq!(ctx.user_delta_t_days(), None);
        ctx.set_iers_delta_t(Some(true));
        assert_eq!(ctx.iers_delta_t_override(), Some(true));
        ctx.set_iers_delta_t(None);
        assert_eq!(ctx.iers_delta_t_override(), None);
        // Tidal-acceleration and lapse-rate legs: in-place setters match
        // the builders, sentinels clear back to the effective defaults.
        ctx.set_tidal_acceleration(Some(-25.85));
        assert_eq!(ctx.tidal_acceleration_override(), Some(-25.85));
        assert_eq!(ctx.effective_tidal_acceleration(), -25.85);
        ctx.set_tidal_acceleration(Some(TIDAL_AUTOMATIC as f64));
        assert_eq!(ctx.tidal_acceleration_override(), None);
        assert_eq!(ctx.effective_tidal_acceleration(), TIDAL_DEFAULT);
        ctx.set_tidal_acceleration(Some(-25.936));
        assert_eq!(ctx.tidal_acceleration_override(), Some(-25.936));
        ctx.set_tidal_acceleration(None);
        assert_eq!(ctx.tidal_acceleration_override(), None);
        ctx.set_lapse_rate(Some(0.005));
        assert_eq!(ctx.lapse_rate_override(), Some(0.005));
        assert_eq!(ctx.effective_lapse_rate(), 0.005);
        ctx.set_lapse_rate(None);
        assert_eq!(ctx.lapse_rate_override(), None);
        assert_eq!(ctx.effective_lapse_rate(), LAPSE_RATE_DEFAULT);
        // Network-policy legs: programmatic override and release-name
        // form; a refused spelling leaves the stored value untouched.
        let mut ctx = EphemerisContext::new();
        assert_eq!(ctx.network_policy_override(), None);
        ctx.set_network_policy(Some(NetworkPolicy::Sealed));
        assert_eq!(ctx.network_policy_override(), Some(NetworkPolicy::Sealed));
        ctx.set_network_policy(None);
        assert_eq!(ctx.network_policy_override(), None);
        ctx.set_network_policy_name(Some(" SEALED ")).unwrap();
        assert_eq!(ctx.network_policy_override(), Some(NetworkPolicy::Sealed));
        assert!(ctx.set_network_policy_name(Some("opne")).is_err());
        assert_eq!(ctx.network_policy_override(), Some(NetworkPolicy::Sealed));
        ctx.set_network_policy_name(None).unwrap();
        assert_eq!(ctx.network_policy_override(), None);
    }

    #[test]
    fn releasing_leb_images_preserves_selection_and_retained_snapshot() {
        let path = PathBuf::from("/data/retained.leb2");
        let site = Observer::new(12.5, 41.9, 35.0).unwrap();
        let mut ctx = EphemerisContext::new()
            .with_leb_file(Some(path.clone()))
            .with_observer(site)
            .with_mode(ProviderMode::SealedLeb)
            .with_network_policy(Some(NetworkPolicy::Sealed));
        ctx.attach_leb_file(bound_leb_file());
        let retained = ctx.clone();
        ctx.release_leb_images();
        assert!(ctx.leb_files().is_empty());
        assert_eq!(ctx.leb_file(), Some(&path));
        assert_eq!(ctx.observer(), Some(site));
        assert_eq!(ctx.mode(), ProviderMode::SealedLeb);
        assert_eq!(ctx.network_policy_override(), Some(NetworkPolicy::Sealed));
        assert_eq!(retained.leb_files().len(), 1);
        assert!(!retained.leb_files()[0].is_closed());
        ctx.release_leb_images();
        assert_eq!(ctx.leb_file(), Some(&path));
    }

    #[test]
    fn leb_file_setter_is_idempotent_and_invalidates_on_change() {
        let mut ctx = EphemerisContext::new();
        // Fresh context: clearing the unconfigured selection is a no-op.
        assert!(!ctx.set_leb_file(None));
        assert_eq!(ctx.leb_file(), None);

        // Configuring a file reports the change.
        assert!(ctx.set_leb_file(Some(PathBuf::from("/data/eph.leb"))));
        assert_eq!(ctx.leb_file(), Some(&PathBuf::from("/data/eph.leb")));
        // Re-setting the same file is an idempotent no-op.
        assert!(!ctx.set_leb_file(Some(PathBuf::from("/data/eph.leb"))));

        // Bind an image resolved from that file.
        ctx.attach_leb_file(bound_leb_file());
        assert_eq!(ctx.leb_files().len(), 1);

        // The unchanged selection keeps the bound store warm.
        assert!(!ctx.set_leb_file(Some(PathBuf::from("/data/eph.leb"))));
        assert_eq!(ctx.leb_files().len(), 1);

        // A changed selection drops the stores so stale data can never be
        // served; the new value is stored.
        assert!(ctx.set_leb_file(Some(PathBuf::from("/data/other.leb"))));
        assert!(ctx.leb_files().is_empty());
        assert_eq!(ctx.leb_file(), Some(&PathBuf::from("/data/other.leb")));

        // Disabling is a change; disabling again is not.
        assert!(ctx.set_leb_file(None));
        assert_eq!(ctx.leb_file(), None);
        assert!(!ctx.set_leb_file(None));

        // Builder routes the same logic.
        let ctx = EphemerisContext::new().with_leb_file(Some(PathBuf::from("/data/b.leb")));
        assert_eq!(ctx.leb_file(), Some(&PathBuf::from("/data/b.leb")));

        // close releases the bound images but preserves the file
        // selection.
        let mut ctx = EphemerisContext::new().with_leb_file(Some(PathBuf::from("/data/b.leb")));
        ctx.attach_leb_file(bound_leb_file());
        ctx.close();
        assert!(ctx.leb_files().is_empty());
        assert_eq!(ctx.leb_file(), Some(&PathBuf::from("/data/b.leb")));
    }

    #[test]
    fn ephemeris_file_override_yields_to_tier_and_close() {
        use crate::positions::select_de_kernel;
        let mut ctx = EphemerisContext::new();
        assert_eq!(ctx.ephemeris_file_override(), None);

        // Explicit name stored, no validation.
        ctx.set_ephemeris_file(Some("de441.bsp"));
        assert_eq!(ctx.ephemeris_file_override(), Some("de441.bsp"));
        // The override feeds the explicit leg of kernel selection ahead
        // of the tier kernel.
        assert_eq!(
            select_de_kernel(ctx.ephemeris_file_override(), None, None, Tier::Base),
            "de441.bsp"
        );

        // A tier change clears the override: tier-based selection resumes.
        ctx.set_tier_name("extended").unwrap();
        assert_eq!(ctx.ephemeris_file_override(), None);
        assert_eq!(
            select_de_kernel(ctx.ephemeris_file_override(), None, None, Tier::Base),
            "de440s.bsp"
        );

        // A refused tier name leaves a stored override untouched.
        ctx.set_ephemeris_file(Some("de441.bsp"));
        assert!(ctx.set_tier_name("ultra").is_err());
        assert_eq!(ctx.ephemeris_file_override(), Some("de441.bsp"));

        // close reverts to tier-based selection.
        ctx.close();
        assert_eq!(ctx.ephemeris_file_override(), None);

        // Clearing explicitly also resumes tier selection.
        ctx.set_ephemeris_file(Some("de441.bsp"));
        ctx.set_ephemeris_file(None);
        assert_eq!(ctx.ephemeris_file_override(), None);

        // Builder routes the same logic.
        let ctx = EphemerisContext::new().with_ephemeris_file(Some("de440s.bsp"));
        assert_eq!(ctx.ephemeris_file_override(), Some("de440s.bsp"));
    }

    #[test]
    fn spk_cache_dir_is_a_plain_store_cleared_by_close() {
        use crate::spk::resolve_spk_cache_dir;
        let mut ctx = EphemerisContext::new();
        assert_eq!(ctx.spk_cache_dir(), None);

        // Plain store: no validation, no invalidation, no changed-flag.
        ctx.set_spk_cache_dir(Some(PathBuf::from("/data/spk")));
        assert_eq!(ctx.spk_cache_dir(), Some(&PathBuf::from("/data/spk")));
        // The override feeds the explicit leg of the cache-dir chain.
        assert_eq!(
            resolve_spk_cache_dir(
                ctx.spk_cache_dir().and_then(|p| p.to_str()),
                Some("/env"),
                Some("/toml")
            ),
            Some("/data/spk".to_string())
        );

        // Re-setting the same value is a plain re-store, still visible.
        ctx.set_spk_cache_dir(Some(PathBuf::from("/data/spk")));
        assert_eq!(ctx.spk_cache_dir(), Some(&PathBuf::from("/data/spk")));

        // Clearing reverts to the environment/TOML/default chain.
        ctx.set_spk_cache_dir(None);
        assert_eq!(ctx.spk_cache_dir(), None);
        assert_eq!(
            resolve_spk_cache_dir(None, Some("/env"), Some("/toml")),
            Some("/env".to_string())
        );

        // close resets the override.
        ctx.set_spk_cache_dir(Some(PathBuf::from("/data/spk")));
        ctx.close();
        assert_eq!(ctx.spk_cache_dir(), None);

        // Builder routes the same logic.
        let ctx = EphemerisContext::new().with_spk_cache_dir(Some(PathBuf::from("/data/spk")));
        assert_eq!(ctx.spk_cache_dir(), Some(&PathBuf::from("/data/spk")));
    }

    #[test]
    fn spk_serving_toggles_default_to_release_values() {
        let ctx = EphemerisContext::new();
        assert_eq!(ctx.auto_spk_download_override(), None);
        assert!(ctx.effective_auto_spk_download());
        assert_eq!(ctx.spk_date_padding_days(), 0);
        assert_eq!(ctx.strict_precision_override(), None);
        assert!(ctx.effective_strict_precision());
        assert!(!ctx.effective_iers_delta_t());
    }

    #[test]
    fn spk_serving_toggles_store_and_clear() {
        let mut ctx = EphemerisContext::new();
        ctx.set_auto_spk_download(Some(false));
        assert_eq!(ctx.auto_spk_download_override(), Some(false));
        assert!(!ctx.effective_auto_spk_download());
        ctx.set_auto_spk_download(None);
        assert_eq!(ctx.auto_spk_download_override(), None);
        assert!(ctx.effective_auto_spk_download());

        ctx.set_spk_date_padding_days(365);
        assert_eq!(ctx.spk_date_padding_days(), 365);

        ctx.set_strict_precision(Some(false));
        assert_eq!(ctx.strict_precision_override(), Some(false));
        assert!(!ctx.effective_strict_precision());
        ctx.set_strict_precision(None);
        assert!(ctx.effective_strict_precision());

        // Builders route the same stores.
        let ctx = EphemerisContext::new()
            .with_auto_spk_download(Some(true))
            .with_spk_date_padding_days(30)
            .with_strict_precision(Some(false));
        assert_eq!(ctx.auto_spk_download_override(), Some(true));
        assert_eq!(ctx.spk_date_padding_days(), 30);
        assert!(!ctx.effective_strict_precision());
    }

    #[test]
    fn spk_toggle_resolvers_mirror_the_release_chains() {
        // Programmatic value wins over everything.
        assert!(resolve_auto_spk_download(
            Some(true),
            Some("0"),
            Some(false)
        ));
        assert!(!resolve_auto_spk_download(
            Some(false),
            Some("1"),
            Some(true)
        ));
        assert!(resolve_strict_precision(
            Some(true),
            Some("off"),
            Some(false)
        ));
        assert!(!resolve_strict_precision(
            Some(false),
            Some("on"),
            Some(true)
        ));
        // Every release env spelling, case/whitespace-insensitive; unknown
        // or blank words fall through to TOML, then to the `true` default.
        for word in ["1", "true", "yes", "on", "enabled", " TRUE ", "Yes"] {
            assert!(
                resolve_auto_spk_download(None, Some(word), None),
                "env {word:?}"
            );
            assert!(
                resolve_auto_spk_download(None, Some(word), Some(false)),
                "env {word:?}"
            );
            assert!(
                resolve_strict_precision(None, Some(word), None),
                "env {word:?}"
            );
        }
        for word in ["0", "false", "no", "off", "disabled", " OFF ", "No"] {
            assert!(
                !resolve_auto_spk_download(None, Some(word), None),
                "env {word:?}"
            );
            assert!(
                !resolve_auto_spk_download(None, Some(word), Some(true)),
                "env {word:?}"
            );
            assert!(
                !resolve_strict_precision(None, Some(word), None),
                "env {word:?}"
            );
        }
        for word in ["", "  ", "maybe", "2"] {
            assert!(
                resolve_auto_spk_download(None, Some(word), Some(true)),
                "env {word:?}"
            );
            assert!(
                !resolve_auto_spk_download(None, Some(word), Some(false)),
                "env {word:?}"
            );
            assert!(
                resolve_auto_spk_download(None, Some(word), None),
                "env {word:?}"
            );
            assert!(
                resolve_strict_precision(None, Some(word), None),
                "env {word:?}"
            );
        }
        assert!(resolve_auto_spk_download(None, None, Some(true)));
        assert!(!resolve_auto_spk_download(None, None, Some(false)));
        assert!(resolve_auto_spk_download(None, None, None));
        assert!(resolve_strict_precision(None, None, None));
    }

    #[test]
    fn iers_cache_dir_is_a_plain_store_preserved_by_close() {
        let mut ctx = EphemerisContext::new();
        assert_eq!(ctx.iers_cache_dir(), None);

        // Plain store: no validation, no directory creation, no
        // invalidation — any path (including the empty string) lands
        // exactly like the release setter.
        ctx.set_iers_cache_dir(Some(PathBuf::from("/data/iers")));
        assert_eq!(ctx.iers_cache_dir(), Some(&PathBuf::from("/data/iers")));
        ctx.set_iers_cache_dir(Some(PathBuf::from("")));
        assert_eq!(ctx.iers_cache_dir(), Some(&PathBuf::from("")));

        // Clearing reverts to the release default cache location.
        ctx.set_iers_cache_dir(None);
        assert_eq!(ctx.iers_cache_dir(), None);

        // reset_session leaves configuration untouched.
        ctx.set_iers_cache_dir(Some(PathBuf::from("/data/iers")));
        ctx.reset_session();
        assert_eq!(ctx.iers_cache_dir(), Some(&PathBuf::from("/data/iers")));

        // close preserves the override — unlike the SPK cache-dir leg, which
        // close resets.
        ctx.set_spk_cache_dir(Some(PathBuf::from("/data/spk")));
        ctx.close();
        assert_eq!(ctx.iers_cache_dir(), Some(&PathBuf::from("/data/iers")));
        assert_eq!(ctx.spk_cache_dir(), None);

        // Builder routes the same logic.
        let ctx = EphemerisContext::new().with_iers_cache_dir(Some(PathBuf::from("/data/iers")));
        assert_eq!(ctx.iers_cache_dir(), Some(&PathBuf::from("/data/iers")));
    }

    #[test]
    fn iers_auto_download_defaults_off_and_survives_close() {
        let mut ctx = EphemerisContext::new();
        assert_eq!(ctx.iers_auto_download_override(), None);
        assert!(!ctx.effective_iers_auto_download());

        ctx.set_iers_auto_download(Some(true));
        assert_eq!(ctx.iers_auto_download_override(), Some(true));
        assert!(ctx.effective_iers_auto_download());
        ctx.set_iers_auto_download(None);
        assert_eq!(ctx.iers_auto_download_override(), None);
        assert!(!ctx.effective_iers_auto_download());

        // Builders route the same store.
        let ctx = EphemerisContext::new().with_iers_auto_download(Some(false));
        assert_eq!(ctx.iers_auto_download_override(), Some(false));

        // close preserves the toggle — unlike the SPK toggle, which
        // close defers back to the chain.
        let mut ctx = EphemerisContext::new().with_iers_auto_download(Some(true));
        ctx.set_auto_spk_download(Some(false));
        ctx.close();
        assert_eq!(ctx.iers_auto_download_override(), Some(true));
        assert_eq!(ctx.auto_spk_download_override(), None);
    }

    #[test]
    fn iers_auto_download_resolver_mirrors_the_release_chain() {
        // Programmatic value wins over everything.
        assert!(resolve_iers_auto_download(
            Some(true),
            Some("0"),
            Some(false)
        ));
        assert!(!resolve_iers_auto_download(
            Some(false),
            Some("1"),
            Some(true)
        ));
        // Every release env spelling, case/whitespace-insensitive; unknown
        // or blank words fall through to TOML, then to the `false`
        // default.
        for word in ["1", "true", "yes", "on", "enabled", " TRUE ", "Yes"] {
            assert!(
                resolve_iers_auto_download(None, Some(word), None),
                "env {word:?}"
            );
            assert!(
                resolve_iers_auto_download(None, Some(word), Some(false)),
                "env {word:?}"
            );
        }
        for word in ["0", "false", "no", "off", "disabled", " OFF ", "No"] {
            assert!(
                !resolve_iers_auto_download(None, Some(word), None),
                "env {word:?}"
            );
            assert!(
                !resolve_iers_auto_download(None, Some(word), Some(true)),
                "env {word:?}"
            );
        }
        for word in ["", "  ", "maybe", "2"] {
            assert!(
                resolve_iers_auto_download(None, Some(word), Some(true)),
                "env {word:?}"
            );
            assert!(
                !resolve_iers_auto_download(None, Some(word), Some(false)),
                "env {word:?}"
            );
            assert!(
                !resolve_iers_auto_download(None, Some(word), None),
                "env {word:?}"
            );
        }
        assert!(resolve_iers_auto_download(None, None, Some(true)));
        assert!(!resolve_iers_auto_download(None, None, Some(false)));
        assert!(!resolve_iers_auto_download(None, None, None));
    }

    #[test]
    fn tier_table_words_match_the_release_tiers() {
        assert_eq!(TIER_TABLE.len(), 3);
        let rows = [
            (
                Tier::Base,
                "base",
                "de440s.bsp",
                "1850-01-01",
                "2150-01-01",
                "Modern usage (1850-2150), ~31 MB",
            ),
            (
                Tier::Medium,
                "medium",
                "de440.bsp",
                "1900-01-01",
                "2100-01-01",
                "General purpose (1550-2650), ~114 MB",
            ),
            (
                Tier::Extended,
                "extended",
                "de441.bsp",
                "1600-01-01",
                "2500-01-01",
                "Extended range (-13200 to +17191), ~3.1 GB",
            ),
        ];
        for (tier, name, kernel, start, end, description) in rows {
            let row = tier_row(tier);
            assert_eq!(
                (
                    row.name,
                    row.ephemeris_file,
                    row.spk_start,
                    row.spk_end,
                    row.description
                ),
                (name, kernel, start, end, description)
            );
            assert_eq!(
                [
                    tier_field(tier, 0),
                    tier_field(tier, 1),
                    tier_field(tier, 2),
                    tier_field(tier, 3),
                    tier_field(tier, 4)
                ],
                [name, kernel, start, end, description]
            );
        }
        // Unknown fields read as empty, never an error.
        for field in [-1, 5, 99] {
            assert_eq!(tier_field(Tier::Base, field), "");
        }
    }

    #[test]
    fn spk_date_ranges_are_cumulative_like_the_release() {
        // Base is its own window; medium unions with base (its literal
        // 1900-2100 is narrower than base 1850-2150); extended unions all.
        assert_eq!(
            spk_date_range_for_tier(Tier::Base),
            ("1850-01-01", "2150-01-01")
        );
        assert_eq!(
            spk_date_range_for_tier(Tier::Medium),
            ("1850-01-01", "2150-01-01")
        );
        assert_eq!(
            spk_date_range_for_tier(Tier::Extended),
            ("1600-01-01", "2500-01-01")
        );
    }

    #[test]
    fn close_resets_auto_spk_and_padding_but_preserves_strict() {
        let mut ctx = EphemerisContext::new();
        ctx.set_auto_spk_download(Some(false));
        ctx.set_spk_date_padding_days(30);
        ctx.set_strict_precision(Some(false));
        ctx.close();
        assert_eq!(ctx.auto_spk_download_override(), None);
        assert!(ctx.effective_auto_spk_download());
        assert_eq!(ctx.spk_date_padding_days(), 0);
        assert_eq!(ctx.strict_precision_override(), Some(false));
        // Reset preserves all three (session state only).
        ctx.set_auto_spk_download(Some(false));
        ctx.set_spk_date_padding_days(7);
        ctx.reset_session();
        assert_eq!(ctx.auto_spk_download_override(), Some(false));
        assert_eq!(ctx.spk_date_padding_days(), 7);
        assert_eq!(ctx.strict_precision_override(), Some(false));
    }

    #[test]
    fn iers_toggle_resolution_truth_table() {
        // Programmatic value wins over everything.
        assert!(resolve_iers_enabled(Some(true), Some("0"), Some(false)));
        assert!(!resolve_iers_enabled(Some(false), Some("1"), Some(true)));
        // Every release env spelling, case/whitespace-insensitive.
        for word in ["1", "true", "yes", "on", "enabled", " TRUE ", "Yes"] {
            assert!(resolve_iers_enabled(None, Some(word), None), "env {word:?}");
            assert!(
                resolve_iers_enabled(None, Some(word), Some(false)),
                "env {word:?} over TOML false"
            );
        }
        for word in ["0", "false", "no", "off", "disabled", " FALSE ", "No"] {
            assert!(
                !resolve_iers_enabled(None, Some(word), None),
                "env {word:?}"
            );
            assert!(
                !resolve_iers_enabled(None, Some(word), Some(true)),
                "env {word:?} over TOML true"
            );
        }
        // Blank or unrecognized env falls through to TOML, then `false`.
        assert!(resolve_iers_enabled(None, Some(""), Some(true)));
        assert!(!resolve_iers_enabled(None, Some(""), Some(false)));
        assert!(resolve_iers_enabled(None, Some("maybe"), Some(true)));
        assert!(!resolve_iers_enabled(None, Some("maybe"), None));
        assert!(!resolve_iers_enabled(None, None, None));
        // Context override plumbing.
        let ctx = EphemerisContext::new().with_iers_delta_t(Some(true));
        assert_eq!(ctx.iers_delta_t_override(), Some(true));
        assert!(resolve_iers_enabled(
            ctx.iers_delta_t_override(),
            None,
            None
        ));
    }

    #[test]
    fn lapse_rate_default_and_reset() {
        // Bit-pinned release default while unset.
        assert_eq!(LAPSE_RATE_DEFAULT.to_bits(), 0x3f7a_9fbe_76c8_b439);
        assert_eq!(
            EphemerisContext::new().effective_lapse_rate(),
            LAPSE_RATE_DEFAULT
        );
        // Explicit values store; None resets to the default.
        let custom = EphemerisContext::new().with_lapse_rate(Some(0.005));
        assert_eq!(custom.lapse_rate_override(), Some(0.005));
        assert_eq!(custom.effective_lapse_rate(), 0.005);
        let reset = custom.with_lapse_rate(None);
        assert_eq!(reset.lapse_rate_override(), None);
        assert_eq!(reset.effective_lapse_rate(), LAPSE_RATE_DEFAULT);
    }

    #[test]
    fn delta_t_config_bridge_carries_context_values() {
        use crate::time::{DeltaTModel, deltat_with_config};
        // The bridge forwards the override, the effective tidal
        // acceleration, and the caller-resolved model/IERS selections; the
        // vendored SMH-2016 provider is attached (DT-1) while IERS stays
        // `None` until table acquisition is wired.
        let ctx = EphemerisContext::new()
            .with_user_delta_t(Some(0.001))
            .with_tidal_acceleration(Some(-25.85));
        let cfg = ctx.delta_t_config(DeltaTModel::EspenakMeeus, false);
        assert_eq!(cfg.user_override_days, Some(0.001));
        assert_eq!(cfg.tidal_acceleration, -25.85);
        assert_eq!(cfg.model, DeltaTModel::EspenakMeeus);
        assert!(!cfg.iers_enabled);
        assert!(cfg.iers_lookup_seconds.is_none());
        assert!(cfg.smh2016_lookup_seconds.is_some());
        // End to end through the core: the override short-circuits
        // computation even for a non-finite epoch.
        assert_eq!(deltat_with_config(f64::NAN, &cfg).unwrap(), 0.001);
        // Without an override the effective (default) tidal acceleration
        // flows into the computed Espenak-Meeus leg.
        let plain = EphemerisContext::new();
        let cfg2 = plain.delta_t_config(DeltaTModel::EspenakMeeus, false);
        assert_eq!(cfg2.tidal_acceleration, TIDAL_DEFAULT);
        let expected = crate::time::deltat_espenak_meeus(2000.0 + (1.0 - 0.5) / 12.0);
        let via_ctx = deltat_with_config(2451545.0, &cfg2).unwrap();
        assert_eq!(via_ctx, (expected + 0.0) / 86400.0);
        // The default model serves end to end through the bridge (J2000
        // via the release construction, adjustment-free at the default
        // tidal acceleration).
        let cfg3 = plain.delta_t_config(DeltaTModel::Smh2016, false);
        let via_default = deltat_with_config(2451545.0, &cfg3).unwrap();
        assert_eq!(
            via_default.to_bits(),
            (63.82896039639765 / 86400.0f64).to_bits()
        );
    }

    #[test]
    fn new_error_messages_mirror_release_wording() {
        let model = Error::InvalidDeltaTModel {
            value: "smh2017".to_string(),
        };
        assert!(model.to_string().starts_with("Unknown delta_t model"));
        let mode = Error::InvalidCalcMode {
            value: "lebb".to_string(),
            source: CALC_MODE_ENV_VAR,
        };
        assert!(
            mode.to_string()
                .starts_with("Invalid calculation mode from LIBEPHEMERIS_MODE")
        );
        let tier = Error::InvalidTier {
            value: "ultra".to_string(),
        };
        assert!(tier.to_string().starts_with("Invalid tier:"));
    }

    #[test]
    fn config_overrides_stay_context_local() {
        // Builders are consuming: sibling contexts never observe each
        // other's Delta-T / refraction overrides (no globals).
        let a = EphemerisContext::new().with_user_delta_t(Some(0.002));
        let b = EphemerisContext::new()
            .with_tidal_acceleration(Some(-25.85))
            .with_lapse_rate(Some(0.008));
        assert_eq!(a.user_delta_t_days(), Some(0.002));
        assert_eq!(b.user_delta_t_days(), None);
        assert!(!a.tidal_acceleration_is_set());
        assert_eq!(b.effective_tidal_acceleration(), -25.85);
        assert_eq!(a.effective_lapse_rate(), LAPSE_RATE_DEFAULT);
        assert_eq!(b.effective_lapse_rate(), 0.008);
    }

    #[test]
    fn angles_cache_copies_in_and_out() {
        use std::collections::HashMap;
        // Release `get_angles_cache`/`set_angles_cache`: the setter copies
        // so external mutation cannot leak into the context.
        let mut ctx = EphemerisContext::new();
        assert!(ctx.angles_cache().is_empty());
        let mut external = HashMap::new();
        external.insert("asc".to_string(), 123.456);
        ctx.set_angles_cache(&external);
        assert_eq!(ctx.angles_cache().get("asc"), Some(&123.456));
        external.insert("asc".to_string(), 0.0);
        assert_eq!(ctx.angles_cache().get("asc"), Some(&123.456));
        ctx.clear_angles_cache();
        assert!(ctx.angles_cache().is_empty());
        // Clearing an empty cache is a no-op, never an error.
        ctx.clear_angles_cache();
    }

    #[test]
    fn reset_session_clears_session_state_and_preserves_configuration() {
        use crate::time::DeltaTModel;
        use std::collections::HashMap;
        // Release `reset_session`: per-calculation state goes, reusable
        // resources and configuration stay.
        let mut ctx = EphemerisContext::new()
            .with_observer(Observer::new(12.5, 41.9, 0.0).unwrap())
            .with_tier(Tier::Extended)
            .with_mode(ProviderMode::Auto)
            .with_paths(DataPaths {
                leb_dir: Some(PathBuf::from("/data/leb")),
                jpl_dir: Some(PathBuf::from("/data/jpl")),
            })
            .with_user_delta_t(Some(0.002))
            .with_tidal_acceleration(Some(-25.85))
            .with_delta_t_model(Some(DeltaTModel::EspenakMeeus))
            .with_iers_delta_t(Some(true))
            .with_lapse_rate(Some(0.008));
        let mut angles = HashMap::new();
        angles.insert("mc".to_string(), 10.0);
        ctx.set_angles_cache(&angles);
        assert!(ctx.observer().is_some());

        ctx.reset_session();

        // Session state cleared.
        assert_eq!(ctx.observer(), None);
        assert_eq!(ctx.sidereal(), SiderealConfig::Tropical);
        assert!(ctx.angles_cache().is_empty());
        // Configuration and overrides preserved.
        assert_eq!(ctx.tier(), Tier::Extended);
        assert_eq!(ctx.mode(), ProviderMode::Auto);
        assert_eq!(ctx.paths().leb_dir, Some(PathBuf::from("/data/leb")));
        assert_eq!(ctx.paths().jpl_dir, Some(PathBuf::from("/data/jpl")));
        assert_eq!(ctx.user_delta_t_days(), Some(0.002));
        assert_eq!(ctx.tidal_acceleration_override(), Some(-25.85));
        assert_eq!(
            ctx.delta_t_model_override(),
            Some(DeltaTModel::EspenakMeeus)
        );
        assert_eq!(ctx.iers_delta_t_override(), Some(true));
        assert_eq!(ctx.lapse_rate_override(), Some(0.008));
        // Reset is idempotent.
        ctx.reset_session();
        assert_eq!(ctx.tier(), Tier::Extended);
    }

    #[test]
    fn close_releases_session_state_and_keeps_process_configuration() {
        use crate::time::DeltaTModel;
        use std::collections::HashMap;
        // Release `_close_inner` preserved/omitted set: mode, LEB
        // selection, Delta-T model, and user Delta-T survive; session
        // state, tier, JPL path, and derived overrides reset.
        let mut ctx = EphemerisContext::new()
            .with_observer(Observer::new(12.5, 41.9, 0.0).unwrap())
            .with_tier(Tier::Extended)
            .with_mode(ProviderMode::HorizonsRemote)
            .with_paths(DataPaths {
                leb_dir: Some(PathBuf::from("/data/leb")),
                jpl_dir: Some(PathBuf::from("/data/jpl")),
            })
            .with_user_delta_t(Some(0.002))
            .with_tidal_acceleration(Some(-25.85))
            .with_delta_t_model(Some(DeltaTModel::EspenakMeeus))
            .with_iers_delta_t(Some(true))
            .with_lapse_rate(Some(0.008));
        let mut angles = HashMap::new();
        angles.insert("mc".to_string(), 10.0);
        ctx.set_angles_cache(&angles);

        ctx.close();

        // Session state released.
        assert_eq!(ctx.observer(), None);
        assert_eq!(ctx.sidereal(), SiderealConfig::Tropical);
        assert!(ctx.angles_cache().is_empty());
        // Tier, JPL path, and derived overrides reset to automatic.
        assert_eq!(ctx.tier(), Tier::Medium);
        assert_eq!(ctx.paths().jpl_dir, None);
        assert_eq!(ctx.tidal_acceleration_override(), None);
        assert_eq!(ctx.iers_delta_t_override(), None);
        assert_eq!(ctx.lapse_rate_override(), None);
        assert_eq!(ctx.effective_lapse_rate(), LAPSE_RATE_DEFAULT);
        // Process-level configuration preserved.
        assert_eq!(ctx.mode(), ProviderMode::HorizonsRemote);
        assert_eq!(ctx.paths().leb_dir, Some(PathBuf::from("/data/leb")));
        assert_eq!(
            ctx.delta_t_model_override(),
            Some(DeltaTModel::EspenakMeeus)
        );
        assert_eq!(ctx.user_delta_t_days(), Some(0.002));
    }

    #[test]
    fn close_is_context_local() {
        // No globals: closing one context never touches a sibling.
        let mut a = EphemerisContext::new()
            .with_observer(Observer::new(0.0, 51.5, 11.0).unwrap())
            .with_tier(Tier::Base);
        let b = EphemerisContext::new().with_tier(Tier::Extended);
        a.close();
        assert_eq!(a.observer(), None);
        assert_eq!(a.tier(), Tier::Medium);
        assert!(b.observer().is_none());
        assert_eq!(b.tier(), Tier::Extended);
    }

    // -- fixtures: tiny hand-packed images (decode-validity only;
    // evaluation is covered by the LEB/SPK family tests) -----------------

    fn push_u32(v: &mut Vec<u8>, x: u32) {
        v.extend_from_slice(&x.to_le_bytes());
    }

    fn push_i32(v: &mut Vec<u8>, x: i32) {
        v.extend_from_slice(&x.to_le_bytes());
    }

    fn push_u64(v: &mut Vec<u8>, x: u64) {
        v.extend_from_slice(&x.to_le_bytes());
    }

    fn push_f64(v: &mut Vec<u8>, x: f64) {
        v.extend_from_slice(&x.to_le_bytes());
    }

    /// Minimal decodable one-body LEB1 image, mirroring
    /// `leb::tests::minimal_leb1`: 64-byte header, one 24-byte section
    /// entry, one 52-byte body entry (id 0, ICRS-barycentric, 2 segments
    /// over JD 0..2), and a 96-byte zero blob.
    fn minimal_leb1_bytes() -> Vec<u8> {
        use crate::leb::{
            COORD_ICRS_BARY, HEADER_SIZE, LEB1_BODY_ENTRY_SIZE, MAGIC_LEB1, SECTION_DIR_SIZE,
        };
        let mut file = Vec::with_capacity(HEADER_SIZE);
        file.extend_from_slice(&MAGIC_LEB1);
        push_u32(&mut file, 1); // version
        push_u32(&mut file, 1); // section count
        push_u32(&mut file, 1); // body count
        push_f64(&mut file, 0.0); // jd_start
        push_f64(&mut file, 2.0); // jd_end
        push_f64(&mut file, 0.0); // generation epoch
        push_u32(&mut file, 0); // flags
        file.extend_from_slice(&[0u8; 20]);
        debug_assert_eq!(file.len(), HEADER_SIZE);
        let body_at = file.len() + SECTION_DIR_SIZE;
        let blob_at = body_at + LEB1_BODY_ENTRY_SIZE;
        push_u32(&mut file, 0); // section id
        push_u32(&mut file, 0); // reserved
        push_u64(&mut file, body_at as u64);
        push_u64(&mut file, LEB1_BODY_ENTRY_SIZE as u64);
        push_i32(&mut file, 0); // body id
        push_u32(&mut file, COORD_ICRS_BARY);
        push_u32(&mut file, 2); // segment count
        push_f64(&mut file, 0.0);
        push_f64(&mut file, 2.0);
        push_f64(&mut file, 1.0); // interval
        push_u32(&mut file, 1); // degree
        push_u32(&mut file, 3); // components
        push_u64(&mut file, blob_at as u64);
        file.extend_from_slice(&[0u8; 96]);
        debug_assert_eq!(file.len(), blob_at + 96);
        file
    }

    /// Minimal decodable SPK image: one 1024-byte DAF file record with an
    /// empty summary chain (`FWARD = BWARD = 0`). Layout mirrors
    /// `spk::tests::build_image` with no segment groups.
    fn empty_daf_bytes() -> Vec<u8> {
        use crate::spk::{RECORD_BYTES, SPK_FILE_ID, SPK_ND, SPK_NI};
        let mut img = vec![0u8; RECORD_BYTES];
        img[0..8].copy_from_slice(SPK_FILE_ID);
        img[8..12].copy_from_slice(&SPK_ND.to_le_bytes());
        img[12..16].copy_from_slice(&SPK_NI.to_le_bytes());
        img[16..76].fill(b' ');
        // Bytes 76..84 (FWARD/BWARD) stay zero: no summary records.
        let free = (RECORD_BYTES / 8 + 1) as i32;
        img[84..88].copy_from_slice(&free.to_le_bytes());
        img[88..96].copy_from_slice(b"LTL-IEEE");
        img
    }

    fn bound_leb_file() -> LebFile {
        LebFile::from_bytes("mem://ctx04.leb1".to_string(), minimal_leb1_bytes())
            .expect("hand fixture must decode")
    }

    #[test]
    fn leb_path_setter_is_idempotent_and_invalidates_on_change() {
        let mut ctx = EphemerisContext::new();
        // Fresh context: clearing the unconfigured path is a no-op.
        assert!(!ctx.set_leb_dir(None));
        assert!(ctx.leb_files().is_empty());

        // Configuring a path reports the change.
        assert!(ctx.set_leb_dir(Some(PathBuf::from("/data/leb"))));
        // Re-setting the same path is an idempotent no-op.
        assert!(!ctx.set_leb_dir(Some(PathBuf::from("/data/leb"))));

        // Bind an image resolved from that path.
        ctx.attach_leb_file(bound_leb_file());
        assert_eq!(ctx.leb_files().len(), 1);

        // The unchanged path keeps the bound store warm.
        assert!(!ctx.set_leb_dir(Some(PathBuf::from("/data/leb"))));
        assert_eq!(ctx.leb_files().len(), 1);
        assert!(!ctx.leb_files()[0].is_closed());

        // A changed path drops the LEB stores so stale data can never be
        // served; the new value is stored.
        assert!(ctx.set_leb_dir(Some(PathBuf::from("/data/leb2"))));
        assert!(ctx.leb_files().is_empty());
        assert_eq!(ctx.paths().leb_dir, Some(PathBuf::from("/data/leb2")));

        // Clearing back to unconfigured is a change; clearing again is not.
        assert!(ctx.set_leb_dir(None));
        assert!(ctx.paths().leb_dir.is_none());
        assert!(!ctx.set_leb_dir(None));
    }

    #[test]
    fn jpl_path_setter_is_idempotent_and_invalidates_on_change() {
        let mut ctx = EphemerisContext::new();
        assert!(!ctx.set_jpl_dir(None));
        assert!(ctx.set_jpl_dir(Some(PathBuf::from("/data/jpl"))));
        assert!(!ctx.set_jpl_dir(Some(PathBuf::from("/data/jpl"))));

        ctx.spk_registry_mut()
            .register_image("mem://ctx04.bsp", empty_daf_bytes())
            .expect("hand fixture must decode");
        assert_eq!(ctx.spk_registry().kernel_count(), 1);

        // The unchanged path keeps the kernel store warm.
        assert!(!ctx.set_jpl_dir(Some(PathBuf::from("/data/jpl"))));
        assert_eq!(ctx.spk_registry().kernel_count(), 1);

        // A changed path releases kernels and registrations together.
        assert!(ctx.set_jpl_dir(Some(PathBuf::from("/data/jpl2"))));
        assert!(ctx.spk_registry().is_empty());
        assert_eq!(ctx.spk_registry().kernel_count(), 0);
        assert_eq!(ctx.paths().jpl_dir, Some(PathBuf::from("/data/jpl2")));
    }

    #[test]
    fn set_paths_splits_invalidation_per_side() {
        let both = || DataPaths {
            leb_dir: Some(PathBuf::from("/data/leb")),
            jpl_dir: Some(PathBuf::from("/data/jpl")),
        };
        let mut ctx = EphemerisContext::new();
        // Fresh context: configuring both paths reports a change.
        assert!(ctx.set_paths(both()));
        // The unchanged pair is a no-op.
        assert!(!ctx.set_paths(both()));

        ctx.attach_leb_file(bound_leb_file());
        ctx.spk_registry_mut()
            .register_image("mem://ctx04.bsp", empty_daf_bytes())
            .expect("hand fixture must decode");

        // Unchanged pair: both stores stay warm.
        assert!(!ctx.set_paths(both()));
        assert_eq!(ctx.leb_files().len(), 1);
        assert_eq!(ctx.spk_registry().kernel_count(), 1);

        // Changing only the LEB side drops only the LEB stores.
        assert!(ctx.set_paths(DataPaths {
            leb_dir: Some(PathBuf::from("/data/other")),
            jpl_dir: Some(PathBuf::from("/data/jpl")),
        }));
        assert!(ctx.leb_files().is_empty());
        assert_eq!(ctx.spk_registry().kernel_count(), 1);

        // The consuming builder follows the same rule: an unchanged pair
        // keeps stores, a changed pair drops exactly its side.
        ctx.attach_leb_file(bound_leb_file());
        let ctx = ctx.with_paths(DataPaths {
            leb_dir: Some(PathBuf::from("/data/other")),
            jpl_dir: Some(PathBuf::from("/data/jpl")),
        });
        assert_eq!(ctx.leb_files().len(), 1);
        assert_eq!(ctx.spk_registry().kernel_count(), 1);
        let ctx = ctx.with_paths(DataPaths {
            leb_dir: Some(PathBuf::from("/data/other")),
            jpl_dir: Some(PathBuf::from("/data/elsewhere")),
        });
        assert_eq!(ctx.leb_files().len(), 1);
        assert!(ctx.spk_registry().is_empty());
    }

    #[test]
    fn reset_session_preserves_bound_stores_and_close_releases_them() {
        use std::collections::HashMap;
        let mut ctx = EphemerisContext::new()
            .with_observer(Observer::new(12.5, 41.9, 0.0).unwrap())
            .with_paths(DataPaths {
                leb_dir: Some(PathBuf::from("/data/leb")),
                jpl_dir: Some(PathBuf::from("/data/jpl")),
            });
        ctx.attach_leb_file(bound_leb_file());
        ctx.spk_registry_mut()
            .register_image("mem://ctx04.bsp", empty_daf_bytes())
            .expect("hand fixture must decode");
        let mut angles = HashMap::new();
        angles.insert("mc".to_string(), 10.0);
        ctx.set_angles_cache(&angles);

        // reset_session clears session state but keeps every bound store.
        ctx.reset_session();
        assert_eq!(ctx.observer(), None);
        assert!(ctx.angles_cache().is_empty());
        assert_eq!(ctx.leb_files().len(), 1);
        assert!(!ctx.leb_files()[0].is_closed());
        assert_eq!(ctx.spk_registry().kernel_count(), 1);

        // Repeated sessions over unchanged paths are no-ops: the stores
        // stay warm and usable, never churned.
        assert!(!ctx.set_leb_dir(Some(PathBuf::from("/data/leb"))));
        assert!(!ctx.set_jpl_dir(Some(PathBuf::from("/data/jpl"))));
        assert!(!ctx.set_paths(DataPaths {
            leb_dir: Some(PathBuf::from("/data/leb")),
            jpl_dir: Some(PathBuf::from("/data/jpl")),
        }));
        assert_eq!(ctx.leb_files().len(), 1);
        assert_eq!(ctx.spk_registry().kernel_count(), 1);

        // close releases session state, derived overrides, and every
        // bound store, while preserving the process-level configuration
        // (mode, LEB selection, Delta-T model/user overrides).
        ctx.close();
        assert!(ctx.leb_files().is_empty());
        assert!(ctx.spk_registry().is_empty());
        assert_eq!(ctx.paths().leb_dir, Some(PathBuf::from("/data/leb")));
        assert_eq!(ctx.paths().jpl_dir, None);
        // Close is idempotent over the stores: a second close is safe.
        ctx.close();
        assert!(ctx.leb_files().is_empty());
        assert!(ctx.spk_registry().is_empty());
    }

    #[test]
    fn path_invalidation_is_context_local() {
        // No globals: path changes and close on one context never touch
        // a sibling, even one sharing the same path values.
        let mut a = EphemerisContext::new().with_paths(DataPaths {
            leb_dir: Some(PathBuf::from("/data/leb")),
            jpl_dir: Some(PathBuf::from("/data/jpl")),
        });
        a.attach_leb_file(bound_leb_file());
        a.spk_registry_mut()
            .register_image("mem://ctx04.bsp", empty_daf_bytes())
            .expect("hand fixture must decode");
        let mut b = EphemerisContext::new().with_paths(DataPaths {
            leb_dir: Some(PathBuf::from("/data/leb")),
            jpl_dir: Some(PathBuf::from("/data/jpl")),
        });
        b.attach_leb_file(bound_leb_file());
        b.spk_registry_mut()
            .register_image("mem://ctx04-b.bsp", empty_daf_bytes())
            .expect("hand fixture must decode");

        assert!(a.set_leb_dir(Some(PathBuf::from("/data/other"))));
        assert!(a.leb_files().is_empty());
        assert_eq!(b.leb_files().len(), 1);
        assert_eq!(b.paths().leb_dir, Some(PathBuf::from("/data/leb")));

        // close releases a's stores while b keeps serving from its own.
        a.close();
        assert!(a.spk_registry().is_empty());
        assert_eq!(b.spk_registry().kernel_count(), 1);
        assert_eq!(b.leb_files().len(), 1);
        a.close();
        assert_eq!(b.spk_registry().kernel_count(), 1);
    }

    #[test]
    fn bound_stores_clone_deeply_and_compare_by_value() {
        // Contexts are independent owned values: cloning a context with
        // bound stores deep-copies the images, and equality compares the
        // stored bytes (a path change on the clone breaks equality).
        let mut ctx = EphemerisContext::new().with_paths(DataPaths {
            leb_dir: Some(PathBuf::from("/data/leb")),
            jpl_dir: None,
        });
        ctx.attach_leb_file(bound_leb_file());
        let clone = ctx.clone();
        assert_eq!(ctx, clone);
        let mut changed = clone.clone();
        assert!(changed.set_leb_dir(Some(PathBuf::from("/data/other"))));
        assert_ne!(ctx, changed);
    }

    #[test]
    fn cold_context_forks_share_provider_initialization_without_sharing_invalidation() {
        use crate::positions::TierPolicy;
        use std::sync::Barrier;

        let mut context = EphemerisContext::new();
        context.attach_leb_file(bound_leb_file());
        let forks: Vec<_> = (0..8).map(|_| context.clone()).collect();
        let barrier = Barrier::new(forks.len());
        let indexes = std::thread::scope(|scope| {
            let workers: Vec<_> = forks
                .iter()
                .map(|fork| {
                    let barrier = &barrier;
                    scope.spawn(move || {
                        barrier.wait();
                        let index = fork.leb_composite(TierPolicy::AttachOrder, false).unwrap();
                        assert_eq!(index.body_coverage(0), Some((0.0, 2.0)));
                        std::ptr::from_ref(index) as usize
                    })
                })
                .collect();
            workers
                .into_iter()
                .map(|worker| worker.join().unwrap())
                .collect::<Vec<_>>()
        });
        let original = context
            .leb_composite(TierPolicy::AttachOrder, false)
            .unwrap();
        assert!(
            indexes
                .iter()
                .all(|address| *address == std::ptr::from_ref(original) as usize)
        );

        let mut changed = forks[0].clone().with_tier(Tier::Extended);
        let rebuilt = changed
            .leb_composite(TierPolicy::AttachOrder, false)
            .unwrap();
        assert!(!std::ptr::eq(original, rebuilt));
        assert!(changed.set_leb_dir(Some(PathBuf::from("/data/changed"))));
        assert!(
            changed
                .leb_composite(TierPolicy::AttachOrder, false)
                .is_err()
        );
        assert_eq!(original.body_coverage(0), Some((0.0, 2.0)));
        assert_eq!(context.leb_files().len(), 1);
    }

    #[test]
    fn horizons_client_is_absent_by_default_and_round_trips() {
        use crate::horizons::{HORIZONS_DEFAULT_MAX_CACHE_SIZE, HORIZONS_DEFAULT_TIMEOUT_SECS};
        use std::num::NonZeroUsize;

        // Sealed-safe default: no client, no shared cache.
        assert!(EphemerisContext::new().horizons_client().is_none());

        let config = HorizonsConfig::new(
            5.0,
            NonZeroUsize::new(16).unwrap(),
            NonZeroUsize::new(2).unwrap(),
        )
        .unwrap();
        let ctx = EphemerisContext::new()
            .with_horizons_config(config)
            .unwrap();
        let client = ctx.horizons_client().expect("attached client");
        assert_eq!(client.timeout_secs(), 5.0);
        assert_eq!(client.max_cache_size(), 16);
        assert_eq!(client.max_workers(), 2);

        // Replacing the config drops the previous client and its cache.
        let ctx = ctx.with_horizons_config(HorizonsConfig::default()).unwrap();
        let client = ctx.horizons_client().expect("replaced client");
        assert_eq!(client.timeout_secs(), HORIZONS_DEFAULT_TIMEOUT_SECS);
        assert_eq!(client.max_cache_size(), HORIZONS_DEFAULT_MAX_CACHE_SIZE);
        assert_eq!(client.cache_len(), 0);

        // Detaching restores the default no-shared-cache behavior.
        assert!(ctx.without_horizons_client().horizons_client().is_none());
    }

    #[test]
    fn horizons_config_refusal_aborts_without_attachment() {
        use std::num::NonZeroUsize;

        // Timeout validation lives in the config value: no `HorizonsConfig`
        // exists for a bad timeout, so no builder call can attach a partial
        // client. A refused chain leaves the previously attached client
        // (here: none) untouched.
        for bad in [0.0, -2.5, f64::NAN, f64::INFINITY] {
            assert!(
                HorizonsConfig::new(
                    bad,
                    NonZeroUsize::new(16).unwrap(),
                    NonZeroUsize::new(2).unwrap(),
                )
                .is_err(),
                "{bad} must refuse"
            );
        }
        assert!(EphemerisContext::new().horizons_client().is_none());
    }

    #[test]
    fn horizons_snapshots_share_owner_and_close_retires_only_the_callers_cache() {
        use crate::horizons::{HorizonsTimeType, StateVector};
        use std::num::NonZeroUsize;

        fn attach_vector(ctx: &EphemerisContext) {
            let client = ctx.horizons_client().expect("attached client");
            let vector = StateVector::new(1.0, 2.0, 3.0, 0.1, 0.2, 0.3).unwrap();
            client
                .get_or_fetch_with("301", 2451545.0, "500@399", HorizonsTimeType::Tdb, || {
                    Ok(vector)
                })
                .unwrap();
        }

        let config = HorizonsConfig::new(
            5.0,
            NonZeroUsize::new(16).unwrap(),
            NonZeroUsize::new(2).unwrap(),
        )
        .unwrap();
        let ctx = EphemerisContext::new()
            .with_horizons_config(config)
            .unwrap();

        // Calculation snapshots explicitly retain the same bounded cache.
        let fork = ctx.clone();
        let retained = ctx.shared_horizons_client().unwrap();
        assert!(std::sync::Arc::ptr_eq(
            &retained,
            &fork.shared_horizons_client().unwrap()
        ));
        assert_eq!(ctx, fork);
        attach_vector(&ctx);
        assert_eq!(ctx.horizons_client().unwrap().cache_len(), 1);
        assert_eq!(fork.horizons_client().unwrap().cache_len(), 1);

        let mut ctx = ctx;
        ctx.reset_session();
        assert!(std::sync::Arc::ptr_eq(
            &retained,
            &ctx.shared_horizons_client().unwrap()
        ));
        assert_eq!(ctx.horizons_client().unwrap().cache_len(), 1);

        // Teardown publishes an empty cache without invalidating active users.
        ctx.close();
        let client = ctx.horizons_client().expect("config survives close");
        assert_eq!(client.timeout_secs(), 5.0);
        assert_eq!(client.cache_len(), 0);
        assert!(!std::sync::Arc::ptr_eq(
            &retained,
            &ctx.shared_horizons_client().unwrap()
        ));
        assert_eq!(retained.cache_len(), 1);
        assert_eq!(fork.horizons_client().unwrap().cache_len(), 1);
    }

    #[test]
    fn shared_reads_are_deterministic_across_scoped_threads() {
        // One inhabited context shared by reference across
        use crate::time::DeltaTModel;
        use std::collections::HashMap;
        let mut angles = HashMap::new();
        angles.insert("mc".to_string(), 10.0);
        angles.insert("asc".to_string(), 20.0);
        let mut ctx = EphemerisContext::new()
            .with_observer(Observer::new(12.5, 41.9, 0.0).unwrap())
            .with_tier(Tier::Extended)
            .with_mode(ProviderMode::Auto)
            .with_user_delta_t(Some(0.002))
            .with_delta_t_model(Some(DeltaTModel::EspenakMeeus))
            .with_iers_delta_t(Some(true));
        ctx.set_angles_cache(&angles);
        fn snapshot(ctx: &EphemerisContext) -> String {
            format!(
                "{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}|{:?}",
                ctx.observer(),
                ctx.sidereal(),
                ctx.sidereal_mode_full(),
                ctx.sidereal_bits(),
                ctx.tier(),
                ctx.mode(),
                ctx.paths(),
                ctx.angles_cache(),
                ctx.user_delta_t_days(),
                ctx.delta_t_model_override(),
                ctx.iers_delta_t_override(),
                ctx.network_policy_override(),
            )
        }
        let expected = snapshot(&ctx);
        std::thread::scope(|scope| {
            let handles: Vec<_> = (0..8)
                .map(|_| {
                    scope.spawn(|| {
                        let mut seen = String::new();
                        for _ in 0..50 {
                            seen = snapshot(&ctx);
                            assert_eq!(seen, expected);
                        }
                        seen
                    })
                })
                .collect();
            for handle in handles {
                assert_eq!(handle.join().unwrap(), expected);
            }
        });
    }
}

#[cfg(test)]
mod main_data_directory_contract_tests {
    use super::*;

    #[test]
    fn resolved_main_data_directory_is_store_only_and_survives_lifecycle() {
        let path = PathBuf::from("/not-created-by-context/explicit-data");
        let mut ctx = EphemerisContext::new().with_jpl_data_dir(Some(path.clone()));
        assert_eq!(ctx.jpl_data_dir(), Some(&path));
        ctx.reset_session();
        assert_eq!(ctx.jpl_data_dir(), Some(&path));
        ctx.close();
        assert_eq!(ctx.jpl_data_dir(), Some(&path));
        ctx.set_jpl_data_dir(None);
        assert_eq!(ctx.jpl_data_dir(), None);
    }

    #[test]
    fn resolved_data_directory_does_not_replace_supplied_registry() {
        let mut ctx = EphemerisContext::new();
        ctx.spk_registry_mut()
            .register_image("owned", crate::spk::moon_test_kernel_bytes())
            .unwrap();
        let registry = ctx.spk_registry().clone();
        ctx.set_jpl_data_dir(Some(PathBuf::from("/another-explicit-data")));
        assert_eq!(ctx.spk_registry(), &registry);
    }
}

#[cfg(test)]
#[path = "iers_context_tests.rs"]
mod iers_owned_tests;

#[cfg(test)]
mod calc_mode_normalization_tests {
    use super::*;

    #[test]
    fn unicode_whitespace_and_lowercase_preserve_four_names() {
        for whitespace in ['\t', '\n', '\u{0085}', '\u{00a0}', '\u{2007}', '\u{3000}'] {
            for (name, mode) in [
                ("AUTO", ProviderMode::Auto),
                ("SKYFIELD", ProviderMode::NativeJpl),
                ("LEB", ProviderMode::SealedLeb),
                ("HORIZONS", ProviderMode::HorizonsRemote),
            ] {
                let raw = format!("{whitespace}{name}{whitespace}");
                assert_eq!(ProviderMode::from_calc_mode_name(&raw), Some(mode));
                assert_eq!(
                    EphemerisContext::new().with_mode_name(&raw).unwrap().mode(),
                    mode
                );
            }
        }
        assert_eq!(
            ProviderMode::from_calc_mode_name("SKYFIELD"),
            Some(ProviderMode::NativeJpl)
        );
        for unsupported in [
            "ſkyfield",
            "SKYFİELD",
            "native_jpl",
            "jpl",
            "moşeph",
            "AUTO\u{200b}",
        ] {
            assert_eq!(ProviderMode::from_calc_mode_name(unsupported), None);
        }
    }

    #[test]
    fn unicode_blank_values_keep_priority_and_toml_fallback() {
        for separator in [' ', '\t', '\u{00a0}', '\u{3000}'] {
            let blank = separator.to_string();
            assert_eq!(
                resolve_calc_mode(None, Some(&blank), Some("LEB")),
                Ok(ProviderMode::SealedLeb)
            );
            assert_eq!(
                resolve_calc_mode(None, Some(&blank), Some(&blank)),
                Ok(ProviderMode::Auto)
            );
        }
        assert_eq!(
            resolve_calc_mode(Some(ProviderMode::SealedLeb), Some("bad"), Some("bad")),
            Ok(ProviderMode::SealedLeb)
        );
        assert_eq!(
            resolve_calc_mode(None, Some("SKYFIELD"), Some("LEB")),
            Ok(ProviderMode::NativeJpl)
        );
    }

    #[test]
    fn non_whitespace_controls_are_not_mode_separators() {
        for control in ['\u{001c}', '\u{001d}', '\u{001e}', '\u{001f}'] {
            assert!(
                ProviderMode::from_calc_mode_name(&format!("{control}auto{control}")).is_none()
            );
            assert!(resolve_calc_mode(None, Some(&control.to_string()), Some("leb")).is_err());
        }
    }

    #[test]
    fn invalid_name_keeps_raw_owned_error_fields_and_source() {
        let raw = "\u{001c}INVALID\u{001f}";
        let error = resolve_calc_mode(None, Some(raw), Some("LEB")).unwrap_err();
        assert!(
            matches!(error, Error::InvalidCalcMode { value, source } if value == raw && source == CALC_MODE_ENV_VAR)
        );
        let error = resolve_calc_mode(None, None, Some(raw)).unwrap_err();
        assert!(
            matches!(error, Error::InvalidCalcMode { value, source } if value == raw && source == "TOML key 'mode'")
        );
    }
}

#[cfg(test)]
mod live_clock_rate_tests {
    use super::*;
    #[test]
    fn native_clock_constant_override_precedes_missing_iers_and_nonfinite_epoch() {
        for seconds in [-600.0, 0.0, 37.125, 123.0] {
            let ctx = EphemerisContext::new().with_user_delta_t(Some(seconds / 86400.0));
            for jd in [f64::NAN, 2451545.0] {
                let value = ctx
                    .delta_t_value_and_rate_days(jd, DeltaTModel::Smh2016, true)
                    .unwrap();
                assert_eq!(value, (seconds / 86400.0, 0.0));
                assert!(ctx.iers_source().loaded_row_snapshot().is_none());
            }
        }
    }
    #[test]
    fn native_clock_value_is_exact_and_month_model_has_only_tidal_rate() {
        for model in [DeltaTModel::EspenakMeeus, DeltaTModel::Smh2016] {
            for jd in [2000000.25, 2430000.5, 2451545.125, 2480000.75] {
                let ctx = EphemerisContext::new().with_tidal_acceleration(Some(-26.0));
                let (value, rate) = ctx.delta_t_value_and_rate_days(jd, model, false).unwrap();
                assert_eq!(
                    value.to_bits(),
                    crate::time::deltat_with_config(jd, &ctx.delta_t_config(model, false))
                        .unwrap()
                        .to_bits()
                );
                assert!(rate.is_finite());
                if model == DeltaTModel::EspenakMeeus {
                    let expected = if jd >= crate::time::TID_ACC_REFERENCE_JD {
                        0.0
                    } else {
                        2.0 * crate::time::tid_acc_adjustment_coefficient()
                            * (-26.0 - crate::constants::TIDAL_DEFAULT)
                            * (jd - crate::time::TID_ACC_REFERENCE_JD)
                            / (36525.0 * 36525.0)
                            / 86400.0
                    };
                    assert_eq!(rate, expected);
                }
            }
        }
    }
    #[test]
    fn native_clock_iers_value_and_slope_share_loaded_owner() {
        let mut ctx = EphemerisContext::new();
        let rows = crate::iers::IersDataSnapshot::from_texts(
            None,
            None,
            Some("2000 1 1 63.0\n2000 2 1 66.1\n2000 3 1 71.9\n"),
        );
        let source = ctx.iers_source().with_rows(rows);
        ctx.attach_iers_source(source);
        for (jd, expected) in [(2451545.0, 0.1), (2451575.5, 0.2), (2451604.5, 0.2)] {
            let (value, rate) = ctx
                .delta_t_value_and_rate_days(jd, DeltaTModel::EspenakMeeus, true)
                .unwrap();
            assert_eq!(
                value.to_bits(),
                crate::time::deltat_with_config(
                    jd,
                    &ctx.delta_t_config(DeltaTModel::EspenakMeeus, true)
                )
                .unwrap()
                .to_bits()
            );
            assert!((rate - expected / 86400.0).abs() < 1e-19);
        }
    }
}
