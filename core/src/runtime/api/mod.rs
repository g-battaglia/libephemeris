//! Public calculation sessions over retained native providers.
//!
//! [`Session`] combines an immutable context snapshot, selected resources,
//! clocks and frame configuration. Provider errors retain native categories
//! and causes. Results use typed records, binary64 values and documented
//! units; bindings convert those values without supplying calculations.
//!
//! Start with [`SessionBuilder::isolated`] for explicit local data, or
//! [`crate::config::Runtime`] for mutable settings and lazy resource discovery.
//! A session has no mutable application settings: cloning it shares immutable
//! data owners, while [`Session::to_builder`] creates independently editable
//! configuration. Calculation methods return [`Result`];
//! use `?` to propagate an [`Error`] to the caller.
//!
//! # Data-independent example
//!
//! Ordinary tropical house geometry needs no planetary file. The epoch below
//! is a UT1 Julian day, not a UTC timestamp; site coordinates are north/east
//! degrees. `b'P'` is the Rust byte literal selecting Placidus.
//!
//! ```
//! use libephemeris_core::api::{SessionBuilder, Error};
//!
//! fn main() -> Result<(), Error> {
//!     let session = SessionBuilder::new().build()?;
//!     let houses = session.houses(2451545.0, 41.9, 12.5, b'P')?;
//!     println!("Ascendant: {} degrees", houses.ascendant());
//!     Ok(())
//! }
//! ```
//!
//! Position results and house results have different layouts; see [`Position`]
//! and [`Houses`]. In particular, position rates are per TT day, whereas house
//! rates differentiate the house query's UT1 Julian day.
pub use crate::constants::*;
use crate::context::{EphemerisContext, NetworkPolicy, Observer, ProviderMode};
pub use crate::errors::{Error, ErrorCategory};
use crate::time::DeltaTModel;
pub use crate::time::{Calendar, CalendarDate, revjul};
mod pure_calendar;
pub use crate::utils::{deg_midp, degnorm, difdeg2n, difdegn, difrad2n, rad_midp, radnorm};
pub use pure_calendar::{
    DateConversion, JulianDay, UtcDateTime, Weekday, date_conversion, day_of_week, julday,
    utc_time_zone,
};

mod advisory;
pub use advisory::{Advisory, AdvisoryDetail};
mod angles;
mod arabic_parts;
pub mod contrib;
pub use arabic_parts::ArabicPartsError;
mod assist_minor;
pub mod eclipse_visibility;
pub mod module_positions;
mod nbody_minor;
pub mod numeric_frames;
pub use angles::{cotrans, cotrans_sp};
pub use assist_minor::{AssistMinorEphemeris, MinorAssistAcquisition, assist_minor_sample};
pub use nbody_minor::{MinorNbodyDiagnostics, NativeMinorAssist, automatic_spk_interval};
mod besselian_utilities;
mod current_file;
pub use current_file::CurrentFileMetadata;
mod calendar;
pub use besselian_utilities::EclipseCycleAnchors;
mod eclipse_contact_utilities;
mod eclipse_remaining;
pub use eclipse_contact_utilities::SolarDurationError;
mod gauquelin;
mod houses;
mod public_houses;
pub use public_houses::{HouseBoundaryError, HouseLayoutWithFallback, HouseLayoutWithRates};
mod local_time;
mod lunar_eclipse;
mod lunar_global_search;
mod lunar_local_search;
pub use lunar_local_search::LunarLocalResult;
mod lunar_eclipse_scalars;
mod mean_lunar;
mod minor_bodies;
mod mutable_star_rows;
pub use mutable_star_rows::{
    AcquiredFixedStar, AcquiredStarRow, NativeStarInputs, StarRowAcquisition, StarRowBackend,
};
mod native_heliacal;
mod node_result;
mod orbital;
pub use node_result::NodesApsides;
mod pctr_context;
mod phenomena;
mod planet_metadata;
mod planetary_moons;
pub use planet_metadata::get_planet_name;
mod public_ayanamsha;
mod public_lunar;
mod vislim;
mod vislim_acquisition;
pub use vislim::{
    VisibilityObject, VisibilityObjectSelector, VisibilityPhenoAcquisition,
    VisibilityStarAcquisition,
};
pub use vislim_acquisition::{VisibilityAtmosphericAcquisition, VisibilityOpticalModel};
mod main_image_lifecycle;
mod positions;
pub use main_image_lifecycle::MainImageLifecycle;
pub use positions::PositionAcquisitions;
mod rise;
mod solar_circumstances;
mod solar_contacts;
mod solar_path;
pub use solar_path::EclipsePath;
mod solar_global_search;
mod solar_local_search;
pub use solar_local_search::SolarLocalResult;
mod star_rise;
mod stars;
pub use phenomena::{
    PhenoCallbacks, PhenoObservationUt1, PhenoPosition, PhenoPositionRetry, Phenomena,
};
pub use stars::{
    FixedStar, PreparedFixedStar, PreparedStarBatch, RetainedStarResult, prepare_fixed_star,
    prepare_star_batch,
};
mod stations;
pub use stations::{
    StationInformation, StationKind, StationMotionKind, StationOutput, StationQuery,
    StationRequest, station_identity, station_request,
};
mod crossings;
pub use crossings::{LongitudeCrossing, NodeCrossing, SearchDirection};
mod iers_lifecycle;
mod observer;
pub use observer::{Celestial, ECL2HOR, EQU2HOR, HOR2ECL, HOR2EQU, Horizontal};

/// Owned calculation components, returned flags and non-fatal advisories.
///
/// For an ordinary body, [`Session::calc`] and [`Session::calc_ut`] return:
///
/// | Indices | Default spherical output | With `FLG_EQUATORIAL` | With `FLG_XYZ` |
/// | --- | --- | --- | --- |
/// | 0, 1, 2 | Longitude (degrees), latitude (degrees), range (AU) | Right ascension (degrees), declination (degrees), range (AU) | x, y, z (AU) |
/// | 3, 4, 5 | Longitude, latitude and range rates | Right ascension, declination and range rates | vx, vy, vz (AU/TT day) |
///
/// Angular rates are degrees per TT day and range rates AU per TT day, even
/// for a UT1-input calculation. `FLG_RADIANS` converts spherical angles and
/// angular rates to radians; Cartesian output takes precedence over radians.
/// Request `FLG_SPEED` or `FLG_SPEED3` for ordinary position rates; without
/// either flag their three slots are zero. Flags select the origin and frame;
/// the default is an apparent geocentric ecliptic place of date. Accessors such
/// as [`Self::longitude`] read slots only: they do not convert units or frames.
///
/// **Exception:** `ECL_NUT` returns `[true_obliquity, mean_obliquity,
/// nutation_longitude, nutation_obliquity, 0, 0]`, in degrees. These are frame
/// angles, not a body position or rates; XYZ/RADIANS and speed requests do not
/// reshape that tuple. Its UT1 entry uses the independent bundled frame clock,
/// not a configured body/model Delta-T override.
///
/// Returned flags describe request normalization, not an authenticated data
/// identity. In particular, `FLG_SWIEPH` does not identify a foreign engine;
/// the configured provider selects the native calculation and its resources.
#[derive(Debug, Clone, PartialEq)]
pub struct Position {
    /// Six components in the layout and units documented on [`Position`].
    pub values: [f64; 6],
    /// Normalized returned request flags, not a foreign-engine data identity claim.
    pub returned_flags: i32,
    /// Configuration/calculation advisories; empty when none were emitted.
    pub diagnostic: String,
    /// Physical advisories emitted during this calculation.
    pub warnings: Vec<Advisory>,
    /// Native integrator messages, separate from configuration diagnostics.
    pub messages: Vec<String>,
}

/// Twelve house cusps and eight principal angles, without a padding element.
///
/// `cusps[0]` is house one and `cusps[11]` house twelve. Central positions are
/// degrees, or radians when the producing extended method selects
/// `FLG_RADIANS`. [`Session::houses`] always returns degrees. See
/// [`Session::houses_ex`] for flags, frames and the fixed-layout limitations.
#[derive(Debug, Clone, PartialEq)]
pub struct Houses {
    /// Ecliptic cusp longitudes, in the producing method's degree/radian units.
    pub cusps: [f64; 12],
    /// Slots 0–7: Ascendant, MC, ARMC, Vertex, equatorial Ascendant,
    /// Koch co-Ascendant, Munkasey co-Ascendant and Munkasey polar Ascendant.
    /// Units match the central cusp positions; ARMC is right ascension of the
    /// local meridian, not an ecliptic longitude.
    pub angles: [f64; 8],
}

/// House positions and their rates per UT1 day, without a padding element.
///
/// Central positions use degrees or requested radians, but **all rates remain
/// degrees per UT1 day**, including when `FLG_RADIANS` is selected. Layout and
/// angle meanings match [`Houses`]; see [`Session::houses_ex2`].
#[derive(Debug, Clone, PartialEq)]
pub struct HousesWithSpeeds {
    /// Twelve cusp longitudes, in degrees or requested radians.
    pub cusps: [f64; 12],
    /// Eight principal angles, matching [`Houses::angles`].
    pub angles: [f64; 8],
    /// Cusp rates in degrees per UT1 day.
    pub cusp_speeds: [f64; 12],
    /// Angle rates in degrees per UT1 day.
    pub angle_speeds: [f64; 8],
}

impl Position {
    /// Component zero: ecliptic longitude in degrees for default output.
    pub fn longitude(&self) -> f64 {
        self.values[0]
    }
    /// Component one: ecliptic latitude in degrees for default output.
    pub fn latitude(&self) -> f64 {
        self.values[1]
    }
    /// Component two: distance in AU for default output.
    pub fn distance(&self) -> f64 {
        self.values[2]
    }
}
impl Houses {
    /// Ascendant in the producing method's degree/radian units.
    pub fn ascendant(&self) -> f64 {
        self.angles[0]
    }
    /// Medium Coeli in the producing method's degree/radian units.
    pub fn mc(&self) -> f64 {
        self.angles[1]
    }
}
impl HousesWithSpeeds {
    /// Ascendant in the producing method's degree/radian units.
    pub fn ascendant(&self) -> f64 {
        self.angles[0]
    }
    /// Medium Coeli in the producing method's degree/radian units.
    pub fn mc(&self) -> f64 {
        self.angles[1]
    }
}

/// Thirty-six Gauquelin sectors and eight principal angles, without padding.
///
/// Central positions are degrees or requested radians; all rates are degrees
/// per UT1 day. See [`Session::houses_gauquelin`] for the producing operation.
#[derive(Debug, Clone, PartialEq)]
pub struct GauquelinSectors {
    /// Sector longitudes in degrees or requested radians; slot zero is sector one.
    pub sectors: [f64; 36],
    /// Principal angles in the same units as sectors, matching [`Houses::angles`].
    pub angles: [f64; 8],
    /// Sector rates in degrees per UT1 day.
    pub sector_speeds: [f64; 36],
    /// Angle rates in degrees per UT1 day.
    pub angle_speeds: [f64; 8],
}

/// Tropical-sidereal offset with flag and diagnostic projection.
#[derive(Debug, Clone, PartialEq)]
pub struct Ayanamsha {
    /// Offset in degrees.
    pub value: f64,
    /// Returned flags of the extended entry.
    pub returned_flags: i32,
    /// Non-fatal producer advisory, when present.
    pub diagnostic: String,
}

/// Immutable calculation configuration over shared images; safe to clone/share.
/// Advanced provider and configuration methods remain in the native SDK.
#[derive(Debug, Clone)]
pub struct Session {
    ctx: EphemerisContext,
    // Clones share the same immutable context. Every independent build gets a
    // fresh non-reusable token, with no owner retained in the numerical memo.
    owned_apparent_id: Option<std::num::NonZeroU64>,
    diagnostics: Vec<String>,
    star_main: std::sync::Arc<crate::positions::star_main::NativeStarMainCache>,
}

/// Reusable session builder with independent settings and shared data owners.
///
/// Methods taking `&mut self` update this builder and return it for chaining;
/// [`Self::build`] borrows it, so it can be reused. Numeric builder settings
/// are checked at build time. [`Self::leb_file`] instead opens and checks its
/// file immediately; building does not reopen it or prove coverage for every
/// later body/epoch query.
#[derive(Debug, Clone)]
pub struct SessionBuilder {
    ctx: EphemerisContext,
    topo: Option<(f64, f64, f64)>,
    sidereal: Option<(i32, f64, f64)>,
    diagnostics: Vec<String>,
}

impl Default for SessionBuilder {
    fn default() -> Self {
        Self::isolated()
    }
}
impl SessionBuilder {
    /// Start with the same explicit local defaults as [`Self::isolated`].
    pub fn new() -> Self {
        Self::isolated()
    }

    /// Start with explicit local settings, without resolving application config.
    ///
    /// Selects sealed LEB mode, sealed networking, SMH2016 body/model Delta-T
    /// with observed IERS disabled, an observer at longitude/latitude/height
    /// `(0°, 0°, 0 m)`, Fagan/Bradley sidereal selection and automatic tidal
    /// acceleration/Delta-T. Construction opens no file. Bind files explicitly
    /// for operations that need body coefficients; model-only operations such
    /// as ordinary tropical houses and `ECL_NUT` need no planetary file.
    pub fn isolated() -> Self {
        Self {
            ctx: EphemerisContext::new()
                .with_mode(ProviderMode::SealedLeb)
                .with_network_policy(Some(NetworkPolicy::Sealed))
                .with_delta_t_model(Some(DeltaTModel::Smh2016))
                .with_iers_delta_t(Some(false)),
            topo: Some((0.0, 0.0, 0.0)),
            sidereal: Some((SIDM_FAGAN_BRADLEY, 0.0, 0.0)),
            diagnostics: Vec::new(),
        }
    }

    /// Retain an explicitly configured context, sharing its immutable images.
    /// Environment/TOML resolution stays with the caller, never calculation.
    pub fn from_context(ctx: EphemerisContext) -> Self {
        Self {
            ctx,
            topo: None,
            sidereal: None,
            diagnostics: Vec::new(),
        }
    }

    /// Bind an existing reader/cache without re-reading its pathname.
    pub fn leb_image(&mut self, image: crate::leb::LebFile) -> &mut Self {
        self.ctx.attach_leb_file(image);
        self
    }

    /// Open and attach one explicit LEB1/LEB2 file immediately.
    ///
    /// `path` is a local pathname. The context's LEB read policy controls its
    /// storage and resource limits. A recognized filename supplies the data
    /// tier; an extension alone does not. Later queries still check actual body
    /// and epoch coverage. Attached immutable storage is shared, not copied.
    ///
    /// # Errors
    ///
    /// Returns file-access, malformed-container or resource-policy errors from
    /// the reader. This method does not download a replacement file.
    pub fn leb_file(&mut self, path: &std::path::Path) -> Result<&mut Self, Error> {
        let file = self.ctx.open_leb_file(path)?;
        if let Some(tier) =
            crate::leb::file_tier(path.file_name().and_then(|s| s.to_str()).unwrap_or(""))
        {
            self.ctx = self.ctx.clone().with_tier(match tier {
                crate::leb::LebTier::Base => crate::context::Tier::Base,
                crate::leb::LebTier::Medium => crate::context::Tier::Medium,
                crate::leb::LebTier::Extended => crate::context::Tier::Extended,
            });
        }
        Ok(self.leb_image(file))
    }

    /// Observer: east longitude/latitude degrees, height metres.
    pub fn topo(&mut self, lon: f64, lat: f64, height: f64) -> &mut Self {
        self.topo = Some((lon, lat, height));
        self
    }
    /// Set the sidereal mode and USER defining epoch/offset.
    ///
    /// `mode` combines a `SIDM_*` identifier with supported `SIDBIT_*` flags.
    /// `t0` is a Julian day, TT unless `SIDBIT_USER_UT` selects UT1;
    /// `ayan_t0` is degrees. These values are used by USER mode. Validation and
    /// normalization advisories are deferred to [`Self::build`].
    pub fn sid_mode(&mut self, mode: i32, t0: f64, ayan_t0: f64) -> &mut Self {
        self.sidereal = Some((mode, t0, ayan_t0));
        self.diagnostics.clear();
        self
    }
    /// Override lunar tidal acceleration in arcseconds per Julian century squared.
    /// `TIDAL_AUTOMATIC` restores computed selection; nonfinite stored values
    /// are rejected by [`Self::build`].
    pub fn tid_acc(&mut self, acc: f64) -> &mut Self {
        self.ctx.set_tidal_acceleration(Some(acc));
        self
    }
    /// Override model TT−UT1 in days, not seconds.
    /// `DELTAT_AUTOMATIC` resumes computed values. Independent frame/observer
    /// clocks are not replaced by this setting. Nonfinite stored values are
    /// rejected by [`Self::build`].
    pub fn delta_t_override(&mut self, days: f64) -> &mut Self {
        self.ctx.set_user_delta_t(Some(days));
        self
    }
    /// Resume configured computed Delta-T without changing its model.
    pub fn clear_delta_t_override(&mut self) -> &mut Self {
        self.ctx.set_user_delta_t(None);
        self
    }
    /// Atmospheric lapse rate in K/m.
    pub fn lapse_rate(&mut self, rate: f64) -> &mut Self {
        self.ctx.set_lapse_rate(Some(rate));
        self
    }

    /// Validate builder settings and create a session sharing attached images.
    ///
    /// Performs no file I/O or position calculation and does not consume this
    /// builder. Sidereal normalization advisories remain available through
    /// [`Session::configuration_diagnostics`]. No body coverage is certified.
    ///
    /// # Errors
    ///
    /// Returns [`Error::NonFiniteArgument`] for nonfinite Delta-T, tidal,
    /// lapse-rate or explicitly supplied sidereal epoch/offset settings, and
    /// [`Error::InvalidObserver`] for invalid explicitly supplied site values.
    /// Already configured contexts can contain settings not checked by this
    /// builder; later consuming operations retain their own validation.
    pub fn build(&self) -> Result<Session, Error> {
        let mut ctx = self.ctx.clone();
        let mut diagnostics = self.diagnostics.clone();
        for value in [
            ctx.user_delta_t_days(),
            ctx.tidal_acceleration_override(),
            ctx.lapse_rate_override(),
        ]
        .into_iter()
        .flatten()
        {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        if let Some((lon, lat, height)) = self.topo {
            ctx.set_observer(Observer::new(lon, lat, height)?);
        }
        if let Some((mode, t0, ayan)) = self.sidereal {
            for value in [t0, ayan] {
                if !value.is_finite() {
                    return Err(Error::NonFiniteArgument { value });
                }
            }
            diagnostics.extend(
                ctx.set_sidereal_mode(mode, t0, ayan)
                    .into_iter()
                    .map(|w| w.to_string()),
            );
        }
        Ok(Session {
            ctx,
            diagnostics,
            star_main: Default::default(),
            owned_apparent_id: crate::pure_cache::new_owned_apparent_identity(),
        })
    }
}

impl Session {
    /// Own an already configured context without file I/O or data duplication.
    pub fn from_context(ctx: EphemerisContext) -> Self {
        Self {
            ctx,
            diagnostics: Vec::new(),
            star_main: Default::default(),
            owned_apparent_id: crate::pure_cache::new_owned_apparent_identity(),
        }
    }
    /// Reconfigure independently, retaining explicitly bound image/cache owners.
    /// The new Session starts a separate lazy stellar main selection; snapshots
    /// already acquired by this Session and its clones remain valid.
    pub fn to_builder(&self) -> SessionBuilder {
        let mut builder = SessionBuilder::from_context(self.ctx.clone());
        builder.diagnostics.clone_from(&self.diagnostics);
        builder
    }
    /// Configuration normalization advisories; calculations do not erase them.
    pub fn configuration_diagnostics(&self) -> &[String] {
        &self.diagnostics
    }
    /// Bound immutable context, for explicit producer metadata/extensions.
    pub fn context(&self) -> &EphemerisContext {
        &self.ctx
    }
    /// Resolve the context's explicit model, defaulting without environment I/O.
    fn delta_t_model(&self) -> DeltaTModel {
        self.ctx
            .delta_t_model_override()
            .unwrap_or(DeltaTModel::Smh2016)
    }

    /// Borrow the bound data with the same time policy used by `deltat`.
    fn with_serving_stack<'c, R>(
        &'c self,
        f: impl FnOnce(&crate::positions::context_serving::ServingStack<'_, 'c>) -> Result<R, Error>,
    ) -> Result<R, Error> {
        crate::positions::context_serving::with_sealed_serving_stack(
            &self.ctx,
            self.delta_t_model(),
            self.ctx.effective_iers_delta_t(),
            false,
            f,
        )
    }

    /// TT−UT in days at UT; explicit context model/overrides, no environment read.
    pub fn deltat(&self, jd_ut: f64) -> Result<f64, Error> {
        crate::time::deltat_with_config(
            jd_ut,
            &self
                .ctx
                .delta_t_config(self.delta_t_model(), self.ctx.effective_iers_delta_t()),
        )
    }
}

pub use solar_path::CapturedSolarPathWidth;

mod numeric_occultation;
pub use numeric_occultation::numeric_occultation_target;

mod numeric_occult_global;

mod numeric_occult_local;
pub use numeric_occult_local::NumericOccultLocalResult;

mod stellar_occultation;
pub use stellar_occultation::{
    StellarOccultationAcquisition, StellarOccultationHorizontal, StellarOccultationPair,
    StellarOccultationRise,
};

#[cfg(test)]
mod owned_session_identity_tests {
    use super::*;

    #[test]
    fn clone_shares_immutable_producer_but_independent_builds_do_not() {
        let first = SessionBuilder::new().build().unwrap();
        let cloned = first.clone();
        assert_eq!(cloned.owned_apparent_id, first.owned_apparent_id);
        let rebuilt = first.to_builder().build().unwrap();
        assert_ne!(rebuilt.owned_apparent_id, first.owned_apparent_id);
        let from_context = Session::from_context(first.context().clone());
        assert_ne!(from_context.owned_apparent_id, first.owned_apparent_id);
    }
}

#[cfg(test)]
mod owned_apparent_session_tests;
