//! Fixed-star catalog lookup, astrometric propagation and apparent positions.
//!
//! Catalog positions use ICRS axes at J2000.0. Right ascension and declination
//! are degrees, proper motions arcseconds/year (RA includes cos(declination)),
//! parallax milliarcseconds and radial velocity km/s. The generated catalog
//! contains the astrometric rows; this module does not modify them.
//!
//! Names use native Unicode whitespace and case mapping. The two public search
//! families retain their catalog name/nomenclature order, aliases, Bayer and
//! Flamsteed forms, and wildcard rules. Exact names and aliases precede the
//! leading Unicode decimal index grammar. Sequential indices use bounded
//! `u64` arithmetic; an overflow is an ordinary catalog miss.
//!
//! Position calculations use binary64 at TT Julian days. Solar-system states
//! and observer geometry come from explicit providers; returned angles are
//! degrees, distances AU and velocities the corresponding units per day.
//! Missing catalog entries, unavailable data and numerical failures are typed
//! errors. Apparent-place formulas and flag domains are documented per entry.
//!
//! Private children separate catalog motion (`motion`), observer-relative
//! corrections and sample rates (`apparent`), dynamic sidereal anchor
//! ownership/sampling (`anchors`), and output-frame/flag projection (`flags`).
//! Public records, acquisition contracts and request serving remain here;
//! child entry points retain facade re-exports.

//! Catalog data and lookup grammar are generated from the self-contained
//! `tools/generators/data/star_catalog.json` snapshot. Its lineage identifies
//! the original project release blobs and separately qualified photometry;
//! see `docs/sources.md` for regeneration and `docs/limitations.md` for
//! magnitude limits.
//!
//! # Example: resolve a star without acquiring ephemeris data
//!
//! ```
//! use libephemeris_core::stars;
//!
//! let star = stars::resolve_star_ref("Spica").expect("built-in catalogue entry");
//! let photometry = stars::star_photometry(star.star_id).expect("catalogue photometry");
//! println!("{}: {} mag in {:?}", star.canonical, photometry.value, photometry.band);
//! // This is a static catalogue measurement, not the brightness at an epoch.
//! ```

use crate::errors::Error;
use crate::star_catalog_tables::STAR_PHOTOMETRY;
#[cfg(test)]
use crate::star_catalog_tables::{NOMEN_SORTED_IDS, STAR_ALIASES};

#[path = "stars/catalog.rs"]
mod catalog;
pub use catalog::{
    canonical_star_name, fixed_star_data, fixstar_mag, fixstar2_mag, hip_for_star_name,
    list_fixed_stars, parse_bayer_designation, parse_flamsteed_designation, resolve_star_name,
    resolve_star_ref, resolve_star2, star_aliases, star_entry_by_id, star_magnitude,
    star_name_to_hip, star_photometry, star_search_type,
};
#[cfg(test)]
use catalog::{nomen_stripped, normalize_phonetic, ref_star_key};

/// Fixed-star astrometric data.
///
/// ICRS J2000.0 epoch. `parallax_mas` of `0.0` means unknown parallax;
/// radial velocity is positive for receding stars. Instances come only from
/// the generated catalog; the serializer checks that measurements are finite.
///
/// An unknown parallax is represented at [`STAR_PARALLAX_FLOOR_MAS`]
/// (one gigaparsec); distance handling lives with the constant, not
/// this struct.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StarData {
    ra_j2000_deg: f64,
    dec_j2000_deg: f64,
    pm_ra_arcsec_yr: f64,
    pm_dec_arcsec_yr: f64,
    parallax_mas: f64,
    radial_velocity_km_s: f64,
}

impl StarData {
    /// Row constructor for the generated table (field order is `STAR_ROWS` tail layout).
    pub(crate) const fn from_row(
        ra_j2000_deg: f64,
        dec_j2000_deg: f64,
        pm_ra_arcsec_yr: f64,
        pm_dec_arcsec_yr: f64,
        parallax_mas: f64,
        radial_velocity_km_s: f64,
    ) -> Self {
        Self {
            ra_j2000_deg,
            dec_j2000_deg,
            pm_ra_arcsec_yr,
            pm_dec_arcsec_yr,
            parallax_mas,
            radial_velocity_km_s,
        }
    }

    /// Right Ascension at J2000.0, degrees.
    pub fn ra_j2000_deg(&self) -> f64 {
        self.ra_j2000_deg
    }
    /// Declination at J2000.0, degrees.
    pub fn dec_j2000_deg(&self) -> f64 {
        self.dec_j2000_deg
    }
    /// Proper motion in RA, arcsec/yr (includes the cos(dec) factor).
    pub fn pm_ra_arcsec_yr(&self) -> f64 {
        self.pm_ra_arcsec_yr
    }
    /// Proper motion in Dec, arcsec/yr.
    pub fn pm_dec_arcsec_yr(&self) -> f64 {
        self.pm_dec_arcsec_yr
    }
    /// Trigonometric parallax, milliarcseconds (`0.0` = unknown).
    pub fn parallax_mas(&self) -> f64 {
        self.parallax_mas
    }
    /// Catalog parallax in arcseconds; a zero parallax remains zero.
    pub fn parallax_arcsec(&self) -> f64 {
        self.parallax_mas / 1000.0
    }

    /// Radial velocity, km/s (positive = receding).
    pub fn radial_velocity_km_s(&self) -> f64 {
        self.radial_velocity_km_s
    }
}

/// Passband of a catalogue magnitude; no passband conversion is implied.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhotometricBand {
    /// The measurement's passband has not been established.
    Unspecified,
    /// Johnson V, distinct from the Hipparcos instrumental passband.
    JohnsonV,
    /// Hipparcos Hp instrumental passband.
    HipparcosHp,
}

/// Published catalogue field identifying a magnitude's provenance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhotometricCatalogueField {
    /// An inherited scalar without an individually established source field.
    Unqualified,
    /// ESA (1997), I/239/hip_main, H5 (`Vmag`, Johnson V).
    HipparcosMainVmag,
}

/// Reduction statistic explicitly documented for a magnitude.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhotometricStatistic {
    /// No mean, median or other reduction statistic has been established.
    Unspecified,
    /// A source explicitly identifies the reported value as a median.
    Median,
}

/// Source-reported variability classification, independent of the magnitude.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PhotometricVariability {
    /// No variability classification has been established; this is not constant.
    Unknown,
    /// The source reports no detected variability.
    Constant,
    /// The source reports periodic variability.
    Periodic,
    /// The source reports variability without a solved light curve.
    Unsolved,
}

/// A static catalogue magnitude and its explicitly established qualifications.
///
/// Values are magnitudes in `band`, not instantaneous brightness estimates.
/// Unspecified fields must not be inferred from the astrometric epoch, star
/// identity or numerical value. Only individually annotated catalogue entries
/// carry a qualified source; the other entries remain unqualified.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StarPhotometry {
    /// Catalogue magnitude in `band`, without light-curve evaluation.
    pub value: f64,
    /// Measurement passband, or explicitly unspecified.
    pub band: PhotometricBand,
    /// Published catalogue field, or an unqualified inherited scalar.
    pub catalogue_field: PhotometricCatalogueField,
    /// Source-defined reduction statistic; H5 does not specify one.
    pub statistic: PhotometricStatistic,
    /// Verbatim source epoch label, if supplied; no time-scale conversion is made.
    ///
    /// The currently qualified H5 entries have no established measurement epoch.
    pub epoch: Option<&'static str>,
    /// Source-reported variability, independently of any epoch or passband.
    pub variability: PhotometricVariability,
    /// Superseded catalogue scalar with no established band, statistic or epoch.
    ///
    /// This records history only and is never used in calculations.
    pub historical_unqualified_value: Option<f64>,
}

impl StarPhotometry {
    /// Preserve a supplied magnitude without asserting a band, source or epoch.
    ///
    /// This value-record constructor performs no validation or conversion;
    /// variability remains unknown, and no historical measurement is inferred.
    pub const fn unqualified(value: f64) -> Self {
        Self {
            value,
            band: PhotometricBand::Unspecified,
            catalogue_field: PhotometricCatalogueField::Unqualified,
            statistic: PhotometricStatistic::Unspecified,
            epoch: None,
            variability: PhotometricVariability::Unknown,
            historical_unqualified_value: None,
        }
    }
}

/// Extended catalog entry.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StarCatalogEntry {
    id: i32,
    name: &'static str,
    nomenclature: &'static str,
    hip_number: i32,
    data: StarData,
    magnitude: f64,
}

impl StarCatalogEntry {
    /// Row constructor for the generated table.
    #[allow(clippy::too_many_arguments)]
    pub(crate) const fn from_row(
        id: i32,
        name: &'static str,
        nomenclature: &'static str,
        hip_number: i32,
        ra_j2000_deg: f64,
        dec_j2000_deg: f64,
        pm_ra_arcsec_yr: f64,
        pm_dec_arcsec_yr: f64,
        parallax_mas: f64,
        radial_velocity_km_s: f64,
        magnitude: f64,
    ) -> Self {
        Self {
            id,
            name,
            nomenclature,
            hip_number,
            data: StarData::from_row(
                ra_j2000_deg,
                dec_j2000_deg,
                pm_ra_arcsec_yr,
                pm_dec_arcsec_yr,
                parallax_mas,
                radial_velocity_km_s,
            ),
            magnitude,
        }
    }

    /// Internal star ID (a `FIXSTAR_OFFSET`-based constant).
    pub fn id(&self) -> i32 {
        self.id
    }
    /// Traditional star name (e.g. `"Regulus"`).
    pub fn name(&self) -> &'static str {
        self.name
    }
    /// Bayer/Flamsteed designation (e.g. `"alLeo"`).
    pub fn nomenclature(&self) -> &'static str {
        self.nomenclature
    }
    /// Hipparcos catalog number.
    pub fn hip_number(&self) -> i32 {
        self.hip_number
    }
    /// Astrometric data for position calculation.
    pub fn data(&self) -> StarData {
        self.data
    }
    /// Static magnitude in the passband described by [`Self::photometry`].
    ///
    /// No light curve or passband conversion is applied. Per-entry measurement
    /// source, band, statistic and epoch are not authenticated for every row;
    /// see `docs/limitations.md` for known qualifications.
    pub fn magnitude(&self) -> f64 {
        self.magnitude
    }

    /// Magnitude and sparse per-entry provenance, without epoch interpolation.
    ///
    /// Unannotated entries explicitly retain an unqualified band and source;
    /// a catalogue's astrometric epoch does not establish its photometric epoch.
    /// An exact copy of a canonical row retains its annotation. Reusing its ID
    /// with different identity fields or measurement bits does not qualify the
    /// modified record, which keeps its own magnitude as an unqualified value.
    pub fn photometry(&self) -> StarPhotometry {
        let Some(value) = STAR_PHOTOMETRY
            .iter()
            .find_map(|&(id, value)| (id == self.id).then_some(value))
        else {
            return StarPhotometry::unqualified(self.magnitude);
        };
        let Some(canonical) = star_entry_by_id(self.id) else {
            return StarPhotometry::unqualified(self.magnitude);
        };
        // Provenance belongs to the complete immutable catalogue row, not an
        // identifier copied onto another value record. Compare IEEE words so
        // a different signed zero cannot inherit an annotation either.
        let same_identity = self.name == canonical.name
            && self.nomenclature == canonical.nomenclature
            && self.hip_number == canonical.hip_number;
        let words = |entry: &Self| {
            [
                entry.data.ra_j2000_deg,
                entry.data.dec_j2000_deg,
                entry.data.pm_ra_arcsec_yr,
                entry.data.pm_dec_arcsec_yr,
                entry.data.parallax_mas,
                entry.data.radial_velocity_km_s,
                entry.magnitude,
            ]
            .map(f64::to_bits)
        };
        if same_identity && words(self) == words(canonical) {
            value
        } else {
            StarPhotometry::unqualified(self.magnitude)
        }
    }
}

/// Phonetic replacements in the calculation order.
const PHONETIC_REPLACEMENTS: &[(&str, &str)] = &[
    ("LL", "L"),
    ("TT", "T"),
    ("SS", "S"),
    ("RR", "R"),
    ("PP", "P"),
    ("CC", "C"),
    ("NN", "N"),
    ("MM", "M"),
    ("FF", "F"),
    ("EU", "U"),
    ("EI", "I"),
    ("AU", "A"),
    ("AE", "E"),
    ("OE", "E"),
    ("OU", "O"),
    ("UE$", "U"),
    ("E$", ""),
    ("GH", "G"),
    ("PH", "F"),
    ("CK", "K"),
    ("CH", "K"),
    ("SCH", "SH"),
    ("X", "KS"),
    ("EE", "E"),
    ("AA", "A"),
    ("OO", "O"),
    ("II", "I"),
    ("UU", "U"),
];

/// Consonant set for the model inter-consonantal-`I` rule.
const PHONETIC_CONSONANTS: &[char] = &[
    'B', 'C', 'D', 'F', 'G', 'H', 'J', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X',
    'Y', 'Z',
];

/// Which resolver grammar classifies a search string.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StarFamily {
    /// `fixstar`/`fixstar_ut`/`fixstar_mag` grammar.
    V1,
    /// `fixstar2`/`fixstar2_ut`/`fixstar2_mag` grammar.
    V2,
}

/// A successful catalog resolution: the star id plus the canonical
/// `"Name,nomenclature"` form the entries return.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StarMatch {
    /// Internal star ID.
    pub star_id: i32,
    /// Canonical `"Name,nomenclature"` form.
    pub canonical: String,
}

/// A catalog lookup miss, carrying its explanation and classified search form.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StarFailure {
    /// Explanation of the catalog miss.
    pub message: String,
    /// Classified search form (`"empty"`, `"wildcard"`, `"nomenclature"`,
    /// `"sequential"`, `"name"`).
    pub search_type: &'static str,
}

impl StarFailure {
    /// Convert the query and search form into a typed native error.
    /// The explanation remains available in [`Self::message`].
    pub fn to_error(&self, star_name: &str) -> Error {
        Error::StarNotFound {
            star_id: star_name.to_string(),
            search_type: Some(self.search_type),
        }
    }
}

/// A catalog lookup miss or a structured validation failure.
#[derive(Clone, Debug, PartialEq)]
pub enum StarResolveError {
    /// No catalog entry matches the query.
    Miss(StarFailure),
    /// A native validation failure.
    Refused(Error),
}

// Catalog space motion and ICRS/J2000 frame geometry.
//
// Trigonometric parallax is milliarcseconds; nonpositive values use the
// one-gigaparsec convention (1e-6 mas). Distance is 1/sin(parallax radians) AU.
// The perspective factor 1/(1-radial_velocity/c) scales tangential and radial
// motion. Linear 3D propagation uses the J2000 position plus velocity times
// elapsed days. J2000 ecliptic output applies frame bias and fixed obliquity;
// of-date apparent-place corrections are composed by the serving layer.
// Epochs and vectors must be finite at checked public boundaries.

/// Parallax floor for unknown distances: a
/// non-positive catalog parallax is represented at one gigaparsec, i.e.
/// one nanoarcsecond (1e-6 milliarcseconds), by the parsec definition.
pub const STAR_PARALLAX_FLOOR_MAS: f64 = 1.0e-6;

/// Arcsecond-to-radian factor.
pub const STAR_ASEC2RAD: f64 = std::f64::consts::PI / 180.0 / 3600.0;

/// Speed of light in m/s.
pub const STAR_C_M_PER_S: f64 = 299792458.0;

/// Astronomical unit in km.
pub const STAR_AU_KM: f64 = 149597870.7;

/// Seconds per day.
pub const STAR_SECONDS_PER_DAY: f64 = 86400.0;

/// Cosine-of-declination floor below which right ascension is undefined
/// and kept unchanged.
pub const STAR_POLE_COS_DEC_FLOOR: f64 = 1e-10;

// Motion is a leaf: no observer selection, apparent reductions or output
// serving. Keep both public and crate-private entry paths on this facade.
#[path = "stars/motion.rs"]
mod motion;
pub(crate) use motion::star_icrs_initial_state_catalog_angles;
pub use motion::{
    propagate_proper_motion, proper_motion_ra_increment, star_geocentric_icrs,
    star_geometric_j2000_by_id, star_geometric_j2000_ecliptic, star_icrs_state,
    star_icrs_state_catalog_angles, star_initial_state_from_native, star_initial_state_from_stored,
};

/// Acquisition channel of one live catalog row.
/// Stored positions select the row before observer reads; the native vector
/// channel selects it after its original Time and observer have been acquired.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StarRowBackend {
    /// Scalar-TT LEB space-motion reduction.
    Stored,
    /// Original-Time vector reduction, also used by topocentric requests.
    Native,
}

/// Owned numerical inputs for catalog space motion.
///
/// Right ascension uses hours, declination degrees, proper motions mas/year,
/// parallax mas and radial velocity km/s. A catalog provider returns this whole
/// value once per sample; the native kernel performs all subsequent arithmetic.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NativeStarInputs {
    /// Right ascension, hours.
    pub ra_hours: f64,
    /// Declination, degrees.
    pub dec_degrees: f64,
    /// Proper motion in RA including cos(dec), mas/year.
    pub ra_mas_per_year: f64,
    /// Proper motion in declination, mas/year.
    pub dec_mas_per_year: f64,
    /// Measured parallax, mas; nonpositive values use the constructor floor.
    pub parallax_mas: f64,
    /// Radial velocity, km/s.
    pub radial_km_per_s: f64,
}

/// Owned numerical catalog row retained for one position sample.
///
/// Implementations return one complete NativeStarInputs value. They must not
/// hold resource or catalog locks while performing external acquisition. All
/// space-motion arithmetic remains in the native kernel.
pub trait AcquiredStarRow: Send + Sync {
    /// Acquire the complete native keyword tuple after observer acquisition.
    fn native_inputs(&self) -> Result<NativeStarInputs, Error>;
}

/// Native catalog selection with an independent owned row per sample.
///
/// The numerical star ID selects the row; selection and typed acquisition
/// errors propagate. A scalar velocity call checks membership before selecting
/// its position samples, while an already-resolved batch can omit that check.
pub trait StarRowAcquisition: Send + Sync {
    /// Consume only scalar velocity's membership check, without selecting a row.
    fn velocity_admission(&self, star_id: i32) -> Result<(), Error>;
    /// Select a fresh row at the backend's actual acquisition boundary.
    fn acquire_row(
        &self,
        star_id: i32,
        backend: StarRowBackend,
    ) -> Result<std::sync::Arc<dyn AcquiredStarRow>, Error>;
    /// Classify a typed acquisition error for an explicitly skipping batch.
    /// A classification error aborts the call; numerical failures are not silently
    /// converted to missing stars.
    fn is_skippable(&self, error: &Error) -> Result<bool, Error>;
    /// Decide whether a stored-data acquisition error can retry the native vector
    /// provider. The context distinguishes scalar and batch diagnostics. The
    /// standalone default retains the original error without a retry.
    fn retry_stored(&self, _error: &Error, _context: &'static str) -> Result<bool, Error> {
        Ok(false)
    }
}

// Apparent stellar positions from catalog motion and observer states.
//
// The pipeline forms the observer-relative ray, applies light travel time,
// then independently selected PPN deflection and Lorentz aberration, followed
// by the requested J2000 or true-of-date frame. Helio/barycentric requests
// suppress both apparent corrections at the entry layer. Topocentric states
// include the physical site and its observer-clock/frame composition.
// Inputs are finite TT epochs and AU/AU-day states; acquisition and frame
// errors propagate without replacing the selected model.

/// Barycentric ICRS observer state (AU, AU/day) for the star pipeline.
///
/// Selected by center: Earth for geocentric,
/// Sun for `FLG_HELCTR`, the SSB origin (zero position and velocity)
/// for `FLG_BARYCTR`. Private fields keep the aberration driver tied
/// to the deflection geometry: both vectors feed the same correction
/// chain.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StarObserver {
    pos_au: [f64; 3],
    vel_au_day: [f64; 3],
    site_gcrs_au: Option<[f64; 3]>,
}

impl StarObserver {
    /// Checked constructor: both vectors must be finite.
    pub fn new(pos_au: [f64; 3], vel_au_day: [f64; 3]) -> Result<Self, Error> {
        for v in pos_au.into_iter().chain(vel_au_day) {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            pos_au,
            vel_au_day,
            site_gcrs_au: None,
        })
    }

    /// Checked topocentric constructor: the BCRS state plus the
    /// geocentric GCRS site vector driving the Earth deflection term
    /// . All vectors must be finite.
    pub fn new_topo(
        pos_au: [f64; 3],
        vel_au_day: [f64; 3],
        site_gcrs_au: [f64; 3],
    ) -> Result<Self, Error> {
        for v in pos_au.into_iter().chain(vel_au_day).chain(site_gcrs_au) {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            pos_au,
            vel_au_day,
            site_gcrs_au: Some(site_gcrs_au),
        })
    }

    /// Barycentric ICRS observer position (AU).
    pub fn pos_au(&self) -> [f64; 3] {
        self.pos_au
    }

    /// Barycentric ICRS observer velocity (AU/day, aberration driver).
    pub fn vel_au_day(&self) -> [f64; 3] {
        self.vel_au_day
    }

    /// Geocentric GCRS site vector (AU) for a topocentric observer;
    /// `None` for Earth/Sun/origin observers.
    pub fn site_gcrs_au(&self) -> Option<[f64; 3]> {
        self.site_gcrs_au
    }
}

/// Of-date frame for the star pipeline: the biased ICRS-to-true-equator
/// matrix plus the true obliquity. J2000 requests ignore it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StarFrame {
    pn: [[f64; 3]; 3],
    eps_true_rad: f64,
}

impl StarFrame {
    /// Checked constructor: every entry and angle must be finite.
    pub fn new(pn: [[f64; 3]; 3], eps_true_rad: f64) -> Result<Self, Error> {
        for row in pn {
            for v in row {
                if !v.is_finite() {
                    return Err(Error::NonFiniteArgument { value: v });
                }
            }
        }
        if !eps_true_rad.is_finite() {
            return Err(Error::NonFiniteArgument {
                value: eps_true_rad,
            });
        }
        Ok(Self { pn, eps_true_rad })
    }

    /// Biased ICRS-to-true-equator-of-date matrix.
    pub fn pn(&self) -> [[f64; 3]; 3] {
        self.pn
    }

    /// True obliquity of date in radians.
    pub fn eps_true_rad(&self) -> f64 {
        self.eps_true_rad
    }
}

/// Remove `FLG_MOSEPH` and replace `FLG_SPEED3` with `FLG_SPEED`.
///
/// Other ephemeris bits are retained. This preprocessing does not validate the
/// flag combination or determine the returned ephemeris-selection bit.
pub fn preprocess_star_flags(iflag: i32) -> i32 {
    use crate::constants::{FLG_MOSEPH, FLG_SPEED, FLG_SPEED3};
    let mut flags = iflag & !FLG_MOSEPH;
    if flags & FLG_SPEED3 != 0 {
        flags = (flags & !FLG_SPEED3) | FLG_SPEED;
    }
    flags
}

/// Effective helio/barycentric center for a star request: `FLG_TOPOCTR` keeps the
/// Earth-based observer (returns 0 — the diurnal offset is not a
/// helio/bary re-centering), else `FLG_BARYCTR` wins over `FLG_HELCTR`,
/// else geocentric 0.
pub fn fixstar_observation_center(flags: i32) -> i32 {
    use crate::constants::{FLG_BARYCTR, FLG_HELCTR, FLG_TOPOCTR};
    if flags & FLG_TOPOCTR != 0 {
        return 0;
    }
    if flags & FLG_BARYCTR != 0 {
        return FLG_BARYCTR;
    }
    if flags & FLG_HELCTR != 0 {
        return FLG_HELCTR;
    }
    0
}

// Apparent samples consume explicit observers and motion states. Only the
// two sample-producing helpers need parent visibility for serving composition.
#[path = "stars/apparent.rs"]
mod apparent;
pub use apparent::{
    sealed_star_position, sealed_star_velocity, star_apparent_ecliptic_topo_via,
    star_apparent_ecliptic_via, star_velocity_via, topo_star_epoch_days,
};
pub(crate) use apparent::{
    star_apparent_ecliptic_initial_via, star_geometric_range_state, star_velocity_samples_via,
};
use apparent::{star_apparent_sample_via, star_apparent_topo_sample_via};

/// One reported stellar place and its physical geometric radial rate.
#[derive(Clone, Copy, Debug)]
pub(crate) struct StarSample {
    /// Reported longitude/latitude degrees and astrometric range AU.
    pub(crate) position: [f64; 3],
    /// Analytic derivative of the astrometric range, AU per TT day.
    pub(crate) radial_rate: f64,
}

/// Half-day central-difference step for geocentric star velocities
pub const STAR_GEOCENTRIC_VELOCITY_STEP_DAYS: f64 = 0.5;

// Anchor serving keeps live row/cache ownership and request-level error
// promotion together. Its internal geometry and sampling helpers stay private.
#[path = "stars/anchors.rs"]
mod anchors;

/// Live dynamic-ayanamsha anchor provider over geometric state sources.
///
/// Standalone rows come from [`crate::sidereal::ANCHOR_STARS`]; supplied live
/// acquisitions retain the actual defining STARS row and per-node cache. The borrowed
/// provider supplies barycentric Earth/deflector states in AU and AU/day at
/// TT epochs, whether backed by LEB or local SPK. Catalog propagation uses
/// the TDB-plus-observer-projection epoch and angle units;
/// apparent-place reductions use PPN deflection and relativistic aberration.
/// Mean-ecliptic longitudes use biased Vondrák 2011 precession and mean
/// obliquity, without nutation. Deflection-free longitudes interpolate three
/// 8-day samples; active deflection always evaluates at the requested date.
/// NOABERR removes the astrometric aberration displacement from the apparent
/// vector. Nutation and mode-specific target subtraction are downstream.
///
/// Formula modes Mardyks and Valens never route here. Missing Earth states,
/// interpolation-node coverage faults, frame faults and reduction errors
/// propagate; no defining-pair substitution or extrapolation is performed.
/// The basic constructor retains the stored provider's optional deflector
/// policy. [`Self::with_native_acquisition`] selects strict scalar-state reads;
/// [`Self::with_native_ray`] binds fresh source-Time native rays. Both preserve
/// mean-frame/node arithmetic and wrap consumed ranges at the outer request.
pub struct SealedAnchorProvider<'a, P: crate::positions::StateProvider> {
    provider: &'a P,
    native_state_at: Option<&'a AnchorStateAt<'a>>,
    native_ray_at: Option<&'a AnchorRayAt<'a>>,
    range_error_at: Option<&'a AnchorRangeErrorAt<'a>>,
    rows: Option<&'a dyn AnchorRowAcquisition>,
    acquired_ray_at: Option<&'a AcquiredAnchorRayAt<'a>>,
    dynamic_override: Option<&'a dyn crate::sidereal::DynamicAnchorProvider>,
}

/// Full scalar cache-key options for a live sidereal anchor sample.
/// The released mean-ecliptic routes use `nonut=true` and zero obliquity.
/// A geometric pole has `geometric=true` with the original false correction
/// flags; its numerical reduction still omits both physical corrections.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnchorSampleOptions {
    /// True obliquity, degrees; zero when the non-nutated path does not use it.
    pub eps_true: f64,
    /// Select the mean ecliptic instead of the true ecliptic.
    pub nonut: bool,
    /// Select the astrometric pole direction without physical corrections.
    pub geometric: bool,
    /// Omit observer-velocity aberration from a light-source anchor.
    pub noaberr: bool,
    /// Omit gravitational deflection from a light-source anchor.
    pub nogdefl: bool,
}

/// Result of the actual six-field, per-node cache lookup.
/// A hit owns a completed longitude; a miss owns the retained constructor
/// operands until its numerical work and success publication have finished.
pub enum AnchorSample {
    /// A previously completed longitude in degrees, with no new provider read.
    Cached(f64),
    /// A miss that must construct and calculate a new scalar sample.
    Pending(std::sync::Arc<dyn AcquiredAnchorSample>),
}

/// One cache miss's owned constructor and successful publication boundary.
/// Attached implementations preserve actual constructor/operator errors and
/// must not hold a cache or resource mutex across callbacks or Rust work.
pub trait AcquiredAnchorSample: Send + Sync {
    /// Evaluate the six Star keyword expressions after cache miss, in order.
    /// RA is hours, Dec degrees, proper motions mas/year, parallax mas and
    /// radial velocity km/s; the source `parallax` field is in arcseconds.
    fn native_inputs(&self) -> Result<NativeStarInputs, Error>;
    /// Publish only the complete, successfully reduced longitude in degrees.
    /// A refused publication remains an error; no partial entry is fabricated.
    fn publish(&self, longitude: f64) -> Result<(), Error>;
}

/// Actual `STARS[key]` row held by one anchor helper invocation.
/// This is distinct from FIXED_STARS and its integer body IDs. The same row
/// supplies new field reads at every exact-date/interpolation-node cache call.
pub trait AcquiredAnchorRow: Send + Sync {
    /// Read the original six cache fields and perform the actual scalar lookup.
    /// Errors in field access, hashing, lookup or retained factory propagate.
    fn sample(&self, tt: f64, options: AnchorSampleOptions) -> Result<AnchorSample, Error>;
    /// Acquire only the four Galactic Center RA constructor expressions.
    /// This constructor deliberately omits parallax and radial velocity:
    /// return their defaults without touching those fields.
    fn galactic_inputs(&self) -> Result<NativeStarInputs, Error>;
}

/// Live defining-namespace row lookup at the consumed stellar anchor stage.
/// Each helper reacquires its actual table/key; callbacks retain row owners,
/// not static catalog substitutions. A refusal never uses ANCHOR_STARS.
pub trait AnchorRowAcquisition: Send + Sync {
    /// Select the actual symbolic `STARS` row (SPICA, REVATI, GAL_CENTER, etc.).
    fn acquire_row(&self, key: &str) -> Result<std::sync::Arc<dyn AcquiredAnchorRow>, Error>;
    /// Classify a retained callback failure at the public fast-mean catch.
    /// `Some(true)` means the actual KeyError/ValueError was handled and its
    /// carrier retired; the caller retries the non-fast route with these same
    /// live acquisitions. `Some(false)` retains a terminal actual failure.
    /// `None` delegates only non-callback errors to the existing typed catch.
    /// Classification failures propagate instead of changing catch families.
    fn fast_mean_miss(&self, _error: &Error) -> Result<Option<bool>, Error> {
        Ok(None)
    }
}

/// Distinct main-image/Time acquisition order of the two anchor helpers.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AnchorRayRole {
    /// The uncached ecliptic helper acquires its ephemeris before scalar Time.
    EclipticLongitude,
    /// Galactic Center RA constructs scalar Time before ephemeris acquisition.
    GalacticRightAscension,
}

/// Source-unit constructor already reduced to ICRS position AU/velocity AU/day.
/// The callback receives the untouched scalar TT and actual helper role. It
/// owns source-Time/observer and apparent reduction, without reacquiring rows,
/// mean-frame projection, interpolation or an outer range wrapper.
pub type AcquiredAnchorRayAt<'a> =
    dyn Fn(([f64; 3], [f64; 3]), f64, bool, bool, AnchorRayRole) -> Result<[f64; 3], Error> + 'a;

/// Caller-owned native anchor acquisition at scalar TT Julian days.
/// The body is the Earth or a Sun/Jupiter/Saturn deflector library id.
/// Return barycentric ICRS position AU and velocity AU/day. Selected native
/// center/barycenter aliases and raw observation/closest-approach errors belong
/// to the caller. Every consumed error propagates, including a range miss;
/// no optional-deflector or substitute-provider policy is applied.
pub type AnchorStateAt<'a> = dyn Fn(i32, f64) -> Result<([f64; 3], [f64; 3]), Error> + 'a;

/// Caller-owned native catalog-anchor reduction at a fresh scalar TT Time.
/// `data` supplies ICRS J2000 catalog angles, proper motion, parallax and
/// radial velocity. Return the geocentric reduced ICRS ray in AU: geometric
/// when both toggles are true, aberration only under `nogdefl`, full apparent
/// minus the astrometric aberration displacement under `noaberr`, or full
/// apparent otherwise. Required Earth and deflector reads remain strict at
/// their actual observation/closest-approach Times; consumed errors propagate.
/// This synchronous borrowed callback supplies no mean-frame projection,
/// node interpolation, outer range promotion or escaped provider reference.
pub type AnchorRayAt<'a> = dyn Fn(&StarData, f64, bool, bool) -> Result<[f64; 3], Error> + 'a;

/// Wrap a consumed raw native range failure at the original anchor request TT.
/// The error retains its actual failed target, ephemeris seconds and bounds
/// when passed to this callback. The caller supplies the active main-image
/// public range envelope and absent star-body identity. Missing aliases and
/// all non-range errors bypass this callback unchanged. The callback is
/// synchronous and returns an owned error; no provider borrow may escape.
pub type AnchorRangeErrorAt<'a> = dyn Fn(f64, Error) -> Error + 'a;

#[cfg(test)]
#[path = "stars/live_anchor_row_tests.rs"]
mod live_anchor_row_tests;

#[cfg(test)]
#[path = "native_anchor_contract_tests.rs"]
mod native_anchor_contract_tests;

// Fixed-star public entry and output transformations.
//
// The apparent-place kernel supplies the selected observer-relative ray.
// This layer applies equatorial, sidereal, defining-epoch, XYZ and radian
// output flags and reports implied return bits. UT1 entries obtain TT from
// the configured clock. A dynamic ayanamsha requires its anchor provider;
// there is no substitution for unavailable data. Site position, motion,
// precession and nutation remain part of the native observer model.

/// Half-step (days) the flag layer reconstructs raw endpoints with
pub const STAR_FLAG_LAYER_STEP_DAYS: f64 = 0.5;

/// Central-difference half-step for topocentric star velocities
pub const STAR_TOPOCENTRIC_VELOCITY_STEP_DAYS: f64 = 0.005;

/// J2000 equatorial obliquity for the J2000 frame arm.
pub const STAR_J2000_OBLIQUITY_DEG: f64 = 84381.406 / 3600.0;

// Projection consumes the acquired raw place through an explicit environment.
// Keep the environment contract public here; only its borrowed serving adapter
// needs parent visibility for composition at the existing call sites.
#[path = "stars/flags.rs"]
mod flags;
use flags::SealedFlagEnv;
pub use flags::{
    apply_fixstar_flags_via, fixstar_positional_transform_via, fixstar_ret_flags,
    precess_ecliptic_to_j2000,
};

/// External inputs the flag layer needs beyond the raw six-slot place
pub trait FixstarFlagEnvironment {
    /// Backward precession to the J2000 ecliptic.
    fn precess_ecliptic_to_j2000(
        &self,
        lon_deg: f64,
        lat_deg: f64,
        jd_tt: f64,
    ) -> Result<(f64, f64), Error>;
    /// Nutation in longitude in radians.
    fn nutation_dpsi_rad(&self, jd_tt: f64) -> Result<f64, Error>;
    /// Mean obliquity of date in degrees.
    fn mean_obliquity_deg(&self, jd_tt: f64) -> Result<f64, Error>;
    /// True obliquity of date in degrees.
    fn true_obliquity_deg(&self, jd_tt: f64) -> Result<f64, Error>;
    /// UT1 epoch for a TT epoch.
    fn ut1_for_tt(&self, jd_tt: f64) -> Result<f64, Error>;
    /// Mean ayanamsha in degrees at a UT epoch under explicit ex flags
    fn mean_ayanamsha_deg_ut(&self, jd_ut: f64, ex_flags: i32) -> Result<f64, Error>;
}

/// Explicit WGS 84 star observer site: east-positive longitude and latitude
/// in degrees, height in metres. The serving borrows an owned scalar value;
/// this type stores no context, provider or mutable global configuration.
/// Construction admits finite components only. Geographic range policy belongs
/// to the caller's observer/configuration layer, not to this value type.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StarTopoSite {
    lon_deg: f64,
    lat_deg: f64,
    height_m: f64,
}

impl StarTopoSite {
    /// Store finite longitude/latitude degrees and height metres.
    /// Returns NonFiniteArgument for the first nonfinite component in that
    /// order. No geographic range, provider or clock is acquired here;
    /// consumed angles retain their trigonometric meaning in the observer.
    pub fn new(lon_deg: f64, lat_deg: f64, height_m: f64) -> Result<Self, Error> {
        for v in [lon_deg, lat_deg, height_m] {
            if !v.is_finite() {
                return Err(Error::NonFiniteArgument { value: v });
            }
        }
        Ok(Self {
            lon_deg,
            lat_deg,
            height_m,
        })
    }

    /// Geographic longitude in degrees (East positive).
    pub fn lon_deg(&self) -> f64 {
        self.lon_deg
    }

    /// Geographic latitude in degrees (North positive).
    pub fn lat_deg(&self) -> f64 {
        self.lat_deg
    }

    /// Altitude above sea level in meters.
    pub fn height_m(&self) -> f64 {
        self.height_m
    }
}

/// One owned star result: six-slot place, canonical name and returned flags.
///
/// Spherical slots are `[longitude, latitude, distance, longitude_rate,
/// latitude_rate, radial_rate]`, in degrees, AU and corresponding units per
/// day. `FLG_EQUATORIAL` replaces longitude/latitude with right ascension and
/// declination, still in degrees. `FLG_XYZ` selects Cartesian AU/AU-day slots
/// and takes precedence over `FLG_RADIANS`, which otherwise changes angular
/// slots to radians and radians/day. Without SPEED the rate slots are zero.
/// The producing entry defines the observation clock, frame and flag echo;
/// in particular TT entries echo input flags, whereas UT1 entries normalize them.
pub type WrappedStarEntry = ([f64; 6], String, i32);

/// Injected native catalog-star place at a fresh scalar TT epoch.
///
/// Arguments are catalog id, TT JD, NOABERR, NOGDEFL, J2000-frame selection,
/// effective observing center, optional WGS 84 site and speed selection.
/// Return longitude/latitude in degrees, distance in AU and corresponding
/// daily rates. The callback owns source vector acquisition and the reported
/// position stencil; errors propagate before the common flag layer runs.
/// It must preserve the independent bundled frame clock and must not replace
/// native SPK geometry with LEB or caller-model TT.
pub type StarRawPlaceAt<'a> = dyn Fn(i32, f64, bool, bool, bool, i32, Option<StarTopoSite>, bool) -> Result<[f64; 6], Error>
    + 'a;

/// Native public UT batch callback over one shared Time.
///
/// Inputs are [`StarFamily::V1`] names in order, original UT1 JD, raw flags, optional site and
/// skip-errors selection. The callback resolves all names before shared
/// geocentric acquisition, retains that Time's split TT/TDB and returns owned
/// six-slot/name/flag entries or absent slots. Topocentric, fixed-epoch and
/// active SIDBIT requests delegate to the wrapped scalar pipeline instead.
/// Shared observer failures remain terminal even when skip-errors is enabled.
pub type StarBatchAt<'a> = dyn Fn(
        &[&str],
        f64,
        i32,
        Option<StarTopoSite>,
        bool,
    ) -> Result<Vec<Option<WrappedStarEntry>>, Error>
    + 'a;

/// Shared batch callback for already acquired catalog identities.
///
/// Arguments are resolved matches, parallel per-star row suppliers, UT1 Julian
/// day, raw flags, optional site and skip-errors selection. The first two slices
/// must have the same length, with a row supplier for each resolved match;
/// empty match slots represent misses already owned by the caller's resolver.
/// Results retain input order and use [`WrappedStarEntry`] units and layout.
/// Shared acquisition failures remain errors; only failures classified as
/// skippable by the bound implementation become absent result slots.
/// No name lookup or catalog row is frozen by this numerical batch interface.
pub type StarResolvedBatchAt<'a> = dyn Fn(
        &[Option<StarMatch>],
        &[Option<std::sync::Arc<dyn StarRowAcquisition>>],
        f64,
        i32,
        Option<StarTopoSite>,
        bool,
    ) -> Result<Vec<Option<WrappedStarEntry>>, Error>
    + 'a;

/// Shared batch callback with a separately borrowed dynamic anchor provider.
/// Arguments and units follow [`StarResolvedBatchAt`]. The final optional
/// provider changes only sidereal shaping; shared observer/Time, fixed-star
/// rows and catches retain their existing owners. It is synchronous and must
/// not retain the supplied provider after return. None keeps static behavior.
pub type StarResolvedAnchorBatchAt<'a> = dyn Fn(
        &[Option<StarMatch>],
        &[Option<std::sync::Arc<dyn StarRowAcquisition>>],
        f64,
        i32,
        Option<StarTopoSite>,
        bool,
        Option<&dyn crate::sidereal::DynamicAnchorProvider>,
    ) -> Result<Vec<Option<WrappedStarEntry>>, Error>
    + 'a;

/// Fixed-star apparent-place serving over a borrowed sealed state provider.
///
/// Composes catalog propagation, observer selection, aberration/deflection,
/// mean/true/fixed-epoch frames, coordinate/output flags and batch entries.
/// Earth, Sun, SSB and explicitly supplied topocentric sites are distinct
/// observing centers; helio/barycentric requests suppress aberration and
/// deflection. Results contain longitude/latitude in degrees, distance in AU
/// and per-day rates unless Cartesian or radian flags select another layout.
///
/// [`Self::new`] retains caller-supplied model Delta-T conversion callbacks.
/// Select [`Self::with_default_timescale`] to acquire star TT and observer UT1
/// through the independent bundled enhanced timescale instead, including
/// evaluation TT for SIDBIT zero points. Model Delta-T remains separate for
/// user/Valens anchors and ordinary sidereal-model terms. Custom timescales/live
/// IERS are not supplied by that fixed-data selector. Sealed coverage misses propagate,
/// without a network fallback or an invented catalog/anchor position.
///
/// Resolve names before calling through [`resolve_star_ref`] or
/// [`resolve_star2`] and [`StarMatch`]. Magnitudes belong to the separate
/// magnitude layer. Model and finite comparison agreement do not certify
/// physical stellar-position accuracy.
pub struct SealedStarServing<'a, P: crate::positions::StateProvider> {
    provider: &'a P,
    nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64), Error>,
    deltat_days_at: &'a dyn Fn(f64) -> Result<f64, Error>,
    use_default_timescale: bool,
    raw_place_at: Option<&'a StarRawPlaceAt<'a>>,
    batch_at: Option<&'a StarBatchAt<'a>>,
    resolved_batch_at: Option<&'a StarResolvedBatchAt<'a>>,
    resolved_anchor_batch_at: Option<&'a StarResolvedAnchorBatchAt<'a>>,
    rows: Option<&'a dyn StarRowAcquisition>,
    stored_vector_error_at: Option<&'a dyn Fn(Error) -> Error>,
    anchor_at: Option<(&'a AnchorRayAt<'a>, &'a AnchorRangeErrorAt<'a>)>,
    anchor_rows: Option<&'a dyn AnchorRowAcquisition>,
    acquired_anchor_at: Option<(&'a AcquiredAnchorRayAt<'a>, &'a AnchorRangeErrorAt<'a>)>,
    dynamic_anchor_provider: Option<&'a dyn crate::sidereal::DynamicAnchorProvider>,
    sidereal: crate::sidereal::SiderealSelection,
}

impl<'a, P: crate::positions::StateProvider> SealedStarServing<'a, P> {
    /// Provenance of the provider actually bound to this calculation.
    ///
    /// The label describes retained numerical resources, not the requested
    /// mode. Automatic mode can therefore report either LEB or native JPL.
    /// Unknown application-provided providers retain their own label.
    pub fn source_tag(&self) -> &'static str {
        match self.provider.provider_name() {
            "sealed-leb/single" | "sealed-leb/composite" | "sealed-leb/tiered" => "LEB",
            "native-jpl/local-spk" => "NativeJpl",
            name => name,
        }
    }

    /// Borrow a state provider, nutation callback, model Delta-T callback and
    /// sidereal selection for the lifetime of this serving.
    ///
    /// Nutation returns radians at TT; Delta-T returns days at the requested
    /// epoch. By default that callback supplies the UT1→TT addition
    /// and scalar TT→UT1 subtraction, and model-based sidereal terms.
    /// [`Self::with_default_timescale`] changes time acquisition, including
    /// SIDBIT zero-point evaluation TT, retaining the model callback for
    /// user/Valens anchors and ordinary sidereal terms. Dynamic star/galactic
    /// anchors evaluate through [`Self::anchor_provider`] over the bound states
    /// or separately selected native ray callback. Consumed errors propagate
    /// without a catalog, state or formula substitute.
    pub fn new(
        provider: &'a P,
        nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64), Error>,
        deltat_days_at: &'a dyn Fn(f64) -> Result<f64, Error>,
        sidereal: crate::sidereal::SiderealSelection,
    ) -> Self {
        Self {
            provider,
            nutation_at,
            deltat_days_at,
            use_default_timescale: false,
            raw_place_at: None,
            batch_at: None,
            resolved_batch_at: None,
            resolved_anchor_batch_at: None,
            rows: None,
            stored_vector_error_at: None,
            anchor_at: None,
            anchor_rows: None,
            acquired_anchor_at: None,
            dynamic_anchor_provider: None,
            sidereal,
        }
    }

    /// Select the fixed bundled timescale for star epochs and observer UT1.
    ///
    /// UT entries acquire TT independently of the bound model Delta-T;
    /// TT observer epochs use the matching fractional-day inverse. SIDBIT
    /// projection zero points acquire their evaluation TT through this clock
    /// too. The model callback remains available for user-defined UT anchors,
    /// the Valens defining epoch and ordinary sidereal-model terms. The original
    /// [`Self::new`] callback-based conversion is retained unless this builder
    /// is selected.
    /// Custom time tables/live IERS are not represented by this fixed-data route.
    pub fn with_default_timescale(mut self) -> Self {
        self.use_default_timescale = true;
        self
    }

    /// Bind native raw acquisition while retaining catalog/flag/wrapper rules.
    ///
    /// Both callbacks are borrowed for this serving's lifetime. `raw_place_at`
    /// owns fresh scalar-TT star acquisition and reported-position rates;
    /// `batch_at` owns the original-UT shared-vector batch. Missing callbacks
    /// retain the existing sealed pipeline without changing its arithmetic.
    /// No callback or provider reference is returned through an entry result.
    pub fn with_native_places(
        mut self,
        raw_place_at: &'a StarRawPlaceAt<'a>,
        batch_at: Option<&'a StarBatchAt<'a>>,
    ) -> Self {
        self.raw_place_at = Some(raw_place_at);
        self.batch_at = batch_at;
        self
    }

    /// Borrow a live catalog callback without selecting a row.
    /// Each consumed sample owns its separate lease. Leaving this unbound keeps
    /// the original immutable catalog path and its accepted arithmetic.
    pub fn with_live_rows(mut self, rows: &'a dyn StarRowAcquisition) -> Self {
        self.rows = Some(rows);
        self
    }

    /// Retain the selected reader's vector-range presentation constructor.
    /// Only live stored-backed vector helpers consume this mapper; native SPK
    /// and static error channels retain their own acquisition paths.
    pub fn with_stored_vector_error(mut self, mapper: &'a dyn Fn(Error) -> Error) -> Self {
        self.stored_vector_error_at = Some(mapper);
        self
    }

    /// Bind a shared native batch that consumes already resolved identities.
    /// Scalar delegated batches still use the ordinary wrapped entry, while
    /// this callback retains original UT1 Time and shared observer acquisition.
    pub fn with_resolved_native_batch(mut self, batch: &'a StarResolvedBatchAt<'a>) -> Self {
        self.resolved_batch_at = Some(batch);
        self
    }

    /// Bind a resolved batch that forwards its current independent anchors.
    /// Unlike a captured static delegate, this callback receives the actual
    /// borrowed override at each batch invocation. The old callback remains a
    /// compatible fallback when this additive hook is absent.
    pub fn with_resolved_native_anchor_batch(
        mut self,
        batch: &'a StarResolvedAnchorBatchAt<'a>,
    ) -> Self {
        self.resolved_anchor_batch_at = Some(batch);
        self
    }

    /// Bind strict native dynamic-anchor acquisition and range promotion.
    /// These callbacks borrow the same immutable provider as catalog places.
    /// Anchor evaluation acquires a fresh scalar TT Time at the exact date or
    /// each interpolation node; the enclosing anchor provider promotes only
    /// raw SPK range failures at its original request TT. Defining model epochs
    /// and the selected bundled evaluation clock remain separate. The binding
    /// lasts only for this serving's synchronous borrowed lifetime.
    pub fn with_native_anchors(
        mut self,
        ray_at: &'a AnchorRayAt<'a>,
        range_error_at: &'a AnchorRangeErrorAt<'a>,
    ) -> Self {
        self.anchor_at = Some((ray_at, range_error_at));
        self
    }

    /// Bind actual defining sidereal STARS rows, independently of FIXED_STARS.
    /// No lookup is performed until a dynamic anchor is consumed. Each helper
    /// owns its row; every node retains the actual cache/constructor boundary.
    pub fn with_live_anchor_rows(mut self, rows: &'a dyn AnchorRowAcquisition) -> Self {
        self.anchor_rows = Some(rows);
        self
    }

    /// Borrow the same serving callbacks with a call-local anchor supplier.
    ///
    /// Only the sidereal catalog acquisition changes. Provider, fixed-star
    /// rows, observer/Time callbacks, range mappers, model clock and selection
    /// retain their original owners and stages. No callback or row is consumed
    /// by constructing this view, and it cannot outlive either borrowed owner.
    pub fn borrowing_live_anchor_rows<'b>(
        &'b self,
        rows: &'b dyn AnchorRowAcquisition,
    ) -> SealedStarServing<'b, P> {
        let mut serving = self.borrowing_view();
        serving.anchor_rows = Some(rows);
        serving
    }

    /// Borrow independent dynamic anchors without changing the fixed-star feed.
    ///
    /// The fixed-star position and observer retain their selected provider and
    /// Time. Only nested sidereal anchor requests use `provider`; this permits
    /// an Auto stored body to consume the separately selected Native main.
    /// Cache and callback failures remain with the supplied owner. Construction
    /// performs no row, provider or clock operation and cannot outlive either
    /// borrowed owner. Other flag/frame operations retain their old callbacks.
    pub fn borrowing_dynamic_anchors<'b>(
        &'b self,
        provider: &'b dyn crate::sidereal::DynamicAnchorProvider,
    ) -> SealedStarServing<'b, P> {
        let mut serving = self.borrowing_view();
        serving.dynamic_anchor_provider = Some(provider);
        serving
    }

    /// Shorten the lifetime of the existing immutable callback view.
    fn borrowing_view(&self) -> SealedStarServing<'_, P> {
        SealedStarServing {
            provider: self.provider,
            nutation_at: self.nutation_at,
            deltat_days_at: self.deltat_days_at,
            use_default_timescale: self.use_default_timescale,
            raw_place_at: self.raw_place_at,
            batch_at: self.batch_at,
            resolved_batch_at: self.resolved_batch_at,
            resolved_anchor_batch_at: self.resolved_anchor_batch_at,
            rows: self.rows,
            stored_vector_error_at: self.stored_vector_error_at,
            anchor_at: self.anchor_at,
            anchor_rows: self.anchor_rows,
            acquired_anchor_at: self.acquired_anchor_at,
            dynamic_anchor_provider: self.dynamic_anchor_provider,
            sidereal: self.sidereal,
        }
    }

    /// Bind source-Time rays for a live anchor's already constructed vectors.
    /// The actual helper role controls ephemeris/Time order. Outer projection,
    /// cache publication and range promotion remain with the anchor provider.
    pub fn with_acquired_anchors(
        mut self,
        ray_at: &'a AcquiredAnchorRayAt<'a>,
        range_error_at: &'a AnchorRangeErrorAt<'a>,
    ) -> Self {
        self.acquired_anchor_at = Some((ray_at, range_error_at));
        self
    }

    /// Shape an already acquired native batch place through the common flags.
    ///
    /// `raw` is longitude/latitude in degrees, distance AU and daily rates at
    /// scalar `jd_tt`; `flags` have undergone star preprocessing. `j2000_native`
    /// identifies the raw fixed frame. Nutation, coordinate conversions,
    /// sidereal model/live anchors and XYZ/radian shaping retain their ordinary
    /// error/clock rules. This method acquires no observer or replacement state.
    pub fn shape_native_place(
        &self,
        raw: [f64; 6],
        jd_tt: f64,
        flags: i32,
        j2000_native: bool,
    ) -> Result<[f64; 6], Error> {
        apply_fixstar_flags_via(
            raw,
            jd_tt,
            flags,
            j2000_native,
            &SealedFlagEnv { serving: self },
        )
    }

    /// Dynamic star/galactic anchors on the mean Vondrák ecliptic at TT.
    /// A native callback supplies fresh source-Time Earth/catalog geometry and
    /// strict apparent deflectors. The common provider preserves exact-date
    /// evaluation with active deflection, quadratic eight-day nodes without it,
    /// galactic hour-circle projection and original-request range promotion.
    /// LEB retains its bound stored geometric states and optional deflector
    /// policy. Missing aliases and consumed frame/time/provider errors propagate;
    /// no catalog, clock or numerical fallback is supplied by this method.
    pub fn anchor_provider(&self) -> SealedAnchorProvider<'_, P> {
        let mut anchors = SealedAnchorProvider::new(self.provider);
        if let Some((ray_at, range_error_at)) = self.anchor_at {
            anchors = anchors.with_native_ray(ray_at, range_error_at);
        }
        if let Some((ray_at, range_error_at)) = self.acquired_anchor_at {
            anchors = anchors.with_acquired_ray(ray_at, range_error_at);
        }
        if let Some(rows) = self.anchor_rows {
            anchors = anchors.with_live_rows(rows);
        }
        if let Some(provider) = self.dynamic_anchor_provider {
            anchors = anchors.with_dynamic_override(provider);
        }
        anchors
    }

    /// TT Julian Day corresponding to a UT1 Julian Day.
    ///
    /// Uses either the explicitly selected bundled timescale or the bound
    /// callback's days offset added to UT1. The bundled route is independent
    /// of model/tidal/user Delta-T settings; errors propagate.
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64, Error> {
        if self.use_default_timescale {
            return crate::time::default_timescale_tt_from_ut1(jd_ut);
        }
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_days_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// Observer UT1 Julian Day corresponding to TT Julian Day.
    ///
    /// The bundled route subtracts its TT-evaluated offset from the fractional
    /// day before restoring the whole day. Legacy conversion subtracts the
    /// bound callback's offset from scalar TT. Neither guarantees exact
    /// round-trip identity; non-finite inputs/offsets are refused.
    pub fn ut1_for_tt(&self, jd_tt: f64) -> Result<f64, Error> {
        if self.use_default_timescale {
            return crate::time::default_timescale_ut1_from_tt(jd_tt);
        }
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let dt = (self.deltat_days_at)(jd_tt)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_tt - dt)
    }

    /// Stored-provider observer at one scalar TT epoch, in ICRS AU and AU/day.
    /// HELCTR reads Sun; BARYCTR is the SSB origin; TOPOCTR requires an explicit
    /// site before reading Earth, then adds its bundled-clock WGS 84 offset.
    /// Otherwise read Earth. Site absence, state, bundled-time and frame failures
    /// propagate. The separately bound native raw callback owns its retained-Time
    /// observer route; this helper does not replace it.
    fn observer_at(
        &self,
        center_flag: i32,
        topo: Option<StarTopoSite>,
        jd_tt: f64,
    ) -> Result<StarObserver, Error> {
        use crate::constants::{EARTH, FLG_BARYCTR, FLG_HELCTR, FLG_TOPOCTR, SUN};
        if center_flag == FLG_HELCTR {
            let state = self.provider.geometric_state(SUN, jd_tt)?;
            return StarObserver::new(state.pos(), state.vel());
        }
        if center_flag == FLG_BARYCTR {
            return StarObserver::new([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        }
        if center_flag == FLG_TOPOCTR {
            // The site resolves before any ephemeris work: a missing site refuses, never silent geocentric.
            let site = topo.ok_or(Error::ConfigurationMissing {
                missing_config: "observer_location",
                suggestion: "topocentric star places need an explicit site",
            })?;
            // equinox-based observer: LEB BCRS Earth state plus the
            // `itrs` topocentric offset, not the LEB-direct `c2t` composition
            // the planet arm serves. Frame time comes from the bundled
            // timescale, like the reference bundled timescale, and not from
            // a caller-selected Delta-T model.
            let earth = self.provider.geometric_state(EARTH, jd_tt)?;
            let delta_t_seconds = crate::smh2016::bundled_deltat_seconds(jd_tt)?;
            let jd_tdb = crate::time::tt_to_tdb_jd(jd_tt)?;
            let (pos, vel, gcrs) = crate::positions::reference_frames::barycentric_site_state(
                earth.pos(),
                earth.vel(),
                site.lon_deg(),
                site.lat_deg(),
                site.height_m(),
                jd_tt,
                delta_t_seconds,
                jd_tdb,
            )?;
            return StarObserver::new_topo(pos, vel, gcrs);
        }
        let earth = self.provider.geometric_state(EARTH, jd_tt)?;
        StarObserver::new(earth.pos(), earth.vel())
    }

    /// Acquire a six-slot place at scalar TT before public output shaping.
    /// A bound native callback owns its raw vector channel. Otherwise use the
    /// stored provider, selected observer and apparent/frame reduction. Output is
    /// longitude/latitude degrees, AU distance and daily rates. SPEED evaluates
    /// central, previous and next reported positions with h=0.005 days for Topos
    /// or h=0.5 otherwise; absent SPEED returns positive-zero rate values.
    /// Unknown catalog id, clock, observer, correction and frame errors propagate
    /// in consumption order; no other provider or stencil is substituted.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn raw_place(
        &self,
        star_id: i32,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
        j2000_frame: bool,
        center_flag: i32,
        topo: Option<StarTopoSite>,
        want_speed: bool,
    ) -> Result<[f64; 6], Error> {
        self.raw_place_with_admission(
            star_id,
            jd_tt,
            noaberr,
            nogdefl,
            j2000_frame,
            center_flag,
            topo,
            want_speed,
            true,
        )
    }

    /// Raw scalar samples with the enclosing caller's velocity-admission role.
    /// The shared stored batch directly evaluates three positions and therefore
    /// does not consume the scalar velocity helper's extra membership check.
    #[allow(clippy::too_many_arguments)]
    fn raw_place_with_admission(
        &self,
        star_id: i32,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
        j2000_frame: bool,
        center_flag: i32,
        topo: Option<StarTopoSite>,
        want_speed: bool,
        velocity_admission: bool,
    ) -> Result<[f64; 6], Error> {
        if let Some(raw_at) = self.raw_place_at {
            return raw_at(
                star_id,
                jd_tt,
                noaberr,
                nogdefl,
                j2000_frame,
                center_flag,
                topo,
                want_speed,
            );
        }
        use crate::constants::FLG_TOPOCTR;
        if let Some(rows) = self.rows {
            if want_speed && velocity_admission {
                rows.velocity_admission(star_id)?;
            }
            let position_at = |jd| {
                self.live_stored_position(
                    rows,
                    star_id,
                    jd,
                    noaberr,
                    nogdefl,
                    j2000_frame,
                    center_flag,
                    topo,
                )
            };
            if want_speed {
                let h = if center_flag == FLG_TOPOCTR {
                    STAR_TOPOCENTRIC_VELOCITY_STEP_DAYS
                } else {
                    STAR_GEOCENTRIC_VELOCITY_STEP_DAYS
                };
                return star_velocity_samples_via(position_at, jd_tt, h);
            }
            let p = position_at(jd_tt)?.position;
            return Ok([p[0], p[1], p[2], 0.0, 0.0, 0.0]);
        }
        let data = fixed_star_data(star_id).ok_or_else(|| Error::StarNotFound {
            star_id: star_id.to_string(),
            search_type: Some("sequential"),
        })?;
        let h = if center_flag == FLG_TOPOCTR {
            STAR_TOPOCENTRIC_VELOCITY_STEP_DAYS
        } else {
            STAR_GEOCENTRIC_VELOCITY_STEP_DAYS
        };
        let position_at = |jd: f64| -> Result<StarSample, Error> {
            // UT1 acquisition errors propagate; the topocentric observer below
            // takes its frame time from the bundled timescale instead.
            let _jd_ut = self.ut1_for_tt(jd)?;
            let observer = self.observer_at(center_flag, topo, jd)?;
            // Topocentric requests follow the topocentric apparent-place
            // composition (star epoch, light time and Lorentz aberration,
            // Earth deflection, equinox-based frames); every other center
            // keeps the LEB-direct core the model serves there.
            if center_flag == FLG_TOPOCTR {
                star_apparent_topo_sample_via(
                    &data,
                    jd,
                    &observer,
                    noaberr,
                    nogdefl,
                    j2000_frame,
                    |id, jd| Ok(self.provider.geometric_state(id, jd)?.pos()),
                )
            } else {
                let frame = if j2000_frame {
                    None
                } else {
                    let (dpsi, deps) = (self.nutation_at)(jd)?;
                    let frames = crate::positions::frame_matrices_for_date(jd, dpsi, deps)?;
                    Some(StarFrame::new(frames.pn(), frames.eps_true_rad())?)
                };
                star_apparent_sample_via(
                    &data,
                    jd,
                    &observer,
                    noaberr,
                    nogdefl,
                    j2000_frame,
                    |id, jd| Ok(self.provider.geometric_state(id, jd)?.pos()),
                    frame,
                )
            }
        };
        if want_speed {
            star_velocity_samples_via(position_at, jd_tt, h)
        } else {
            let triple = position_at(jd_tt)?.position;
            Ok([triple[0], triple[1], triple[2], 0.0, 0.0, 0.0])
        }
    }

    /// Live stored sample with its row selected before the observer read.
    /// Topocentric samples instead use the native constructor after their
    /// independently acquired Time/observer. Both tails retain the row through
    /// corrections/frame errors; no lease survives into the next sample.
    #[allow(clippy::too_many_arguments)]
    fn live_stored_position(
        &self,
        rows: &dyn StarRowAcquisition,
        id: i32,
        jd: f64,
        noaberr: bool,
        nogdefl: bool,
        j2000: bool,
        center: i32,
        topo: Option<StarTopoSite>,
    ) -> Result<StarSample, Error> {
        match self.live_stored_position_attempt(rows, id, jd, noaberr, nogdefl, j2000, center, topo)
        {
            Err(error) if center != crate::constants::FLG_TOPOCTR => {
                if rows.retry_stored(&error, "star")? {
                    crate::positions::star_serving::live_stored_vector_place_via(
                        self.provider,
                        self.stored_vector_error_at,
                        rows,
                        id,
                        jd,
                        noaberr,
                        nogdefl,
                        j2000,
                        center,
                    )
                } else {
                    Err(error)
                }
            }
            result => result,
        }
    }

    /// One stored attempt, with no catch crossing the subsequent vector retry.
    #[allow(clippy::too_many_arguments)]
    fn live_stored_position_attempt(
        &self,
        rows: &dyn StarRowAcquisition,
        id: i32,
        jd: f64,
        noaberr: bool,
        nogdefl: bool,
        j2000: bool,
        center: i32,
        topo: Option<StarTopoSite>,
    ) -> Result<StarSample, Error> {
        use crate::constants::FLG_TOPOCTR;
        let row = if center != FLG_TOPOCTR {
            Some(rows.acquire_row(id, StarRowBackend::Stored)?)
        } else {
            None
        };
        let observer = self.observer_at(center, topo, jd)?;
        if center == FLG_TOPOCTR {
            return crate::positions::star_serving::live_topocentric_place_via(
                self.provider,
                self.stored_vector_error_at,
                rows,
                id,
                jd,
                &observer,
                noaberr,
                nogdefl,
                j2000,
            );
        }
        let row = row.expect("stored row selected before observer");
        let initial = star_initial_state_from_stored(row.as_ref())?;
        star_apparent_ecliptic_initial_via(
            initial,
            jd,
            &observer,
            noaberr,
            nogdefl,
            j2000,
            |body, epoch| Ok(self.provider.geometric_state(body, epoch)?.pos()),
            || {
                if j2000 {
                    return Ok(None);
                }
                let (dpsi, deps) = (self.nutation_at)(jd)?;
                let frames = crate::positions::frame_matrices_for_date(jd, dpsi, deps)?;
                Ok(Some(StarFrame::new(frames.pn(), frames.eps_true_rad())?))
            },
        )
    }

    /// Reduce an already-resolved star at scalar TT to place and returned flags.
    /// Preprocess calculation flags, force centered/TRUEPOS correction selection,
    /// acquire the raw place, then apply coordinate/sidereal/output shaping.
    /// `implied_echo` chooses UT normalization versus TT echo. Positions
    /// and rates have the units selected by the output flags; consumed catalog,
    /// observer, provider, model and frame errors propagate.
    fn entry_core(
        &self,
        star_id: i32,
        jd_tt: f64,
        flags_in: i32,
        implied_echo: bool,
        topo: Option<StarTopoSite>,
    ) -> Result<([f64; 6], i32), Error> {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NOABERR, FLG_NOGDEFL, FLG_SPEED, FLG_TOPOCTR,
            FLG_TRUEPOS,
        };
        let ret_flags = fixstar_ret_flags(flags_in, implied_echo);
        let flags = preprocess_star_flags(flags_in);
        let mut noaberr = flags & FLG_NOABERR != 0 || flags & FLG_TRUEPOS != 0;
        let mut nogdefl = flags & FLG_NOGDEFL != 0 || flags & FLG_TRUEPOS != 0;
        let center = fixstar_observation_center(flags);
        if center != 0 {
            noaberr = true;
            nogdefl = true;
        }
        let use_j2000 = flags & FLG_J2000 != 0;
        let eff_center = if flags & FLG_TOPOCTR != 0 {
            FLG_TOPOCTR
        } else if center == FLG_BARYCTR {
            FLG_BARYCTR
        } else if center == FLG_HELCTR {
            FLG_HELCTR
        } else {
            0
        };
        let raw = self.raw_place(
            star_id,
            jd_tt,
            noaberr,
            nogdefl,
            use_j2000,
            eff_center,
            topo,
            flags & FLG_SPEED != 0,
        )?;
        let env = SealedFlagEnv { serving: self };
        let place = apply_fixstar_flags_via(raw, jd_tt, flags, use_j2000, &env)?;
        Ok((place, ret_flags))
    }

    /// Apparent place of an already-resolved star at a UT1 Julian Day.
    ///
    /// Acquire TT through the selected timescale, then apply observer/frame
    /// and output flags. Return the owned [`WrappedStarEntry`] layout, with
    /// normalized flags and applicable implied bits. The default spherical
    /// frame is the true ecliptic/equinox of date; flags select equatorial,
    /// J2000, mean or sidereal coordinates and Cartesian/radian representation.
    /// Use [`Self::fixstar_ut_wrapped`] for fixed-epoch or active sidereal-plane
    /// dispatch. Provider, clock, configuration and coordinate errors propagate;
    /// a topocentric request requires an explicit `topo` site.
    pub fn fixstar_ut(
        &self,
        star: &StarMatch,
        tjdut: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
    ) -> Result<([f64; 6], String, i32), Error> {
        let jd_tt = self.jd_tt_for_ut(tjdut)?;
        let (place, ret) = self.entry_core(star.star_id, jd_tt, flags, true, topo)?;
        Ok((place, star.canonical.clone(), ret))
    }

    /// Place of an already-resolved star at finite scalar TT Julian day.
    /// Use the common entry core without fixed-epoch/projection prologues and echo
    /// the caller flags. Output owns six components, canonical name and
    /// retflags; spherical units are degrees/AU and their daily rates before
    /// RADIANS/XYZ shaping. Nonfinite time and consumed acquisition errors refuse.
    pub fn fixstar_tt(
        &self,
        star: &StarMatch,
        jd_tt: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
    ) -> Result<([f64; 6], String, i32), Error> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let (place, ret) = self.entry_core(star.star_id, jd_tt, flags, false, topo)?;
        Ok((place, star.canonical.clone(), ret))
    }

    /// Shape an explicit correction recipe without entry-level forcing.
    /// The supplied correction pair, observing center and raw J2000 selection
    /// are taken literally. Acquire position or the geocentric rate stencil,
    /// then apply flags using `j2000_native` to identify the raw frame. No site
    /// is supplied by this recipe. Unknown ids, provider, frame and model errors
    /// propagate; public spherical/Cartesian units are selected by `flags`.
    #[allow(clippy::too_many_arguments)]
    pub fn fixstar_flags_entry(
        &self,
        star_id: i32,
        jd_tt: f64,
        flags: i32,
        j2000_native: bool,
        noaberr: bool,
        nogdefl: bool,
        j2000_frame: bool,
        center_flag: i32,
    ) -> Result<[f64; 6], Error> {
        use crate::constants::FLG_SPEED;
        let raw = if flags & FLG_SPEED != 0 {
            self.sealed_raw_velocity(star_id, jd_tt, noaberr, nogdefl, j2000_frame, center_flag)?
        } else {
            let triple = sealed_star_position(
                self.provider,
                self.nutation_at,
                star_id,
                jd_tt,
                noaberr,
                nogdefl,
                j2000_frame,
                center_flag,
            )?;
            [triple[0], triple[1], triple[2], 0.0, 0.0, 0.0]
        };
        let env = SealedFlagEnv { serving: self };
        apply_fixstar_flags_via(raw, jd_tt, flags, j2000_native, &env)
    }

    /// Stored six-slot place with the geocentric h=0.5-day rate stencil.
    /// The recipe carries no surface site, so no topocentric step is consumed.
    /// Observer/position failures propagate; rates differentiate reported
    /// longitude/latitude degrees and distance AU with the wrapped longitude fold.
    fn sealed_raw_velocity(
        &self,
        star_id: i32,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
        j2000_frame: bool,
        center_flag: i32,
    ) -> Result<[f64; 6], Error> {
        sealed_star_velocity(
            self.provider,
            self.nutation_at,
            star_id,
            jd_tt,
            noaberr,
            nogdefl,
            j2000_frame,
            center_flag,
        )
    }

    /// Dispatch an already-resolved TT/UT1 entry by frame/projection priority.
    /// Fixed-epoch sidereal frames precede active ECL_T 0/SSY_PLANE projections;
    /// responsive live modes suppress those planes. Other calls use the ordinary
    /// entry. `is_ut` selects clock acquisition; `verbatim_echo` selects the TT
    /// wrapper flag echo. The canonical identifier remains owned and unchanged.
    fn wrapped_entry(
        &self,
        star: &StarMatch,
        jd_in: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
        is_ut: bool,
        verbatim_echo: bool,
    ) -> Result<([f64; 6], String, i32), Error> {
        use crate::constants::{FLG_SIDEREAL, SIDBIT_ECL_T0, SIDBIT_SSY_PLANE};
        if (flags & FLG_SIDEREAL) != 0 {
            let mode = self.sidereal.mode();
            if crate::sidereal::is_fixed_epoch_request(flags, mode) {
                return self.fixed_epoch_entry(star, jd_in, flags, topo, is_ut, verbatim_echo);
            }
            let bits = self.sidereal.bits();
            if (bits & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE)) != 0
                && !crate::sidereal::SIDBIT_PROJECTION_SUPPRESS_MODES.contains(&mode)
            {
                return self.sidbit_entry(
                    star,
                    jd_in,
                    flags,
                    topo,
                    is_ut,
                    verbatim_echo,
                    mode,
                    bits,
                );
            }
        }
        if is_ut {
            self.fixstar_ut(star, jd_in, flags, topo)
        } else {
            self.fixstar_tt(star, jd_in, flags, topo)
        }
    }

    /// Evaluate a catalog place in the selected fixed mean epoch frame.
    /// Re-evaluate with internal J2000|NONUT and no SIDEREAL recursion, then
    /// rotate the complete position/rate result to the mode's defining frame.
    /// UT entries echo the caller flags with NONUT; TT entries echo.
    /// Acquisition and transform errors propagate before any result is returned.
    fn fixed_epoch_entry(
        &self,
        star: &StarMatch,
        jd_in: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
        is_ut: bool,
        verbatim_echo: bool,
    ) -> Result<([f64; 6], String, i32), Error> {
        let sub_flags = crate::sidereal::fixed_epoch_request_flags(flags);
        let (sub_place, _, sub_rf) = if is_ut {
            self.fixstar_ut(star, jd_in, sub_flags, topo)
        } else {
            self.fixstar_tt(star, jd_in, sub_flags, topo)
        }?;
        let xx_t0 =
            crate::sidereal::transform_fixed_epoch_result(sub_place, flags, self.sidereal.mode())?;
        let ret = if verbatim_echo {
            flags
        } else {
            crate::sidereal::fixed_epoch_retflag(sub_rf, flags)
        };
        Ok((xx_t0, star.canonical.clone(), ret))
    }

    /// Project an already-resolved star onto a configured sidereal plane.
    /// ECL_T 0 takes precedence over SSY_PLANE. Start from a J2000|NONUT raw place
    /// without SIDEREAL recursion. Equatorial output uses the mean equator and
    /// equinox of t 0 (J2000 for SSY_PLANE) with no zero-point subtraction.
    /// Ecliptic output projects through the zero-point-folded plane matrix.
    /// Responsive live modes are excluded by the enclosing dispatcher.
    ///
    /// ECL_T 0 resolves its defining epoch from the configured sidereal selection;
    /// user-UT and Valens epoch conversion uses the model Delta-T callback.
    /// Ecliptic zero-point evaluation treats t 0 as UT1 and acquires evaluation TT
    /// through the selected star clock. Equatorial output consumes no zero point.
    /// Frame, anchor, model and provider errors propagate; raw anchor ranges use
    /// their original request TT envelope. Units and flag echoes match the wrapped
    /// entry; no fallback plane or anchor is manufactured.
    #[allow(clippy::too_many_arguments)]
    fn sidbit_entry(
        &self,
        star: &StarMatch,
        jd_in: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
        is_ut: bool,
        verbatim_echo: bool,
        mode: i32,
        bits: i32,
    ) -> Result<([f64; 6], String, i32), Error> {
        use crate::constants::{FLG_EQUATORIAL, SIDBIT_ECL_T0};
        use crate::constants::{J2000, SIDBIT_SSY_PLANE};
        debug_assert!((bits & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE)) != 0);
        let sub_flags = crate::sidereal::fixed_epoch_request_flags(flags);
        let (sub_place, _, sub_rf) = if is_ut {
            self.fixstar_ut(star, jd_in, sub_flags, topo)
        } else {
            self.fixstar_tt(star, jd_in, sub_flags, topo)
        }?;
        let anchors = self.anchor_provider();
        let aya = crate::sidereal::SealedAyanamshaServing::new(
            self.sidereal,
            self.deltat_days_at,
            Some(&anchors),
        );
        // Flag precedence mirrors the calc path: `SIDBIT_ECL_T0` selects
        // the mode's defining plane, `SIDBIT_SSY_PLANE` alone the J2000
        // invariable plane.
        let t0 = if (bits & SIDBIT_ECL_T0) != 0 {
            let user_t0 = if mode == crate::constants::SIDM_USER {
                aya.user_pair_tt()?.0
            } else {
                0.0
            };
            let valens = if mode == crate::constants::SIDM_VALENS_MOON {
                aya.valens_epoch_tt()?
            } else {
                0.0
            };
            crate::sidereal::ecl_t0_epoch_jd(mode, user_t0, valens)
        } else {
            J2000
        };
        let xx_p = if (flags & FLG_EQUATORIAL) != 0 {
            crate::sidereal::transform_equatorial_epoch_result(sub_place, flags, t0)?
        } else {
            // `_calc_ayanamsa(t0_jd, sidm)` acquires its evaluation TT
            // through get_timescale.ut 1_jd, independently of deltat.
            // Keep deltat for the user/Valens model anchors inside aya.
            let zero_point_tt = self.jd_tt_for_ut(t0)?;
            let (ayan, _) = aya.mean_tt_with(zero_point_tt, false, false)?;
            let zero = if (bits & SIDBIT_ECL_T0) != 0 {
                ayan
            } else {
                crate::sidereal::ssy_plane_zero_point_deg(ayan)?
            };
            let m_ecl = crate::sidereal::sidbit_ecliptic_matrix(bits, t0, zero)?.ok_or(
                Error::ConfigurationMissing {
                    missing_config: "sidereal_projection",
                    suggestion: "a SIDBIT projection bit is set (guard above)",
                },
            )?;
            crate::sidereal::transform_sidbit_result(sub_place, flags, m_ecl)?
        };
        let ret = if verbatim_echo {
            flags
        } else {
            crate::sidereal::fixed_epoch_retflag(sub_rf, flags)
        };
        Ok((xx_p, star.canonical.clone(), ret))
    }

    /// UT1 entry with fixed-epoch and active-plane dispatch for a resolved star.
    /// Clock conversion and output units follow [`Self::fixstar_ut`]; wrapper
    /// priority follows the common wrapped dispatcher. UT-normalized returned flags and
    /// the supplied owned canonical identifier accompany the six components.
    pub fn fixstar_ut_wrapped(
        &self,
        star: &StarMatch,
        tjdut: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
    ) -> Result<([f64; 6], String, i32), Error> {
        self.wrapped_entry(star, tjdut, flags, topo, true, false)
    }

    /// Scalar TT entry with fixed-epoch and active-plane dispatch.
    /// The supplied catalog match is already resolved. Use the selected frame
    /// prologue or [`Self::fixstar_tt`]; return six components, owned canonical
    /// identifier and the input flags. Consumed errors propagate.
    pub fn fixstar_tt_wrapped(
        &self,
        star: &StarMatch,
        jd_tt: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
    ) -> Result<([f64; 6], String, i32), Error> {
        self.wrapped_entry(star, jd_tt, flags, topo, false, true)
    }

    /// UT1 wrapped entry for a match already resolved with V 2 grammar.
    /// The caller owns flexible lookup; calculation, units, frame priority and
    /// UT echo follow [`Self::fixstar_ut_wrapped`]. No name is reinterpreted here.
    pub fn fixstar2_ut_wrapped(
        &self,
        star: &StarMatch,
        tjdut: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
    ) -> Result<([f64; 6], String, i32), Error> {
        self.wrapped_entry(star, tjdut, flags, topo, true, false)
    }

    /// Scalar TT wrapped entry for a match already resolved with V 2 grammar.
    /// Calculation, units, frame priority and flag echo follow
    /// [`Self::fixstar_tt_wrapped`]. The canonical name identifies the V 2 match.
    pub fn fixstar2_tt_wrapped(
        &self,
        star: &StarMatch,
        jd_tt: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
    ) -> Result<([f64; 6], String, i32), Error> {
        self.wrapped_entry(star, jd_tt, flags, topo, false, true)
    }

    /// UT1 batch places in input order, with optional skipped-error slots.
    ///
    /// A bound [`StarBatchAt`] callback owns the complete batch operation. Its
    /// shared native route resolves names before shared observer acquisition;
    /// shared clock/observer failures remain terminal, even when skipping is
    /// enabled. Topocentric, fixed-epoch and active-plane requests use its
    /// wrapped scalar delegation policy.
    ///
    /// Without that callback, this method checks `tjdut` for finiteness, then
    /// resolves and calculates each name sequentially using V1 grammar. With
    /// `skip_errors=true`, lookup misses and every error returned by that star's
    /// wrapped calculation become `None`, including its clock/provider errors.
    /// Resolver `Refused` errors and the initial nonfinite epoch remain terminal.
    /// With skipping disabled the first consumed error aborts, with no partial
    /// result. Successful slots use [`WrappedStarEntry`] units and layout.
    pub fn batch_fixstars_ut(
        &self,
        stars: &[&str],
        tjdut: f64,
        flags: i32,
        topo: Option<StarTopoSite>,
        skip_errors: bool,
    ) -> Result<Vec<Option<WrappedStarEntry>>, Error> {
        if let Some(batch_at) = self.batch_at {
            return batch_at(stars, tjdut, flags, topo, skip_errors);
        }
        if !tjdut.is_finite() {
            return Err(Error::NonFiniteArgument { value: tjdut });
        }
        let mut out = Vec::with_capacity(stars.len());
        for key in stars {
            let slot = match resolve_star_ref(key) {
                Ok(star) => match self.fixstar_ut_wrapped(&star, tjdut, flags, topo) {
                    Ok(entry) => Some(entry),
                    Err(e) => {
                        if skip_errors {
                            None
                        } else {
                            return Err(e);
                        }
                    }
                },
                Err(StarResolveError::Miss(failure)) => {
                    if skip_errors {
                        None
                    } else {
                        return Err(failure.to_error(key));
                    }
                }
                // raises here instead of
                // returning a miss, escaping every `skip_errors` arm.
                Err(StarResolveError::Refused(e)) => return Err(e),
            };
            out.push(slot);
        }
        Ok(out)
    }

    /// Reduce a genuine already-resolved batch with live per-sample rows.
    /// Shared native acquisition belongs to the bound callback. The stored fast
    /// arm acquires UT1→TT once and evaluates center/minus/plus without scalar
    /// velocity admission. The actual stored catch/log callback selects a whole
    /// vector retry on the same reader; it preserves original Time and partial
    /// stored results. Corruption and unmatched failures remain terminal.
    /// The Session dispatcher handles delegated topocentric/fixed-epoch/plane
    /// requests before selecting this shared fast arm. `row_owners` must have
    /// the same length as `resolved`, with an owner for each resolved star.
    /// This low-level entry trusts that alignment; the Session validates it
    /// before dispatch. Mismatched slices can truncate the paired iteration.
    ///
    /// # Panics
    ///
    /// In the stored route, a resolved star paired with an absent owner panics.
    /// Bound batch callbacks retain their own documented preconditions.
    pub fn batch_resolved_ut(
        &self,
        resolved: &[Option<StarMatch>],
        row_owners: &[Option<std::sync::Arc<dyn StarRowAcquisition>>],
        ut: f64,
        flags_in: i32,
        topo: Option<StarTopoSite>,
        skip: bool,
    ) -> Result<Vec<Option<WrappedStarEntry>>, Error> {
        if resolved.iter().all(Option::is_none) {
            return Ok(vec![None; resolved.len()]);
        }
        if let Some(batch) = self.resolved_anchor_batch_at {
            return batch(
                resolved,
                row_owners,
                ut,
                flags_in,
                topo,
                skip,
                self.dynamic_anchor_provider,
            );
        }
        if let Some(batch) = self.resolved_batch_at {
            return batch(resolved, row_owners, ut, flags_in, topo, skip);
        }
        let ret = fixstar_ret_flags(flags_in, true);
        let flags = preprocess_star_flags(flags_in);
        let time = if self.use_default_timescale {
            let (hi, lo) = crate::smh2016::timescale_parts_from_ut1(ut)?;
            crate::types::TwoPartEpoch { hi, lo }
        } else {
            let jd = self.jd_tt_for_ut(ut)?;
            crate::types::TwoPartEpoch {
                hi: jd.floor(),
                lo: jd - jd.floor(),
            }
        };
        let jd = time.value();
        let center = fixstar_observation_center(flags);
        let noaberr = flags & (crate::constants::FLG_NOABERR | crate::constants::FLG_TRUEPOS) != 0
            || center != 0;
        let nogdefl = flags & (crate::constants::FLG_NOGDEFL | crate::constants::FLG_TRUEPOS) != 0
            || center != 0;
        let j2000 = flags & crate::constants::FLG_J2000 != 0;
        let mut output = Vec::with_capacity(resolved.len());
        for (star, rows) in resolved.iter().zip(row_owners) {
            let Some(star) = star else {
                output.push(None);
                continue;
            };
            let rows = rows.as_ref().expect("validated occurrence owner");
            let serving = SealedStarServing {
                provider: self.provider,
                nutation_at: self.nutation_at,
                deltat_days_at: self.deltat_days_at,
                use_default_timescale: self.use_default_timescale,
                raw_place_at: None,
                batch_at: None,
                resolved_batch_at: None,
                resolved_anchor_batch_at: None,
                rows: Some(rows.as_ref()),
                stored_vector_error_at: self.stored_vector_error_at,
                anchor_at: self.anchor_at,
                anchor_rows: self.anchor_rows,
                acquired_anchor_at: self.acquired_anchor_at,
                dynamic_anchor_provider: self.dynamic_anchor_provider,
                sidereal: self.sidereal,
            };
            // The shared stored arm directly samples positions, without
            // either scalar velocity admission or per-sample fallback.
            let position_at = |epoch| {
                serving.live_stored_position_attempt(
                    rows.as_ref(),
                    star.star_id,
                    epoch,
                    noaberr,
                    nogdefl,
                    j2000,
                    center,
                    topo,
                )
            };
            let result = (|| {
                let raw = if flags & crate::constants::FLG_SPEED != 0 {
                    star_velocity_samples_via(position_at, jd, STAR_GEOCENTRIC_VELOCITY_STEP_DAYS)?
                } else {
                    let p = position_at(jd)?.position;
                    [p[0], p[1], p[2], 0.0, 0.0, 0.0]
                };
                let shaped = serving.shape_native_place(raw, jd, flags, j2000)?;
                Ok((shaped, star.canonical.clone(), ret))
            })();
            match result {
                Ok(entry) => output.push(Some(entry)),
                Err(error) => {
                    if rows.retry_stored(&error, "star batch")? {
                        // Source keeps already assigned stored slots. A
                        // skipped vector retry leaves that earlier value.
                        output.resize(resolved.len(), None);
                        return crate::positions::star_serving::live_stored_vector_batch_via(
                            self.provider,
                            self.stored_vector_error_at,
                            self,
                            resolved,
                            row_owners,
                            output,
                            time,
                            flags_in,
                            skip,
                        );
                    }
                    return Err(error);
                }
            }
        }
        Ok(output)
    }
}

#[cfg(test)]
mod live_row_tests {
    use super::*;
    use std::sync::{
        Arc, Mutex,
        atomic::{AtomicUsize, Ordering},
    };

    struct Row {
        trace: Arc<Mutex<Vec<String>>>,
        serial: usize,
    }

    impl Drop for Row {
        fn drop(&mut self) {
            self.trace
                .lock()
                .unwrap()
                .push(format!("drop:{}", self.serial));
        }
    }

    impl AcquiredStarRow for Row {
        fn native_inputs(&self) -> Result<NativeStarInputs, Error> {
            self.trace
                .lock()
                .unwrap()
                .push(format!("native:{}", self.serial));
            Ok(NativeStarInputs {
                ra_hours: 2.0 * self.serial as f64,
                dec_degrees: 2.0,
                ra_mas_per_year: 0.0,
                dec_mas_per_year: 0.0,
                parallax_mas: 500.0,
                radial_km_per_s: 0.0,
            })
        }
    }

    struct Rows {
        trace: Arc<Mutex<Vec<String>>>,
        next: AtomicUsize,
    }

    impl StarRowAcquisition for Rows {
        fn velocity_admission(&self, _: i32) -> Result<(), Error> {
            self.trace
                .lock()
                .unwrap()
                .push("velocity-membership".into());
            Ok(())
        }
        fn acquire_row(
            &self,
            _: i32,
            backend: StarRowBackend,
        ) -> Result<Arc<dyn AcquiredStarRow>, Error> {
            let serial = self.next.fetch_add(1, Ordering::SeqCst) + 1;
            self.trace
                .lock()
                .unwrap()
                .push(format!("select:{backend:?}:{serial}"));
            Ok(Arc::new(Row {
                trace: Arc::clone(&self.trace),
                serial,
            }))
        }
        fn is_skippable(&self, _: &Error) -> Result<bool, Error> {
            Ok(false)
        }
    }

    struct RefusesObserver(Arc<Mutex<Vec<String>>>);
    impl crate::positions::StateProvider for RefusesObserver {
        fn provider_name(&self) -> &'static str {
            "owned-star-operator-control"
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
        fn eval_raw(&self, body_id: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
            Err(Error::LebBodyNotFound { body_id })
        }
        fn eval_raw_two_part(
            &self,
            body_id: i32,
            _: crate::types::TwoPartEpoch,
        ) -> Result<(Vec<f64>, Vec<f64>), Error> {
            Err(Error::LebBodyNotFound { body_id })
        }
        fn geometric_state(
            &self,
            body_id: i32,
            _: f64,
        ) -> Result<crate::positions::GeometricState, Error> {
            self.0.lock().unwrap().push("observer".into());
            Err(Error::LebBodyNotFound { body_id })
        }
        fn serving_entry_snapshot(
            &self,
            body_id: i32,
        ) -> Result<crate::positions::EntrySnapshot, Error> {
            Err(Error::LebBodyNotFound { body_id })
        }
    }

    #[test]
    fn stored_row_is_selected_before_observer_and_fields_wait_for_success() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let rows = Rows {
            trace: Arc::clone(&trace),
            next: AtomicUsize::new(0),
        };
        let provider = RefusesObserver(Arc::clone(&trace));
        let nutation = |_| panic!("unused frame");
        let delta = |_| Ok(0.0);
        let serving = SealedStarServing::new(
            &provider,
            &nutation,
            &delta,
            crate::sidereal::SiderealSelection::unset(),
        )
        .with_live_rows(&rows);
        let result = serving.raw_place(
            90000001,
            crate::constants::J2000,
            true,
            true,
            true,
            0,
            None,
            false,
        );
        assert!(matches!(
            result,
            Err(Error::LebBodyNotFound {
                body_id: crate::constants::EARTH
            })
        ));
        assert_eq!(
            *trace.lock().unwrap(),
            ["select:Stored:1", "observer", "drop:1"]
        );
    }

    #[test]
    fn scalar_speed_admits_once_and_drops_each_row_before_reselection() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let rows = Rows {
            trace: Arc::clone(&trace),
            next: AtomicUsize::new(0),
        };
        let provider = RefusesObserver(Arc::clone(&trace));
        let nutation = |_| panic!("unused J2000 frame");
        let delta = |_| Ok(0.0);
        let serving = SealedStarServing::new(
            &provider,
            &nutation,
            &delta,
            crate::sidereal::SiderealSelection::unset(),
        )
        .with_live_rows(&rows);
        let raw = serving
            .raw_place(
                90000001,
                crate::constants::J2000,
                true,
                true,
                true,
                crate::constants::FLG_BARYCTR,
                None,
                true,
            )
            .unwrap();
        assert!(
            raw[3].abs() > 1.0,
            "changed rows must affect the actual rate"
        );
        let trace = trace.lock().unwrap();
        let stages = trace
            .iter()
            .filter(|stage| !stage.starts_with("field:"))
            .map(String::as_str)
            .collect::<Vec<_>>();
        assert_eq!(
            stages,
            [
                "velocity-membership",
                "select:Stored:1",
                "native:1",
                "drop:1",
                "select:Stored:2",
                "native:2",
                "drop:2",
                "select:Stored:3",
                "native:3",
                "drop:3"
            ]
        );
    }

    #[test]
    fn stored_shared_speed_has_fresh_rows_without_scalar_membership() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let rows = Arc::new(Rows {
            trace: Arc::clone(&trace),
            next: AtomicUsize::new(0),
        }) as Arc<dyn StarRowAcquisition>;
        let provider = RefusesObserver(Arc::clone(&trace));
        let nutation = |_| panic!("unused J2000 frame");
        let delta = |_| Ok(0.0);
        let serving = SealedStarServing::new(
            &provider,
            &nutation,
            &delta,
            crate::sidereal::SiderealSelection::unset(),
        );
        let resolved = [Some(StarMatch {
            star_id: 90000001,
            canonical: "actual".into(),
        })];
        let output = serving
            .batch_resolved_ut(
                &resolved,
                &[Some(rows)],
                crate::constants::J2000,
                crate::constants::FLG_J2000
                    | crate::constants::FLG_TRUEPOS
                    | crate::constants::FLG_BARYCTR
                    | crate::constants::FLG_SPEED,
                None,
                false,
            )
            .unwrap();
        assert!(output[0].as_ref().unwrap().0[3].abs() > 1.0);
        let trace = trace.lock().unwrap();
        let stages = trace
            .iter()
            .filter(|stage| !stage.starts_with("field:"))
            .map(String::as_str)
            .collect::<Vec<_>>();
        assert_eq!(
            stages,
            [
                "select:Stored:1",
                "native:1",
                "drop:1",
                "select:Stored:2",
                "native:2",
                "drop:2",
                "select:Stored:3",
                "native:3",
                "drop:3"
            ]
        );
    }

    #[test]
    fn native_keyword_units_preserve_original_constructor_words() {
        let data = fixed_star_data(crate::constants::SIRIUS).unwrap();
        let inputs = NativeStarInputs {
            ra_hours: data.ra_j2000_deg() / 15.0,
            dec_degrees: data.dec_j2000_deg(),
            ra_mas_per_year: data.pm_ra_arcsec_yr() * 1000.0,
            dec_mas_per_year: data.pm_dec_arcsec_yr() * 1000.0,
            parallax_mas: data.parallax_mas(),
            radial_km_per_s: data.radial_velocity_km_s(),
        };
        let (actual_pos, actual_vel) = star_initial_state_from_native(inputs);
        let (old_pos, old_vel) = star_icrs_initial_state_catalog_angles(&data);
        for (actual, old) in actual_pos
            .into_iter()
            .chain(actual_vel)
            .zip(old_pos.into_iter().chain(old_vel))
        {
            assert_eq!(actual.to_bits(), old.to_bits());
        }
    }

    struct RefusedRow;
    impl AcquiredStarRow for RefusedRow {
        fn native_inputs(&self) -> Result<NativeStarInputs, Error> {
            Err(Error::UnknownBody { body_id: 77 })
        }
    }
    struct RetryRows {
        trace: Arc<Mutex<Vec<String>>>,
    }
    impl StarRowAcquisition for RetryRows {
        fn velocity_admission(&self, _: i32) -> Result<(), Error> {
            panic!("unused shared velocity admission")
        }
        fn acquire_row(
            &self,
            _: i32,
            backend: StarRowBackend,
        ) -> Result<Arc<dyn AcquiredStarRow>, Error> {
            self.trace
                .lock()
                .unwrap()
                .push(format!("retry-select:{backend:?}"));
            if backend == StarRowBackend::Stored {
                return Ok(Arc::new(RefusedRow));
            }
            Ok(Arc::new(Row {
                trace: Arc::clone(&self.trace),
                serial: 2,
            }))
        }
        fn is_skippable(&self, _: &Error) -> Result<bool, Error> {
            Ok(false)
        }
        fn retry_stored(&self, error: &Error, context: &'static str) -> Result<bool, Error> {
            assert!(matches!(error, Error::UnknownBody { body_id: 77 }));
            self.trace.lock().unwrap().push(format!("log:{context}"));
            Ok(true)
        }
    }

    #[test]
    fn stored_scalar_retry_reacquires_native_row_without_jpl_or_inverse_time() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let rows = RetryRows {
            trace: Arc::clone(&trace),
        };
        let provider = RefusesObserver(Arc::clone(&trace));
        let unused = |_| panic!("no model/inverse clock in BARYCTR J2000 constructor");
        let nutation = |_| panic!("unused frame");
        let serving = SealedStarServing::new(
            &provider,
            &nutation,
            &unused,
            crate::sidereal::SiderealSelection::unset(),
        )
        .with_default_timescale();
        let result = serving
            .live_stored_position(
                &rows,
                90000001,
                crate::constants::J2000,
                true,
                true,
                true,
                crate::constants::FLG_BARYCTR,
                None,
            )
            .unwrap();
        assert!(result.position[2] > 0.0);
        assert_eq!(
            *trace.lock().unwrap(),
            [
                "retry-select:Stored",
                "log:star",
                "retry-select:Native",
                "native:2",
                "drop:2"
            ]
        );
    }

    #[test]
    fn huge_finite_stored_epoch_reaches_row_refusal_without_unused_inverse() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let rows = RetryRows {
            trace: Arc::clone(&trace),
        };
        let provider = RefusesObserver(Arc::clone(&trace));
        let unused = |_| panic!("unused clock");
        let nutation = |_| panic!("unused frame");
        let serving = SealedStarServing::new(
            &provider,
            &nutation,
            &unused,
            crate::sidereal::SiderealSelection::unset(),
        )
        .with_default_timescale();
        let result = serving.live_stored_position_attempt(
            &rows,
            90000001,
            1e200,
            true,
            true,
            true,
            crate::constants::FLG_BARYCTR,
            None,
        );
        assert!(matches!(result, Err(Error::UnknownBody { body_id: 77 })));
    }

    #[test]
    fn stored_batch_retry_recomputes_successes_and_retains_original_time_once() {
        let trace = Arc::new(Mutex::new(Vec::new()));
        let first = Arc::new(Rows {
            trace: Arc::clone(&trace),
            next: AtomicUsize::new(0),
        }) as Arc<dyn StarRowAcquisition>;
        let second = Arc::new(RetryRows {
            trace: Arc::clone(&trace),
        }) as Arc<dyn StarRowAcquisition>;
        let provider = RefusesObserver(Arc::clone(&trace));
        let delta_calls = std::cell::Cell::new(0);
        let delta = |_| {
            delta_calls.set(delta_calls.get() + 1);
            Ok(0.0)
        };
        let nutation = |_| panic!("unused frame");
        let serving = SealedStarServing::new(
            &provider,
            &nutation,
            &delta,
            crate::sidereal::SiderealSelection::unset(),
        );
        let resolved = ["first", "second"].map(|canonical| {
            Some(StarMatch {
                star_id: 90000001,
                canonical: canonical.into(),
            })
        });
        let output = serving
            .batch_resolved_ut(
                &resolved,
                &[Some(first), Some(second)],
                crate::constants::J2000,
                crate::constants::FLG_J2000
                    | crate::constants::FLG_TRUEPOS
                    | crate::constants::FLG_BARYCTR,
                None,
                false,
            )
            .unwrap();
        assert!(output.iter().all(Option::is_some));
        assert_eq!(
            delta_calls.get(),
            1,
            "retry must reuse the already acquired batch Time"
        );
        let first = &output[0].as_ref().unwrap().0;
        let second = &output[1].as_ref().unwrap().0;
        assert_eq!(
            first, second,
            "both retried rows use the same native catalog data"
        );
        assert!(first.iter().all(|value| value.is_finite()));
        assert!(first[2] > 0.0);
        assert_eq!(&first[3..], &[0.0; 3]);
        let events = trace.lock().unwrap();
        let snapshots = events
            .iter()
            .filter(|event| event.starts_with("native:"))
            .count();
        let retired = events
            .iter()
            .filter(|event| event.starts_with("drop:"))
            .count();
        assert_eq!(
            snapshots, 3,
            "one stored snapshot and two replacement snapshots"
        );
        assert_eq!(retired, snapshots, "all owned numerical rows are released");
        assert_eq!(
            events
                .iter()
                .filter(|event| *event == "log:star batch")
                .count(),
            1
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::star_catalog_tables::{REF_SORTED_IDS, STAR_CATALOG};

    #[test]
    fn catalog_ids_stay_in_fixed_star_namespace_and_preserve_entries() {
        for entry in STAR_CATALOG {
            assert!(
                entry.id() >= crate::constants::FIXSTAR_OFFSET,
                "catalog identifier {} is outside the fixed-star namespace",
                entry.id()
            );
            let resolved = star_entry_by_id(entry.id()).expect("generated catalog row");
            assert_eq!(resolved, entry, "catalog row {}", entry.id());
        }
    }

    #[test]
    fn catalog_id_domain_guard_preserves_unknown_and_nonstellar_misses() {
        use crate::constants::{AST_OFFSET, EARTH, FIXSTAR_OFFSET, SUN};

        let beyond_catalog = STAR_CATALOG
            .iter()
            .map(StarCatalogEntry::id)
            .max()
            .expect("nonempty catalog")
            .checked_add(1)
            .expect("catalog identifiers leave an upper-domain miss");
        for star_id in [
            i32::MIN,
            -1,
            SUN,
            EARTH,
            AST_OFFSET,
            FIXSTAR_OFFSET - 1,
            FIXSTAR_OFFSET,
            beyond_catalog,
            i32::MAX,
        ] {
            let original_match = STAR_CATALOG.iter().find(|entry| entry.id() == star_id);
            if star_id < FIXSTAR_OFFSET || star_id >= beyond_catalog {
                assert!(original_match.is_none(), "expected catalog miss {star_id}");
            }
            assert_eq!(star_entry_by_id(star_id), original_match, "id {star_id}");
            assert_eq!(
                fixed_star_data(star_id),
                original_match.map(StarCatalogEntry::data),
                "id {star_id}"
            );
            assert_eq!(
                star_magnitude(star_id),
                original_match.map(StarCatalogEntry::magnitude),
                "id {star_id}"
            );
            assert_eq!(
                canonical_star_name(star_id),
                original_match.map(StarCatalogEntry::name),
                "id {star_id}"
            );
        }
    }

    #[test]
    fn catalog_row_count_and_identity() {
        // The catalog contains 1447 rows.
        assert_eq!(STAR_CATALOG.len(), 1447);
        assert_eq!(list_fixed_stars().len(), 1447);
        let regulus = star_entry_by_id(1_000_001).expect("Regulus row");
        assert_eq!(regulus.name(), "Regulus");
        assert_eq!(regulus.nomenclature(), "alLeo");
        assert_eq!(regulus.hip_number(), 49669);
        assert_eq!(canonical_star_name(1_000_001), Some("Regulus"));
        assert_eq!(canonical_star_name(999), None);
    }

    #[test]
    fn table_ids_unique_and_alias_targets_resolve() {
        let mut ids: Vec<i32> = STAR_CATALOG.iter().map(|e| e.id()).collect();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), STAR_CATALOG.len());
        for (_, id) in STAR_ALIASES {
            assert!(
                star_entry_by_id(*id).is_some(),
                "alias target {id} missing from catalog"
            );
        }
        for id in REF_SORTED_IDS.iter().chain(NOMEN_SORTED_IDS.iter()) {
            assert!(star_entry_by_id(*id).is_some());
        }
        let mut sorted = ids.clone();
        sorted.sort_unstable();
        assert_eq!(sorted.len(), STAR_CATALOG.len());
    }

    #[test]
    fn emitted_sort_orders_match_rust_key_functions() {
        // The generator computed these with the model key functions; the
        // Rust keys must reproduce them over every row.
        let mut by_ref: Vec<i32> = STAR_CATALOG.iter().map(|e| e.id()).collect();
        by_ref.sort_by_key(|id| {
            ref_star_key(star_entry_by_id(*id).expect("sorted id in catalog").name())
        });
        assert_eq!(by_ref, REF_SORTED_IDS);
        let mut by_nomen: Vec<i32> = STAR_CATALOG.iter().map(|e| e.id()).collect();
        by_nomen.sort_by_key(|id| {
            nomen_stripped(
                star_entry_by_id(*id)
                    .expect("sorted id in catalog")
                    .nomenclature(),
            )
        });
        assert_eq!(by_nomen, NOMEN_SORTED_IDS);
    }

    #[test]
    fn phonetic_normalization_pins() {
        // Pinned against the model code itself (extracted function,
        // scratch only): the docstring
        // examples (`BETLGUS`, `FOMALHAT`) are stale — the code yields the
        // values below, which the port reproduces exactly.
        assert_eq!(normalize_phonetic("Betelgeuse"), "BTLGS");
        assert_eq!(normalize_phonetic("Betelgeux"), "BTLGKS");
        assert_eq!(normalize_phonetic("Fomalhaut"), "FMLHT");
        assert_eq!(normalize_phonetic("Formalhaut"), "FRMLHT");
        assert_eq!(normalize_phonetic("Aldebaran"), "ALDBRN");
    }

    #[test]
    fn bayer_and_flamsteed_docstring_pins() {
        assert_eq!(
            parse_bayer_designation("Alpha Leonis").as_deref(),
            Some("alLeo")
        );
        assert_eq!(
            parse_bayer_designation("Beta Persei").as_deref(),
            Some("bePer")
        );
        assert_eq!(
            parse_bayer_designation("Gamma Virginis").as_deref(),
            Some("gaVir")
        );
        assert_eq!(parse_bayer_designation("Invalid Name"), None);
        assert_eq!(parse_bayer_designation(""), None);
        assert_eq!(
            parse_flamsteed_designation("32 Leonis").as_deref(),
            Some("32 LEO")
        );
        assert_eq!(
            parse_flamsteed_designation("87 Virginis").as_deref(),
            Some("87 VIR")
        );
        assert_eq!(
            parse_flamsteed_designation("21 Tauri").as_deref(),
            Some("21 TAU")
        );
        assert_eq!(parse_flamsteed_designation("Invalid Name"), None);
    }

    #[test]
    fn library_resolve_name_tiers() {
        assert_eq!(resolve_star_name("Regulus"), Some(1_000_001));
        assert_eq!(resolve_star_name("SIRIUS"), Some(crate::constants::SIRIUS));
        assert_eq!(resolve_star_name("Alpha Leonis"), Some(1_000_001));
        assert_eq!(resolve_star_name("32 Leonis"), Some(1_000_001));
        assert_eq!(resolve_star_name(""), None);
        assert_eq!(resolve_star_name(","), None);
        // Misspellings reach the fuzzy tier like the model.
        assert_eq!(
            resolve_star_name("Betelgeux"),
            resolve_star_name("Betelgeuse")
        );
        assert!(resolve_star_name("NoSuchStarXyz").is_none());
    }

    fn miss_message(result: Result<StarMatch, StarResolveError>) -> String {
        match result {
            Err(StarResolveError::Miss(failure)) => failure.message,
            _ => panic!("expected a catalog miss"),
        }
    }

    #[test]
    fn v1_resolver_pins() {
        let regulus = resolve_star_ref("Regulus").expect("Regulus");
        assert_eq!(regulus.star_id, 1_000_001);
        assert_eq!(regulus.canonical, "Regulus,alLeo");
        // Shared designation resolves to the primary (brightest) component.
        let al_cen = resolve_star_ref(",alCen").expect(",alCen");
        assert!(al_cen.canonical.starts_with("Rigil Kentaurus,"));
        // Name-before-comma wins; unknown name fails with the original string.
        assert_eq!(
            resolve_star_ref("Regulus,zzZzz")
                .expect("name wins")
                .star_id,
            1_000_001
        );
        assert_eq!(
            miss_message(resolve_star_ref("Nosuch,alTau")),
            "star Nosuch,alTau not found"
        );
        // `%` is an ordinary name character in v 1.
        assert_eq!(
            miss_message(resolve_star_ref("Alde%")),
            "star Alde% not found"
        );
        assert_eq!(miss_message(resolve_star_ref("")), "star name empty");
        // Sequential numbers index the name-sorted catalog.
        let first = resolve_star_ref("1").expect("seq 1");
        assert_eq!(first.star_id, REF_SORTED_IDS[0]);
        assert!(resolve_star_ref("0").is_err());
        assert!(resolve_star_ref("1448").is_err());
    }

    #[test]
    fn v2_resolver_pins() {
        let regulus = resolve_star2("Regulus").expect("Regulus");
        assert_eq!(regulus.star_id, 1_000_001);
        // Every comma form keys on the nomenclature: Regulus,alVir -> Spica.
        let spica = resolve_star2("Regulus,alVir").expect("nomenclature wins");
        assert!(spica.canonical.starts_with("Spica,"));
        // Bare nomenclature refuses.
        assert_eq!(
            miss_message(resolve_star2("alCen")),
            "could not find star name alcen"
        );
        // Wildcard grammar.
        assert!(resolve_star2("Alde%").is_ok());
        assert_eq!(
            miss_message(resolve_star2("a%b")),
            "invalid search string a%b"
        );
        assert_eq!(
            miss_message(resolve_star2("0")),
            "could not find star name 0"
        );
        assert_eq!(
            miss_message(resolve_star2("99999")),
            "sequential fixed star number 99999 is not available"
        );
        // Curated HIP-keyed correction: Alaraph -> beta Vir (Zavijava).
        let alaraph = resolve_star2("Alaraph").expect("Alaraph");
        let zavijava = STAR_CATALOG.iter().find(|e| e.hip_number == 57757).unwrap();
        assert_eq!(alaraph.star_id, zavijava.id());
        assert_eq!(miss_message(resolve_star2("")), "star name empty");
    }

    #[test]
    fn search_type_taxonomy() {
        assert_eq!(star_search_type("", StarFamily::V1), "empty");
        assert_eq!(star_search_type("Alde%", StarFamily::V1), "name");
        assert_eq!(star_search_type("Alde%", StarFamily::V2), "wildcard");
        assert_eq!(star_search_type(",alCen", StarFamily::V1), "nomenclature");
        assert_eq!(star_search_type("12", StarFamily::V2), "sequential");
        assert_eq!(star_search_type("Regulus", StarFamily::V1), "name");
        assert_eq!(star_search_type("10CMa", StarFamily::V1), "name");
        assert_eq!(star_search_type("109Her", StarFamily::V2), "name");
    }

    #[test]
    fn every_catalog_name_resolves_to_its_own_physical_entry() {
        for entry in STAR_CATALOG {
            for resolver in [resolve_star_ref, resolve_star2] {
                let actual = resolver(entry.name).expect("catalog names must be reachable");
                assert_eq!(actual.star_id, entry.id, "{}", entry.name);
            }
        }
        // Unknown numeric suffixes retain their public sequential grammar.
        assert_eq!(
            resolve_star_ref("10unknown").unwrap().star_id,
            REF_SORTED_IDS[9]
        );
        assert_eq!(
            resolve_star2("10unknown").unwrap().star_id,
            NOMEN_SORTED_IDS[9]
        );
    }

    #[test]
    fn sequential_indices_use_unicode_decimal_digits_and_bounded_arithmetic() {
        assert_eq!(resolve_star_ref("١").unwrap().star_id, REF_SORTED_IDS[0]);
        assert_eq!(resolve_star2("１").unwrap().star_id, NOMEN_SORTED_IDS[0]);
        let leading_zeros = format!("{}1", "0".repeat(20_000));
        assert_eq!(
            resolve_star_ref(&leading_zeros).unwrap().star_id,
            REF_SORTED_IDS[0]
        );
        for query in ["²", "99999999999999999999999"] {
            assert!(matches!(
                resolve_star_ref(query),
                Err(StarResolveError::Miss(_))
            ));
            assert!(matches!(
                resolve_star2(query),
                Err(StarResolveError::Miss(_))
            ));
        }
        assert_eq!(star_search_type("²", StarFamily::V1), "name");
        assert_eq!(resolve_star2("1²").unwrap().star_id, NOMEN_SORTED_IDS[0]);
        assert_eq!(parse_flamsteed_designation("² Leonis"), None);
    }

    #[test]
    fn hip_lookup_pins() {
        assert_eq!(hip_for_star_name("Regulus"), Some(49669));
        assert_eq!(hip_for_star_name("alpha leonis"), Some(49669));
        assert_eq!(hip_for_star_name("Unknown Star"), None);
        assert_eq!(hip_for_star_name(""), None);
        // Mapped star without a HIP number keeps the model -1 sentinel.
        assert_eq!(hip_for_star_name("ABSOLUTNO"), Some(-1));
    }

    #[test]
    fn propagate_pm_sirius_docstring_pin() {
        // `propagate_proper_motion` docstring example (Sirius at
        // J 1991.25 -> J2000), values captured from an ephemeral pinned-tag
        // probe (Rust-owned scratch only, never stored).
        let (ra, dec) = propagate_proper_motion(
            101.289,
            -16.713,
            -0.54601,
            -1.22307,
            crate::constants::J1991_25,
            crate::constants::J2000,
        );
        assert_eq!(ra.to_bits(), 101.28761435910722f64.to_bits());
        assert_eq!(dec.to_bits(), (-16.715972739583336f64).to_bits());
    }

    #[test]
    fn propagate_pm_pole_clamp_and_wrap_edges() {
        use crate::constants::J2000;
        // Pole: RA undefined, kept unchanged (same ephemeral probe).
        assert_eq!(
            propagate_proper_motion(10.0, 90.0, 1.0, 0.5, J2000, J2000 + 36525.0),
            (10.0, 90.0)
        );
        // Declination clamps at +90 (same ephemeral probe).
        assert_eq!(
            propagate_proper_motion(0.0, 89.0, 0.0, 36000.0, J2000, J2000 + 3652.5),
            (0.0, 90.0)
        );
        // Right ascension wraps at 360 degrees.
        assert_eq!(
            propagate_proper_motion(359.0, 0.0, 3.6, 0.0, J2000, J2000 + 36525.0),
            (359.1, 0.0)
        );
        assert_eq!(
            propagate_proper_motion(1.0, 0.0, -3.6, 0.0, J2000, J2000 + 36525.0),
            (0.9, 0.0)
        );
        // Southern clamp arm.
        let (_, dec) = propagate_proper_motion(0.0, -89.0, 0.0, -36000.0, J2000, J2000 + 3652.5);
        assert_eq!(dec, -90.0);
    }

    #[test]
    fn icrs_state_linearity_and_parallax_distance() {
        use crate::constants::J2000;
        let regulus = fixed_star_data(1_000_001).expect("Regulus");
        let (pos0, vel) = star_icrs_state(&regulus, J2000);
        // At J2000 the propagation term vanishes: the state is the
        // parallax distance on the (ra, dec) direction.
        let dist = (pos0[0] * pos0[0] + pos0[1] * pos0[1] + pos0[2] * pos0[2]).sqrt();
        let expect = 1.0 / (regulus.parallax_mas() * 1.0e-3 * STAR_ASEC2RAD).sin();
        assert_eq!(dist.to_bits(), expect.to_bits());
        // Linearity: one Julian year later the state is pos 0 + vel, with
        // no t^2 term (the pinned evaluated path is first-order in 3D).
        let year = crate::constants::DAYS_PER_JULIAN_YEAR;
        let (pos1, vel1) = star_icrs_state(&regulus, J2000 + year);
        for i in 0..3 {
            assert_eq!(pos1[i].to_bits(), (pos0[i] + vel[i] * year).to_bits());
            assert_eq!(vel1[i].to_bits(), vel[i].to_bits());
        }
    }

    #[test]
    fn parallax_floor_is_one_gigaparsec() {
        use crate::constants::J2000;
        // theta-1 Ori A carries a non-positive parallax at the pin
        // (unknown distance); the model `> 0.0` gate floors it.
        let star = fixed_star_data(1_000_416).expect("th01OriA");
        assert!(star.parallax_mas() <= 0.0);
        let (pos, _vel) = star_icrs_state(&star, J2000);
        let dist = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
        let expect = 1.0 / (STAR_PARALLAX_FLOOR_MAS * 1.0e-3 * STAR_ASEC2RAD).sin();
        assert_eq!(dist.to_bits(), expect.to_bits());
    }

    #[test]
    fn geocentric_zero_coincidence_has_no_light_time() {
        use crate::constants::J2000;
        let regulus = fixed_star_data(1_000_001).expect("Regulus");
        let (pos, _vel) = star_icrs_state(&regulus, J2000);
        let (geo, dist, lt) = star_geocentric_icrs(&regulus, J2000, pos);
        assert_eq!(geo, [0.0, 0.0, 0.0]);
        assert_eq!(dist, 0.0);
        assert_eq!(lt, 0.0);
    }

    #[test]
    fn geometric_j2000_by_id_refuses_unknown_stars() {
        use crate::constants::J2000;
        match star_geometric_j2000_by_id(0, J2000, [0.0, 0.0, 0.0]) {
            Err(Error::StarNotFound { .. }) => {}
            other => panic!("expected StarNotFound, got {other:?}"),
        }
        // Every catalog id resolves (the sweep pins the values).
        for id in 1_000_001..1_000_448 {
            assert!(
                star_geometric_j2000_by_id(id, J2000, [1.0, 0.0, 0.0]).is_ok(),
                "id {id}"
            );
        }
    }

    #[test]
    fn star_flag_and_center_helpers_match_release() {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_MOSEPH, FLG_SPEED, FLG_SPEED3, FLG_SWIEPH, FLG_TOPOCTR,
        };
        // `_preprocess_flags`: MOSEPH stripped, SPEED 3 -> SPEED, TOPOCTR kept.
        assert_eq!(preprocess_star_flags(FLG_SWIEPH), FLG_SWIEPH);
        assert_eq!(preprocess_star_flags(FLG_MOSEPH), 0);
        assert_eq!(preprocess_star_flags(FLG_SWIEPH | FLG_MOSEPH), FLG_SWIEPH);
        assert_eq!(preprocess_star_flags(FLG_SPEED3), FLG_SPEED);
        assert_eq!(preprocess_star_flags(FLG_SPEED3 | FLG_SPEED), FLG_SPEED);
        assert_eq!(
            preprocess_star_flags(FLG_TOPOCTR | FLG_SWIEPH),
            FLG_TOPOCTR | FLG_SWIEPH
        );
        // `_fixstar_observation_center`: TOPOCTR > BARYCTR > HELCTR > geo.
        assert_eq!(fixstar_observation_center(0), 0);
        assert_eq!(fixstar_observation_center(FLG_HELCTR), FLG_HELCTR);
        assert_eq!(fixstar_observation_center(FLG_BARYCTR), FLG_BARYCTR);
        assert_eq!(
            fixstar_observation_center(FLG_HELCTR | FLG_BARYCTR),
            FLG_BARYCTR
        );
        assert_eq!(
            fixstar_observation_center(FLG_TOPOCTR | FLG_BARYCTR | FLG_HELCTR),
            0
        );
    }

    #[test]
    fn apparent_geometric_branch_equals_star03_geometric() {
        use crate::constants::J2000;
        // noaberr + nogdefl + J2000 is the geometric path: the
        // shared core must reproduce it exactly (same operation order,
        // zero observer velocity keeps aberration an identity anyway).
        let regulus = fixed_star_data(1_000_001).expect("Regulus");
        let observer = StarObserver::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]).expect("observer");
        let no_deflect = |_: i32, _: f64| -> Result<[f64; 3], Error> {
            panic!("no deflectors consulted on the geometric branch")
        };
        for jd in [J2000, J2000 + 1000.0] {
            let apparent = star_apparent_ecliptic_via(
                &regulus, jd, &observer, true, true, true, no_deflect, None,
            )
            .expect("geometric branch");
            let geometric = star_geometric_j2000_ecliptic(&regulus, jd, [1.0, 0.0, 0.0])
                .expect("STAR-03 geometric");
            assert_eq!(apparent, geometric, "jd {jd}");
        }
    }

    #[test]
    fn apparent_correction_branches_are_independent() {
        use crate::constants::J2000;
        let regulus = fixed_star_data(1_000_001).expect("Regulus");
        let observer = StarObserver::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]).expect("observer");
        // A skipped deflector never
        // perturbs the branch that skips deflection.
        let missing = |_: i32, _: f64| -> Result<[f64; 3], Error> {
            Err(Error::LebBodyNotFound { body_id: 0 })
        };
        let skipped =
            star_apparent_ecliptic_via(&regulus, J2000, &observer, true, true, true, missing, None)
                .expect("skipped deflectors continue");
        let geometric = star_geometric_j2000_ecliptic(&regulus, J2000, [1.0, 0.0, 0.0])
            .expect("STAR-03 geometric");
        assert_eq!(skipped, geometric);
        // A fatal deflector error propagates when deflection runs and is
        // ignored when `nogdefl` skips it (never a silent substitution).
        let fatal = |_: i32, _: f64| -> Result<[f64; 3], Error> {
            Err(Error::NonFiniteArgument { value: f64::NAN })
        };
        assert!(matches!(
            star_apparent_ecliptic_via(&regulus, J2000, &observer, true, false, true, fatal, None),
            Err(Error::NonFiniteArgument { .. })
        ));
        assert!(
            star_apparent_ecliptic_via(&regulus, J2000, &observer, true, true, true, fatal, None)
                .is_ok()
        );
        // Zero observer velocity is the aberration identity: toggling
        // `noaberr` changes nothing there, pinning the independent gate.
        let with_aber = star_apparent_ecliptic_via(
            &regulus, J2000, &observer, false, true, true, missing, None,
        )
        .expect("aberration identity");
        assert_eq!(with_aber, geometric);
    }

    #[test]
    fn direction_corrections_preserve_geometric_range_and_radial_rate() {
        let initial = ([1e7, 2e7, 3e7], [0.001, -0.002, 0.003]);
        let observer = StarObserver::new([1.0, 0.0, 0.0], [0.0, 0.017, 0.0]).unwrap();
        let place = |noaberr, nogdefl| {
            star_apparent_ecliptic_initial_via(
                initial,
                crate::constants::J2000,
                &observer,
                noaberr,
                nogdefl,
                true,
                |_, _| Ok([0.0; 3]),
                || Ok(None),
            )
            .unwrap()
        };
        let geometric = place(true, true);
        for (noaberr, nogdefl) in [(true, false), (false, true), (false, false)] {
            let apparent = place(noaberr, nogdefl);
            assert_eq!(apparent.position[2], geometric.position[2]);
            assert_eq!(apparent.radial_rate, geometric.radial_rate);
            assert_ne!(apparent.position[..2], geometric.position[..2]);
        }
    }

    #[test]
    fn apparent_ofdate_frame_wiring_and_missing_frame_refusal() {
        use crate::constants::J2000;
        let regulus = fixed_star_data(1_000_001).expect("Regulus");
        let observer = StarObserver::new([1.0, 0.0, 0.0], [0.0, 0.0, 0.0]).expect("observer");
        let missing = |_: i32, _: f64| -> Result<[f64; 3], Error> {
            Err(Error::LebBodyNotFound { body_id: 0 })
        };
        // Of-date without a frame refuses (never a silent J2000).
        assert!(matches!(
            star_apparent_ecliptic_via(
                &regulus, J2000, &observer, true, true, false, missing, None
            ),
            Err(Error::ConfigurationMissing { .. })
        ));
        // Identity-matrix wiring: pn = I reproduces the manual
        // mat 3 + equatorial-to-ecliptic composition exactly.
        let eps = 23.43928f64.to_radians();
        let frame = StarFrame::new([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], eps)
            .expect("identity frame");
        let got = star_apparent_ecliptic_via(
            &regulus,
            J2000,
            &observer,
            true,
            true,
            false,
            missing,
            Some(frame),
        )
        .expect("of-date wiring");
        let (geo, _, _) = star_geocentric_icrs(&regulus, J2000, [1.0, 0.0, 0.0]);
        let expected = crate::positions::cartesian_to_spherical_deg(
            crate::positions::rotate_equatorial_to_ecliptic(geo, eps).expect("rotate"),
        )
        .expect("shape");
        assert_eq!(got, expected);
        // Checked constructors refuse non-finite inputs.
        assert!(StarObserver::new([f64::NAN, 0.0, 0.0], [0.0; 3]).is_err());
        assert!(StarFrame::new([[1.0, 0.0, 0.0]; 3], f64::INFINITY).is_err());
    }

    #[test]
    fn star_velocity_wraps_and_scales() {
        // Longitude straddling 0/360 takes the shortest arc.
        let crossing = |jd: f64| -> Result<[f64; 3], Error> {
            Ok(if jd < 0.0 {
                [359.0, 0.0, 1.0]
            } else {
                [1.0, 0.0, 1.0]
            })
        };
        let v = star_velocity_via(crossing, 0.0, 0.5).expect("wrap");
        assert_eq!(v[0].to_bits(), 1.0f64.to_bits());
        assert_eq!(v[3].to_bits(), 2.0f64.to_bits());
        assert_eq!(v[4].to_bits(), 0.0f64.to_bits());
        // Linear drift differentiates exactly.
        let drift =
            |jd: f64| -> Result<[f64; 3], Error> { Ok([10.0 + 0.1 * jd, 5.0 - 0.2 * jd, 3.0]) };
        let v = star_velocity_via(drift, 100.0, 0.5).expect("drift");
        assert!((v[3] - 0.1).abs() < 1e-12);
        assert!((v[4] + 0.2).abs() < 1e-12);
        assert_eq!(v[5].to_bits(), 0.0f64.to_bits());
        // Non-positive steps refuse (never a silent divide-by-zero).
        assert!(star_velocity_via(drift, 100.0, 0.0).is_err());
        assert!(star_velocity_via(drift, f64::NAN, 0.5).is_err());
    }

    #[test]
    fn flag_ret_flags_pins() {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_MOSEPH, FLG_NOABERR, FLG_NOGDEFL, FLG_NONUT,
            FLG_SIDEREAL, FLG_SWIEPH, FLG_TOPOCTR,
        };
        // TT entries echo.
        assert_eq!(fixstar_ret_flags(0, false), 0);
        assert_eq!(
            fixstar_ret_flags(FLG_MOSEPH | FLG_SIDEREAL, false),
            FLG_MOSEPH | FLG_SIDEREAL
        );
        // UT entries default the ephemeris bit.
        assert_eq!(fixstar_ret_flags(0, true), FLG_SWIEPH);
        assert_eq!(fixstar_ret_flags(FLG_SWIEPH, true), FLG_SWIEPH);
        // MOSEPH is accepted-for-compatibility and echoed when alone.
        assert_eq!(fixstar_ret_flags(FLG_MOSEPH, true), FLG_MOSEPH);
        // J2000/SIDEREAL imply NONUT on the echo only.
        assert_eq!(
            fixstar_ret_flags(FLG_SWIEPH | FLG_J2000, true),
            FLG_SWIEPH | FLG_J2000 | FLG_NONUT
        );
        assert_eq!(
            fixstar_ret_flags(FLG_SWIEPH | FLG_SIDEREAL, true),
            FLG_SWIEPH | FLG_SIDEREAL | FLG_NONUT
        );
        // Helio/bary imply NOGDEFL|NOABERR on the echo only.
        assert_eq!(
            fixstar_ret_flags(FLG_SWIEPH | FLG_HELCTR, true),
            FLG_SWIEPH | FLG_HELCTR | FLG_NOGDEFL | FLG_NOABERR
        );
        // Center priority strips losers from the echo too.
        assert_eq!(
            fixstar_ret_flags(FLG_SWIEPH | FLG_TOPOCTR | FLG_BARYCTR | FLG_HELCTR, true)
                & (FLG_HELCTR | FLG_BARYCTR),
            0
        );
    }

    #[test]
    fn precess_ecliptic_identity_and_edges() {
        use crate::constants::J2000;
        // Inside the model `|from - to| < 1e-6` window the input is
        // returned untouched (exactly).
        assert_eq!(
            precess_ecliptic_to_j2000(123.456, -20.5, J2000).expect("identity"),
            (123.456, -20.5)
        );
        assert_eq!(
            precess_ecliptic_to_j2000(123.456, -20.5, J2000 + 5e-7).expect("window"),
            (123.456, -20.5)
        );
        // Outside the window the frame moves (24 years carry ~0.34 deg
        // of precession) while latitude stays bounded.
        let (lon, lat) = precess_ecliptic_to_j2000(100.0, 20.0, J2000 + 8766.0).expect("precess");
        assert!((lon - 100.0).abs() > 0.1 && (lon - 100.0).abs() < 1.0);
        assert!((lat - 20.0).abs() < 0.5);
        assert!((0.0..360.0).contains(&lon));
        // Non-finite inputs refuse.
        assert!(precess_ecliptic_to_j2000(f64::NAN, 0.0, J2000).is_err());
        assert!(precess_ecliptic_to_j2000(0.0, 0.0, f64::INFINITY).is_err());
    }

    /// Stub flag environment: every provider returns a canned value and
    /// records which arms fired, so each `_positional_transform` branch is
    /// pinned independently of any ephemeris.
    struct StubEnv {
        precess_to: Option<(f64, f64)>,
        dpsi_rad: f64,
        mean_eps: f64,
        true_eps: f64,
        ut1_shift: f64,
        ayan: f64,
        seen_ex_flags: std::cell::Cell<i32>,
        seen_precess: std::cell::Cell<bool>,
    }

    impl StubEnv {
        fn plain() -> Self {
            Self {
                precess_to: None,
                dpsi_rad: 0.0,
                mean_eps: 23.44,
                true_eps: 23.45,
                ut1_shift: -0.0008,
                ayan: 25.0,
                seen_ex_flags: std::cell::Cell::new(0),
                seen_precess: std::cell::Cell::new(false),
            }
        }
    }

    impl FixstarFlagEnvironment for StubEnv {
        fn precess_ecliptic_to_j2000(
            &self,
            lon_deg: f64,
            lat_deg: f64,
            _jd_tt: f64,
        ) -> Result<(f64, f64), Error> {
            self.seen_precess.set(true);
            Ok(self.precess_to.unwrap_or((lon_deg, lat_deg)))
        }

        fn nutation_dpsi_rad(&self, _jd_tt: f64) -> Result<f64, Error> {
            Ok(self.dpsi_rad)
        }

        fn mean_obliquity_deg(&self, _jd_tt: f64) -> Result<f64, Error> {
            Ok(self.mean_eps)
        }

        fn true_obliquity_deg(&self, _jd_tt: f64) -> Result<f64, Error> {
            Ok(self.true_eps)
        }

        fn ut1_for_tt(&self, jd_tt: f64) -> Result<f64, Error> {
            Ok(jd_tt + self.ut1_shift)
        }

        fn mean_ayanamsha_deg_ut(&self, _jd_ut: f64, ex_flags: i32) -> Result<f64, Error> {
            self.seen_ex_flags.set(ex_flags);
            Ok(self.ayan)
        }
    }

    #[test]
    fn flag_layer_uses_public_spherical_cotrans() {
        use crate::constants::{FLG_EQUATORIAL, FLG_J2000, J2000};
        let env = StubEnv::plain();
        // Independently acquired from the pinned public spherical rotation.
        // Normalizing an intermediate Cartesian vector changes latitude.
        let got = fixstar_positional_transform_via(
            152.1,
            0.5,
            3.0,
            J2000,
            FLG_EQUATORIAL | FLG_J2000,
            true,
            &env,
        )
        .unwrap();
        assert_eq!(
            got.map(f64::to_bits),
            [0x406348ba2e5bfff3, 0x4026634f9dfc0fca, 0x4008000000000000]
        );
        let zero = fixstar_positional_transform_via(
            -0.0,
            0.0,
            3.0,
            J2000,
            FLG_EQUATORIAL | FLG_J2000,
            true,
            &env,
        )
        .unwrap();
        assert_eq!(zero[0].to_bits(), 0);
        let mut env = StubEnv::plain();
        env.true_eps = f64::NAN;
        assert!(matches!(
            fixstar_positional_transform_via(100.0, 20.0, 3.0, J2000, FLG_EQUATORIAL, false, &env,),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn flag_layer_bare_request_is_identity() {
        use crate::constants::J2000;
        // No output flags: the position passes through, speeds are zero,
        // and no provider beyond the trait boundary is consulted.
        let env = StubEnv::plain();
        let raw = [100.0, 20.0, 3.0, 0.5, -0.1, 0.01];
        let got = apply_fixstar_flags_via(raw, J2000, 0, false, &env).expect("bare");
        assert_eq!(got, [100.0, 20.0, 3.0, 0.0, 0.0, 0.0]);
    }

    #[test]
    fn flag_layer_frame_and_equatorial_arms() {
        use crate::constants::J2000;
        use crate::constants::{FLG_EQUATORIAL, FLG_J2000, FLG_NONUT, FLG_SIDEREAL};
        use crate::constants::{FLG_NOABERR, FLG_TRUEPOS};
        // J2000 native skips the backward precession entirely.
        let env = StubEnv::plain();
        let got = fixstar_positional_transform_via(100.0, 20.0, 3.0, J2000, FLG_J2000, true, &env)
            .expect("j2000 native");
        assert_eq!(got, [100.0, 20.0, 3.0]);
        assert!(!env.seen_precess.get());
        // J2000 non-native consults the precession provider.
        let env = StubEnv {
            precess_to: Some((99.5, 20.1)),
            ..StubEnv::plain()
        };
        let got = fixstar_positional_transform_via(100.0, 20.0, 3.0, J2000, FLG_J2000, false, &env)
            .expect("j2000 precess");
        assert_eq!(got, [99.5, 20.1, 3.0]);
        assert!(env.seen_precess.get());
        // NONUT subtracts the stub nutation in longitude.
        let env = StubEnv {
            dpsi_rad: 0.01,
            ..StubEnv::plain()
        };
        let got = fixstar_positional_transform_via(100.0, 20.0, 3.0, J2000, FLG_NONUT, false, &env)
            .expect("nonut");
        assert_eq!(
            got[0].to_bits(),
            (100.0 - 0.01f64.to_degrees()).rem_euclid(360.0).to_bits()
        );
        assert_eq!(got[1].to_bits(), 20.0f64.to_bits());
        // EQUATORIAL without frame bits uses the true-obliquity arm
        // through the shared `cotrans`.
        let env = StubEnv::plain();
        let got =
            fixstar_positional_transform_via(100.0, 20.0, 3.0, J2000, FLG_EQUATORIAL, false, &env)
                .expect("equatorial true");
        let ra = f64::from_bits(0x4059b922d8260765);
        let dec = f64::from_bits(0x4045806b57ca2b88);
        assert_eq!(got, [ra, dec, 3.0]);
        // EQUATORIAL with NONUT uses the mean arm instead (stub dpsi
        // is 0, so the reduction is the identity before the rotation).
        let env = StubEnv::plain();
        let got = fixstar_positional_transform_via(
            100.0,
            20.0,
            3.0,
            J2000,
            FLG_EQUATORIAL | FLG_NONUT,
            false,
            &env,
        )
        .expect("equatorial mean");
        // Pinned public cotrans values, not the fast-calc normalized helper.
        let expect = [
            f64::from_bits(0x4059b900c2471c1e),
            f64::from_bits(0x40457f2bec2f965d),
            3.0,
        ];
        assert_eq!(got, expect);
        // SIDEREAL subtracts the stub mean ayanamsha and carries exactly
        // the request's TRUEPOS/NOABERR/NOGDEFL bits into the ex flags
        // (NONUT always set) — no forced bits are added there.
        let env = StubEnv::plain();
        let got = fixstar_positional_transform_via(
            100.0,
            20.0,
            3.0,
            J2000,
            FLG_SIDEREAL | FLG_TRUEPOS,
            false,
            &env,
        )
        .expect("sidereal");
        // SIDEREAL implies the mean-ecliptic reduction first (stub dpsi
        // is 0 here), then the 25-degree subtraction.
        assert_eq!(
            got[0].to_bits(),
            (100.0f64 - 25.0).rem_euclid(360.0).to_bits()
        );
        assert_eq!(env.seen_ex_flags.get(), FLG_NONUT | FLG_TRUEPOS);
        // An explicit NOABERR bit travels into the ex flags too.
        let env = StubEnv::plain();
        fixstar_positional_transform_via(
            100.0,
            20.0,
            3.0,
            J2000,
            FLG_SIDEREAL | FLG_NOABERR,
            false,
            &env,
        )
        .expect("sidereal noaberr");
        assert_eq!(env.seen_ex_flags.get(), FLG_NONUT | FLG_NOABERR);
        // SIDEREAL equatorial subtracts from right ascension, not longitude.
        let env = StubEnv::plain();
        let eq_only = fixstar_positional_transform_via(
            100.0,
            20.0,
            3.0,
            J2000,
            FLG_EQUATORIAL | FLG_SIDEREAL,
            false,
            &env,
        )
        .expect("sidereal equatorial");
        let (ra, _) = crate::positions::cotrans(
            (100.0f64 - 0.0f64.to_degrees()).rem_euclid(360.0),
            20.0,
            -23.44,
        )
        .expect("cotrans");
        assert_eq!(
            eq_only[0].to_bits(),
            (ra - 25.0).rem_euclid(360.0).to_bits()
        );
    }

    #[test]
    fn flag_layer_speed_differentiates_reported_positions() {
        use crate::constants::J2000;
        use crate::constants::{FLG_NONUT, FLG_SPEED};
        // The stub nutation drifts linearly with epoch; the raw speeds
        // are zero, so the reported speed must equal exactly the drift
        // rate — a frozen-frame rotation would report 0.
        struct DriftEnv {
            rate_deg_day: f64,
        }
        impl FixstarFlagEnvironment for DriftEnv {
            fn precess_ecliptic_to_j2000(
                &self,
                lon_deg: f64,
                lat_deg: f64,
                _jd_tt: f64,
            ) -> Result<(f64, f64), Error> {
                Ok((lon_deg, lat_deg))
            }

            fn nutation_dpsi_rad(&self, jd_tt: f64) -> Result<f64, Error> {
                Ok((self.rate_deg_day * (jd_tt - crate::constants::J2000)).to_radians())
            }

            fn mean_obliquity_deg(&self, _jd_tt: f64) -> Result<f64, Error> {
                Ok(23.44)
            }

            fn true_obliquity_deg(&self, _jd_tt: f64) -> Result<f64, Error> {
                Ok(23.44)
            }

            fn ut1_for_tt(&self, jd_tt: f64) -> Result<f64, Error> {
                Ok(jd_tt)
            }

            fn mean_ayanamsha_deg_ut(&self, _jd_ut: f64, _ex_flags: i32) -> Result<f64, Error> {
                Ok(0.0)
            }
        }
        let env = DriftEnv { rate_deg_day: 0.01 };
        let raw = [100.0, 20.0, 3.0, 0.0, 0.0, 0.0];
        let got =
            apply_fixstar_flags_via(raw, J2000, FLG_SPEED | FLG_NONUT, false, &env).expect("speed");
        // Reported position at the request epoch carries the drift offset
        // (zero there by construction); the speed carries minus the drift
        // rate (the layer *subtracts* the growing nutation).
        assert_eq!(got[0].to_bits(), 100.0f64.to_bits());
        assert!((got[3] + 0.01).abs() < 1e-12);
        assert_eq!(got[4].to_bits(), 0.0f64.to_bits());
        assert_eq!(got[5].to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn flag_layer_xyz_and_radians_shaping() {
        use crate::constants::J2000;
        use crate::constants::{FLG_RADIANS, FLG_SPEED, FLG_XYZ};
        let env = StubEnv::plain();
        let raw = [100.0, 20.0, 3.0, 0.02, -0.01, 0.001];
        let xyz =
            apply_fixstar_flags_via(raw, J2000, FLG_SPEED | FLG_XYZ, false, &env).expect("xyz");
        // The flag-layer speed reconstruction rounds through the
        // midpoint ± speed·h association (like the model), so the
        // comparison is tight-absolute, not exactly.
        let expect = crate::positions::spherical_to_cartesian_with_velocity(
            100.0, 20.0, 3.0, 0.02, -0.01, 0.001,
        );
        for (a, b) in xyz.iter().zip(expect.iter()) {
            assert!((a - b).abs() < 1e-12, "{a} vs {b}");
        }
        // XYZ wins over RADIANS when both are set.
        // (No FLG_SPEED here, so the shaped speeds are zero.)
        let both =
            apply_fixstar_flags_via(raw, J2000, FLG_XYZ | FLG_RADIANS, false, &env).expect("both");
        let expect_nospeed =
            crate::positions::spherical_to_cartesian_with_velocity(100.0, 20.0, 3.0, 0.0, 0.0, 0.0);
        assert_eq!(both, expect_nospeed);
        let rad = apply_fixstar_flags_via(raw, J2000, FLG_RADIANS, false, &env).expect("radians");
        assert_eq!(
            rad,
            [
                100.0f64.to_radians(),
                20.0f64.to_radians(),
                3.0,
                0.0,
                0.0,
                3.0 - 3.0,
            ]
        );
        assert_eq!(rad[5].to_bits(), 0.0f64.to_bits());
    }

    #[test]
    fn angular_flag_transforms_preserve_tiny_stellar_radial_rates() {
        use crate::constants::{
            FLG_EQUATORIAL, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_SPEED, J2000,
        };
        let env = StubEnv::plain();
        for range in [1e7, 1e8, 1e14] {
            for rate in [
                -0.002_021_395_546_499_799_6,
                0.000_577_550_253_864_272,
                1e-12,
            ] {
                for flags in [
                    0,
                    FLG_EQUATORIAL,
                    FLG_J2000,
                    FLG_NONUT,
                    FLG_SIDEREAL,
                    FLG_RADIANS,
                ] {
                    let raw = [100.0, 20.0, range, 0.02, -0.01, rate];
                    let reported =
                        apply_fixstar_flags_via(raw, J2000, flags | FLG_SPEED, false, &env)
                            .unwrap();
                    assert_eq!(reported[2], range);
                    assert_eq!(reported[5], rate);
                }
            }
        }
    }

    #[test]
    fn topo_site_constructor_and_missing_site_refusal() {
        // Checked constructor pins the components and refuses
        // non-finite input.
        let site = StarTopoSite::new(10.0, 47.0, 500.0).expect("site");
        assert_eq!(
            (site.lon_deg(), site.lat_deg(), site.height_m()),
            (10.0, 47.0, 500.0)
        );
        assert!(StarTopoSite::new(f64::NAN, 0.0, 0.0).is_err());
        // TOPOCTR without a site refuses instead of silently
        // falling back to geocentric (needs a provider only to reach
        // the observer branch — the refusal fires first).
        struct NoProvider;
        impl crate::positions::StateProvider for NoProvider {
            fn provider_name(&self) -> &'static str {
                "test-stub"
            }

            fn has_body(&self, _body_id: i32) -> bool {
                false
            }

            fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
                None
            }

            fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
                None
            }

            fn eval_raw(&self, _body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
                Err(Error::LebBodyNotFound { body_id: 0 })
            }

            fn eval_raw_two_part(
                &self,
                _body_id: i32,
                _epoch: crate::types::TwoPartEpoch,
            ) -> Result<(Vec<f64>, Vec<f64>), Error> {
                Err(Error::LebBodyNotFound { body_id: 0 })
            }

            fn serving_entry_snapshot(
                &self,
                body_id: i32,
            ) -> Result<crate::positions::EntrySnapshot, Error> {
                Err(Error::LebBodyNotFound { body_id })
            }
        }
        let nutation_at = |_jd: f64| -> Result<(f64, f64), Error> { Ok((0.0, 0.0)) };
        let deltat_at = |_jd: f64| -> Result<f64, Error> { Ok(0.0008) };
        let serving = SealedStarServing::new(
            &NoProvider,
            &nutation_at,
            &deltat_at,
            crate::sidereal::SiderealSelection::unset(),
        );
        let regulus = resolve_star_ref("Regulus").expect("Regulus");
        match serving.fixstar_tt(
            &regulus,
            crate::constants::J2000,
            crate::constants::FLG_TOPOCTR,
            None,
        ) {
            Err(Error::ConfigurationMissing { .. }) => {}
            other => panic!("expected ConfigurationMissing, got {other:?}"),
        }
    }

    // Vega under J 1900 mode 19 at 2460463.0 TT from Greenwich. The pinned
    // reference uses the same hash-checked base coefficients and site.
    #[test]
    fn topocentric_vega_matches_source_within_precision_contract() {
        use crate::positions::FRAME_LIGHT_SPEED_AU_DAY;
        assert_eq!(
            FRAME_LIGHT_SPEED_AU_DAY.to_bits(),
            0x4065a4a0d4b3b973,
            "light speed in AU per day"
        );
        let Some(bytes) = anchor_test_fixture() else {
            eprintln!(
                "topocentric_vega_matches_source_within_precision_contract: SKIP (missing LEB2 core fixture)"
            );
            return;
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let provider = crate::positions::SingleLeb::new(&index, &bytes);
        let nutation_at = |jd: f64| index.eval_nutation(&bytes, jd);
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
        let mut selection = crate::sidereal::SiderealSelection::unset();
        let _warnings = selection.set_sidereal_mode(19, 0.0, 0.0);
        let serving = SealedStarServing::new(&provider, &nutation_at, &deltat_at, selection)
            .with_default_timescale();
        // Time conversion keeps the one-microsecond PLAN limit.
        for (jd, word) in [
            (2460463.0, 0x4142c59780000016),
            (2451545.0, 0x4142b42c7ffffffe),
            (2415020.0, 0x41426cd5ffffffff),
            (2500000.0, 0x414312cfffffffdf),
        ] {
            assert!(
                (crate::time::tt_to_tdb_jd(jd).expect("tdb") - f64::from_bits(word)).abs()
                    * 86400.0
                    <= 1e-6,
                "TDB conversion at {jd}"
            );
        }
        let vega = resolve_star_ref("Vega").expect("Vega");
        let site = StarTopoSite::new(0.0, 51.4769, 46.0).expect("site");
        let (place, name, ret) = serving
            .fixstar_tt_wrapped(&vega, 2460463.0, 98304, Some(site))
            .expect("topocentric Vega");
        assert_eq!(name, "Vega,alLyr");
        assert_eq!(ret, 98304);
        for (slot, expected) in [0x4071bf0f9c83c704, 0x404edf3e2c088a77, 0x41382a7faaea6e0e]
            .into_iter()
            .enumerate()
        {
            let expected = f64::from_bits(expected);
            let tolerance = if slot < 2 {
                0.001 / 3600.0
            } else {
                expected.abs() * 1e-8
            };
            assert!(
                (place[slot] - expected).abs() <= tolerance,
                "Vega slot {slot}: {} vs {expected}",
                place[slot]
            );
        }
        assert_eq!(place[3..], [0.0; 3], "rates are zero without a speed flag");
    }

    // Regression cases for observer-state acquisition, coordinate frames and stellar rates.
    #[test]
    fn topocentric_rotation_preserves_the_geodetic_site_radius() {
        let Some(bytes) = anchor_test_fixture() else {
            eprintln!(
                "topocentric_rotation_preserves_the_geodetic_site_radius: SKIP (missing LEB2 core fixture)"
            );
            return;
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let provider = crate::positions::SingleLeb::new(&index, &bytes);
        let nutation_at = |jd: f64| index.eval_nutation(&bytes, jd);
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
        let selection = crate::sidereal::SiderealSelection::unset();
        let serving = SealedStarServing::new(&provider, &nutation_at, &deltat_at, selection)
            .with_default_timescale();
        let site = StarTopoSite::new(0.0, 51.4769, 46.0).expect("site");
        let obs = serving
            .observer_at(crate::constants::FLG_TOPOCTR, Some(site), 2460463.0)
            .expect("observer");
        let site = obs.site_gcrs_au().expect("site vector");
        let norm = |vector: [f64; 3]| vector.iter().map(|x| x * x).sum::<f64>().sqrt();
        let terrestrial = crate::positions::geodetic_to_itrf_m_in(
            0.0,
            51.4769_f64.to_radians(),
            46.0,
            crate::positions::SiteEllipsoid::Wgs84,
        )
        .unwrap();
        assert!((norm(site) * crate::positions::AU_IN_M - norm(terrestrial)).abs() < 1e-6);
        assert!(obs.pos_au().iter().all(|v| v.is_finite()));
        assert!(obs.vel_au_day().iter().all(|v| v.is_finite()));
    }

    // Regulus with the default Fagan/Bradley selection at a common nominal
    // UT/TT epoch, using the hash-checked base coefficients. Distinct clock
    // conversions preserve the separate pinned reference inputs.
    #[test]
    fn fagan_default_ayanamsha_matches_source_within_precision_contract() {
        let Some(bytes) = anchor_test_fixture() else {
            eprintln!(
                "fagan_default_ayanamsha_matches_source_within_precision_contract: SKIP (missing LEB2 core fixture)"
            );
            return;
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let provider = crate::positions::SingleLeb::new(&index, &bytes);
        let nutation_at = |jd: f64| index.eval_nutation(&bytes, jd);
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        let deltat_at = |jd_ut: f64| crate::time::deltat_with_config(jd_ut, &cfg);
        let selection = crate::sidereal::SiderealSelection::unset();
        let serving = SealedStarServing::new(&provider, &nutation_at, &deltat_at, selection)
            .with_default_timescale();
        let regulus = resolve_star_ref("Regulus").expect("Regulus");
        let assert_reference = |slot: usize, actual: f64, word: u64| {
            let expected = f64::from_bits(word);
            let tolerance = match slot {
                0 | 1 => 0.001 / 3600.0,
                2 => expected.abs() * 1e-8,
                _ => 0.0, // No speed flag: the physical result has zero rates.
            };
            assert!(
                (actual - expected).abs() <= tolerance,
                "Regulus slot {slot}: {actual} vs {expected} (limit {tolerance})"
            );
        };
        // UT entry.
        let (place, name, ret) = serving
            .fixstar_ut_wrapped(&regulus, 2460463.0, 65536, None)
            .expect("Fagan-default Regulus UT");
        assert_eq!(name, "Regulus,alLeo");
        assert_eq!(ret, 65602);
        for (i, (g, w)) in place
            .iter()
            .zip([
                0x405f458288d376b8u64,
                0x3fddcdbc6cee295eu64,
                0x415321709dbbd452u64,
                0x0000000000000000u64,
                0x0000000000000000u64,
                0x0000000000000000u64,
            ])
            .enumerate()
        {
            assert_reference(i, *g, w);
        }
        // TT entry at the same nominal epoch.
        let (place, name, ret) = serving
            .fixstar_tt_wrapped(&regulus, 2460463.0, 65536, None)
            .expect("Fagan-default Regulus TT");
        assert_eq!(name, "Regulus,alLeo");
        assert_eq!(ret, 65536);
        for (i, (g, w)) in place
            .iter()
            .zip([
                0x405f458289242a18u64,
                0x3fddcdbc6cfd7d4au64,
                0x415321709dbb90dcu64,
                0x0000000000000000u64,
                0x0000000000000000u64,
                0x0000000000000000u64,
            ])
            .enumerate()
        {
            assert_reference(i, *g, w);
        }
        // The batch entry has the same reference UT longitude.
        let batch = serving
            .batch_fixstars_ut(&["Regulus"], 2460463.0, 65536, None, false)
            .expect("batch")[0]
            .clone()
            .expect("slot");
        assert_reference(0, batch.0[0], 0x405f458288d376b8u64);
    }

    #[test]
    fn catalog_initial_states_are_finite_and_nonzero() {
        for (name, id) in [("Vega", 1000007), ("Regulus", 1000001), ("Sirius", 1000004)] {
            let data = fixed_star_data(id).expect("catalog star is present");
            let (position, velocity) =
                star_icrs_state_catalog_angles(&data, crate::constants::J2000);
            // A catalog direction with a finite positive distance defines a
            // nonzero position. This physical invariant needs no captured
            // implementation words or fitted numerical tolerance.
            assert!(
                position.iter().all(|value| value.is_finite()),
                "{name} position"
            );
            assert!(
                velocity.iter().all(|value| value.is_finite()),
                "{name} velocity"
            );
            assert!(
                position.iter().any(|value| *value != 0.0),
                "{name} distance"
            );
        }
    }

    // Observer-state acquisition, coordinate frames and stellar rates.

    #[test]
    fn bundled_star_timescale_is_independent_of_model_delta_t() {
        struct NoProvider;
        impl crate::positions::StateProvider for NoProvider {
            fn provider_name(&self) -> &'static str {
                "clock-only"
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
            fn eval_raw(&self, body_id: i32, _: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
                Err(Error::LebBodyNotFound { body_id })
            }
            fn eval_raw_two_part(
                &self,
                body_id: i32,
                _: crate::types::TwoPartEpoch,
            ) -> Result<(Vec<f64>, Vec<f64>), Error> {
                Err(Error::LebBodyNotFound { body_id })
            }
            fn serving_entry_snapshot(
                &self,
                body_id: i32,
            ) -> Result<crate::positions::EntrySnapshot, Error> {
                Err(Error::LebBodyNotFound { body_id })
            }
        }
        let calls = std::cell::Cell::new(0);
        let model = |_: f64| {
            calls.set(calls.get() + 1);
            Ok(0.001)
        };
        let nutation = |_: f64| Ok((0.0, 0.0));
        let serving = SealedStarServing::new(
            &NoProvider,
            &nutation,
            &model,
            crate::sidereal::SiderealSelection::unset(),
        );
        let ut = 2460000.5;
        let legacy = serving.jd_tt_for_ut(ut).unwrap();
        assert_eq!(legacy.to_bits(), (ut + 0.001).to_bits());
        assert_eq!(calls.get(), 1);
        let serving = serving.with_default_timescale();
        let tt = serving.jd_tt_for_ut(ut).unwrap();
        assert_eq!(
            tt.to_bits(),
            crate::time::default_timescale_tt_from_ut1(ut)
                .unwrap()
                .to_bits()
        );
        assert_ne!(tt.to_bits(), legacy.to_bits());
        assert_eq!(
            serving.ut1_for_tt(tt).unwrap().to_bits(),
            crate::time::default_timescale_ut1_from_tt(tt)
                .unwrap()
                .to_bits()
        );
        assert_eq!(calls.get(), 1, "timescale must not consume model Delta-T");
        assert_eq!((serving.deltat_days_at)(ut).unwrap(), 0.001);
        assert!(serving.jd_tt_for_ut(f64::NAN).is_err());
        assert!(serving.ut1_for_tt(f64::INFINITY).is_err());

        // Classical SIDBIT zero points must use the selected independent
        // clock as well, not a second UT->TT leg through the model callback.
        // Barycentric carriers need no fixture/provider and isolate that leg.
        use crate::constants::{FLG_BARYCTR, FLG_EQUATORIAL, FLG_SIDEREAL, FLG_XYZ};
        let model_a = |_: f64| Ok(0.001);
        let model_b = |_: f64| Ok(0.002);
        let star = resolve_star_ref("Spica").unwrap();
        for mode in [257, 769, 256, 513] {
            let mut selection = crate::sidereal::SiderealSelection::unset();
            selection.set_sidereal_mode(mode, 0.0, 0.0);
            let a = SealedStarServing::new(&NoProvider, &nutation, &model_a, selection)
                .with_default_timescale();
            let b = SealedStarServing::new(&NoProvider, &nutation, &model_b, selection)
                .with_default_timescale();
            for epoch in [2451545.0, 2460463.0] {
                for representation in [FLG_XYZ, FLG_XYZ | FLG_EQUATORIAL] {
                    let flags = FLG_SIDEREAL | FLG_BARYCTR | representation;
                    let (xx_a, name_a, rf_a) =
                        a.fixstar_tt_wrapped(&star, epoch, flags, None).unwrap();
                    let (xx_b, name_b, rf_b) =
                        b.fixstar_tt_wrapped(&star, epoch, flags, None).unwrap();
                    assert_eq!(
                        xx_a.map(f64::to_bits),
                        xx_b.map(f64::to_bits),
                        "mode={mode}, epoch={epoch}, flags={flags}"
                    );
                    assert_eq!((name_a, rf_a), (name_b, rf_b));
                }
            }
        }
    }

    #[test]
    fn wrapped_entries_match_plain_without_sidereal_request() {
        // Without `FLG_SIDEREAL` neither prologue fires: the wrapped
        // entries are the plain entries (same value or same failure).
        // `NoProvider` fails every computation identically, so this pins
        // the dispatch wiring, not the ephemeris.
        struct NoProvider;
        impl crate::positions::StateProvider for NoProvider {
            fn provider_name(&self) -> &'static str {
                "test-stub"
            }

            fn has_body(&self, _body_id: i32) -> bool {
                false
            }

            fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
                None
            }

            fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
                None
            }

            fn eval_raw(&self, _body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
                Err(Error::LebBodyNotFound { body_id: 0 })
            }

            fn eval_raw_two_part(
                &self,
                _body_id: i32,
                _epoch: crate::types::TwoPartEpoch,
            ) -> Result<(Vec<f64>, Vec<f64>), Error> {
                Err(Error::LebBodyNotFound { body_id: 0 })
            }

            fn serving_entry_snapshot(
                &self,
                body_id: i32,
            ) -> Result<crate::positions::EntrySnapshot, Error> {
                Err(Error::LebBodyNotFound { body_id })
            }
        }
        let nutation_at = |_jd: f64| -> Result<(f64, f64), Error> { Ok((0.0, 0.0)) };
        let deltat_at = |_jd: f64| -> Result<f64, Error> { Ok(0.0008) };
        let serving = SealedStarServing::new(
            &NoProvider,
            &nutation_at,
            &deltat_at,
            crate::sidereal::SiderealSelection::unset(),
        );
        let regulus = resolve_star_ref("Regulus").expect("Regulus");
        let jd = crate::constants::J2000;
        // Plain tropical and sidereal-under-default-mode requests agree
        // between wrapped and plain spellings (same provider failure).
        for flags in [
            0,
            crate::constants::FLG_SPEED,
            crate::constants::FLG_SIDEREAL,
            crate::constants::FLG_SIDEREAL | crate::constants::FLG_EQUATORIAL,
        ] {
            let plain_ut = serving.fixstar_ut(&regulus, jd, flags, None);
            let wrapped_ut = serving.fixstar_ut_wrapped(&regulus, jd, flags, None);
            assert_eq!(
                format!("{plain_ut:?}"),
                format!("{wrapped_ut:?}"),
                "UT flags {flags}"
            );
            let plain_tt = serving.fixstar_tt(&regulus, jd, flags, None);
            let wrapped_tt = serving.fixstar_tt_wrapped(&regulus, jd, flags, None);
            assert_eq!(
                format!("{plain_tt:?}"),
                format!("{wrapped_tt:?}"),
                "TT flags {flags}"
            );
        }
        // The v 2 spellings share the same wrapped core.
        let regulus2 = resolve_star2("Regulus").expect("Regulus v2");
        let w1 = serving.fixstar2_ut_wrapped(&regulus2, jd, 0, None);
        let w2 = serving.fixstar_ut_wrapped(&regulus, jd, 0, None);
        assert_eq!(format!("{w1:?}"), format!("{w2:?}"));
    }

    #[test]
    fn batch_slots_skip_errors_and_fail_fast_identity() {
        // Batch resolution/slot/fail-fast semantics over a failing
        // provider: misses and computation failures become `None` slots
        // under `skip_errors`, input order is preserved, and without it
        // the first element's error aborts the batch.
        struct NoProvider;
        impl crate::positions::StateProvider for NoProvider {
            fn provider_name(&self) -> &'static str {
                "test-stub"
            }

            fn has_body(&self, _body_id: i32) -> bool {
                false
            }

            fn body_coverage(&self, _body_id: i32) -> Option<(f64, f64)> {
                None
            }

            fn serving_coord_type(&self, _body_id: i32) -> Option<u32> {
                None
            }

            fn eval_raw(&self, _body_id: i32, _jd: f64) -> Result<(Vec<f64>, Vec<f64>), Error> {
                Err(Error::LebBodyNotFound { body_id: 0 })
            }

            fn eval_raw_two_part(
                &self,
                _body_id: i32,
                _epoch: crate::types::TwoPartEpoch,
            ) -> Result<(Vec<f64>, Vec<f64>), Error> {
                Err(Error::LebBodyNotFound { body_id: 0 })
            }

            fn serving_entry_snapshot(
                &self,
                body_id: i32,
            ) -> Result<crate::positions::EntrySnapshot, Error> {
                Err(Error::LebBodyNotFound { body_id })
            }
        }
        let nutation_at = |_jd: f64| -> Result<(f64, f64), Error> { Ok((0.0, 0.0)) };
        let deltat_at = |_jd: f64| -> Result<f64, Error> { Ok(0.0008) };
        let serving = SealedStarServing::new(
            &NoProvider,
            &nutation_at,
            &deltat_at,
            crate::sidereal::SiderealSelection::unset(),
        );
        let jd = crate::constants::J2000;
        // skip_errors: every slot is None (miss + computation failure
        // alike), order and length preserved.
        let slots = serving
            .batch_fixstars_ut(&["NoSuchStar_xyz", "Regulus", ""], jd, 0, None, true)
            .expect("skip_errors batch");
        assert_eq!(slots.len(), 3);
        assert!(slots.iter().all(|s| s.is_none()));
        // Empty input is an empty success, never an error.
        let empty = serving
            .batch_fixstars_ut(&[], jd, 0, None, true)
            .expect("empty batch");
        assert!(empty.is_empty());
        // Fail-fast: the resolution miss aborts before the later
        // computation failure is reached.
        match serving.batch_fixstars_ut(&["NoSuchStar_xyz", "Regulus"], jd, 0, None, false) {
            Err(Error::StarNotFound { .. }) => {}
            other => panic!("expected StarNotFound, got {other:?}"),
        }
        // Fail-fast: a resolvable first element surfaces the
        // computation failure instead.
        match serving.batch_fixstars_ut(&["Regulus", "NoSuchStar_xyz"], jd, 0, None, false) {
            Err(Error::LebBodyNotFound { .. }) => {}
            other => panic!("expected LebBodyNotFound, got {other:?}"),
        }
        // Nondecimal numeric characters follow ordinary missing-name policy.
        match serving.batch_fixstars_ut(&["²"], jd, 0, None, false) {
            Err(Error::StarNotFound {
                search_type: Some("name"),
                ..
            }) => {}
            other => panic!("expected StarNotFound/name, got {other:?}"),
        }
        let skipped = serving
            .batch_fixstars_ut(&["Regulus", "²"], jd, 0, None, true)
            .expect("both computation and missing-name failures can be skipped");
        assert_eq!(skipped.len(), 2);
        assert!(skipped.iter().all(Option::is_none));
        // A non-finite epoch refuses before any element runs.
        match serving.batch_fixstars_ut(&["Regulus"], f64::NAN, 0, None, true) {
            Err(Error::NonFiniteArgument { .. }) => {}
            other => panic!("expected NonFiniteArgument, got {other:?}"),
        }
    }

    #[test]
    fn photometry_qualifies_only_reviewed_johnson_v_entries() {
        for (name, hip, expected, historical) in [
            ("Spica", 65474, 0.98_f64, Some(1.04_f64)),
            ("Antares", 80763, 1.06_f64, None),
        ] {
            let id = resolve_star_name(name).expect("catalogued name");
            let entry = star_entry_by_id(id).expect("catalogued ID");
            let value = star_photometry(id).expect("photometry for every catalogued ID");
            assert_eq!(entry.hip_number(), hip);
            assert_eq!(value, entry.photometry());
            assert_eq!(value.value.to_bits(), expected.to_bits());
            assert_eq!(entry.magnitude().to_bits(), expected.to_bits());
            assert_eq!(value.band, PhotometricBand::JohnsonV);
            assert_eq!(
                value.catalogue_field,
                PhotometricCatalogueField::HipparcosMainVmag
            );
            assert_eq!(value.statistic, PhotometricStatistic::Unspecified);
            assert_eq!(value.epoch, None);
            assert_eq!(value.variability, PhotometricVariability::Unknown);
            assert_eq!(
                value.historical_unqualified_value.map(f64::to_bits),
                historical.map(f64::to_bits)
            );
        }
        for entry in list_fixed_stars()
            .iter()
            .filter(|entry| !matches!(entry.hip_number(), 65474 | 80763))
        {
            let value = entry.photometry();
            assert_eq!(value.value.to_bits(), entry.magnitude().to_bits());
            assert_eq!(value.band, PhotometricBand::Unspecified);
            assert_eq!(
                value.catalogue_field,
                PhotometricCatalogueField::Unqualified
            );
            assert_eq!(value.statistic, PhotometricStatistic::Unspecified);
            assert_eq!(value.epoch, None);
            assert_eq!(value.variability, PhotometricVariability::Unknown);
            assert_eq!(value.historical_unqualified_value, None);
        }
        assert_ne!(
            PhotometricVariability::Unknown,
            PhotometricVariability::Constant
        );
        assert!(star_photometry(-1).is_none());
        assert!(star_photometry(crate::constants::FIXSTAR_OFFSET).is_none());
    }

    #[test]
    fn copied_ids_do_not_qualify_modified_photometric_rows() {
        let id = resolve_star_name("Spica").expect("catalogued name");
        let canonical = star_entry_by_id(id).expect("catalogued ID");
        let copied = *canonical;
        assert_eq!(copied.photometry(), canonical.photometry());
        let different_magnitude = StarCatalogEntry::from_row(
            canonical.id,
            canonical.name,
            canonical.nomenclature,
            canonical.hip_number,
            canonical.data.ra_j2000_deg,
            canonical.data.dec_j2000_deg,
            canonical.data.pm_ra_arcsec_yr,
            canonical.data.pm_dec_arcsec_yr,
            canonical.data.parallax_mas,
            canonical.data.radial_velocity_km_s,
            7.0,
        );
        let mut different_astrometry = copied;
        different_astrometry.data.ra_j2000_deg =
            f64::from_bits(canonical.data.ra_j2000_deg.to_bits() ^ 1);
        let mut different_identity = copied;
        different_identity.hip_number = 80763;
        for entry in [
            different_magnitude,
            different_astrometry,
            different_identity,
        ] {
            let value = entry.photometry();
            assert_eq!(value.value.to_bits(), entry.magnitude().to_bits());
            assert_eq!(value, StarPhotometry::unqualified(entry.magnitude()));
        }
    }

    #[test]
    fn qualified_photometry_preserves_magnitude_aliases_and_missing_name_errors() {
        for (name, designation, expected) in [
            ("Spica", ",alVir", 0.98_f64),
            ("Antares", ",alSco", 1.06_f64),
        ] {
            let id = resolve_star_name(name).expect("catalogued name");
            let photometry = star_photometry(id).expect("catalogued photometry");
            for spelling in [name, designation] {
                for result in [fixstar_mag(spelling), fixstar2_mag(spelling)] {
                    let (value, canonical) = result.expect("name and designation");
                    assert_eq!(value.to_bits(), expected.to_bits());
                    assert_eq!(value.to_bits(), photometry.value.to_bits());
                    assert!(canonical.starts_with(name));
                }
            }
        }
        for result in [
            fixstar_mag("missing-photometry"),
            fixstar2_mag("missing-photometry"),
        ] {
            assert!(matches!(result, Err(Error::StarNotFound { .. })));
        }
    }

    #[test]
    fn fixstar_magnitude_entries() {
        // Regulus pins: the model docstring example (1.40) and the
        // "Name,nomenclature" spelling, exactly on the stored value.
        let (mag, name) = fixstar_mag("Regulus").expect("Regulus v1");
        assert_eq!(mag.to_bits(), 1.4f64.to_bits());
        assert_eq!(name, "Regulus,alLeo");
        let (mag2, name2) = fixstar2_mag("Regulus").expect("Regulus v2");
        assert_eq!(mag2.to_bits(), mag.to_bits());
        assert_eq!(name2, name);
        // Nomenclature forms key the other part on each side (v 1 ignores
        // the name, v 2 ignores everything but the nomenclature).
        let (mag_n, name_n) = fixstar_mag(",alLeo").expect(",alLeo v1");
        assert_eq!(
            (mag_n.to_bits(), name_n.as_str()),
            (mag.to_bits(), name.as_str())
        );
        let (mag_n2, name_n2) = fixstar2_mag("Regulus,alLeo").expect("v2 nomenclature");
        assert_eq!(
            (mag_n2.to_bits(), name_n2.as_str()),
            (mag.to_bits(), name.as_str())
        );
        // The v 1 implicit prefix and the v 2 trailing-% wildcard agree.
        let (mag_p, _) = fixstar_mag("Regu").expect("v1 prefix");
        assert_eq!(mag_p.to_bits(), mag.to_bits());
        let (mag_w, _) = fixstar2_mag("Regu%").expect("v2 wildcard");
        assert_eq!(mag_w.to_bits(), mag.to_bits());
        // A '%' is an ordinary name character in the v 1 family.
        match fixstar_mag("Regu%") {
            Err(Error::StarNotFound { .. }) => {}
            other => panic!("expected StarNotFound, got {other:?}"),
        }
        // Sequential numbers index the family sort orders; both spellings
        // of star #12 agree with the resolver triple.
        let v1seq = resolve_star_ref("12").expect("v1 seq 12");
        let (mag_s, name_s) = fixstar_mag("12").expect("v1 seq mag");
        assert_eq!(
            mag_s.to_bits(),
            star_magnitude(v1seq.star_id).expect("v1 seq row").to_bits()
        );
        assert_eq!(name_s, v1seq.canonical);
        let v2seq = resolve_star2("12").expect("v2 seq 12");
        let (mag_s2, name_s2) = fixstar2_mag("12").expect("v2 seq mag");
        assert_eq!(
            mag_s2.to_bits(),
            star_magnitude(v2seq.star_id).expect("v2 seq row").to_bits()
        );
        assert_eq!(name_s2, v2seq.canonical);
        // Unknown keys refuse with the model search form on both sides.
        match fixstar_mag("NoSuchStar_xyz") {
            Err(Error::StarNotFound {
                search_type: Some("name"),
                ..
            }) => {}
            other => panic!("expected StarNotFound/name, got {other:?}"),
        }
        match fixstar2_mag("") {
            Err(Error::StarNotFound {
                search_type: Some("empty"),
                ..
            }) => {}
            other => panic!("expected StarNotFound/empty, got {other:?}"),
        }
        match fixstar2_mag("NoSuchStar_xyz%") {
            Err(Error::StarNotFound {
                search_type: Some("wildcard"),
                ..
            }) => {}
            other => panic!("expected StarNotFound/wildcard, got {other:?}"),
        }
        // Nondecimal numeric characters are ordinary unresolved catalog names.
        match fixstar_mag("²") {
            Err(Error::StarNotFound {
                search_type: Some("name"),
                ..
            }) => {}
            other => panic!("expected StarNotFound/name, got {other:?}"),
        }
        match fixstar2_mag("²") {
            Err(Error::StarNotFound {
                search_type: Some("name"),
                ..
            }) => {}
            other => panic!("expected StarNotFound/name, got {other:?}"),
        }
    }

    #[test]
    fn every_catalog_row_carries_a_magnitude() {
        // A missing magnitude is reported before the `StarMagnitudeUnavailable`
        // arm is unreachable over the swept resolvers because every
        // catalog id has a magnitude row.
        for entry in list_fixed_stars() {
            assert!(
                star_magnitude(entry.id()).is_some(),
                "no magnitude for {}",
                entry.id()
            );
        }
        assert_eq!(list_fixed_stars().len(), 1447);
    }

    #[test]
    fn sirius_place_matches_source_within_precision_at_shared_tt() {
        // Sirius under the same TT epoch and immutable coefficient image as
        // the pinned reference. The separate model-clock check also retains
        // the distinction between explicit Delta-T and the bundled UT clock.
        use crate::constants::EARTH;
        use crate::positions::StateProvider;
        let Some(bytes) = anchor_test_fixture() else {
            eprintln!(
                "sirius_place_matches_source_within_precision_at_shared_tt: SKIP (missing released core fixture)"
            );
            return;
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let composite = crate::leb::CompositeIndex::new(vec![(index.clone(), &bytes[..])])
            .expect("one reader must build");
        let cfg = crate::time::DeltaTConfig {
            model: crate::time::DeltaTModel::EspenakMeeus,
            user_override_days: None,
            tidal_acceleration: crate::constants::TIDAL_DEFAULT,
            iers_enabled: false,
            iers_lookup_seconds: None,
            smh2016_lookup_seconds: None,
        };
        // The native chain converts 2460000.5 UT through the
        // Espenak-Meeus polynomial (never the observed table).
        let jd_ut = 2460000.5;
        let jd_tt = jd_ut + crate::time::deltat_with_config(jd_ut, &cfg).expect("dt serves");
        assert_eq!(jd_tt.to_bits(), 0x4142c4b0401bd273);
        // Bundled-clock TT for the same UT epoch, independently of the
        // Espenak-Meeus polynomial used above.
        let jd_rel = f64::from_bits(0x4142c4b0401a3e92);
        assert!(
            (jd_tt - jd_rel).abs() < 1e-4,
            "native/observed TT gap exceeds the characterized envelope"
        );
        let nutation_at = |jd: f64| composite.eval_nutation(jd);
        let data = crate::stars::fixed_star_data(1000004).expect("sirius row");
        let earth = composite
            .geometric_state(EARTH, jd_rel)
            .expect("earth serves");
        let (dpsi, deps) = nutation_at(jd_rel).expect("nut serves");
        let frames =
            crate::positions::frame_matrices_for_date(jd_rel, dpsi, deps).expect("frames serve");
        let observer = crate::stars::StarObserver::new(earth.pos(), earth.vel()).expect("observer");
        let frame =
            crate::stars::StarFrame::new(frames.pn(), frames.eps_true_rad()).expect("frame");
        let triple = crate::stars::star_apparent_ecliptic_via(
            &data,
            jd_rel,
            &observer,
            false,
            false,
            false,
            |id, jd| Ok(composite.geometric_state(id, jd)?.pos()),
            Some(frame),
        )
        .expect("apparent serves");
        // Original reference values, tested by their physical observables.
        for (slot, word) in [0x405a19c6bba3b1a4, 0xc043ce85c46ce2c6, 0x4120994335085ce4]
            .into_iter()
            .enumerate()
        {
            let reference = f64::from_bits(word);
            let tolerance = if slot < 2 {
                0.001 / 3600.0
            } else {
                reference.abs() * 1e-8
            };
            assert!(
                (triple[slot] - reference).abs() <= tolerance,
                "Sirius slot {slot}: {} vs {reference}",
                triple[slot]
            );
        }
    }

    /// Released LEB 2 core for the anchor-provider gate: the
    /// `LIBEPHEMERIS_TEST_LEB2` override or the hash-verified
    /// current core in the configured test data directory. Admitted hashes
    /// are checked by `crate::test_data::admitted_core_hash`.
    /// Absence skips; a present file with the wrong hash fails (never
    /// silent evaluation on unrelated coefficients).
    fn anchor_test_fixture() -> Option<Vec<u8>> {
        use sha2::{Digest, Sha256};
        let path = std::env::var_os("LIBEPHEMERIS_TEST_LEB2")
            .map(std::path::PathBuf::from)
            .unwrap_or_else(crate::test_data::core_leb_path);
        let bytes = match std::fs::read(&path) {
            Ok(bytes) => bytes,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return None,
            Err(error) => panic!("cannot read fixture {}: {error}", path.display()),
        };
        assert!(
            crate::test_data::admitted_core_hash(
                &Sha256::digest(&bytes)
                    .iter()
                    .map(|byte| format!("{byte:02x}"))
                    .collect::<String>()
            ),
            "fixture {} is not the released base_core",
            path.display()
        );
        Some(bytes)
    }

    #[test]
    fn sealed_anchor_provider_serves_all_live_anchors() {
        // Fixed-star public entry and output transformations.
        //
        // The apparent-place kernel supplies the selected observer-relative ray.
        // This layer applies equatorial, sidereal, defining-epoch, XYZ and radian
        // output flags and reports implied return bits. UT1 entries obtain TT from
        // the configured clock. A dynamic ayanamsha requires its anchor provider;
        // there is no substitution for unavailable data. Site position, motion,
        // precession and nutation remain part of the native observer model.
        use crate::sidereal::DynamicAnchor;
        use crate::sidereal::DynamicAnchorProvider;
        let Some(bytes) = anchor_test_fixture() else {
            eprintln!(
                "sealed_anchor_provider_serves_all_live_anchors: SKIP (missing LEB2 core fixture)"
            );
            return;
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let provider = crate::positions::SingleLeb::new(&index, &bytes);
        let anchors = SealedAnchorProvider::new(&provider);
        let all = [
            DynamicAnchor::Spica,
            DynamicAnchor::Revati,
            DynamicAnchor::Pushya,
            DynamicAnchor::Mula,
            DynamicAnchor::Aldebaran,
            DynamicAnchor::GalCenter,
            DynamicAnchor::GalNorthPole,
            DynamicAnchor::MulaWilhelm,
        ];
        for epoch in [crate::constants::J2000, crate::constants::J2000 + 36525.0] {
            for anchor in all {
                let mut first = None;
                for (noaberr, nogdefl) in
                    [(false, false), (true, false), (false, true), (true, true)]
                {
                    let lon = anchors
                        .anchor_longitude_deg(anchor, epoch, noaberr, nogdefl)
                        .unwrap_or_else(|e| panic!("{anchor:?} @ {epoch}: {e:?}"));
                    assert!(
                        lon.is_finite() && (0.0..360.0).contains(&lon),
                        "{anchor:?} @ {epoch}: unreduced {lon}"
                    );
                    // Determinism: the same call repeats identical.
                    let again = anchors
                        .anchor_longitude_deg(anchor, epoch, noaberr, nogdefl)
                        .expect("repeat serves");
                    assert_eq!(
                        lon.to_bits(),
                        again.to_bits(),
                        "{anchor:?} nondeterministic"
                    );
                    if anchor == DynamicAnchor::GalNorthPole {
                        match first {
                            None => first = Some(lon.to_bits()),
                            Some(word) => assert_eq!(
                                word,
                                lon.to_bits(),
                                "pole must ignore aberration/deflection toggles"
                            ),
                        }
                    }
                }
            }
        }
        // A non-finite epoch refuses before any provider state is read.
        match anchors.anchor_longitude_deg(DynamicAnchor::Spica, f64::NAN, false, false) {
            Err(crate::errors::Error::NonFiniteArgument { .. }) => {}
            other => panic!("expected NonFiniteArgument, got {other:?}"),
        }
    }

    #[test]
    fn sealed_anchor_provider_feeds_all_live_dispatch_arms() {
        // Fixed-star public entry and output transformations.
        //
        // The apparent-place kernel supplies the selected observer-relative ray.
        // This layer applies equatorial, sidereal, defining-epoch, XYZ and radian
        // output flags and reports implied return bits. UT1 entries obtain TT from
        // the configured clock. A dynamic ayanamsha requires its anchor provider;
        // there is no substitution for unavailable data. Site position, motion,
        // precession and nutation remain part of the native observer model.
        use crate::constants::{
            SIDM_ALDEBARAN_15TAU, SIDM_GALCENT_0SAG, SIDM_GALCENT_COCHRANE,
            SIDM_GALCENT_MULA_WILHELM, SIDM_GALCENT_RGILBRAND, SIDM_GALEQU_IAU1958,
            SIDM_GALEQU_MULA, SIDM_GALEQU_TRUE, SIDM_TRUE_CITRA, SIDM_TRUE_MULA, SIDM_TRUE_PUSHYA,
            SIDM_TRUE_REVATI, SIDM_TRUE_SHEORAN,
        };
        let Some(bytes) = anchor_test_fixture() else {
            eprintln!(
                "sealed_anchor_provider_feeds_all_live_dispatch_arms: SKIP (missing LEB2 core fixture)"
            );
            return;
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let provider = crate::positions::SingleLeb::new(&index, &bytes);
        let anchors = SealedAnchorProvider::new(&provider);
        let live = [
            SIDM_TRUE_CITRA,
            SIDM_TRUE_REVATI,
            SIDM_TRUE_PUSHYA,
            SIDM_TRUE_MULA,
            SIDM_TRUE_SHEORAN,
            SIDM_ALDEBARAN_15TAU,
            SIDM_GALCENT_0SAG,
            SIDM_GALCENT_RGILBRAND,
            SIDM_GALCENT_COCHRANE,
            SIDM_GALCENT_MULA_WILHELM,
            SIDM_GALEQU_IAU1958,
            SIDM_GALEQU_TRUE,
            SIDM_GALEQU_MULA,
        ];
        for mode in live {
            let (value, warning) = crate::sidereal::calc_mean_ayanamsha_deg(
                crate::constants::J2000,
                mode,
                (0.0, 0.0),
                0.0,
                false,
                false,
                false,
                Some(&anchors),
            )
            .unwrap_or_else(|e| panic!("mode {mode}: {e:?}"));
            assert!(
                value.is_finite() && (0.0..360.0).contains(&value),
                "mode {mode}: unreduced {value}"
            );
            assert!(warning.is_none(), "mode {mode}: unexpected {warning:?}");
        }
    }

    #[test]
    fn sealed_serving_serves_live_anchor_sidereal_star_place() {
        // A TRUE_CITRA stellar place consumes its live-anchor provider.
        // Check finite reduced longitude, deterministic reuse and a genuine
        // provider-versus-formula difference. Deflection-free anchors retain
        // their sampled nodes; deflected anchors use the actual request date.
        use crate::constants::{FLG_SIDEREAL, J2000, SIDM_TRUE_CITRA};
        let Some(bytes) = anchor_test_fixture() else {
            eprintln!(
                "sealed_serving_serves_live_anchor_sidereal_star_place: SKIP (missing LEB2 core fixture)"
            );
            return;
        };
        let index = crate::leb::FileIndex::decode(&bytes).expect("core LEB decodes");
        let composite = crate::leb::CompositeIndex::new(vec![(index.clone(), &bytes[..])])
            .expect("one reader must build");
        let nutation_at = |jd: f64| composite.eval_nutation(jd);
        let deltat_at = |jd: f64| {
            let cfg = crate::time::DeltaTConfig {
                model: crate::time::DeltaTModel::EspenakMeeus,
                user_override_days: None,
                tidal_acceleration: crate::constants::TIDAL_DEFAULT,
                iers_enabled: false,
                iers_lookup_seconds: None,
                smh2016_lookup_seconds: None,
            };
            crate::time::deltat_with_config(jd, &cfg)
        };
        let mut live = crate::sidereal::SiderealSelection::unset();
        let warnings = live.set_sidereal_mode(SIDM_TRUE_CITRA, 0.0, 0.0);
        assert!(
            warnings.is_empty(),
            "unexpected sidereal warnings: {warnings:?}"
        );
        let serving = SealedStarServing::new(&composite, &nutation_at, &deltat_at, live);
        let regulus = resolve_star_ref("Regulus").expect("Regulus resolves");
        let (place, name, _) = serving
            .fixstar_tt(&regulus, J2000, FLG_SIDEREAL, None)
            .unwrap_or_else(|e| panic!("live-anchor sidereal place refused: {e:?}"));
        assert_eq!(name, regulus.canonical);
        for (i, v) in place.iter().enumerate() {
            assert!(v.is_finite(), "slot {i}: non-finite {v}");
        }
        assert!(
            (0.0..360.0).contains(&place[0]),
            "lon unreduced: {}",
            place[0]
        );
        // Determinism: the same call repeats identical.
        let (again, _, _) = serving
            .fixstar_tt(&regulus, J2000, FLG_SIDEREAL, None)
            .expect("repeat serves");
        assert_eq!(place, again, "live-anchor serving nondeterministic");
        // Provider-fed, not skipped: the provider-free default selection
        // (Fagan/Bradley formula arm) serves a different longitude.
        let plain = SealedStarServing::new(
            &composite,
            &nutation_at,
            &deltat_at,
            crate::sidereal::SiderealSelection::unset(),
        );
        let (formula, _, _) = plain
            .fixstar_tt(&regulus, J2000, FLG_SIDEREAL, None)
            .expect("formula arm serves");
        assert_ne!(
            place[0].to_bits(),
            formula[0].to_bits(),
            "live anchor indistinguishable from the formula arm"
        );
    }

    #[test]
    fn catalog_grammar_uses_native_unicode_whitespace() {
        assert_eq!(
            parse_bayer_designation("\u{a0}Alpha\u{2003}Leonis\u{3000}").as_deref(),
            Some("alLeo")
        );
        assert_eq!(
            resolve_star_ref("\u{a0}Regulus\u{3000}").unwrap().star_id,
            1_000_001
        );
        assert_eq!(ref_star_key("Rigil Kentaurus"), "rigilkentaurus");
        assert_eq!(ref_star_key("Regulus,alLeo"), "regulus,alLeo");
        assert_eq!(ref_star_key("\u{1c}Regulus"), "\u{1c}regulus");
        assert!(matches!(
            resolve_star_ref("\u{1c}Regulus"),
            Err(StarResolveError::Miss(_))
        ));
    }
}

#[cfg(test)]
mod astropy_core_time_tests {
    use super::*;

    #[test]
    fn astropy_core_time_raw_star_constructor_keeps_zero_sign_before_propagation() {
        let row = StarData::from_row(-0.0, -0.0, 0.0, 0.0, 0.0, 0.0);
        let (initial, velocity) = star_icrs_initial_state_catalog_angles(&row);
        assert_eq!(initial[1].to_bits(), (-0.0_f64).to_bits());
        assert_eq!(initial[2].to_bits(), (-0.0_f64).to_bits());
        // The unchanged catalog propagation adds velocity*0 at J2000.
        // An arbitrary star must not take that addition before its own
        // TDB-plus-projection propagation.
        let (catalog, catalog_velocity) =
            star_icrs_state_catalog_angles(&row, crate::constants::J2000);
        assert_eq!(catalog[1].to_bits(), 0.0_f64.to_bits());
        assert_eq!(catalog[2].to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            velocity.map(f64::to_bits),
            catalog_velocity.map(f64::to_bits)
        );
    }
}
