//! House cusps, auxiliary angles, house positions and Gauquelin sectors.
//!
//! The kernels implement ecliptic, equatorial, semi-arc, horizon and vector
//! constructions. Inputs are degrees unless an entry explicitly uses radians;
//! epoch entries use UT1 Julian days, not UTC, and
//! the configured clock/frame. Ordinary systems return twelve cusps and eight
//! angles `[Asc, MC, ARMC, Vertex, EquAsc, CoAscKoch, CoAscMunkasey, PolarAsc]`.
//! Gauquelin `G` returns thirty-six sectors. Longitudes are normalized to
//! `[0, 360)`.
//!
//! Placidus, Koch and Gauquelin semi-arcs can be undefined in the polar
//! circle; their typed errors retain the requested system and latitude.
//! Explicit fallback entries select the requested alternate system and report
//! that choice. Speed entries differentiate the actual cusp/angle calculation
//! with the documented discontinuity and sign-boundary rules.
//!
//! | Code | `house_name` | Construction family | Polar-circle behavior | Implementation family |
//! | --- | --- | --- | --- | --- |
//! | `P` | Placidus | Time-based trisection of diurnal arcs | Raises (`PolarCircle`) | quadrant/time-based layer |
//! | `K` | Koch | Time-based (Koch birthplace arc + rising points) | Raises | Koch-kernel layer |
//! | `R` | Regiomontanus | Equatorial trisection, MC/ARMC flip when MC below horizon | Computes (latitude sign flip) | quadrant/time-based layer |
//! | `C` | Campanus | Prime-vertical trisection, MC/ARMC flip | Computes (latitude sign flip) | quadrant/time-based layer |
//! | `E`/`A` | equal | Ascendant + 30 deg steps (`A` uses the same algorithm) | Stable at all latitudes | simple-ecliptic layer |
//! | `W` | equal/ whole sign | Zodiac-sign floor of the Ascendant | Sign-pinned speeds; stable | simple-ecliptic layer |
//! | `O` | Porphyry | Quadrant trisection between ASC and MC | Stable; default fallback | simple-ecliptic layer |
//! | `B` | Alcabitius | Time-based semi-arc trisection | Computes | quadrant/time-based layer |
//! | `T` | Polich/Page | Topocentric (true topocentric Ascendant/MC) + MC flip | Computes | topocentric/specialized layer |
//! | `M` | Morinus | Equatorial division from ARMC (no latitude dependence) | Stable | quadrant/time-based layer |
//! | `X` | axial rotation system/Meridian houses | Meridian-plane division | Stable | quadrant/time-based layer |
//! | `V` | equal/Vehlow | Ascendant − 15 deg + 30 deg steps | Stable | simple-ecliptic layer |
//! | `H` | horizon/azimut | Horizon/azimuthal division | Computes; southern coasc2 branch at lat 0 | topocentric/specialized layer |
//! | `Y` | APC houses | Ascendant-parallel-circle vector construction | Computes; own MC-below-horizon post-flip | topocentric/specialized layer |
//! | `F` | Carter poli-equ. | Poli-equatorial construction | Computes | topocentric/specialized layer |
//! | `U` | Krusinski-Pisa-Goelzer | Ascendant–zenith great-circle vector construction | Computes | topocentric/specialized layer |
//! | `N` | equal/1=Aries | Natural-gradient (whole-sign from Aries 0) | Sign-pinned speeds; stable | topocentric/specialized layer |
//! | `G` | Gauquelin sectors | 36 diurnal sectors (not 12 cusps) | Raises | Gauquelin-sector layer |
//! | `g` | Gauquelin sectors | Same computation; **first 12 of the 36 sectors** | Raises in the polar circle | Gauquelin-sector layer |
//! | `S` | Sripati | Porphyry-like quadrant from ASC/MC | Computes | simple-ecliptic layer |
//! | `L` | Pullen SD | Sinusoidal-delta (neo-Porphyry) quadrant weighting | Computes | simple-ecliptic layer |
//! | `Q` | Pullen SR | Sinusoidal-ratio quadrant weighting | Computes | simple-ecliptic layer |
//! | `D` | equal (MC) | MC + 30 deg steps | Computes | simple-ecliptic layer |
//! | `I` | Sunshine | Treindl diurnal-arc division (needs Sun declination via `calc_ut`) | MC re-sync from cusp 10 | topocentric/specialized layer |
//! | `i` | Sunshine/alt. | Makransky spherical-trigonometry construction | MC re-sync from cusp 10 | topocentric/specialized layer |
//! | `J` | Savard-A | Horizon north/south-point position circles | Computes (latitude sign flip) | topocentric/specialized layer |
//! | other | `""` | Unknown selector: Placidus default computation | Raises like `P` in the polar circle | quadrant/time-based layer (default arm) |
//!
//! # Example: request an explicit polar-safe fallback
//!
//! ```
//! use libephemeris_core::houses;
//!
//! // ARMC, latitude, obliquity and the Sunshine-only declination are degrees.
//! let result = houses::houses_armc_with_fallback(
//!     123.45, 70.0, 23.44, 'P', 'O', true, 0.0,
//! ).expect("Porphyry fallback is available");
//! println!("House 1: {} degrees", result.cusps[0]);
//! println!("Fallback used: {}, note: {:?}", result.used_fallback, result.warning);
//! // Rust index 0 is house 1; there is no unused leading cusp slot.
//! ```

use crate::errors::Error;

/// Result type for house-inventory operations.
pub type Result<T> = std::result::Result<T, Error>;

/// Selector inventory: uppercase ASCII `A`–`Y` plus lowercase `i`.
///
/// Lowercase `g` is deliberately absent from this inventory but is an explicit
/// dispatch exception: Gauquelin geometry returns its first twelve sectors,
/// rather than the 36 returned by `G`. Both refuse in the polar circle. Other
/// selectors outside the inventory reach the Placidus-default computation.
pub const KNOWN_HSYS_CODES: &str = "ABCDEFGHIJKLMNOPQRSTUVWXYi";

/// Replacement character for an integer selector outside the Unicode
/// scalar range (`_UNKNOWN_HSYS_CHAR`): never a house-system code, so it
/// reaches the same Placidus default as every other unmapped value
/// instead of leaking a conversion error.
pub const UNKNOWN_HSYS_CHAR: char = '\u{FFFD}';

/// Latitude magnitude at or beyond which many systems lose accuracy
/// (`EXTREME_LATITUDE_THRESHOLD = 80.0`). Distinct from the polar
/// circle: an advisory warning domain, not a refusal.
pub const EXTREME_LATITUDE_THRESHOLD_DEG: f64 = 80.0;

/// Number of cusp longitudes for ordinary systems (houses 1–12).
pub const N_CUSPS: usize = 12;
/// Number of ASCMC auxiliary angles.
pub const N_ASCMC: usize = 8;
/// Number of Gauquelin sector longitudes (`G` return shape).
pub const N_GAUQUELIN_SECTORS: usize = 36;

/// Internal fixed-capacity cusp storage; public vector APIs allocate only at
/// their output boundary. Iteration visits the original 12/36 words in order.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct CuspBuffer {
    values: [f64; N_GAUQUELIN_SECTORS],
    len: usize,
}

impl CuspBuffer {
    pub(crate) fn as_slice(&self) -> &[f64] {
        &self.values[..self.len]
    }

    pub(crate) fn into_vec(self) -> Vec<f64> {
        self.into_iter().collect()
    }
}

impl AsRef<[f64]> for CuspBuffer {
    fn as_ref(&self) -> &[f64] {
        self.as_slice()
    }
}

impl std::ops::Deref for CuspBuffer {
    type Target = [f64];
    fn deref(&self) -> &Self::Target {
        self.as_slice()
    }
}

impl std::ops::DerefMut for CuspBuffer {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.values[..self.len]
    }
}

impl<const N: usize> From<[f64; N]> for CuspBuffer {
    fn from(values: [f64; N]) -> Self {
        values.into_iter().collect()
    }
}

impl FromIterator<f64> for CuspBuffer {
    fn from_iter<T: IntoIterator<Item = f64>>(iter: T) -> Self {
        let mut result = Self {
            values: [0.0; N_GAUQUELIN_SECTORS],
            len: 0,
        };
        for value in iter {
            assert!(result.len < N_GAUQUELIN_SECTORS, "internal cusp capacity");
            result.values[result.len] = value;
            result.len += 1;
        }
        result
    }
}

impl IntoIterator for CuspBuffer {
    type Item = f64;
    type IntoIter = std::iter::Take<std::array::IntoIter<f64, N_GAUQUELIN_SECTORS>>;
    fn into_iter(self) -> Self::IntoIter {
        self.values.into_iter().take(self.len)
    }
}

/// ASCMC slot indices.
pub const ASCMC_ASC: usize = 0;
/// ASCMC slot: Medium Coeli.
pub const ASCMC_MC: usize = 1;
/// ASCMC slot: ARMC (apparent sidereal time + longitude).
pub const ASCMC_ARMC: usize = 2;
/// ASCMC slot: Vertex.
pub const ASCMC_VERTEX: usize = 3;
/// ASCMC slot: Equatorial Ascendant (east point).
pub const ASCMC_EQUASC: usize = 4;
/// ASCMC slot: Co-Ascendant W. Koch (index 5).
pub const ASCMC_COASC_KOCH: usize = 5;
/// ASCMC slot: Co-Ascendant M. Munkasey (index 6).
pub const ASCMC_COASC_MUNKASEY: usize = 6;
/// ASCMC slot: Polar Ascendant M. Munkasey (index 7).
pub const ASCMC_POLAR_ASC: usize = 7;

/// Systems whose cusp speeds are sign-pinned (`_SIGN_PINNED_HSYS`):
/// Whole Sign `W` and Natural Gradient `N`. Used by the ARMC/time-stencil
/// speed layer, recorded here so the rule travels with the inventory.
pub const SIGN_PINNED_HSYS: [char; 2] = ['W', 'N'];

/// Systems that fail outright inside the polar circle
/// (`get_extreme_latitude_info` affected list).
pub const POLAR_AFFECTED_HSYS: [char; 3] = ['P', 'K', 'G'];
/// Systems that compute inside the polar circle but may be numerically
/// unstable at extreme latitudes (`>= 80°`; unstable list).
pub const EXTREME_UNSTABLE_HSYS: [char; 8] = ['C', 'R', 'T', 'B', 'H', 'U', 'Y', 'F'];
/// Systems reliable at all latitudes (stable list).
pub const STABLE_ALL_LATITUDE_HSYS: [char; 7] = ['E', 'W', 'O', 'M', 'X', 'V', 'N'];
/// Default fallback system for `houses_with_fallback` (Porphyry).
pub const DEFAULT_FALLBACK_HSYS: char = 'O';

/// A single cusp-sanity defect (`_validate_cusps` failure arms).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CuspDefect {
    /// Cusp `index` (0-based) is NaN.
    Nan {
        /// 0-based cusp position.
        index: usize,
    },
    /// Cusp `index` (0-based) is infinite.
    Infinite {
        /// 0-based cusp position.
        index: usize,
    },
    /// Cusp `index` (0-based) lies outside `[0, 360)`.
    OutOfRange {
        /// 0-based cusp position.
        index: usize,
        /// Offending value.
        value: f64,
    },
}

/// Whole-sign sign-cusp nudge (`_houses_whole_sign` ): `1e-9` deg
/// added before the `floor(asc / 30)` so an Ascendant within one float
/// ULP below an exact sign cusp (e.g. `180 − 1.4e-13` at cardinal ARMC)
/// resolves to the new sign, matching the reference.
pub const WHOLE_SIGN_NUDGE_DEG: f64 = 1e-9;

/// Vehlow offset: the Ascendant sits 15° into house 1.
pub const VEHLOW_OFFSET_DEG: f64 = 15.0;

/// Pullen SR equal-quadrant band (`calc_house_sizes_ratio` ): a
/// quadrant within `0.0001`° of 90° takes three exact 30° houses without
/// running the Newton solve.
pub const PULLEN_SR_EQUAL_BAND_DEG: f64 = 0.0001;

/// Pullen SR Newton-solve constants : at most 20 iterations,
/// `1e-12` derivative floor, `0.01` positivity clamp, `1e-10`
/// fixed-point stop, `1000.0` degenerate target for an empty quadrant,
/// `0.0001` divisor guard for the `x` recovery.
pub const PULLEN_SR_MAX_ITER: u32 = 20;
/// Pullen SR Newton derivative floor guarding the correction quotient.
pub const PULLEN_SR_DERIVATIVE_FLOOR: f64 = 1e-12;
/// Pullen SR positivity clamp on the solved radius.
pub const PULLEN_SR_R_MIN: f64 = 0.01;
/// Pullen SR fixed-point stop on successive radius iterates.
pub const PULLEN_SR_CONVERGENCE: f64 = 1e-10;
/// Pullen SR degenerate target for an empty quadrant.
pub const PULLEN_SR_EMPTY_QUADRANT_TARGET: f64 = 1000.0;
/// Pullen SR divisor guard for the `x` recovery.
pub const PULLEN_SR_DIVISOR_GUARD: f64 = 0.0001;

// Quadrant and time-based cusp kernels.

/// Placidus fixed-point convergence threshold in degrees
/// (`CONVERGENCE_THRESHOLD = 1e-7`, about 0.36 milliarcsec).
pub const PLACIDUS_CONVERGENCE_DEG: f64 = 1e-7;

/// Placidus fixed-point iteration budget (`MAX_ITERATIONS = 50`).
pub const PLACIDUS_MAX_ITERATIONS: u32 = 50;

/// Placidus bisection refinement budget (80 halvings, ).
pub const PLACIDUS_BISECTION_ITERATIONS: u32 = 80;

/// Placidus divergence warm-up (`iteration > 5`, ): the
/// runaway guard only arms after the sixth iteration.
pub const PLACIDUS_DIVERGENCE_WARMUP_ITERATIONS: u32 = 5;

/// Placidus divergence factor : bailing out when the wrapped
/// step grows beyond 1.5× the previous step.
pub const PLACIDUS_DIVERGENCE_FACTOR: f64 = 1.5;

/// Near-zero band for the RA→ecliptic conversion
/// (`_NEAR_ZERO = 1e-10`, ): polar pole heights and cardinal
/// snaps.
pub const RA_LON_NEAR_ZERO_DEG: f64 = 1e-10;

/// Singularity floor for the RA→ecliptic conversion : a
/// numerator and denominator both within `1e-12` take the continuous
/// right-ascension limit instead of `atan2` noise.
pub const RA_LON_SINGULARITY_FLOOR: f64 = 1e-12;

// Topocentric and specialized cusp kernels.

/// Three-component Cartesian vector in the equatorial frame of date
/// (`house_constructions.py` tuple convention).
type EqVec = [f64; 3];

/// Degenerate-anchor guard for house-circle construction
/// (`_TINY` ).
pub const HOUSE_CIRCLE_TINY: f64 = 1e-15;

/// Horizontal (azimuthal) near-zero guard (`_NEAR_ZERO`, local
/// to `_houses_horizontal` in the release).
pub const HORIZONTAL_NEAR_ZERO_DEG: f64 = 1e-10;

/// Sunshine/Treindl near-zero guard (`NEAR_ZERO`, local to
/// `_houses_sunshine` in the release).
pub const SUNSHINE_NEAR_ZERO_DEG: f64 = 1e-10;

// Gauquelin 36-sector cusps.

/// Gauquelin pole-height near-zero band (`NEAR_ZERO = 1e-10`, ):
/// guards the `tan ε` division and the equatorial-cusp branches.
pub const GAUQUELIN_NEAR_ZERO: f64 = 1e-10;

/// Gauquelin iterative-refinement convergence threshold
/// (`CONVERGENCE_THRESHOLD = 1e-8` degrees).
pub const GAUQUELIN_CONVERGENCE_DEG: f64 = 1e-8;

/// Gauquelin iterative-refinement iteration budget (`niter_max = 100`).
pub const GAUQUELIN_MAX_ITERATIONS: u32 = 100;

// Cusp/ASCMC speeds (`houses_armc_ex2`, `houses_ex2` stencils).

/// ARMC rotation rate per mean solar day (Meeus ch. 12,
/// `houses_armc_ex2` spelling): one extra rotation relative to the Sun.
pub const SIDEREAL_RATE_DEG_PER_DAY: f64 = 360.98564736629;

/// Centered-difference ARMC half-step: one sidereal second in degrees
/// ( `_SIDEREAL_RATE / 86400.0`).
pub const ARMC_SPEED_STEP_DEG: f64 = SIDEREAL_RATE_DEG_PER_DAY / 86400.0;

/// JD-based speed stencil half-step ( `_DT_DAYS = 1.0 / 4096.0`):
/// an exact binary fraction, so `tjdut ± h` stays symmetric at
/// JD ~2.45e6 where a one-second step would lose ~5 digits to roundoff.
pub const EX2_DT_DAYS: f64 = 1.0 / 4096.0;

/// JD-based Richardson fine half-step ( `_DT_HALF = 1.0 / 8192.0`).
pub const EX2_DT_HALF_DAYS: f64 = 1.0 / 8192.0;

/// Wrapped change across a speed stencil beyond which no finite
/// derivative exists ( `90.0` in both `fd_speed` and `_rate`).
pub const SPEED_DISCONTINUITY_DEG: f64 = 90.0;

/// Sign-pinned wheel step ( `30.0`): the `W`/`N` cusp-function
/// discontinuity arm.
pub const SIGN_PINNED_STEP_DEG: f64 = 30.0;

/// Sign-pinned step match tolerance ( `1e-6`).
pub const SIGN_PINNED_STEP_TOL_DEG: f64 = 1e-6;

/// ARMC-direct speed solution: cusps + ASCMC angles + cusp speeds +
/// ASCMC speeds in deg/day (the `houses_armc_ex2` / `houses_ex2` shape;
/// the cusp vector is length 12, length 36 for Gauquelin `G`).
pub type HouseSpeedSolution = (Vec<f64>, [f64; N_ASCMC], Vec<f64>, [f64; N_ASCMC]);

// Decimal house positions.

/// Near-zero threshold of the house-position kernels.
pub const HPOS_NEAR_ZERO_DEG: f64 = 1e-10;

/// Koch valid-fraction tolerance:
/// fractions infinitesimally outside the closed `[0, 2]` interval clamp
/// back into it; genuinely outside fractions stay undefined (`0.0`).
pub const KOCH_FRACTION_TOL: f64 = 1e-9;

/// Fractional clamp bound inside one house:
/// a body stays strictly inside its containing house.
pub const HOUSE_FRACTION_CAP: f64 = 0.9999999999;

/// Cusp-boundary offset added before the final house-circle wrap
pub const CUSP_BOUNDARY_OFFSET: f64 = 0.0;

/// Latitude classification for house calculations: whether `lat_deg` is
/// extreme (`|lat| >= 80°`, project convention), inside the polar
/// circle (`|lat| > 90° − obliquity`, where Placidus/Koch/Gauquelin
/// fail), the polar threshold, the extreme threshold, the hemisphere
/// letter, the failing/affected systems, the numerically unstable
/// systems and the systems reliable at all latitudes.
#[derive(Debug, Clone, PartialEq)]
pub struct ExtremeLatitudeInfo {
    /// Input latitude in degrees.
    pub latitude: f64,
    /// `|lat| >= 80°`.
    pub is_extreme: bool,
    /// `|lat| > 90° − obliquity`.
    pub is_polar_circle: bool,
    /// `90° − obliquity`.
    pub polar_threshold: f64,
    /// `80°` ([`EXTREME_LATITUDE_THRESHOLD_DEG`]).
    pub extreme_threshold: f64,
    /// `'N'` for `lat >= 0` (NaN reads the southern arm, like the
    /// release comparison), else `'S'`.
    pub hemisphere: char,
    /// `['P', 'K', 'G']` inside the polar circle, else empty.
    pub affected_systems: Vec<char>,
    /// The eight quadrant systems fragile past 80°, else empty.
    pub unstable_systems: Vec<char>,
    /// The seven systems reliable at all latitudes, always.
    pub stable_systems: Vec<char>,
}

/// Outcome of a fallback entry: cusps, ASCMC, whether the fallback
/// system served the result, and the informational warning.
#[derive(Debug, Clone, PartialEq)]
pub struct HouseFallbackOutcome {
    /// 12 cusps (36 for `G`), degnorm-shaped like the entry.
    pub cusps: Vec<f64>,
    /// 8 ASCMC angles.
    pub ascmc: [f64; N_ASCMC],
    /// True when the fallback system served the result.
    pub used_fallback: bool,
    /// Informational message when the fallback served or an extreme
    /// latitude behaved unstably; `None` otherwise.
    pub warning: Option<String>,
}

/// Borrowed synchronous providers for Julian-day house acquisition.
///
/// `frame_at` maps UT1 JD to `(armc0_deg, eps_true_deg)`; the sealed binding
/// obtains these through [`crate::time::house_armc_obliquity_deg`].
/// `sun_declination_at` maps `(UT1 JD, ephemeris_flags)` to equatorial solar
/// declination in degrees for Sunshine houses. `deltat_days_at` supplies
/// TT minus UT1 in days for the analytical Sun fallback. Unused callbacks
/// remain lazy; consumed errors and callback panics are not hidden.
pub struct HousesJdProviders<'a> {
    /// `(armc0_deg, eps_true_deg)` at a UT Julian Day.
    pub frame_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    /// Equatorial Sun declination at `(tjdut, ephemeris_flags)`.
    pub sun_declination_at: &'a dyn Fn(f64, i32) -> Result<f64>,
    /// TT−UT in days at a UT Julian Day.
    pub deltat_days_at: &'a dyn Fn(f64) -> Result<f64>,
}

/// Borrowed providers for sidereal house post-processing.
///
/// `ayanamsha_ex_at` maps `(UT1 JD, masked_flags)` to `(returned_flags,
/// ayanamsha_deg)`; the house reduction uses the value. The sealed binding
/// uses [`crate::sidereal::SealedAyanamshaServing::ex_ut`].
/// `mean_ut_plain_at` supplies mean ayanamsha in degrees at UT1 JD for plane
/// zero points. `nutation_at` supplies `(dpsi, deps)` in radians at TT JD.
/// `user_t0_tt` and `valens_epoch_tt` retain TT anchors for their respective
/// ECL_T0 projections. Callbacks are borrowed for the synchronous reduction
/// and can have caller-controlled state; they are not assumed pure.
pub struct HousesSiderealProviders<'a> {
    /// `(retflag, ayanamsha_deg)` at `(tjdut, masked_flags)`.
    pub ayanamsha_ex_at: &'a dyn Fn(f64, i32) -> Result<(i32, f64)>,
    /// Plain mean ayanamsha at a UT epoch (zero points).
    pub mean_ut_plain_at: &'a dyn Fn(f64) -> Result<f64>,
    /// `(dpsi_rad, deps_rad)` nutation at a TT epoch.
    pub nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    /// `SIDM_USER` epoch on TT (`SIDBIT_USER_UT` already applied).
    pub user_t0_tt: f64,
    /// Valens defining epoch on TT.
    pub valens_epoch_tt: f64,
}

/// Sidereal house feeds with defining epochs acquired only when consumed.
///
/// The first three feeds have the same units and flag contract as
/// [`HousesSiderealProviders`]. The epoch callbacks return Julian Day TT:
/// `user_t0_tt_at` applies the selected USER_UT model conversion and
/// `valens_epoch_tt_at` applies the Valens UT defining-epoch conversion.
/// Neither callback is used by tropical, fixed-epoch, ordinary ayanamsha,
/// or SSY-plane geometry. ECL_T0 reads only the selected USER or Valens
/// callback after the tropical seed and observation-frame feeds succeed.
/// All callbacks and any borrowed image snapshots must outlive the call.
pub struct HousesSiderealCallbacks<'a> {
    /// `(retflag, ayanamsha_deg)` at `(tjdut, masked_flags)`.
    pub ayanamsha_ex_at: &'a dyn Fn(f64, i32) -> Result<(i32, f64)>,
    /// Plain mean ayanamsha at a UT epoch, with default anchor toggles.
    pub mean_ut_plain_at: &'a dyn Fn(f64) -> Result<f64>,
    /// `(dpsi_rad, deps_rad)` at a TT observation epoch.
    pub nutation_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    /// Selected USER defining epoch in TT, evaluated on demand.
    pub user_t0_tt_at: &'a dyn Fn() -> Result<f64>,
    /// Selected Valens defining epoch in TT, evaluated on demand.
    pub valens_epoch_tt_at: &'a dyn Fn() -> Result<f64>,
}

impl HousesSiderealCallbacks<'_> {
    /// ECL_T0 plane epoch in TT for the selected base mode.
    ///
    /// USER and Valens consume their corresponding model callback. Other
    /// modes read the classical plane table, then the defining-epoch table,
    /// then J2000; they do not evaluate either unused model callback.
    /// A consumed provider error propagates before plane reconstruction.
    fn defining_epoch_tt(&self, mode: i32) -> Result<f64> {
        use crate::constants::{SIDM_USER, SIDM_VALENS_MOON};
        use crate::sidereal::{AYANAMSHA_DEFINING, ECL_T0_CLASSICAL_EPOCHS, J2000_TT};
        match mode {
            SIDM_USER => (self.user_t0_tt_at)(),
            SIDM_VALENS_MOON => (self.valens_epoch_tt_at)(),
            _ => Ok(ECL_T0_CLASSICAL_EPOCHS
                .iter()
                .find(|(id, _)| *id == mode)
                .map(|(_, epoch)| *epoch)
                .or_else(|| {
                    AYANAMSHA_DEFINING
                        .iter()
                        .find(|(id, _, _)| *id == mode)
                        .map(|(_, _, epoch)| *epoch)
                })
                .unwrap_or(J2000_TT)),
        }
    }
}

// Gauquelin-sector JD orchestration (methods 0-1)

/// Rise/set-time query for the Gauquelin methods 2-5 leg.
/// `body_id` is the provider-resolved integer body (star names stay with
/// the unserved star-geometry leg — the provider resolves integer ids only); `rsmi` carries
/// `CALC_RISE`/`CALC_SET` plus the release `BIT_DISC_CENTER` /
/// `BIT_NO_REFRACTION` spelling for the method.
pub struct GauquelinRiseSetQuery {
    /// Search start as a UT Julian Day.
    pub jd_start_ut: f64,
    /// Provider-resolved integer body id.
    pub body_id: i32,
    /// Event kind plus modifier bits.
    pub rsmi: i32,
    /// Observer longitude East, degrees.
    pub lon_deg_east: f64,
    /// Observer latitude, degrees.
    pub lat_deg: f64,
    /// Observer altitude above sea level, metres.
    pub height_m: f64,
    /// Atmospheric pressure, mbar (already defaulted: the release
    /// `gauquelin_sector` wrapper maps `0.0` to 1013.25 before the
    /// `rise_trans` barometric estimate could see it).
    pub atpress_mbar: f64,
    /// Atmospheric temperature, C (already defaulted: `0.0` to 15.0).
    pub attemp_c: f64,
    /// Ephemeris flags for the rise/set acquisition.
    pub flags: i32,
}

/// Gauquelin-sector providers: the methods 0-1 hour-angle legs plus the methods 2-5
/// rise/set-time leg (the rise/set sector-path leg).
pub struct GauquelinSectorProviders<'a> {
    /// UT Julian Day → `(armc0_deg, eps_true_deg)`.
    pub frame_at: &'a dyn Fn(f64) -> Result<(f64, f64)>,
    /// `(jd_ut, body_id, flags | FLG_SPEED)` → ecliptic `(lon, lat_body)`
    /// degrees.
    pub body_ecliptic_at: &'a dyn Fn(f64, i32, i32) -> Result<(f64, f64)>,
    /// Rise/set-time search → `(res, tret0)` with the release
    /// `rise_trans` contract (`res = 0` carries the event JD in UT,
    /// `res = -2` is the circumpolar/exhausted verdict). `None` keeps
    /// the sector-entry refusal: methods 2-5 report
    /// [`Error::GauquelinMethodUnserved`] instead of a silent
    /// method-0 substitution.
    pub rise_set_at: Option<&'a dyn Fn(GauquelinRiseSetQuery) -> Result<(i32, f64)>>,
}

// Private implementation boundaries; public logical paths remain unchanged.
mod geometry;
pub use geometry::angular_diff;
pub use geometry::fold_hsys_case;
pub use geometry::house_name;
pub use geometry::hsys_code_from_bytes;
pub use geometry::hsys_code_from_int;
pub use geometry::hsys_code_from_str;
pub use geometry::hsys_from_bytes;
pub use geometry::hsys_from_int;
pub use geometry::hsys_from_str;
pub use geometry::is_extreme_latitude;
pub use geometry::is_polar_circle;
pub use geometry::is_sign_pinned;
pub use geometry::polar_latitude_threshold_deg;
pub use geometry::polar_raise_applies;
pub use geometry::validate_cusps;
use geometry::wrap360;
mod systems;
pub use systems::calc_ascendant;
use systems::clamp_unit;
pub use systems::degnorm_deg;
pub use systems::difdeg2n_deg;
use systems::horizon_frame;
use systems::hour_circle_lon_deg;
use systems::house_circle_lon_deg;
pub use systems::houses_alcabitius;
pub use systems::houses_campanus;
pub use systems::houses_carter;
pub use systems::houses_equal;
pub use systems::houses_equal_mc;
pub use systems::houses_meridian;
pub use systems::houses_morinus;
pub use systems::houses_placidus;
pub use systems::houses_polich_page;
pub use systems::houses_porphyry;
pub use systems::houses_pullen_sd;
pub use systems::houses_pullen_sr;
pub use systems::houses_quadrant;
pub use systems::houses_regiomontanus;
pub use systems::houses_sripati;
pub use systems::houses_vehlow;
pub use systems::houses_whole_sign;
use systems::init_cardinal_cusps_13;
#[cfg(test)]
use systems::pullen_sr_quadrant_sizes;
pub use systems::ra_to_ecliptic_longitude;
use systems::semiarc_cusp_hour_angles;
use systems::set_opposite_cusps_13;
use systems::sin_cos;
use systems::strip_index0;
use systems::unit_radec;
use systems::vec_cross;
use systems::vec_dot;
use systems::vec_normalized;
mod specialized_systems;
pub use specialized_systems::apc_cusp;
pub use specialized_systems::apc_cusp_checked;
pub use specialized_systems::armc_fd_speed_deg_per_day;
pub use specialized_systems::asc_from_armc;
pub use specialized_systems::calc_lat_for_entry;
pub use specialized_systems::calc_vertex;
pub use specialized_systems::coasc_koch_deg;
pub use specialized_systems::coasc_munkasey_deg;
pub use specialized_systems::compute_ascmc;
pub use specialized_systems::equatorial_ascendant_deg;
pub use specialized_systems::ex2_angle_rate_richardson_deg_per_day;
pub use specialized_systems::ex2_cusp_rate_deg_per_day;
pub use specialized_systems::houses_apc;
pub use specialized_systems::houses_armc;
pub(crate) use specialized_systems::houses_armc_buffer;
pub use specialized_systems::houses_armc_ex2;
pub use specialized_systems::houses_armc_from_normalized;
pub use specialized_systems::houses_ex2_via;
pub use specialized_systems::houses_gauquelin;
pub use specialized_systems::houses_gauquelin_12;
pub use specialized_systems::houses_horizontal;
pub use specialized_systems::houses_koch;
pub use specialized_systems::houses_krusinski;
pub use specialized_systems::houses_natural_gradient;
pub use specialized_systems::houses_savard_a;
pub use specialized_systems::houses_specialized;
pub use specialized_systems::houses_sunshine;
pub use specialized_systems::houses_sunshine_makransky;
pub use specialized_systems::mc_from_armc;
pub use specialized_systems::polar_ascendant_deg;
pub use specialized_systems::resolve_armc_active;
pub use specialized_systems::sun_declination_analytic;
pub use specialized_systems::validate_houses_armc_latitude;
#[cfg(test)]
use specialized_systems::wrap_speed_diff_deg;
mod serving;
use serving::above_horizon_test;
use serving::acos_deg;
use serving::asc_diff_saturated;
pub use serving::ascendant_on_eastern_horizon;
use serving::asin_deg;
use serving::atan_deg;
use serving::cos_deg;
pub use serving::extreme_latitude_info;
pub use serving::fallback_system_name;
pub use serving::houses_armc_with_fallback;
pub use serving::houses_ex_via;
pub use serving::houses_ex_via_callbacks;
pub(crate) use serving::houses_ex_via_callbacks_model_clock;
pub use serving::houses_ex_with_context;
pub use serving::houses_fallback_via;
pub use serving::houses_jd_via;
pub(crate) use serving::houses_jd_via_model_clock;
use serving::rotate_frame;
pub use serving::sealed_sun_declination_deg;
pub use serving::sidereal_selection_for_context;
use serving::sin_deg;
use serving::tan_deg;
pub use serving::validate_ecliptic_longitude;
pub use serving::validate_houses_coordinates;
pub(crate) use serving::{
    houses_ex_buffer_via, houses_ex_buffer_via_callbacks, houses_ex2_buffer_via,
};
mod house_position;
pub use house_position::gauquelin_sector_from_rise_set_via;
pub use house_position::gauquelin_sector_via;
pub use house_position::house_pos;
pub use house_position::house_pos_from_normalized;
use house_position::transpose3_local;

#[cfg(test)]
#[path = "tests/tests.rs"]
mod tests;

#[cfg(test)]
#[path = "native_callbacks_tests.rs"]
mod native_callbacks_tests;
