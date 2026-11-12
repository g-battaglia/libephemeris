//! Sidereal modes, defining epochs, ayanamshas and position projections.
//!
//! Predefined selections combine defining-epoch pairs, long-term precession,
//! fixed-epoch frames, and live stellar or galactic anchors. A context owns its
//! [`SiderealSelection`]; dynamic modes require an explicit anchor provider.
//! Projection flags choose an ecliptic, equatorial or invariable reference plane.
//!
//! Evaluation epochs are TT Julian days; UT entries use UT1 and apply their
//! configured Delta-T conversion. A caller's defining epoch also uses UT1 when
//! `SIDBIT_USER_UT` is set. Ayanamshas and longitude offsets are degrees, and state rates are
//! per TT day. Nutation follows the selected actual model or bound polynomial;
//! value and rate calculations retain that same frame. Fixed-epoch projections
//! and of-date corrections are documented at their entries.
//!
//! The scalar ayanamsha dispatcher warns for unknown modes and uses the
//! Fagan/Bradley defining pair. Missing required anchors and checked nonfinite
//! evaluation epochs return typed errors; storing a selection alone does not
//! validate its epochs. These conventional zodiac definitions do not replace
//! the underlying position's precision. Bound callbacks are synchronous and
//! may perform I/O or panic; adapters do not imply callback purity.

use crate::constants::{
    FLG_BARYCTR, FLG_EQUATORIAL, FLG_HELCTR, FLG_J2000, FLG_JPLEPH, FLG_MOSEPH, FLG_NOABERR,
    FLG_NOGDEFL, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_SPEED, FLG_SWIEPH, FLG_TRUEPOS, FLG_XYZ,
    SIDBIT_ECL_DATE, SIDBIT_ECL_T0, SIDBIT_NO_PREC_OFFSET, SIDBIT_PREC_ORIG, SIDBIT_SSY_PLANE,
    SIDBIT_USER_UT, SIDM_ALDEBARAN_15TAU, SIDM_ARYABHATA, SIDM_ARYABHATA_522, SIDM_ARYABHATA_MSUN,
    SIDM_B1950, SIDM_BABYL_BRITTON, SIDM_BABYL_ETPSC, SIDM_BABYL_HUBER, SIDM_BABYL_KUGLER1,
    SIDM_BABYL_KUGLER2, SIDM_BABYL_KUGLER3, SIDM_DELUCE, SIDM_DJWHAL_KHUL, SIDM_FAGAN_BRADLEY,
    SIDM_GALALIGN_MARDYKS, SIDM_GALCENT_0SAG, SIDM_GALCENT_COCHRANE, SIDM_GALCENT_MULA_WILHELM,
    SIDM_GALCENT_RGILBRAND, SIDM_GALEQU_FIORENZA, SIDM_GALEQU_IAU1958, SIDM_GALEQU_MULA,
    SIDM_GALEQU_TRUE, SIDM_HIPPARCHOS, SIDM_J1900, SIDM_J2000, SIDM_JN_BHASIN, SIDM_KRISHNAMURTI,
    SIDM_KRISHNAMURTI_VP291, SIDM_LAHIRI, SIDM_LAHIRI_1940, SIDM_LAHIRI_ICRC, SIDM_LAHIRI_VP285,
    SIDM_RAMAN, SIDM_SASSANIAN, SIDM_SS_CITRA, SIDM_SS_REVATI, SIDM_SURYASIDDHANTA,
    SIDM_SURYASIDDHANTA_MSUN, SIDM_TRUE_CITRA, SIDM_TRUE_MULA, SIDM_TRUE_PUSHYA, SIDM_TRUE_REVATI,
    SIDM_TRUE_SHEORAN, SIDM_USER, SIDM_USHASHASHI, SIDM_VALENS_MOON, SIDM_YUKTESHWAR,
};
use crate::errors::{Error, Warning};
use crate::stars::StarData;

/// Mask isolating the base ayanamsha mode from a `SIDBIT`-combined ID.
///
/// the model applies `mode & 0xFF` (planets `_calc_ayanamsa` reducer,
/// `get_ayanamsa_name` lookup); high bits are the projection flags.
pub const BASE_MODE_MASK: i32 = 0xFF;

/// Highest predefined base mode ID (`NSIDM_PREDEF - 1`); base IDs `0..=46`
/// are all assigned, `47..=254` are unassigned.
pub const MAX_PREDEFINED_MODE: i32 = 46;

/// SIDBIT projection flags with a distinct realization on the sidereal
/// path: frame projections
/// (`SIDBIT_ECL_T0`/`SIDBIT_SSY_PLANE`, handled by the frame-projection layer), the user-epoch
/// timescale flag (`SIDBIT_USER_UT`) and the ecliptic-of-date flag
/// (`SIDBIT_ECL_DATE`).
pub const APPLIED_SIDBITS: i32 =
    SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE | SIDBIT_USER_UT | SIDBIT_ECL_DATE;

/// SIDBIT flags accepted silently for compatibility but currently reducing
/// to the base ayanamsha value.
/// Accepted silently because both are no-ops at this layer: `PREC_ORIG`
/// changes results only below the ~0.1″ ayanamsha floor, and
/// `NO_PREC_OFFSET` (a six-mode time-constant correction up to ~0.83″)
/// is not reconstructible without fitting the reference output.
pub const ACCEPTED_SIDBITS: i32 = SIDBIT_NO_PREC_OFFSET | SIDBIT_PREC_ORIG;

/// Union retained by `set_sid_mode`.
/// Bits outside this union are forwarded to the state-level setter, which
/// strips and reports genuinely unknown high bits.
pub const IMPLEMENTED_SIDBITS: i32 = APPLIED_SIDBITS | ACCEPTED_SIDBITS;

/// Predefined epoch-pair modes: `(mode, ayanamsha at defining epoch
/// [deg], defining epoch [TT JD])`.
///
/// Each pair is anchored at its defining TT epoch, where accumulated
/// precession vanishes. Definitions, historical ambiguities and source-status
/// labels are documented in `beginners-manual/en/11-sidereal-zodiac.md`.
/// Per-entry literature provenance: Fagan/Bradley SVP;
/// Lahiri (Calendar Reform Committee); Babylonian modes (van der
/// Waerden/Huber tabulation); Siddhanta modes (Burgess/Clark readings);
/// Britton (formula derivation); galactic-center modes (Reid & Brunthaler
/// position); Sheoran (2017 anchor rationale).
pub const AYANAMSHA_DEFINING: [(i32, f64, f64); 32] = [
    (SIDM_FAGAN_BRADLEY, 24.042044444444457, 2433282.42345905),
    (SIDM_LAHIRI, 23.85709166666667, 2451545.0),
    (SIDM_DELUCE, 26.413055555555555, 2415020.31352),
    (SIDM_RAMAN, 21.014166666666668, 2415020.0),
    (SIDM_USHASHASHI, 0.0, 1925312.0),
    (SIDM_KRISHNAMURTI, 22.36388888888889, 2415020.0),
    (SIDM_DJWHAL_KHUL, 30.0, 2494279.0),
    (SIDM_YUKTESHWAR, 20.91, 2412908.1187773715),
    (SIDM_JN_BHASIN, 0.0, 1854089.0),
    (SIDM_BABYL_KUGLER1, 0.0, 1834000.0),
    (SIDM_BABYL_KUGLER2, 0.0, 1797110.0),
    (SIDM_BABYL_KUGLER3, 0.0, 1774464.0),
    (SIDM_BABYL_HUBER, -4.466666666666667, 1684613.0),
    (SIDM_BABYL_ETPSC, 0.0, 1807702.0),
    (SIDM_HIPPARCHOS, 0.0, 1920199.0),
    (SIDM_SASSANIAN, 0.0, 1927136.0),
    (SIDM_J2000, 0.0, 2451545.0),
    (SIDM_J1900, 0.0, 2415020.0),
    (SIDM_B1950, 0.0, 2433282.42345905),
    (SIDM_SURYASIDDHANTA, 0.0, 1903395.1338202194),
    (
        SIDM_SURYASIDDHANTA_MSUN,
        -0.21510829031467438,
        1903396.8775657457,
    ),
    (SIDM_ARYABHATA, 0.0, 1903395.1338202194),
    (
        SIDM_ARYABHATA_MSUN,
        -0.23810671619139612,
        1903396.8542324193,
    ),
    (SIDM_SS_REVATI, -0.7922573227712064, 1903395.1338202194),
    (SIDM_SS_CITRA, 2.1105983371992068, 1903395.1338202194),
    (SIDM_ARYABHATA_522, 0.0, 1911798.0),
    (SIDM_BABYL_BRITTON, 0.0, 1805145.0),
    (SIDM_GALEQU_FIORENZA, 25.0, 2451544.5),
    (SIDM_LAHIRI_1940, 0.0, 1825599.0),
    (SIDM_LAHIRI_VP285, 0.0, 1825233.3060523525),
    (SIDM_KRISHNAMURTI_VP291, 0.0, 1827424.7590899938),
    (SIDM_LAHIRI_ICRC, 23.245339462821466, 2435553.5),
];

/// Modes evaluated from live catalog/frame geometry in `planets.py`
/// : no defining pair above.
pub const DYNAMIC_AYANAMSHA_MODES: [i32; 15] = [
    SIDM_ALDEBARAN_15TAU,
    SIDM_GALCENT_0SAG,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_GALALIGN_MARDYKS,
    SIDM_TRUE_MULA,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_TRUE_SHEORAN,
    SIDM_GALCENT_COCHRANE,
    SIDM_VALENS_MOON,
];

/// Catalog anchor each dynamic mode is evaluated from.
///
/// `TRUE_SHEORAN` shares the `PUSHYA` (delta Cancri) entry with
/// `TRUE_PUSHYA`: the model evaluates both from `STARS["PUSHYA"]` with
/// different target longitudes.
pub const DYNAMIC_MODE_ANCHORS: [(i32, &str); 15] = [
    (SIDM_TRUE_CITRA, "SPICA"),
    (SIDM_TRUE_REVATI, "REVATI"),
    (SIDM_TRUE_PUSHYA, "PUSHYA"),
    (SIDM_TRUE_MULA, "MULA"),
    (SIDM_TRUE_SHEORAN, "PUSHYA"),
    (SIDM_ALDEBARAN_15TAU, "ALDEBARAN"),
    (SIDM_GALCENT_0SAG, "GAL_CENTER"),
    (SIDM_GALCENT_RGILBRAND, "GAL_CENTER"),
    (SIDM_GALCENT_COCHRANE, "GAL_CENTER"),
    (SIDM_GALCENT_MULA_WILHELM, "GAL_CENTER"),
    (SIDM_GALEQU_IAU1958, "GAL_NORTH_POLE"),
    (SIDM_GALEQU_TRUE, "GAL_NORTH_POLE"),
    (SIDM_GALEQU_MULA, "GAL_NORTH_POLE"),
    (SIDM_GALALIGN_MARDYKS, "IAU2006_PRECESSION_FORMULA"),
    (SIDM_VALENS_MOON, "EPOCH_PAIR"),
];

/// Galactic-center modes with their published sidereal target longitude of
/// Sgr A* in degrees.
pub const GALCENT_TARGET_LON: [(i32, f64); 3] = [
    (SIDM_GALCENT_0SAG, 240.0),
    (SIDM_GALCENT_RGILBRAND, 244.37694101250946),
    (SIDM_GALCENT_COCHRANE, 270.0),
];

/// Galactic-equator modes evaluated from the live galactic north pole,
/// with their ascending-node targets in degrees.
pub const GALEQU_TARGET_LON: [(i32, f64); 3] = [
    (SIDM_GALEQU_IAU1958, 240.0),
    (SIDM_GALEQU_TRUE, 239.94708),
    (SIDM_GALEQU_MULA, 246.6137),
];

/// Wilhelm's Mula target: Sgr A* hour-circle (dhruva) projection held at
/// the middle of Mula, 246°40'.
pub const MULA_WILHELM_TARGET_LON: f64 = 246.66666666666666;

/// Aldebaran anchor: held at 15° Taurus, 45° absolute.
pub const ALDEBARAN_TARGET_LON: f64 = 45.0;

/// Sheoran anchor: delta Cancri held at 103°29'32.9375" sidereal
pub const SHEORAN_TARGET_LON: f64 = 103.49248263888889;

/// Baked J2000 anchor astrometry for the live star/galactic ayanamsha
/// arms.
///
/// Each entry carries ICRS J2000.0 right ascension and declination in
/// degrees, proper motion in arcsec/yr (`pm_ra` is the great-circle rate
/// mu_alpha*, i.e. it already includes the cos(dec) factor, exactly as
/// the model documents), trigonometric parallax in milliarcseconds
/// and radial velocity in km/s (positive receding). the model stores
/// parallax in arcseconds and converts with `parallax * 1000.0` where
/// consumed; every entry below was checked to be exactly identical to
/// that product, so the milliarcsecond literals carry the model values
/// exactly.
///
/// Provenance per entry: Spica (alpha Virginis,
/// HIP 65474, Hipparcos), Revati (zeta Piscium A, HIP 5737, Gaia DR 3),
/// Aldebaran (alpha Tauri, HIP 21421, Hipparcos/van Leeuwen 2007),
/// Pushya (delta Cancri, HIP 42911, Gaia DR 3; shared by the Sheoran
/// mode), Mula (lambda Scorpii, HIP 85927, Hipparcos 2007 — too bright
/// for Gaia DR 3), the Galactic Center (Sgr A*, Reid & Brunthaler 2004
/// position with the 2020 apparent proper motion) and the Galactic
/// North Pole (J2000, fixed).
///
/// The Galactic Center `pm_ra` needs care: the tabulated
/// -3.151 mas/yr is the RA *coordinate* rate, so the model stores
/// the great-circle rate `-0.003151 * cos(radians(-29.00781056))`
/// (-2.7557 mas/yr); the other entries already store catalog mu_alpha*
/// values. Entries without a parallax or radial velocity carry
/// `0.0` (unknown parallax is represented at the one-gigaparsec floor
/// by the propagation layer, not here).
///
/// This table is data only: the apparent-place reduction (proper-motion
/// propagation, Earth vector, aberration/deflection, precession, node
/// cache) is fixed-star-pipeline scope, so the 13 live arms keep
/// refusing with [`Error::AyanamshaAnchorUnserved`] until a
/// [`DynamicAnchorProvider`] is served.
pub const ANCHOR_STARS: [(&str, StarData); 7] = [
    (
        "SPICA",
        StarData::from_row(201.2982475, -11.1613219, -0.04250, -0.03173, 12.44, 1.0),
    ),
    (
        "REVATI",
        StarData::from_row(
            18.4328583349,
            7.5753601597,
            0.142693,
            -0.053051,
            24.4595,
            15.0,
        ),
    ),
    (
        "ALDEBARAN",
        StarData::from_row(68.98016279, 16.509302351, 0.06345, -0.18894, 48.94, 54.26),
    ),
    (
        "PUSHYA",
        StarData::from_row(
            131.1712460977,
            18.1543080691,
            -0.018435,
            -0.227813,
            23.8271,
            17.14,
        ),
    ),
    (
        "MULA",
        StarData::from_row(263.40216717, -37.10382356, -0.00853, -0.03080, 5.71, -3.00),
    ),
    (
        "GAL_CENTER",
        // -0.003151 * cos(radians(-29.00781056)): the model expression
        // evaluated once; the literal below is its exact value.
        StarData::from_row(
            266.41683708,
            -29.00781056,
            -0.002755718424389021,
            -0.005547,
            0.0,
            0.0,
        ),
    ),
    (
        "GAL_NORTH_POLE",
        StarData::from_row(192.85948, 27.12825, 0.0, 0.0, 0.0, 0.0),
    ),
];

/// Look up a baked anchor entry by its `STARS` key (see
/// [`DYNAMIC_MODE_ANCHORS`]). Returns `None` for unknown keys,
/// including the two formula markers (`"IAU2006_PRECESSION_FORMULA"`,
/// `"EPOCH_PAIR"`), which name no catalog entry by design.
pub fn anchor_star(name: &str) -> Option<StarData> {
    ANCHOR_STARS
        .iter()
        .find(|(key, _)| *key == name)
        .map(|(_, star)| *star)
}

/// Mardyks defining pair: ayanamsha exactly 30° at the September equinox
/// 1998.
pub const MARDYKS_DEFINING: (f64, f64) = (30.0, 2451079.733114197);

/// Valens lunar-epoch mode: defining UT epoch (JD) and ayanamsha offset
/// in degrees.
pub const VALENS_MOON_T0_UT: f64 = 1775845.5;
/// Valens defining ayanamsha at [`VALENS_MOON_T0_UT`], in degrees.
pub const VALENS_MOON_AYAN_T0_DEG: f64 = -2.9422;

/// Fixed-epoch sidereal modes: positions are expressed on the mean ecliptic/equinox of a
/// reference epoch rather than by a scalar longitude offset. Frame math
/// lives in the precession/projection layers, not in this table.
pub const FIXED_EPOCH_MODES: [i32; 4] = [SIDM_J2000, SIDM_J1900, SIDM_B1950, SIDM_GALALIGN_MARDYKS];

/// Constant ecliptic-longitude offset of a fixed-epoch frame in degrees,
/// subtracted from the ecliptic channels only.
pub fn fixed_epoch_lon_offset_deg(mode: i32) -> f64 {
    if mode == SIDM_GALALIGN_MARDYKS {
        MARDYKS_DEFINING.0
    } else {
        0.0
    }
}

/// True if `mode` is a fixed-epoch request mode (base ID; caller strips
/// SIDBITs first).
pub fn is_fixed_epoch_mode(mode: i32) -> bool {
    FIXED_EPOCH_MODES.contains(&mode)
}

/// Live star / galactic-center modes whose anchor carries annual aberration
/// and solar deflection in its default apparent-place realization, so the
/// `FLG_TRUEPOS`/`FLG_NOABERR`/`FLG_NOGDEFL` request flags (and
/// helio/barycentric centers) shift the subtracted zero point. Inert for every other mode.
pub const ABERRANT_ANCHOR_MODES: [i32; 9] = [
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_TRUE_MULA,
    SIDM_TRUE_SHEORAN,
    SIDM_GALCENT_0SAG,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALCENT_COCHRANE,
    SIDM_GALCENT_MULA_WILHELM,
];

/// Star/galactic "calculated" modes for which `get_ayanamsa_ex[_ut]`
/// drops `FLG_NONUT` from the echoed return flag.
pub const AYANAMSA_EX_NONUT_DROP_MODES: [i32; 12] = [
    SIDM_GALCENT_0SAG,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_TRUE_MULA,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_TRUE_SHEORAN,
    SIDM_GALCENT_COCHRANE,
];

/// Solar-system invariable plane relative to the J2000 mean ecliptic:
/// inclination and ascending-node longitude in degrees from Souami &
/// Souchay (2012), A&A 543, A 133.
pub const INVARIABLE_PLANE_INCL_DEG: f64 = 1.578690;
/// Node of the invariable plane in degrees; see [`INVARIABLE_PLANE_INCL_DEG`].
pub const INVARIABLE_PLANE_NODE_DEG: f64 = 107.582222;

/// Predefined sidereal-mode names in canonical ID order.
const AYANAMSHA_NAMES: [(i32, &str); 47] = [
    (SIDM_FAGAN_BRADLEY, "Fagan/Bradley"),
    (SIDM_LAHIRI, "Lahiri"),
    (SIDM_DELUCE, "De Luce"),
    (SIDM_RAMAN, "Raman"),
    (SIDM_USHASHASHI, "Usha/Shashi"),
    (SIDM_KRISHNAMURTI, "Krishnamurti"),
    (SIDM_DJWHAL_KHUL, "Djwhal Khul"),
    (SIDM_YUKTESHWAR, "Yukteshwar"),
    (SIDM_JN_BHASIN, "J.N. Bhasin"),
    (SIDM_BABYL_KUGLER1, "Babylonian/Kugler 1"),
    (SIDM_BABYL_KUGLER2, "Babylonian/Kugler 2"),
    (SIDM_BABYL_KUGLER3, "Babylonian/Kugler 3"),
    (SIDM_BABYL_HUBER, "Babylonian/Huber"),
    (SIDM_BABYL_ETPSC, "Babylonian/Eta Piscium"),
    (SIDM_ALDEBARAN_15TAU, "Babylonian/Aldebaran = 15 Tau"),
    (SIDM_HIPPARCHOS, "Hipparchos"),
    (SIDM_SASSANIAN, "Sassanian"),
    (SIDM_GALCENT_0SAG, "Galact. Center = 0 Sag"),
    (SIDM_J2000, "J2000"),
    (SIDM_J1900, "J1900"),
    (SIDM_B1950, "B1950"),
    (SIDM_SURYASIDDHANTA, "Suryasiddhanta"),
    (SIDM_SURYASIDDHANTA_MSUN, "Suryasiddhanta, mean Sun"),
    (SIDM_ARYABHATA, "Aryabhata"),
    (SIDM_ARYABHATA_MSUN, "Aryabhata, mean Sun"),
    (SIDM_SS_REVATI, "SS Revati"),
    (SIDM_SS_CITRA, "SS Citra"),
    (SIDM_TRUE_CITRA, "True Citra"),
    (SIDM_TRUE_REVATI, "True Revati"),
    (SIDM_TRUE_PUSHYA, "True Pushya (PVRN Rao)"),
    (SIDM_GALCENT_RGILBRAND, "Galactic Center (Gil Brand)"),
    (SIDM_GALEQU_IAU1958, "Galactic Equator (IAU1958)"),
    (SIDM_GALEQU_TRUE, "Galactic Equator"),
    (SIDM_GALEQU_MULA, "Galactic Equator mid-Mula"),
    (SIDM_GALALIGN_MARDYKS, "Skydram (Mardyks)"),
    (SIDM_TRUE_MULA, "True Mula (Chandra Hari)"),
    (
        SIDM_GALCENT_MULA_WILHELM,
        "Dhruva/Gal.Center/Mula (Wilhelm)",
    ),
    (SIDM_ARYABHATA_522, "Aryabhata 522"),
    (SIDM_BABYL_BRITTON, "Babylonian/Britton"),
    (SIDM_TRUE_SHEORAN, "\"Vedic\"/Sheoran"),
    (SIDM_GALCENT_COCHRANE, "Cochrane (Gal.Center = 0 Cap)"),
    (SIDM_GALEQU_FIORENZA, "Galactic Equator (Fiorenza)"),
    (SIDM_VALENS_MOON, "Vettius Valens"),
    (SIDM_LAHIRI_1940, "Lahiri 1940"),
    (SIDM_LAHIRI_VP285, "Lahiri VP285"),
    (SIDM_KRISHNAMURTI_VP291, "Krishnamurti-Senthilathiban"),
    (SIDM_LAHIRI_ICRC, "Lahiri ICRC"),
];

/// Name of a sidereal mode.
///
/// Compatibility contract: only the predefined modes `0..=46` carry a
/// name; every other ID (the unassigned block, `SIDM_USER`, anything
/// wider) yields `""`. The `SIDBIT_*` projection flags occupy bits >= 8,
/// so the lookup is driven by the low byte only, and the mask is a plain
/// two's-complement `& 0xFF`, identical for negative IDs in Rust and
/// for either sign.
pub fn ayanamsha_name(mode: i32) -> &'static str {
    let base = mode & BASE_MODE_MASK;
    for (id, name) in AYANAMSHA_NAMES {
        if id == base {
            return name;
        }
    }
    ""
}

/// Strip a possibly `SIDBIT`-combined mode to its base ayanamsha ID, with
/// the reducer rule: non-negative IDs mask to the
/// low byte, negative IDs pass through unchanged so the invalid-mode
/// fallback (Fagan/Bradley, applied by the value layer) applies instead of wrapping onto
/// `SIDM_USER`.
pub fn strip_to_base_mode(mode: i32) -> i32 {
    if mode >= 0 {
        mode & BASE_MODE_MASK
    } else {
        mode
    }
}

/// True when `mode` names a computable ayanamsha: a non-negative ID whose
/// base is a predefined mode `0..=46` (defining pair or dynamic anchor) or
/// `SIDM_USER` (caller-supplied epoch pair). Anything else reaches the
/// the model unknown-mode fallback (warn + Fagan/Bradley) in the value layer.
pub fn is_registered_mode(mode: i32) -> bool {
    if mode < 0 {
        return false;
    }
    let base = mode & BASE_MODE_MASK;
    base <= MAX_PREDEFINED_MODE || base == SIDM_USER
}

/// Context-local sidereal selection: the mode ID, the `SIDM_USER`
/// reference epoch/ayanamsha pair and the retained `SIDBIT` flags.
///
/// This owned value does not mutate process-global configuration. Negative
/// and unassigned IDs can be stored; evaluation applies its documented warning
/// and fallback. Defining `t0` is TT JD, or UT1 JD with `SIDBIT_USER_UT`;
/// `ayan_t0` is degrees. Both scalars are stored without finiteness validation.
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct SiderealSelection {
    mode: i32,
    t0: f64,
    ayan_t0: f64,
    bits: i32,
}

impl SiderealSelection {
    /// Never-set state: reads back as Fagan/Bradley with a zero epoch
    /// pair and no projection flags.
    pub fn unset() -> Self {
        Self {
            mode: SIDM_FAGAN_BRADLEY,
            t0: 0.0,
            ayan_t0: 0.0,
            bits: 0,
        }
    }

    /// Store a combined mode/flag selection and return its projection warnings.
    ///
    /// For nonnegative IDs, retains [`IMPLEMENTED_SIDBITS`] and strips unknown
    /// high bits from the remaining base, reporting
    /// [`Warning::UnsupportedSiderealProjection`]. Negative IDs are stored
    /// unchanged with no retained flags. Unknown base-mode fallback is deferred
    /// to evaluation; no astronomical value or callback is consumed here.
    ///
    /// `SIDM_USER` reads `t0` as TT Julian days, or UT1 Julian days if `mode`
    /// includes `SIDBIT_USER_UT`, and `ayan_t0` as degrees. These values are
    /// stored unchanged, including nonfinite values, and predefined modes
    /// generally ignore the pair. For example, `(SIDM_J2000, 0.0, 0.0)` selects
    /// a predefined frame; `(SIDM_USER, t0, ayan_t0)` supplies a defining pair.
    pub fn set_sidereal_mode(&mut self, mode: i32, t0: f64, ayan_t0: f64) -> Vec<Warning> {
        let (base, bits) = if mode >= 0 {
            (mode & !IMPLEMENTED_SIDBITS, mode & IMPLEMENTED_SIDBITS)
        } else {
            (mode, 0)
        };
        let mut warnings = Vec::new();
        let stored = if base >= 0 {
            let unknown = base & !BASE_MODE_MASK;
            if unknown != 0 {
                warnings.push(Warning::UnsupportedSiderealProjection { bits: unknown });
                base & BASE_MODE_MASK
            } else {
                base
            }
        } else {
            base
        };
        self.mode = stored;
        self.t0 = t0;
        self.ayan_t0 = ayan_t0;
        self.bits = bits;
        warnings
    }

    /// Store the implemented projection and epoch-timescale flags independently.
    ///
    /// This field update preserves the base mode and the raw defining epoch
    /// and angular-offset values, including negative modes and non-finite
    /// values. `bits` is an integer mask; bits outside
    /// [`IMPLEMENTED_SIDBITS`] are discarded without diagnostics. Callers own
    /// any staged warning delivery before assigning the remaining fields.
    pub fn set_retained_bits(&mut self, bits: i32) {
        self.bits = bits & IMPLEMENTED_SIDBITS;
    }

    /// Stored mode ID: the base ayanamsha mode after flag stripping, or a
    /// negative/unassigned ID unchanged.
    pub fn mode(self) -> i32 {
        self.mode
    }

    /// Stored `(base_mode, defining_epoch_jd, defining_ayanamsha_deg)` triple.
    /// The epoch is TT unless [`Self::bits`] contains `SIDBIT_USER_UT` (UT1).
    /// Projection bits are returned separately by [`Self::bits`].
    pub fn full(self) -> (i32, f64, f64) {
        (self.mode, self.t0, self.ayan_t0)
    }

    /// Retained implemented `SIDBIT` flags.
    pub fn bits(self) -> i32 {
        self.bits
    }
}

// Sidereal defining epochs, reference angles and frame projections.
//
// A defining epoch and ayanamsha identify the zero point. The Vondrák pole
// geometry carries that definition to the requested TT epoch. Projection
// flags choose the ecliptic of the defining epoch or the invariable plane;
// the actual frame matrices are shared with the position pipeline.

/// J2000.0 epoch in TT Julian Days.
pub const J2000_TT: f64 = 2451545.0;

/// Accumulated precession in ecliptic longitude from `t0_tt` to `jd_tt`,
/// in degrees, measured along the fixed mean ecliptic of `t0_tt`.
///
/// The sidereal reference direction is fixed on the **mean ecliptic of
/// `t0_tt`**; the returned value is the amount by which the vernal equinox
/// of `jd_tt` has regressed along that fixed ecliptic: the equinox of
/// `jd_tt` is projected into the ecliptic-of-`t0_tt` frame, giving its
/// longitude `lambda`, and the accumulation is `-lambda`. NOT reduced mod
/// 360 and non-additive in `t0_tt` (see the section docs). Both TT arguments
/// must be finite. Finite extrapolation may overflow the pole polynomials or
/// normalization; their IEEE components are consumed in the ordered dot
/// products and `atan2`, so an undefined projection returns NaN. This does
/// not extend the Vondrák model's physical ±200,000-year validity, and it does
/// not change the stricter checked pole/frame accessors.
pub fn method_b_accumulated_precession(jd_tt: f64, t0_tt: f64) -> Result<f64, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !t0_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: t0_tt });
    }
    let e0 = crate::erfa::ltp_ecliptic_frame_ieee(t0_tt)?;
    let equinox_t = crate::erfa::ltp_ecliptic_frame_ieee(jd_tt)?[3];
    let wx = e0[0][0] * equinox_t[0] + e0[0][1] * equinox_t[1] + e0[0][2] * equinox_t[2];
    let wy = e0[1][0] * equinox_t[0] + e0[1][1] * equinox_t[1] + e0[1][2] * equinox_t[2];
    let lam_deg = wy.atan2(wx).to_degrees();
    Ok(-lam_deg)
}

/// Longitude shift (deg) referring an epoch-pair mode to the mean ecliptic
/// of date. the model LEB fast path
/// (`_calc_ayanamsa_from_leb`) delegates to the same function so both
/// backends stay numerically identical; the Rust fast-path mapping lives
/// outside this band.
///
/// The fixed sidereal zero point `Z` sits on the mean ecliptic of `t0_tt`
/// at ecliptic longitude `defining_value_deg` (latitude 0). Returns the
/// difference between the tropical longitude of `Z` measured on the mean
/// ecliptic of date `tjd_tt` and on the mean ecliptic of `t0_tt` (the
/// Method-B baseline), unwrapped to `(-180, 180]`. Zero at
/// `tjd_tt == t0_tt`; a few milli-arcseconds for modern epochs, up to
/// ~0.5-0.7" for modes whose defining epoch is a millennium or more from
/// the request. Live star/galactic modes never route here (their value is
/// already an of-date longitude).
pub fn ecl_date_ayanamsha_delta_deg(
    defining_value_deg: f64,
    t0_tt: f64,
    tjd_tt: f64,
) -> Result<f64, Error> {
    if !defining_value_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: defining_value_deg,
        });
    }
    if !t0_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: t0_tt });
    }
    if !tjd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: tjd_tt });
    }
    // Sidereal zero point Z on the mean ecliptic of t 0, in the GCRS basis:
    // (ecliptic-of-t 0 basis rows)^T @ [cos, sin, 0]. Separate cos/sin calls
    // preserve the model evaluation order.
    let e0 = crate::erfa::ltp_ecliptic_frame(t0_tt)?;
    let lon0 = defining_value_deg.to_radians();
    let cx = lon0.cos();
    let sx = lon0.sin();
    let zx = e0[0][0] * cx + e0[1][0] * sx;
    let zy = e0[0][1] * cx + e0[1][1] * sx;
    let zz = e0[0][2] * cx + e0[1][2] * sx;
    // Tropical longitude of Z on the mean ecliptic and equinox of date.
    let ed = crate::erfa::ltp_ecliptic_frame(tjd_tt)?;
    let wx = ed[0][0] * zx + ed[0][1] * zy + ed[0][2] * zz;
    let wy = ed[1][0] * zx + ed[1][1] * zy + ed[1][2] * zz;
    let lon_date = wy.atan2(wx).to_degrees();
    let method_b = defining_value_deg + method_b_accumulated_precession(tjd_tt, t0_tt)?;
    Ok((lon_date - method_b + 180.0).rem_euclid(360.0) - 180.0)
}

/// Epoch-pair ayanamsha scalar at `jd_tt` (TT) in degrees, reduced to
/// `[0, 360)`.
///
/// `value + method_b_accumulated_precession(jd_tt, defining_epoch_tt)`,
/// plus the [`ecl_date_ayanamsha_delta_deg`] shift when `ecl_date` (the
/// `SIDBIT_ECL_DATE` bit) is set. `SIDM_USER` callers pass their stored
/// `(t0, ayan_t0)` pair (already on TT; the `SIDBIT_USER_UT` epoch still
/// needs its UT→TT conversion in the value layer); Valens callers pass `(VALENS_MOON_AYAN_T 0_DEG, defining TT
/// epoch)`. UT/TT conversion, mode dispatch, dynamic anchors and the
/// unknown-mode fallback are the value layer's job, not this scalar kernel's.
pub fn epoch_pair_ayanamsha_deg(
    jd_tt: f64,
    defining_value_deg: f64,
    defining_epoch_tt: f64,
    ecl_date: bool,
) -> Result<f64, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    if !defining_value_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: defining_value_deg,
        });
    }
    if !defining_epoch_tt.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: defining_epoch_tt,
        });
    }
    let mut value = defining_value_deg + method_b_accumulated_precession(jd_tt, defining_epoch_tt)?;
    if ecl_date {
        value += ecl_date_ayanamsha_delta_deg(defining_value_deg, defining_epoch_tt, jd_tt)?;
    }
    Ok(value.rem_euclid(360.0))
}

// Sidereal defining epochs, reference angles and frame projections.
//
// A defining epoch and ayanamsha identify the zero point. The Vondrák pole
// geometry carries that definition to the requested TT epoch. Projection
// flags choose the ecliptic of the defining epoch or the invariable plane;
// the actual frame matrices are shared with the position pipeline.

/// Mardyks ayanamsha at `jd_tt` (TT) in degrees, reduced to `[0, 360)`
///
/// Defining 30° at the September equinox 1998
/// ([`MARDYKS_DEFINING`]), propagated with the IAU 2006
/// general-precession angle difference. Operation order is preserved exactly
/// (`value + (P(t) − P(t0))`, then Euclidean modulo via `rem_euclid`).
pub fn mardyks_ayanamsha_deg(jd_tt: f64) -> Result<f64, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let (defining_value, defining_epoch) = MARDYKS_DEFINING;
    let now = crate::positions::general_precession_deg(jd_tt)?;
    let then = crate::positions::general_precession_deg(defining_epoch)?;
    Ok((defining_value + (now - then)).rem_euclid(360.0))
}

/// Nutation in longitude at `jd_tt` (TT) in degrees (the true-ayanamsha
/// increment of: IAU 2006/2000 A dpsi via the
/// native [`crate::erfa::nut06a`], converted with `to_degrees` like the
/// the model `math.degrees`).
pub fn nutation_in_longitude_deg(jd_tt: f64) -> Result<f64, Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let (dpsi_rad, _) = crate::erfa::nut06a(jd_tt)?;
    if !dpsi_rad.is_finite() {
        return Err(Error::NonFiniteArgument { value: dpsi_rad });
    }
    Ok(dpsi_rad.to_degrees())
}

/// True ayanamsha in degrees: the mean
/// value (already reduced to `[0, 360)` by the dispatch) plus the signed
/// nutation in longitude. The sum is deliberately NOT renormalized —
/// near a zero crossing it may read slightly negative or above 360,
/// matching the model convention.
pub fn true_ayanamsha_deg(mean_deg: f64, jd_tt: f64) -> Result<f64, Error> {
    if !mean_deg.is_finite() {
        return Err(Error::NonFiniteArgument { value: mean_deg });
    }
    Ok(mean_deg + nutation_in_longitude_deg(jd_tt)?)
}

/// Echoed return flag for `get_ayanamsa_ex[_ut]`.
///
/// Compatibility contract, preserved exactly: the echo carries only the resolved
/// ephemeris-selection bit plus (conditionally) `FLG_NONUT` — every other
/// input bit (`FLG_SPEED`, `FLG_SIDEREAL`, spurious bits) is dropped. The
/// ephemeris bit is single-winner with priority
/// `JPLEPH > SWIEPH > MOSEPH`, defaulting to `SWIEPH` when none is set
/// (`0 → 2`, `258 → 2`, `6 → 2`, `3 → 1`). `FLG_NONUT` is echoed only
/// when set AND the base mode is outside
/// [`AYANAMSA_EX_NONUT_DROP_MODES`].
pub fn ayanamsha_ex_retflag(flags: i32, base_mode: i32) -> i32 {
    let mut retflag = if flags & FLG_JPLEPH != 0 {
        FLG_JPLEPH
    } else if flags & FLG_SWIEPH != 0 {
        FLG_SWIEPH
    } else if flags & FLG_MOSEPH != 0 {
        FLG_MOSEPH
    } else {
        FLG_SWIEPH
    };
    if flags & FLG_NONUT != 0 && !AYANAMSA_EX_NONUT_DROP_MODES.contains(&base_mode) {
        retflag |= FLG_NONUT;
    }
    retflag
}

/// Flag resolution for the ayanamsha value.
///
/// Returns `(noaberr, nogdefl, use_true)`: a heliocentric or
/// barycentric request evaluates the live anchor without annual
/// aberration and without solar deflection (there is no Earth to be
/// aberrated by — the compatibility contract measured up to ~21" on the
/// aberrant modes); `FLG_TRUEPOS`/`FLG_NOABERR`/`FLG_NOGDEFL` set the
/// same toggles individually. `FLG_NONUT` or `FLG_J2000` selects the mean
/// ayanamsha (a J2000 ecliptic carries no nutation component, so the
/// true value would inject a spurious ~9–17" offset); otherwise the true
/// (mean + nutation) value applies.
pub fn resolve_ayanamsha_flags(iflag: i32) -> (bool, bool, bool) {
    let centered = iflag & (FLG_HELCTR | FLG_BARYCTR) != 0;
    let noaberr = iflag & (FLG_TRUEPOS | FLG_NOABERR) != 0 || centered;
    let nogdefl = iflag & (FLG_TRUEPOS | FLG_NOGDEFL) != 0 || centered;
    let use_true = iflag & FLG_NONUT == 0 && iflag & FLG_J2000 == 0;
    (noaberr, nogdefl, use_true)
}

/// Live anchor geometry a dynamic ayanamsha mode is evaluated from.
///
/// Variants name `STARS` entry (or derived direction) the
/// corresponding `_calc_ayanamsa` arm consumes; frame/formula modes
/// (`Mardyks`, `Valens`) never route here. Catalog propagation and the
/// apparent-place reduction belong to the fixed-star pipeline; this enum only fixes the
/// dispatch contract so the geometry can plug in without changing the
/// orchestration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DynamicAnchor {
    /// `STARS["SPICA"]` apparent mean-ecliptic longitude (True Citra).
    Spica,
    /// `STARS["REVATI"]` apparent mean-ecliptic longitude (True Revati).
    Revati,
    /// `STARS["PUSHYA"]` apparent mean-ecliptic longitude (True Pushya,
    /// True Sheoran — the targets differ, the anchor is shared).
    Pushya,
    /// `STARS["MULA"]` apparent mean-ecliptic longitude (True Mula).
    Mula,
    /// `STARS["ALDEBARAN"]` apparent mean-ecliptic longitude.
    Aldebaran,
    /// Sgr A* apparent mean-ecliptic longitude (`STARS["GAL_CENTER"]`).
    GalCenter,
    /// `STARS["GAL_NORTH_POLE"]` geometric (deflection/aberration-free)
    /// mean-ecliptic longitude (galactic-equator modes).
    GalNorthPole,
    /// Sgr A* hour-circle projection onto the mean ecliptic of date
    MulaWilhelm,
}

/// Provider of live anchor longitudes for the dynamic ayanamsha arms.
///
/// Returns the anchor's longitude on the mean ecliptic of date at
/// `jd_tt` (TT) in degrees; `noaberr`/`nogdefl` carry the
/// `FLG_NOABERR`/`FLG_NOGDEFL` anchor convention (geometric
/// realization for `FLG_TRUEPOS`, inert for the pole). The
/// `GalNorthPole` anchor is geometric by definition and ignores both
/// toggles. Returned longitudes need not be pre-reduced; the dispatch
/// applies the model target subtraction and `% 360`.
pub trait DynamicAnchorProvider {
    /// Whether every served anchor has complete owned source/clock kinematics.
    /// This read-free capability never infers purity from scalar callbacks.
    fn has_analytic_anchor_rates(&self) -> bool {
        false
    }

    /// Optional complete mean-anchor derivative in degrees per input TT day.
    /// The default declares opaque scalar providers; it performs no reads.
    fn anchor_rate_deg_day(
        &self,
        _anchor: DynamicAnchor,
        _jd_tt: f64,
        _noaberr: bool,
        _nogdefl: bool,
    ) -> Result<Option<f64>, Error> {
        Ok(None)
    }

    /// Anchor longitude in degrees (any finite value; the dispatch
    /// reduces mod 360 after the model target subtraction).
    fn anchor_longitude_deg(
        &self,
        anchor: DynamicAnchor,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<f64, Error>;
}

/// Mean ayanamsha at `jd_tt` (TT) in degrees, reduced to `[0, 360)`
///
/// - `base_mode` is the stripped base ID ([`strip_to_base_mode`]
///   contract: non-negative IDs mask to the low byte, negatives fall
///   through to the fallback).
/// - `user_pair` is the `(t0_tt, ayan_t0)` epoch pair, consumed only by
///   `SIDM_USER` (already on TT — the `SIDBIT_USER_UT` conversion lives
///   in the serving).
/// - `valens_epoch_tt` is the Valens defining TT epoch
///   (`VALENS_MOON_T0_UT + ΔT`, baked by the serving like the model).
/// - `ecl_date` is the `SIDBIT_ECL_DATE` bit: epoch-pair arms add the
///   ecliptic-of-date delta; live star/galactic arms are naturally inert (their
///   value is already an of-date longitude).
/// - `noaberr`/`nogdefl` are pre-gated to the aberrant modes by the
///   caller; epoch arms
///   ignore them.
/// - `dynamic` serves the 13 live star/galactic arms; `None` refuses
///   them with [`Error::AyanamshaAnchorUnserved`] (live-anchor geometry is
///   fixed-star-pipeline scope), never
///   a silent epoch substitution.
/// - Unknown IDs (negatives, unassigned positives) warn
///   ([`Warning::UnknownSiderealMode`]) and fall back to the
///   Fagan/Bradley defining pair — the model code contract, not the
///   stale docstring mention of Lahiri. The returned warning is `None`
///   for every recognized mode.
#[allow(clippy::too_many_arguments)] // positional arity mirrors the model dispatch
pub fn calc_mean_ayanamsha_deg(
    jd_tt: f64,
    base_mode: i32,
    user_pair: (f64, f64),
    valens_epoch_tt: f64,
    ecl_date: bool,
    noaberr: bool,
    nogdefl: bool,
    dynamic: Option<&dyn DynamicAnchorProvider>,
) -> Result<(f64, Option<Warning>), Error> {
    calc_mean_ayanamsha_with_epochs(
        jd_tt,
        base_mode,
        &|| Ok(user_pair),
        &|| Ok(valens_epoch_tt),
        ecl_date,
        noaberr,
        nogdefl,
        dynamic,
        None,
    )
}

/// Dispatch the mean angle while acquiring only the selected defining epoch.
/// `user_pair_at` returns a TT JD and degrees; `valens_epoch_at` returns a TT
/// JD. Each callback is consumed once in its owning mode, after the finite
/// evaluation-TT guard. Other modes consume neither callback. Formula, anchor,
/// warning and selected provider errors follow the scalar dispatcher.
#[allow(clippy::too_many_arguments)]
fn calc_mean_ayanamsha_with_epochs(
    jd_tt: f64,
    base_mode: i32,
    user_pair_at: &dyn Fn() -> Result<(f64, f64), Error>,
    valens_epoch_at: &dyn Fn() -> Result<f64, Error>,
    ecl_date: bool,
    noaberr: bool,
    nogdefl: bool,
    dynamic: Option<&dyn DynamicAnchorProvider>,
    unknown_mode_at: Option<&dyn Fn(i32) -> Result<(), Error>>,
) -> Result<(f64, Option<Warning>), Error> {
    if !jd_tt.is_finite() {
        return Err(Error::NonFiniteArgument { value: jd_tt });
    }
    let mode = strip_to_base_mode(base_mode);
    // The aberration/deflection-free anchor applies only to the live
    // star/galactic-center modes; epoch arms
    // ignore both toggles regardless of the request flags.
    let aberrant = ABERRANT_ANCHOR_MODES.contains(&mode);
    let eff_noaberr = noaberr && aberrant;
    let eff_nogdefl = nogdefl && aberrant;

    if mode == SIDM_USER {
        let (t0_tt, ayan_t0) = user_pair_at()?;
        return Ok((
            epoch_pair_ayanamsha_deg(jd_tt, ayan_t0, t0_tt, ecl_date)?,
            None,
        ));
    }

    if DYNAMIC_AYANAMSHA_MODES.contains(&mode) {
        // Mardyks and Valens are formula/epoch arms that need no catalog.
        if mode == SIDM_GALALIGN_MARDYKS {
            return Ok((mardyks_ayanamsha_deg(jd_tt)?, None));
        }
        if mode == SIDM_VALENS_MOON {
            let valens_epoch_tt = valens_epoch_at()?;
            let mut value =
                VALENS_MOON_AYAN_T0_DEG + method_b_accumulated_precession(jd_tt, valens_epoch_tt)?;
            if ecl_date {
                value +=
                    ecl_date_ayanamsha_delta_deg(VALENS_MOON_AYAN_T0_DEG, valens_epoch_tt, jd_tt)?;
            }
            return Ok((value.rem_euclid(360.0), None));
        }
        // The 13 live star/galactic arms need catalog geometry.
        let provider = dynamic.ok_or(Error::AyanamshaAnchorUnserved { mode })?;
        let anchor_lon = |anchor: DynamicAnchor| -> Result<f64, Error> {
            let lon = provider.anchor_longitude_deg(anchor, jd_tt, eff_noaberr, eff_nogdefl)?;
            if !lon.is_finite() {
                return Err(Error::NonFiniteArgument { value: lon });
            }
            Ok(lon)
        };
        let value = if mode == SIDM_TRUE_CITRA {
            anchor_lon(DynamicAnchor::Spica)? - 180.0
        } else if mode == SIDM_TRUE_REVATI {
            // Revati (zeta Piscium) at 29°50' Pisces: star longitude +
            // 10 arcmin, exactly as released.
            anchor_lon(DynamicAnchor::Revati)? + 10.0 / 60.0
        } else if mode == SIDM_TRUE_PUSHYA {
            // delta Cancri held at 16° Cancer (106° absolute), exactly as released.
            anchor_lon(DynamicAnchor::Pushya)? - 106.0
        } else if mode == SIDM_TRUE_MULA {
            // lambda Scorpii held at 0° Sagittarius (240°), exactly as released.
            anchor_lon(DynamicAnchor::Mula)? - 240.0
        } else if mode == SIDM_TRUE_SHEORAN {
            anchor_lon(DynamicAnchor::Pushya)? - SHEORAN_TARGET_LON
        } else if mode == SIDM_ALDEBARAN_15TAU {
            // Aldebaran held at 15° Taurus (45° absolute). This fixed
            // catalog anchor is absent from `_ABERRANT_ANCHOR_MODES`, so
            // the toggles above are already inert here by construction.
            anchor_lon(DynamicAnchor::Aldebaran)? - ALDEBARAN_TARGET_LON
        } else if mode == SIDM_GALCENT_0SAG
            || mode == SIDM_GALCENT_RGILBRAND
            || mode == SIDM_GALCENT_COCHRANE
        {
            let target = GALCENT_TARGET_LON
                .iter()
                .find(|(id, _)| *id == mode)
                .map(|(_, lon)| *lon)
                .ok_or(Error::AyanamshaAnchorUnserved { mode })?;
            anchor_lon(DynamicAnchor::GalCenter)? - target
        } else if mode == SIDM_GALEQU_IAU1958
            || mode == SIDM_GALEQU_TRUE
            || mode == SIDM_GALEQU_MULA
        {
            let target = GALEQU_TARGET_LON
                .iter()
                .find(|(id, _)| *id == mode)
                .map(|(_, lon)| *lon)
                .ok_or(Error::AyanamshaAnchorUnserved { mode })?;
            let pole = anchor_lon(DynamicAnchor::GalNorthPole)?;
            (pole + 90.0).rem_euclid(360.0) - target
        } else if mode == SIDM_GALCENT_MULA_WILHELM {
            anchor_lon(DynamicAnchor::MulaWilhelm)? - MULA_WILHELM_TARGET_LON
        } else {
            return Err(Error::AyanamshaAnchorUnserved { mode });
        };
        return Ok((value.rem_euclid(360.0), None));
    }

    // Registered defining pair, or the unknown-mode fallback (warn +
    // Fagan/Bradley, the model code contract).
    let mut warning = None;
    let (defining_value, defining_epoch) =
        match AYANAMSHA_DEFINING.iter().find(|(id, _, _)| *id == mode) {
            Some((_, value, epoch)) => (*value, *epoch),
            None => {
                if let Some(observe) = unknown_mode_at {
                    observe(mode)?;
                }
                warning = Some(Warning::UnknownSiderealMode { mode });
                let (id0, value0, epoch0) = AYANAMSHA_DEFINING[0];
                debug_assert_eq!(id0, SIDM_FAGAN_BRADLEY);
                (value0, epoch0)
            }
        };
    Ok((
        epoch_pair_ayanamsha_deg(jd_tt, defining_value, defining_epoch, ecl_date)?,
        warning,
    ))
}

/// A context-bound ayanamsha provider.
///
/// Stores a selection by value and borrows a `UT1 JD -> TT−UT1 days` callback.
/// Dynamic stellar/galactic modes additionally require [`DynamicAnchorProvider`].
/// Scalar ayanamshas do not apply defining-ecliptic or invariable-plane
/// projections; the coordinate serving layer owns those transforms. The adapter
/// performs no direct I/O or global mutation, but consumed callbacks may perform
/// I/O, fail or panic. All borrowed providers must outlive this serving value.
pub struct SealedAyanamshaServing<'a> {
    selection: SiderealSelection,
    deltat_days_at: &'a dyn Fn(f64) -> Result<f64, Error>,
    dynamic: Option<&'a dyn DynamicAnchorProvider>,
    unknown_mode_at: Option<&'a dyn Fn(i32) -> Result<(), Error>>,
}

impl<'a> SealedAyanamshaServing<'a> {
    /// Store the selection and borrow providers without invoking them.
    ///
    /// This infallible constructor validates neither epochs nor provider
    /// coverage. Missing anchors and invalid selected values fail only when
    /// their evaluation path is consumed.
    pub fn new(
        selection: SiderealSelection,
        deltat_days_at: &'a dyn Fn(f64) -> Result<f64, Error>,
        dynamic: Option<&'a dyn DynamicAnchorProvider>,
    ) -> Self {
        Self {
            selection,
            deltat_days_at,
            dynamic,
            unknown_mode_at: None,
        }
    }

    /// Observe an actually consumed unknown-mode fallback before its defining
    /// pair is evaluated. No callback runs for recognized or unselected modes.
    /// Observer errors propagate before fallback arithmetic and true nutation;
    /// returned diagnostics remain intact for ordinary Rust callers.
    pub fn with_unknown_mode_observer(
        mut self,
        observe: Option<&'a dyn Fn(i32) -> Result<(), Error>>,
    ) -> Self {
        self.unknown_mode_at = observe;
        self
    }

    /// Complete live mean rate from the selected physical anchor, when owned.
    pub(crate) fn live_mean_rate_deg_day(
        &self,
        tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<Option<f64>, Error> {
        let mode = strip_to_base_mode(self.selection.mode());
        let Some((_, key)) = DYNAMIC_MODE_ANCHORS.iter().find(|(id, _)| *id == mode) else {
            return Ok(None);
        };
        let anchor = if mode == SIDM_GALCENT_MULA_WILHELM {
            DynamicAnchor::MulaWilhelm
        } else {
            match *key {
                "SPICA" => DynamicAnchor::Spica,
                "REVATI" => DynamicAnchor::Revati,
                "PUSHYA" => DynamicAnchor::Pushya,
                "MULA" => DynamicAnchor::Mula,
                "ALDEBARAN" => DynamicAnchor::Aldebaran,
                "GAL_CENTER" => DynamicAnchor::GalCenter,
                "GAL_NORTH_POLE" => DynamicAnchor::GalNorthPole,
                _ => return Ok(None),
            }
        };
        let aberrant = ABERRANT_ANCHOR_MODES.contains(&mode);
        self.dynamic
            .map(|provider| {
                provider.anchor_rate_deg_day(anchor, tt, noaberr && aberrant, nogdefl && aberrant)
            })
            .transpose()
            .map(Option::flatten)
    }

    /// Bound sidereal selection (mode, epoch pair, retained bits).
    pub fn selection(&self) -> SiderealSelection {
        self.selection
    }

    /// Whether the `SIDBIT_ECL_DATE` projection bit is retained.
    fn ecl_date(&self) -> bool {
        self.selection.bits() & SIDBIT_ECL_DATE != 0
    }

    /// Whether the `SIDBIT_USER_UT` timescale bit is retained.
    fn user_ut(&self) -> bool {
        self.selection.bits() & SIDBIT_USER_UT != 0
    }

    /// TT epoch for a UT epoch through the bound ΔT provider
    /// (`jd_ut + ΔT(jd_ut)`, the calculation order).
    pub fn jd_tt_for_ut(&self, jd_ut: f64) -> Result<f64, Error> {
        if !jd_ut.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_ut });
        }
        let dt = (self.deltat_days_at)(jd_ut)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_ut + dt)
    }

    /// UT1 epoch for a TT epoch using the bound one-step inverse.
    /// The callback supplies TT-minus-UT1 days at the requested TT Julian day;
    /// errors propagate without substituting another clock.
    pub fn ut1_for_tt(&self, jd_tt: f64) -> Result<f64, Error> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let dt = (self.deltat_days_at)(jd_tt)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(jd_tt - dt)
    }

    /// Valens defining TT epoch (`VALENS_MOON_T0_UT + ΔT`, the model
    /// association inside `_calc_ayanamsa`).
    pub fn valens_epoch_tt(&self) -> Result<f64, Error> {
        let dt = (self.deltat_days_at)(VALENS_MOON_T0_UT)?;
        if !dt.is_finite() {
            return Err(Error::NonFiniteArgument { value: dt });
        }
        Ok(VALENS_MOON_T0_UT + dt)
    }

    /// `SIDM_USER` epoch pair on TT: with `SIDBIT_USER_UT` the stored
    /// `t0` is a UT date and moves through the ΔT provider; otherwise it is already TT.
    ///
    /// Crate-visible for the sealed projection serving, which anchors the
    /// `SIDBIT_ECL_T0` projection plane on the same TT epoch the
    /// ayanamsha dispatch uses.
    pub(crate) fn user_pair_tt(&self) -> Result<(f64, f64), Error> {
        let (_, t0, ayan_t0) = self.selection.full();
        if self.user_ut() {
            if !t0.is_finite() {
                return Err(Error::NonFiniteArgument { value: t0 });
            }
            let dt = (self.deltat_days_at)(t0)?;
            if !dt.is_finite() {
                return Err(Error::NonFiniteArgument { value: dt });
            }
            Ok((t0 + dt, ayan_t0))
        } else {
            Ok((t0, ayan_t0))
        }
    }

    /// Mean ayanamsha at a TT epoch with explicit anchor toggles
    /// (the TT-layer dispatch; UT conversion is the caller's job). Only the
    /// selected user-UT or Valens defining epoch consumes the Delta-T provider;
    /// classical and live-anchor modes do not acquire unused epoch clocks.
    pub fn mean_tt_with(
        &self,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<(f64, Option<Warning>), Error> {
        calc_mean_ayanamsha_with_epochs(
            jd_tt,
            self.selection.mode(),
            &|| self.user_pair_tt(),
            &|| self.valens_epoch_tt(),
            self.ecl_date(),
            noaberr,
            nogdefl,
            self.dynamic,
            self.unknown_mode_at,
        )
    }

    /// Mean ayanamsha at a UT epoch.
    pub fn mean_ut(&self, jd_ut: f64) -> Result<(f64, Option<Warning>), Error> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.mean_tt_with(jd_tt, false, false)
    }

    /// Mean ayanamsha at a TT epoch.
    pub fn mean_tt(&self, jd_tt: f64) -> Result<(f64, Option<Warning>), Error> {
        let jd_ut = self.ut1_for_tt(jd_tt)?;
        self.mean_ut(jd_ut)
    }

    /// True ayanamsha at a TT epoch with explicit anchor toggles
    pub fn true_tt_with(
        &self,
        jd_tt: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<(f64, Option<Warning>), Error> {
        let (mean, warning) = self.mean_tt_with(jd_tt, noaberr, nogdefl)?;
        Ok((true_ayanamsha_deg(mean, jd_tt)?, warning))
    }

    /// True ayanamsha at a UT epoch with explicit anchor toggles.
    pub fn true_ut_with(
        &self,
        jd_ut: f64,
        noaberr: bool,
        nogdefl: bool,
    ) -> Result<(f64, Option<Warning>), Error> {
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        self.true_tt_with(jd_tt, noaberr, nogdefl)
    }

    /// Flag-aware ayanamsha at a UT epoch: `FLG_NONUT`/`FLG_J2000` select the
    /// mean value, otherwise the true value; helio/barycentric centers
    /// evaluate the live anchor without aberration or deflection.
    pub fn for_flags_ut(&self, jd_ut: f64, iflag: i32) -> Result<(f64, Option<Warning>), Error> {
        let (noaberr, nogdefl, use_true) = resolve_ayanamsha_flags(iflag);
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        if use_true {
            self.true_tt_with(jd_tt, noaberr, nogdefl)
        } else {
            self.mean_tt_with(jd_tt, noaberr, nogdefl)
        }
    }

    /// Extended entry at a UT epoch:
    /// the echoed retflag plus the true value, or the mean value under
    /// `FLG_NONUT`. Anchor toggles come from the flags alone (no
    /// helio/bary centering — the contract, distinct from `for_flags_ut`).
    pub fn ex_ut(&self, jd_ut: f64, flags: i32) -> Result<(i32, f64, Option<Warning>), Error> {
        let base = strip_to_base_mode(self.selection.mode());
        let retflag = ayanamsha_ex_retflag(flags, base);
        let noaberr = flags & (FLG_TRUEPOS | FLG_NOABERR) != 0;
        let nogdefl = flags & (FLG_TRUEPOS | FLG_NOGDEFL) != 0;
        let jd_tt = self.jd_tt_for_ut(jd_ut)?;
        let (mean, warning) = self.mean_tt_with(jd_tt, noaberr, nogdefl)?;
        let value = if flags & FLG_NONUT != 0 {
            mean
        } else {
            true_ayanamsha_deg(mean, jd_tt)?
        };
        Ok((retflag, value, warning))
    }

    /// Extended entry at a TT epoch:
    /// TT→UT1 through the native inversion, then the UT extended entry.
    pub fn ex_tt(&self, jd_tt: f64, flags: i32) -> Result<(i32, f64, Option<Warning>), Error> {
        let jd_ut = self.ut1_for_tt(jd_tt)?;
        self.ex_ut(jd_ut, flags)
    }
}

// Scalar ayanamsha evaluation.
//
// Epoch-defined modes use their defining TT epochs and angles. Mardyks uses
// the IAU 2006 general-precession angle; true ayanamsha adds IAU 2006/2000 A
// nutation in longitude without renormalizing the final sum. Stellar and
// galactic definitions require an explicit DynamicAnchorProvider; a missing
// anchor returns AyanamshaAnchorUnserved rather than an epoch substitution.
// TT/UT1 conversion uses the supplied Delta-T model.

/// Central-difference half-step (days) for the observer-frame ayanamsha
/// drift removal.
pub const SIDEREAL_SPEED_HALF_STEP_DAYS: f64 = 1.0 / 86400.0;

/// Subtract the ayanamsha from an ecliptic longitude.
///
/// A non-finite `lon_deg` (an undefined lunar point) passes through
/// unchanged, using the same physical model propagation; a non-finite ayanamsha
/// is a provider fault and is refused.
pub fn apply_sidereal_offset_deg(lon_deg: f64, ayanamsha_deg: f64) -> Result<f64, Error> {
    if !ayanamsha_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: ayanamsha_deg,
        });
    }
    if !lon_deg.is_finite() {
        return Ok(lon_deg);
    }
    Ok((lon_deg - ayanamsha_deg).rem_euclid(360.0))
}

/// Ayanamsha-drift rate in deg/day from two 1-second samples (speed half and the fast-path star-based
/// branch: the mod-360 samples are unwrapped into `(-180, 180]` before
/// dividing, so a 0/360 straddle cannot spike the speed).
pub fn sidereal_rate_correction_deg_day(
    ayan_prev_deg: f64,
    ayan_next_deg: f64,
) -> Result<f64, Error> {
    if !ayan_prev_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: ayan_prev_deg,
        });
    }
    if !ayan_next_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: ayan_next_deg,
        });
    }
    let delta = (ayan_next_deg - ayan_prev_deg + 180.0).rem_euclid(360.0) - 180.0;
    Ok(delta / (2.0 * SIDEREAL_SPEED_HALF_STEP_DAYS))
}

/// observer-frame sidereal correction of one ecliptic longitude and its
/// rate: offset
/// subtraction plus, under `with_speed` (the `FLG_SPEED` bit), removal
/// of the sampled ayanamsha drift. Non-finite longitudes pass through
/// with their rate untouched. `ayan_prev_deg`/`ayan_next_deg` must be the
/// dispatch values at ∓/± `SIDEREAL_SPEED_HALF_STEP_DAYS` around the same
/// epoch as the offset value; the factored signature takes precomputed
/// samples, so the caller owns the association. The finite-`lon_deg` /
/// finite-`dlon` gate on the speed leg differs from the model, which
/// subtracts drift from the rate even when `lon` is non-finite (Rust
/// leaves the rate untouched).
pub fn apply_sidereal_correction(
    lon_deg: f64,
    dlon_deg_day: f64,
    ayanamsha_deg: f64,
    ayan_prev_deg: f64,
    ayan_next_deg: f64,
    with_speed: bool,
) -> Result<(f64, f64), Error> {
    let lon = apply_sidereal_offset_deg(lon_deg, ayanamsha_deg)?;
    if with_speed && lon_deg.is_finite() && dlon_deg_day.is_finite() {
        let drift = sidereal_rate_correction_deg_day(ayan_prev_deg, ayan_next_deg)?;
        Ok((lon, dlon_deg_day - drift))
    } else {
        Ok((lon, dlon_deg_day))
    }
}

/// Degenerate-origin sentinel: the exact all-zero
/// 6-tuple is the coordinate origin, not a direction — every rotation
/// fixes it and no ayanamsha is ever subtracted from it. `-0.0` counts
/// as zero and NaN is non-degenerate (`NaN == 0.0` is false), matching
/// the zero-distance guard; NaN remains an undefined nonzero input.
pub fn is_degenerate_origin_6tuple(xx: &[f64; 6]) -> bool {
    xx.iter().all(|&v| v == 0.0)
}

/// True when a `FLG_SIDEREAL` request runs under a fixed-epoch mode
/// : a frame request,
/// not an ayanamsha offset. `base_mode` is the stripped base ID.
pub fn is_fixed_epoch_request(iflag: i32, base_mode: i32) -> bool {
    (iflag & FLG_SIDEREAL) != 0 && is_fixed_epoch_mode(base_mode)
}

/// Rewrite a fixed-epoch sidereal request to its equivalent frame
/// request:
/// strip the sidereal projection and the radians shaping, force the
/// J2000 mean frame without nutation.
pub fn fixed_epoch_request_flags(iflag: i32) -> i32 {
    (iflag & !(FLG_SIDEREAL | FLG_RADIANS)) | FLG_J2000 | FLG_NONUT
}

/// Echo flags for a fixed-epoch sidereal call: the rewritten sub-request's
/// retflag (`sub_retflag`) with its backend-added bits kept (ephemeris
/// echo bits, implied bits) minus the internal `FLG_J2000` rewrite, plus
/// the caller's flags and `FLG_NONUT` (the reference frame is mean); the
/// internal `FLG_J2000` rewrite is not exposed unless the caller passed
/// it explicitly.
pub fn fixed_epoch_retflag(sub_retflag: i32, iflag: i32) -> i32 {
    (sub_retflag & !FLG_J2000) | (iflag & (FLG_SIDEREAL | FLG_RADIANS | FLG_J2000)) | FLG_NONUT
}

/// Reference epoch (TT JD) of a fixed-epoch sidereal mode.
///
/// J2000 is the literal TT value; J 1900/B 1950 go through the native
/// `epj2jd`/`epb2jd` epoch converters (native fused form); Mardyks
/// uses its defining September-equinox-1998 epoch. `None` for any other
/// mode (the caller routes only fixed-epoch requests here).
pub fn fixed_epoch_t0_jd(mode: i32) -> Option<f64> {
    match mode {
        SIDM_J2000 => Some(J2000_TT),
        SIDM_J1900 => {
            let (d1, d2) = crate::erfa::epj2jd(1900.0).ok()?;
            Some(d1 + d2)
        }
        SIDM_B1950 => {
            let (d1, d2) = crate::erfa::epb2jd(1950.0).ok()?;
            Some(d1 + d2)
        }
        SIDM_GALALIGN_MARDYKS => Some(MARDYKS_DEFINING.1),
        _ => None,
    }
}

/// Passive x-axis rotation.
fn rot_x_passive(angle_rad: f64) -> [[f64; 3]; 3] {
    let (c, s) = (angle_rad.cos(), angle_rad.sin());
    [[1.0, 0.0, 0.0], [0.0, c, s], [0.0, -s, c]]
}

/// Passive z-axis rotation.
fn rot_z_passive(angle_rad: f64) -> [[f64; 3]; 3] {
    let (c, s) = (angle_rad.cos(), angle_rad.sin());
    [[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]]
}

/// Transpose of a 3×3 direction-cosine matrix (exact index permutation;
/// an involution that preserves orthogonality), factored once for every
/// `V(t0)·V(J2000)ᵀ` chain in this module where the model transposes inline.
fn transpose3(m: [[f64; 3]; 3]) -> [[f64; 3]; 3] {
    [
        [m[0][0], m[1][0], m[2][0]],
        [m[0][1], m[1][1], m[2][1]],
        [m[0][2], m[1][2], m[2][2]],
    ]
}

/// Paired t 0-frame rotations of [`fixed_epoch_matrices`].
pub type EpochMatrixPair = ([[f64; 3]; 3], [[f64; 3]; 3]);

/// `(M_eq, M_ecl)` rotating J2000 frames to the t 0 frames of a
/// fixed-epoch mode.
///
/// `M_eq`: mean equatorial J2000 → mean equatorial of t 0 (Vondrak 2011
/// chain; the ICRS frame bias cancels in the product, so the biased
/// realization is used on both sides like the model default).
/// `M_ecl`: mean ecliptic J2000 → mean ecliptic of t 0 (IAU 2006 mean
/// obliquity sandwich), additionally rotated by the defining 30° for
/// Mardyks. Identity pair when `t0 == J2000_TT` and `lon_offset == 0.0`
///
/// The IAU 2006 obliquity is evaluated directly in radians with fused Horner
/// accumulation, without a degree-polynomial roundtrip. Frame composition
/// uses ordinary binary64 matrix multiplication. The precession product uses
/// the transposed J2000 frame; subsequent obliquity rotations compose the
/// defining equatorial and ecliptic coordinate planes.
///
/// # Errors
///
/// Returns [`Error::AyanamshaAnchorUnserved`] when the mode has no served
/// defining epoch (the `None` arm), where the model raises the
/// partial-lookup miss. Clock/polynomial and matrix evaluation errors propagate.
pub fn fixed_epoch_matrices(mode: i32) -> Result<EpochMatrixPair, Error> {
    let t0 = fixed_epoch_t0_jd(mode).ok_or(Error::AyanamshaAnchorUnserved { mode })?;
    let lon_offset = fixed_epoch_lon_offset_deg(mode);
    if t0 == J2000_TT && lon_offset == 0.0 {
        let eye = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        return Ok((eye, eye));
    }
    let v_t0 = crate::positions::vondrak_precession_matrix(t0, true)?;
    let v_j2k = crate::positions::vondrak_precession_matrix(J2000_TT, true)?;
    let m_eq = sidereal_frame_product_transposed_right(v_t0, v_j2k)?;
    let eps_j2k = crate::erfa::pfw06(J2000_TT)?.3;
    let eps_t0 = crate::erfa::pfw06(t0)?.3;
    let mut m_ecl = sidereal_frame_product(
        sidereal_frame_product(rot_x_passive(eps_t0), m_eq)?,
        rot_x_passive(-eps_j2k),
    )?;
    if lon_offset != 0.0 {
        m_ecl = sidereal_frame_product(rot_z_passive(lon_offset.to_radians()), m_ecl)?;
    }
    Ok((m_eq, m_ecl))
}

/// Compose finite row-major 3×3 frame matrices, `C = A B`.
///
/// The shared scalar matrix kernel evaluates all products in binary64.
/// Non-finite inputs are refused; finite arithmetic follows IEEE overflow.
pub(crate) fn sidereal_frame_product(
    a: [[f64; 3]; 3],
    b: [[f64; 3]; 3],
) -> Result<[[f64; 3]; 3], Error> {
    crate::positions::mat3_mul(a, b)
}

/// Compose finite row-major frames as `C = A Bᵀ`.
/// Transposition changes indices without changing coefficients. Units and
/// finite-domain checks match [`sidereal_frame_product`].
pub(crate) fn sidereal_frame_product_transposed_right(
    a: [[f64; 3]; 3],
    b: [[f64; 3]; 3],
) -> Result<[[f64; 3]; 3], Error> {
    crate::positions::mat3_mul(a, transpose3(b))
}

/// Rotate a spherical `(lon, lat, dist, dlon, dlat, ddist)` state by a
/// time-independent matrix (same-level counterpart
/// `sidereal_epoch._rotate_spherical`).
///
/// Position and velocity rotate with the same matrix; the spherical
/// rates are re-derived from the rotated Cartesian pair. The
/// matrix-vector products and the radial-rate dot use binary64 arithmetic.
/// The
/// degenerate-origin sentinel passes through untouched (rotating it
/// through a unit-radius substitution manufactured a longitude out of
/// nothing in the model history). A non-finite longitude or latitude
/// passes through untouched — no frame conversion is defined for it
/// (the `_maybe_equatorial_convert` NaN contract), where the model
/// would propagate NaN through trigonometry. A zero distance keeps the direction
/// meaningful through a unit-radius substitution, restored to zero
/// afterwards; fully-static zero states keep zero rates.
pub fn rotate_spherical_6tuple(xx: [f64; 6], m: [[f64; 3]; 3]) -> Result<[f64; 6], Error> {
    if is_degenerate_origin_6tuple(&xx) {
        return Ok([0.0, 0.0, 0.0, 0.0, 0.0, 0.0]);
    }
    let [lon, lat, r, dlon, dlat, dr] = xx;
    if !lon.is_finite() || !lat.is_finite() {
        // Undefined lunar point (the `_maybe_equatorial_convert` NaN
        // contract): no frame conversion is defined, pass through.
        return Ok(xx);
    }
    let lon_r = lon.to_radians();
    let lat_r = lat.to_radians();
    let dlon_r = dlon.to_radians();
    let dlat_r = dlat.to_radians();
    let (cl, sl) = (lon_r.cos(), lon_r.sin());
    let (cb, sb) = (lat_r.cos(), lat_r.sin());
    let (r_eff, dr_eff) = if r == 0.0 { (1.0, 0.0) } else { (r, dr) };
    let pos = [r_eff * cb * cl, r_eff * cb * sl, r_eff * sb];
    let vel = [
        dr_eff * cb * cl - r_eff * sb * cl * dlat_r - r_eff * cb * sl * dlon_r,
        dr_eff * cb * sl - r_eff * sb * sl * dlat_r + r_eff * cb * cl * dlon_r,
        dr_eff * sb + r_eff * cb * dlat_r,
    ];
    let p = crate::positions::mat3_vec3(m, pos)?;
    let v = crate::positions::mat3_vec3(m, vel)?;
    let rho2 = p[0] * p[0] + p[1] * p[1];
    let rho = rho2.sqrt();
    let rn = (rho2 + p[2] * p[2]).sqrt();
    let lon2 = (p[1].atan2(p[0])).rem_euclid(2.0 * std::f64::consts::PI);
    let lat2 = (p[2] / rn).clamp(-1.0, 1.0).asin();
    let mut dlon2 = if rho2 > 0.0 {
        (p[0] * v[1] - p[1] * v[0]) / rho2
    } else {
        0.0
    };
    let mut dr2 = (0.0 + p[0] * v[0] + p[1] * v[1] + p[2] * v[2]) / rn;
    let mut dlat2 = if rho > 0.0 {
        (v[2] - dr2 * p[2] / rn) / rho
    } else {
        0.0
    };
    let mut rn_out = rn;
    if r == 0.0 {
        rn_out = 0.0;
        dr2 = 0.0;
        if dlon == 0.0 && dlat == 0.0 && dr == 0.0 {
            dlon2 = 0.0;
            dlat2 = 0.0;
        }
    }
    Ok([
        lon2.to_degrees(),
        lat2.to_degrees(),
        rn_out,
        dlon2.to_degrees(),
        dlat2.to_degrees(),
        dr2,
    ])
}

/// J2000 mean ecliptic → solar-system invariable-plane frame: node onto the x-axis,
/// then tilt by the inclination.
pub fn invariable_plane_matrix() -> [[f64; 3]; 3] {
    // Closed-form composition of two passive rotations over published
    // constants; both factors are finite by construction.
    let node = rot_z_passive(INVARIABLE_PLANE_NODE_DEG.to_radians());
    let incl = rot_x_passive(INVARIABLE_PLANE_INCL_DEG.to_radians());
    // `mat3_mul` refuses non-finite rows only; rotation rows are finite,
    // so the `Rx @ Rz` order holds without a fallible seam.
    let mut out = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = incl[i][0] * node[0][j] + incl[i][1] * node[1][j] + incl[i][2] * node[2][j];
        }
    }
    out
}

/// J2000 mean ecliptic → mean ecliptic and equinox of `t0_jd`.
///
/// Built from the Vondrak 2011 long-term chain with the Vondrak mean
/// obliquity at each epoch (frame self-consistency for ancient ECL_T 0
/// epochs far outside the IAU 2006 optimum) — deliberately NOT the
/// IAU 2006 `obl06` sandwich the fixed-epoch matrices use. The matrices are
/// composed in binary64 with a shared scalar kernel.
///
/// Formal provenance: Vondrák et al. 2011 (A&A 534, A 22), valid over
/// ±200 millennia against the IAU 2006 1800–2200 optimum; the ICRS
/// frame bias cancels in the `v_t0 @ v_j2kᵀ` product (the same
/// realization as the other epoch rotations).
pub fn ecliptic_of_t0_matrix(t0_jd: f64) -> Result<[[f64; 3]; 3], Error> {
    if !t0_jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: t0_jd });
    }
    let v_t0 = crate::positions::vondrak_precession_matrix(t0_jd, true)?;
    let v_j2k = crate::positions::vondrak_precession_matrix(J2000_TT, true)?;
    let m_eq = sidereal_frame_product_transposed_right(v_t0, v_j2k)?;
    let eps0 = crate::positions::vondrak_mean_obliquity_deg(J2000_TT)?.to_radians();
    let epst = crate::positions::vondrak_mean_obliquity_deg(t0_jd)?.to_radians();
    sidereal_frame_product(
        sidereal_frame_product(rot_x_passive(epst), m_eq)?,
        rot_x_passive(-eps0),
    )
}

/// Ecliptic-frame rotation for a SIDBIT projection, or `None`.
///
/// Maps a J2000 mean-ecliptic Cartesian state to the projection frame
/// with the sidereal zero point rotated onto the x-axis, so reading
/// spherical longitude yields the sidereal longitude directly.
/// `zero_point_deg` is the zero-point longitude *in the projection
/// frame*. `SIDBIT_ECL_T0` takes precedence with both bits set;
/// `SIDBIT_SSY_PLANE` alone selects the invariable plane. The final zero-point
/// rotation uses the shared binary64 matrix product.
pub fn sidbit_ecliptic_matrix(
    sid_bits: i32,
    t0_jd: f64,
    zero_point_deg: f64,
) -> Result<Option<[[f64; 3]; 3]>, Error> {
    if !t0_jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: t0_jd });
    }
    if !zero_point_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: zero_point_deg,
        });
    }
    let base = if (sid_bits & SIDBIT_ECL_T0) != 0 {
        ecliptic_of_t0_matrix(t0_jd)?
    } else if (sid_bits & SIDBIT_SSY_PLANE) != 0 {
        invariable_plane_matrix()
    } else {
        return Ok(None);
    };
    sidereal_frame_product(rot_z_passive(zero_point_deg.to_radians()), base).map(Some)
}

/// Invariable-plane longitude of the J2000 ayanamsha zero-point
/// direction: the sidereal zero point at J2000 mean-ecliptic longitude
/// `ayan_j2000_deg` projected onto the invariable plane. The finite degree
/// angle forms the unit ray `(cos a, sin a, 0)`; the fixed plane matrix and
/// `atan2` yield longitude in `[0, 360)`. Nonfinite arguments refuse before
/// trigonometry. This geometric projection does not establish a model-valid
/// epoch or acquire a clock or body provider.
pub fn ssy_plane_zero_point_deg(ayan_j2000_deg: f64) -> Result<f64, Error> {
    if !ayan_j2000_deg.is_finite() {
        return Err(Error::NonFiniteArgument {
            value: ayan_j2000_deg,
        });
    }
    ssy_plane_zero_point_deg_ieee(ayan_j2000_deg)
}

/// Project an internally computed J2000 mean-ecliptic zero-point longitude
/// into the invariable plane. Angles are degrees; the returned longitude is
/// in `[0, 360)` when defined. The fixed dimensionless matrix acts on the
/// unit ray `(cos a, sin a, 0)` with row sums starting at positive zero.
/// A derived NaN passes through trigonometry, row products and `atan2`;
/// it is not interpreted as a provider/configuration failure. Infinite angles
/// refuse before trigonometry. The checked public entry continues to reject
/// every nonfinite input, including NaN. No caller epoch is inferred here.
pub(crate) fn ssy_plane_zero_point_deg_ieee(ayan_j2000_deg: f64) -> Result<f64, Error> {
    if ayan_j2000_deg.is_infinite() {
        return Err(Error::MathDomain {
            value: ayan_j2000_deg,
        });
    }
    let m = invariable_plane_matrix();
    let a = ayan_j2000_deg.to_radians();
    let v = [a.cos(), a.sin(), 0.0];
    let p = m.map(|row| 0.0 + row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
    Ok(p[1].atan2(p[0]).to_degrees().rem_euclid(360.0))
}

/// Project a J2000|NONUT *ecliptic* result onto a SIDBIT frame.
///
/// Only ecliptic output is routed here (spherical or `FLG_XYZ`);
/// equatorial requests take the T 0-mean-equator arm instead. `m_ecl`
/// already folds in the zero-point rotation. `FLG_RADIANS` is re-applied
/// for spherical output only.
pub fn transform_sidbit_result(
    xx: [f64; 6],
    iflag: i32,
    m_ecl: [[f64; 3]; 3],
) -> Result<[f64; 6], Error> {
    if (iflag & FLG_XYZ) != 0 {
        let pos = crate::positions::mat3_vec3(m_ecl, [xx[0], xx[1], xx[2]])?;
        let vel = crate::positions::mat3_vec3(m_ecl, [xx[3], xx[4], xx[5]])?;
        return Ok([pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]]);
    }
    let mut out = rotate_spherical_6tuple(xx, m_ecl)?;
    if (iflag & FLG_RADIANS) != 0 {
        out[0] = out[0].to_radians();
        out[1] = out[1].to_radians();
        out[3] = out[3].to_radians();
        out[4] = out[4].to_radians();
    }
    Ok(out)
}

/// J2000 mean equatorial → mean equatorial and equinox of `t0_jd`
/// : the fixed-epoch
/// chain for an arbitrary epoch; identity at J2000. The defining frame and
/// transposed J2000 frame are composed by the shared binary64 matrix kernel.
pub fn equatorial_epoch_matrix(t0_jd: f64) -> Result<[[f64; 3]; 3], Error> {
    if !t0_jd.is_finite() {
        return Err(Error::NonFiniteArgument { value: t0_jd });
    }
    if t0_jd == J2000_TT {
        return Ok([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]);
    }
    let v_t0 = crate::positions::vondrak_precession_matrix(t0_jd, true)?;
    let v_j2k = crate::positions::vondrak_precession_matrix(J2000_TT, true)?;
    sidereal_frame_product_transposed_right(v_t0, v_j2k)
}

/// Map a FLG_J2000|FLG_NONUT *equatorial* result onto the mean equator
/// and equinox of `t0_jd`.
///
/// The `SIDBIT_ECL_T0` equatorial arm: no ayanamsha subtraction —
/// position and velocity in the t 0 mean frame. Handles spherical and
/// `FLG_XYZ`, re-applying `FLG_RADIANS` for spherical output.
pub fn transform_equatorial_epoch_result(
    xx: [f64; 6],
    iflag: i32,
    t0_jd: f64,
) -> Result<[f64; 6], Error> {
    let m_eq = equatorial_epoch_matrix(t0_jd)?;
    if (iflag & FLG_XYZ) != 0 {
        let pos = crate::positions::mat3_vec3(m_eq, [xx[0], xx[1], xx[2]])?;
        let vel = crate::positions::mat3_vec3(m_eq, [xx[3], xx[4], xx[5]])?;
        return Ok([pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]]);
    }
    let mut out = rotate_spherical_6tuple(xx, m_eq)?;
    if (iflag & FLG_RADIANS) != 0 {
        out[0] = out[0].to_radians();
        out[1] = out[1].to_radians();
        out[3] = out[3].to_radians();
        out[4] = out[4].to_radians();
    }
    Ok(out)
}

/// Map a FLG_J2000|FLG_NONUT result onto the t 0 frame of a fixed-epoch
/// mode.
///
/// Ecliptic-spherical, equatorial-spherical and XYZ representations;
/// `SIDM_J2000` is the identity (still re-applying `FLG_RADIANS` for
/// spherical output). XYZ selects the equatorial or ecliptic matrix by
/// `FLG_EQUATORIAL`, like the spherical branch.
pub fn transform_fixed_epoch_result(
    xx: [f64; 6],
    iflag: i32,
    sid_mode: i32,
) -> Result<[f64; 6], Error> {
    let mut out = xx;
    if sid_mode != SIDM_J2000 {
        let (m_eq, m_ecl) = fixed_epoch_matrices(sid_mode)?;
        let m = if (iflag & FLG_EQUATORIAL) != 0 {
            m_eq
        } else {
            m_ecl
        };
        if (iflag & FLG_XYZ) != 0 {
            let pos = crate::positions::mat3_vec3(m, [xx[0], xx[1], xx[2]])?;
            let vel = crate::positions::mat3_vec3(m, [xx[3], xx[4], xx[5]])?;
            out = [pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]];
        } else {
            out = rotate_spherical_6tuple(xx, m)?;
        }
    }
    if (iflag & FLG_RADIANS) != 0 && (iflag & FLG_XYZ) == 0 {
        out[0] = out[0].to_radians();
        out[1] = out[1].to_radians();
        out[3] = out[3].to_radians();
        out[4] = out[4].to_radians();
    }
    Ok(out)
}

/// Star/galactic "true" modes whose zero point is defined by a live
/// ecliptic-of-date direction, for which the SIDBIT_ECL_T 0/SSY_PLANE
/// re-projection is inert.
///
/// The projection is a no-op for exactly these modes (the un-projected
/// baseline is reproduced); epoch-anchored mean modes, Aldebaran = 15
/// Tau (an absolute-longitude anchor) and Valens (a Method-B epoch
/// mode) do accept the projection. Mardyks reaches the same result
/// through the fixed-epoch frame path and needs no listing. This is the
/// same "calculated" star/galactic class that drops `FLG_NONUT` from the
/// extended entry — see [`AYANAMSA_EX_NONUT_DROP_MODES`].
pub const SIDBIT_PROJECTION_SUPPRESS_MODES: [i32; 12] = [
    SIDM_GALCENT_0SAG,
    SIDM_TRUE_CITRA,
    SIDM_TRUE_REVATI,
    SIDM_TRUE_PUSHYA,
    SIDM_GALCENT_RGILBRAND,
    SIDM_GALEQU_IAU1958,
    SIDM_GALEQU_TRUE,
    SIDM_GALEQU_MULA,
    SIDM_TRUE_MULA,
    SIDM_GALCENT_MULA_WILHELM,
    SIDM_TRUE_SHEORAN,
    SIDM_GALCENT_COCHRANE,
];

/// Classical defining epochs (TT JD) for the ECL_T 0 projection plane
///
/// Where a mode's classical defining epoch — the mean ecliptic its
/// authors referred the sidereal zero to — differs from the VALUE
/// anchor epoch in [`AYANAMSHA_DEFINING`], the projection plane follows
/// the classical epoch. Entries without a row here resolve from the
/// defining table (or J2000 when unlisted there).
///
/// Entry provenance (published defining statement + JD rationale):
/// Lahiri — Calendar Reform Committee 1955; Kugler triplet —
/// Huber 1958 / Fagan 1950; Hipparchos — Mercier; De Luce —
/// De Luce 1963; Djwhal Khul — Lindsay 2006 (1900.0 standard epoch);
/// Britton — Britton 2010; Lahiri-1940 — Lahiri 1939/1940.
pub const ECL_T0_CLASSICAL_EPOCHS: [(i32, f64); 9] = [
    (SIDM_LAHIRI, 2435553.5),
    (SIDM_BABYL_KUGLER1, 1684532.5),
    (SIDM_BABYL_KUGLER2, 1684532.5),
    (SIDM_BABYL_KUGLER3, 1684532.5),
    (SIDM_HIPPARCHOS, 1674484.0),
    (SIDM_DELUCE, 1721057.5),
    (SIDM_DJWHAL_KHUL, 2415020.0),
    (SIDM_BABYL_BRITTON, 1721057.5),
    (SIDM_LAHIRI_1940, 2415020.0),
];

/// Defining epoch (JD TT) of a sidereal mode for the ECL_T 0 projection
///
/// `SIDM_USER` takes the caller's stored epoch on TT (`user_t0_tt`,
/// already through the `SIDBIT_USER_UT` conversion when set —
/// mirroring the `_calc_ayanamsa` anchoring convention); Valens takes
/// its UT-anchored defining epoch through the same ΔT step
/// (`valens_epoch_tt`); classical modes take the table above; the rest
/// resolve from [`AYANAMSHA_DEFINING`], defaulting to J2000. The plane
/// is independent of the stored ayanamsha value (rotating the zero
/// point within the plane moves longitudes, never latitudes). The plane
/// is the stored t 0 taken literally (`t0 = 0.0` means JD 0): the
/// latitude shift vanishes when t 0 equals the computation date and
/// grows ~0.36″ per year of |date − t 0| — the Valens special case
/// exists because silent J2000 fallback cost up to ~4′ at high
/// ecliptic latitudes.
pub fn ecl_t0_epoch_jd(mode: i32, user_t0_tt: f64, valens_epoch_tt: f64) -> f64 {
    if mode == SIDM_USER {
        return user_t0_tt;
    }
    if mode == SIDM_VALENS_MOON {
        return valens_epoch_tt;
    }
    if let Some((_, epoch)) = ECL_T0_CLASSICAL_EPOCHS.iter().find(|(id, _)| *id == mode) {
        return *epoch;
    }
    AYANAMSHA_DEFINING
        .iter()
        .find(|(id, _, _)| *id == mode)
        .map(|(_, _, epoch)| *epoch)
        .unwrap_or(J2000_TT)
}

/// Sidereal position route for one request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SiderealRoute {
    /// Plain ecliptic correction: subtract the flag-aware ayanamsha
    /// (mean or true) with the sealed speed-drift rule.
    BaseCorrection,
    /// Fixed-epoch frame request (`SIDM_J2000/J1900/B1950/Mardyks`):
    /// rewrite to J2000|NONUT, sub-request, rotate to the t 0 frame.
    FixedEpoch,
    /// `SIDBIT_ECL_T0`/`SIDBIT_SSY_PLANE` ecliptic projection onto the
    /// plane with the zero point rotated out.
    SidbitEcliptic,
    /// `SIDBIT_ECL_T0` equatorial arm: the t 0 mean equator/equinox with
    /// no ayanamsha subtraction.
    SidbitEquatorialT0,
}

/// Resolve the sidereal route, or `None` for a tropical request
pub fn resolve_sidereal_route(iflag: i32, base_mode: i32, sid_bits: i32) -> Option<SiderealRoute> {
    if (iflag & FLG_SIDEREAL) == 0 {
        return None;
    }
    if is_fixed_epoch_request(iflag, base_mode) {
        return Some(SiderealRoute::FixedEpoch);
    }
    if (sid_bits & (SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE)) != 0
        && ((iflag & FLG_EQUATORIAL) == 0 || (sid_bits & SIDBIT_ECL_T0) != 0)
        && !SIDBIT_PROJECTION_SUPPRESS_MODES.contains(&base_mode)
    {
        if (iflag & FLG_EQUATORIAL) != 0 {
            return Some(SiderealRoute::SidbitEquatorialT0);
        }
        return Some(SiderealRoute::SidbitEcliptic);
    }
    Some(SiderealRoute::BaseCorrection)
}

/// Echo for the SIDBIT ecliptic arm: drop the internal J2000
/// rewrite bit, keep the caller's representation bits — including an
/// explicitly requested `FLG_J2000` — and retain `FLG_NONUT` (the
/// projection is defined on a mean ecliptic). The caller echoes this
/// through the entry `_echo_request_bits` (SPEED 3/MOSEPH) and, on the
/// TT path, the ephemeris-bit convention.
pub fn sidbit_ecliptic_retflag(sub_retflag: i32, flags: i32) -> i32 {
    (sub_retflag & !FLG_J2000) | (flags & (FLG_SIDEREAL | FLG_RADIANS | FLG_XYZ | FLG_J2000))
}

type LiveAnchorRate<'a> = &'a dyn Fn(f64, bool, bool) -> Result<Option<f64>, Error>;

type EpochPairRate<'a> = &'a dyn Fn(f64, f64, f64, bool) -> Result<f64, Error>;

/// One native tropical result and its optional same-invocation word motion.
pub(crate) struct TropicalWordResult {
    pub(crate) result: (crate::positions::PositionResult, i32),
    pub(crate) motion: Option<crate::positions::derivative::OutputWordMotion>,
}
pub(crate) type TropicalWordsAt<'a> =
    dyn Fn(f64, i32, i32) -> Result<TropicalWordResult, Error> + 'a;

pub(crate) type ProjectionWordsAt<'a> = dyn Fn(
        SiderealRoute,
        f64,
        bool,
    ) -> Result<Option<crate::positions::sidereal_frame_words::WordMatrix>, Error>
    + 'a;

/// Sealed-TT sidereal serving over the native pipeline (frame projection
/// over the tropical core; scalar dispatch stays in the value layer).
///
/// Binds the tropical sealed core (`tropical_compute_tt`: TT epoch,
/// body, stripped calc flags → unshaped spherical degrees + computation
/// retflag, i.e. the `CalcEngine::compute_tt` contract) with the
/// flag-aware mean-ayanamsha provider at TT, the plain UT-layer mean
/// provider for ECL_T 0/SSY_PLANE zero points, the LEB-frame nutation
/// providers (dpsi in degrees, `_VEL_H` central-difference rate) and
/// the TT-anchored `SIDM_USER`/Valens epochs. This adapter reads its bound
/// selection, not process globals. Construction consumes no callback; later
/// evaluations can perform I/O or panic through their bound providers.
///
/// Out of scope with rationale: live star/galactic anchor geometry
/// is out of scope here — the mean provider refuses those modes without one, never
/// a silent epoch substitution); the `calc_ut` UT→TT timescale skin
/// (TT epochs in, TT epochs out).
pub struct SiderealFastServing<'a> {
    selection: SiderealSelection,
    tropical_compute_tt:
        &'a dyn Fn(f64, i32, i32) -> Result<(crate::positions::PositionResult, i32), Error>,
    mean_ayanamsha_tt: &'a dyn Fn(f64, bool, bool) -> Result<f64, Error>,
    mean_ayanamsha_ut_plain: &'a dyn Fn(f64) -> Result<f64, Error>,
    nutation_dpsi_deg: &'a dyn Fn(f64) -> Result<f64, Error>,
    nutation_rate_deg_day: &'a dyn Fn(f64) -> Result<f64, Error>,
    user_t0_tt: f64,
    valens_epoch_tt: f64,
    epoch_pair_rate: Option<EpochPairRate<'a>>,
    projection_words: Option<&'a ProjectionWordsAt<'a>>,
    tropical_words: Option<&'a TropicalWordsAt<'a>>,
    live_anchor_rate: Option<LiveAnchorRate<'a>>,
}

impl<'a> SiderealFastServing<'a> {
    /// Store the selection, defining TT epochs and borrowed providers.
    ///
    /// This passive, infallible constructor invokes no callback and performs
    /// no I/O or epoch/coverage validation. Tropical results use degrees/AU;
    /// nutation values use degrees and rates degrees per TT day. Errors and
    /// callback panics belong to the later consumed evaluation path.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        selection: SiderealSelection,
        tropical_compute_tt: &'a dyn Fn(
            f64,
            i32,
            i32,
        )
            -> Result<(crate::positions::PositionResult, i32), Error>,
        mean_ayanamsha_tt: &'a dyn Fn(f64, bool, bool) -> Result<f64, Error>,
        mean_ayanamsha_ut_plain: &'a dyn Fn(f64) -> Result<f64, Error>,
        nutation_dpsi_deg: &'a dyn Fn(f64) -> Result<f64, Error>,
        nutation_rate_deg_day: &'a dyn Fn(f64) -> Result<f64, Error>,
        user_t0_tt: f64,
        valens_epoch_tt: f64,
    ) -> Self {
        Self {
            selection,
            tropical_compute_tt,
            mean_ayanamsha_tt,
            mean_ayanamsha_ut_plain,
            nutation_dpsi_deg,
            nutation_rate_deg_day,
            user_t0_tt,
            valens_epoch_tt,
            epoch_pair_rate: None,
            projection_words: None,
            tropical_words: None,
            live_anchor_rate: None,
        }
    }

    /// Bind the complete derivative of the same native epoch-pair scalar law.
    /// Only the immutable native context installs this capability. The callback
    /// receives observation TT, defining angle in degrees, defining TT and
    /// ECL_DATE; its result is degrees per TT day. It consumes no clock provider.
    pub(crate) fn with_epoch_pair_rate(mut self, rate: EpochPairRate<'a>) -> Self {
        self.epoch_pair_rate = Some(rate);
        self
    }

    /// Bind a complete native live mean-anchor derivative, degrees per TT day.
    /// None explicitly retains an opaque provider's established sampling law.
    pub(crate) fn with_live_anchor_rate(mut self, rate: LiveAnchorRate<'a>) -> Self {
        self.live_anchor_rate = Some(rate);
        self
    }

    /// Select the already acquired defining operands; live anchors own their
    /// existing sampled law and therefore never invoke the epoch-pair callback.
    fn native_epoch_pair_rate(&self, jd_tt: f64) -> Result<Option<f64>, Error> {
        let Some(rate) = self.epoch_pair_rate else {
            return Ok(None);
        };
        let mode = self.base_mode();
        let (value, epoch) = if mode == SIDM_USER {
            (self.selection.full().2, self.user_t0_tt)
        } else if mode == SIDM_VALENS_MOON {
            (VALENS_MOON_AYAN_T0_DEG, self.valens_epoch_tt)
        } else if DYNAMIC_AYANAMSHA_MODES.contains(&mode) {
            return Ok(None);
        } else {
            // Unknown modes already emitted their diagnostic in the scalar
            // provider; its selected Fagan/Bradley pair is the derivative law.
            let (_, value, epoch) = AYANAMSHA_DEFINING
                .iter()
                .find(|(id, _, _)| *id == mode)
                .copied()
                .unwrap_or(AYANAMSHA_DEFINING[0]);
            (value, epoch)
        };
        rate(
            jd_tt,
            value,
            epoch,
            self.selection.bits() & SIDBIT_ECL_DATE != 0,
        )
        .map(Some)
    }

    /// Bind the same trusted tropical invocation with its optional word motion.
    pub(crate) fn with_tropical_words(mut self, tropical: &'a TropicalWordsAt<'a>) -> Self {
        self.tropical_words = Some(tropical);
        self
    }
    fn projection_subrequest(
        &self,
        tt: f64,
        body: i32,
        flags: i32,
    ) -> Result<TropicalWordResult, Error> {
        match self.tropical_words {
            Some(calculate) => calculate(tt, body, flags),
            None => (self.tropical_compute_tt)(tt, body, flags).map(|result| TropicalWordResult {
                result,
                motion: None,
            }),
        }
    }

    /// Native constant-frame arithmetic over the same actual scalar call.
    /// Public/custom callback constructors keep their existing tuple route.
    pub(crate) fn with_projection_words(mut self, projection: &'a ProjectionWordsAt<'a>) -> Self {
        self.projection_words = Some(projection);
        self
    }

    fn refine_projection_rates(
        &self,
        route: SiderealRoute,
        tt: f64,
        sub: [f64; 6],
        motion: Option<crate::positions::derivative::OutputWordMotion>,
        mut output: [f64; 6],
        flags: i32,
    ) -> Result<[f64; 6], Error> {
        if flags & FLG_SPEED == 0 {
            return Ok(output);
        }
        if let Some(matrix_at) = self.projection_words {
            if let Some(matrix) = matrix_at(route, tt, flags & FLG_EQUATORIAL != 0)? {
                let rate = match motion {
                    Some(motion) => {
                        let zero = crate::positions::derivative::AngularWord::from(0.0);
                        let rates = crate::positions::derivative::expanded_output_rates(
                            matrix,
                            [[zero; 3]; 3],
                            motion.position,
                            motion.velocity,
                            false,
                        )?;
                        Some([rates.0, rates.1, rates.2])
                    }
                    None => crate::positions::sidereal_frame_words::spherical_rates(sub, matrix)?,
                };
                if let Some(rate) = rate {
                    output[3..].copy_from_slice(&rate);
                }
            }
        }
        Ok(output)
    }

    /// Bound sidereal selection (mode, epoch pair, retained bits).
    pub fn selection(&self) -> SiderealSelection {
        self.selection
    }

    /// Stripped base mode driving the dispatch.
    fn base_mode(&self) -> i32 {
        strip_to_base_mode(self.selection.mode())
    }

    /// ECL_T 0 projection epoch for the bound mode (TT).
    fn projection_t0_tt(&self) -> f64 {
        ecl_t0_epoch_jd(self.base_mode(), self.user_t0_tt, self.valens_epoch_tt)
    }

    /// ECL_T 0/SSY_PLANE zero point in the projection frame: the plain mean ayanamsha at
    /// the defining epoch through the UT-layer entry — the convention with default toggles.
    fn projection_zero_point_deg(&self, route: SiderealRoute) -> Result<f64, Error> {
        match route {
            SiderealRoute::SidbitEcliptic if (self.selection.bits() & SIDBIT_ECL_T0) != 0 => {
                (self.mean_ayanamsha_ut_plain)(self.projection_t0_tt())
            }
            SiderealRoute::SidbitEcliptic => {
                (self.mean_ayanamsha_ut_plain)(J2000_TT).and_then(ssy_plane_zero_point_deg)
            }
            _ => Err(Error::ConfigurationMissing {
                missing_config: "sidereal_projection_zero_point",
                suggestion: "zero points exist only on the SIDBIT ecliptic arms",
            }),
        }
    }

    /// Sealed-TT base correction of one tropical core tuple.
    ///
    /// Equatorial output and zero-distance states pass through (the
    /// equator carries no longitude subtraction; the origin has no
    /// defined longitude). Otherwise the flag-aware mean-or-true
    /// ayanamsha is subtracted; under `FLG_SPEED` the drift is removed
    /// by the complete selected epoch-pair derivative when that native
    /// capability is bound, plus the actual nutation rate for true frames.
    /// Arbitrary callbacks keep their established general-precession or
    /// actual-provider one-second difference policy.
    /// Provider-error contract: a nutation-provider failure propagates
    /// as a typed error (`?` on `nutation_dpsi_deg`/`nutation_rate_deg_day`) —
    /// the model instead warns and continues with mean-only ayanamsha;
    /// whether a warn-and-continue numerical comparison fallback is wanted is an
    /// integration decision, not taken here.
    fn base_correct(
        &self,
        jd_tt: f64,
        mut xx: crate::positions::PositionResult,
        calc_flags: i32,
    ) -> Result<crate::positions::PositionResult, Error> {
        if (calc_flags & FLG_EQUATORIAL) != 0 || xx[2] == 0.0 {
            return Ok(xx);
        }
        let (noaberr, nogdefl, use_true) = resolve_ayanamsha_flags(calc_flags);
        let mean = (self.mean_ayanamsha_tt)(jd_tt, noaberr, nogdefl)?;
        let aya = if use_true {
            mean + (self.nutation_dpsi_deg)(jd_tt)?
        } else {
            mean
        };
        xx[0] = apply_sidereal_offset_deg(xx[0], aya)?;
        if (calc_flags & FLG_SPEED) != 0 && xx[0].is_finite() {
            let mode = self.base_mode();
            if let Some(rate) = self.native_epoch_pair_rate(jd_tt)? {
                xx[3] -= rate;
                if use_true {
                    xx[3] -= (self.nutation_rate_deg_day)(jd_tt)?;
                }
            } else if DYNAMIC_AYANAMSHA_MODES.contains(&mode) {
                if let Some(rate) = self
                    .live_anchor_rate
                    .map(|rate| rate(jd_tt, noaberr, nogdefl))
                    .transpose()?
                    .flatten()
                {
                    xx[3] -= rate;
                    if use_true {
                        xx[3] -= (self.nutation_rate_deg_day)(jd_tt)?;
                    }
                    return Ok(xx);
                }
                // Star-anchored drift: central-difference the actual
                // provider (mean or true per the active frame).
                let sample = |jd: f64| -> Result<f64, Error> {
                    let m = (self.mean_ayanamsha_tt)(jd, noaberr, nogdefl)?;
                    if use_true {
                        Ok(m + (self.nutation_dpsi_deg)(jd)?)
                    } else {
                        Ok(m)
                    }
                };
                let prev = sample(jd_tt - SIDEREAL_SPEED_HALF_STEP_DAYS)?;
                let next = sample(jd_tt + SIDEREAL_SPEED_HALF_STEP_DAYS)?;
                xx[3] -= sidereal_rate_correction_deg_day(prev, next)?;
            } else {
                xx[3] -= crate::positions::general_precession_rate_deg_day(jd_tt)?;
                if use_true {
                    xx[3] -= (self.nutation_rate_deg_day)(jd_tt)?;
                }
            }
        }
        Ok(xx)
    }

    /// Sidereal `compute_tt` over a TT epoch.
    ///
    /// `calc_flags` are the stripped computation flags (no XYZ/RADIANS
    /// shaping — the entry shapes exactly once). Returns the unshaped
    /// spherical core with its computation retflag; the entry adds the
    /// representation echo and implied bits. Each projection arm prepares
    /// its sub-request as
    /// `strip_output_flags(fixed_epoch_request_flags(calc_flags))` — a
    /// sealed-core adaptation with no counterpart in the model body
    /// (which passes `fixed_epoch_request_flags(flags)` straight into the
    /// full entry that strips/shapes itself). The TT ephemeris-bit echo
    /// duty lives with the entry-timescale skins: `tt_echo`
    /// parameter applies `_calc_tt_epheflag_echo` on the `calc` path
    /// only, while this core carries no such parameter. Tropical requests
    /// (route `None`) delegate unchanged.
    pub fn compute_tt(
        &self,
        jd_tt: f64,
        body_id: i32,
        calc_flags: i32,
    ) -> Result<(crate::positions::PositionResult, i32), Error> {
        if !jd_tt.is_finite() {
            return Err(Error::NonFiniteArgument { value: jd_tt });
        }
        let mode = self.base_mode();
        let bits = self.selection.bits();
        let Some(route) = resolve_sidereal_route(calc_flags, mode, bits) else {
            return (self.tropical_compute_tt)(jd_tt, body_id, calc_flags);
        };
        match route {
            SiderealRoute::FixedEpoch => {
                let sub_flags =
                    crate::positions::strip_output_flags(fixed_epoch_request_flags(calc_flags));
                let TropicalWordResult {
                    result: (sub_xx, sub_rf),
                    motion,
                } = self.projection_subrequest(jd_tt, body_id, sub_flags)?;
                let out = transform_fixed_epoch_result(sub_xx, calc_flags, mode)?;
                let out = self.refine_projection_rates(
                    route,
                    fixed_epoch_t0_jd(mode).unwrap(),
                    sub_xx,
                    motion,
                    out,
                    calc_flags,
                )?;
                Ok((out, fixed_epoch_retflag(sub_rf, calc_flags)))
            }
            SiderealRoute::SidbitEquatorialT0 => {
                let sub_flags =
                    crate::positions::strip_output_flags(fixed_epoch_request_flags(calc_flags));
                let TropicalWordResult {
                    result: (sub_xx, sub_rf),
                    motion,
                } = self.projection_subrequest(jd_tt, body_id, sub_flags)?;
                let out =
                    transform_equatorial_epoch_result(sub_xx, calc_flags, self.projection_t0_tt())?;
                let out = self.refine_projection_rates(
                    route,
                    self.projection_t0_tt(),
                    sub_xx,
                    motion,
                    out,
                    calc_flags,
                )?;
                let eq_rf = crate::positions::echo_request_bits(
                    fixed_epoch_retflag(sub_rf, calc_flags),
                    calc_flags,
                );
                Ok((out, eq_rf))
            }
            SiderealRoute::SidbitEcliptic => {
                let t0 = self.projection_t0_tt();
                let zero = self.projection_zero_point_deg(route)?;
                let m_ecl = sidbit_ecliptic_matrix(bits, t0, zero)?.ok_or(
                    Error::ConfigurationMissing {
                        missing_config: "sidereal_projection_matrix",
                        suggestion: "the route guarantees a projection bit; report a routing defect",
                    },
                )?;
                let sub_flags =
                    crate::positions::strip_output_flags(fixed_epoch_request_flags(calc_flags));
                let TropicalWordResult {
                    result: (sub_xx, sub_rf),
                    motion,
                } = self.projection_subrequest(jd_tt, body_id, sub_flags)?;
                let out = transform_sidbit_result(sub_xx, calc_flags, m_ecl)?;
                let out =
                    self.refine_projection_rates(route, t0, sub_xx, motion, out, calc_flags)?;
                let rf = crate::positions::echo_request_bits(
                    sidbit_ecliptic_retflag(sub_rf, calc_flags),
                    calc_flags,
                );
                Ok((out, rf))
            }
            SiderealRoute::BaseCorrection => {
                let (trop_xx, trop_rf) = (self.tropical_compute_tt)(jd_tt, body_id, calc_flags)?;
                Ok((self.base_correct(jd_tt, trop_xx, calc_flags)?, trop_rf))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{NSIDM_PREDEF, SIDBITS};

    #[test]
    fn defining_table_shape_and_locks() {
        // 32 static epoch pairs; every entry is a predefined non-dynamic
        // base mode, keys unique and ascending (transcription lock).
        assert_eq!(AYANAMSHA_DEFINING.len(), 32);
        let mut prev = -1;
        for (mode, value, epoch) in AYANAMSHA_DEFINING {
            assert!(mode > prev, "defining table out of order at {mode}");
            prev = mode;
            assert!((0..=MAX_PREDEFINED_MODE).contains(&mode));
            assert!(!DYNAMIC_AYANAMSHA_MODES.contains(&mode));
            assert!(value.is_finite() && epoch.is_finite());
        }
        // Spot pins against the model literals.
        assert_eq!(
            AYANAMSHA_DEFINING[0],
            (0, 24.042044444444457, 2433282.42345905)
        );
        assert_eq!(AYANAMSHA_DEFINING[1], (1, 23.85709166666667, 2451545.0));
        assert_eq!(AYANAMSHA_DEFINING[12], (12, -4.466666666666667, 1684613.0));
    }

    #[test]
    fn dynamic_set_shape_and_partition() {
        // 15 dynamic modes; together with the 32 defining pairs they cover
        // every predefined base mode exactly once.
        assert_eq!(DYNAMIC_AYANAMSHA_MODES.len(), 15);
        let mut covered = [false; 47];
        for (mode, _, _) in AYANAMSHA_DEFINING {
            assert!(!covered[mode as usize]);
            covered[mode as usize] = true;
        }
        for mode in DYNAMIC_AYANAMSHA_MODES {
            assert!(!covered[mode as usize], "mode {mode} in both tables");
            covered[mode as usize] = true;
        }
        assert!(covered.iter().all(|c| *c), "predefined mode uncovered");
        assert_eq!(NSIDM_PREDEF, 47);
        assert_eq!(DYNAMIC_MODE_ANCHORS.len(), 15);
    }

    #[test]
    fn anchor_catalog_shape_and_lookup() {
        // Seven baked anchor entries in catalog row order; keys unique,
        // every entry finite, unknown keys (incl. the two formula
        // markers) refuse with None.
        assert_eq!(ANCHOR_STARS.len(), 7);
        let mut seen = std::collections::HashSet::new();
        for (key, star) in ANCHOR_STARS {
            assert!(seen.insert(key), "duplicate anchor key {key}");
            assert!(star.ra_j2000_deg().is_finite());
            assert!(star.dec_j2000_deg().is_finite());
            assert!(star.pm_ra_arcsec_yr().is_finite());
            assert!(star.pm_dec_arcsec_yr().is_finite());
            assert!(star.parallax_mas().is_finite());
            assert!(star.radial_velocity_km_s().is_finite());
            assert_eq!(anchor_star(key), Some(star));
        }
        assert_eq!(anchor_star("SPICA").unwrap().ra_j2000_deg(), 201.2982475);
        assert_eq!(anchor_star("SIRIUS"), None);
        assert_eq!(anchor_star("IAU2006_PRECESSION_FORMULA"), None);
        assert_eq!(anchor_star("EPOCH_PAIR"), None);
        assert_eq!(anchor_star("spica"), None);
    }

    #[test]
    fn anchor_catalog_matches_release_words() {
        // exactly transcription lock against `STARS` values
        // (parallax in milliarcseconds, checked identical to the model
        // `arcsec * 1000.0` product; GAL_CENTER pm_ra is the model
        // cos(dec) expression evaluated once).
        let words = |name: &str| {
            let star = anchor_star(name).unwrap();
            [
                star.ra_j2000_deg().to_bits(),
                star.dec_j2000_deg().to_bits(),
                star.pm_ra_arcsec_yr().to_bits(),
                star.pm_dec_arcsec_yr().to_bits(),
                star.parallax_mas().to_bits(),
                star.radial_velocity_km_s().to_bits(),
            ]
        };
        assert_eq!(
            words("SPICA"),
            [
                0x4069298b3e5753a4,
                0xc0265298c8b941d6,
                0xbfa5c28f5c28f5c3,
                0xbfa03eea209aaa3b,
                0x4028e147ae147ae1,
                0x3ff0000000000000,
            ]
        );
        assert_eq!(
            words("REVATI"),
            [
                0x40326ecfcdc8324f,
                0x401e4d2b36b554d3,
                0x3fc243c3a42f1ed1,
                0xbfab29802c0a4a06,
                0x403875a1cac08312,
                0x402e000000000000,
            ]
        );
        assert_eq!(
            words("ALDEBARAN"),
            [
                0x40513ebafcb5f397,
                0x40308261a38d5227,
                0x3fb03e425aee6320,
                0xbfc82f2f9873ffac,
                0x40487851eb851eb8,
                0x404b2147ae147ae1,
            ]
        );
        assert_eq!(
            words("PUSHYA"),
            [
                0x4060657ad918a60d,
                0x40322780bbce4b1d,
                0xbf92e09fe86833c6,
                0xbfcd28f9f44d4456,
                0x4037d3bcd35a8588,
                0x403123d70a3d70a4,
            ]
        );
        assert_eq!(
            words("MULA"),
            [
                0x4070766f46d7aacc,
                0xc0428d4a1725608d,
                0xbf81782d38476f2a,
                0xbf9f8a0902de00d2,
                0x4016d70a3d70a3d7,
                0xc008000000000000,
            ]
        );
        assert_eq!(
            words("GAL_CENTER"),
            [
                0x4070a6ab5d5ba5c3,
                0xc03d01ffdf73c371,
                0xbf669329104fc0d0,
                0xbf76b87379746024,
                0x0,
                0x0,
            ]
        );
        assert_eq!(
            words("GAL_NORTH_POLE"),
            [0x40681b80dc33721d, 0x403b20d4fdf3b646, 0x0, 0x0, 0x0, 0x0,]
        );
    }

    #[test]
    fn anchor_catalog_covers_dynamic_star_modes() {
        // Every catalog-backed anchor named by DYNAMIC_MODE_ANCHORS
        // resolves in the baked table; the Mardyks/Valens formula
        // markers name no entry by design.
        for (_, anchor) in DYNAMIC_MODE_ANCHORS {
            if anchor == "IAU2006_PRECESSION_FORMULA" || anchor == "EPOCH_PAIR" {
                assert_eq!(anchor_star(anchor), None);
            } else {
                assert!(
                    anchor_star(anchor).is_some(),
                    "dynamic anchor {anchor} missing from catalog"
                );
            }
        }
    }

    #[test]
    fn sidbit_sets_match_release_bits() {
        assert_eq!(SIDBITS, 256);
        assert_eq!(APPLIED_SIDBITS, 256 + 512 + 1024 + 2048);
        assert_eq!(ACCEPTED_SIDBITS, 4096 + 8192);
        assert_eq!(IMPLEMENTED_SIDBITS, APPLIED_SIDBITS | ACCEPTED_SIDBITS);
        assert_eq!(IMPLEMENTED_SIDBITS, 256 + 512 + 1024 + 2048 + 4096 + 8192);
    }

    #[test]
    fn names_cover_all_predefined_modes() {
        assert_eq!(AYANAMSHA_NAMES.len(), 47);
        for mode in 0..=MAX_PREDEFINED_MODE {
            assert!(!ayanamsha_name(mode).is_empty(), "mode {mode} unnamed");
        }
        // the model string pins (incl. quoting and abbreviations).
        assert_eq!(ayanamsha_name(0), "Fagan/Bradley");
        assert_eq!(ayanamsha_name(39), "\"Vedic\"/Sheoran");
        assert_eq!(ayanamsha_name(34), "Skydram (Mardyks)");
        // Unassigned block, SIDM_USER and negatives carry no name.
        assert_eq!(ayanamsha_name(47), "");
        assert_eq!(ayanamsha_name(254), "");
        assert_eq!(ayanamsha_name(SIDM_USER), "");
        assert_eq!(ayanamsha_name(-1), "");
        // Low-byte masking for SIDBIT composites.
        assert_eq!(ayanamsha_name(256 | 1), "Lahiri");
        assert_eq!(
            ayanamsha_name(256 + 512 + 1024 + 2048 + 4096 + 8192 + 1),
            "Lahiri"
        );
    }

    #[test]
    fn names_verbatim_all_47() {
        // the full `get_ayanamsa_name` table, transcribed
        // value-for-value (ephemeral probe 2026-09-23 against peeled
        // `37ffb05...`; scratch only, nothing stored). Guards transcription
        // slips the non-empty check in `names_cover_all_predefined_modes`
        // cannot see.
        const EXPECTED: [(i32, &str); 47] = [
            (0, "Fagan/Bradley"),
            (1, "Lahiri"),
            (2, "De Luce"),
            (3, "Raman"),
            (4, "Usha/Shashi"),
            (5, "Krishnamurti"),
            (6, "Djwhal Khul"),
            (7, "Yukteshwar"),
            (8, "J.N. Bhasin"),
            (9, "Babylonian/Kugler 1"),
            (10, "Babylonian/Kugler 2"),
            (11, "Babylonian/Kugler 3"),
            (12, "Babylonian/Huber"),
            (13, "Babylonian/Eta Piscium"),
            (14, "Babylonian/Aldebaran = 15 Tau"),
            (15, "Hipparchos"),
            (16, "Sassanian"),
            (17, "Galact. Center = 0 Sag"),
            (18, "J2000"),
            (19, "J1900"),
            (20, "B1950"),
            (21, "Suryasiddhanta"),
            (22, "Suryasiddhanta, mean Sun"),
            (23, "Aryabhata"),
            (24, "Aryabhata, mean Sun"),
            (25, "SS Revati"),
            (26, "SS Citra"),
            (27, "True Citra"),
            (28, "True Revati"),
            (29, "True Pushya (PVRN Rao)"),
            (30, "Galactic Center (Gil Brand)"),
            (31, "Galactic Equator (IAU1958)"),
            (32, "Galactic Equator"),
            (33, "Galactic Equator mid-Mula"),
            (34, "Skydram (Mardyks)"),
            (35, "True Mula (Chandra Hari)"),
            (36, "Dhruva/Gal.Center/Mula (Wilhelm)"),
            (37, "Aryabhata 522"),
            (38, "Babylonian/Britton"),
            (39, "\"Vedic\"/Sheoran"),
            (40, "Cochrane (Gal.Center = 0 Cap)"),
            (41, "Galactic Equator (Fiorenza)"),
            (42, "Vettius Valens"),
            (43, "Lahiri 1940"),
            (44, "Lahiri VP285"),
            (45, "Krishnamurti-Senthilathiban"),
            (46, "Lahiri ICRC"),
        ];
        assert_eq!(AYANAMSHA_NAMES.len(), EXPECTED.len());
        for (id, want) in EXPECTED {
            assert_eq!(ayanamsha_name(id), want, "mode {id}");
        }
    }

    #[test]
    fn base_mode_and_registration_rules() {
        assert_eq!(strip_to_base_mode(256 | 5), 5);
        assert_eq!(strip_to_base_mode(3), 3);
        assert_eq!(strip_to_base_mode(-1), -1);
        assert!(is_registered_mode(0) && is_registered_mode(46));
        assert!(is_registered_mode(SIDM_USER));
        assert!(!is_registered_mode(47) && !is_registered_mode(254));
        assert!(!is_registered_mode(-1));
    }

    #[test]
    fn set_get_roundtrip_plain_modes() {
        let mut sel = SiderealSelection::unset();
        assert_eq!(sel.full(), (0, 0.0, 0.0));
        assert_eq!(sel.bits(), 0);
        // Plain base mode: stored, no diagnostics.
        assert!(sel.set_sidereal_mode(1, 0.0, 0.0).is_empty());
        assert_eq!(sel.full(), (1, 0.0, 0.0));
        assert_eq!(sel.bits(), 0);
        // SIDM_USER keeps its epoch pair.
        assert!(sel.set_sidereal_mode(SIDM_USER, 2451545.0, 24.0).is_empty());
        assert_eq!(sel.full(), (SIDM_USER, 2451545.0, 24.0));
    }

    #[test]
    fn implemented_bits_retained_silently() {
        let mut sel = SiderealSelection::unset();
        let warnings = sel.set_sidereal_mode(256 | 512 | 1024 | 1, 0.0, 0.0);
        assert!(warnings.is_empty());
        assert_eq!(sel.mode(), 1);
        assert_eq!(sel.bits(), 256 | 512 | 1024);
        // Accepted-but-unrealized bits are likewise silent.
        let warnings = sel.set_sidereal_mode(4096 | 8192 | 2, 0.0, 0.0);
        assert!(warnings.is_empty());
        assert_eq!(sel.mode(), 2);
        assert_eq!(sel.bits(), 4096 | 8192);
    }

    #[test]
    fn unknown_high_bits_stripped_with_diagnostic() {
        let mut sel = SiderealSelection::unset();
        let warnings = sel.set_sidereal_mode(0x4000 | 1, 0.0, 0.0);
        assert_eq!(
            warnings,
            vec![Warning::UnsupportedSiderealProjection { bits: 0x4000 }]
        );
        assert_eq!(sel.mode(), 1);
        assert_eq!(sel.bits(), 0);
        // Implemented bits survive while unknown ones are reported.
        let warnings = sel.set_sidereal_mode(256 | 0x8000 | 5, 1.5, 2.5);
        assert_eq!(
            warnings,
            vec![Warning::UnsupportedSiderealProjection { bits: 0x8000 }]
        );
        assert_eq!(sel.full(), (5, 1.5, 2.5));
        assert_eq!(sel.bits(), 256);
    }

    #[test]
    fn negative_modes_pass_through_verbatim() {
        let mut sel = SiderealSelection::unset();
        // Scalar ayanamsha evaluation.
        //
        // Epoch-defined modes use their defining TT epochs and angles. Mardyks uses
        // the IAU 2006 general-precession angle; true ayanamsha adds IAU 2006/2000 A
        // nutation in longitude without renormalizing the final sum. Stellar and
        // galactic definitions require an explicit DynamicAnchorProvider; a missing
        // anchor returns AyanamshaAnchorUnserved rather than an epoch substitution.
        // TT/UT1 conversion uses the supplied Delta-T model.
        let warnings = sel.set_sidereal_mode(-1, 0.0, 0.0);
        assert!(warnings.is_empty());
        assert_eq!(sel.mode(), -1);
        assert_eq!(sel.bits(), 0);
    }

    #[test]
    fn set_sidereal_mode_i32_boundaries() {
        // the two-stage rule at the `i32` edges. `i32::MAX` is
        // non-negative, so every implemented SIDBIT is retained silently,
        // the remaining high bits are reported once with their exact value,
        // and the base reduces to the low byte (`0xFF` = `SIDM_USER`).
        let mut sel = SiderealSelection::unset();
        let warnings = sel.set_sidereal_mode(i32::MAX, 1.5, 2.5);
        assert_eq!(
            warnings,
            vec![Warning::UnsupportedSiderealProjection {
                bits: i32::MAX & !IMPLEMENTED_SIDBITS & !BASE_MODE_MASK
            }]
        );
        assert_eq!(sel.mode(), SIDM_USER);
        assert_eq!(sel.full(), (SIDM_USER, 1.5, 2.5));
        assert_eq!(sel.bits(), IMPLEMENTED_SIDBITS);
        // `i32::MIN` is negative: an invalid ID, not a composite — stored
        // with no retained flags and no diagnostics.
        let mut sel = SiderealSelection::unset();
        let warnings = sel.set_sidereal_mode(i32::MIN, -0.0, f64::NAN);
        assert!(warnings.is_empty());
        assert_eq!(sel.mode(), i32::MIN);
        assert_eq!(sel.bits(), 0);
        assert!(sel.full().1.is_sign_negative());
        assert!(sel.full().2.is_nan());
        // Registration at the edges: `i32::MAX` masks to `SIDM_USER`
        // (registered), `i32::MIN` is negative (unregistered).
        assert!(is_registered_mode(i32::MAX));
        assert!(!is_registered_mode(i32::MIN));
        assert_eq!(strip_to_base_mode(i32::MIN), i32::MIN);
    }

    #[test]
    fn fixed_epoch_inventory() {
        assert_eq!(FIXED_EPOCH_MODES, [18, 19, 20, 34]);
        for mode in FIXED_EPOCH_MODES {
            assert!(is_fixed_epoch_mode(mode));
            // Mode 34 is deliberately both: the dynamic branch evaluates
            // its scalar ayanamsha from the IAU 2006 formula while the
            // fixed-epoch path projects positions onto its t 0 frame.
            if mode == SIDM_GALALIGN_MARDYKS {
                assert!(DYNAMIC_AYANAMSHA_MODES.contains(&mode));
            } else {
                assert!(!DYNAMIC_AYANAMSHA_MODES.contains(&mode));
            }
        }
        assert!(!is_fixed_epoch_mode(1));
        assert_eq!(fixed_epoch_lon_offset_deg(34), 30.0);
        assert_eq!(fixed_epoch_lon_offset_deg(18), 0.0);
    }

    // Method-B kernel: analytic properties (no the model values).
    #[test]
    fn method_b_zero_at_own_epoch() {
        // The equinox of t 0 projected into its own ecliptic frame has
        // longitude ~0, so the accumulation vanishes (rounding-level only).
        for t0 in [J2000_TT, 2415020.0, 1903395.1338202194, 1684613.0] {
            let acc = method_b_accumulated_precession(t0, t0).unwrap();
            assert!(acc.abs() < 1e-9, "t0={t0}: acc={acc}");
        }
    }

    #[test]
    fn method_b_sign_and_century_scale() {
        // The equinox regresses: forward accumulation is positive and close
        // to the general precession (~5029"/cy ~ 1.397 deg/cy); backward is
        // the near-mirror image. Bounds are a-priori (paper scale), not
        // fitted: the exact values come from the model sweep.
        let fwd = method_b_accumulated_precession(J2000_TT + 36525.0, J2000_TT).unwrap();
        assert!(fwd > 1.39 && fwd < 1.41, "fwd={fwd}");
        let bwd = method_b_accumulated_precession(J2000_TT - 36525.0, J2000_TT).unwrap();
        assert!(bwd < -1.39 && bwd > -1.41, "bwd={bwd}");
        // Deliberately unreduced: a 400-year accumulation exceeds 360 only
        // far outside this window; here it must equal the raw sum scale.
        let long = method_b_accumulated_precession(J2000_TT + 365250.0, J2000_TT).unwrap();
        assert!(long > 13.9 && long < 14.1, "long={long}");
    }

    #[test]
    fn method_b_non_additive_in_t0() {
        // The documented non-additivity: splitting across an intermediate
        // epoch changes the result because each leg refers to a different
        // fixed ecliptic. The residual is second-order in the precession
        // angles — small, but far above rounding over century legs.
        let t0 = 2415020.0; // 1900
        let t1 = 2451545.0; // 2000
        let t2 = 2488070.0; // 2100
        let direct = method_b_accumulated_precession(t2, t0).unwrap();
        let split = method_b_accumulated_precession(t2, t1).unwrap()
            + method_b_accumulated_precession(t1, t0).unwrap();
        let residual = (direct - split).abs();
        assert!(
            residual > 1e-12,
            "expected non-additive residual, got {residual}"
        );
        assert!(residual < 1e-3, "residual implausibly large: {residual}");
    }

    #[test]
    fn ecl_date_delta_zero_at_own_epoch() {
        // Referring the zero point to its own epoch's ecliptic is a no-op.
        for (value, epoch) in [
            (23.85709166666667, 2451545.0),
            (0.0, 1903395.1338202194),
            (-4.466666666666667, 1684613.0),
        ] {
            let d = ecl_date_ayanamsha_delta_deg(value, epoch, epoch).unwrap();
            assert!(d.abs() < 1e-9, "value={value}: delta={d}");
            // And the epoch-pair scalar reproduces the defining pair there,
            // up to the 0/360 wrap: value 0.0 with a -4 e-16 accumulation
            // reduces to exactly 360.0 on both sides (an ulp of 360 is
            // 5.7 e-14, so 360 - 4 e-16 rounds to 360.0), hence the
            // circular distance below rather than a linear one.
            let circ = |a: f64| (a - value + 540.0).rem_euclid(360.0) - 180.0;
            let a = epoch_pair_ayanamsha_deg(epoch, value, epoch, false).unwrap();
            assert!(circ(a).abs() < 1e-9, "a={a}");
            let ae = epoch_pair_ayanamsha_deg(epoch, value, epoch, true).unwrap();
            assert!(circ(ae).abs() < 1e-9, "ae={ae}");
        }
    }

    #[test]
    fn epoch_pair_j2000_lahiri_pin() {
        // At its defining epoch the Lahiri scalar is the table value.
        let a = epoch_pair_ayanamsha_deg(2451545.0, 23.85709166666667, 2451545.0, false).unwrap();
        assert!((a - 23.85709166666667).abs() < 1e-9, "a={a}");
    }

    #[test]
    fn method_b_hostile_inputs_refused() {
        let bad = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
        for v in bad {
            assert!(method_b_accumulated_precession(v, J2000_TT).is_err());
            assert!(method_b_accumulated_precession(J2000_TT, v).is_err());
            assert!(ecl_date_ayanamsha_delta_deg(v, J2000_TT, J2000_TT).is_err());
            assert!(ecl_date_ayanamsha_delta_deg(0.0, v, J2000_TT).is_err());
            assert!(ecl_date_ayanamsha_delta_deg(0.0, J2000_TT, v).is_err());
            assert!(epoch_pair_ayanamsha_deg(v, 0.0, J2000_TT, false).is_err());
            assert!(epoch_pair_ayanamsha_deg(J2000_TT, v, J2000_TT, false).is_err());
            assert!(epoch_pair_ayanamsha_deg(J2000_TT, 0.0, v, true).is_err());
        }
    }

    // value layer: analytic properties (no the model values).

    /// Fixed stub anchor provider: returns scripted longitudes and
    /// records the received toggles so the `eff_` gating is observable.
    struct StubAnchors {
        lon: f64,
        seen: std::cell::RefCell<Vec<(DynamicAnchor, bool, bool)>>,
    }

    impl StubAnchors {
        fn new(lon: f64) -> Self {
            Self {
                lon,
                seen: std::cell::RefCell::new(Vec::new()),
            }
        }
    }

    impl DynamicAnchorProvider for StubAnchors {
        fn anchor_longitude_deg(
            &self,
            anchor: DynamicAnchor,
            jd_tt: f64,
            noaberr: bool,
            nogdefl: bool,
        ) -> Result<f64, Error> {
            assert!(jd_tt.is_finite());
            self.seen.borrow_mut().push((anchor, noaberr, nogdefl));
            Ok(self.lon)
        }
    }

    #[test]
    fn mardyks_defining_pair_and_scale() {
        // Exactly 30° at the September-equinox-1998 defining epoch, and
        // forward accumulation at the general-precession rate (~1.4
        // deg/cy — a-priori paper scale, not fitted).
        let (_, epoch) = MARDYKS_DEFINING;
        let a0 = mardyks_ayanamsha_deg(epoch).unwrap();
        assert!((a0 - 30.0).abs() < 1e-9, "a0={a0}");
        let a1 = mardyks_ayanamsha_deg(epoch + 36525.0).unwrap();
        let d = (a1 - a0).rem_euclid(360.0);
        assert!(d > 1.39 && d < 1.41, "d={d}");
        assert!(mardyks_ayanamsha_deg(f64::NAN).is_err());
    }

    #[test]
    fn true_value_adds_signed_nutation_without_renormalizing() {
        // At J2000 the IAU 2006 dpsi is negative (~-0.0039°), so the
        // true value of a zero mean reads slightly negative — the
        // the model no-renormalization contract, not a wrap to 359.99.
        let dpsi = nutation_in_longitude_deg(J2000_TT).unwrap();
        assert!(dpsi < -0.003 && dpsi > -0.005, "dpsi={dpsi}");
        let t = true_ayanamsha_deg(0.0, J2000_TT).unwrap();
        assert!((t - dpsi).abs() < 1e-15, "t={t}");
        assert!(t < 0.0, "expected no renormalization, got {t}");
        // Composition identity: true − mean is exactly the increment.
        let m = 23.85709166666667;
        let t2 = true_ayanamsha_deg(m, 2460000.5).unwrap();
        let d2 = nutation_in_longitude_deg(2460000.5).unwrap();
        assert!((t2 - m - d2).abs() < 1e-15);
        assert!(true_ayanamsha_deg(f64::INFINITY, J2000_TT).is_err());
        assert!(true_ayanamsha_deg(0.0, f64::NAN).is_err());
    }

    #[test]
    fn ex_retflag_single_winner_and_nonut_drop() {
        use crate::constants::{FLG_JPLEPH, FLG_MOSEPH, FLG_SWIEPH};
        // Ephemeris-bit priority with the SWIEPH default; every other
        // input bit is dropped.
        assert_eq!(ayanamsha_ex_retflag(0, 1), FLG_SWIEPH);
        assert_eq!(ayanamsha_ex_retflag(FLG_JPLEPH, 1), FLG_JPLEPH);
        assert_eq!(ayanamsha_ex_retflag(FLG_MOSEPH, 1), FLG_MOSEPH);
        assert_eq!(ayanamsha_ex_retflag(256, 1), FLG_SWIEPH);
        assert_eq!(ayanamsha_ex_retflag(FLG_SWIEPH | FLG_MOSEPH, 1), FLG_SWIEPH);
        assert_eq!(ayanamsha_ex_retflag(FLG_JPLEPH | FLG_SWIEPH, 1), FLG_JPLEPH);
        // FLG_NONUT echoes for ordinary modes (0 → 2, NONUT → 66,
        // JPLEPH|NONUT → 65 in the model spellings).
        assert_eq!(ayanamsha_ex_retflag(FLG_NONUT, 1), FLG_SWIEPH | FLG_NONUT);
        assert_eq!(
            ayanamsha_ex_retflag(FLG_JPLEPH | FLG_NONUT, 1),
            FLG_JPLEPH | FLG_NONUT
        );
        // ... but drops for the 12 star/galactic calculated modes.
        assert_eq!(AYANAMSA_EX_NONUT_DROP_MODES.len(), 12);
        for mode in AYANAMSA_EX_NONUT_DROP_MODES {
            assert_eq!(ayanamsha_ex_retflag(FLG_NONUT, mode), FLG_SWIEPH);
            assert_eq!(
                ayanamsha_ex_retflag(FLG_JPLEPH | FLG_NONUT, mode),
                FLG_JPLEPH
            );
        }
    }

    #[test]
    fn flag_resolution_mean_vs_true_and_center_toggles() {
        use crate::constants::{
            FLG_BARYCTR, FLG_HELCTR, FLG_J2000, FLG_NOABERR, FLG_NOGDEFL, FLG_NONUT, FLG_TRUEPOS,
        };
        // Default: true value, no anchor toggles.
        assert_eq!(resolve_ayanamsha_flags(0), (false, false, true));
        assert_eq!(resolve_ayanamsha_flags(FLG_NONUT), (false, false, false));
        assert_eq!(resolve_ayanamsha_flags(FLG_J2000), (false, false, false));
        assert_eq!(
            resolve_ayanamsha_flags(FLG_NONUT | FLG_J2000),
            (false, false, false)
        );
        // Individual toggles; TRUEPOS sets both.
        assert_eq!(resolve_ayanamsha_flags(FLG_NOABERR), (true, false, true));
        assert_eq!(resolve_ayanamsha_flags(FLG_NOGDEFL), (false, true, true));
        assert_eq!(resolve_ayanamsha_flags(FLG_TRUEPOS), (true, true, true));
        // Helio/barycentric centers set both toggles (no Earth to be
        // aberrated by), keeping the true/mean selection orthogonal.
        assert_eq!(resolve_ayanamsha_flags(FLG_HELCTR), (true, true, true));
        assert_eq!(
            resolve_ayanamsha_flags(FLG_BARYCTR | FLG_NONUT),
            (true, true, false)
        );
    }

    #[test]
    fn dispatch_epoch_modes_reproduce_defining_pairs() {
        // Sidereal defining epochs, reference angles and frame projections.
        //
        // A defining epoch and ayanamsha identify the zero point. The Vondrák pole
        // geometry carries that definition to the requested TT epoch. Projection
        // flags choose the ecliptic of the defining epoch or the invariable plane;
        // the actual frame matrices are shared with the position pipeline.
        for (mode, value, epoch) in AYANAMSHA_DEFINING {
            // Mardyks (34) is dynamic-formula, not a defining pair.
            if mode == SIDM_GALALIGN_MARDYKS {
                continue;
            }
            let (a, w) =
                calc_mean_ayanamsha_deg(epoch, mode, (0.0, 0.0), 0.0, false, false, false, None)
                    .unwrap();
            assert!(w.is_none(), "mode={mode}");
            let circ = (a - value + 540.0).rem_euclid(360.0) - 180.0;
            assert!(circ.abs() < 1e-9, "mode={mode}: a={a}");
            // ECL_DATE is a no-op at the defining epoch itself.
            let (ae, _) =
                calc_mean_ayanamsha_deg(epoch, mode, (0.0, 0.0), 0.0, true, false, false, None)
                    .unwrap();
            let circe = (ae - value + 540.0).rem_euclid(360.0) - 180.0;
            assert!(circe.abs() < 1e-9, "mode={mode}: ae={ae}");
        }
    }

    #[test]
    fn consumed_unknown_mode_observer_precedes_fallback_and_is_fallible() {
        use std::cell::RefCell;
        let seen = RefCell::new(Vec::new());
        let observe = |mode| {
            seen.borrow_mut().push(mode);
            Ok(())
        };
        let delta = |_| Ok(0.0);
        for mode in [SIDM_FAGAN_BRADLEY, 99, -1] {
            let mut selection = SiderealSelection::unset();
            selection.set_sidereal_mode(mode, 0., 0.);
            let aya = SealedAyanamshaServing::new(selection, &delta, None)
                .with_unknown_mode_observer(Some(&observe));
            let (value, warning) = aya.mean_tt_with(J2000_TT, false, false).unwrap();
            let ordinary = SealedAyanamshaServing::new(selection, &delta, None)
                .mean_tt_with(J2000_TT, false, false)
                .unwrap();
            assert_eq!(value.to_bits(), ordinary.0.to_bits());
            assert_eq!(warning, ordinary.1);
        }
        assert_eq!(*seen.borrow(), [99, -1]);
        let mut selection = SiderealSelection::unset();
        selection.set_sidereal_mode(99, 0., 0.);
        let refuse = |_| Err(Error::UnknownBody { body_id: 983 });
        let aya = SealedAyanamshaServing::new(selection, &delta, None)
            .with_unknown_mode_observer(Some(&refuse));
        assert!(matches!(
            aya.ex_ut(J2000_TT, 0),
            Err(Error::UnknownBody { body_id: 983 })
        ));
        assert!(matches!(
            aya.mean_tt_with(1.0e308, false, false),
            Err(Error::UnknownBody { body_id: 983 })
        ));
        assert!(matches!(
            aya.mean_tt_with(f64::NAN, false, false),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn dispatch_unknown_mode_warns_and_falls_back_to_fagan() {
        // Unassigned positives and negatives warn and evaluate the
        // Fagan/Bradley pair.
        for bad in [47, 100, 254, -1, -7] {
            let (a, w) =
                calc_mean_ayanamsha_deg(J2000_TT, bad, (0.0, 0.0), 0.0, false, false, false, None)
                    .unwrap();
            assert_eq!(w, Some(Warning::UnknownSiderealMode { mode: bad }));
            let (f, wf) =
                calc_mean_ayanamsha_deg(J2000_TT, 0, (0.0, 0.0), 0.0, false, false, false, None)
                    .unwrap();
            assert!(wf.is_none());
            assert!((a - f).abs() < 1e-15, "bad={bad}: a={a} f={f}");
        }
        // SIDM_USER never warns, even with a wild stored pair.
        let (_, w) = calc_mean_ayanamsha_deg(
            J2000_TT,
            SIDM_USER,
            (0.0, 0.0),
            0.0,
            false,
            false,
            false,
            None,
        )
        .unwrap();
        assert!(w.is_none());
    }

    #[test]
    fn dispatch_user_and_valens_epoch_pairs() {
        // SIDM_USER is the stored pair plus Method-B accumulation.
        let (a, _) = calc_mean_ayanamsha_deg(
            2451545.0,
            SIDM_USER,
            (2451545.0, 24.0),
            0.0,
            false,
            false,
            false,
            None,
        )
        .unwrap();
        let circ = (a - 24.0 + 540.0).rem_euclid(360.0) - 180.0;
        assert!(circ.abs() < 1e-9, "a={a}");
        // Valens reads its defining offset at its defining TT epoch.
        let vt = VALENS_MOON_T0_UT + 69.0 / 86400.0; // ~69 s ΔT stand-in
        let (v, _) = calc_mean_ayanamsha_deg(
            vt,
            SIDM_VALENS_MOON,
            (0.0, 0.0),
            vt,
            false,
            false,
            false,
            None,
        )
        .unwrap();
        let circv = (v - VALENS_MOON_AYAN_T0_DEG + 540.0).rem_euclid(360.0) - 180.0;
        assert!(circv.abs() < 1e-9, "v={v}");
    }

    #[test]
    fn dispatch_dynamic_arms_use_verbatim_targets() {
        // Stub anchor at 200°: each arm subtracts its target.
        let stub = StubAnchors::new(200.0);
        let some = Some(&stub as &dyn DynamicAnchorProvider);
        let cases = [
            (SIDM_TRUE_CITRA, 20.0, DynamicAnchor::Spica),
            (SIDM_TRUE_REVATI, 200.0 + 10.0 / 60.0, DynamicAnchor::Revati),
            (SIDM_TRUE_PUSHYA, 94.0, DynamicAnchor::Pushya),
            (SIDM_TRUE_MULA, 320.0, DynamicAnchor::Mula),
            (
                SIDM_TRUE_SHEORAN,
                200.0 - SHEORAN_TARGET_LON,
                DynamicAnchor::Pushya,
            ),
            (SIDM_ALDEBARAN_15TAU, 155.0, DynamicAnchor::Aldebaran),
            (SIDM_GALCENT_0SAG, 320.0, DynamicAnchor::GalCenter),
            (
                SIDM_GALCENT_RGILBRAND,
                200.0 - 244.37694101250946,
                DynamicAnchor::GalCenter,
            ),
            (SIDM_GALCENT_COCHRANE, 290.0, DynamicAnchor::GalCenter),
            (
                SIDM_GALCENT_MULA_WILHELM,
                200.0 - MULA_WILHELM_TARGET_LON,
                DynamicAnchor::MulaWilhelm,
            ),
        ];
        for (mode, expected, anchor) in cases {
            stub.seen.borrow_mut().clear();
            let (a, w) =
                calc_mean_ayanamsha_deg(J2000_TT, mode, (0.0, 0.0), 0.0, false, true, true, some)
                    .unwrap();
            assert!(w.is_none(), "mode={mode}");
            let circ = (a - expected.rem_euclid(360.0) + 540.0).rem_euclid(360.0) - 180.0;
            assert!(circ.abs() < 1e-9, "mode={mode}: a={a}");
            // The requested anchor was consumed exactly once.
            let seen = stub.seen.borrow();
            assert_eq!(seen.len(), 1, "mode={mode}");
            assert_eq!(seen[0].0, anchor, "mode={mode}");
        }
        // Galactic-equator arms: ascending node = pole + 90, minus target.
        for (mode, target) in GALEQU_TARGET_LON {
            stub.seen.borrow_mut().clear();
            let (a, _) =
                calc_mean_ayanamsha_deg(J2000_TT, mode, (0.0, 0.0), 0.0, false, false, false, some)
                    .unwrap();
            let expected = (200.0f64 + 90.0).rem_euclid(360.0) - target;
            let circ = (a - expected.rem_euclid(360.0) + 540.0).rem_euclid(360.0) - 180.0;
            assert!(circ.abs() < 1e-9, "mode={mode}: a={a}");
            // The geometric pole ignores both toggles even when set.
            let seen = stub.seen.borrow();
            assert_eq!(
                seen[0],
                (DynamicAnchor::GalNorthPole, false, false),
                "mode={mode}"
            );
        }
    }

    #[test]
    fn dispatch_dynamic_gating_and_refusals() {
        let stub = StubAnchors::new(200.0);
        let some = Some(&stub as &dyn DynamicAnchorProvider);
        // Without a provider every live arm refuses (never a silent
        // epoch substitution); formula arms stay servable.
        for mode in DYNAMIC_AYANAMSHA_MODES {
            let out = calc_mean_ayanamsha_deg(
                J2000_TT,
                mode,
                (0.0, 0.0),
                J2000_TT,
                false,
                false,
                false,
                None,
            );
            if mode == SIDM_GALALIGN_MARDYKS || mode == SIDM_VALENS_MOON {
                assert!(out.is_ok(), "mode={mode}");
            } else {
                assert_eq!(out, Err(Error::AyanamshaAnchorUnserved { mode }));
            }
        }
        // Aberrant Citra forwards the toggles; non-aberrant Aldebaran
        // receives them cleared (its anchor is catalog-fixed).
        stub.seen.borrow_mut().clear();
        calc_mean_ayanamsha_deg(
            J2000_TT,
            SIDM_TRUE_CITRA,
            (0.0, 0.0),
            0.0,
            false,
            true,
            true,
            some,
        )
        .unwrap();
        assert_eq!(stub.seen.borrow()[0], (DynamicAnchor::Spica, true, true));
        stub.seen.borrow_mut().clear();
        calc_mean_ayanamsha_deg(
            J2000_TT,
            SIDM_ALDEBARAN_15TAU,
            (0.0, 0.0),
            0.0,
            false,
            true,
            true,
            some,
        )
        .unwrap();
        assert_eq!(
            stub.seen.borrow()[0],
            (DynamicAnchor::Aldebaran, false, false)
        );
        // Non-finite provider longitudes are refused, not wrapped.
        struct Bad;
        impl DynamicAnchorProvider for Bad {
            fn anchor_longitude_deg(
                &self,
                _: DynamicAnchor,
                _: f64,
                _: bool,
                _: bool,
            ) -> Result<f64, Error> {
                Ok(f64::NAN)
            }
        }
        let bad = Bad;
        assert!(
            calc_mean_ayanamsha_deg(
                J2000_TT,
                SIDM_TRUE_CITRA,
                (0.0, 0.0),
                0.0,
                false,
                false,
                false,
                Some(&bad as &dyn DynamicAnchorProvider),
            )
            .is_err()
        );
        assert!(
            calc_mean_ayanamsha_deg(f64::NAN, 1, (0.0, 0.0), 0.0, false, false, false, None,)
                .is_err()
        );
    }

    #[test]
    fn provider_less_serving_refuses_live_anchor_ayanamsha() {
        // Scalar-routes disposition: the production provider-less serving
        // construction (scalar comparison arms, scalar C-ABI entry) refuses
        // every live star/galactic arm with `AyanamshaAnchorUnserved` —
        // never a silent epoch substitution — while the formula/epoch arms
        // stay servable. Provider-bound sealed servings evaluate the live
        // arms; the swept ayanamsha manifests exclude live modes by
        // generator construction, so this refusal class never meets a
        // served the model value in comparison.
        let dt = |_: f64| -> Result<f64, Error> { Ok(69.0 / 86400.0) };
        for mode in DYNAMIC_AYANAMSHA_MODES {
            let mut sel = SiderealSelection::unset();
            sel.set_sidereal_mode(mode, 0.0, 0.0);
            let bound = SealedAyanamshaServing::new(sel, &dt, None);
            let out = bound.mean_ut(J2000_TT);
            if mode == SIDM_GALALIGN_MARDYKS || mode == SIDM_VALENS_MOON {
                assert!(out.is_ok(), "mode={mode}");
            } else {
                assert_eq!(out, Err(Error::AyanamshaAnchorUnserved { mode }));
            }
        }
        // Epoch-arm contrast on the same provider-less serving.
        let mut sel = SiderealSelection::unset();
        sel.set_sidereal_mode(1, 0.0, 0.0);
        let bound = SealedAyanamshaServing::new(sel, &dt, None);
        assert!(bound.mean_ut(J2000_TT).is_ok());
    }

    #[test]
    fn serving_timescale_and_user_ut_rules() {
        // ΔT stand-in: constant 69 s, so jd_tt = jd_ut + 69/86400.
        let dt = |_: f64| -> Result<f64, Error> { Ok(69.0 / 86400.0) };
        let mut sel = SiderealSelection::unset();
        sel.set_sidereal_mode(1, 0.0, 0.0);
        let bound = SealedAyanamshaServing::new(sel, &dt, None);
        assert!(
            (bound.jd_tt_for_ut(2451545.0).unwrap() - (2451545.0 + 69.0 / 86400.0)).abs() < 1e-15
        );
        assert!(
            (bound.ut1_for_tt(2451545.0).unwrap() - (2451545.0 - 69.0 / 86400.0)).abs() < 1e-15
        );
        assert!(
            (bound.valens_epoch_tt().unwrap() - (VALENS_MOON_T0_UT + 69.0 / 86400.0)).abs() < 1e-15
        );
        // mean_ut converts then dispatches: Lahiri at J2000 TT reads
        // the defining value up to the 69 s accumulation (~1 e-7 deg).
        let (a, w) = bound.mean_ut(2451545.0 - 69.0 / 86400.0).unwrap();
        assert!(w.is_none());
        assert!((a - 23.85709166666667).abs() < 1e-6, "a={a}");
        // SIDM_USER + SIDBIT_USER_UT interprets t 0 as a UT date.
        let mut sel_ut = SiderealSelection::unset();
        sel_ut.set_sidereal_mode(SIDM_USER | SIDBIT_USER_UT, 2451545.0, 24.0);
        let bound_ut = SealedAyanamshaServing::new(sel_ut, &dt, None);
        let (t0_tt, _) = bound_ut.user_pair_tt().unwrap();
        assert!((t0_tt - (2451545.0 + 69.0 / 86400.0)).abs() < 1e-15);
        // ... while without the bit the same values stay TT.
        let mut sel_tt = SiderealSelection::unset();
        sel_tt.set_sidereal_mode(SIDM_USER, 2451545.0, 24.0);
        let bound_tt = SealedAyanamshaServing::new(sel_tt, &dt, None);
        let (t0_raw, _) = bound_tt.user_pair_tt().unwrap();
        assert_eq!(t0_raw, 2451545.0);
        // ex entries echo the retflag and select true vs mean.
        let (ret, v_true, _) = bound.ex_ut(2451545.0, 0).unwrap();
        assert_eq!(ret, crate::constants::FLG_SWIEPH);
        let (ret_n, v_mean, _) = bound.ex_ut(2451545.0, crate::constants::FLG_NONUT).unwrap();
        assert_eq!(
            ret_n,
            crate::constants::FLG_SWIEPH | crate::constants::FLG_NONUT
        );
        let d = (v_true - v_mean).abs();
        assert!(d < 0.01 && d > 1e-6, "d={d}");
        // for_flags_ut honors the NONUT-mean rule.
        let (f_mean, _) = bound
            .for_flags_ut(2451545.0, crate::constants::FLG_NONUT)
            .unwrap();
        assert!((f_mean - v_mean).abs() < 1e-15);
        // Hostile inputs refused; ΔT failures propagate.
        assert!(bound.mean_ut(f64::NAN).is_err());
        assert!(bound.ut1_for_tt(f64::INFINITY).is_err());
        let bad_dt = |_: f64| -> Result<f64, Error> { Ok(f64::NAN) };
        let bad_bound = SealedAyanamshaServing::new(sel, &bad_dt, None);
        assert!(bad_bound.mean_ut(2451545.0).is_err());
    }

    #[test]
    fn mean_tt_dispatch_keeps_unused_epoch_providers_lazy() {
        let refuse = |_: f64| -> Result<f64, Error> {
            panic!("unused defining epoch must not consume Delta-T")
        };
        for mode in [
            SIDM_LAHIRI | SIDBIT_ECL_T0 | SIDBIT_USER_UT,
            SIDM_GALALIGN_MARDYKS | SIDBIT_USER_UT,
            999 | SIDBIT_USER_UT,
        ] {
            let mut selection = SiderealSelection::unset();
            selection.set_sidereal_mode(mode, 1e308, 0.0);
            let serving = SealedAyanamshaServing::new(selection, &refuse, None);
            assert!(serving.mean_tt_with(J2000_TT, false, false).is_ok());
        }
        let mut live = SiderealSelection::unset();
        live.set_sidereal_mode(SIDM_TRUE_CITRA | SIDBIT_USER_UT, 1e308, 0.0);
        let serving = SealedAyanamshaServing::new(live, &refuse, None);
        assert_eq!(
            serving.mean_tt_with(J2000_TT, false, false),
            Err(Error::AyanamshaAnchorUnserved {
                mode: SIDM_TRUE_CITRA
            })
        );
        assert!(matches!(
            serving.mean_tt_with(f64::NAN, false, false),
            Err(Error::NonFiniteArgument { .. })
        ));
    }

    #[test]
    fn mean_tt_dispatch_acquires_only_the_selected_epoch_clock() {
        for (mode, expected) in [
            (SIDM_USER | SIDBIT_USER_UT, vec![J2000_TT]),
            (SIDM_USER, vec![]),
            (SIDM_VALENS_MOON | SIDBIT_USER_UT, vec![VALENS_MOON_T0_UT]),
        ] {
            let calls = std::cell::RefCell::new(Vec::new());
            let provider = |ut| {
                calls.borrow_mut().push(ut);
                Ok(69.0 / 86400.0)
            };
            let mut selection = SiderealSelection::unset();
            selection.set_sidereal_mode(mode, J2000_TT, 24.0);
            let serving = SealedAyanamshaServing::new(selection, &provider, None);
            assert!(serving.mean_tt_with(J2000_TT, false, false).is_ok());
            assert_eq!(*calls.borrow(), expected);
        }
        let provider = |_: f64| Ok(f64::INFINITY);
        for mode in [SIDM_USER | SIDBIT_USER_UT, SIDM_VALENS_MOON] {
            let mut selection = SiderealSelection::unset();
            selection.set_sidereal_mode(mode, J2000_TT, 24.0);
            let serving = SealedAyanamshaServing::new(selection, &provider, None);
            assert!(matches!(
                serving.mean_tt_with(J2000_TT, false, false),
                Err(Error::NonFiniteArgument { value }) if value.is_infinite()
            ));
        }
    }

    #[test]
    fn sidereal_frame_product_layout_identity_and_refusals() {
        let identity = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        let a = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 10.0]];
        let b = [[2.0, 0.0, 1.0], [3.0, 4.0, 0.0], [0.0, 5.0, 6.0]];
        assert_eq!(sidereal_frame_product(a, identity).unwrap(), a);
        assert_eq!(sidereal_frame_product(identity, a).unwrap(), a);
        assert_eq!(
            sidereal_frame_product(a, b).unwrap(),
            [[8.0, 23.0, 19.0], [23.0, 50.0, 40.0], [38.0, 82.0, 67.0]]
        );
        assert_eq!(
            sidereal_frame_product_transposed_right(a, identity).unwrap(),
            a
        );
        assert_eq!(
            sidereal_frame_product_transposed_right(identity, b).unwrap(),
            transpose3(b)
        );
        assert_eq!(
            sidereal_frame_product_transposed_right(a, b).unwrap(),
            [[5.0, 11.0, 28.0], [14.0, 32.0, 61.0], [24.0, 53.0, 100.0]]
        );
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            let mut bad = a;
            bad[1][2] = value;
            assert!(sidereal_frame_product(bad, b).is_err());
            assert!(sidereal_frame_product(a, bad).is_err());
            assert!(sidereal_frame_product_transposed_right(bad, b).is_err());
            assert!(sidereal_frame_product_transposed_right(a, bad).is_err());
        }
        let zero = sidereal_frame_product(a, [[0.0; 3]; 3]).unwrap();
        assert_eq!(zero, [[0.0; 3]; 3]);
        let zero = sidereal_frame_product_transposed_right(a, [[0.0; 3]; 3]).unwrap();
        assert_eq!(zero, [[0.0; 3]; 3]);
    }

    #[test]
    fn native_epoch_frames_preserve_vector_norm_and_roundtrip() {
        let vector = [0.3, -0.4, 0.5];
        let square_norm = |value: [f64; 3]| value.iter().map(|word| word * word).sum::<f64>();
        for epoch in [
            1500000.5,
            2378496.5,
            2433282.423,
            J2000_TT,
            2460310.5,
            3000000.5,
        ] {
            for frame in [
                equatorial_epoch_matrix(epoch).unwrap(),
                ecliptic_of_t0_matrix(epoch).unwrap(),
            ] {
                let projected = crate::positions::mat3_vec3(frame, vector).unwrap();
                assert!((square_norm(projected) - square_norm(vector)).abs() < 2e-13);
                let restored = crate::positions::mat3_vec3(transpose3(frame), projected).unwrap();
                for (expected, actual) in vector.into_iter().zip(restored) {
                    assert!(
                        (actual - expected).abs() < 2e-13,
                        "epoch={epoch}: {actual} != {expected}"
                    );
                }
            }
        }
    }

    // Position projection: analytic properties.

    #[test]
    fn base_offset_wrap_and_refusals() {
        // Lahiri-like subtraction wraps into [0, 360).
        let lon = apply_sidereal_offset_deg(10.0, 23.85709166666667).unwrap();
        let want = (10.0 - 23.85709166666667_f64).rem_euclid(360.0);
        assert!((lon - want).abs() < 1e-12, "lon={lon}");
        let lon2 = apply_sidereal_offset_deg(5.0, 10.0).unwrap();
        assert!((lon2 - 355.0).abs() < 1e-12, "lon2={lon2}");
        // Undefined lunar points pass through; provider faults refuse.
        assert!(apply_sidereal_offset_deg(f64::NAN, 23.0).unwrap().is_nan());
        assert!(apply_sidereal_offset_deg(10.0, f64::NAN).is_err());
        assert!(apply_sidereal_offset_deg(10.0, f64::INFINITY).is_err());
        // 0/360 straddle unwraps instead of spiking the speed.
        let drift = sidereal_rate_correction_deg_day(359.9999, 0.0001).unwrap();
        assert!(drift.abs() < 100.0, "drift={drift}");
        assert!(drift > 0.0, "drift={drift}");
        assert!(sidereal_rate_correction_deg_day(f64::NAN, 0.0).is_err());
        // Full observer-frame correction composes both halves.
        let (l, dl) =
            apply_sidereal_correction(10.0, 1.0, 23.85709166666667, 23.857, 23.8572, true).unwrap();
        assert!((l - want).abs() < 1e-12, "l={l}");
        // The 0.0002° sample drift over 2 s removes ~8.64°/day from dlon.
        let want_dl = 1.0 - 0.0002 / (2.0 * SIDEREAL_SPEED_HALF_STEP_DAYS);
        assert!((dl - want_dl).abs() < 1e-9, "dl={dl}");
        let (ln, _) = apply_sidereal_correction(f64::NAN, 0.0, 1.0, 1.0, 1.0, true).unwrap();
        assert!(ln.is_nan());
    }

    #[test]
    fn degenerate_sentinel_contract() {
        assert!(is_degenerate_origin_6tuple(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]));
        assert!(is_degenerate_origin_6tuple(&[
            -0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        ]));
        assert!(!is_degenerate_origin_6tuple(&[
            1.0, 0.0, 0.0, 0.0, 0.0, 0.0
        ]));
        assert!(!is_degenerate_origin_6tuple(&[
            0.0, 0.0, 1e-300, 0.0, 0.0, 0.0
        ]));
    }

    #[test]
    fn fixed_epoch_flag_helpers() {
        use crate::constants::{
            FLG_EQUATORIAL, FLG_J2000, FLG_NONUT, FLG_RADIANS, FLG_SIDEREAL, FLG_SPEED, FLG_SWIEPH,
        };
        // Request detection needs both the bit and a fixed mode.
        assert!(is_fixed_epoch_request(FLG_SIDEREAL, SIDM_J2000));
        assert!(is_fixed_epoch_request(FLG_SIDEREAL, SIDM_GALALIGN_MARDYKS));
        assert!(!is_fixed_epoch_request(FLG_SIDEREAL, SIDM_LAHIRI));
        assert!(!is_fixed_epoch_request(0, SIDM_J2000));
        // Rewrite strips SIDEREAL/RADIANS, forces J2000|NONUT.
        let rw = fixed_epoch_request_flags(FLG_SIDEREAL | FLG_RADIANS | FLG_SPEED | FLG_SWIEPH);
        assert_eq!(rw, FLG_SPEED | FLG_SWIEPH | FLG_J2000 | FLG_NONUT);
        // Echo keeps the caller bits plus NONUT, drops the rewrite J2000.
        let sub = FLG_J2000 | FLG_NONUT | FLG_SPEED | FLG_SWIEPH;
        let req = FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH;
        assert_eq!(
            fixed_epoch_retflag(sub, req),
            FLG_NONUT | FLG_SPEED | FLG_SWIEPH | FLG_SIDEREAL
        );
        // An explicitly requested J2000 survives the echo.
        let req_j = req | FLG_J2000;
        assert_eq!(fixed_epoch_retflag(sub, req_j) & FLG_J2000, FLG_J2000);
        // Equatorial bit is inert on the flag helpers (routing lives in
        // the transforms).
        assert!(is_fixed_epoch_request(
            FLG_SIDEREAL | FLG_EQUATORIAL,
            SIDM_J1900
        ));
    }

    #[test]
    fn fixed_epoch_t0_table() {
        assert_eq!(fixed_epoch_t0_jd(SIDM_J2000), Some(J2000_TT));
        // J 1900/B 1950 go through the fused native converters, so
        // the values are identical to `epj2jd`/`epb2jd`
        // sums (ephemeral probe 2026-09-23 against;
        // scratch only). The former 1 e-6 bound is tightened to exactly.
        assert_eq!(
            fixed_epoch_t0_jd(SIDM_J1900).map(f64::to_bits),
            Some(2415020.0f64.to_bits())
        );
        assert_eq!(
            fixed_epoch_t0_jd(SIDM_B1950).map(f64::to_bits),
            Some(0x4142_9081_3633_e7f9u64)
        );
        assert_eq!(
            fixed_epoch_t0_jd(SIDM_GALALIGN_MARDYKS),
            Some(MARDYKS_DEFINING.1)
        );
        assert_eq!(fixed_epoch_t0_jd(SIDM_LAHIRI), None);
        assert_eq!(fixed_epoch_t0_jd(-1), None);
    }

    #[test]
    fn epoch_matrices_identity_and_orthogonality() {
        // J2000 is the identity pair.
        let (eq, ecl) = fixed_epoch_matrices(SIDM_J2000).unwrap();
        let eye = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        assert_eq!((eq, ecl), (eye, eye));
        // Other modes give proper rotations (orthogonal, det +1).
        for mode in [SIDM_J1900, SIDM_B1950, SIDM_GALALIGN_MARDYKS] {
            let (m_eq, m_ecl) = fixed_epoch_matrices(mode).unwrap();
            for m in [m_eq, m_ecl] {
                let mt = transpose3(m);
                let prod = crate::positions::mat3_mul(m, mt).unwrap();
                for (i, row) in prod.iter().enumerate() {
                    for (j, got) in row.iter().enumerate() {
                        let want = if i == j { 1.0 } else { 0.0 };
                        assert!(
                            (got - want).abs() < 1e-12,
                            "mode={mode}: prod[{i}][{j}]={got}"
                        );
                    }
                }
                let det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
                    - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
                    + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
                assert!((det - 1.0).abs() < 1e-12, "mode={mode}: det={det}");
            }
        }
        // Non-fixed modes refuse with a typed error, never a silent frame.
        assert!(fixed_epoch_matrices(SIDM_LAHIRI).is_err());
    }

    #[test]
    fn rotate_spherical_identity_sentinel_and_z180() {
        let eye = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        // Identity round-trips a generic state to ~1 ulp.
        let xx = [123.456, -20.25, 1.5, 0.75, -0.01, 0.002];
        let back = rotate_spherical_6tuple(xx, eye).unwrap();
        for (a, b) in xx.iter().zip(back.iter()) {
            assert!((a - b).abs() < 1e-12, "a={a} b={b}");
        }
        // The sentinel passes through untouched (no manufactured longitude).
        assert_eq!(
            rotate_spherical_6tuple([0.0; 6], rot_z_passive(1.0)).unwrap(),
            [0.0; 6]
        );
        // Half-turn about z adds 180° of longitude, keeps latitude/rates.
        let z180 = rot_z_passive(std::f64::consts::PI);
        let out = rotate_spherical_6tuple([10.0, 20.0, 2.0, 1.5, 0.1, 0.01], z180).unwrap();
        assert!((out[0] - 190.0).abs() < 1e-9, "lon={}", out[0]);
        assert!((out[1] - 20.0).abs() < 1e-9, "lat={}", out[1]);
        assert!((out[2] - 2.0).abs() < 1e-12, "r={}", out[2]);
        assert!((out[3] - 1.5).abs() < 1e-9, "dlon={}", out[3]);
        // Zero-distance direction keeps its angles, restores zero radius.
        let zd = rotate_spherical_6tuple([45.0, 30.0, 0.0, 0.0, 0.0, 0.0], eye).unwrap();
        assert!((zd[0] - 45.0).abs() < 1e-9 && (zd[1] - 30.0).abs() < 1e-9);
        assert_eq!(zd[2], 0.0);
        // Undefined lunar points pass through.
        let nan = rotate_spherical_6tuple([f64::NAN, 0.0, 1.0, 0.0, 0.0, 0.0], eye).unwrap();
        assert!(nan[0].is_nan());
    }

    #[test]
    fn invariable_plane_and_ssy_zero_point_properties() {
        // Proper rotation.
        let m = invariable_plane_matrix();
        let prod = crate::positions::mat3_mul(m, transpose3(m)).unwrap();
        for (i, row) in prod.iter().enumerate() {
            for (j, got) in row.iter().enumerate() {
                let want = if i == j { 1.0 } else { 0.0 };
                assert!((got - want).abs() < 1e-12, "prod[{i}][{j}]={got}");
            }
        }
        // Zero point lands in [0, 360) and moves monotonically locally.
        let z0 = ssy_plane_zero_point_deg(23.85709166666667).unwrap();
        assert!((0.0..360.0).contains(&z0), "z0={z0}");
        let z1 = ssy_plane_zero_point_deg(24.85709166666667).unwrap();
        let dz = (z1 - z0).rem_euclid(360.0);
        assert!(dz > 0.9 && dz < 1.1, "dz={dz}");
        assert!(ssy_plane_zero_point_deg(f64::NAN).is_err());
    }

    #[test]
    fn ecl_t0_epoch_resolution() {
        // USER/Valens passthrough; classical table; defining table; J2000.
        assert_eq!(ecl_t0_epoch_jd(SIDM_USER, 2400000.5, 0.0), 2400000.5);
        assert_eq!(ecl_t0_epoch_jd(SIDM_VALENS_MOON, 0.0, 1775846.0), 1775846.0);
        assert_eq!(ecl_t0_epoch_jd(SIDM_LAHIRI, 0.0, 0.0), 2435553.5);
        assert_eq!(ecl_t0_epoch_jd(SIDM_HIPPARCHOS, 0.0, 0.0), 1674484.0);
        assert_eq!(
            ecl_t0_epoch_jd(SIDM_FAGAN_BRADLEY, 0.0, 0.0),
            2433282.42345905
        );
        assert_eq!(ecl_t0_epoch_jd(47, 0.0, 0.0), J2000_TT);
        assert_eq!(ECL_T0_CLASSICAL_EPOCHS.len(), 9);
    }

    #[test]
    fn ecl_t0_classical_table_verbatim() {
        // every row of the table pinned exactly (read file-scoped ).
        // The triplet shares one projection plane (one Babylonian zodiac).
        const EXPECTED: [(i32, f64); 9] = [
            (SIDM_LAHIRI, 2435553.5),
            (SIDM_BABYL_KUGLER1, 1684532.5),
            (SIDM_BABYL_KUGLER2, 1684532.5),
            (SIDM_BABYL_KUGLER3, 1684532.5),
            (SIDM_HIPPARCHOS, 1674484.0),
            (SIDM_DELUCE, 1721057.5),
            (SIDM_DJWHAL_KHUL, 2415020.0),
            (SIDM_BABYL_BRITTON, 1721057.5),
            (SIDM_LAHIRI_1940, 2415020.0),
        ];
        assert_eq!(ECL_T0_CLASSICAL_EPOCHS, EXPECTED);
        for (mode, epoch) in EXPECTED {
            assert_eq!(ecl_t0_epoch_jd(mode, 0.0, 0.0), epoch);
        }
    }

    #[test]
    fn ecl_t0_epoch_boundary_defaults() {
        // modes with no classical row and no defining row — the
        // unassigned block, `SIDM_USER`-adjacent 254, and negatives — fall
        // back to J2000 like the model `.get` chain. USER/Valens never
        // reach that chain (caller epochs win, even pathological ones:
        // `t0 = 0.0` really means JD 0, NaN passes through).
        for mode in [47, 254, -1, i32::MIN, i32::MAX] {
            assert_eq!(ecl_t0_epoch_jd(mode, 0.0, 0.0), J2000_TT, "mode={mode}");
        }
        assert_eq!(ecl_t0_epoch_jd(SIDM_USER, 0.0, 0.0), 0.0);
        assert!(ecl_t0_epoch_jd(SIDM_USER, f64::NAN, 1775846.0).is_nan());
    }

    #[test]
    fn route_resolution_matches_release_branches() {
        use crate::constants::{FLG_EQUATORIAL, FLG_SIDEREAL};
        // Tropical requests never route.
        assert_eq!(resolve_sidereal_route(0, SIDM_LAHIRI, 0), None);
        // Fixed-epoch modes win over SIDBIT bits.
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL, SIDM_J1900, SIDBIT_ECL_T0),
            Some(SiderealRoute::FixedEpoch)
        );
        // SIDBIT ecliptic arm; ECL_T 0 precedence is in the matrix builder.
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL, SIDM_LAHIRI, SIDBIT_ECL_T0),
            Some(SiderealRoute::SidbitEcliptic)
        );
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL, SIDM_LAHIRI, SIDBIT_SSY_PLANE),
            Some(SiderealRoute::SidbitEcliptic)
        );
        // Equatorial goes to the T 0 arm only under ECL_T 0...
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL | FLG_EQUATORIAL, SIDM_LAHIRI, SIDBIT_ECL_T0),
            Some(SiderealRoute::SidbitEquatorialT0)
        );
        // ... SSY_PLANE alone leaves the equator on the base path...
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL | FLG_EQUATORIAL, SIDM_LAHIRI, SIDBIT_SSY_PLANE),
            Some(SiderealRoute::BaseCorrection)
        );
        // ... and the suppress-listed live-direction modes stay on base.
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL, SIDM_TRUE_CITRA, SIDBIT_ECL_T0),
            Some(SiderealRoute::BaseCorrection)
        );
        // Aldebaran (absolute anchor) and Valens (epoch mode) accept it.
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL, SIDM_ALDEBARAN_15TAU, SIDBIT_ECL_T0),
            Some(SiderealRoute::SidbitEcliptic)
        );
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL, SIDM_VALENS_MOON, SIDBIT_SSY_PLANE),
            Some(SiderealRoute::SidbitEcliptic)
        );
        // Plain epoch modes correct the scalar.
        assert_eq!(
            resolve_sidereal_route(FLG_SIDEREAL, SIDM_LAHIRI, 0),
            Some(SiderealRoute::BaseCorrection)
        );
    }

    #[test]
    fn sidbit_matrix_selection_and_precedence() {
        // No bit → no matrix.
        assert!(sidbit_ecliptic_matrix(0, J2000_TT, 0.0).unwrap().is_none());
        // Both bits → ECL_T 0 wins (identical to ECL_T 0 alone).
        let t0 = 2435553.5;
        let both = sidbit_ecliptic_matrix(SIDBIT_ECL_T0 | SIDBIT_SSY_PLANE, t0, 23.5)
            .unwrap()
            .unwrap();
        let ecl = sidbit_ecliptic_matrix(SIDBIT_ECL_T0, t0, 23.5)
            .unwrap()
            .unwrap();
        assert_eq!(both, ecl);
        // At the J2000 epoch the ECL_T 0 base is the identity up to the
        // Vondrak-chain rounding (~1 e-16 per element through the
        // obliquity sandwich), so the matrix is the zero-point
        // z-rotation within rotation-rounding class.
        let m = sidbit_ecliptic_matrix(SIDBIT_ECL_T0, J2000_TT, 30.0)
            .unwrap()
            .unwrap();
        let want = rot_z_passive(30.0_f64.to_radians());
        for (row_m, row_w) in m.iter().zip(want.iter()) {
            for (got, w) in row_m.iter().zip(row_w.iter()) {
                assert!((got - w).abs() < 1e-12, "got={got} want={w}");
            }
        }
        assert!(sidbit_ecliptic_matrix(SIDBIT_ECL_T0, f64::NAN, 0.0).is_err());
    }

    /// Stub tropical core: echoes a scripted tuple and records the flags.
    struct StubTropical {
        xx: [f64; 6],
        seen: std::cell::RefCell<Vec<i32>>,
    }

    impl StubTropical {
        fn call(
            &self,
            _jd: f64,
            _body: i32,
            flags: i32,
        ) -> Result<(crate::positions::PositionResult, i32), Error> {
            self.seen.borrow_mut().push(flags);
            Ok((self.xx, flags))
        }
    }

    /// Bind a [`SiderealFastServing`] over stub closures.
    ///
    /// The serving borrows its providers, so the closures live in the
    /// calling test's scope and are passed in by reference.
    #[allow(clippy::too_many_arguments)]
    fn stub_serve<'a>(
        sel: SiderealSelection,
        tropical: &'a dyn Fn(
            f64,
            i32,
            i32,
        ) -> Result<(crate::positions::PositionResult, i32), Error>,
        mean_tt: &'a dyn Fn(f64, bool, bool) -> Result<f64, Error>,
        mean_ut: &'a dyn Fn(f64) -> Result<f64, Error>,
        dpsi: &'a dyn Fn(f64) -> Result<f64, Error>,
        rate: &'a dyn Fn(f64) -> Result<f64, Error>,
    ) -> SiderealFastServing<'a> {
        SiderealFastServing::new(
            sel, tropical, mean_tt, mean_ut, dpsi, rate, 2451545.0, 1775846.0,
        )
    }

    #[test]
    fn native_epoch_pair_rate_preserves_acquisition_order_and_opaque_defaults() {
        use std::cell::RefCell;
        let calls = RefCell::new(Vec::new());
        let tropical = |_: f64, _: i32, flags: i32| {
            calls.borrow_mut().push("tropical");
            Ok(([100.0, 5.0, 1.2, 0.9, 0.01, 0.001], flags))
        };
        let mean = |_: f64, _: bool, _: bool| {
            calls.borrow_mut().push("mean");
            Ok(23.0)
        };
        let nutation = |_: f64| {
            calls.borrow_mut().push("nutation");
            Ok(0.004)
        };
        let nutation_rate = |_: f64| {
            calls.borrow_mut().push("nutation_rate");
            Ok(0.002)
        };
        let epoch_rate = |tt: f64, value: f64, epoch: f64, ecl_date: bool| {
            calls.borrow_mut().push("epoch_rate");
            assert_eq!(
                (tt, value, epoch, ecl_date),
                (J2000_TT, AYANAMSHA_DEFINING[1].1, J2000_TT, false)
            );
            Ok(0.003)
        };
        let mut selection = SiderealSelection::unset();
        selection.set_sidereal_mode(SIDM_LAHIRI, 0.0, 0.0);
        let flags = FLG_SIDEREAL | FLG_SPEED;
        let opaque = stub_serve(
            selection,
            &tropical,
            &mean,
            &|_| Ok(0.0),
            &nutation,
            &nutation_rate,
        );
        let expected = opaque.compute_tt(J2000_TT, 4, flags).unwrap();
        assert_eq!(
            *calls.borrow(),
            vec!["tropical", "mean", "nutation", "nutation_rate"]
        );
        assert_eq!(
            expected.0[3],
            0.9 - crate::positions::general_precession_rate_deg_day(J2000_TT).unwrap() - 0.002
        );
        calls.borrow_mut().clear();
        let native = opaque.with_epoch_pair_rate(&epoch_rate);
        let actual = native.compute_tt(J2000_TT, 4, flags).unwrap();
        assert_eq!(
            *calls.borrow(),
            vec![
                "tropical",
                "mean",
                "nutation",
                "epoch_rate",
                "nutation_rate"
            ]
        );
        assert_eq!(actual.0[..3], expected.0[..3]);
        assert_eq!(actual.0[4..], expected.0[4..]);
        assert_eq!(actual.1, expected.1);
        assert_eq!(actual.0[3], 0.9 - 0.003 - 0.002);
        calls.borrow_mut().clear();
        native.compute_tt(J2000_TT, 4, flags & !FLG_SPEED).unwrap();
        assert_eq!(*calls.borrow(), vec!["tropical", "mean", "nutation"]);
        calls.borrow_mut().clear();
        native
            .compute_tt(J2000_TT, 4, flags | FLG_EQUATORIAL)
            .unwrap();
        assert_eq!(*calls.borrow(), vec!["tropical"]);
        calls.borrow_mut().clear();
        native.compute_tt(J2000_TT, 4, flags | FLG_NONUT).unwrap();
        assert_eq!(*calls.borrow(), vec!["tropical", "mean", "epoch_rate"]);
    }

    #[test]
    fn native_epoch_pair_rate_keeps_value_failures_ahead_of_rate() {
        use std::cell::RefCell;
        let calls = RefCell::new(Vec::new());
        let tropical = |_: f64, _: i32, flags: i32| Ok(([100.0, 0.0, 1.0, 1.0, 0.0, 0.0], flags));
        let epoch_rate = |_: f64, _: f64, _: f64, _: bool| {
            calls.borrow_mut().push("epoch_rate");
            Err(Error::MathDomain { value: 3.0 })
        };
        let nutation_rate = |_: f64| {
            calls.borrow_mut().push("nutation_rate");
            Ok(0.0)
        };
        let mean_error = |_: f64, _: bool, _: bool| Err(Error::MathDomain { value: 1.0 });
        let nutation_error = |_: f64| Err(Error::MathDomain { value: 2.0 });
        let flags = FLG_SIDEREAL | FLG_SPEED;
        let selection = SiderealSelection::unset();
        let first = stub_serve(
            selection,
            &tropical,
            &mean_error,
            &|_| Ok(0.0),
            &nutation_error,
            &nutation_rate,
        )
        .with_epoch_pair_rate(&epoch_rate);
        assert!(matches!(
            first.compute_tt(J2000_TT, 4, flags),
            Err(Error::MathDomain { value: 1.0 })
        ));
        assert!(calls.borrow().is_empty());
        let second = stub_serve(
            selection,
            &tropical,
            &|_, _, _| Ok(1.0),
            &|_| Ok(0.0),
            &nutation_error,
            &nutation_rate,
        )
        .with_epoch_pair_rate(&epoch_rate);
        assert!(matches!(
            second.compute_tt(J2000_TT, 4, flags),
            Err(Error::MathDomain { value: 2.0 })
        ));
        assert!(calls.borrow().is_empty());
        let third = stub_serve(
            selection,
            &tropical,
            &|_, _, _| Ok(1.0),
            &|_| Ok(0.0),
            &|_| Ok(0.0),
            &nutation_rate,
        )
        .with_epoch_pair_rate(&epoch_rate);
        assert!(matches!(
            third.compute_tt(J2000_TT, 4, flags),
            Err(Error::MathDomain { value: 3.0 })
        ));
        assert_eq!(*calls.borrow(), vec!["epoch_rate"]);
    }

    #[test]
    fn native_epoch_pair_rate_uses_acquired_user_and_valens_clocks() {
        use std::cell::RefCell;
        let calls = RefCell::new(Vec::new());
        let epoch_rate = |tt: f64, value: f64, epoch: f64, ecl_date: bool| {
            calls.borrow_mut().push((tt, value, epoch, ecl_date));
            Ok(0.0)
        };
        let tropical = |_: f64, _: i32, flags: i32| Ok(([100.0, 0.0, 1.0, 1.0, 0.0, 0.0], flags));
        for mode in [
            SIDM_USER | SIDBIT_USER_UT | SIDBIT_ECL_DATE,
            SIDM_VALENS_MOON,
        ] {
            let mut selection = SiderealSelection::unset();
            selection.set_sidereal_mode(mode, 123.0, 24.5);
            let serving = SiderealFastServing::new(
                selection,
                &tropical,
                &|_, _, _| Ok(2.0),
                &|_| Ok(0.0),
                &|_| Ok(0.0),
                &|_| Ok(0.0),
                124.25,
                1775845.501,
            )
            .with_epoch_pair_rate(&epoch_rate);
            serving
                .compute_tt(J2000_TT, 4, FLG_SIDEREAL | FLG_SPEED)
                .unwrap();
        }
        assert_eq!(
            *calls.borrow(),
            vec![
                (J2000_TT, 24.5, 124.25, true),
                (J2000_TT, VALENS_MOON_AYAN_T0_DEG, 1775845.501, false)
            ]
        );
    }

    #[test]
    fn native_epoch_pair_rate_never_replaces_live_anchor_sampling() {
        use std::cell::RefCell;
        let calls = RefCell::new(Vec::new());
        let tropical = |_: f64, _: i32, flags: i32| Ok(([100.0, 0.0, 1.0, 1.0, 0.0, 0.0], flags));
        let mean = |tt: f64, _: bool, _: bool| {
            calls.borrow_mut().push(tt.to_bits());
            Ok(20.0 + (tt - J2000_TT) * 0.001)
        };
        let forbidden = |_: f64, _: f64, _: f64, _: bool| -> Result<f64, Error> {
            panic!("live anchors do not expose an epoch-pair derivative")
        };
        for mode in DYNAMIC_AYANAMSHA_MODES {
            if mode == SIDM_VALENS_MOON || mode == SIDM_GALALIGN_MARDYKS {
                continue;
            }
            let mut selection = SiderealSelection::unset();
            selection.set_sidereal_mode(mode, 0.0, 0.0);
            let base = stub_serve(
                selection,
                &tropical,
                &mean,
                &|_| Ok(0.0),
                &|_| Ok(0.0),
                &|_| Ok(0.0),
            );
            calls.borrow_mut().clear();
            let expected = base
                .compute_tt(J2000_TT, 4, FLG_SIDEREAL | FLG_SPEED)
                .unwrap();
            let expected_calls = calls.borrow().clone();
            calls.borrow_mut().clear();
            let actual = base
                .with_epoch_pair_rate(&forbidden)
                .compute_tt(J2000_TT, 4, FLG_SIDEREAL | FLG_SPEED)
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(*calls.borrow(), expected_calls);
            assert_eq!(expected_calls.len(), 3);
        }
    }

    #[test]
    fn serving_base_correction_gates_and_mean_true_selection() {
        use crate::constants::{
            FLG_EQUATORIAL, FLG_J2000, FLG_NONUT, FLG_SIDEREAL, FLG_SPEED, FLG_SWIEPH,
        };
        let trop = StubTropical {
            xx: [100.0, 5.0, 1.2, 0.9, 0.01, 0.001],
            seen: std::cell::RefCell::new(Vec::new()),
        };
        let mut sel = SiderealSelection::unset();
        sel.set_sidereal_mode(SIDM_LAHIRI, 0.0, 0.0);
        let tropical = |jd: f64, body: i32, flags: i32| trop.call(jd, body, flags);
        let mean_tt = |_jd: f64, _nb: bool, _ng: bool| -> Result<f64, Error> { Ok(23.857) };
        let mean_ut = |_jd: f64| -> Result<f64, Error> { Ok(23.857) };
        let dpsi = |_jd: f64| -> Result<f64, Error> { Ok(0.0042) };
        let rate = |_jd: f64| -> Result<f64, Error> { Ok(1e-6) };
        let svc = stub_serve(sel, &tropical, &mean_tt, &mean_ut, &dpsi, &rate);
        // Of-date frame: true = mean + stub dpsi.
        let (out, rf) = svc
            .compute_tt(J2000_TT, 4, FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH)
            .unwrap();
        assert!(
            (out[0] - (100.0 - 23.857 - 0.0042)).abs() < 1e-9,
            "lon={}",
            out[0]
        );
        assert_eq!(rf, FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH);
        // Tropical sub-call keeps the original flags (SIDEREAL included
        // for the mean-equator selection downstream).
        assert_eq!(
            *trop.seen.borrow(),
            vec![FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH]
        );
        // J2000 frame: mean only, and the analytic drift still applies.
        let (out_j, _) = svc
            .compute_tt(
                J2000_TT,
                4,
                FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH | FLG_J2000,
            )
            .unwrap();
        assert!(
            (out_j[0] - (100.0 - 23.857)).abs() < 1e-9,
            "lon={}",
            out_j[0]
        );
        // NONUT likewise selects the mean.
        let (out_n, _) = svc
            .compute_tt(J2000_TT, 4, FLG_SIDEREAL | FLG_SWIEPH | FLG_NONUT)
            .unwrap();
        assert!((out_n[0] - (100.0 - 23.857)).abs() < 1e-9);
        // Equatorial output never subtracts; zero distance never subtracts.
        let (out_eq, _) = svc
            .compute_tt(
                J2000_TT,
                4,
                FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH | FLG_EQUATORIAL,
            )
            .unwrap();
        assert_eq!(out_eq, trop.xx);
        let trop0 = StubTropical {
            xx: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            seen: std::cell::RefCell::new(Vec::new()),
        };
        let tropical0 = |jd: f64, body: i32, flags: i32| trop0.call(jd, body, flags);
        let svc0 = stub_serve(sel, &tropical0, &mean_tt, &mean_ut, &dpsi, &rate);
        let (out0, _) = svc0
            .compute_tt(J2000_TT, 14, FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH)
            .unwrap();
        assert_eq!(out0, [0.0; 6]);
    }

    #[test]
    fn serving_fixed_epoch_rewrites_and_echoes() {
        use crate::constants::{FLG_J2000, FLG_NONUT, FLG_SIDEREAL, FLG_SPEED, FLG_SWIEPH};
        let trop = StubTropical {
            xx: [100.0, 5.0, 1.2, 0.9, 0.01, 0.001],
            seen: std::cell::RefCell::new(Vec::new()),
        };
        let mut sel = SiderealSelection::unset();
        sel.set_sidereal_mode(SIDM_J1900, 0.0, 0.0);
        let tropical = |jd: f64, body: i32, flags: i32| trop.call(jd, body, flags);
        let mean_tt = |_jd: f64, _nb: bool, _ng: bool| -> Result<f64, Error> { Ok(0.0) };
        let mean_ut = |_jd: f64| -> Result<f64, Error> { Ok(0.0) };
        let dpsi = |_jd: f64| -> Result<f64, Error> { Ok(0.0042) };
        let rate = |_jd: f64| -> Result<f64, Error> { Ok(1e-6) };
        let svc = stub_serve(sel, &tropical, &mean_tt, &mean_ut, &dpsi, &rate);
        let flags = FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH;
        let (out, rf) = svc.compute_tt(J2000_TT, 4, flags).unwrap();
        // Sub-request ran rewritten (SIDEREAL stripped, J2000|NONUT set).
        assert_eq!(
            *trop.seen.borrow(),
            vec![FLG_SPEED | FLG_SWIEPH | FLG_J2000 | FLG_NONUT]
        );
        // Output moved to the t 0 frame and echoes SIDEREAL+NONUT.
        assert!((out[0] - 100.0).abs() > 0.01, "lon={}", out[0]);
        assert_eq!(rf & FLG_SIDEREAL, FLG_SIDEREAL);
        assert_eq!(rf & FLG_NONUT, FLG_NONUT);
        assert_eq!(rf & FLG_J2000, 0);
        // J2000 mode is the near-identity (input values preserved).
        let mut sel2 = SiderealSelection::unset();
        sel2.set_sidereal_mode(SIDM_J2000, 0.0, 0.0);
        let svc2 = stub_serve(sel2, &tropical, &mean_tt, &mean_ut, &dpsi, &rate);
        let (out2, _) = svc2.compute_tt(J2000_TT, 4, flags).unwrap();
        for (a, b) in trop.xx.iter().zip(out2.iter()) {
            assert!((a - b).abs() < 1e-9, "a={a} b={b}");
        }
    }

    #[test]
    fn serving_sidbit_arms_route_and_compose() {
        use crate::constants::{FLG_EQUATORIAL, FLG_NONUT, FLG_SIDEREAL, FLG_SPEED, FLG_SWIEPH};
        let trop = StubTropical {
            xx: [100.0, 5.0, 1.2, 0.9, 0.01, 0.001],
            seen: std::cell::RefCell::new(Vec::new()),
        };
        let mut sel = SiderealSelection::unset();
        sel.set_sidereal_mode(SIDM_LAHIRI | SIDBIT_ECL_T0, 0.0, 0.0);
        assert_eq!(sel.bits(), SIDBIT_ECL_T0);
        let tropical = |jd: f64, body: i32, flags: i32| trop.call(jd, body, flags);
        let mean_tt = |_jd: f64, _nb: bool, _ng: bool| -> Result<f64, Error> { Ok(23.857) };
        let mean_ut = |_jd: f64| -> Result<f64, Error> { Ok(23.857) };
        let dpsi = |_jd: f64| -> Result<f64, Error> { Ok(0.0042) };
        let rate = |_jd: f64| -> Result<f64, Error> { Ok(1e-6) };
        let svc = stub_serve(sel, &tropical, &mean_tt, &mean_ut, &dpsi, &rate);
        let flags = FLG_SIDEREAL | FLG_SPEED | FLG_SWIEPH;
        let (out, rf) = svc.compute_tt(J2000_TT, 4, flags).unwrap();
        // Ecliptic arm keeps SIDEREAL+NONUT, drops the rewrite J2000.
        assert_eq!(rf & FLG_SIDEREAL, FLG_SIDEREAL);
        assert_eq!(rf & FLG_NONUT, FLG_NONUT);
        assert_eq!(rf & crate::constants::FLG_J2000, 0);
        assert!((out[0] - 100.0).abs() > 1e-6, "lon={}", out[0]);
        // Equatorial ECL_T 0 arm lands in the t 0 mean frame.
        let (out_eq, _) = svc.compute_tt(J2000_TT, 4, flags | FLG_EQUATORIAL).unwrap();
        assert!((out_eq[0] - 100.0).abs() > 1e-6);
        // SSY_PLANE equatorial stays on the base path (tropical values).
        let mut sel_s = SiderealSelection::unset();
        sel_s.set_sidereal_mode(SIDM_LAHIRI | SIDBIT_SSY_PLANE, 0.0, 0.0);
        let svc_s = stub_serve(sel_s, &tropical, &mean_tt, &mean_ut, &dpsi, &rate);
        let (out_s, _) = svc_s
            .compute_tt(J2000_TT, 4, flags | FLG_EQUATORIAL)
            .unwrap();
        assert_eq!(out_s, trop.xx);
        // Zero-point helper refuses off-arm routes instead of guessing.
        assert!(
            svc.projection_zero_point_deg(SiderealRoute::BaseCorrection)
                .is_err()
        );
        // Non-finite epochs refuse before any provider contact.
        assert!(svc.compute_tt(f64::NAN, 4, flags).is_err());
    }
}

#[cfg(test)]
mod live_anchor_rate_dispatch_tests {
    use super::*;
    struct Toggles;
    impl DynamicAnchorProvider for Toggles {
        fn anchor_longitude_deg(
            &self,
            _: DynamicAnchor,
            _: f64,
            _: bool,
            _: bool,
        ) -> Result<f64, Error> {
            panic!("rate dispatch must not sample scalar anchors")
        }
        fn anchor_rate_deg_day(
            &self,
            anchor: DynamicAnchor,
            _: f64,
            a: bool,
            b: bool,
        ) -> Result<Option<f64>, Error> {
            if anchor == DynamicAnchor::Aldebaran {
                assert!(!a && !b)
            } else {
                assert!(a && b)
            }
            Ok(Some(0.125))
        }
    }
    #[test]
    fn live_anchor_rate_keeps_aldedebaran_toggle_neutrality() {
        let delta = |_| Ok(0.0);
        let anchors = Toggles;
        for mode in [SIDM_ALDEBARAN_15TAU, SIDM_TRUE_CITRA] {
            let mut selection = SiderealSelection::unset();
            selection.set_sidereal_mode(mode, 0.0, 0.0);
            let serving = SealedAyanamshaServing::new(selection, &delta, Some(&anchors));
            assert_eq!(
                serving
                    .live_mean_rate_deg_day(2451545.0, true, true)
                    .unwrap(),
                Some(0.125)
            );
        }
    }
}
