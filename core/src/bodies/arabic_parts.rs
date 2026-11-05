//! Arabic lots, zodiac partitions, dignities, aspects and civil string helpers.
//!
//! Longitudes and geographic angles use degrees. Pure kernels evaluate in
//! binary64 and do not load ephemeris files. The four lots combine ascendant
//! and planetary longitudes, with the documented day/night sect reversal.
//! Provider-dependent compositions live in `api::arabic_parts` and
//! `api::contrib`. ISO timestamps use proleptic Gregorian civil dates and
//! retain explicit UTC offsets and subsecond precision for durations.
//!
//! Domain and parsing failures are typed [`ArabicError`] values. Each entry
//! documents any finite-input, latitude or calendar restriction.

use std::fmt;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// Typed failures for the Arabic-parts / contrib pure kernels.
#[derive(Debug, Clone, PartialEq)]
pub enum ArabicError {
    /// A floating-point input was NaN or infinite. `field` names the
    /// argument (`"asc"`, `"sun"`, `"lon"`,...).
    NonFinite {
        /// Name of the rejected floating-point argument.
        field: &'static str,
        /// Non-finite input value, in the units of the named argument.
        value: f64,
    },
    /// Geographic latitude outside `[-90, 90]`.
    InvalidLatitude {
        /// Rejected geographic latitude in degrees.
        value: f64,
    },
    /// Geographic longitude outside `[-180, 360]`.
    InvalidLongitude {
        /// Rejected geographic longitude in degrees East.
        value: f64,
    },
    /// The 3-D high-latitude day/night arm was requested (`|lat| > 60°`
    /// with JD + observer supplied) but the `azalt`/`ECL2HOR` geometric
    /// pipeline is outside this lane. the model flags the 2-D answer
    /// as unreliable there, so this fails instead of guessing.
    HighLatitude3dUnavailable {
        /// Observer latitude in degrees requiring the unavailable 3-D sect calculation.
        lat: f64,
        /// Requested observer epoch as a UT1 Julian day.
        jd: f64,
    },
    /// A coordinate/time string did not match the model-admitted shape.
    /// `what` names the expected shape (`"HH:MM[:SS]"`, `"DDNMM"`,...).
    Parse {
        /// Expected textual format or component description.
        what: &'static str,
        /// Text that could not be parsed.
        input: String,
    },
    /// A parsed date/time component is out of range.
    OutOfRange {
        /// Name of the civil date/time component outside its accepted range.
        field: &'static str,
        /// Rejected integer component, in the units of the named field.
        value: i32,
    },
    /// A duration cannot mix civil timestamps with and without UTC offsets.
    MixedTimezoneAwareness,
}

impl fmt::Display for ArabicError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite { field, value } => {
                write!(f, "arabic_parts: non-finite {field}: {value}")
            }
            Self::InvalidLatitude { value } => {
                write!(f, "arabic_parts: latitude {value} outside [-90, 90]")
            }
            Self::InvalidLongitude { value } => {
                write!(f, "arabic_parts: longitude {value} outside [-180, 360]")
            }
            Self::HighLatitude3dUnavailable { lat, jd } => write!(
                f,
                "arabic_parts: 3-D sect at |lat| {lat} > 60 (jd {jd}) needs the azalt pipeline, which this lane does not provide"
            ),
            Self::Parse { what, input } => {
                write!(f, "arabic_parts: cannot parse {input:?} as {what}")
            }
            Self::MixedTimezoneAwareness => {
                f.write_str("cannot subtract timestamps with mixed timezone awareness")
            }
            Self::OutOfRange { field, value } => {
                write!(f, "arabic_parts: {field} = {value} out of range")
            }
        }
    }
}

impl std::error::Error for ArabicError {}

/// Crate-local result alias (this module must not touch `errors.rs`).
pub type ArabicResult<T> = std::result::Result<T, ArabicError>;

/// Reject non-finite inputs with the field-tagged failure the module doc advertises.
fn require_finite(field: &'static str, value: f64) -> ArabicResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(ArabicError::NonFinite { field, value })
    }
}

/// Normalize a finite angle in degrees to the half-open circle `[0, 360)`.
/// A negative subnormal remainder can round to 360; that circular endpoint
/// is the same direction as zero and is represented canonically as zero.
#[inline]
pub fn wrap360(x: f64) -> f64 {
    let value = x.rem_euclid(360.0);
    if value == 360.0 { 0.0 } else { value }
}

/// Zero-based index of the 108 equal zodiac quarters, from a normalized angle.
/// A quarter spans exactly 10/3 degrees. One native division/floor supplies
/// both lunar-mansion quarters and the equal navamsa partitions.
fn zodiac_quarter(normalized_longitude: f64) -> u8 {
    (normalized_longitude * 3.0 / 10.0).floor() as u8
}

/// Latitude threshold above which the model switches to the 3-D
/// horizontal-coordinate day/night arm. Value and switching role preserved.
pub const EXTREME_LATITUDE_THRESHOLD: f64 = 60.0;

// ---------------------------------------------------------------------------
// Arabic parts
// ---------------------------------------------------------------------------

/// Part of Fortune (Pars Fortunae): day `ASC + Moon − Sun`,
/// night `ASC + Sun − Moon`.
/// Inputs accept any finite ecliptic longitude; the result is an ecliptic longitude in degrees wrapped to [0, 360).
pub fn part_of_fortune(asc: f64, sun: f64, moon: f64, diurnal: bool) -> ArabicResult<f64> {
    let asc = require_finite("asc", asc)?;
    let sun = require_finite("sun", sun)?;
    let moon = require_finite("moon", moon)?;
    Ok(if diurnal {
        wrap360(asc + moon - sun)
    } else {
        wrap360(asc + sun - moon)
    })
}

/// Part of Spirit (Pars Spiritus / Daimon): the sect mirror of Fortune.
/// Inputs accept any finite ecliptic longitude; the result is an ecliptic longitude in degrees wrapped to [0, 360).
pub fn part_of_spirit(asc: f64, sun: f64, moon: f64, diurnal: bool) -> ArabicResult<f64> {
    let asc = require_finite("asc", asc)?;
    let sun = require_finite("sun", sun)?;
    let moon = require_finite("moon", moon)?;
    Ok(if diurnal {
        wrap360(asc + sun - moon)
    } else {
        wrap360(asc + moon - sun)
    })
}

/// Part of Love (Pars Amoris): `ASC + Venus − Sun` (sect-independent
/// project convention at the pin, not the §476 lot).
/// Inputs accept any finite ecliptic longitude; the result is an ecliptic longitude in degrees wrapped to [0, 360).
pub fn part_of_love(asc: f64, venus: f64, sun: f64) -> ArabicResult<f64> {
    let asc = require_finite("asc", asc)?;
    let venus = require_finite("venus", venus)?;
    let sun = require_finite("sun", sun)?;
    Ok(wrap360(asc + venus - sun))
}

/// Part of Faith (Pars Fidei): `ASC + Mercury − Moon` (project
/// convention at the pin).
/// Inputs accept any finite ecliptic longitude; the result is an ecliptic longitude in degrees wrapped to [0, 360).
pub fn part_of_faith(asc: f64, mercury: f64, moon: f64) -> ArabicResult<f64> {
    let asc = require_finite("asc", asc)?;
    let mercury = require_finite("mercury", mercury)?;
    let moon = require_finite("moon", moon)?;
    Ok(wrap360(asc + mercury - moon))
}

/// 2-D ecliptic sect rule: `sun_rel = (sun − asc) % 360`; night iff
/// `0 < sun_rel < 180`. Boundary exactness preserved: the Sun exactly
/// on the Ascendant or Descendant counts as day.
pub fn is_day_chart_2d(sun_lon: f64, asc: f64) -> ArabicResult<bool> {
    let sun_lon = require_finite("sun_lon", sun_lon)?;
    let asc = require_finite("asc", asc)?;
    let sun_rel = wrap360(sun_lon - asc);
    Ok(!(0.0 < sun_rel && sun_rel < 180.0))
}

/// Checked observer context for the sect decision.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GeoContext {
    /// Julian Day (UT), forwarded to the 3-D pipeline when required.
    pub jd: f64,
    /// Geographic latitude in degrees, `[-90, 90]`.
    pub lat: f64,
    /// Geographic longitude in degrees East, `[-180, 360]`.
    pub lon: f64,
}

impl GeoContext {
    /// Checked constructor: rejects non-finite and out-of-range components.
    pub fn new(jd: f64, lat: f64, lon: f64) -> ArabicResult<Self> {
        let jd = require_finite("jd", jd)?;
        let lat = require_finite("geo_lat", lat)?;
        let lon = require_finite("geo_lon", lon)?;
        if !(-90.0..=90.0).contains(&lat) {
            return Err(ArabicError::InvalidLatitude { value: lat });
        }
        if !(-180.0..=360.0).contains(&lon) {
            return Err(ArabicError::InvalidLongitude { value: lon });
        }
        Ok(Self { jd, lat, lon })
    }
}

/// Sect decision mirroring `is_day_chart` at the pin: the 2-D rule for
/// `|lat| <= 60°` (or no observer supplied); an explicit
/// [`ArabicError::HighLatitude3dUnavailable`] refusal when the 3-D arm
/// would be required, since that pipeline is outside this lane.
/// `sun_lat` defaults to `0.0` at the pin and is carried for signature
/// numerical comparison (it only feeds the 3-D arm).
pub fn is_day_chart(
    sun_lon: f64,
    asc: f64,
    sun_lat: f64,
    observer: Option<&GeoContext>,
) -> ArabicResult<bool> {
    require_finite("sun_lat", sun_lat)?;
    match observer {
        Some(ctx) if ctx.lat.abs() > EXTREME_LATITUDE_THRESHOLD => {
            Err(ArabicError::HighLatitude3dUnavailable {
                lat: ctx.lat,
                jd: ctx.jd,
            })
        }
        _ => is_day_chart_2d(sun_lon, asc),
    }
}

/// Typed equivalent of the model five-key position map
/// (`Asc`, `Sun`, `Moon`, `Mercury`, `Venus`; `Sun_lat` optional,
/// default `0.0`). A struct makes a missing key a compile-time error
/// instead of the model's runtime `KeyError`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArabicInputs {
    /// Ascendant ecliptic longitude in degrees.
    pub asc: f64,
    /// Solar ecliptic longitude in degrees, used for sect and lot formulas.
    pub sun: f64,
    /// Lunar ecliptic longitude in degrees.
    pub moon: f64,
    /// Mercury's ecliptic longitude in degrees, used by the Part of Faith.
    pub mercury: f64,
    /// Venus's ecliptic longitude in degrees, used by the Part of Love.
    pub venus: f64,
    /// Solar ecliptic latitude in degrees; [`Self::ecliptic`] sets this to zero.
    pub sun_lat: f64,
}

impl ArabicInputs {
    /// Checked constructor (`sun_lat` defaults to `0.0` via [`Self::ecliptic`]).
    pub fn new(
        asc: f64,
        sun: f64,
        moon: f64,
        mercury: f64,
        venus: f64,
        sun_lat: f64,
    ) -> ArabicResult<Self> {
        require_finite("asc", asc)?;
        require_finite("sun", sun)?;
        require_finite("moon", moon)?;
        require_finite("mercury", mercury)?;
        require_finite("venus", venus)?;
        require_finite("sun_lat", sun_lat)?;
        Ok(Self {
            asc,
            sun,
            moon,
            mercury,
            venus,
            sun_lat,
        })
    }

    /// Constructor with `sun_lat = 0.0`.
    pub fn ecliptic(asc: f64, sun: f64, moon: f64, mercury: f64, venus: f64) -> ArabicResult<Self> {
        Self::new(asc, sun, moon, mercury, venus, 0.0)
    }
}

/// All four standard lots, keyed by their dictionary names
/// (`Pars_Fortunae`, `Pars_Spiritus`, `Pars_Amoris`, `Pars_Fidei`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ArabicParts {
    /// Part of Fortune ecliptic longitude, degrees in `[0, 360)`.
    pub pars_fortunae: f64,
    /// Part of Spirit ecliptic longitude, degrees in `[0, 360)`.
    pub pars_spiritus: f64,
    /// Part of Love ecliptic longitude, degrees in `[0, 360)`.
    pub pars_amoris: f64,
    /// Part of Faith ecliptic longitude, degrees in `[0, 360)`.
    pub pars_fidei: f64,
}

/// Evaluate all four lots with one shared sect decision.
pub fn calc_all_parts(
    inputs: &ArabicInputs,
    observer: Option<&GeoContext>,
) -> ArabicResult<ArabicParts> {
    let diurnal = is_day_chart(inputs.sun, inputs.asc, inputs.sun_lat, observer)?;
    Ok(ArabicParts {
        pars_fortunae: part_of_fortune(inputs.asc, inputs.sun, inputs.moon, diurnal)?,
        pars_spiritus: part_of_spirit(inputs.asc, inputs.sun, inputs.moon, diurnal)?,
        pars_amoris: part_of_love(inputs.asc, inputs.venus, inputs.sun)?,
        pars_fidei: part_of_faith(inputs.asc, inputs.mercury, inputs.moon)?,
    })
}

// ---------------------------------------------------------------------------
// Contrib constants
// ---------------------------------------------------------------------------

/// Zodiac sign index, 0 = Aries through 11 = Pisces.
pub const ARIES: i64 = 0;
/// Zodiac sign index 1 = Taurus.
pub const TAURUS: i64 = 1;
/// Zodiac sign index 2 = Gemini.
pub const GEMINI: i64 = 2;
/// Zodiac sign index 3 = Cancer.
pub const CANCER: i64 = 3;
/// Zodiac sign index 4 = Leo.
pub const LEO: i64 = 4;
/// Zodiac sign index 5 = Virgo.
pub const VIRGO: i64 = 5;
/// Zodiac sign index 6 = Libra.
pub const LIBRA: i64 = 6;
/// Zodiac sign index 7 = Scorpio.
pub const SCORPIO: i64 = 7;
/// Zodiac sign index 8 = Sagittarius.
pub const SAGITTARIUS: i64 = 8;
/// Zodiac sign index 9 = Capricorn.
pub const CAPRICORN: i64 = 9;
/// Zodiac sign index 10 = Aquarius.
pub const AQUARIUS: i64 = 10;
/// Zodiac sign index 11 = Pisces.
pub const PISCES: i64 = 11;

/// Sanskrit rasi alias of the English sign index above.
pub const MESHA: i64 = ARIES;
/// Sanskrit rasi alias of `TAURUS` (sign index 1).
pub const VRISHABA: i64 = TAURUS;
/// Sanskrit rasi alias of `GEMINI` (sign index 2).
pub const MITHUNA: i64 = GEMINI;
/// Sanskrit rasi alias of `CANCER` (sign index 3).
pub const KATAKA: i64 = CANCER;
/// Sanskrit rasi alias of `LEO` (sign index 4).
pub const SIMHA: i64 = LEO;
/// Sanskrit rasi alias of `VIRGO` (sign index 5).
pub const KANYA: i64 = VIRGO;
/// Sanskrit rasi alias of `LIBRA` (sign index 6).
pub const THULA: i64 = LIBRA;
/// Sanskrit rasi alias of `SCORPIO` (sign index 7).
pub const VRISHIKA: i64 = SCORPIO;
/// Sanskrit rasi alias of `SAGITTARIUS` (sign index 8).
pub const DHANUS: i64 = SAGITTARIUS;
/// Sanskrit rasi alias of `CAPRICORN` (sign index 9).
pub const MAKARA: i64 = CAPRICORN;
/// Sanskrit rasi alias of `AQUARIUS` (sign index 10).
pub const KUMBHA: i64 = AQUARIUS;
/// Sanskrit rasi alias of `PISCES` (sign index 11).
pub const MEENA: i64 = PISCES;

/// Vedic planet IDs.
pub const SURYA: i32 = 0;
/// Alternate Sun name, id 0 (same body as `SURYA`).
pub const RAVI: i32 = 0;
/// Moon id 1.
pub const CHANDRA: i32 = 1;
/// Alternate Moon name, id 1 (same body as `CHANDRA`).
pub const SOMA: i32 = 1;
/// Mercury id 2.
pub const BUDHA: i32 = 2;
/// Alternate Mercury name, id 2 (same body as `BUDHA`).
pub const SOUMYA: i32 = 2;
/// Alternate Mercury name carrying id 10 (per reference contrib).
pub const THAMA: i32 = 10;
/// Venus id 3.
pub const SUKRA: i32 = 3;
/// Alternate Venus name, id 3 (same body as `SUKRA`).
pub const BHARGAVA: i32 = 3;
/// Mars id 4.
pub const KUJA: i32 = 4;
/// Alternate Mars name, id 4 (same body as `KUJA`).
pub const ANGARAKA: i32 = 4;
/// Jupiter id 5.
pub const GURU: i32 = 5;
/// Alternate Jupiter name, id 5 (same body as `GURU`).
pub const BRIHASPATI: i32 = 5;
/// Saturn id 6.
pub const SANI: i32 = 6;
/// Alternate Saturn name, id 6 (same body as `SANI`).
pub const MANDA: i32 = 6;
/// Ascending lunar node id 10 (shares the value with `THAMA`).
pub const RAHU: i32 = 10;
/// Descending lunar node id −10.
pub const KETU: i32 = -10;
/// Alternate descending-node name, id −10 (same body as `KETU`).
pub const SIKHI: i32 = -10;

/// Classical aspect angles in degrees (integer family).
pub const CONJUNCTION: f64 = 0.0;
/// Semisextile aspect angle, 30 degrees.
pub const SEMISEXTILE: f64 = 30.0;
/// Sextile aspect angle, 60 degrees.
pub const SEXTILE: f64 = 60.0;
/// Square aspect angle, 90 degrees.
pub const SQUARE: f64 = 90.0;
/// Trine aspect angle, 120 degrees.
pub const TRINE: f64 = 120.0;
/// Quincunx aspect angle, 150 degrees.
pub const QUINCUNX: f64 = 150.0;
/// Opposition aspect angle, 180 degrees.
pub const OPPOSITION: f64 = 180.0;
/// Semisquare aspect angle, 45 degrees.
pub const SEMISQUARE: f64 = 45.0;
/// Sesquisquare aspect angle, 135 degrees.
pub const SESQUISQUARE: f64 = 135.0;
/// Quintile aspect angle, 72 degrees.
pub const QUINTILE: f64 = 72.0;
/// Biquintile aspect angle, 144 degrees.
pub const BIQUINTILE: f64 = 144.0;
/// Semiquintile aspect angle, 36 degrees.
pub const SEMIQUINTILE: f64 = 36.0;
/// Septile family `360/7`, `720/7`, `1080/7`.
pub const SEPTILE: f64 = 360.0 / 7.0;
/// Biseptile aspect angle, degrees (`720/7` float division).
pub const BISEPTILE: f64 = 720.0 / 7.0;
/// Triseptile aspect angle, degrees (`1080/7` float division).
pub const TRISEPTILE: f64 = 1080.0 / 7.0;
/// Novile aspect angle, 40 degrees (`360/9` family).
pub const NOVILE: f64 = 40.0;
/// Seminovile aspect angle, 20 degrees (`360/9` family).
pub const SEMINOVILE: f64 = 20.0;
/// Binovile aspect angle, 80 degrees (`360/9` family).
pub const BINOVILE: f64 = 80.0;
/// Quatronovile aspect angle, 160 degrees (`360/9` family).
pub const QUATRONOVILE: f64 = 160.0;
/// Undecile aspect angle, degrees (`360/11` float division).
pub const UNDECILE: f64 = 360.0 / 11.0;
/// Biundecile aspect angle, degrees (`720/11` float division).
pub const BIUNDECILE: f64 = 720.0 / 11.0;
/// Triundecile aspect angle, degrees (`1080/11` float division).
pub const TRIUNDECILE: f64 = 1080.0 / 11.0;
/// Quadundecile aspect angle, degrees (`1440/11` float division).
pub const QUADUNDECILE: f64 = 1440.0 / 11.0;
/// Quinundecile aspect angle, degrees (`1800/11` float division).
pub const QUINUNDECILE: f64 = 1800.0 / 11.0;
/// 15-degree-family sextile variant, 15 degrees.
pub const SQUISEXTILE: f64 = 15.0;
/// 15-degree-family square variant, 22.5 degrees.
pub const SQUISQUARE: f64 = 22.5;

/// English sign names in index order; display helper for sign indices.
const RASI_NAMES: [&str; 12] = [
    "Aries",
    "Taurus",
    "Gemini",
    "Cancer",
    "Leo",
    "Virgo",
    "Libra",
    "Scorpio",
    "Sagittarius",
    "Capricorn",
    "Aquarius",
    "Pisces",
];

/// Twenty-seven lunar-mansion names in index order; display helper for nakshatra indices.
const NAKSHATRA_NAMES: [&str; 27] = [
    "Aswini",
    "Bharani",
    "Kritikha",
    "Rohini",
    "Mrigasira",
    "Aridra",
    "Punarvasu",
    "Pushyami",
    "Aslesha",
    "Makha",
    "Pubba",
    "Uttara",
    "Hasta",
    "Chitta",
    "Swathi",
    "Vishaka",
    "Anuradha",
    "Jyesta",
    "Moola",
    "Poorvashada",
    "Uttarashada",
    "Sravana",
    "Dhanishta",
    "Satabisha",
    "Poorvabhadra",
    "Uttarabhadra",
    "Revathi",
];

/// Nakshatra width `360/27` and derived pada / navamsa widths.
/// All widths are in degrees; names are public here while the model keeps them private helpers.
pub const NAKSHATRA_WIDTH: f64 = 360.0 / 27.0;
/// Pada (quarter-nakshatra) width in degrees (`NAKSHATRA_WIDTH / 4`).
pub const PADA_WIDTH: f64 = NAKSHATRA_WIDTH / 4.0;
/// Navamsa (ninth of a 30-degree sign) width in degrees (`30/9`).
pub const NAVAMSA_WIDTH: f64 = 30.0 / 9.0;

// ---------------------------------------------------------------------------
// Contrib: Vedic conversion helpers
// ---------------------------------------------------------------------------

/// Rasi (sign `0 = Aries .. 11 = Pisces`) for a longitude.
pub fn long2rasi(deg: f64) -> ArabicResult<u8> {
    let deg = require_finite("lon", deg)?;
    Ok((wrap360(deg) / 30.0).floor() as u8)
}

/// `(nakshatra_id 0..26, pada 0..3)` for a longitude.
pub fn long2nakshatra(deg: f64) -> ArabicResult<(u8, u8)> {
    let norm = wrap360(require_finite("lon", deg)?);
    let quarter = zodiac_quarter(norm);
    Ok((quarter / 4, quarter % 4))
}

/// Navamsa `0..11` for a longitude (fire-from-Aries / earth-from-
/// Capricorn / air-from-Libra / water-from-Cancer element rule).
pub fn long2navamsa(deg: f64) -> ArabicResult<u8> {
    let norm = wrap360(require_finite("lon", deg)?);
    // The element-based starts 0,9,6,3 are the same continuous sequence
    // of equal 10/3-degree subdivisions modulo the twelve signs.
    Ok(zodiac_quarter(norm) % 12)
}

/// Normalize a rasi number to `[0, 11]`.
pub fn rasinorm(rasi: i64) -> u8 {
    rasi.rem_euclid(12) as u8
}

/// Forward rasi distance `r1 → r2` in `[0, 11]`.
pub fn rasi_dif(r1: i64, r2: i64) -> u8 {
    rasinorm(r2 - r1)
}

/// Signed rasi difference `r1 → r2` in `(−6, 6]`.
pub fn rasi_dif2(r1: i64, r2: i64) -> i8 {
    let d = rasi_dif(r1, r2) as i8;
    if d > 6 { d - 12 } else { d }
}

/// Antiscion and contra-antiscion as six-component position/speed arrays.
/// Reflect the longitude about the Cancer/Capricorn axis (`180-lon`), then
/// take its antipode. Unspecified latitude, distance and speed slots are zero.
pub fn antiscion(lon: f64) -> ArabicResult<([f64; 6], [f64; 6])> {
    let norm = wrap360(require_finite("lon", lon)?);
    let ant = wrap360(180.0 - norm);
    let contra = wrap360(ant + 180.0);
    Ok((
        [ant, 0.0, 0.0, 0.0, 0.0, 0.0],
        [contra, 0.0, 0.0, 0.0, 0.0, 0.0],
    ))
}

/// Sexagesimal split `(degrees-within-sign, sign, minutes, seconds)`
/// with the model carry rule (rounded `60` propagates through
/// seconds → minutes → degrees → sign). Half-second ties round to the
/// nearest even second before carries are applied.
pub fn degsplit(deg: f64) -> ArabicResult<(u8, u8, u8, u8)> {
    let norm = wrap360(require_finite("lon", deg)?);
    let mut sign = (norm / 30.0).floor() as u8;
    let within = norm - f64::from(sign) * 30.0;
    let mut d = within as u8;
    let rem_min = (within - f64::from(d)) * 60.0;
    let mut m = rem_min as u8;
    let mut s = ((rem_min - f64::from(m)) * 60.0).round_ties_even() as i64;
    if s >= 60 {
        s -= 60;
        m += 1;
    }
    if m >= 60 {
        m -= 60;
        d += 1;
    }
    if d >= 30 {
        d -= 30;
        sign = (sign + 1) % 12;
    }
    #[allow(clippy::cast_sign_loss)]
    Ok((d, sign, m, s as u8))
}

/// English zodiac name for a sign number (normalized `mod 12`).
pub fn sign_name(sign: i64) -> &'static str {
    RASI_NAMES[sign.rem_euclid(12) as usize]
}

/// Sanskrit nakshatra name (normalized `mod 27`).
pub fn nakshatra_name(nak: i64) -> &'static str {
    NAKSHATRA_NAMES[nak.rem_euclid(27) as usize]
}

/// Human-readable house-system name for a letter code; unknown codes
/// yield `"Unknown"`.
pub fn house_system_name(code: char) -> &'static str {
    match code {
        'A' => "Equal",
        'E' => "Equal",
        'B' => "Alcabitius",
        'C' => "Campanus",
        'D' => "Equal MC",
        'F' => "Carter poli-equatorial",
        'G' => "Gauquelin sectors",
        'H' => "Horizontal/Azimuthal",
        'I' => "Sunshine",
        'i' => "Sunshine alternative",
        'K' => "Koch",
        'L' => "Pullen S-Delta",
        'M' => "Morinus",
        'N' => "Whole sign",
        'O' => "Porphyry",
        'P' => "Placidus",
        'Q' => "Pullen S-Ratio",
        'R' => "Regiomontanus",
        'S' => "Sripati",
        'T' => "Polich/Page (topocentric)",
        'U' => "Krusinski-Pisa-Goelzer",
        'V' => "Vehlow equal",
        'W' => "Whole sign",
        'X' => "Axial rotation/Meridian",
        'Y' => "APC",
        _ => "Unknown",
    }
}

// ---------------------------------------------------------------------------
// Contrib: dignities / strengths / relations
// ---------------------------------------------------------------------------

// Per the model: Aries Mars, Taurus Venus, Gemini Mercury, Cancer Moon,
// Leo Sun, Virgo Mercury, Libra Venus, Scorpio Mars, Sagittarius
// Jupiter, Capricorn Saturn, Aquarius Saturn, Pisces Jupiter.
/// Domicile lords indexed by rasi 0..11: Aries Mars, Taurus Venus, Gemini Mercury, Cancer Moon, Leo Sun, Virgo Mercury, Libra Venus, Scorpio Mars, Sagittarius Jupiter, Capricorn Saturn, Aquarius Saturn, Pisces Jupiter.
const DOMICILE: [i32; 12] = [
    KUJA, SUKRA, BUDHA, CHANDRA, SURYA, BUDHA, SUKRA, KUJA, GURU, SANI, SANI, GURU,
];

/// Vedic domicile lord of a rasi (normalized `mod 12`).
pub fn lord(rasi: i64) -> i32 {
    DOMICILE[rasinorm(rasi) as usize]
}

/// Exaltation table: planet → `(sign, degree-within-sign)`.
fn exaltation(planet: i32) -> Option<(f64, f64)> {
    match planet {
        SURYA => Some((0.0, 10.0)),
        CHANDRA => Some((1.0, 3.0)),
        KUJA => Some((9.0, 28.0)),
        BUDHA => Some((5.0, 15.0)),
        GURU => Some((3.0, 5.0)),
        SUKRA => Some((11.0, 27.0)),
        SANI => Some((6.0, 20.0)),
        _ => None,
    }
}

/// Ochchabala (exaltation strength) `0` at debilitation → `60` at
/// exaltation, linear over the `180°` between them. Nodes and unknown
/// IDs yield `0.0`.
pub fn ochchabala(planet: i32, longitude: f64) -> ArabicResult<f64> {
    let longitude = require_finite("lon", longitude)?;
    let Some((sign, deg)) = exaltation(planet) else {
        return Ok(0.0);
    };
    let ex_long = sign * 30.0 + deg;
    let mut diff = wrap360(longitude - ex_long);
    if diff > 180.0 {
        diff = 360.0 - diff;
    }
    Ok(60.0 * (1.0 - diff / 180.0))
}

/// Residential strength: linear `1.0 → 0.0` across the sign
pub fn residential_strength(deg: f64) -> ArabicResult<f64> {
    let within = wrap360(require_finite("lon", deg)?) % 30.0;
    Ok(1.0 - within / 30.0)
}

/// Natural-relation row for a planet ID: column order Sun..Saturn (IDs 0..6); 2 friend, 1 neutral, 0 enemy.
fn naisargika_row(planet: i32) -> Option<[u8; 7]> {
    // Column order Sun..Saturn (IDs 0..6); 2 friend / 1 neutral / 0 enemy.
    match planet {
        SURYA => Some([1, 2, 1, 0, 2, 2, 0]),
        CHANDRA => Some([2, 1, 2, 1, 1, 1, 1]),
        KUJA => Some([2, 2, 0, 1, 1, 2, 1]),
        BUDHA => Some([2, 0, 1, 2, 1, 1, 1]),
        GURU => Some([2, 2, 0, 0, 2, 1, 1]),
        SUKRA => Some([0, 0, 2, 1, 1, 1, 2]),
        SANI => Some([0, 0, 2, 2, 0, 1, 1]),
        _ => None,
    }
}

/// Natural (`naisargika`) relation `0` enemy / `1` neutral / `2` friend.
/// Unknown IDs yield neutral `1`.
///
/// Column order is planet ID `0..=6` (Sun, Moon, Mercury, Venus, Mars,
/// Jupiter, Saturn); spot checks against the pin: `BUDHA → SUKRA = 2`,
/// `GURU → BUDHA = 0`, `GURU → KUJA = 2`, `SANI → KUJA = 0`.
pub fn naisargika_relation(p1: i32, p2: i32) -> u8 {
    if !(0..=6).contains(&p2) {
        return 1;
    }
    match naisargika_row(p1) {
        Some(row) => row[p2 as usize],
        None => 1,
    }
}

/// Temporal (`tatkalika`) friendship from the house separation
/// `p1 → p2`: same house neutral `1`,
/// houses 2, 3, 4, 10, 11, 12 friends `2`, else enemies `0`.
pub fn tatkalika_relation(p1_house: i64, p2_house: i64) -> u8 {
    let d = (p2_house - p1_house).rem_euclid(12);
    if d == 0 {
        1
    } else if matches!(d, 1 | 2 | 3 | 9 | 10 | 11) {
        2
    } else {
        0
    }
}

// ---------------------------------------------------------------------------
// Contrib: aspect predicates
// ---------------------------------------------------------------------------

/// Shortest angular separation of two longitudes, `0..=180`.
pub fn angular_separation(a_lon: f64, b_lon: f64) -> ArabicResult<f64> {
    let a = require_finite("a_lon", a_lon)?;
    let b = require_finite("b_lon", b_lon)?;
    let d = (wrap360(a) - wrap360(b)).abs();
    Ok(if d <= 180.0 { d } else { 360.0 - d })
}

/// Aspect predicate: `(is_match, deviation)` where `is_match` holds iff
/// the separation is within `orb` of `aspect`. the model planet-ID
/// arguments are dropped — they never enter the computation.
pub fn match_aspect(a_lon: f64, b_lon: f64, aspect: f64, orb: f64) -> ArabicResult<(bool, f64)> {
    let aspect = require_finite("aspect", aspect)?;
    let orb = require_finite("orb", orb)?;
    let delta = (angular_separation(a_lon, b_lon)? - aspect).abs();
    Ok((delta <= orb, delta))
}

/// [`match_aspect`] plus the signed `b − a` arc in `(−180, 180]`.
pub fn match_aspect_signed(
    a_lon: f64,
    b_lon: f64,
    aspect: f64,
    orb: f64,
) -> ArabicResult<(bool, f64, f64)> {
    let a = require_finite("a_lon", a_lon)?;
    let b = require_finite("b_lon", b_lon)?;
    let (hit, delta) = match_aspect(a, b, aspect, orb)?;
    let mut signed = wrap360(b - a);
    if signed > 180.0 {
        signed -= 360.0;
    }
    Ok((hit, delta, signed))
}

/// [`match_aspect_signed`] plus the absolute separation.
pub fn match_aspect_full(
    a_lon: f64,
    b_lon: f64,
    aspect: f64,
    orb: f64,
) -> ArabicResult<(bool, f64, f64, f64)> {
    let (hit, delta, signed) = match_aspect_signed(a_lon, b_lon, aspect, orb)?;
    Ok((hit, delta, signed, angular_separation(a_lon, b_lon)?))
}

// ---------------------------------------------------------------------------
// Contrib: pure time / coordinate helpers
// ---------------------------------------------------------------------------

/// Duration between two Julian Days in fractional days.
pub fn jduration(jd_start: f64, jd_end: f64) -> ArabicResult<f64> {
    let a = require_finite("jd_start", jd_start)?;
    let b = require_finite("jd_end", jd_end)?;
    Ok(b - a)
}

/// Parse one time/date field as `i32`, mapping failure to the field-tagged `ArabicError::Parse`.
fn parse_int(field: &'static str, piece: &str) -> ArabicResult<i32> {
    piece.trim().parse::<i32>().map_err(|_| ArabicError::Parse {
        what: field,
        input: piece.to_owned(),
    })
}

/// Parse `HH[:MM[:SS]]` into `(hour, minute, second)`.
/// Missing minute/second components are zero; malformed integers are errors.
pub fn parse_hms(time_string: &str) -> ArabicResult<(i32, i32, i32)> {
    let parts: Vec<&str> = time_string.trim().split(':').take(3).collect();
    if parts.iter().any(|p| p.trim().is_empty()) {
        return Err(ArabicError::Parse {
            what: "HH:MM[:SS]",
            input: time_string.to_owned(),
        });
    }
    let h = parse_int("HH:MM[:SS]", parts[0])?;
    let m = if parts.len() > 1 {
        parse_int("HH:MM[:SS]", parts[1])?
    } else {
        0
    };
    let s = if parts.len() > 2 {
        parse_int("HH:MM[:SS]", parts[2])?
    } else {
        0
    };
    Ok((h, m, s))
}

/// Gregorian civil fields of an ISO timestamp, before UTC-offset adjustment.
/// Fractional seconds are retained by the internal parser for durations but
/// omitted from this integral civil-field record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DateTimeParts {
    /// Gregorian year, 1 through 9999.
    pub year: i32,
    /// Calendar month, 1 through 12.
    pub month: u8,
    /// Day of the selected month.
    pub day: u8,
    /// Civil hour, 0 through 23.
    pub hour: u8,
    /// Civil minute, 0 through 59.
    pub minute: u8,
    /// Civil second, 0 through 59.
    pub second: u8,
}

struct IsoMoment {
    local: chrono::NaiveDateTime,
    utc: chrono::NaiveDateTime,
    aware: bool,
}

/// Read ISO calendar dates/times with optional fractional seconds and offset.
/// Naive input remains a civil time, without an inferred system timezone.
fn parse_iso(date_string: &str) -> ArabicResult<IsoMoment> {
    use chrono::{DateTime, Datelike, NaiveDate, NaiveDateTime, Timelike};
    let mut input = date_string.replace(' ', "T");
    if input.as_bytes().get(10) == Some(&b'T')
        && input.len() >= 13
        && matches!(input.as_bytes().get(13), None | Some(b'+' | b'-' | b'Z'))
    {
        input.insert_str(13, ":00");
    }
    let invalid = || ArabicError::Parse {
        what: "ISO Gregorian timestamp in years 1..9999",
        input: date_string.to_owned(),
    };
    let aware = DateTime::parse_from_rfc3339(&input).ok().or_else(|| {
        [
            "%Y-%m-%dT%H:%M%:z",
            "%Y-%m-%dT%H:%M%z",
            "%Y-%m-%dT%H:%M:%S%.f%z",
        ]
        .iter()
        .find_map(|format| DateTime::parse_from_str(&input, format).ok())
    });
    let parsed = if let Some(date) = aware {
        IsoMoment {
            local: date.naive_local(),
            utc: date.naive_utc(),
            aware: true,
        }
    } else {
        let date = ["%Y-%m-%dT%H:%M:%S%.f", "%Y-%m-%dT%H:%M", "%Y-%m-%dT%H"]
            .iter()
            .find_map(|format| NaiveDateTime::parse_from_str(&input, format).ok())
            .or_else(|| {
                NaiveDate::parse_from_str(&input, "%Y-%m-%d")
                    .ok()
                    .and_then(|date| date.and_hms_opt(0, 0, 0))
            })
            .ok_or_else(invalid)?;
        IsoMoment {
            local: date,
            utc: date,
            aware: false,
        }
    };
    if !(1..=9999).contains(&parsed.local.year()) || parsed.local.nanosecond() >= 1_000_000_000 {
        return Err(invalid());
    }
    Ok(parsed)
}

/// Parse an ISO date or timestamp into its six integral civil fields.
/// Fractional seconds and UTC offsets are accepted; malformed dates, leap
/// seconds and years outside 1..9999 return a parse error.
pub fn parse_datetime(date_string: &str) -> ArabicResult<DateTimeParts> {
    use chrono::{Datelike, Timelike};
    let date = parse_iso(date_string)?.local;
    Ok(DateTimeParts {
        year: date.year(),
        month: date.month() as u8,
        day: date.day() as u8,
        hour: date.hour() as u8,
        minute: date.minute() as u8,
        second: date.second() as u8,
    })
}

/// Decimal Julian years between two ISO timestamps, `b - a`.
/// One year is 365.25 days. Explicit UTC offsets are normalized before
/// subtraction; fractional seconds remain in the native duration. Mixing
/// timestamps with and without offsets returns `MixedTimezoneAwareness`.
pub fn years_diff(date_a: &str, date_b: &str) -> ArabicResult<f64> {
    let a = parse_iso(date_a)?;
    let b = parse_iso(date_b)?;
    if a.aware != b.aware {
        return Err(ArabicError::MixedTimezoneAwareness);
    }
    let duration = b.utc.signed_duration_since(a.utc);
    let seconds = duration.num_seconds() as f64 + f64::from(duration.subsec_nanos()) / 1e9;
    Ok(seconds / (365.25 * 86400.0))
}

/// Parse `12N34` / `12S34` / `12E34` / `12W34` into signed decimal
/// degrees (N/E positive, S/W negative); a bare number passes through.
/// With a hemisphere letter, a missing degree or minute component is zero.
/// Non-finite hemisphere components are rejected.
pub fn geoc2d(coord: &str) -> ArabicResult<f64> {
    let s = coord.trim().to_uppercase();
    for letter in ['N', 'S', 'W', 'E'] {
        if let Some(idx) = s.find(letter) {
            let (deg_part, min_part) = s.split_at(idx);
            let min_part = &min_part[letter.len_utf8()..];
            let deg: f64 = if deg_part.trim().is_empty() {
                0.0
            } else {
                deg_part.trim().parse().map_err(|_| ArabicError::Parse {
                    what: "DD[H]MM",
                    input: coord.to_owned(),
                })?
            };
            let mins: f64 = if min_part.trim().is_empty() {
                0.0
            } else {
                min_part.trim().parse().map_err(|_| ArabicError::Parse {
                    what: "DD[H]MM",
                    input: coord.to_owned(),
                })?
            };
            require_finite("coord", deg)?;
            require_finite("coord", mins)?;
            let sign = if letter == 'S' || letter == 'W' {
                -1.0
            } else {
                1.0
            };
            return Ok(sign * (deg + mins / 60.0));
        }
    }
    s.parse::<f64>().map_err(|_| ArabicError::Parse {
        what: "DD[H]MM",
        input: coord.to_owned(),
    })
}

/// Shared `DD[H]MM` formatter: hemisphere letter by sign, minutes rounded with `60` carried into degrees.
fn degmin_string(value: f64, pos: char, neg: char) -> ArabicResult<String> {
    let value = require_finite("coord", value)?;
    let hem = if value >= 0.0 { pos } else { neg };
    let abs = value.abs();
    let mut d = abs as i64;
    let mut m = ((abs - d as f64) * 60.0).round_ties_even() as i64;
    if m >= 60 {
        m -= 60;
        d += 1;
    }
    Ok(format!("{d}{hem}{m:02}"))
}

/// Decimal latitude → `DDNMM` / `DDSMM` (rounded `60` carried, as at
/// the pin).
pub fn geolat2c(lat: f64) -> ArabicResult<String> {
    degmin_string(lat, 'N', 'S')
}

/// Decimal longitude → `DDDEMM` / `DDDWMM` (rounded `60` carried).
pub fn geolon2c(lon: f64) -> ArabicResult<String> {
    degmin_string(lon, 'E', 'W')
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= EPS
    }

    // -- Fortune / Spirit sect rules -------------------------------------

    #[test]
    fn fortune_day_rule() {
        // ASC 15.5, Moon 240, Sun 120 → 15.5 + 240 − 120 = 135.5.
        assert!(close(
            part_of_fortune(15.5, 120.0, 240.0, true).unwrap(),
            135.5
        ));
    }

    #[test]
    fn fortune_night_reverses() {
        assert!(close(
            part_of_fortune(15.5, 120.0, 240.0, false).unwrap(),
            255.5
        ));
    }

    #[test]
    fn spirit_is_sect_mirror() {
        assert!(close(
            part_of_spirit(15.5, 120.0, 240.0, true).unwrap(),
            255.5
        ));
        assert!(close(
            part_of_spirit(15.5, 120.0, 240.0, false).unwrap(),
            135.5
        ));
    }

    #[test]
    fn fortune_wraps_negative_arc() {
        // 10 + 5 − 20 = −5 → 355.
        assert!(close(
            part_of_fortune(10.0, 20.0, 5.0, true).unwrap(),
            355.0
        ));
    }

    #[test]
    fn fortune_wraps_full_circle() {
        // 350 + 300 − 290 = 360 → 0.
        assert!(close(
            part_of_fortune(350.0, 290.0, 300.0, true).unwrap(),
            0.0
        ));
    }

    #[test]
    fn love_and_faith_formulas() {
        // ASC + Venus − Sun; ASC + Mercury − Moon.
        assert!(close(part_of_love(15.5, 110.0, 120.0).unwrap(), 5.5));
        assert!(close(part_of_faith(15.5, 130.0, 240.0).unwrap(), 265.5));
    }

    #[test]
    fn parts_reject_non_finite() {
        match part_of_fortune(f64::NAN, 0.0, 0.0, true) {
            Err(ArabicError::NonFinite { field, .. }) => assert_eq!(field, "asc"),
            other => panic!("expected NonFinite, got {:?}", other),
        }
        assert!(part_of_love(0.0, f64::INFINITY, 0.0).is_err());
        assert!(part_of_faith(0.0, 0.0, f64::NEG_INFINITY).is_err());
        assert!(part_of_spirit(0.0, f64::NAN, 0.0, false).is_err());
    }

    // -- Sect decision ----------------------------------------------------

    #[test]
    fn day_night_2d_nominal() {
        // Sun 120 past ASC 15.5: rel 104.5 < 180 → night.
        assert!(!is_day_chart_2d(120.0, 15.5).unwrap());
        // Sun 300: rel 284.5 → day.
        assert!(is_day_chart_2d(300.0, 15.5).unwrap());
    }

    #[test]
    fn day_night_boundaries_count_as_day() {
        // Sunrise (rel 0) and sunset (rel 180) are day charts.
        assert!(is_day_chart_2d(15.5, 15.5).unwrap());
        assert!(is_day_chart_2d(195.5, 15.5).unwrap());
        // Just inside the lower hemisphere → night.
        assert!(!is_day_chart_2d(15.5 + 1e-6, 15.5).unwrap());
        assert!(!is_day_chart_2d(195.5 - 1e-6, 15.5).unwrap());
        // Just outside → day.
        assert!(is_day_chart_2d(15.5 - 1e-6, 15.5).unwrap());
        assert!(is_day_chart_2d(195.5 + 1e-6, 15.5).unwrap());
    }

    #[test]
    fn sect_with_moderate_observer_uses_2d() {
        let ctx = GeoContext::new(2451545.0, 45.0, 10.0).unwrap();
        assert_eq!(
            is_day_chart(120.0, 15.5, 0.0, Some(&ctx)),
            is_day_chart_2d(120.0, 15.5)
        );
        assert_eq!(
            is_day_chart(120.0, 15.5, 0.0, None),
            is_day_chart_2d(120.0, 15.5)
        );
    }

    #[test]
    fn sect_high_latitude_refuses() {
        let ctx = GeoContext::new(2451545.0, 70.0, 10.0).unwrap();
        assert_eq!(
            is_day_chart(120.0, 15.5, 0.0, Some(&ctx)),
            Err(ArabicError::HighLatitude3dUnavailable {
                lat: 70.0,
                jd: 2451545.0
            })
        );
        // Exactly 60° still takes the 2-D arm.
        let edge = GeoContext::new(2451545.0, 60.0, 10.0).unwrap();
        assert!(is_day_chart(300.0, 15.5, 0.0, Some(&edge)).is_ok());
    }

    #[test]
    fn geo_context_validation() {
        assert!(GeoContext::new(2451545.0, 91.0, 0.0).is_err());
        assert!(GeoContext::new(2451545.0, 0.0, 361.0).is_err());
        assert!(GeoContext::new(f64::NAN, 0.0, 0.0).is_err());
        assert!(GeoContext::new(2451545.0, -90.0, -180.0).is_ok());
    }

    #[test]
    fn calc_all_parts_doc_example() {
        let inputs = ArabicInputs::ecliptic(15.5, 120.0, 240.0, 130.0, 110.0).unwrap();
        let parts = calc_all_parts(&inputs, None).unwrap();
        // Sun rel 104.5 → night: Fortune = ASC + Sun − Moon = 255.5.
        assert!(close(parts.pars_fortunae, 255.5));
        assert!(close(parts.pars_spiritus, 135.5));
        assert!(close(parts.pars_amoris, 5.5));
        assert!(close(parts.pars_fidei, 265.5));
    }

    // -- Vedic partitions -------------------------------------------------

    #[test]
    fn rasi_boundaries() {
        assert_eq!(long2rasi(0.0).unwrap(), 0);
        assert_eq!(long2rasi(29.999).unwrap(), 0);
        assert_eq!(long2rasi(30.0).unwrap(), 1);
        assert_eq!(long2rasi(359.9).unwrap(), 11);
        assert_eq!(long2rasi(360.0).unwrap(), 0);
        assert_eq!(long2rasi(-0.1).unwrap(), 11);
        assert!(long2rasi(f64::NAN).is_err());
    }

    #[test]
    fn nakshatra_pada() {
        assert_eq!(long2nakshatra(0.0).unwrap(), (0, 0));
        assert_eq!(long2nakshatra(360.0 - 1e-9).unwrap(), (26, 3));
        // 13°20′ = one nakshatra width → (1, 0).
        let (nak, pada) = long2nakshatra(40.0 / 3.0).unwrap();
        assert_eq!((nak, pada), (1, 0));
        // Mid first nakshatra: 6.6667° → (0, 2).
        assert_eq!(long2nakshatra(20.0 / 3.0).unwrap(), (0, 2));
    }

    #[test]
    fn equal_zodiac_partitions_share_the_same_quarter_index() {
        for degree in 0..360 {
            // Integer-degree reference: floor(3*longitude/10), evaluated
            // exactly with small integers rather than rounded width division.
            let quarter = ((3 * degree) / 10) as u8;
            assert_eq!(
                long2nakshatra(f64::from(degree)).unwrap(),
                (quarter / 4, quarter % 4)
            );
            assert_eq!(long2navamsa(f64::from(degree)).unwrap(), quarter % 12);
            assert_eq!(
                long2nakshatra(f64::from(degree) - 720.0).unwrap(),
                (quarter / 4, quarter % 4)
            );
        }
        assert_eq!(long2nakshatra(120.0).unwrap(), (9, 0));
        assert_eq!(long2nakshatra(240.0).unwrap(), (18, 0));
        assert!(long2nakshatra(f64::INFINITY).is_err());
        assert!(long2navamsa(f64::NAN).is_err());
    }

    #[test]
    fn wrapped_subnormal_angles_keep_a_valid_half_open_partition() {
        for angle in [-f64::from_bits(1), -0.0, 0.0, 360.0, f64::MAX] {
            let norm = wrap360(angle);
            assert!((0.0..360.0).contains(&norm));
            assert!(long2rasi(angle).unwrap() < 12);
            let (mansion, quarter) = long2nakshatra(angle).unwrap();
            assert!(mansion < 27 && quarter < 4);
            assert!(long2navamsa(angle).unwrap() < 12);
        }
    }

    #[test]
    fn navamsa_element_starts() {
        assert_eq!(long2navamsa(0.0).unwrap(), 0); // Aries from Aries
        assert_eq!(long2navamsa(30.0).unwrap(), 9); // Taurus from Capricorn
        assert_eq!(long2navamsa(60.0).unwrap(), 6); // Gemini from Libra
        assert_eq!(long2navamsa(90.0).unwrap(), 3); // Cancer from Cancer
        // 10 degrees is exactly the fourth equal 10/3-degree subdivision.
        assert_eq!(long2navamsa(10.0).unwrap(), 3);
    }

    #[test]
    fn rasi_distances() {
        assert_eq!(rasinorm(12), 0);
        assert_eq!(rasinorm(-1), 11);
        assert_eq!(rasi_dif(0, 3), 3);
        assert_eq!(rasi_dif(10, 2), 4);
        assert_eq!(rasi_dif2(0, 3), 3);
        assert_eq!(rasi_dif2(0, 9), -3);
        assert_eq!(rasi_dif2(0, 6), 6);
        assert_eq!(rasi_dif2(0, 7), -5);
    }

    #[test]
    fn antiscion_values() {
        let (ant, contra) = antiscion(0.0).unwrap();
        assert!(close(ant[0], 180.0));
        assert!(close(contra[0], 0.0));
        assert_eq!(ant[3], 0.0);
        let (ant, contra) = antiscion(10.0).unwrap();
        assert!(close(ant[0], 170.0));
        assert!(close(contra[0], 350.0));
        assert_eq!(ant[3], 0.0);
        assert_eq!(antiscion(ant[0]).unwrap().0[0], 10.0);
    }

    #[test]
    fn degsplit_nominal_and_carry() {
        assert_eq!(degsplit(15.5).unwrap(), (15, 0, 30, 0));
        assert_eq!(degsplit(30.0).unwrap(), (0, 1, 0, 0));
        // 359.9999999999 rounds seconds to 60 → carries into Pisces→Aries.
        assert_eq!(degsplit(359.9999999999).unwrap(), (0, 0, 0, 0));
    }

    #[test]
    fn degsplit_half_seconds_use_even_rounding() {
        // These binary64 degrees produce exact half-second residuals after
        // the two successive sexagesimal multiplications.
        assert_eq!(degsplit(0.5 / 3600.0).unwrap(), (0, 0, 0, 0));
        assert_eq!(degsplit(1.5 / 3600.0).unwrap(), (0, 0, 0, 2));
        assert_eq!(degsplit(4.5 / 3600.0).unwrap(), (0, 0, 0, 4));
        assert_eq!(degsplit(3.5 / 3600.0).unwrap(), (0, 0, 0, 4));
    }

    #[test]
    fn names_and_house_systems() {
        assert_eq!(sign_name(0), "Aries");
        assert_eq!(sign_name(12), "Aries");
        assert_eq!(sign_name(-1), "Pisces");
        assert_eq!(nakshatra_name(0), "Aswini");
        assert_eq!(nakshatra_name(27), "Aswini");
        assert_eq!(nakshatra_name(26), "Revathi");
        assert_eq!(house_system_name('P'), "Placidus");
        assert_eq!(house_system_name('S'), "Sripati");
        assert_eq!(house_system_name('Z'), "Unknown");
    }

    #[test]
    fn lords() {
        assert_eq!(lord(0), KUJA);
        assert_eq!(lord(1), SUKRA);
        assert_eq!(lord(4), SURYA);
        assert_eq!(lord(9), SANI);
        assert_eq!(lord(12), KUJA);
    }

    #[test]
    fn ochchabala_endpoints_and_unknown() {
        // Sun exalted at 10° Aries → 60; debilitated 10° Libra → 0.
        assert!(close(ochchabala(SURYA, 10.0).unwrap(), 60.0));
        assert!(close(ochchabala(SURYA, 190.0).unwrap(), 0.0));
        assert!(close(ochchabala(SURYA, 100.0).unwrap(), 30.0));
        assert_eq!(ochchabala(RAHU, 10.0).unwrap(), 0.0);
        assert_eq!(ochchabala(99, 10.0).unwrap(), 0.0);
    }

    #[test]
    fn residential_strength_endpoints() {
        assert!(close(residential_strength(0.0).unwrap(), 1.0));
        assert!(close(residential_strength(15.0).unwrap(), 0.5));
        assert!(residential_strength(29.999).unwrap() < 0.001);
    }

    #[test]
    fn relations() {
        assert_eq!(naisargika_relation(SURYA, SUKRA), 0);
        assert_eq!(naisargika_relation(SURYA, KUJA), 2);
        assert_eq!(naisargika_relation(BUDHA, SUKRA), 2);
        assert_eq!(naisargika_relation(GURU, BUDHA), 0);
        assert_eq!(naisargika_relation(GURU, KUJA), 2);
        assert_eq!(naisargika_relation(SANI, KUJA), 0);
        assert_eq!(naisargika_relation(CHANDRA, SOMA), 1);
        assert_eq!(naisargika_relation(99, SURYA), 1);
        assert_eq!(naisargika_relation(SURYA, 99), 1);
        assert_eq!(tatkalika_relation(0, 0), 1);
        assert_eq!(tatkalika_relation(0, 1), 2);
        assert_eq!(tatkalika_relation(0, 3), 2);
        assert_eq!(tatkalika_relation(0, 9), 2);
        assert_eq!(tatkalika_relation(0, 4), 0);
        assert_eq!(tatkalika_relation(0, 6), 0);
    }

    // -- Aspects ----------------------------------------------------------

    #[test]
    fn separation_and_match() {
        assert!(close(angular_separation(10.0, 190.0).unwrap(), 180.0));
        assert!(close(angular_separation(10.0, 350.0).unwrap(), 20.0));
        let (hit, delta) = match_aspect(0.0, 90.5, SQUARE, 1.0).unwrap();
        assert!(hit);
        assert!(close(delta, 0.5));
        let (hit, _) = match_aspect(0.0, 92.5, SQUARE, 1.0).unwrap();
        assert!(!hit);
        let (hit, _, signed, sep) = match_aspect_full(10.0, 350.0, CONJUNCTION, 25.0).unwrap();
        assert!(hit);
        assert!(close(signed, -20.0));
        assert!(close(sep, 20.0));
    }

    // -- Time / coordinate helpers ----------------------------------------

    #[test]
    fn time_helpers() {
        assert!(close(jduration(2451545.0, 2451546.5).unwrap(), 1.5));
        assert_eq!(parse_hms("12:34:56").unwrap(), (12, 34, 56));
        assert_eq!(parse_hms("3:04").unwrap(), (3, 4, 0));
        assert_eq!(parse_hms("12").unwrap(), (12, 0, 0));
        assert!(parse_hms("ab:cd").is_err());
        assert!(close(
            years_diff("2000-01-01", "2001-01-01").unwrap(),
            366.0 / 365.25
        ));
        assert!(close(
            years_diff("2001-01-01", "2000-01-01").unwrap(),
            -366.0 / 365.25
        ));
        assert!(parse_datetime("2020-02-30").is_err());
        assert!(parse_datetime("2020-13-01").is_err());
        assert_eq!(parse_datetime("2024-03-20T12").unwrap().hour, 12);
        assert_eq!(parse_datetime("2024-03-20T12:34+02:00").unwrap().minute, 34);
        let civil = parse_datetime("2024-03-20T12:34:56.123456+02:00").unwrap();
        assert_eq!((civil.hour, civil.minute, civil.second), (12, 34, 56));
        assert_eq!(
            years_diff("2024-01-01T12:00:00+02:00", "2024-01-01T10:00:00Z").unwrap(),
            0.0
        );
        assert!(
            (years_diff("2024-01-01T00:00:00.1", "2024-01-01T00:00:00.2").unwrap()
                * 365.25
                * 86400.0
                - 0.1)
                .abs()
                < 1e-15
        );
        assert_eq!(
            years_diff("2024-01-01", "2024-01-01T00:00:00Z"),
            Err(ArabicError::MixedTimezoneAwareness)
        );
    }

    #[test]
    fn coord_helpers() {
        assert!(close(geoc2d("12N34").unwrap(), 12.0 + 34.0 / 60.0));
        assert!(close(geoc2d("12S34").unwrap(), -(12.0 + 34.0 / 60.0)));
        assert!(close(geoc2d("10W30").unwrap(), -(10.5)));
        assert!(close(geoc2d("10.5").unwrap(), 10.5));
        assert!(geoc2d("abc").is_err());
        assert_eq!(geoc2d("N30").unwrap(), 0.5);
        assert_eq!(geoc2d("W").unwrap(), -0.0);
        assert_eq!(geolat2c(12.5).unwrap(), "12N30");
        assert_eq!(geolat2c(-12.5).unwrap(), "12S30");
        assert_eq!(geolon2c(10.5).unwrap(), "10E30");
        assert_eq!(geolon2c(-10.5).unwrap(), "10W30");
        // Rounded 60′ carries into degrees.
        assert_eq!(geolat2c(12.999999999).unwrap(), "13N00");
    }

    #[test]
    fn error_display_is_stable() {
        let e = ArabicError::HighLatitude3dUnavailable {
            lat: 70.0,
            jd: 2451545.0,
        };
        assert!(format!("{e}").contains("azalt"));
    }
}
