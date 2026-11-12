//! Thin bindings for native convenience calculations.
use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::arabic_parts as native;
use libephemeris_core::errors::Result;
use libephemeris_core::python_support::{Session, contrib};
use pyo3::exceptions::{PyException, PyNotImplementedError, PyTypeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

pyo3::create_exception!(
    _contrib,
    Error,
    PyException,
    "Compatibility exception for libephemeris.contrib, derived from Exception.\n\nNumerical helpers can instead raise ValueError, TypeError or native engine\nexceptions; this class is not a catch-all base for every convenience failure."
);
/// Empty extensible attribute container for application convenience data.
///
/// Data() takes no arguments and creates an independent __dict__. No chart
/// fields, units, defaults, validation or ephemeris resources are installed;
/// assign application-owned attributes explicitly.
#[pyclass(module = "libephemeris.contrib", dict)]
struct Data;
#[pymethods]
impl Data {
    /// Create an empty Data attribute container; takes no arguments.
    #[new]
    fn new() -> Self {
        Self
    }
}
/// Empty extensible attribute container for application-owned user information.
///
/// User() takes no arguments and creates an independent __dict__. It performs
/// no database access, authentication, chart calculation or field validation.
/// There are no predefined attributes or implicit observing settings.
#[pyclass(module = "libephemeris.contrib", dict)]
struct User;
#[pymethods]
impl User {
    /// Create an empty User attribute container; takes no arguments.
    #[new]
    fn new() -> Self {
        Self
    }
}

fn pure<T: Send>(
    py: Python<'_>,
    operation: impl FnOnce() -> native::ArabicResult<T> + Send,
) -> PyResult<T> {
    py.detach(operation)
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
fn call<T: Send>(
    py: Python<'_>,
    operation: impl FnOnce(&Session) -> Result<T> + Send,
) -> PyResult<T> {
    py.detach(|| operation(runtime().shared_session()?.as_ref()))
        .map_err(|error| to_pyerr(py, error))
}
macro_rules! longitude {
    ($name:ident,$result:ty,$doc:literal) => {
        #[doc = $doc]
        #[pyfunction]
        fn $name(py: Python<'_>, deg: f64) -> PyResult<$result> {
            pure(py, || native::$name(deg))
        }
    };
}
longitude!(
    long2rasi,
    u8,
    "Return the zodiac sign index for finite longitude deg in degrees.\n\nNormalizes modulo 360; output is 0 (Aries) through 11 (Pisces). Nonfinite\ninput raises ValueError. No clock, chart state or ephemeris is consulted."
);
longitude!(
    long2nakshatra,
    (u8, u8),
    "Return (nakshatra_index, pada_index) for finite longitude deg in degrees.\n\nLongitude is normalized modulo 360; mansion indices are 0..26 and quarters\n0..3, not one-based. Nonfinite input raises ValueError. No sidereal correction\nis inferred: supply longitude in the zodiac convention you intend."
);
longitude!(
    long2navamsa,
    u8,
    "Return the navamsa sign index 0..11 for finite longitude deg in degrees.\n\nUses equal 10/3-degree partitions modulo the zodiac, with element-based\nstarting signs. Nonfinite input raises ValueError. No ayanamsha or clock\nis acquired; the caller chooses the longitude's zodiac convention."
);
longitude!(
    degsplit,
    (u8, u8, u8, u8),
    "Split finite longitude deg into (degrees_within_sign, sign, minutes, seconds).\n\nInput is degrees, normalized modulo 360. sign is 0..11; seconds round to\nthe nearest integer with ties to even and carries through minutes/degrees/\nsign. Nonfinite input raises ValueError. This is not the split_deg layout."
);
longitude!(
    residential_strength,
    f64,
    "Return dimensionless residential strength for finite longitude deg in degrees.\n\nAfter modulo-360 normalization, strength decreases linearly from 1 towards 0\nacross each 30-degree sign, resetting to 1 at its boundary. Nonfinite input\nraises ValueError. This convention does not acquire or choose an ayanamsha."
);
/// Normalize integer sign index rasi modulo 12, returning an int in 0..11.
///
/// Negative indices wrap normally. Input must fit a signed 64-bit integer;
/// extraction overflow raises OverflowError. No chart, units or data are read.
#[pyfunction]
fn rasinorm(py: Python<'_>, rasi: i64) -> u8 {
    py.detach(|| native::rasinorm(rasi))
}
/// Return forward sign distance r1 to r2 modulo 12, as an int in 0..11.
///
/// Inputs are integer sign indices, not longitude degrees. Use ordinary sign
/// indices: native r2-r1 must fit signed 64 bits; overflowing differences are
/// outside this helper's domain and debug builds can panic. No chart is read.
#[pyfunction]
fn rasi_dif(py: Python<'_>, r1: i64, r2: i64) -> u8 {
    py.detach(|| native::rasi_dif(r1, r2))
}
/// Return signed sign distance r1 to r2 in (-6, 6], with a half-turn = +6.
///
/// Inputs are integer sign indices, not degrees. As with rasi_dif, r2-r1 must
/// fit signed 64 bits; extreme overflowing differences are outside the domain
/// and debug builds can panic. No chart, clock or ephemeris is consulted.
#[pyfunction]
fn rasi_dif2(py: Python<'_>, r1: i64, r2: i64) -> i8 {
    py.detach(|| native::rasi_dif2(r1, r2))
}
/// Return (antiscion, contra_antiscion) for finite longitude lon in degrees.
///
/// Each result is a six-float tuple: longitude followed by five zeros. The
/// longitudes are (180-lon) modulo 360 and its antipode. Latitude, distance
/// and rates are placeholders, not physical states. Nonfinite lon raises
/// ValueError. No chart frame, clock or ephemeris is acquired.
#[pyfunction]
fn antiscion<'py>(
    py: Python<'py>,
    lon: f64,
) -> PyResult<(Bound<'py, PyTuple>, Bound<'py, PyTuple>)> {
    let (a, b) = pure(py, || native::antiscion(lon))?;
    Ok((PyTuple::new(py, a)?, PyTuple::new(py, b)?))
}
/// Return the English zodiac name for integer sign, normalized modulo 12.
///
/// 0 is Aries and 11 Pisces. Negative indices wrap; extraction requires a
/// signed 64-bit integer. No chart or ephemeris is consulted.
#[pyfunction]
fn signtostr(py: Python<'_>, sign: i64) -> &'static str {
    py.detach(|| native::sign_name(sign))
}
/// Return the catalog mansion name for integer nak, normalized modulo 27.
///
/// Indices are zero-based (0 = Aswini); negative values wrap. Input must fit
/// signed 64 bits. This naming lookup does not calculate a stellar position.
#[pyfunction]
fn get_nakshatra_name(py: Python<'_>, nak: i64) -> &'static str {
    py.detach(|| native::nakshatra_name(nak))
}
/// Return the convenience-module display name for the complete string code.
///
/// code must match a one-character entry exactly; there is no ASCII case
/// folding here, and 'I'/'i' have distinct names. Unknown, empty or multi-
/// character strings return 'Unknown', not the main house_name empty string.
/// This lookup does not calculate a house or check its polar domain.
#[pyfunction]
fn house_system_name(py: Python<'_>, code: &str) -> &'static str {
    py.detach(|| contrib::house_name(code))
}
/// Return the classical domicile ruler's contrib planet ID for sign rasi.
///
/// Integer rasi is normalized modulo 12 (0 = Aries). The result uses this
/// module's Vedic IDs, not a newly acquired planetary position. Extraction
/// requires a signed 64-bit input; no clock/data or additional validation.
#[pyfunction]
fn lord(py: Python<'_>, rasi: i64) -> i32 {
    py.detach(|| native::lord(rasi))
}
/// Return conventional exaltation strength, 0 at debilitation to 60 at exaltation.
///
/// planet is a contrib classical planet ID; longitude is finite ecliptic
/// degrees in the intended zodiac. Strength is linear over 180 degrees.
/// Nodes/unknown IDs return 0, but nonfinite longitude still raises ValueError.
/// No ayanamsha, observer or planetary state is acquired.
#[pyfunction]
fn ochchabala(py: Python<'_>, planet: i32, longitude: f64) -> PyResult<f64> {
    pure(py, || native::ochchabala(planet, longitude))
}
/// Return the natural-relation code from planet p1 to planet p2.
///
/// Uses contrib classical IDs 0..6: 0 enemy, 1 neutral, 2 friend. Unknown
/// IDs return neutral (1). This asymmetric table is a convention, not a
/// calculation from orbital positions; no clock or resources are consulted.
#[pyfunction]
fn naisargika_relation(py: Python<'_>, p1: i32, p2: i32) -> u8 {
    py.detach(|| native::naisargika_relation(p1, p2))
}
/// Return temporal friendship from integer house indices p1_house to p2_house.
///
/// Same house is neutral (1); forward differences 1,2,3,9,10,11 modulo 12
/// are friendly (2), others enemy (0). Consistent zero- or one-based indices
/// give the same result. The signed 64-bit subtraction must be representable;
/// extreme overflowing inputs are outside the domain and can panic in debug
/// builds. No positions, clock or ephemeris are acquired.
#[pyfunction]
fn tatkalika_relation(py: Python<'_>, p1_house: i64, p2_house: i64) -> u8 {
    py.detach(|| native::tatkalika_relation(p1_house, p2_house))
}
/// Return twelve Sripati/Raman cusp longitudes in degrees at UT1 JD jd.
///
/// geolon is east longitude [-180,360] and geolat north latitude [-90,90],
/// both finite degrees. Returns an unpadded tuple: index 0 is house 1, not
/// a leading-zero list or a (cusps, angles) pair. Uses configured house clocks;
/// coordinate, clock and construction failures propagate. No alternate
/// system is automatically selected and caller observer state is unchanged.
#[pyfunction]
fn raman_houses<'py>(
    py: Python<'py>,
    jd: f64,
    geolon: f64,
    geolat: f64,
) -> PyResult<Bound<'py, PyTuple>> {
    PyTuple::new(
        py,
        call(py, |session| {
            session.contrib_raman_houses(jd, geolon, geolat)
        })?,
    )
}
/// Return apparent longitudes of Aldebaran, Regulus, Antares and Fomalhaut.
///
/// jd is UT1 Julian day; the four-float tuple is ecliptic degrees in that
/// order, calculated through the configured fixed-star route. An unavailable
/// individual place becomes NaN with a native debug log, not None or a
/// shortened tuple. Runtime/session acquisition errors can still propagate
/// before the four lookups. The name does not calculate Saturn's position.
#[pyfunction]
fn saturn_4_stars(py: Python<'_>, jd: f64) -> PyResult<(f64, f64, f64, f64)> {
    call(py, |session| Ok(session.contrib_saturn_four_stars(jd))).map(|[a, b, c, d]| (a, b, c, d))
}
/// Return the operating system's CURRENT UTC civil Julian day as a float.
///
/// Uses Unix-epoch wall-clock seconds; it is not TT or a UT1 measurement and
/// applies no leap-second/Earth-orientation correction. Takes no arguments,
/// reads no ephemeris and does not provide a monotonic duration clock.
#[pyfunction]
fn jdnow(py: Python<'_>) -> f64 {
    py.detach(contrib::current_utc_jd)
}
/// Return jd_end minus jd_start in fractional days.
///
/// Both inputs must be finite Julian days in the SAME caller-chosen scale;
/// no UTC/TT/UT1 conversion is performed. Nonfinite inputs raise ValueError.
/// Negative results are allowed; there is no system-clock or ephemeris read.
#[pyfunction]
fn jduration(py: Python<'_>, jd_start: f64, jd_end: f64) -> PyResult<f64> {
    pure(py, || native::jduration(jd_start, jd_end))
}
/// Parse `HH[:MM[:SS]]` into an integer (hour, minute, second) tuple.
///
/// time_string may omit minute/second fields, which become zero. Surrounding
/// field whitespace is accepted; only the first three colon fields are read.
/// Invalid/empty integer fields raise ValueError. This parser does NOT enforce
/// 24/60/60 bounds, select a timezone, or validate a leap-second label.
#[pyfunction]
fn t2i(py: Python<'_>, time_string: &str) -> PyResult<(i32, i32, i32)> {
    pure(py, || native::parse_hms(time_string))
}
/// Parse an ISO Gregorian label to (year, month, day, hour, minute, second).
///
/// date_string may be a date or timestamp with optional fractional seconds
/// and explicit UTC offset. Returns the LOCAL integral civil fields; fractional
/// seconds are omitted and the offset is not applied to these returned fields.
/// Missing time is midnight. Malformed dates, leap seconds and years outside
/// 1..9999 raise ValueError. No machine timezone or astronomical clock is used.
#[pyfunction]
fn dt2i(py: Python<'_>, date_string: &str) -> PyResult<(i32, u8, u8, u8, u8, u8)> {
    pure(py, || native::parse_datetime(date_string)).map(|value| {
        (
            value.year,
            value.month,
            value.day,
            value.hour,
            value.minute,
            value.second,
        )
    })
}
/// Format finite Gregorian Julian day jd as 'YYYY-MM-DDTHH:MM:SS'.
///
/// Rounds seconds to nearest with ties to even and carries midnight. Returns
/// no timezone suffix and infers/converts no time scale; the label has the
/// caller's JD scale. Nonfinite input or dates outside years 1..9999 raise
/// ValueError. No system clock, leap table or ephemeris is read.
#[pyfunction]
fn jd2isostr(py: Python<'_>, jd: f64) -> PyResult<String> {
    pure(py, || contrib::iso_timestamp(jd))
}
/// Return signed fractional Julian years from ISO date_a to date_b (b-a).
///
/// One year is 365.25 days. Explicit UTC offsets are normalized for subtraction
/// and fractional seconds retained; two naive labels use civil subtraction,
/// without a machine timezone. Mixing offset-aware and naive labels raises
/// TypeError. Malformed/out-of-domain dates or leap-second labels raise
/// ValueError. This is not a calendar anniversary count or TT/UT1 conversion.
#[pyfunction]
fn years_diff(py: Python<'_>, date_a: &str, date_b: &str) -> PyResult<f64> {
    py.detach(|| native::years_diff(date_a, date_b))
        .map_err(|error| match error {
            native::ArabicError::MixedTimezoneAwareness => PyTypeError::new_err(error.to_string()),
            _ => PyValueError::new_err(error.to_string()),
        })
}
/// Parse '12N34', '12S34', '12E34' or '12W34' into signed decimal degrees.
///
/// coord is case-insensitive with outer whitespace stripped. N/E are positive,
/// S/W negative; missing degrees/minutes around a hemisphere letter are zero.
/// A bare float string is accepted unchanged, including nonfinite spellings.
/// Invalid text or nonfinite hemisphere components raises ValueError. This
/// parser does not enforce latitude/longitude or minute bounds.
#[pyfunction]
fn geoc2d(py: Python<'_>, coord: &str) -> PyResult<f64> {
    pure(py, || native::geoc2d(coord))
}
/// Format geographic latitude lat in degrees as 'DDNMM' or 'DDSMM'.
///
/// Minutes round to an integer with ties to even, carrying 60 into degrees;
/// zero is north. Nonfinite lat raises ValueError. Geographic bounds are not
/// enforced; use ordinary latitude magnitudes, since integer-degree/carry
/// overflow is outside this formatter's domain and can panic in debug builds.
#[pyfunction]
fn geolat2c(py: Python<'_>, lat: f64) -> PyResult<String> {
    pure(py, || native::geolat2c(lat))
}
/// Format east-positive longitude lon in degrees as 'DDDEMM' or 'DDDWMM'.
///
/// Minutes round to an integer with ties to even and carry 60 into degrees;
/// zero is east. Nonfinite lon raises ValueError. No longitude normalization
/// or range guard is imposed; extreme integer-degree/carry overflow is outside
/// the domain and can panic in debug builds. No observer setting is changed.
#[pyfunction]
fn geolon2c(py: Python<'_>, lon: f64) -> PyResult<String> {
    pure(py, || native::geolon2c(lon))
}
/// Test an aspect from supplied degree longitudes, returning (hit, deviation).
///
/// a_lon, b_lon, aspect and orb are finite degrees; a_id/b_id are accepted
/// but UNUSED. hit is integer 1 iff abs(shortest_separation-aspect) <= orb,
/// otherwise 0; deviation is that nonnegative absolute residual in degrees.
/// Negative orb is not rejected but cannot match. Nonfinite angle inputs
/// raise ValueError; no body lookup, chart, clock or position is acquired.
#[pyfunction]
fn match_aspect(
    py: Python<'_>,
    a_id: i32,
    a_lon: f64,
    aspect: f64,
    b_id: i32,
    b_lon: f64,
    orb: f64,
) -> PyResult<(i32, f64)> {
    let _ = (a_id, b_id);
    pure(py, || native::match_aspect(a_lon, b_lon, aspect, orb))
        .map(|(hit, delta)| (i32::from(hit), delta))
}
/// Return (integer_hit, deviation_deg) for a_lon/b_lon, aspect and orb in degrees.
///
/// Uses the same shortest-separation/absolute-residual test as match_aspect,
/// without unused body-ID arguments. All four angles must be finite or
/// ValueError is raised. Negative orb yields no match; no chart/data is read.
#[pyfunction]
fn match_aspect2(
    py: Python<'_>,
    a_lon: f64,
    b_lon: f64,
    aspect: f64,
    orb: f64,
) -> PyResult<(i32, f64)> {
    pure(py, || native::match_aspect(a_lon, b_lon, aspect, orb))
        .map(|(hit, delta)| (i32::from(hit), delta))
}
/// Return (integer_hit, deviation_deg, signed_arc_deg) for a supplied aspect.
///
/// a_lon/b_lon, aspect and orb are finite degrees; a_id/b_id are UNUSED.
/// Match/deviation follow match_aspect. signed_arc is b_lon-a_lon wrapped
/// to (-180,180], not a signed aspect residual. Nonfinite angle inputs raise
/// ValueError. No clock, planetary position or chart state is consulted.
#[pyfunction]
fn match_aspect3(
    py: Python<'_>,
    a_id: i32,
    a_lon: f64,
    aspect: f64,
    b_id: i32,
    b_lon: f64,
    orb: f64,
) -> PyResult<(i32, f64, f64)> {
    let _ = (a_id, b_id);
    pure(py, || {
        native::match_aspect_signed(a_lon, b_lon, aspect, orb)
    })
    .map(|(hit, delta, signed)| (i32::from(hit), delta, signed))
}
/// Return (integer_hit, deviation_deg, signed_arc_deg, separation_deg).
///
/// a_lon/b_lon, aspect and orb are finite degrees; a_id/b_id are UNUSED.
/// Match and residual follow match_aspect; signed_arc is b-a in (-180,180],
/// separation is the shortest arc in `0..=180` degrees. Nonfinite angle inputs raise
/// ValueError. This scalar predicate acquires no bodies, time or chart.
#[pyfunction]
fn match_aspect4(
    py: Python<'_>,
    a_id: i32,
    a_lon: f64,
    aspect: f64,
    b_id: i32,
    b_lon: f64,
    orb: f64,
) -> PyResult<(i32, f64, f64, f64)> {
    let _ = (a_id, b_id);
    pure(py, || native::match_aspect_full(a_lon, b_lon, aspect, orb))
        .map(|(hit, delta, signed, sep)| (i32::from(hit), delta, signed, sep))
}
/// Find the first half-day sample matching an aspect to a FIXED longitude.
///
/// body is a public body ID; aspect, target_lon and orb are degrees (orb
/// default 1); jd_start is UT1 JD. Samples begin at jd_start and advance
/// 0.5 day over less than three Julian years using configured ecliptic places.
/// Returns (sample_jd_ut1, absolute_residual_deg), or (NaN, +infinity) for a
/// completed miss, not None. This coarse convenience search is NOT an exact
/// crossing solver and can miss narrow windows. Nonfinite scalar arguments
/// raise ValueError; data, coverage, configuration and calculation errors
/// propagate rather than becoming a search miss.
#[pyfunction(signature=(body,aspect,target_lon,jd_start,orb=1.0))]
fn next_aspect(
    py: Python<'_>,
    body: i32,
    aspect: f64,
    target_lon: f64,
    jd_start: f64,
    orb: f64,
) -> PyResult<(f64, f64)> {
    call(py, |session| {
        session.contrib_next_aspect(
            body,
            contrib::AspectTarget::Longitude(target_lon),
            aspect,
            jd_start,
            orb,
        )
    })
}
/// Search as next_aspect with a fixed one-degree orb and fixed target_lon.
///
/// body is a public ID, aspect/target_lon degrees and jd_start UT1 JD.
/// Returns (sample_jd_ut1, absolute_residual_deg), or (NaN, +infinity) after
/// the same three-year half-day grid. Not an exact crossing; finite-input
/// validation and provider/configuration/coverage errors follow next_aspect.
#[pyfunction]
fn next_aspect2(
    py: Python<'_>,
    body: i32,
    aspect: f64,
    target_lon: f64,
    jd_start: f64,
) -> PyResult<(f64, f64)> {
    call(py, |session| {
        session.contrib_next_aspect(
            body,
            contrib::AspectTarget::Longitude(target_lon),
            aspect,
            jd_start,
            1.0,
        )
    })
}
/// Search for a sampled aspect to the FIXED supplied cusp longitude.
///
/// body is a public ID; aspect/cusp_lon/orb are degrees, orb defaults to 1;
/// jd_start is UT1 JD. No house/cusp is recomputed. Returns (sample_jd_ut1,
/// absolute_residual_deg), or (NaN, +infinity), on next_aspect's half-day grid
/// over less than three years. Finite-input/provider errors propagate. This
/// is a coarse convenience search, not a moving-cusp or exact event solver.
#[pyfunction(signature=(body,aspect,cusp_lon,jd_start,orb=1.0))]
fn next_aspect_cusp(
    py: Python<'_>,
    body: i32,
    aspect: f64,
    cusp_lon: f64,
    jd_start: f64,
    orb: f64,
) -> PyResult<(f64, f64)> {
    call(py, |session| {
        session.contrib_next_aspect(
            body,
            contrib::AspectTarget::Longitude(cusp_lon),
            aspect,
            jd_start,
            orb,
        )
    })
}
/// Search next_aspect_cusp's fixed cusp with a fixed one-degree orb.
///
/// body is a public ID; aspect/cusp_lon are degrees and jd_start UT1 JD.
/// Returns (sample_jd_ut1, absolute_residual_deg), or (NaN, +infinity), over
/// the same three-year half-day grid. No house is recalculated; nonfinite
/// inputs and provider/configuration/coverage failures raise exceptions.
#[pyfunction]
fn next_aspect_cusp2(
    py: Python<'_>,
    body: i32,
    aspect: f64,
    cusp_lon: f64,
    jd_start: f64,
) -> PyResult<(f64, f64)> {
    call(py, |session| {
        session.contrib_next_aspect(
            body,
            contrib::AspectTarget::Longitude(cusp_lon),
            aspect,
            jd_start,
            1.0,
        )
    })
}
/// Find the first half-day sampled aspect between TWO moving bodies.
///
/// body_a/body_b are public IDs; aspect and orb are degrees (orb default 1);
/// jd_start is UT1 JD. Both configured ecliptic places are evaluated at each
/// sample. Returns (sample_jd_ut1, absolute_residual_deg), or (NaN,+infinity)
/// for a completed three-year-grid miss. No exact crossing refinement occurs;
/// narrow windows can be missed. Nonfinite scalar and provider/coverage/
/// configuration failures propagate, not a silent fixed-target fallback.
#[pyfunction(signature=(body_a,body_b,aspect,jd_start,orb=1.0))]
fn next_aspect_with(
    py: Python<'_>,
    body_a: i32,
    body_b: i32,
    aspect: f64,
    jd_start: f64,
    orb: f64,
) -> PyResult<(f64, f64)> {
    call(py, |session| {
        session.contrib_next_aspect(
            body_a,
            contrib::AspectTarget::Body(body_b),
            aspect,
            jd_start,
            orb,
        )
    })
}
/// Search next_aspect_with's two moving bodies with a fixed one-degree orb.
///
/// body_a/body_b are public IDs, aspect degrees and jd_start UT1 JD. Returns
/// (sample_jd_ut1, absolute_residual_deg) or (NaN,+infinity) on the three-year
/// half-day grid. It is not an exact crossing; nonfinite scalar and required
/// provider/configuration/coverage errors propagate.
#[pyfunction]
fn next_aspect_with2(
    py: Python<'_>,
    body_a: i32,
    body_b: i32,
    aspect: f64,
    jd_start: f64,
) -> PyResult<(f64, f64)> {
    call(py, |session| {
        session.contrib_next_aspect(
            body_a,
            contrib::AspectTarget::Body(body_b),
            aspect,
            jd_start,
            1.0,
        )
    })
}
/// Find the first half-day sample changing to negative longitude speed.
///
/// body is a public ID; jd_start is UT1 JD; max_days defaults to 1095.75
/// days. The previous sample must have nonnegative speed, so an initially
/// retrograde body is not itself a transition. Returns sampled UT1 JD or
/// NaN for no detected change (also for max_days <= 0), not None. No exact
/// station is refined. Nonfinite start/horizon raises ValueError; actual
/// data, coverage, configuration and calculation failures propagate.
#[pyfunction(signature=(body,jd_start,max_days=1095.75))]
fn next_retro(py: Python<'_>, body: i32, jd_start: f64, max_days: f64) -> PyResult<f64> {
    call(py, |session| {
        session.contrib_next_retro(body, jd_start, max_days)
    })
}

macro_rules! unsupported {
    ($name:ident,$doc:literal) => {
        #[doc = $doc]
        #[pyfunction(signature=(*_args,**_kwargs))]
        fn $name(_args: &Bound<'_, PyTuple>, _kwargs: Option<&Bound<'_, PyDict>>) -> PyResult<()> {
            Err(PyNotImplementedError::new_err(concat!(
                "libephemeris.contrib.",
                stringify!($name),
                " is not implemented"
            )))
        }
    };
}
unsupported!(
    atlas_connect,
    "Unsupported atlas connection compatibility entry.\n\nAccepts arbitrary positional/keyword arguments and always raises\nNotImplementedError. Opens no connection, files or network resources;\nthere is no successful return value."
);
unsupported!(
    atlas_close,
    "Unsupported atlas close compatibility entry.\n\nAccepts arbitrary positional/keyword arguments and always raises\nNotImplementedError; it is not a no-op successful close and retires\nno runtime or atlas resources."
);
unsupported!(
    atlas_search,
    "Unsupported atlas search compatibility entry.\n\nAccepts arbitrary positional/keyword arguments and always raises\nNotImplementedError, not None or an empty result. Performs no atlas\nlookup, database query, file access or download."
);
unsupported!(
    atlas_countries_list,
    "Unsupported atlas country-list compatibility entry.\n\nAccepts arbitrary positional/keyword arguments and always raises\nNotImplementedError, not an empty list. No atlas/database is acquired."
);
unsupported!(
    db_connect,
    "Unsupported database connection compatibility entry.\n\nAccepts arbitrary positional/keyword arguments and always raises\nNotImplementedError. Opens no database and has no successful return value."
);
unsupported!(
    db_close,
    "Unsupported database close compatibility entry.\n\nAccepts arbitrary positional/keyword arguments and always raises\nNotImplementedError; no database, context or resource cache is closed."
);
unsupported!(
    tzabbr_find,
    "Unsupported timezone-abbreviation lookup compatibility entry.\n\nAccepts arbitrary positional/keyword arguments and always raises\nNotImplementedError, not None. No timezone database or network is consulted."
);
unsupported!(
    tzabbr_list,
    "Unsupported timezone-abbreviation list compatibility entry.\n\nAccepts arbitrary positional/keyword arguments and always raises\nNotImplementedError, not an empty list. No timezone data is acquired."
);

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module
        .py()
        .get_type::<Error>()
        .setattr("__module__", "libephemeris.contrib")?;
    module.add("Error", module.py().get_type::<Error>())?;
    module.add_class::<Data>()?;
    module.add_class::<User>()?;
    module.add_function(wrap_pyfunction!(atlas_connect, module)?)?;
    module.add_function(wrap_pyfunction!(atlas_close, module)?)?;
    module.add_function(wrap_pyfunction!(atlas_search, module)?)?;
    module.add_function(wrap_pyfunction!(atlas_countries_list, module)?)?;
    module.add_function(wrap_pyfunction!(db_connect, module)?)?;
    module.add_function(wrap_pyfunction!(db_close, module)?)?;
    module.add_function(wrap_pyfunction!(tzabbr_find, module)?)?;
    module.add_function(wrap_pyfunction!(tzabbr_list, module)?)?;

    module.add_function(wrap_pyfunction!(long2rasi, module)?)?;
    module.add_function(wrap_pyfunction!(long2nakshatra, module)?)?;
    module.add_function(wrap_pyfunction!(long2navamsa, module)?)?;
    module.add_function(wrap_pyfunction!(rasinorm, module)?)?;
    module.add_function(wrap_pyfunction!(rasi_dif, module)?)?;
    module.add_function(wrap_pyfunction!(rasi_dif2, module)?)?;
    module.add_function(wrap_pyfunction!(antiscion, module)?)?;
    module.add_function(wrap_pyfunction!(degsplit, module)?)?;
    module.add_function(wrap_pyfunction!(signtostr, module)?)?;
    module.add_function(wrap_pyfunction!(get_nakshatra_name, module)?)?;
    module.add_function(wrap_pyfunction!(house_system_name, module)?)?;
    module.add_function(wrap_pyfunction!(lord, module)?)?;
    module.add_function(wrap_pyfunction!(ochchabala, module)?)?;
    module.add_function(wrap_pyfunction!(residential_strength, module)?)?;
    module.add_function(wrap_pyfunction!(naisargika_relation, module)?)?;
    module.add_function(wrap_pyfunction!(tatkalika_relation, module)?)?;
    module.add_function(wrap_pyfunction!(raman_houses, module)?)?;
    module.add_function(wrap_pyfunction!(saturn_4_stars, module)?)?;
    module.add_function(wrap_pyfunction!(jdnow, module)?)?;
    module.add_function(wrap_pyfunction!(jduration, module)?)?;
    module.add_function(wrap_pyfunction!(t2i, module)?)?;
    module.add_function(wrap_pyfunction!(dt2i, module)?)?;
    module.add_function(wrap_pyfunction!(jd2isostr, module)?)?;
    module.add_function(wrap_pyfunction!(years_diff, module)?)?;
    module.add_function(wrap_pyfunction!(geoc2d, module)?)?;
    module.add_function(wrap_pyfunction!(geolat2c, module)?)?;
    module.add_function(wrap_pyfunction!(geolon2c, module)?)?;
    module.add_function(wrap_pyfunction!(match_aspect, module)?)?;
    module.add_function(wrap_pyfunction!(match_aspect2, module)?)?;
    module.add_function(wrap_pyfunction!(match_aspect3, module)?)?;
    module.add_function(wrap_pyfunction!(match_aspect4, module)?)?;
    module.add_function(wrap_pyfunction!(next_aspect, module)?)?;
    module.add_function(wrap_pyfunction!(next_aspect2, module)?)?;
    module.add_function(wrap_pyfunction!(next_aspect_cusp, module)?)?;
    module.add_function(wrap_pyfunction!(next_aspect_cusp2, module)?)?;
    module.add_function(wrap_pyfunction!(next_aspect_with, module)?)?;
    module.add_function(wrap_pyfunction!(next_aspect_with2, module)?)?;
    module.add_function(wrap_pyfunction!(next_retro, module)?)?;
    module.add("ARIES", native::ARIES)?;
    module.add("TAURUS", native::TAURUS)?;
    module.add("GEMINI", native::GEMINI)?;
    module.add("CANCER", native::CANCER)?;
    module.add("LEO", native::LEO)?;
    module.add("VIRGO", native::VIRGO)?;
    module.add("LIBRA", native::LIBRA)?;
    module.add("SCORPIO", native::SCORPIO)?;
    module.add("SAGITTARIUS", native::SAGITTARIUS)?;
    module.add("CAPRICORN", native::CAPRICORN)?;
    module.add("AQUARIUS", native::AQUARIUS)?;
    module.add("PISCES", native::PISCES)?;
    module.add("MESHA", native::MESHA)?;
    module.add("VRISHABA", native::VRISHABA)?;
    module.add("MITHUNA", native::MITHUNA)?;
    module.add("KATAKA", native::KATAKA)?;
    module.add("SIMHA", native::SIMHA)?;
    module.add("KANYA", native::KANYA)?;
    module.add("THULA", native::THULA)?;
    module.add("VRISHIKA", native::VRISHIKA)?;
    module.add("DHANUS", native::DHANUS)?;
    module.add("MAKARA", native::MAKARA)?;
    module.add("KUMBHA", native::KUMBHA)?;
    module.add("MEENA", native::MEENA)?;
    module.add("SURYA", native::SURYA)?;
    module.add("RAVI", native::RAVI)?;
    module.add("CHANDRA", native::CHANDRA)?;
    module.add("SOMA", native::SOMA)?;
    module.add("BUDHA", native::BUDHA)?;
    module.add("SOUMYA", native::SOUMYA)?;
    module.add("THAMA", native::THAMA)?;
    module.add("SUKRA", native::SUKRA)?;
    module.add("BHARGAVA", native::BHARGAVA)?;
    module.add("KUJA", native::KUJA)?;
    module.add("ANGARAKA", native::ANGARAKA)?;
    module.add("GURU", native::GURU)?;
    module.add("BRIHASPATI", native::BRIHASPATI)?;
    module.add("SANI", native::SANI)?;
    module.add("MANDA", native::MANDA)?;
    module.add("RAHU", native::RAHU)?;
    module.add("KETU", native::KETU)?;
    module.add("SIKHI", native::SIKHI)?;
    module.add("CONJUNCTION", native::CONJUNCTION as i64)?;
    module.add("SEMISEXTILE", native::SEMISEXTILE as i64)?;
    module.add("SEXTILE", native::SEXTILE as i64)?;
    module.add("SQUARE", native::SQUARE as i64)?;
    module.add("TRINE", native::TRINE as i64)?;
    module.add("QUINCUNX", native::QUINCUNX as i64)?;
    module.add("OPPOSITION", native::OPPOSITION as i64)?;
    module.add("SEMISQUARE", native::SEMISQUARE as i64)?;
    module.add("SESQUISQUARE", native::SESQUISQUARE as i64)?;
    module.add("QUINTILE", native::QUINTILE as i64)?;
    module.add("BIQUINTILE", native::BIQUINTILE as i64)?;
    module.add("SEMIQUINTILE", native::SEMIQUINTILE as i64)?;
    module.add("SEPTILE", native::SEPTILE)?;
    module.add("BISEPTILE", native::BISEPTILE)?;
    module.add("TRISEPTILE", native::TRISEPTILE)?;
    module.add("NOVILE", native::NOVILE as i64)?;
    module.add("SEMINOVILE", native::SEMINOVILE as i64)?;
    module.add("BINOVILE", native::BINOVILE as i64)?;
    module.add("QUATRONOVILE", native::QUATRONOVILE as i64)?;
    module.add("UNDECILE", native::UNDECILE)?;
    module.add("BIUNDECILE", native::BIUNDECILE)?;
    module.add("TRIUNDECILE", native::TRIUNDECILE)?;
    module.add("QUADUNDECILE", native::QUADUNDECILE)?;
    module.add("QUINUNDECILE", native::QUINUNDECILE)?;
    module.add("SQUISEXTILE", native::SQUISEXTILE as i64)?;
    module.add("SQUISQUARE", native::SQUISQUARE)?;
    for &(name, id) in contrib::NAKSHATRA_IDS {
        module.add(name, id)?;
    }
    Ok(())
}
