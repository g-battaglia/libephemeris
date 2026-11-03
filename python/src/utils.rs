//! Thin numerical and coordinate conversions; every formula is in the core.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::{python_support as api, utils};
use pyo3::exceptions::PyUnicodeDecodeError;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyString};

type Triple = (f64, f64, f64);
type Sextuple = (f64, f64, f64, f64, f64, f64);

fn triple(value: [f64; 3]) -> Triple {
    let [a, b, c] = value;
    (a, b, c)
}

/// Decode only the two public text representations at the language boundary.
pub(crate) fn text(value: &Bound<'_, PyAny>) -> PyResult<String> {
    if let Ok(value) = value.cast::<PyString>() {
        return Ok(value.to_str()?.to_owned());
    }
    let value = value.cast::<PyBytes>()?;
    let bytes = value.as_bytes();
    if let Some(index) = bytes.iter().position(|byte| !byte.is_ascii()) {
        return Err(PyUnicodeDecodeError::new_err((
            "ascii",
            value.clone().unbind(),
            index,
            index + 1,
            "ordinal not in range(128)",
        )));
    }
    Ok(String::from_utf8_lossy(bytes).into_owned())
}

macro_rules! float_one {
    ($name:ident, $arg:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyfunction]
        fn $name($arg: f64) -> f64 {
            utils::$name($arg)
        }
    };
}
macro_rules! float_two {
    ($name:ident, $a:ident, $b:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyfunction]
        fn $name($a: f64, $b: f64) -> f64 {
            utils::$name($a, $b)
        }
    };
}
macro_rules! integer_one {
    ($name:ident, $arg:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyfunction]
        fn $name($arg: i64) -> i64 {
            utils::$name($arg)
        }
    };
}
macro_rules! integer_two {
    ($name:ident, $a:ident, $b:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyfunction]
        fn $name($a: i64, $b: i64) -> i64 {
            utils::$name($a, $b)
        }
    };
}

float_one!(
    degnorm,
    x,
    r"Normalize x degrees to [0, 360).

Returns a float; negative angles wrap forward and magnitudes below 1e-13
degrees snap to zero. Nonfinite inputs propagate as NaN."
);
float_one!(
    radnorm,
    x,
    r"Normalize x radians to [0, 2*pi).

Returns a float; negative angles wrap forward and magnitudes below 1e-13
radians snap to zero. Nonfinite inputs propagate as NaN."
);
float_two!(
    difdeg2n,
    p1,
    p2,
    r"Return the signed angular difference p1 - p2 in degrees.

Inputs are degrees; the float result lies in [-180, 180), with exact
half-turns returning -180. Nonfinite arithmetic propagates as NaN."
);
float_two!(
    difdegn,
    p1,
    p2,
    r"Return p1 - p2 reduced to [0, 360) degrees.

Both arguments and the float result are degrees. Unlike difdeg2n, this
unsigned difference wraps negative differences forward."
);
float_two!(
    difrad2n,
    p1,
    p2,
    r"Return the signed angular difference p1 - p2 in radians.

Inputs are radians; the float result lies in [-pi, pi), with exact
half-turns returning -pi. Nonfinite arithmetic propagates as NaN."
);
float_two!(
    deg_midp,
    x1,
    x2,
    r"Return the shorter-arc midpoint of degree angles x1 and x2.

The float result is normalized to [0, 360). Exact oppositions take the
positive half-turn from x1; nonfinite arithmetic propagates as NaN."
);
float_two!(
    rad_midp,
    x,
    y,
    r"Return the shorter-arc midpoint of radian angles x and y.

The float result is normalized to [0, 2*pi). Exact oppositions take the
positive half-turn from x; nonfinite arithmetic propagates as NaN."
);
integer_one!(
    csnorm,
    cs,
    r"Normalize an integer angle in centiseconds to [0, 129600000).

One centisecond is 0.01 arcsecond. The returned integer has the same unit;
negative values wrap forward. Inputs must fit a signed 64-bit integer."
);
integer_one!(
    csroundsec,
    cs,
    r"Round centisecond angle cs to arcsecond display precision.

Returns centiseconds, normally a multiple of 100, with the native
30-degree zodiac-boundary guard. The +50 bias truncates toward zero:
150 becomes 200, while -150 becomes -100. Not Python round()."
);
integer_two!(
    difcs2n,
    p1,
    p2,
    r"Return the signed difference p1 - p2 in angular centiseconds.

One centisecond is 0.01 arcsecond. Returns an integer in
[-64800000, 64800000); exact half-turns select the negative boundary."
);
integer_two!(
    difcsn,
    p1,
    p2,
    r"Return p1 - p2 reduced to [0, 129600000) angular centiseconds.

Inputs and the integer result use 0.01 arcsecond units. Unlike difcs2n,
negative differences wrap forward; inputs must fit signed 64-bit integers."
);

/// Round d to an integer with half-way cases away from zero.
///
/// Unlike Python round(), d2l(1.5) is 2 and d2l(-1.5) is -2. Returns a
/// signed-64-bit-range Python int. NaN raises ValueError; infinity or a
/// rounded result outside that range raises OverflowError.
#[pyfunction]
fn d2l(py: Python<'_>, d: f64) -> PyResult<i64> {
    utils::d2l_checked(d).map_err(|error| to_pyerr(py, error))
}

/// Format angular centiseconds cs as degrees, minutes and seconds.
///
/// One centisecond is 0.01 arcsecond. Returns a space-padded str with a
/// hundredth-of-arcsecond field; the sign belongs to the degree field, so
/// negative sub-degree magnitudes lose their sign. Full turns are not
/// normalized. cs must fit a signed 64-bit integer.
#[pyfunction]
fn cs2degstr(cs: i64) -> String {
    utils::cs2degstr(cs)
}

/// Format angular centiseconds with a positive or negative direction marker.
///
/// cs is an integer in 0.01 arcsecond units; plus and minus are str or ASCII
/// bytes such as 'E'/'W' or 'N'/'S'. Zero uses plus. Returns str with whole
/// arcseconds rounded and zero seconds omitted; rounding carries into minutes
/// and degrees. Non-ASCII bytes raise UnicodeDecodeError.
#[pyfunction]
fn cs2lonlatstr(cs: i64, plus: &Bound<'_, PyAny>, minus: &Bound<'_, PyAny>) -> PyResult<String> {
    Ok(utils::cs2lonlatstr(cs, &text(plus)?, &text(minus)?))
}

/// Format centiseconds of time cs as a 24-hour clock string.
///
/// Here cs means 0.01 second of time, not angular centiseconds. sep is str
/// or ASCII bytes, default ':'. Returns zero-padded HH:MM:SS, with native
/// nearest-second rounding and wrapping at 24 hours, including negative
/// inputs. suppresszero=True omits a zero seconds field. Non-ASCII sep
/// bytes raise UnicodeDecodeError; cs must fit a signed 64-bit integer.
#[pyfunction(signature = (cs, sep = ":".to_owned(), suppresszero = false), text_signature = "(cs, sep=':', suppresszero=False)")]
fn cs2timestr(cs: i64, #[pyo3(from_py_with = text)] sep: String, suppresszero: bool) -> String {
    utils::cs2timestr(cs, &sep, suppresszero)
}

/// Rotate a spherical position between equatorial and ecliptic coordinates.
///
/// coord is a three-number sequence (longitude_or_RA, latitude_or_Dec,
/// distance). Angles and obliquity eps are degrees; RA is not hours.
/// Positive eps converts equatorial to ecliptic; negative eps reverses the
/// rotation. Returns a three-float tuple with unchanged distance/unit. This
/// is a fixed rotation, not precession or apparent-place correction.
/// Nonfinite components or eps raise ValueError.
///
/// Examples
/// --------
/// ```python
/// from libephemeris import cotrans
/// lon, lat, radius = cotrans((120.0, 10.0, 1.0), 23.44)
/// ```
#[pyfunction]
fn cotrans(py: Python<'_>, coord: [f64; 3], eps: f64) -> PyResult<Triple> {
    py.detach(|| api::cotrans(coord[0], coord[1], coord[2], eps))
        .map(triple)
        .map_err(|error| to_pyerr(py, error))
}

/// Rotate a spherical position and its rates using either supported layout.
///
/// Parameters
/// ----------
/// coord : sequence of float
///     Six components (lon, lat, distance, lon_rate, lat_rate, radial_rate),
///     or three position components when rates are supplied separately.
/// eps_or_speed : float or sequence of float
///     Obliquity in degrees for the six-component layout, or a three-rate
///     sequence for the separate layout.
/// eps : float or None, optional
///     Obliquity in degrees for separate position/rate sequences. None selects
///     the six-component layout. Positive obliquity converts equatorial to
///     ecliptic; negative converts back.
///
/// Returns
/// -------
/// tuple
///     Six floats for cotrans_sp(six_components, eps), or (position_tuple,
///     rate_tuple) for cotrans_sp(position, rates, eps). Angles/rates are
///     degrees/degrees per day; distance/radial rate retain their length unit.
///     No physical clock conversion is performed.
///
/// Raises
/// ------
/// ValueError
///     For nonfinite numerical components or incompatible sequence lengths.
/// TypeError
///     If an argument cannot be extracted in the selected layout.
#[pyfunction(signature = (coord, eps_or_speed, eps = None))]
fn cotrans_sp(
    py: Python<'_>,
    coord: &Bound<'_, PyAny>,
    eps_or_speed: &Bound<'_, PyAny>,
    eps: Option<f64>,
) -> PyResult<Py<PyAny>> {
    if let Some(eps) = eps {
        let [a, b, c] = coord.extract::<[f64; 3]>()?;
        let [d, e, f] = eps_or_speed.extract::<[f64; 3]>()?;
        let (position, speed) = py
            .detach(|| api::cotrans_sp(a, b, c, d, e, f, eps))
            .map_err(|error| to_pyerr(py, error))?;
        Ok((triple(position), triple(speed))
            .into_pyobject(py)?
            .into_any()
            .unbind())
    } else {
        let [a, b, c, d, e, f] = coord.extract::<[f64; 6]>()?;
        let eps = eps_or_speed.extract::<f64>()?;
        let ([a, b, c], [d, e, f]) = py
            .detach(|| api::cotrans_sp(a, b, c, d, e, f, eps))
            .map_err(|error| to_pyerr(py, error))?;
        let result: Sextuple = (a, b, c, d, e, f);
        Ok(result.into_pyobject(py)?.into_any().unbind())
    }
}

/// Split degree into degree/minute/second fields and a sign or segment index.
///
/// Returns (degrees, minutes, seconds, fraction, index). With roundflag=0,
/// fields describe abs(degree), fraction is fractional arcseconds, and index
/// is +1 or -1. SPLIT_DEG_ZODIACAL selects a zero-based 30-degree segment;
/// SPLIT_DEG_NAKSHATRA selects 13-degree-20-minute segments for nonnegative
/// inputs. Named SPLIT_DEG_ROUND_* flags choose rounding; KEEP_DEG/KEEP_SIGN
/// prevent the corresponding boundary carry. With rounding, fraction repeats
/// the rounded integer seconds field rather than a subsecond remainder.
/// Angles are not generally normalized before segment indexing. NaN raises
/// ValueError; infinite or out-of-range integer magnitudes raise OverflowError.
#[pyfunction(signature = (degree, roundflag = 0))]
fn split_deg(py: Python<'_>, degree: f64, roundflag: i32) -> PyResult<(i64, i64, i64, f64, i64)> {
    utils::split_deg_checked(degree, roundflag).map_err(|error| to_pyerr(py, error))
}

/// Convert celestial coordinates to local horizontal coordinates at UT1 JD.
///
/// tjdut is UT1 Julian day. flag is ECL2HOR or EQU2HOR; xin is
/// (ecliptic_longitude, latitude, distance) or (RA, Dec, distance), with
/// celestial angles in degrees (RA not hours). Distance is ignored.
/// geopos is (longitude_east_deg, latitude_north_deg, height_m); atpress
/// is mbar (zero estimates pressure from height), and attemp is Celsius.
/// Returns (azimuth, true_altitude, apparent_altitude), all degrees.
/// Azimuth runs from South through West: 0=South, 90=West, 180=North.
/// Apparent altitude includes native refraction; true altitude is geometric.
/// Frame/clock errors propagate as Error subclasses; checked nonfinite
/// projection inputs/results raise ValueError.
#[pyfunction]
fn azalt(
    py: Python<'_>,
    tjdut: f64,
    flag: i32,
    geopos: [f64; 3],
    atpress: f64,
    attemp: f64,
    xin: [f64; 3],
) -> PyResult<Triple> {
    let [lon, lat, alt] = geopos;
    let [ra, dec, _] = xin;
    py.detach(|| {
        runtime()
            .model_session()?
            .azalt(tjdut, flag, lon, lat, alt, atpress, attemp, ra, dec)
    })
    .map(|value| (value.azimuth, value.true_altitude, value.apparent_altitude))
    .map_err(|error| to_pyerr(py, error))
}

/// Convert geometric horizontal coordinates back to celestial coordinates.
///
/// tjdut is UT1 JD; geopos is (longitude_east_deg, latitude_north_deg,
/// height_m), with height accepted but ignored by this inverse kernel.
/// azimuth is degrees South through West and true_altitude is unrefracted
/// degrees. flag=HOR2EQU returns (RA, Dec); HOR2ECL returns ecliptic
/// (longitude, latitude), both in degrees. RA is not hours. This does not
/// remove refraction from an apparent altitude. Frame/clock failures
/// propagate; checked nonfinite numerical inputs raise ValueError.
#[pyfunction]
fn azalt_rev(
    py: Python<'_>,
    tjdut: f64,
    flag: i32,
    geopos: [f64; 3],
    azimuth: f64,
    true_altitude: f64,
) -> PyResult<(f64, f64)> {
    let [lon, lat, alt] = geopos;
    py.detach(|| {
        runtime()
            .model_session()?
            .azalt_rev(tjdut, flag, lon, lat, alt, azimuth, true_altitude)
    })
    .map(|value| (value.longitude_or_ra, value.latitude_or_dec))
    .map_err(|error| to_pyerr(py, error))
}

/// Convert true/apparent altitude with the native closed-form refraction model.
///
/// alt is degrees; atpress is mbar and attemp is Celsius. flag=0 converts
/// true to apparent altitude; any other value selects apparent to true.
/// Returns the converted altitude in degrees, not the refraction increment.
/// Horizon/low-altitude limits and IEEE values follow the model; this scalar
/// helper does not perform a separate physical-input validation pass.
#[pyfunction(signature = (alt, atpress = 1013.25, attemp = 15.0, flag = 0))]
fn refrac(py: Python<'_>, alt: f64, atpress: f64, attemp: f64, flag: i32) -> f64 {
    py.detach(|| runtime().refrac(alt, atpress, attemp, flag))
}

/// Convert altitude with refraction and horizon-dip diagnostic values.
///
/// alt is degrees, geoalt is height in metres, atpress is mbar, attemp is
/// Celsius and lapserate is K/m. None uses the configured lapse rate.
/// flag=0 converts true to apparent; other values select the inverse.
/// Returns (converted_altitude, (true_altitude, apparent_altitude,
/// refraction, horizon_dip)), all angles in degrees. Native low-altitude
/// conventions and IEEE values are retained without a new validation pass.
#[pyfunction(signature = (alt, geoalt, atpress = 1013.25, attemp = 15.0, lapserate = None, flag = 0))]
fn refrac_extended(
    py: Python<'_>,
    alt: f64,
    geoalt: f64,
    atpress: f64,
    attemp: f64,
    lapserate: Option<f64>,
    flag: i32,
) -> (f64, (f64, f64, f64, f64)) {
    py.detach(|| runtime().refrac_extended(alt, geoalt, atpress, attemp, lapserate, flag))
}

/// Calculate and cache eight chart angles and four planetary longitudes.
///
/// jd_ut is UT1 Julian day; lat/lon are geographic degrees north/east.
/// Returns a dict of degree values keyed by Ascendant, MC, ARMC, Vertex,
/// Equatorial_Ascendant, Descendant, IC, AntiVertex, Sun, Moon, Mercury
/// and Venus. Planetary longitudes use ordinary calc_ut calls with flags=0;
/// polar house geometry uses the Porphyry fallback.
///
/// Notes
/// -----
/// This function sets the default runtime's observer to (lon, lat, 0 metres)
/// before calculation and publishes the angle cache on success. It is not
/// the separate eight-angle helper at libephemeris.angles.calc_angles.
/// CoordinateError, coverage/data/clock errors and checked numerical
/// ValueError failures propagate; it is not a side-effect-free geometry call.
#[pyfunction]
fn calc_angles(
    py: Python<'_>,
    jd_ut: f64,
    lat: f64,
    lon: f64,
) -> PyResult<std::collections::HashMap<String, f64>> {
    py.detach(|| runtime().calc_angles(jd_ut, lat, lon))
        .map_err(|error| to_pyerr(py, error))
}

/// Return the great-circle separation of two spherical positions in degrees.
///
/// lon1/lat1 and lon2/lat2 are longitude/latitude pairs in degrees in the
/// same frame; equatorial RA must also be degrees, not hours. Returns a
/// float using the native Vincenty spherical formula, stable near coincidence
/// and opposition. No epoch/frame conversion is performed. Checked nonfinite
/// operands raise ValueError.
#[pyfunction]
fn angular_separation(py: Python<'_>, lon1: f64, lat1: f64, lon2: f64, lat2: f64) -> PyResult<f64> {
    py.detach(|| libephemeris_core::events::vincenty_separation_deg(lon1, lat1, lon2, lat2))
        .map_err(|error| to_pyerr(py, error))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let helpers = PyModule::new(module.py(), "libephemeris._core._utils")?;
    helpers.add_function(wrap_pyfunction!(angular_separation, &helpers)?)?;
    module
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item("libephemeris._core._utils", &helpers)?;
    module.add("_utils", helpers)?;
    macro_rules! add {
        ($($name:ident),+ $(,)?) => { $(module.add_function(wrap_pyfunction!($name, module)?)?;)+ };
    }
    add!(
        degnorm,
        radnorm,
        difdeg2n,
        difdegn,
        difrad2n,
        difcs2n,
        difcsn,
        csnorm,
        csroundsec,
        cs2degstr,
        cs2lonlatstr,
        cs2timestr,
        deg_midp,
        rad_midp,
        d2l,
        split_deg,
        calc_angles,
        cotrans,
        cotrans_sp,
        azalt,
        azalt_rev,
        refrac,
        refrac_extended
    );
    Ok(())
}
