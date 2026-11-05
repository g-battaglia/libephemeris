//! House-system bindings over native geometry and session-owned clocks.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::houses as core_houses;
use pyo3::prelude::*;
use pyo3::types::{PyBytes, PyDict, PyInt, PyString, PyTuple};

type Layout = (Py<PyTuple>, Py<PyTuple>);
type Rates = (Py<PyTuple>, Py<PyTuple>, Py<PyTuple>, Py<PyTuple>);
type Fallback = (Py<PyTuple>, Py<PyTuple>, bool, Option<String>);

fn selector(value: &Bound<'_, PyAny>) -> PyResult<char> {
    let result = if let Ok(text) = value.cast::<PyString>() {
        core_houses::hsys_from_str(text.to_str()?)
    } else if let Ok(bytes) = value.cast::<PyBytes>() {
        core_houses::hsys_from_bytes(bytes.as_bytes())
    } else {
        Ok(core_houses::hsys_from_int(value.extract::<i32>()?))
    };
    result.map_err(|error| to_pyerr(value.py(), error))
}

fn layout(py: Python<'_>, value: (Vec<f64>, [f64; 8])) -> PyResult<Layout> {
    Ok((
        PyTuple::new(py, value.0)?.unbind(),
        PyTuple::new(py, value.1)?.unbind(),
    ))
}

fn rates(py: Python<'_>, value: (Vec<f64>, [f64; 8], Vec<f64>, [f64; 8])) -> PyResult<Rates> {
    Ok((
        PyTuple::new(py, value.0)?.unbind(),
        PyTuple::new(py, value.1)?.unbind(),
        PyTuple::new(py, value.2)?.unbind(),
        PyTuple::new(py, value.3)?.unbind(),
    ))
}

/// Return the rounded polar latitude threshold ``90 - obliquity`` in degrees.
///
/// obliquity is degrees, default 23.44. This helper does not validate inputs
/// or choose a house system; NaN/infinity follow floating-point arithmetic.
/// The displayed threshold is rounded: get_extreme_latitude_info classifies
/// the actual latitude/obliquity pair without relying on this subtraction.
#[pyfunction(signature = (obliquity = 23.44))]
fn get_polar_latitude_threshold(obliquity: f64) -> f64 {
    core_houses::polar_latitude_threshold_deg(obliquity)
}

/// Describe house-system latitude limits without calculating cusps.
///
/// lat is north-positive geographic degrees; obliquity is degrees (23.44 by
/// default). Returns a dict with keys ``latitude``, ``is_extreme`` (abs(lat)
/// >= 80), ``is_polar_circle`` (abs(lat) + obliquity > 90 using the exact
/// supplied floating-point pair), ``polar_threshold``, ``extreme_threshold``,
/// ``hemisphere`` ('N' or 'S'), ``affected_systems``, ``unstable_systems`` and
/// ``stable_systems`` (lists of selector strings). Thresholds are degrees.
/// This advisory helper imposes no finite/range validation and does not select
/// a fallback; exact equality at the polar boundary is not beyond the circle.
#[pyfunction(signature = (lat, obliquity = 23.44))]
fn get_extreme_latitude_info(py: Python<'_>, lat: f64, obliquity: f64) -> PyResult<Py<PyDict>> {
    let value = core_houses::extreme_latitude_info(lat, obliquity);
    let result = PyDict::new(py);
    result.set_item("latitude", value.latitude)?;
    result.set_item("is_extreme", value.is_extreme)?;
    result.set_item("is_polar_circle", value.is_polar_circle)?;
    result.set_item("polar_threshold", value.polar_threshold)?;
    result.set_item("extreme_threshold", value.extreme_threshold)?;
    result.set_item("hemisphere", value.hemisphere.to_string())?;
    result.set_item(
        "affected_systems",
        value
            .affected_systems
            .into_iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>(),
    )?;
    result.set_item(
        "unstable_systems",
        value
            .unstable_systems
            .into_iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>(),
    )?;
    result.set_item(
        "stable_systems",
        value
            .stable_systems
            .into_iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>(),
    )?;
    Ok(result.unbind())
}

/// Calculate tropical house cusps and eight auxiliary angles at a UT1 epoch.
///
/// Parameters
/// ----------
/// tjdut : float
///     UT1 Julian day, not UTC or TT.
/// lat, lon : float
///     Finite geographic latitude [-90, 90] north and longitude [-180, 360]
///     east, in degrees. This query does not change the configured observer.
/// hsys : int, str or bytes
///     Character code, one-character string or one Latin-1 byte; default 'P'
///     (80, Placidus). ASCII lowercase folds to uppercase except 'i' (a
///     distinct Sunshine construction) and 'g'. Unknown selectors use the
///     Placidus construction, including its polar refusal.
/// iflag : int
///     Default 0. Ephemeris bits select Sunshine's Sun-declination source;
///     this tropical entry does not apply SIDEREAL or RADIANS formatting.
///
/// Returns
/// -------
/// (cusps, angles) : tuple[tuple[float, ...], tuple[float, ...]]
///     Ecliptic longitudes in degrees. Cusps have 12 entries, index 0 = house
///     1, without padding; uppercase 'G' instead returns 36 Gauquelin sectors,
///     lowercase 'g' only their first 12. The eight angles are Ascendant, MC,
///     ARMC, Vertex, equatorial Ascendant, Koch co-Ascendant, Munkasey
///     co-Ascendant and polar Ascendant, in that order (not Descendant/IC).
///
/// Raises
/// ------
/// TypeError
///     A string/bytes selector does not contain exactly one character/byte.
/// CoordinateError, PolarCircleError
///     Coordinates are invalid or the selected construction is polar-undefined.
///     Configuration, clock, required data and other construction failures
///     propagate. Sunshine alone may use analytic Sun declination on a range
///     miss or zero declination on selected value-domain failures; other data
///     and configuration errors are not a general fallback.
#[pyfunction(signature = (tjdut, lat, lon, hsys = 'P', iflag = 0), text_signature = "(tjdut, lat, lon, hsys=80, iflag=0)")]
fn houses(
    py: Python<'_>,
    tjdut: f64,
    lat: f64,
    lon: f64,
    #[pyo3(from_py_with = selector)] hsys: char,
    iflag: i32,
) -> PyResult<Layout> {
    let value = py
        .detach(|| {
            runtime()
                .house_session(hsys, iflag)?
                .public_houses(tjdut, lat, lon, hsys, iflag)
        })
        .map_err(|error| to_pyerr(py, error))?;
    layout(py, value)
}

/// Calculate tropical or configured sidereal houses at UT1 Julian day tjdut.
///
/// lat/lon are north/east geographic degrees with the domains of houses.
/// hsys defaults to 'P'/80 and accepts the same int, one-character str or
/// one-byte Latin-1 selector, including distinct 'i' and 'g' behavior and
/// unknown-selector Placidus geometry. flags defaults to 0; FLG_SIDEREAL
/// selects the configured sidereal plane/zero point, and FLG_RADIANS converts
/// central cusp and angle values to radians. Ephemeris bits affect Sunshine.
///
/// Returns (cusps, angles): 12 cusps, or 36 for uppercase 'G' (lowercase 'g'
/// gives the first 12 sectors), with index 0 = house/sector 1 and no padding.
/// Eight angles follow houses: Asc, MC, ARMC, Vertex, equatorial Ascendant,
/// Koch co-Ascendant, Munkasey co-Ascendant, polar Ascendant. Default units
/// are degrees. Sidereal branches retain their own frame/ARMC conventions;
/// this is not simply a scalar shift of every tropical angle.
///
/// Invalid selectors/coordinates, nonfinite epochs, polar-undefined systems
/// and consumed clock, data, anchor or construction failures raise exceptions.
/// This entry does not automatically retry another house system; use
/// houses_with_fallback for that separate tropical contract.
#[pyfunction(signature = (tjdut, lat, lon, hsys = 'P', flags = 0), text_signature = "(tjdut, lat, lon, hsys=80, flags=0)")]
fn houses_ex(
    py: Python<'_>,
    tjdut: f64,
    lat: f64,
    lon: f64,
    #[pyo3(from_py_with = selector)] hsys: char,
    flags: i32,
) -> PyResult<Layout> {
    let value = py
        .detach(|| {
            runtime()
                .house_session(hsys, flags)?
                .public_houses_ex(tjdut, lat, lon, hsys, flags)
        })
        .map_err(|error| to_pyerr(py, error))?;
    layout(py, value)
}

/// Calculate houses_ex positions and their rates at UT1 Julian day tjdut.
///
/// lat/lon are north/east degrees; hsys defaults to 'P'/80 and flags to 0.
/// Selectors, sidereal settings, coordinate domains and errors follow
/// houses_ex, including 'G' = 36 sectors and 'g' = the first 12 sectors.
///
/// Returns (cusps, angles, cusp_speeds, angle_speeds), all tuples. Cusps have
/// no leading zero (index 0 = house 1); angles use the eight-element houses
/// order. FLG_RADIANS changes central cusps/angles only: both speed tuples
/// ALWAYS remain degrees per UT1 day. A speed flag is not required.
///
/// Rates use centered time samples, with Richardson refinement for angles.
/// Samples retain the selected time/frame/ayanamsha models; any failed sample
/// fails the call. Discontinuities over 90 degrees and sign-pinned Whole Sign/
/// Aries cusp jumps use the established zero-rate convention, not a claim of
/// differentiability at a jump. No polar fallback is automatic.
#[pyfunction(signature = (tjdut, lat, lon, hsys = 'P', flags = 0), text_signature = "(tjdut, lat, lon, hsys=80, flags=0)")]
fn houses_ex2(
    py: Python<'_>,
    tjdut: f64,
    lat: f64,
    lon: f64,
    #[pyo3(from_py_with = selector)] hsys: char,
    flags: i32,
) -> PyResult<Rates> {
    let value = py
        .detach(|| {
            runtime()
                .house_session(hsys, flags)?
                .public_houses_ex2(tjdut, lat, lon, hsys, flags)
        })
        .map_err(|error| to_pyerr(py, error))?;
    rates(py, value)
}

/// Calculate tropical houses from ARMC and obliquity without a clock/data read.
///
/// armc is right ascension of the MC in degrees (reduced modulo 360); lat is
/// finite geographic north latitude [-90, 90] degrees; eps is true obliquity
/// in degrees. hsys defaults to 'P'/80 with the selector rules of houses.
/// ascmc9 is Sun declination in degrees for 'I'/'i' only (default 0), not an
/// extra returned angle; other systems ignore it.
///
/// Returns (cusps, angles) in degrees: 12 cusps with index 0 = house 1, or 36
/// for 'G'; 'g' returns their first 12. The eight angle meanings/order are the
/// same as houses. Invalid latitude/selector and polar or geometric failures
/// raise exceptions. Other floating-point inputs are not uniformly guarded;
/// nonphysical or nonfinite geometry can propagate IEEE values. No fallback
/// or sidereal/radian formatting is requested by this API.
#[pyfunction(signature = (armc, lat, eps, hsys = 'P', ascmc9 = 0.0), text_signature = "(armc, lat, eps, hsys=80, ascmc9=0.0)")]
fn houses_armc(
    py: Python<'_>,
    armc: f64,
    lat: f64,
    eps: f64,
    #[pyo3(from_py_with = selector)] hsys: char,
    ascmc9: f64,
) -> PyResult<Layout> {
    let value = py
        .detach(|| core_houses::houses_armc(armc, lat, eps, hsys, ascmc9))
        .map_err(|error| to_pyerr(py, error))?;
    layout(py, value)
}

/// Return ARMC-based cusps, angles and approximate daily angular rates.
///
/// Arguments/defaults and errors follow houses_armc: all angles are degrees;
/// ascmc9 (default 0) is Sun declination only for Sunshine, and hsys defaults
/// to 'P'/80. Returns (cusps, angles, cusp_speeds, angle_speeds), each a tuple,
/// with 12/36 cusp entries and eight angles in the houses order, no padding.
/// 'G' returns 36 sectors; lowercase 'g' returns their first 12.
///
/// Positions are degrees; rates are degrees per mean solar day from ARMC
/// samples separated by one sidereal second. Obliquity and supplied solar
/// declination are held fixed, so these are not full time derivatives.
/// Prefer houses_ex2 when a UT1 epoch is available. Discontinuous/sign-pinned
/// cusp branches retain zero-rate conventions. No files or clocks are read;
/// any central or side-sample construction error fails the call.
#[pyfunction(signature = (armc, lat, eps, hsys = 'P', ascmc9 = 0.0), text_signature = "(armc, lat, eps, hsys=80, ascmc9=0.0)")]
fn houses_armc_ex2(
    py: Python<'_>,
    armc: f64,
    lat: f64,
    eps: f64,
    #[pyo3(from_py_with = selector)] hsys: char,
    ascmc9: f64,
) -> PyResult<Rates> {
    let value = py
        .detach(|| core_houses::houses_armc_ex2(armc, lat, eps, hsys, ascmc9))
        .map_err(|error| to_pyerr(py, error))?;
    rates(py, value)
}

/// Calculate tropical houses with an explicit polar/invalid-cusp fallback.
///
/// jd_ut is UT1 JD; lat/lon are north/east degrees with the houses domains.
/// hsys is required; fallback_hsys defaults to 'O'/79 (Porphyry). Both accept
/// houses selectors, including the distinct G/g layouts. validate_cusps
/// defaults to True and checks primary cusps for finite values in [0, 360).
///
/// Returns (cusps, angles, used_fallback, warning), with degree tuples in the
/// houses layout, a bool and str or None. Index 0 is house/sector 1; the
/// actually served system determines whether there are 12 or 36 cusps.
/// A polar refusal or, when enabled, invalid primary cusps triggers a single
/// fallback calculation. At abs(lat) >= 80 an instability advisory can instead
/// accompany successful primary cusps with used_fallback=False. The JD-based
/// advisory uses nominal obliquity 23.44 degrees. Warning text is returned,
/// not emitted through Python warnings.
///
/// Coordinate, clock, data and non-polar construction errors propagate;
/// fallback failures also propagate, without recursive retries. Disabling
/// validation skips only the cusp check, not coordinate/polar admission.
#[pyfunction(signature = (jd_ut, lat, lon, hsys, fallback_hsys = 'O', validate_cusps = true), text_signature = "(jd_ut, lat, lon, hsys, fallback_hsys=79, validate_cusps=True)")]
fn houses_with_fallback(
    py: Python<'_>,
    jd_ut: f64,
    lat: f64,
    lon: f64,
    #[pyo3(from_py_with = selector)] hsys: char,
    #[pyo3(from_py_with = selector)] fallback_hsys: char,
    validate_cusps: bool,
) -> PyResult<Fallback> {
    let (cusps, angles, used, warning) = py
        .detach(|| {
            runtime().shared_session()?.public_houses_with_fallback(
                jd_ut,
                lat,
                lon,
                hsys,
                fallback_hsys,
                validate_cusps,
            )
        })
        .map_err(|error| to_pyerr(py, error))?;
    let (cusps, angles) = layout(py, (cusps, angles))?;
    Ok((cusps, angles, used, warning))
}

/// Calculate ARMC-based tropical houses with one explicit fallback attempt.
///
/// armc, lat, eps and Sunshine-only ascmc9 are degrees with the houses_armc
/// meanings/domains. hsys is required; fallback_hsys defaults to 'O'/79,
/// validate_cusps to True, ascmc9 to 0. Both selectors follow houses rules.
///
/// Returns (cusps, angles, used_fallback, warning), in degrees, with no cusp
/// padding and the eight-angle houses order. Uppercase 'G' supplies 36 sectors;
/// lowercase 'g' only 12. Polar primary failures or validated invalid cusps
/// trigger the fallback once; extreme-latitude advisories may return the
/// primary unchanged with used_fallback=False. Warning is str or None, not
/// an emitted Python warning. The supplied eps drives advisory classification.
/// All other failures, and any fallback failure, propagate. No ephemeris or
/// time-scale data are read, and this is not a sidereal/radian entry.
#[pyfunction(signature = (armc, lat, eps, hsys, fallback_hsys = 'O', validate_cusps = true, ascmc9 = 0.0), text_signature = "(armc, lat, eps, hsys, fallback_hsys=79, validate_cusps=True, ascmc9=0.0)")]
#[allow(clippy::too_many_arguments)]
fn houses_armc_with_fallback(
    py: Python<'_>,
    armc: f64,
    lat: f64,
    eps: f64,
    #[pyo3(from_py_with = selector)] hsys: char,
    #[pyo3(from_py_with = selector)] fallback_hsys: char,
    validate_cusps: bool,
    ascmc9: f64,
) -> PyResult<Fallback> {
    let value = py
        .detach(|| {
            core_houses::houses_armc_with_fallback(
                armc,
                lat,
                eps,
                hsys,
                fallback_hsys,
                validate_cusps,
                ascmc9,
            )
        })
        .map_err(|error| to_pyerr(py, error))?;
    let (cusps, angles) = layout(py, (value.cusps, value.ascmc))?;
    Ok((cusps, angles, value.used_fallback, value.warning))
}

/// Return the display name for an int, one-character str or one-byte selector.
///
/// ASCII case folds as in houses, except that 'i' retains its distinct name;
/// 'g' is named Gauquelin even though its cusp layout has only 12 entries.
/// Unknown selectors return the empty string. Invalid str/bytes lengths raise
/// TypeError; integer extraction/overflow errors propagate.
#[pyfunction]
fn house_name(#[pyo3(from_py_with = selector)] hsys: char) -> &'static str {
    core_houses::house_name(hsys)
}

/// Locate an ecliptic point within a house, using supplied geometry only.
///
/// Parameters
/// ----------
/// armc, lat, obliquity : float
///     ARMC, geographic north latitude and true obliquity, all degrees.
/// hsys_or_objcoord, lon_or_hsys
///     Use system-first ``house_pos(armc, lat, eps, hsys, longitude,
///     lat_body=0.0)`` or coordinate-first ``house_pos(armc, lat, eps,
///     (longitude, latitude), hsys)``. The latter selector defaults to 'P'
///     when omitted. A system-first call also accepts a two-component
///     coordinate sequence. With a system but no coordinate, longitude is 0.
///     A pair supplies both angles and ignores lat_body. Selectors follow
///     houses; the actual keyword names are those in the signature.
///
/// Returns
/// -------
/// float
///     Decimal house in [1, 13), or sector in [1, 37) for either 'G' or 'g'.
///     The integer part is the one-based house; Koch returns 0.0 when its
///     placement is undefined. No time conversion, refraction or observer
///     mutation occurs. Sunshine placement uses its context-free geometry,
///     not an acquired solar declination.
///
/// Raises
/// ------
/// ValueError
///     A consumed angle is nonfinite. Malformed coordinate/selector inputs
///     raise extraction/TypeError exceptions; cusp-based branches can propagate
///     coordinate, polar or geometric construction errors. Body coordinates
///     must refer to the ecliptic consistent with obliquity and ARMC.
///
/// Examples
/// --------
/// >>> import libephemeris as ephe
/// >>> a = ephe.house_pos(120.0, 41.9, 23.44, 'P', 150.0, 2.0)
/// >>> b = ephe.house_pos(120.0, 41.9, 23.44, (150.0, 2.0), 'P')
/// >>> a == b
/// True
#[pyfunction(signature = (armc, lat, obliquity, hsys_or_objcoord, lon_or_hsys = None, lat_body = 0.0))]
fn house_pos(
    py: Python<'_>,
    armc: f64,
    lat: f64,
    obliquity: f64,
    hsys_or_objcoord: &Bound<'_, PyAny>,
    lon_or_hsys: Option<&Bound<'_, PyAny>>,
    lat_body: f64,
) -> PyResult<f64> {
    let system_input = hsys_or_objcoord.is_instance_of::<PyInt>()
        || hsys_or_objcoord.is_instance_of::<PyString>()
        || hsys_or_objcoord.is_instance_of::<PyBytes>();
    let (system, longitude, latitude) = if system_input {
        let system = selector(hsys_or_objcoord)?;
        if let Some(value) = lon_or_hsys {
            if let Ok([longitude, latitude]) = value.extract::<[f64; 2]>() {
                (system, longitude, latitude)
            } else {
                (system, value.extract::<f64>()?, lat_body)
            }
        } else {
            (system, 0.0, lat_body)
        }
    } else {
        let [longitude, latitude] = hsys_or_objcoord.extract::<[f64; 2]>()?;
        let system = lon_or_hsys.map(selector).transpose()?.unwrap_or('P');
        (system, longitude, latitude)
    };
    py.detach(|| core_houses::house_pos(armc, lat, obliquity, system, longitude, latitude))
        .map_err(|error| to_pyerr(py, error))
}

/// Return the decimal Gauquelin sector [1, 37) for a body or fixed star.
///
/// tjdut is UT1 JD; body is an integer body ID or star-name str. geopos is
/// (east longitude degrees, north latitude degrees, height metres). method 0
/// uses the ecliptic place, 1 forces its latitude to zero; 2/3 use rise/set
/// with disc center and 4/5 use limb events. Methods 2/4 omit refraction,
/// 3/5 include it. Other method integers follow method 0. atpress is mbar
/// and attemp Celsius; zero defaults become 1013.25 mbar and 15 C for events.
/// flags defaults to 32770; methods 0/1 pass the requested flags to the
/// position entry with SPEED added, without undoing XYZ/radian formatting.
/// Supply ordinary ecliptic-degree flags for meaningful sector geometry.
///
/// The returned sector is a scalar, not a 36-cusp tuple. Event methods use
/// actual neighboring rise/set events and fail when they cannot be found,
/// rather than falling back to method 0. Catalog, coordinate, coverage,
/// configuration, clock and calculation exceptions propagate. Calls are
/// synchronous and retain any per-call observer without changing caller state.
#[pyfunction(signature = (tjdut, body, method, geopos, atpress = 0.0, attemp = 0.0, flags = 32770))]
#[allow(clippy::too_many_arguments)]
fn gauquelin_sector(
    py: Python<'_>,
    tjdut: f64,
    body: &Bound<'_, PyAny>,
    method: i32,
    geopos: [f64; 3],
    atpress: f64,
    attemp: f64,
    flags: i32,
) -> PyResult<f64> {
    if let Ok(star) = body.cast::<PyString>() {
        let star = star.to_str()?.to_owned();
        py.detach(|| {
            runtime()
                .shared_session()?
                .star_gauquelin_sector(&star, tjdut, method, geopos, atpress, attemp, flags)
        })
        .map_err(|error| to_pyerr(py, error))
    } else {
        let body = body.extract::<i32>()?;
        py.detach(|| {
            runtime()
                .shared_session()?
                .gauquelin_sector(tjdut, body, method, geopos, atpress, attemp, flags)
        })
        .map_err(|error| to_pyerr(py, error))
    }
}

/// Return ``a - b`` wrapped to (-180, 180] degrees (half-turn = +180).
///
/// Inputs are degree angles; this helper performs no validation/acquisition.
/// Nonfinite arithmetic may return NaN rather than an exception.
#[pyfunction]
fn angular_diff(py: Python<'_>, a: f64, b: f64) -> f64 {
    py.detach(|| core_houses::angular_diff(a, b))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let helpers = PyModule::new(module.py(), "libephemeris._core._houses")?;
    helpers.add_function(wrap_pyfunction!(angular_diff, &helpers)?)?;
    module
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item("libephemeris._core._houses", &helpers)?;
    module.add("_houses", helpers)?;
    macro_rules! add {
        ($($name:ident),+ $(,)?) => { $(module.add_function(wrap_pyfunction!($name, module)?)?;)+ };
    }
    add!(
        houses,
        houses_ex,
        houses_ex2,
        houses_armc,
        houses_armc_ex2,
        houses_with_fallback,
        houses_armc_with_fallback,
        house_name,
        house_pos,
        gauquelin_sector,
        get_polar_latitude_threshold,
        get_extreme_latitude_info
    );
    Ok(())
}
