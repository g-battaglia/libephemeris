//! Thin conversions for the public numeric frame and scalar-cache modules.

use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::python_support::numeric_frames as native;
use libephemeris_core::{erfa, lunar, orbital, positions, python_support as api, sidereal};
use numpy::PyArray2;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};

type Triple = (f64, f64, f64);
type Six = (f64, f64, f64, f64, f64, f64);
type MatrixTuple = (Triple, Triple, Triple);
type FlatMatrixTuple = (f64, f64, f64, f64, f64, f64, f64, f64, f64);
fn triple([x, y, z]: [f64; 3]) -> Triple {
    (x, y, z)
}
fn six([a, b, c, d, e, f]: [f64; 6]) -> Six {
    (a, b, c, d, e, f)
}
fn matrix(value: native::Matrix3) -> MatrixTuple {
    (triple(value[0]), triple(value[1]), triple(value[2]))
}
fn call<T: Send>(
    py: Python<'_>,
    operation: impl FnOnce() -> libephemeris_core::errors::Result<T> + Send,
) -> PyResult<T> {
    py.detach(operation).map_err(|error| to_pyerr(py, error))
}
fn session<T: Send>(
    py: Python<'_>,
    operation: impl FnOnce(&api::Session) -> libephemeris_core::errors::Result<T> + Send,
) -> PyResult<T> {
    call(py, || operation(runtime().shared_session()?.as_ref()))
}
fn array(py: Python<'_>, value: native::Matrix3) -> PyResult<Bound<'_, PyArray2<f64>>> {
    PyArray2::from_vec2(py, &value.map(|row| row.to_vec())).map_err(Into::into)
}
fn info(py: Python<'_>, value: native::CacheInfo) -> PyResult<Bound<'_, PyDict>> {
    let out = PyDict::new(py);
    out.set_item("hits", value.hits)?;
    out.set_item("misses", value.misses)?;
    out.set_item("maxsize", value.maxsize)?;
    out.set_item("currsize", value.currsize)?;
    Ok(out)
}

/// Return IAU 2000B nutation (dpsi, deps) in DEGREES at TT Julian day jd_tt.
///
/// This is the abbreviated 77-term model, not the radian IAU 2006/2000A
/// cache interface. Nonfinite TT raises ValueError. The series is intended
/// for dates within a few centuries of J2000, not arbitrary-epoch accuracy.
#[pyfunction(name="nutation_angles", signature=(jd_tt))]
fn astrometry_nutation_angles(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64)> {
    call(py, || native::nutation_angles(jd_tt))
}
/// Precess J2000 mean-ecliptic longitude/latitude to the mean ecliptic of date.
///
/// lon and lat are degrees; jd_tt is a TT Julian day. Returns (lon, lat) in
/// degrees using the native long-term precession frame, without nutation or
/// apparent-place corrections. Nonfinite inputs and frame evaluation failures
/// raise exceptions; finite extrapolation does not guarantee model accuracy.
#[pyfunction(signature=(lon,lat,jd_tt))]
fn precess_from_j2000(py: Python<'_>, lon: f64, lat: f64, jd_tt: f64) -> PyResult<(f64, f64)> {
    call(py, || native::precess_from_j2000(lon, lat, jd_tt))
}
/// Apply normalized first-order annual aberration to an ICRS Cartesian ray.
///
/// position is a three-component observer-relative vector in AU;
/// earth_velocity is the observer's barycentric ICRS velocity in AU/day.
/// Returns a three-float tuple in AU, preserving radius. This helper does not
/// acquire a velocity or solve light time. A zero ray is returned unchanged.
/// Nonfinite input components raise ValueError; arbitrarily large finite
/// inputs can overflow derived arithmetic and produce IEEE NaN/infinity.
#[pyfunction(signature=(position,earth_velocity))]
fn apply_aberration_to_position(
    py: Python<'_>,
    position: [f64; 3],
    earth_velocity: [f64; 3],
) -> PyResult<Triple> {
    call(py, || {
        native::apply_aberration_to_position(position, earth_velocity)
    })
    .map(triple)
}
/// Return cached IAU 2006/2000A (dpsi, deps) in RADIANS at TT JD jd_tt.
///
/// Uses a process-wide bounded 256-pair cache, keyed by the exact epoch word.
/// Nonfinite TT raises ValueError. No observed free-core nutation is included;
/// the IAU series is intended for epochs within a few centuries of J2000.
#[pyfunction(signature=(jd_tt))]
fn get_cached_nutation(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64)> {
    call(py, || native::cached_nutation(jd_tt))
}
/// Return IAU 2006/2000A (dpsi, deps) in DEGREES at TT Julian day jd_tt.
///
/// Uses the same bounded cache and model as get_cached_nutation, converting
/// its radian pair. Nonfinite TT raises ValueError; no observed nutation feed
/// or planetary ephemeris is read.
#[pyfunction(signature=(jd_tt))]
fn get_nutation_degrees(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64)> {
    call(py, || native::nutation_angles_iau2006_degrees(jd_tt))
}
/// Return cached (mean_obliquity, true_obliquity) in degrees at TT JD jd_tt.
///
/// Mean obliquity is the angle between Vondrak mean poles; true adds the
/// IAU 2006/2000A obliquity nutation. The 256-entry cache is process-wide.
/// Nonfinite TT raises ValueError; long-term precession validity does not
/// extend the near-J2000 accuracy domain of the nutation series.
#[pyfunction(signature=(jd_tt))]
fn get_cached_obliquity(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64)> {
    call(py, || native::cached_obliquity(jd_tt))
}
/// Return Vondrak pole-angle mean obliquity in DEGREES at TT JD jd_tt.
///
/// Uses the mean/true obliquity cache, which also evaluates nutation.
/// Nonfinite TT raises ValueError; no planetary data is required.
#[pyfunction(signature=(jd_tt))]
fn get_mean_obliquity(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    call(py, || native::cached_obliquity(jd_tt).map(|value| value.0))
}
/// Return true obliquity in DEGREES at TT Julian day jd_tt.
///
/// This is Vondrak pole-angle mean obliquity plus IAU 2006/2000A deps, from
/// the bounded shared cache. Nonfinite TT raises ValueError. No measured EOP
/// correction is acquired; the nutation series has its own near-J2000 domain.
#[pyfunction(signature=(jd_tt))]
fn get_true_obliquity(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    call(py, || native::cached_obliquity(jd_tt).map(|value| value.1))
}
/// Convert UT1 Julian day jd to a cached TWO-PART TT Julian day.
///
/// Returns (part1, part2), whose sum is TT, not (TT, UT1). Uses the independent
/// bundled frame clock, not the model Delta-T user override. Invalid epochs
/// or unavailable frame-clock coverage raise native clock/range exceptions.
/// Retains at most 256 scalar pairs; no Python time-scale object is created.
#[pyfunction(signature=(jd))]
fn get_cached_time_ut1(py: Python<'_>, jd: f64) -> PyResult<(f64, f64)> {
    call(py, || native::cached_time_ut1(runtime(), jd))
}
/// Split finite TT Julian day jd into cached (floor(jd), jd - floor(jd)).
///
/// Both components are TT days: this is not a TT/UT1 pair or a clock
/// conversion. Nonfinite input raises ValueError. The shared cache retains
/// at most 256 pairs and performs no ephemeris or observation-file access.
#[pyfunction(signature=(jd))]
fn get_cached_time_tt(py: Python<'_>, jd: f64) -> PyResult<(f64, f64)> {
    call(py, || native::cached_time_tt(jd))
}
/// Clear numeric nutation, obliquity, TT/UT1-split and ERFA model caches.
///
/// Returns None and resets their hit/miss counters process-wide. Active calls
/// retain already copied values. Does not close data files, reset observing
/// configuration or promise reclamation of every engine/OS cache.
#[pyfunction]
fn clear_caches(py: Python<'_>) {
    py.detach(native::clear_caches);
}
/// Return counters for the four numeric caches, without acquiring data.
///
/// The dict keys are 'nutation', 'obliquity', 'time_ut1' and 'time_tt'. Each
/// value is a dict with integer 'hits', 'misses', 'maxsize' (256) and
/// 'currsize'. The separate model-selectable ERFA cache is not included.
/// Counts are records, not a bound on total process memory.
#[pyfunction]
fn get_cache_info(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    let values = py.detach(native::cache_info);
    let out = PyDict::new(py);
    for (name, value) in values {
        out.set_item(name, info(py, value)?)?;
    }
    Ok(out)
}
/// Return True: native ERFA-derived formulas are compiled into this engine.
///
/// This does not import Python ERFA, probe a planetary dataset or certify
/// model accuracy at a requested date. Takes no arguments or external data.
#[pyfunction]
fn has_erfa() -> bool {
    true
}
/// Return IAU 2000A (dpsi, deps) in RADIANS at TT Julian day jd_tt.
///
/// Evaluates the full native luni-solar/planetary series without observed
/// free-core nutation. Nonfinite TT raises ValueError. Model accuracy applies
/// within a few centuries of J2000; this does not read an ephemeris file.
#[pyfunction(signature=(jd_tt))]
fn get_erfa_nutation_nut00a(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64)> {
    call(py, || erfa::nut00a(jd_tt))
}
/// Return IAU 2006/2000A (dpsi, deps) in RADIANS at TT Julian day jd_tt.
///
/// Applies the P03 adjustments to IAU 2000A; no observed free-core nutation
/// is included. Nonfinite TT raises ValueError. Accuracy is that of the IAU
/// series near J2000, not arbitrary-date certification.
#[pyfunction(signature=(jd_tt))]
fn get_erfa_nutation_nut06a(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64)> {
    call(py, || erfa::nut06a(jd_tt))
}
/// Return the IAU 2006/2000A bias-precession-nutation matrix at TT JD jd_tt.
///
/// Returns a dimensionless NumPy float64 array of shape (3, 3), row major,
/// mapping GCRS column vectors to the true equator/equinox of date. No Earth
/// rotation or polar motion is included. Nonfinite TT and frame arithmetic
/// failures raise exceptions; the IAU model's near-J2000 accuracy applies.
#[pyfunction(signature=(jd_tt))]
fn get_erfa_pnm06a_matrix(py: Python<'_>, jd_tt: f64) -> PyResult<Bound<'_, PyArray2<f64>>> {
    array(py, call(py, || erfa::pnm06a(jd_tt))?)
}
/// Return cached (dpsi, deps) in radians for an explicitly selected IAU model.
///
/// jd_tt is TT JD. model is exactly 'nut00a' or 'nut06a' (default 'nut06a'),
/// case-sensitive; other names and nonfinite epochs raise ValueError. The
/// process-wide cache holds at most 256 (epoch, model) records. No external
/// Python ERFA package or observed free-core-nutation data is used.
#[pyfunction(signature=(jd_tt,model="nut06a"))]
fn get_erfa_nutation_cached(py: Python<'_>, jd_tt: f64, model: &str) -> PyResult<(f64, f64)> {
    call(py, || native::erfa_nutation_cached(jd_tt, model))
}
/// Return IAU 2006 polynomial mean obliquity in RADIANS at TT JD jd_tt.
///
/// This differs from the Vondrak pole-angle obliquity helpers. Nonfinite TT
/// raises ValueError. Finite extrapolation is evaluated without a date-window
/// guard, but is not a guarantee of distant-epoch physical accuracy.
#[pyfunction(signature=(jd_tt))]
fn get_erfa_obliquity_iau2006(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    call(py, || native::erfa_obliquity_iau2006(jd_tt))
}
/// Return native nutation model diagnostics at TT Julian day jd_tt.
///
/// The dict contains 'jd_tt' and model dicts 'skyfield_iau2000b',
/// 'skyfield_iau2000a', 'erfa_nut00a', 'erfa_nut06a', each with
/// 'dpsi_arcsec'/'deps_arcsec'. Historical backend labels name mathematical
/// models: every value is native. 'differences_mas' contains dpsi/deps keys
/// prefixed 'skyfield_2000a_vs_erfa_00a_', 'erfa_00a_vs_erfa_06a_' and
/// 'skyfield_2000b_vs_2000a_', in milliarcseconds, left model minus right.
/// The first pair is zero by construction. Nonfinite TT raises ValueError;
/// agreement between these series is not independent accuracy evidence.
#[pyfunction(signature=(jd_tt))]
fn compare_nutation_models(py: Python<'_>, jd_tt: f64) -> PyResult<Bound<'_, PyDict>> {
    let values = call(py, || native::compare_nutation_models(jd_tt))?;
    let out = PyDict::new(py);
    out.set_item("jd_tt", jd_tt)?;
    for (key, value) in values {
        out.set_item(key, value)?;
    }
    Ok(out)
}
/// Clear the process-wide model-selectable ERFA nutation cache and counters.
///
/// Returns None. Other numeric caches, contexts and ephemeris owners remain
/// unchanged; active calculations retain their copied results.
#[pyfunction]
fn clear_erfa_nutation_cache(py: Python<'_>) {
    py.detach(native::clear_erfa_nutation_cache);
}
/// Return ERFA nutation-cache integer counters as a dict.
///
/// Keys are 'hits', 'misses', 'maxsize' (256 records) and 'currsize'. No
/// astronomical data is acquired; these are not process memory measurements.
#[pyfunction]
fn get_erfa_nutation_cache_info(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    info(py, py.detach(native::erfa_nutation_cache_info))
}

/// Return IERS 2003 lunar arguments (D, M, M_prime, F) in RADIANS at TT JD.
///
/// These are mean lunar elongation, solar anomaly, lunar anomaly and lunar
/// argument of latitude, in that order. No coverage/finite-input check is
/// performed; nonfinite/overflowing polynomial arithmetic may yield NaNs.
/// No provider, apparent correction or observation clock is acquired.
#[pyfunction(signature=(jd_tt))]
fn lunar_delaunay_arguments(py: Python<'_>, jd_tt: f64) -> (f64, f64, f64, f64) {
    py.detach(|| lunar::lunar_delaunay_arguments(jd_tt))
}
/// Return (mean_apogee_plane_longitude, mean_node_longitude) in degrees.
///
/// jd_tt is TT JD. These are the analytic mean ecliptic-of-date arguments,
/// not an apparent position or a pair of nutation angles. The configured
/// runtime supplies the main-source coverage policy, potentially acquiring
/// metadata. Nonfinite/uncovered epochs raise EphemerisRangeError; required
/// configuration or resource errors propagate. No extrapolation is promised.
#[pyfunction(signature=(jd_tt))]
fn mean_lunar_apse_plane_node(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64)> {
    session(py, |s| s.numeric_mean_lunar_plane_node(jd_tt))
}
/// Return mean lunar apogee (longitude, latitude, distance) at TT JD jd_tt.
///
/// Longitude/latitude refer to the model's mean ecliptic of date in degrees;
/// distance is the conventional constant mean-apogee radius in AU. No apparent
/// corrections or user output flags are applied. Nonfinite/uncovered dates
/// raise EphemerisRangeError under configured main-source metadata policy;
/// metadata/configuration/resource errors can also propagate.
#[pyfunction(signature=(jd_tt))]
fn mean_lunar_apogee_position(py: Python<'_>, jd_tt: f64) -> PyResult<Triple> {
    session(py, |s| s.numeric_mean_lunar_apogee_position(jd_tt))
}
/// Return analytic mean-apogee spherical position and rates at TT JD jd_tt.
///
/// The six-float tuple is (lon, lat, radius, lon_rate, lat_rate, radius_rate)
/// in mean ecliptic-of-date degrees, AU, degrees/TT-day and AU/TT-day.
/// Radius is conventional and constant; radius_rate is zero. Keyword-only
/// speed3 defaults to False and is accepted but ignored by this model.
/// Nonfinite/uncovered epochs raise EphemerisRangeError under the configured
/// main-source coverage policy; required metadata/resource errors propagate.
#[pyfunction(signature=(jd_tt,*,speed3=false))]
fn mean_lunar_apogee_state(py: Python<'_>, jd_tt: f64, speed3: bool) -> PyResult<Six> {
    session(py, |s| s.numeric_mean_lunar_apogee_state(jd_tt, speed3))
}
/// Return mean lunar node (longitude, 0, conventional_distance_au) at TT JD.
///
/// Longitude is degrees on the model's mean ecliptic of date; latitude is
/// zero and distance is the constant mean-node radius, not a fitted range.
/// Nonfinite/uncovered jd_tt raises EphemerisRangeError under configured
/// main-source metadata policy; required configuration/data errors propagate.
#[pyfunction(signature=(jd_tt))]
fn mean_lunar_node_position(py: Python<'_>, jd_tt: f64) -> PyResult<Triple> {
    session(py, |s| s.numeric_mean_lunar_node_position(jd_tt))
}
/// Return mean-node (lon, lat, radius, lon_rate, lat_rate, radius_rate) at TT.
///
/// jd_tt is TT JD; output is mean ecliptic-of-date degrees, AU and rates per
/// TT day. Latitude, latitude rate and radial rate are zero; radius is the
/// conventional mean-node distance. Keyword-only speed3 defaults to False
/// and is ignored. Nonfinite/uncovered epochs raise EphemerisRangeError under
/// configured main-source policy; metadata/configuration errors propagate.
#[pyfunction(signature=(jd_tt,*,speed3=false))]
fn mean_lunar_node_state(py: Python<'_>, jd_tt: f64, speed3: bool) -> PyResult<Six> {
    session(py, |s| s.numeric_mean_lunar_node_state(jd_tt, speed3))
}
/// Return mean elements (a, e, inclination, ascending_node, perihelion_longitude).
///
/// ipl is a supported Mercury-through-Neptune or Earth body ID; jd_tt is TT
/// JD. a is AU, e dimensionless, and angles are degrees referred to the mean
/// ecliptic/equinox of date. Earth has zero inclination/node; node and
/// perihelion longitude are wrapped to [0, 360). Unsupported bodies (including
/// Pluto) or nonfinite epochs return None, not an exception. The Meeus/Simon
/// polynomial is an approximation; no ephemeris or clock is acquired.
#[pyfunction(signature=(ipl,jd_tt))]
fn mean_orbital_elements(
    py: Python<'_>,
    ipl: i32,
    jd_tt: f64,
) -> Option<(f64, f64, f64, f64, f64)> {
    py.detach(|| orbital::mean_orbital_elements(ipl, jd_tt).map(|[a, b, c, d, e]| (a, b, c, d, e)))
}
/// Return mean obliquity in radians from Vondrak mean poles at TT JD jd_tt.
///
/// Nonfinite epochs raise ValueError. This pole-angle realization, not the
/// separate direct series, keeps the equator/ecliptic frames consistent.
/// The long-term model spans about +/-200,000 years around J2000 with growing
/// physical uncertainty; no accuracy certification or date clamp is implied.
#[pyfunction(signature=(jd_tt))]
fn vondrak_mean_obliquity_rad(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    call(py, || positions::vondrak_mean_obliquity_rad(jd_tt))
}
/// Return Vondrak pole-angle mean obliquity in degrees at TT JD jd_tt.
///
/// Degree form of vondrak_mean_obliquity_rad with the same long-term model
/// limitations. Nonfinite epochs raise ValueError; no data file is acquired.
#[pyfunction(signature=(jd_tt))]
fn vondrak_mean_obliquity_deg(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    call(py, || positions::vondrak_mean_obliquity_deg(jd_tt))
}
/// Return a Vondrak precession matrix to the mean equator/equinox of date.
///
/// jd_tt is TT JD. With frame_bias=True (default), input column vectors are
/// ICRS; with False, they are mean equatorial J2000. Output is a nested tuple
/// of three row tuples, NOT a NumPy array or flat tuple. Entries are
/// dimensionless; matrix multiplication preserves the vector's length unit.
/// No nutation, Earth rotation or polar motion is included. Nonfinite epochs,
/// overflowing or undefined frame geometry raise exceptions. The published
/// long-term model's growing uncertainty remains; no dataset is acquired.
#[pyfunction(signature=(jd_tt,frame_bias=true))]
fn vondrak_precession_matrix(
    py: Python<'_>,
    jd_tt: f64,
    frame_bias: bool,
) -> PyResult<MatrixTuple> {
    call(py, || {
        positions::vondrak_precession_matrix(jd_tt, frame_bias)
    })
    .map(matrix)
}
/// Return (precession_nutation_matrix, true_obliquity) for a TT epoch.
///
/// jd_tt is TT JD; dpsi/deps are caller-supplied longitude/obliquity nutation
/// in RADIANS. frame_bias=True (default) selects ICRS input; False selects
/// mean-equatorial J2000 input. The dimensionless nested row-tuple matrix
/// maps column vectors to the true equator/equinox of date; true_obliquity
/// is radians, using Vondrak pole-angle mean obliquity plus deps. No nutation
/// feed is fetched. Nonfinite inputs and frame arithmetic failures raise
/// exceptions; model validity is that of the selected precession/nutation.
#[pyfunction(signature=(jd_tt,dpsi,deps,frame_bias=true))]
fn vondrak_pn_matrix(
    py: Python<'_>,
    jd_tt: f64,
    dpsi: f64,
    deps: f64,
    frame_bias: bool,
) -> PyResult<(MatrixTuple, f64)> {
    call(py, || {
        positions::vondrak_pn_matrix(jd_tt, dpsi, deps, frame_bias)
    })
    .map(|(m, e)| (matrix(m), e))
}
/// Return Method-B precession accumulated from t0_tt to jd_tt, in degrees.
///
/// Both inputs are TT Julian days. The date equinox is projected onto the
/// fixed mean ecliptic of t0_tt; the result is minus its longitude. It is not
/// reduced modulo 360 and is not additive in the reference epoch. Nonfinite
/// arguments raise ValueError. Extreme finite extrapolation can yield NaN
/// and does not extend the Vondrak model's physical validity.
#[pyfunction(signature=(jd_tt,t0_tt))]
fn method_b_accumulated_precession(py: Python<'_>, jd_tt: f64, t0_tt: f64) -> PyResult<f64> {
    call(py, || {
        sidereal::method_b_accumulated_precession(jd_tt, t0_tt)
    })
}
/// Evaluate the direct Vondrak mean-obliquity series, in radians at TT JD.
///
/// jd_tt must be finite or ValueError is raised. This diagnostic series is
/// distinct from mean_obliquity_rad's angle between the mean frame poles.
/// No data or nutation is acquired; distant-epoch uncertainty remains.
#[pyfunction(signature=(jd_tt))]
fn mean_obliquity_series_rad(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    call(py, || positions::vondrak_mean_obliquity_series_rad(jd_tt))
}
/// Return Vondrak pole-angle mean obliquity in radians at TT JD jd_tt.
///
/// Same native realization as vondrak_mean_obliquity_rad, without the
/// direct-series diagnostic or nutation. Nonfinite input raises ValueError;
/// the long-term model is not an arbitrary-date accuracy guarantee.
#[pyfunction(signature=(jd_tt))]
fn mean_obliquity_rad(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    call(py, || positions::vondrak_mean_obliquity_rad(jd_tt))
}
/// Return Vondrak pole-angle mean obliquity in degrees at TT JD jd_tt.
///
/// Same realization as mean_obliquity_rad with degree output. Nonfinite input
/// raises ValueError; no ephemeris file, nutation feed or clock is acquired.
#[pyfunction(signature=(jd_tt))]
fn mean_obliquity_deg(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    call(py, || positions::vondrak_mean_obliquity_deg(jd_tt))
}
/// Return unbiased J2000-to-date mean-equatorial precession as a FLAT tuple.
///
/// jd_tt is TT JD. The nine floats are row-major (m00,m01,m02,m10,...,m22),
/// acting on mean-equatorial J2000 column vectors to produce mean equator/
/// equinox-of-date components. No ICRS frame bias or nutation is included.
/// Nonfinite epochs, overflowing or undefined Vondrak frame geometry raise
/// exceptions. Entries are dimensionless; no data file is acquired.
#[pyfunction(signature=(jd_tt))]
fn precession_matrix(py: Python<'_>, jd_tt: f64) -> PyResult<FlatMatrixTuple> {
    call(py, || positions::vondrak_precession_matrix_nobias(jd_tt)).map(|m| {
        (
            m[0][0], m[0][1], m[0][2], m[1][0], m[1][1], m[1][2], m[2][0], m[2][1], m[2][2],
        )
    })
}
/// Return Greenwich MEAN sidereal angle in [0, 360) degrees at UT1 JD.
///
/// jd_ut1 is not TT/UTC. Uses the configured session's Delta-T/sidereal
/// model, including its long-term continuity terms. Invalid epochs or
/// required clock/configuration failures propagate; no nutation is added.
#[pyfunction(signature=(jd_ut1))]
fn mean_sidereal_time_deg(py: Python<'_>, jd_ut1: f64) -> PyResult<f64> {
    session(py, |s| s.numeric_mean_sidereal_time_deg(jd_ut1))
}
/// Return LOCAL apparent sidereal angle in [0, 360) degrees at UT1 JD.
///
/// jd_ut1 is UT1 JD; longitude is east geographic degrees, default 0
/// (Greenwich). Required keyword-only dpsi_deg and eps_true_deg are nutation
/// in longitude and true obliquity in DEGREES, supplied by the caller. Adds
/// their equation of equinoxes to the configured mean sidereal model.
/// Consumed nonfinite values or clock/configuration failures raise exceptions;
/// no observed nutation or site is acquired implicitly.
#[pyfunction(signature=(jd_ut1,longitude=0.0,*,dpsi_deg,eps_true_deg))]
fn apparent_sidereal_time_deg(
    py: Python<'_>,
    jd_ut1: f64,
    longitude: f64,
    dpsi_deg: f64,
    eps_true_deg: f64,
) -> PyResult<f64> {
    session(py, |s| {
        s.numeric_apparent_sidereal_time_deg(jd_ut1, longitude, dpsi_deg, eps_true_deg)
    })
}
/// Test whether iflag requests SIDEREAL under a fixed-epoch base sid_mode.
///
/// sid_mode must be the already stripped base ID (J2000, J1900, B1950 or
/// Galactic Alignment/Mardyks); projection bits are not stripped here.
/// Returns bool without calculation, configuration lookup or validation.
#[pyfunction(signature=(iflag,sid_mode))]
fn is_fixed_epoch_request(iflag: i32, sid_mode: i32) -> bool {
    sidereal::is_fixed_epoch_request(iflag, sid_mode)
}
/// Rewrite flags for the internal mean-J2000 fixed-epoch calculation.
///
/// Returns an int with SIDEREAL and RADIANS cleared and J2000 and NONUT set;
/// other bits are retained. This is bit manipulation, not validation or a
/// calculation, and it does not change the runtime's sidereal mode.
#[pyfunction(signature=(iflag))]
fn fixed_epoch_request_flags(iflag: i32) -> i32 {
    sidereal::fixed_epoch_request_flags(iflag)
}
/// Restore public flags after a fixed-epoch sidereal sub-calculation.
///
/// retflag is the sub-request result; iflag is the original request. Returns
/// retflag without internal J2000, plus the original SIDEREAL/RADIANS/J2000
/// bits and NONUT. Other returned bits survive. Performs no validation/I/O.
#[pyfunction(signature=(retflag,iflag))]
fn fixed_epoch_retflag(retflag: i32, iflag: i32) -> i32 {
    sidereal::fixed_epoch_retflag(retflag, iflag)
}
/// Return a J2000 mean-ecliptic to sidereal-projection rotation, or None.
///
/// sid_bits chooses SIDBIT_ECL_T0 (mean ecliptic of TT JD t0_jd), otherwise
/// SIDBIT_SSY_PLANE (invariable plane); ECL_T0 wins if both are set. The
/// zero_point_deg longitude IN that projection plane is rotated to the x
/// axis. With neither bit, returns None. t0_jd and zero_point_deg must still
/// be finite in every branch. The result is a dimensionless NumPy float64
/// (3, 3) row-major matrix acting on column vectors, not a tuple. Nonfinite
/// inputs or invalid frame arithmetic raise exceptions; no clock is inferred.
#[pyfunction(signature=(sid_bits,t0_jd,zero_point_deg))]
fn sidbit_ecliptic_matrix(
    py: Python<'_>,
    sid_bits: i32,
    t0_jd: f64,
    zero_point_deg: f64,
) -> PyResult<Option<Bound<'_, PyArray2<f64>>>> {
    call(py, || {
        sidereal::sidbit_ecliptic_matrix(sid_bits, t0_jd, zero_point_deg)
    })?
    .map(|m| array(py, m))
    .transpose()
}
/// Project a J2000 ecliptic zero-point direction onto the invariable plane.
///
/// ayan_j2000_deg is its finite mean-ecliptic longitude in degrees, latitude
/// zero. Returns invariable-plane longitude in [0, 360) degrees. Nonfinite
/// input raises ValueError. This geometric projection reads no clock, body
/// data or configured ayanamsha and does not establish model validity.
#[pyfunction(signature=(ayan_j2000_deg))]
fn ssy_plane_zero_point_deg(py: Python<'_>, ayan_j2000_deg: f64) -> PyResult<f64> {
    call(py, || sidereal::ssy_plane_zero_point_deg(ayan_j2000_deg))
}
/// Rotate a mean-J2000 ecliptic six-state into a supplied SIDBIT frame.
///
/// xx is a length-six sequence. With FLG_XYZ, it is (x,y,z,vx,vy,vz) in
/// AU and AU/TT-day; otherwise (lon,lat,r,dlon,dlat,dr) in DEGREES, AU and
/// corresponding per-TT-day rates, even if iflag contains FLG_RADIANS.
/// m_ecl is a dimensionless (3,3) row-major matrix from J2000 ecliptic to the
/// target projection, including its zero point. Returns a six-float tuple;
/// RADIANS changes spherical output angles/rates only. Equatorial requests
/// need the separate equatorial helper; this function does not reroute them.
///
/// Position and velocity use the SAME fixed matrix (no matrix time derivative).
/// Consumed nonfinite Cartesian/matrix values raise ValueError. In spherical
/// form an all-zero state or nonfinite input longitude/latitude bypasses the
/// rotation; zero radius retains direction but restores zero radius/radial
/// rate. Shape/extraction errors propagate. No provider or epoch is acquired.
#[pyfunction(signature=(xx,iflag,m_ecl))]
fn transform_sidbit_result(
    py: Python<'_>,
    xx: [f64; 6],
    iflag: i32,
    m_ecl: [[f64; 3]; 3],
) -> PyResult<Six> {
    call(py, || sidereal::transform_sidbit_result(xx, iflag, m_ecl)).map(six)
}
/// Return a mean-equatorial J2000 to mean-equatorial t0_jd rotation.
///
/// t0_jd is finite TT JD. Returns a dimensionless NumPy float64 (3,3)
/// row-major matrix acting on column vectors, identity at J2000. Uses the
/// Vondrak epoch chain with frame bias cancelling between epochs; no nutation
/// or Earth rotation is included. Invalid epochs/frame arithmetic raise
/// exceptions. Matrix entries preserve any consistent vector length units.
#[pyfunction(signature=(t0_jd))]
fn equatorial_epoch_matrix(py: Python<'_>, t0_jd: f64) -> PyResult<Bound<'_, PyArray2<f64>>> {
    array(py, call(py, || sidereal::equatorial_epoch_matrix(t0_jd))?)
}
/// Rotate a mean-J2000 EQUATORIAL six-state to the mean equator of t0_jd.
///
/// t0_jd is TT JD. xx is Cartesian AU/AU per TT day with FLG_XYZ, otherwise
/// (RA,Dec,r,dRA,dDec,dr) in DEGREES, AU and per-TT-day rates. iflag describes
/// the output; input spherical angles are degrees even with FLG_RADIANS.
/// Returns a six-float tuple, converting output angular slots/rates to radians
/// only for spherical RADIANS. No ayanamsha is subtracted; the same fixed
/// epoch matrix rotates position and velocity, without a rotation-rate term.
/// Nonfinite epochs, consumed matrix/Cartesian values or frame failures raise
/// exceptions. Spherical all-zero/undefined-direction states retain the
/// transform_sidbit_result sentinel rules. No provider is acquired.
#[pyfunction(signature=(xx,iflag,t0_jd))]
fn transform_equatorial_epoch_result(
    py: Python<'_>,
    xx: [f64; 6],
    iflag: i32,
    t0_jd: f64,
) -> PyResult<Six> {
    call(py, || {
        sidereal::transform_equatorial_epoch_result(xx, iflag, t0_jd)
    })
    .map(six)
}
/// Rotate a mean-J2000 six-state into a fixed sidereal mode's defining frame.
///
/// sid_mode is a base J2000, J1900, B1950 or Galactic Alignment/Mardyks ID.
/// iflag selects equatorial versus ecliptic and XYZ versus spherical output.
/// xx is AU/AU per TT day for XYZ, or degree/degree/AU and corresponding
/// per-TT-day rates for spherical input, even when RADIANS is requested.
/// Returns a six-float tuple; RADIANS affects spherical output angles/rates
/// only. Mardyks adds its defining ecliptic zero-point rotation; equatorial
/// output has no ayanamsha subtraction. J2000 is an identity apart from units.
///
/// Matrices are fixed in time: position and velocity rotate together without
/// a frame-rate term. Unsupported modes or consumed nonfinite/frame values
/// raise native exceptions; spherical undefined/zero states use the sentinel
/// rules of transform_sidbit_result. No configuration is changed or acquired.
#[pyfunction(signature=(xx,iflag,sid_mode))]
fn transform_fixed_epoch_result(
    py: Python<'_>,
    xx: [f64; 6],
    iflag: i32,
    sid_mode: i32,
) -> PyResult<Six> {
    call(py, || {
        sidereal::transform_fixed_epoch_result(xx, iflag, sid_mode)
    })
    .map(six)
}
/// Return eight named geometric chart angles in degrees at UT1 JD jd_ut.
///
/// lat/lon are finite north/east geographic degrees in [-90,90]/[-180,360].
/// Keys are 'Ascendant', 'MC', 'ARMC', 'Vertex', 'Equatorial_Ascendant',
/// 'Descendant', 'IC' and 'AntiVertex'. The last three are the antipodes of
/// Ascendant/MC/Vertex, NOT the co-Ascendants in a houses angle tuple.
/// Uses Placidus with Porphyry polar/invalid-cusp fallback and configured
/// clocks; coordinate, clock and non-fallback geometry errors propagate.
/// No planetary longitudes are included and caller observer state is unchanged.
#[pyfunction(name="calc_angles",signature=(jd_ut,lat,lon))]
fn angles_calc_angles(
    py: Python<'_>,
    jd_ut: f64,
    lat: f64,
    lon: f64,
) -> PyResult<std::collections::BTreeMap<&'static str, f64>> {
    session(py, |s| s.numeric_angles(jd_ut, lat, lon)).map(native::named_angles)
}
/// Return one named chart-angle value in degrees, or zero for an unknown ID.
///
/// angle_id is a public angle constant; jd_ut is UT1 JD, lat/lon north/east
/// geographic degrees. Uses the same eight-angle geometry/fallback and
/// coordinate/clock errors as calc_angles. Unknown IDs still calculate that
/// geometry first; they do not bypass errors. No planet or observer mutation.
#[pyfunction(signature=(angle_id,jd_ut,lat,lon))]
fn get_angle_value(py: Python<'_>, angle_id: i32, jd_ut: f64, lat: f64, lon: f64) -> PyResult<f64> {
    session(py, |s| s.numeric_angle_value(angle_id, jd_ut, lat, lon))
}

/// Convert an ecliptic Cartesian state to spherical position and derivatives.
///
/// state is (x,y,z,vx,vy,vz) in AU and AU per TT day, all finite, in one
/// ecliptic frame. Returns (lon,lat,r,dlon,dlat,dr) as a six-float tuple in
/// degrees, AU and corresponding per-TT-day units, retaining the actual
/// radial dot-product derivative. It changes representation, not frame/clock.
///
/// Nonfinite components raise ValueError. Zero vector or zero cylindrical
/// radius (including an ecliptic pole, or squared-radius underflow) raises
/// CalculationError because a longitude derivative is undefined. There is no
/// blanket overflow guard for arbitrarily large finite values; derived IEEE
/// values may be nonfinite. No data files or epoch samples are consumed.
#[pyfunction(signature=(state))]
fn state_to_polar(py: Python<'_>, state: [f64; 6]) -> PyResult<Six> {
    py.detach(|| native::interpolated_state_to_polar(state))
        .map_err(|error| {
            if matches!(
                error,
                libephemeris_core::errors::Error::ZeroVectorToSpherical
            ) {
                match py
                    .import("libephemeris._core")
                    .and_then(|m| m.getattr("CalculationError"))
                    .and_then(|class| class.call1((error.to_string(),)))
                {
                    Ok(value) => PyErr::from_value(value),
                    Err(error) => error,
                }
            } else {
                to_pyerr(py, error)
            }
        })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    macro_rules! add { ($($name:ident),* $(,)?) => { $(module.add(concat!("_nf_",stringify!($name)),wrap_pyfunction!($name,module)?)?;)* }; }
    add!(
        state_to_polar,
        astrometry_nutation_angles,
        precess_from_j2000,
        apply_aberration_to_position,
        get_cached_nutation,
        get_nutation_degrees,
        get_cached_obliquity,
        get_mean_obliquity,
        get_true_obliquity,
        get_cached_time_ut1,
        get_cached_time_tt,
        clear_caches,
        get_cache_info,
        has_erfa,
        get_erfa_nutation_nut00a,
        get_erfa_nutation_nut06a,
        get_erfa_pnm06a_matrix,
        get_erfa_nutation_cached,
        get_erfa_obliquity_iau2006,
        compare_nutation_models,
        clear_erfa_nutation_cache,
        get_erfa_nutation_cache_info,
        lunar_delaunay_arguments,
        mean_lunar_apse_plane_node,
        mean_lunar_apogee_position,
        mean_lunar_apogee_state,
        mean_lunar_node_position,
        mean_lunar_node_state,
        mean_orbital_elements,
        vondrak_mean_obliquity_rad,
        vondrak_mean_obliquity_deg,
        vondrak_precession_matrix,
        vondrak_pn_matrix,
        method_b_accumulated_precession,
        mean_obliquity_series_rad,
        mean_obliquity_rad,
        mean_obliquity_deg,
        precession_matrix,
        mean_sidereal_time_deg,
        apparent_sidereal_time_deg,
        is_fixed_epoch_request,
        fixed_epoch_request_flags,
        fixed_epoch_retflag,
        sidbit_ecliptic_matrix,
        ssy_plane_zero_point_deg,
        transform_sidbit_result,
        equatorial_epoch_matrix,
        transform_equatorial_epoch_result,
        transform_fixed_epoch_result,
        angles_calc_angles,
        get_angle_value
    );
    let py = module.py();
    for (name, value) in native::astrometry_constants() {
        module.add(name, value)?;
    }
    module.add(
        "_nf_OBLIQUITY_COEFFS",
        PyTuple::new(py, native::OBLIQUITY_COEFFS)?,
    )?;
    module.add(
        "_nf_ZETA_A_COEFFS",
        PyTuple::new(py, native::ZETA_A_COEFFS)?,
    )?;
    module.add("_nf_Z_A_COEFFS", PyTuple::new(py, native::Z_A_COEFFS)?)?;
    module.add(
        "_nf_THETA_A_COEFFS",
        PyTuple::new(py, native::THETA_A_COEFFS)?,
    )?;
    for (name, value) in native::lunar_constants() {
        module.add(name, value)?;
    }
    module.add(
        "_nf_INVARIABLE_PLANE_INCL_DEG",
        sidereal::INVARIABLE_PLANE_INCL_DEG,
    )?;
    module.add(
        "_nf_INVARIABLE_PLANE_NODE_DEG",
        sidereal::INVARIABLE_PLANE_NODE_DEG,
    )?;
    module.add("_nf_FIXED_EPOCH_T0", native::fixed_epoch_times())?;
    module.add("_nf_FIXED_EPOCH_LON_OFFSET", native::fixed_epoch_offsets())?;
    module.add("_nf_MARDYKS_DEFINING", sidereal::MARDYKS_DEFINING)?;
    module.add(
        "_nf_SUPPORTED_LUNAR_APSIDES",
        (
            libephemeris_core::constants::INTP_APOG,
            libephemeris_core::constants::INTP_PERG,
        ),
    )?;
    Ok(())
}
