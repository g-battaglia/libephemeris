//! Lunar-node and apogee calls over the core's geometric and analytic models.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::lunar;
use pyo3::exceptions::{PyUserWarning, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

pyo3::create_exception!(
    libephemeris,
    MeeusPolynomialWarning,
    PyUserWarning,
    "UserWarning subtype for polynomial-model validity advisories. The legacy optimal interval is +/-2 Julian centuries from J2000, not a universal accuracy bound. Direct ERFA/IERS mean-point helpers do not emit this warning merely because that interval is exceeded."
);
pyo3::create_exception!(
    libephemeris,
    MeeusRangeError,
    PyValueError,
    "Deprecated ValueError subtype retained for lunar polynomial-range compatibility. Modern mean-point helpers enforce their selected main-image coverage through native range errors; these legacy validity bands do not define file coverage or guarantee accuracy."
);

/// Calculate the osculating ascending lunar node at ``jd_tt`` (TT Julian day).
///
/// Return ``(longitude, latitude, distance)`` on the true ecliptic of date,
/// in degrees, degrees and AU; node latitude is zero. Uses geometric Moon/Earth
/// states from the configured ephemeris, not an apparent lunar position.
/// Invalid time, missing data, coverage/frame errors and singular orbital
/// geometry raise the corresponding native exception.
#[pyfunction]
fn calc_true_lunar_node(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64, f64)> {
    py.detach(|| runtime().shared_session()?.true_lunar_point(jd_tt, false))
        .map_err(|error| to_pyerr(py, error))
}

/// Calculate the osculating lunar apogee (True Lilith) at a TT Julian day.
///
/// Return ``(longitude_deg, latitude_deg, distance_au)`` on the true ecliptic
/// of date from the instantaneous geometric eccentricity vector. This is not
/// the next physical apogee event. ``calc_true_lilith_orbital_elements`` is
/// an alias of this function, not an independent model. Invalid time, missing
/// ephemeris/coverage, frame errors and singular conic geometry propagate.
#[pyfunction]
fn calc_true_lilith(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64, f64)> {
    py.detach(|| runtime().shared_session()?.true_lunar_point(jd_tt, true))
        .map_err(|error| to_pyerr(py, error))
}

/// Evaluate the DE440-fitted interpolated apogee at ``jd_tt`` (TT Julian day).
///
/// Return ``(longitude_deg, latitude_deg, distance_au)`` on the true ecliptic
/// of date. Distance is conventional and constant, not the Moon's range at an
/// event. The residual table tapers to zero outside its grid; the underlying
/// harmonic model continues without an accuracy guarantee. No kernel is loaded.
/// NaN/infinite time raises an exception. Example:
///
/// ```python
/// from libephemeris.lunar import calc_interpolated_apogee
/// lon, lat, au = calc_interpolated_apogee(2451545.0)
/// ```
#[pyfunction]
fn calc_interpolated_apogee(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64, f64)> {
    py.detach(|| lunar::calc_interpolated_apogee(jd_tt))
        .map_err(|error| to_pyerr(py, error))
}

/// Evaluate the independent DE440-fitted interpolated perigee at a TT Julian day.
///
/// Return ``(longitude_deg, latitude_deg, distance_au)`` on the true ecliptic
/// of date. This is an abstract apsidal track with constant conventional range,
/// not a perigee/tidal event finder. The residual tapers outside its source grid;
/// harmonic extrapolation is not precision-qualified. Non-finite time raises;
/// no ephemeris file is required.
#[pyfunction]
fn calc_interpolated_perigee(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64, f64)> {
    py.detach(|| lunar::calc_interpolated_perigee(jd_tt))
        .map_err(|error| to_pyerr(py, error))
}

/// Return residual-table metadata for ``INTP_APOG`` or ``INTP_PERG``.
///
/// The tuple is ``(first_tt_jd, last_tt_jd, step_days, taper_days)``; any other
/// ``body_id`` returns ``None``. Residuals interpolate linearly between knots
/// and taper linearly to zero outside the grid, creating derivative breaks.
/// This describes the fitted source, not guaranteed physical coverage; no file
/// is opened and no epoch is evaluated.
#[pyfunction]
fn interpolated_apse_source_grid(py: Python<'_>, body_id: i32) -> Option<(f64, f64, f64, f64)> {
    py.detach(|| lunar::interpolated_apse_source_grid(body_id))
}

/// Return mean ascending-node longitude in degrees, normalized to [0, 360).
///
/// ``jd_tt`` is a TT Julian day. The ERFA/IERS analytic argument is on the mean
/// ecliptic of date, independent of observer and configured Delta-T. Active JPL
/// image bounds may restrict the epoch; ordinary LEB windows do not restrict
/// this analytic model. Invalid time, range and configuration errors raise.
#[pyfunction]
fn calc_mean_lunar_node(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    py.detach(|| runtime().shared_session()?.calc_mean_lunar_node(jd_tt))
        .map_err(|error| to_pyerr(py, error))
}

/// Return mean lunar-apogee longitude in degrees on the mean ecliptic of date.
///
/// ``jd_tt`` is a TT Julian day. Uses the ERFA/IERS eccentricity-plane argument
/// and conventional orbit geometry, not the instantaneous lunar ellipse.
/// Observer/Delta-T settings do not enter; active main-image coverage and
/// acquisition policy still apply. Non-finite time, range and data errors raise.
#[pyfunction]
fn calc_mean_lilith(py: Python<'_>, jd_tt: f64) -> PyResult<f64> {
    py.detach(|| runtime().shared_session()?.calc_mean_lilith(jd_tt))
        .map_err(|error| to_pyerr(py, error))
}

/// Return mean-apogee ``(longitude, latitude)`` in degrees at a TT Julian day.
///
/// The direction lies on the mean ecliptic-of-date coordinate system and uses
/// a conventionally inclined lunar orbit. Coverage, configuration and invalid
/// time errors are the same as ``calc_mean_lilith``; no observer correction is
/// applied.
#[pyfunction]
fn calc_mean_lilith_with_latitude(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64)> {
    py.detach(|| {
        runtime()
            .shared_session()?
            .calc_mean_lilith_with_latitude(jd_tt)
    })
    .map_err(|error| to_pyerr(py, error))
}

/// Return osculating perigee ``(longitude_deg, latitude_deg, distance_au)``.
///
/// ``jd_tt`` is TT Julian day. The configured geometric Moon/Earth state gives
/// the eccentricity-vector direction on the true ecliptic of date and conic
/// range ``p/(1+e)``. This is not a search for the next distance minimum.
/// Invalid time, missing resources, coverage/frame and singular-orbit errors raise.
#[pyfunction]
fn calc_osculating_perigee(py: Python<'_>, jd_tt: f64) -> PyResult<(f64, f64, f64)> {
    py.detach(|| runtime().shared_session()?.calc_osculating_perigee(jd_tt))
        .map_err(|error| to_pyerr(py, error))
}

/// Compare the maintained apogee spellings at ``jd_tt`` (TT Julian day).
///
/// Return a dict with ``eccentricity_vector`` and ``orbital_elements`` triples
/// ``(longitude_deg, latitude_deg, distance_au)``, plus ``lon_diff``/``lat_diff``
/// in degrees and ``e_diff`` in AU (despite its name, it is a distance difference,
/// not eccentricity), on the true ecliptic of date. Both spellings are the same
/// model, so differences are zero
/// by construction, not an independent accuracy test. State/geometry errors raise.
#[pyfunction]
fn compare_true_lilith_methods(py: Python<'_>, jd_tt: f64) -> PyResult<Py<PyDict>> {
    let value = py
        .detach(|| {
            runtime()
                .shared_session()?
                .compare_true_lilith_methods(jd_tt)
        })
        .map_err(|error| to_pyerr(py, error))?;
    let result = PyDict::new(py);
    result.set_item("eccentricity_vector", value.eccentricity_vector)?;
    result.set_item("orbital_elements", value.orbital_elements)?;
    result.set_item("lon_diff", value.lon_diff)?;
    result.set_item("lat_diff", value.lat_diff)?;
    result.set_item("e_diff", value.dist_diff)?;
    Ok(result.unbind())
}

macro_rules! mean_state {
    ($name:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyfunction(signature=(jd_tt, *, speed3=false))]
        fn $name(
            py: Python<'_>,
            jd_tt: f64,
            speed3: bool,
        ) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
            py.detach(|| runtime().shared_session()?.$name(jd_tt, speed3))
                .map_err(|error| to_pyerr(py, error))
        }
    };
}
mean_state!(
    calc_mean_lunar_node_state,
    "Return mean-node (lon, lat, range, dlon, dlat, drange) at a TT Julian day.\n\nCoordinates are degrees/degrees/AU on the mean ecliptic of date and rates are per TT day. Latitude, latitude rate and radial rate are zero; range is conventional. Keyword-only speed3=False is accepted but ignored: both values use the same smooth argument derivative. Active main-image coverage, invalid time and configuration errors raise."
);
mean_state!(
    calc_mean_lilith_state,
    "Return mean-apogee (lon, lat, range, dlon, dlat, drange) at a TT Julian day.\n\nCoordinates are degrees/degrees/AU on the mean ecliptic of date and rates are per TT day. Range is conventional and radial rate is zero. Keyword-only speed3=False is accepted but ignored: the same smooth model derivative is returned. Invalid time, active main-image coverage and configuration errors raise."
);

/// Return an interpolated apsidal state at ``jd_tt`` (TT Julian day).
///
/// ``body_id`` must be ``INTP_APOG`` or ``INTP_PERG``; otherwise ``ValueError``
/// is raised. Layout is ``(lon, lat, range, dlon, dlat, drange)`` in degrees,
/// degrees, AU and corresponding units per TT day, on the true ecliptic of date.
/// ``speed3=False`` is accepted but ignored; both settings use the same centered
/// model stencil. The residual-table taper/extrapolation limits of the position
/// helpers apply; these are not event states. Non-finite time raises.
#[pyfunction(signature=(jd_tt, body_id, *, speed3=false))]
fn calc_interpolated_apse_state(
    py: Python<'_>,
    jd_tt: f64,
    body_id: i32,
    speed3: bool,
) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
    py.detach(|| lunar::calc_interpolated_apse_state(jd_tt, body_id, speed3))
        .map_err(|error| match error {
            libephemeris_core::errors::Error::UnknownBody { .. } => {
                PyValueError::new_err(error.to_string())
            }
            error => to_pyerr(py, error),
        })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "MeeusPolynomialWarning",
        module.py().get_type::<MeeusPolynomialWarning>(),
    )?;
    module.add("MeeusRangeError", module.py().get_type::<MeeusRangeError>())?;
    macro_rules! add {
        ($($name:ident),+ $(,)?) => { $(module.add_function(wrap_pyfunction!($name, module)?)?;)+ };
    }
    add!(
        calc_true_lunar_node,
        calc_true_lilith,
        calc_interpolated_apogee,
        calc_interpolated_perigee,
        interpolated_apse_source_grid,
        calc_mean_lunar_node,
        calc_mean_lilith,
        calc_mean_lilith_with_latitude,
        calc_osculating_perigee,
        compare_true_lilith_methods,
        calc_mean_lunar_node_state,
        calc_mean_lilith_state,
        calc_interpolated_apse_state
    );
    module.add(
        "calc_true_lilith_orbital_elements",
        module.getattr("calc_true_lilith")?,
    )?;
    module.add("MEEUS_OPTIMAL_CENTURIES", lunar::MEEUS_OPTIMAL_CENTURIES)?;
    module.add("MEEUS_VALID_CENTURIES", lunar::MEEUS_VALID_CENTURIES)?;
    module.add("MEEUS_MAX_CENTURIES", lunar::MEEUS_MAX_CENTURIES)?;
    module.add("GM_EARTH_MOON_AU3_DAY2", lunar::GM_EARTH_MOON_AU3_DAY2)?;
    Ok(())
}
