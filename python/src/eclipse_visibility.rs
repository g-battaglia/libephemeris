//! Thin wrappers for the eclipse module's scalar, eight-detail optical APIs.
use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::python_support::eclipse_visibility::EclipseVisibilityError;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

fn error(py: Python<'_>, value: EclipseVisibilityError) -> PyErr {
    match value {
        EclipseVisibilityError::Native(error) => to_pyerr(py, error),
        EclipseVisibilityError::NativeMainDisabled => PyRuntimeError::new_err(value.to_string()),
    }
}
/// Search the eclipse module's compact scalar heliacal model.
///
/// Parameters
/// ----------
/// jd_start : float
///     Starting UT1 Julian day, not TT/UTC.
/// lat, lon, altitude : float
///     Geographic latitude north degrees, longitude east degrees and height m.
/// pressure, temperature, humidity : float
///     Public atmospheric arguments: mbar, Celsius and humidity FRACTION
///     (default 0.5). These arguments do not enter this model's visibility
///     relation; they do not configure the sequence-based retinal model.
/// body : int
///     Planet ID; Sun, Moon and unknown IDs are rejected. The default body=0
///     is therefore not a usable event target.
/// event_type : int
///     1 morning first, 2 evening last, 3 evening first, 4 morning last.
/// flags : int
///     Calculation flags used for photometry, default 2.
///
/// Returns
/// -------
/// tuple
///     (event_ut1_jd, event_type), or (0., -1) after a completed sampled miss.
///
/// Raises
/// ------
/// RuntimeError
///     Forced LEB mode cannot supply this direct JPL/SPK geometry.
/// ValueError, libephemeris.Error
///     Invalid event/body/site or data/calculation failure.
///
/// Notes
/// -----
/// This libephemeris.eclipse API differs from root heliacal_ut's sequence
/// arguments and visibility-window triple. Uses a 400-day three-point daily
/// dawn/dusk scan (evening-last may extend 120 days), simple twilight magnitude
/// and elongation >=10 degrees. Bounded bisection does not guarantee a physical
/// event, convergence or observational visibility. Geometry and configured
/// photometry retain their separate clocks; failures are not completed misses.
#[pyfunction(signature=(jd_start,lat,lon,altitude=0.0,pressure=1013.25,temperature=15.0,humidity=0.5,body=0,event_type=1,flags=2))]
#[allow(clippy::too_many_arguments)]
fn heliacal_ut(
    py: Python<'_>,
    jd_start: f64,
    lat: f64,
    lon: f64,
    altitude: f64,
    pressure: f64,
    temperature: f64,
    humidity: f64,
    body: i32,
    event_type: i32,
    flags: i32,
) -> PyResult<(f64, i32)> {
    py.detach(|| {
        runtime()
            .shared_session()
            .map_err(EclipseVisibilityError::Native)?
            .eclipse_heliacal_ut(
                jd_start,
                lat,
                lon,
                altitude,
                pressure,
                temperature,
                humidity,
                body,
                event_type,
                flags,
            )
    })
    .map_err(|value| error(py, value))
}
/// Return a compact vision flag and EIGHT fixed-instant detail slots.
///
/// Parameters
/// ----------
/// jd : float
///     UT1 Julian day, not TT/UTC.
/// geopos : sequence of float
///     (longitude east degrees, latitude north degrees, height metres),
///     missing words default zero.
/// atmo : sequence of float
///     (pressure mbar, temperature Celsius, humidity PERCENT, met_range),
///     defaults (1013.25, 15, 50, 0). Positive range <1 supplies total
///     extinction mag/airmass; >=1 gives meteorological visibility km.
/// observer : sequence of float
///     Age years and Snellen ratio, default 36/1. Additional optical words
///     are not used by this compact model.
/// objname : str
///     Ordinary planet name, supported numeric-ID string or catalog-star name.
/// flags : int
///     Calculation/visibility mask; DARK and NOMOON suppress the respective
///     background. Default 2. This is not the full retinal optical model.
///
/// Returns
/// -------
/// tuple
///     (vision_flag_as_float, details8): limiting magnitude, object altitude,
///     object azimuth, Sun altitude, Sun azimuth, Moon altitude, Moon azimuth,
///     extinction-adjusted object magnitude. Angles are degrees; planet
///     azimuths start north, star azimuths start south toward west. Below the
///     horizon the flag is HELFLAG_BELOW_HORIZON, limit is 0 and the last slot
///     is intrinsic magnitude. DARK/NOMOON set corresponding altitude to -90.
///
/// Raises
/// ------
/// RuntimeError
///     Direct JPL/SPK is disabled by forced LEB mode.
/// ValueError, libephemeris.Error
///     Invalid object/site or provider/data/calculation failure.
///
/// Notes
/// -----
/// Unlike root vis_limit_mag (ten details), this uses approximate piecewise
/// solar/lunar backgrounds, simple extinction and linear age correction.
/// Native apparent geometry and public photometry/horizon clocks are distinct.
/// No event search, telescope model or observational guarantee is supplied.
#[pyfunction(signature=(jd,geopos,atmo,observer,objname,flags=2))]
fn vis_limit_mag<'py>(
    py: Python<'py>,
    jd: f64,
    geopos: Vec<f64>,
    atmo: Vec<f64>,
    observer: Vec<f64>,
    objname: &str,
    flags: i32,
) -> PyResult<(f64, Bound<'py, PyTuple>)> {
    let (status, values) = py
        .detach(|| {
            runtime()
                .shared_session()
                .map_err(EclipseVisibilityError::Native)?
                .eclipse_vis_limit_mag(jd, &geopos, &atmo, &observer, objname, flags)
        })
        .map_err(|value| error(py, value))?;
    Ok((status, PyTuple::new(py, values)?))
}
pub(crate) fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new(parent.py(), "libephemeris._core._eclipse_visibility")?;
    module.add_function(wrap_pyfunction!(heliacal_ut, &module)?)?;
    module.add_function(wrap_pyfunction!(vis_limit_mag, &module)?)?;
    parent
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item("libephemeris._core._eclipse_visibility", &module)?;
    parent.add("_eclipse_visibility", module)?;
    Ok(())
}
