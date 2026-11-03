//! Map conversion around the native four-lot and geometric sect calculation.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::arabic_parts;
use libephemeris_core::python_support::ArabicPartsError;
use pyo3::exceptions::{PyKeyError, PyValueError};
use pyo3::prelude::*;
use std::collections::HashMap;

/// Calculate four Arabic lots from supplied ecliptic longitudes.
///
/// Parameters
/// ----------
/// positions : dict[str, float]
///     Required keys are ``Asc``, ``Sun``, ``Moon``, ``Mercury`` and ``Venus``;
///     values are degrees in one consistent ecliptic frame. Optional ``Sun_lat``
///     is solar latitude in degrees, default 0. Extra keys are ignored.
/// jd : float or None
///     UT1 Julian day, used only by the high-latitude sect calculation.
/// geo_lat, geo_lon : float or None
///     Geographic north latitude and east longitude in degrees. All three
///     optional values default to None. Only a complete set with
///     ``abs(geo_lat) > 60`` selects geometric, unrefracted solar altitude at
///     zero height; otherwise sect uses the supplied Ascendant/Sun arc.
///
/// Returns
/// -------
/// dict[str, float]
///     ``Pars_Fortunae``, ``Pars_Spiritus``, ``Pars_Amoris`` and ``Pars_Fidei``,
///     each a longitude in [0, 360) degrees. Fortune and Spirit reverse their
///     formulas at night; Love and Faith use the fixed project conventions.
///
/// Raises
/// ------
/// KeyError
///     A required position key is missing.
/// ValueError
///     A supplied required longitude or Sun_lat is nonfinite.
/// libephemeris.Error
///     Native configuration, clock or horizon geometry fails; the boundary may
///     also report its documented built-in exceptions. No planetary positions
///     are acquired: the supplied dictionary defines the chart.
#[pyfunction(signature = (positions, *, jd = None, geo_lat = None, geo_lon = None))]
fn calc_all_arabic_parts(
    py: Python<'_>,
    positions: HashMap<String, f64>,
    jd: Option<f64>,
    geo_lat: Option<f64>,
    geo_lon: Option<f64>,
) -> PyResult<HashMap<String, f64>> {
    let result = py.detach(|| runtime().calc_all_arabic_parts(&positions, jd, geo_lat, geo_lon));
    result.map_err(|error| match error {
        ArabicPartsError::MissingPositions(_) => PyKeyError::new_err(error.to_string()),
        ArabicPartsError::Input(_) => PyValueError::new_err(error.to_string()),
        ArabicPartsError::Native(error) => to_pyerr(py, error),
    })
}

/// Return Fortune: ``asc + moon - sun`` by day, ``asc + sun - moon`` by night.
///
/// asc, sun and moon are finite ecliptic longitudes in degrees in one frame.
/// is_diurnal defaults to True. The result is wrapped to [0, 360) degrees;
/// nonfinite inputs raise ValueError. This scalar formula reads no clock/data.
#[pyfunction(signature=(asc,sun,moon,is_diurnal=true))]
fn calc_arabic_part_of_fortune(
    py: Python<'_>,
    asc: f64,
    sun: f64,
    moon: f64,
    is_diurnal: bool,
) -> PyResult<f64> {
    py.detach(|| arabic_parts::part_of_fortune(asc, sun, moon, is_diurnal))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
/// Return Spirit: ``asc + sun - moon`` by day, ``asc + moon - sun`` by night.
///
/// asc, sun and moon are finite ecliptic longitudes in degrees in one frame.
/// is_diurnal defaults to True. The result is wrapped to [0, 360) degrees;
/// nonfinite inputs raise ValueError. This scalar formula reads no clock/data.
#[pyfunction(signature=(asc,sun,moon,is_diurnal=true))]
fn calc_arabic_part_of_spirit(
    py: Python<'_>,
    asc: f64,
    sun: f64,
    moon: f64,
    is_diurnal: bool,
) -> PyResult<f64> {
    py.detach(|| arabic_parts::part_of_spirit(asc, sun, moon, is_diurnal))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
/// Return Love as ``(asc + venus - sun) % 360`` in degrees.
///
/// All inputs are finite ecliptic longitudes in one frame. This project
/// convention is independent of sect and reads no clock/data. Nonfinite
/// inputs raise ValueError.
#[pyfunction]
fn calc_arabic_part_of_love(py: Python<'_>, asc: f64, venus: f64, sun: f64) -> PyResult<f64> {
    py.detach(|| arabic_parts::part_of_love(asc, venus, sun))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
/// Return Faith as ``(asc + mercury - moon) % 360`` in degrees.
///
/// All inputs are finite ecliptic longitudes in one frame. This project
/// convention is independent of sect and reads no clock/data. Nonfinite
/// inputs raise ValueError.
#[pyfunction]
fn calc_arabic_part_of_faith(py: Python<'_>, asc: f64, mercury: f64, moon: f64) -> PyResult<f64> {
    py.detach(|| arabic_parts::part_of_faith(asc, mercury, moon))
        .map_err(|error| PyValueError::new_err(error.to_string()))
}
/// Determine day/night sect from supplied solar and Ascendant longitudes.
///
/// sun_lon and asc are ecliptic degrees in one consistent frame; sun_lat is
/// solar latitude in degrees, default 0. With jd (UT1 JD), geo_lat (north
/// degrees) and geo_lon (east degrees) all supplied and abs(geo_lat) > 60,
/// use unrefracted solar altitude at zero height: altitude >= 0 is day.
/// Otherwise return False only for ``0 < (sun_lon - asc) % 360 < 180``.
/// Exact horizon/arc boundaries count as day. Incomplete optional inputs
/// select the arc rule rather than requesting a partial observer calculation.
///
/// Returns bool. Consumed nonfinite arc angles raise ValueError; the 3-D
/// branch propagates native coordinate, clock and projection exceptions.
/// Unused optional values and sun_lat in the arc branch are not validated.
#[pyfunction(signature=(sun_lon,asc,*,jd=None,geo_lat=None,geo_lon=None,sun_lat=0.0))]
fn is_day_chart(
    py: Python<'_>,
    sun_lon: f64,
    asc: f64,
    jd: Option<f64>,
    geo_lat: Option<f64>,
    geo_lon: Option<f64>,
    sun_lat: f64,
) -> PyResult<bool> {
    py.detach(|| runtime().is_day_chart(sun_lon, asc, jd, geo_lat, geo_lon, sun_lat))
        .map_err(|error| match error {
            ArabicPartsError::MissingPositions(_) => PyKeyError::new_err(error.to_string()),
            ArabicPartsError::Input(_) => PyValueError::new_err(error.to_string()),
            ArabicPartsError::Native(error) => to_pyerr(py, error),
        })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(calc_all_arabic_parts, module)?)?;
    let helpers = PyModule::new(module.py(), "libephemeris._core._arabic_parts")?;
    helpers.add_function(wrap_pyfunction!(calc_arabic_part_of_fortune, &helpers)?)?;
    helpers.add_function(wrap_pyfunction!(calc_arabic_part_of_spirit, &helpers)?)?;
    helpers.add_function(wrap_pyfunction!(calc_arabic_part_of_love, &helpers)?)?;
    helpers.add_function(wrap_pyfunction!(calc_arabic_part_of_faith, &helpers)?)?;
    helpers.add_function(wrap_pyfunction!(is_day_chart, &helpers)?)?;
    module
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item("libephemeris._core._arabic_parts", &helpers)?;
    module.add("_arabic_parts", helpers)?;
    Ok(())
}
