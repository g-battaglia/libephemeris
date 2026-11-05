//! Typed position calls and process-wide configuration forwarding.

use libephemeris_core::config::{ConfigurationSources, Runtime};
use libephemeris_core::python_support::Position;
use pyo3::prelude::*;
use std::path::PathBuf;
use std::sync::OnceLock;

pub(crate) type PositionResult = ((f64, f64, f64, f64, f64, f64), i32);

pub(crate) fn runtime() -> &'static Runtime {
    static DEFAULT: OnceLock<Runtime> = OnceLock::new();
    DEFAULT.get_or_init(|| Runtime::with_sources(ConfigurationSources::discover()))
}

pub(crate) fn position(py: Python<'_>, result: Position) -> PyResult<PositionResult> {
    for (message, category) in result
        .warnings
        .into_iter()
        .map(|warning| (warning.message, false))
        .chain(result.messages.into_iter().map(|message| (message, true)))
    {
        let message = std::ffi::CString::new(message)
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
        let category = if category {
            py.get_type::<pyo3::exceptions::PyRuntimeWarning>()
        } else {
            py.get_type::<pyo3::exceptions::PyUserWarning>()
        };
        PyErr::warn(py, &category, &message, 1)?;
    }
    let [a, b, c, d, e, f] = result.values;
    Ok(((a, b, c, d, e, f), result.returned_flags))
}

/// Calculate a body's position and motion at a UT1 Julian day.
///
/// Parameters
/// ----------
/// tjdut : float
///     Julian day in UT1, not UTC. Use utc_to_jd for a civil UTC label.
/// planet : int
///     Body identifier, for example SUN, MOON or MARS; ECL_NUT requests
///     obliquity and nutation instead of a body position.
/// flags : int, optional
///     Bitwise OR of FLG_* constants. The default is FLG_SWIEPH | FLG_SPEED.
///     Provider selection follows set_calc_mode, not the compatibility bit's name.
///     FLG_J2000/FLG_NONUT select epoch/mean-frame reductions;
///     FLG_TRUEPOS suppresses apparent corrections, while FLG_NOABERR and
///     FLG_NOGDEFL disable aberration and deflection individually.
///
/// Returns
/// -------
/// tuple
///     (values, returned_flags). Default output is apparent, geocentric and
///     ecliptic-of-date: (longitude, latitude, distance, longitude_rate,
///     latitude_rate, distance_rate), in degrees,
///     AU, degrees/TT-day and AU/TT-day. FLG_EQUATORIAL selects RA/Dec
///     (RA in degrees, not hours); FLG_RADIANS changes angular units to radians.
///     FLG_XYZ instead returns (x, y, z, vx, vy, vz) in AU and AU/TT-day.
///     Rates require a speed flag; returned_flags describes the served request.
///     For ECL_NUT, values is (true_obliquity, mean_obliquity, dpsi, deps,
///     0.0, 0.0), always in degrees; ordinary output shaping and speed flags
///     do not turn these frame quantities into coordinates or rates.
///     Angle/Arabic-part IDs are another special case: degree longitude plus
///     five placeholders, without ordinary XYZ/RADIANS shaping. Angles require
///     the configured observer even without FLG_TOPOCTR. Arabic parts consume
///     cached chart inputs rather than recomputing them at tjdut.
///
/// Raises
/// ------
/// ValueError
///     For nonfinite numerical input on checked calculation paths.
/// Error
///     For unsupported requests or required configuration/data failures;
///     subclasses include UnknownBodyError, ConfigurationError,
///     EphemerisRangeError and NetworkSealedError as applicable.
///
/// Notes
/// -----
/// FLG_TOPOCTR uses the observer configured by set_topo; FLG_SIDEREAL uses
/// set_sid_mode. FLG_HELCTR/FLG_BARYCTR select Sun/barycenter origins.
/// Actual body/support coverage, not a tier label, limits dates.
/// Calls are synchronous, releasing the GIL during native calculation.
///
/// Examples
/// --------
/// With a previously configured local ephemeris:
/// ```python
/// import libephemeris as eph
/// jd_tt, jd_ut1 = eph.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
/// values, returned_flags = eph.calc_ut(jd_ut1, eph.MARS)
/// lon, lat, distance, lon_rate, lat_rate, distance_rate = values
/// ```
#[pyfunction(signature = (tjdut, planet, flags = 258))]
fn calc_ut(py: Python<'_>, tjdut: f64, planet: i32, flags: i32) -> PyResult<PositionResult> {
    py.detach(|| runtime().calc_ut(tjdut, planet, flags))
        .map_err(|error| crate::errors::to_pyerr(py, error))
        .and_then(|value| position(py, value))
}

/// Calculate a body's position and motion at a TT Julian day.
///
/// Parameters
/// ----------
/// tjdet : float
///     Terrestrial Time Julian day, not UT1 or UTC.
/// planet : int
///     Body identifier, or ECL_NUT for obliquity/nutation.
/// flags : int, optional
///     FLG_* bitmask; default FLG_SWIEPH | FLG_SPEED. Observer, frame,
///     sidereal and correction selections follow calc_ut.
///
/// Returns
/// -------
/// tuple
///     (values, returned_flags), with the six-component layout documented
///     by calc_ut: spherical degrees/AU and rates per TT day; equatorial,
///     radians and Cartesian alternatives follow the corresponding flags.
///     ECL_NUT returns (true_obliquity, mean_obliquity, dpsi, deps, 0.0,
///     0.0), always in degrees, even with FLG_RADIANS; these are not rates.
///     Angle/Arabic-part IDs retain the special degree/placeholder layout
///     and observer/chart-cache requirements described by calc_ut.
///
/// Raises
/// ------
/// ValueError
///     For nonfinite numerical input on checked paths.
/// Error
///     For unsupported requests, missing configuration/data or coverage
///     failures, with the same exception subclasses as calc_ut.
#[pyfunction(signature = (tjdet, planet, flags = 258))]
fn calc(py: Python<'_>, tjdet: f64, planet: i32, flags: i32) -> PyResult<PositionResult> {
    py.detach(|| runtime().calc(tjdet, planet, flags))
        .map_err(|error| crate::errors::to_pyerr(py, error))
        .and_then(|value| position(py, value))
}

/// Select the calculation provider for the process-wide default runtime.
///
/// Parameters
/// ----------
/// mode : str or None
///     'leb', 'skyfield', 'auto' or 'horizons'. 'skyfield' names the native
///     Rust JPL/SPK route, not a Python engine. None clears the explicit
///     override so captured environment/dotenv/TOML settings and then the
///     default 'auto' apply. No ephemeris is opened by this setter.
///
/// Raises
/// ------
/// ValueError
///     If mode is unknown; the previous configuration remains selected.
#[pyfunction(signature = (mode))]
fn set_calc_mode(py: Python<'_>, mode: Option<&str>) -> PyResult<()> {
    runtime()
        .set_calc_mode(mode)
        .map_err(|error| crate::errors::to_pyerr(py, error))
}

/// Select the default runtime's 'base', 'medium' or 'extended' precision tier.
///
/// tier is a string naming the data-selection tier, not an accuracy tolerance.
/// Changing the tier clears an explicit JPL filename so tier-based selection
/// applies, but preserves an explicit LEB path. No file opens here; actual
/// per-body coverage still applies. Raises ValueError for an unknown tier.
#[pyfunction]
fn set_precision_tier(py: Python<'_>, tier: &str) -> PyResult<()> {
    runtime()
        .set_precision_tier(tier)
        .map_err(|error| crate::errors::to_pyerr(py, error))
}

/// Select an explicit local LEB file, or clear that selection with None.
///
/// filepath is a string/path-like filename or None. Opening and validation
/// occur when data is consumed, not here; existing retained handles remain
/// valid. None allows ordinary configured discovery, rather than closing
/// every handle or disabling LEB calculation. This does not change calc mode.
#[pyfunction(signature = (filepath))]
fn set_leb_file(filepath: Option<PathBuf>) {
    runtime().set_leb_file(filepath);
}

/// Override whether observed IERS Delta-T data participates in model clocks.
///
/// enabled is True or False; None clears the explicit override and restores
/// the captured configuration/default policy. Data acquisition follows the
/// ordinary IERS/network settings. This does not replace independent bundled
/// observation clocks used by APIs that document their own time conversion.
#[pyfunction(signature = (enabled))]
fn set_iers_delta_t_enabled(enabled: Option<bool>) {
    runtime().set_iers_delta_t_enabled(enabled);
}

/// Override TT minus UT1 with a constant number of days.
///
/// acc is a float in days, despite its historical name. None or
/// DELTAT_AUTOMATIC clears the override and resumes computed Delta-T.
/// Other values are stored as supplied; this setter does not validate them.
/// The override affects configured model clocks, not every independent
/// observation clock. See deltat and set_delta_t_model.
#[pyfunction(signature = (acc))]
fn set_delta_t_userdef(acc: Option<f64>) {
    runtime().set_delta_t_userdef(acc);
}

/// Select the default runtime's analytical Delta-T model.
///
/// model is 'smh2016' or 'espenak_meeus', case/whitespace-insensitive.
/// None clears the explicit selection, allowing captured configuration and
/// then the 'smh2016' default. User Delta-T overrides and enabled observed
/// IERS data take precedence over computed model values. Raises ValueError
/// for an unknown name without replacing the previous selection.
#[pyfunction(signature = (model))]
fn set_delta_t_model(py: Python<'_>, model: Option<&str>) -> PyResult<()> {
    runtime()
        .set_delta_t_model(model)
        .map_err(|error| crate::errors::to_pyerr(py, error))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(calc_ut, module)?)?;
    module.add_function(wrap_pyfunction!(calc, module)?)?;
    module.add_function(wrap_pyfunction!(set_calc_mode, module)?)?;
    module.add_function(wrap_pyfunction!(set_precision_tier, module)?)?;
    module.add_function(wrap_pyfunction!(set_leb_file, module)?)?;
    module.add_function(wrap_pyfunction!(set_iers_delta_t_enabled, module)?)?;
    module.add_function(wrap_pyfunction!(set_delta_t_userdef, module)?)?;
    module.add_function(wrap_pyfunction!(set_delta_t_model, module)?)?;
    Ok(())
}
