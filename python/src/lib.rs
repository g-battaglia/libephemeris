// SPDX-License-Identifier: AGPL-3.0
//! Minimal Python interface to the LibEphemeris Rust engine.

use pyo3::prelude::*;
mod arabic_parts;
mod cli;
mod constants;
mod context;
mod contrib;
mod crossing;
mod dotenv;
mod download;
mod eclipse;
mod eclipse_visibility;
mod errors;
mod extinction;
mod fixed_stars;
mod heliacal;
mod heliacal_model;
mod houses;
mod hypothetical;
mod identity;
mod iers_data;
mod inventory;
mod logging_config;
mod lunar;
mod minor_bodies;
mod net;
mod numeric_frames;
mod physical_utilities;
mod planetary_moons;
mod planets;
mod positions;
mod rebound_integration;
mod spk;
mod spk_auto;
mod spk_auto_extra;
mod state;
mod time_utils;
mod tracing;
mod utils;

/// Astronomical calculations, time scales and data management in native Rust.
///
/// Python entry points are re-exported by libephemeris and its public modules.
/// calc takes TT Julian days; calc_ut takes UT1 Julian days. Neither accepts
/// a civil UTC label directly: utc_to_jd returns (jd_tt, jd_ut1).
/// Ordinary position results are (six_values, returned_flags), with spherical
/// degrees/AU and rates per TT day unless output flags select another layout.
///
/// Global setters configure one process-wide runtime. EphemerisContext provides
/// independently configurable observers and sidereal settings. Configure data
/// and network policy before concurrent work; calculations release the GIL but
/// remain synchronous. Coverage and missing resources are reported as errors.
/// See individual functions for flags, units, result layouts and exceptions.
#[pymodule]
fn _core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    libephemeris_core::identity::set_distribution_version(identity::package_version())
        .map_err(pyo3::exceptions::PyRuntimeError::new_err)?;
    for warning in positions::runtime().sources().warnings() {
        let message = std::ffi::CString::new(warning.to_string())
            .map_err(|error| pyo3::exceptions::PyValueError::new_err(error.to_string()))?;
        PyErr::warn(
            py,
            &py.get_type::<pyo3::exceptions::PyUserWarning>(),
            &message,
            1,
        )?;
    }
    module.add("__version__", identity::package_version())?;
    errors::register(module)?;
    constants::register(module)?;
    positions::register(module)?;
    cli::register(module)?;
    planets::register(module)?;
    utils::register(module)?;
    time_utils::register(module)?;
    state::register(module)?;
    context::register(module)?;
    houses::register(module)?;
    lunar::register(module)?;
    crossing::register(module)?;
    extinction::register(module)?;
    heliacal::register(module)?;
    heliacal_model::register(module)?;
    fixed_stars::register(module)?;
    planetary_moons::register(module)?;
    eclipse::register(module)?;
    eclipse_visibility::register(module)?;
    rebound_integration::register(module)?;
    logging_config::register(module)?;
    tracing::register(module)?;
    spk::register(module)?;
    spk_auto::register(module)?;
    download::register(module)?;
    net::register(module)?;
    iers_data::register(module)?;
    inventory::register(module)?;
    dotenv::register(module)?;
    hypothetical::register(module)?;
    minor_bodies::register(module)?;
    physical_utilities::register(module)?;
    numeric_frames::register(module)?;
    arabic_parts::register(module)?;
    let modules = py.import("sys")?.getattr("modules")?;
    let convenience = PyModule::new(py, "libephemeris._core._contrib")?;
    contrib::register(&convenience)?;
    modules.set_item("libephemeris._core._contrib", &convenience)?;
    module.add("_contrib", &convenience)?;
    module.add("contrib", &convenience)?;
    let automatic = PyModule::new(py, "libephemeris._core._spk_auto")?;
    spk_auto_extra::register(&automatic)?;
    modules.set_item("libephemeris._core._spk_auto", &automatic)?;
    module.add("_spk_auto", &automatic)?;
    module.add("spk_auto", &automatic)?;
    let names: Vec<String> = module
        .dict()
        .keys()
        .iter()
        .filter_map(|key| key.extract::<String>().ok())
        .filter(|name| !name.starts_with('_') || name == "__version__")
        .collect();
    module.add("__all__", names)?;
    Ok(())
}
