// SPDX-License-Identifier: AGPL-3.0
//! Convert native constant values and forward the four SPK identity lookups.

use libephemeris_core::constants as core;
use libephemeris_core::public_constants::{PUBLIC_CONSTANTS, PublicConstantValue};
use libephemeris_core::public_module_constants::{ModuleConstantValue, PUBLIC_MODULE_CONSTANTS};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PySet, PyTuple};

/// Look up the built-in JPL Horizons target string for a library body ID.
///
/// ipl is an integer libephemeris identifier, not a Horizons command.
/// Returns str, or None when the built-in SPK identity map has no entry.
/// This lookup performs no network access or kernel registration.
#[pyfunction]
fn get_horizons_id(py: Python<'_>, ipl: i32) -> Option<&'static str> {
    py.detach(|| core::get_horizons_id(ipl))
}

/// Look up the built-in NAIF/SPICE target ID for a library body ID.
///
/// ipl is an integer libephemeris identifier. Returns int, or None when
/// unmapped; it does not inspect user registrations or acquire an SPK file.
#[pyfunction]
fn get_naif_id_from_ipl(py: Python<'_>, ipl: i32) -> Option<i32> {
    py.detach(|| core::get_naif_id_from_ipl(ipl))
}

/// Return a body's built-in (Horizons target string, NAIF integer ID) pair.
///
/// ipl is a libephemeris body ID. Returns None if absent from the static map.
/// Unlike get_spk_body_info, this does not query registered filenames or
/// imply that a kernel is installed or covers any particular epoch.
#[pyfunction]
fn get_spk_body_info_from_map(py: Python<'_>, ipl: i32) -> Option<(&'static str, i32)> {
    py.detach(|| core::get_spk_body_info_from_map(ipl))
}

/// Test whether a library body ID is in the built-in SPK-download blocklist.
///
/// ipl is an integer libephemeris identifier; the result is bool. False
/// does not guarantee target coverage, network permission or a successful
/// Horizons download. No network request is made by this lookup.
#[pyfunction]
fn is_spk_auto_download_blocked(py: Python<'_>, ipl: i32) -> bool {
    py.detach(|| core::is_spk_auto_download_blocked(ipl))
}

fn module_constant(py: Python<'_>, value: ModuleConstantValue) -> PyResult<Py<PyAny>> {
    Ok(match value {
        ModuleConstantValue::None => py.None(),
        ModuleConstantValue::Bool(value) => value.into_pyobject(py)?.to_owned().into_any().unbind(),
        ModuleConstantValue::Integer(value) => value.into_pyobject(py)?.into_any().unbind(),
        ModuleConstantValue::Float(value) => value.into_pyobject(py)?.into_any().unbind(),
        ModuleConstantValue::Text(value) => value.into_pyobject(py)?.into_any().unbind(),
        ModuleConstantValue::Tuple(values) => {
            let values = values
                .iter()
                .map(|&value| module_constant(py, value))
                .collect::<PyResult<Vec<_>>>()?;
            PyTuple::new(py, values)?.into_any().unbind()
        }
        ModuleConstantValue::List(values) => {
            let values = values
                .iter()
                .map(|&value| module_constant(py, value))
                .collect::<PyResult<Vec<_>>>()?;
            PyList::new(py, values)?.into_any().unbind()
        }
        ModuleConstantValue::Set(values) => {
            let values = values
                .iter()
                .map(|&value| module_constant(py, value))
                .collect::<PyResult<Vec<_>>>()?;
            PySet::new(py, values)?.into_any().unbind()
        }
        ModuleConstantValue::Map(values) => {
            let dictionary = PyDict::new(py);
            for &(key, value) in values {
                dictionary.set_item(module_constant(py, key)?, module_constant(py, value)?)?;
            }
            dictionary.into_any().unbind()
        }
    })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    // Package identity belongs to the Python distribution, not a generated
    // scientific-constant snapshot or the independently versioned Rust core.
    module.add("version", crate::identity::package_version())?;
    module.add(
        "_module_state_DEFAULT_DATA_DIR",
        crate::positions::runtime()
            .default_data_directory()
            .map_err(|error| crate::errors::to_pyerr(py, error))?
            .to_string_lossy()
            .into_owned(),
    )?;
    for &(key, _, _, value) in PUBLIC_MODULE_CONSTANTS {
        module.add(key, module_constant(py, value)?)?;
    }
    for &(name, value) in PUBLIC_CONSTANTS {
        match value {
            PublicConstantValue::Integer(value) => module.add(name, value)?,
            PublicConstantValue::Float(value) => module.add(name, value)?,
            PublicConstantValue::Text(value) => module.add(name, value)?,
            PublicConstantValue::IntegerSet(values) => {
                module.add(name, PySet::new(py, values.iter().copied())?)?;
            }
            PublicConstantValue::IntegerMap(values) => {
                let dictionary = PyDict::new(py);
                for &(key, value) in values {
                    dictionary.set_item(key, value)?;
                }
                module.add(name, dictionary)?;
            }
            PublicConstantValue::IntegerTextMap(values) => {
                let dictionary = PyDict::new(py);
                for &(key, value) in values {
                    dictionary.set_item(key, value)?;
                }
                module.add(name, dictionary)?;
            }
            PublicConstantValue::IntegerFloatMap(values) => {
                let dictionary = PyDict::new(py);
                for &(key, value) in values {
                    dictionary.set_item(key, value)?;
                }
                module.add(name, dictionary)?;
            }
            PublicConstantValue::SpkIdentityMap(values) => {
                let dictionary = PyDict::new(py);
                for &(body, horizons, naif) in values {
                    dictionary.set_item(body, (horizons, naif))?;
                }
                module.add(name, dictionary)?;
            }
            PublicConstantValue::MajorAsteroidMap(values) => {
                let dictionary = PyDict::new(py);
                for &(body, number, horizons, naif, display_name) in values {
                    dictionary.set_item(body, (number, horizons, naif, display_name))?;
                }
                module.add(name, dictionary)?;
            }
        }
    }
    module.add_function(wrap_pyfunction!(get_horizons_id, module)?)?;
    module.add_function(wrap_pyfunction!(get_naif_id_from_ipl, module)?)?;
    module.add_function(wrap_pyfunction!(get_spk_body_info_from_map, module)?)?;
    module.add_function(wrap_pyfunction!(is_spk_auto_download_blocked, module)?)?;
    Ok(())
}
