// SPDX-License-Identifier: AGPL-3.0
//! Standard Python logging adapters backed by native threshold configuration.

use crate::positions::runtime;
use libephemeris_core::logging;
use pyo3::prelude::*;
use pyo3::sync::PyOnceLock;
use pyo3::types::{PyCFunction, PyDict, PyTuple};
use std::ffi::CString;
use std::sync::OnceLock;

static LOGGER: PyOnceLock<Py<PyAny>> = PyOnceLock::new();
static BRIDGE: OnceLock<Option<pyo3_log::ResetHandle>> = OnceLock::new();
static ENABLE_LOGGING_DOC: OnceLock<CString> = OnceLock::new();

fn default_level() -> i32 {
    runtime().sources().log_level()
}

fn logger(py: Python<'_>) -> PyResult<Bound<'_, PyAny>> {
    let value = LOGGER.get_or_try_init(py, || {
        let standard = py.import("logging")?;
        let logger = standard
            .getattr("getLogger")?
            .call1((logging::LOGGER_NAME,))?;
        let level = default_level();
        logging::set_threshold(level);
        logger.call_method1("setLevel", (level,))?;
        let stream = py.import("sys")?.getattr("stderr")?;
        let handler = standard.getattr("StreamHandler")?.call1((stream,))?;
        handler.call_method1("setLevel", (level,))?;
        let formatter = standard
            .getattr("Formatter")?
            .call1(("[libephemeris] %(levelname)s: %(message)s",))?;
        handler.call_method1("setFormatter", (formatter,))?;
        logger.call_method1("addHandler", (handler,))?;
        logger.setattr("propagate", false)?;
        Ok::<_, PyErr>(logger.unbind())
    })?;
    Ok(value.bind(py).clone())
}

/// Return the shared standard-library logging.Logger named "libephemeris".
///
/// The extension initializes it with a stderr StreamHandler, the configured
/// threshold and propagation disabled. This is the live logger, not a copy;
/// Python logging configuration and handler changes affect later messages.
/// Python logging/import errors propagate. No astronomical data is acquired.
#[pyfunction]
fn get_logger(py: Python<'_>) -> PyResult<Bound<'_, PyAny>> {
    logger(py)
}

/// Set the native threshold, Python logger and all its current handler levels.
///
/// level is an integer using Python logging values (DEBUG=10, INFO=20,
/// WARNING=30, ERROR=40, CRITICAL=50). Arbitrary integers are passed through;
/// values above 50 suppress all named native severities. Return None.
/// Errors raised by the Python logger or handlers propagate.
#[pyfunction]
fn set_log_level(py: Python<'_>, level: i32) -> PyResult<()> {
    let value = logger(py)?;
    value.call_method1("setLevel", (level,))?;
    for handler in value.getattr("handlers")?.try_iter()? {
        handler?.call_method1("setLevel", (level,))?;
    }
    logging::set_threshold(level);
    Ok(())
}

/// Suppress named library severities by setting logger/native thresholds to 51.
///
/// Return None; handlers are retained and their levels are not changed.
/// This does not call logging.disable() or disable other applications' loggers.
/// Use enable_logging(level) to restore a chosen threshold. Logger errors
/// propagate and numerical calculations are unaffected.
#[pyfunction]
fn disable_logging(py: Python<'_>) -> PyResult<()> {
    logger(py)?.call_method1("setLevel", (logging::DISABLED_THRESHOLD,))?;
    logging::disable_logging();
    Ok(())
}

/// Enable logging at an integer Python severity; see the exported docstring
/// for the dynamically resolved default and handler effects.
#[pyfunction(signature = (level=default_level()))]
fn enable_logging(py: Python<'_>, level: i32) -> PyResult<()> {
    set_log_level(py, level)
}

/// Format size_bytes with one decimal digit and binary 1024-based units.
///
/// Return a str labeled B, KB, MB, GB, TB or PB, for example
/// format_file_size(1536) == "1.5 KB". Fractional and negative values are
/// accepted, not validated as actual file sizes. NaN/infinities use the native
/// floating-point spelling. No filesystem access or logging occurs.
#[pyfunction]
fn format_file_size(py: Python<'_>, size_bytes: f64) -> String {
    py.detach(|| logging::format_file_size_f64(size_bytes))
}

/// Install the native log bridge and register the five public logging helpers.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    logger(py)?;
    // Logger-only caching keeps direct stdlib Logger.setLevel changes live.
    // Another embedding application may already own Rust's global logger;
    // preserve it rather than replacing the application's sink.
    if BRIDGE.get().is_none() {
        let bridge = pyo3_log::Logger::new(py, pyo3_log::Caching::Loggers)?;
        let _ = BRIDGE.set(bridge.install().ok());
    }
    module.add("_DEFAULT_LOG_LEVEL", default_level())?;
    module.add("_LOGGER_NAME", logging::LOGGER_NAME)?;
    module.add_function(wrap_pyfunction!(get_logger, module)?)?;
    module.add_function(wrap_pyfunction!(set_log_level, module)?)?;
    module.add_function(wrap_pyfunction!(disable_logging, module)?)?;
    let typed_enable = wrap_pyfunction!(enable_logging, module)?.unbind();
    let signature = ENABLE_LOGGING_DOC.get_or_init(|| {
        CString::new(format!(
            "enable_logging(level={})\n--\n\nEnable library logging at an integer Python logging severity.\n\nSet the native threshold, shared libephemeris logger and every current\nhandler to level; return None. DEBUG=10, INFO=20, WARNING=30, ERROR=40\nand CRITICAL=50. Arbitrary integer thresholds are passed through.\n\nWhen omitted, level comes from the runtime's captured environment/config\n(default WARNING), not the previously active level before disable_logging().\nPython logger/handler errors propagate. Example: enable_logging(20).",
            default_level()
        ))
        .expect("integer signature and constant documentation contain no NUL")
    });
    let enable = PyCFunction::new_closure(
        py,
        Some(c"enable_logging"),
        Some(signature.as_c_str()),
        move |args: &Bound<'_, PyTuple>, kwargs: Option<&Bound<'_, PyDict>>| {
            typed_enable
                .bind(args.py())
                .call(args, kwargs)
                .map(Bound::unbind)
        },
    )?;
    enable.setattr("__module__", module.name()?)?;
    module.add_function(enable)?;
    module.add_function(wrap_pyfunction!(format_file_size, module)?)?;
    Ok(())
}
