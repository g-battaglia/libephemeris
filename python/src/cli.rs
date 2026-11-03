//! Console entry point: convert Python argv and call the native application.
use pyo3::prelude::*;

/// Run the native command-line application using this interpreter's sys.argv.
///
/// Parameters
/// ----------
/// None
///     Arguments come from ``sys.argv[1:]``; do not pass an argv list here.
///     The optional generation executable is resolved beside sys.executable,
///     not from another virtual environment's PATH.
///
/// Returns
/// -------
/// int
///     Native command exit status: 0 for success, 1 for operational failure,
///     2 for invalid arguments. Offline LEB generation preserves its child
///     application's normal exit code; a missing executable or abnormal child
///     termination gives 1. Returns the code rather than calling sys.exit.
///
/// Notes
/// -----
/// Commands may print to stdout/stderr, prompt on stdin, read configuration
/// and data, and, when explicitly selected, download, install, generate or
/// replace files. Status inspection can trigger lazy IERS acquisition and
/// directory creation when configuration allows it; it is not universally
/// side-effect-free. Native work releases the GIL but remains synchronous. Python import/attribute/type
/// extraction errors propagate before dispatch; command failures are reported
/// by the returned status and native diagnostics.
#[pyfunction]
fn cli_main(py: Python<'_>) -> PyResult<i32> {
    let system = py.import("sys")?;
    let arguments: Vec<String> = system.getattr("argv")?.extract()?;
    let executable: std::path::PathBuf = system.getattr("executable")?.extract()?;
    // Absolute console invocations need not activate their virtualenv's PATH.
    // Select the helper beside this interpreter, with no cross-env fallback.
    let generator = executable
        .parent()
        .unwrap_or(std::path::Path::new(""))
        .join(if cfg!(windows) {
            "libephemeris-generate.exe"
        } else {
            "libephemeris-generate"
        });
    Ok(py.detach(|| {
        libephemeris_core::cli::run_with_generator(arguments.iter().skip(1), Some(&generator))
    }))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(cli_main, module)?)
}
