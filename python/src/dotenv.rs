//! Native configuration-file acquisition without process-environment writes.
use crate::positions::runtime;
use pyo3::prelude::*;
use std::path::PathBuf;
/// Read a dotenv file into the native runtime's captured configuration.
///
/// path accepts str or os.PathLike. None uses a nonempty LIBEPHEMERIS_ENV_FILE
/// exclusively if set; otherwise it selects the first existing .env in the
/// working directory, then ~/.libephemeris/.env. Missing/unreadable selected
/// files do not fall through to another location. Only LIBEPHEMERIS_ keys are
/// accepted; no variable interpolation occurs. Captured environment values win;
/// override=True permits replacing them, but explicit API settings still win.
/// No os.environ entries are changed and no network request is made.
///
/// Return False for a missing/unreadable UTF-8 file, True once read (even if
/// no usable keys occur). A successful load retires default reader discovery
/// for future calls; already retained calculation/data snapshots stay valid.
#[pyfunction(signature=(path=None,*,r#override=false))]
fn load_dotenv(py: Python<'_>, path: Option<PathBuf>, r#override: bool) -> bool {
    py.detach(|| runtime().load_dotenv(path.as_deref(), r#override))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(load_dotenv, module)?)?;
    Ok(())
}
