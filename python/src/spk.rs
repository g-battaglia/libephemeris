//! Native immutable SPK registration and policy-gated acquisition.
use crate::{errors::to_pyerr, positions::runtime};
use pyo3::exceptions::{PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use std::collections::HashMap;
use std::path::PathBuf;

type Position = (f64, f64, f64, f64, f64, f64);

#[derive(FromPyObject)]
enum NaifId {
    Integer(i32),
    Text(String),
}
impl NaifId {
    fn number(self) -> PyResult<i32> {
        match self {
            Self::Integer(value) => Ok(value),
            Self::Text(value) => value
                .trim()
                .parse()
                .map_err(|error: std::num::ParseIntError| PyValueError::new_err(error.to_string())),
        }
    }
}
/// Acquire a Horizons-generated SPK and return its filename as str.
///
/// body is a Horizons command (for example "1;"), not a library body ID;
/// start/end are calendar strings interpreted by Horizons, typically
/// "YYYY-MM-DD". center is its location code (default solar-system barycenter
/// "500@0"). path is a str/PathLike output directory, not a filename; None
/// selects the configured SPK cache or data-root/spk.
///
/// overwrite=False keeps an existing structurally valid file without checking
/// requested target/epoch coverage; True requests replacement. Network misses
/// create directories and atomically publish validated bytes. timeout is a
/// positive integer number of seconds per transport attempt (not an overall
/// retry deadline). Policy, HTTP, payload and filesystem errors propagate.
/// No body registration or extrapolation is performed.
#[pyfunction(signature=(body,start,end,path=None,center="500@0",overwrite=false,timeout=120))]
#[allow(clippy::too_many_arguments)]
fn download_spk(
    py: Python<'_>,
    body: &str,
    start: &str,
    end: &str,
    path: Option<PathBuf>,
    center: &str,
    overwrite: bool,
    timeout: i64,
) -> PyResult<String> {
    py.detach(|| {
        runtime().download_spk(
            body,
            start,
            end,
            path.as_deref(),
            center,
            overwrite,
            timeout as f64,
        )
    })
    .map(|path| path.to_string_lossy().into_owned())
    .map_err(|error| to_pyerr(py, error))
}
/// Bind library body ID ipl to NAIF target naif_id in a local SPK image.
///
/// spk_file accepts str or os.PathLike; relative names search the configured
/// SPK cache, data root, then the supplied relative path. naif_id accepts an
/// integer or decimal integer string. These two ID systems are not
/// interchangeable. Re-registering ipl replaces its mapping for future calls.
///
/// The kernel is opened eagerly and retained as an immutable owned image;
/// no download occurs. Return None. Invalid integer text raises ValueError;
/// missing files, absent targets and malformed/unsupported SPKs raise native
/// SPK errors. Registration alone does not prove a complete center chain or
/// usable coverage at a particular epoch.
#[pyfunction]
fn register_spk_body(py: Python<'_>, ipl: i32, spk_file: PathBuf, naif_id: NaifId) -> PyResult<()> {
    let target = naif_id.number()?;
    py.detach(|| {
        let path = runtime().resolve_spk_path(&spk_file)?;
        runtime().register_spk_body(ipl, &path, target)
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Remove the default mapping for library body ID ipl; return None.
///
/// Unknown IDs are a no-op. This does not delete the kernel on disk or
/// invalidate immutable images held by active calculations/other snapshots.
#[pyfunction]
fn unregister_spk_body(py: Python<'_>, ipl: i32) {
    py.detach(|| runtime().unregister_spk_body(ipl));
}
/// Return (kernel_path_str, naif_target_id) for library body ID ipl, or None.
///
/// This is copied registration metadata, not an epoch-coverage check; it
/// remains available for mappings whose images await lazy reopening. No I/O.
#[pyfunction]
fn get_spk_body_info(py: Python<'_>, ipl: i32) -> Option<(String, i32)> {
    py.detach(|| runtime().registered_spks().remove(&ipl))
}
/// Return a copied {library_body_id: (kernel_path_str, naif_target_id)} dict.
///
/// Mutating the result does not change registrations. Metadata includes
/// lazily reopenable definitions and performs no file or network access.
#[pyfunction]
fn list_spk_bodies(py: Python<'_>) -> HashMap<i32, (String, i32)> {
    py.detach(|| runtime().registered_spks())
}
/// Return the whole kernel's (start_jd, end_jd) TDB envelope, or None.
///
/// spk_file is str or os.PathLike, resolved like register_spk_body(). Missing
/// or unreadable/corrupt files return None; configuration errors can raise.
/// The local kernel is opened for inspection, without registration/download.
/// This union envelope is not proof of continuous coverage, a particular
/// NAIF target's coverage or a complete target-to-center chain. Calculations
/// still check their actual segments and never extrapolate beyond them.
#[pyfunction]
fn get_spk_coverage(py: Python<'_>, spk_file: PathBuf) -> PyResult<Option<(f64, f64)>> {
    py.detach(|| runtime().spk_coverage(&spk_file))
        .map_err(|error| to_pyerr(py, error))
}
/// Download a Horizons SPK, register ipl, and return the kernel path as str.
///
/// body/start/end/path/center/overwrite/timeout follow download_spk(); ipl is
/// the library body ID, while naif_id is the kernel's integer NAIF target.
/// None infers the target from the downloaded file and body command. A valid
/// cached file is retained only if its whole-file envelope passes the native
/// requested-window check; an insufficient envelope triggers replacement even
/// with overwrite=False. That check is not target-chain proof.
///
/// Download/validation/deduction/registration errors propagate. Completed
/// downloads may remain on disk after a later registration failure; no new
/// body mapping is published on that failure. Existing calculation snapshots
/// retain their images when a mapping/file is replaced. No extrapolation.
#[pyfunction(signature=(body,ipl,start,end,naif_id=None,path=None,center="500@0",overwrite=false,timeout=120))]
#[allow(clippy::too_many_arguments)]
fn download_and_register_spk(
    py: Python<'_>,
    body: &str,
    ipl: i32,
    start: &str,
    end: &str,
    naif_id: Option<i32>,
    path: Option<PathBuf>,
    center: &str,
    overwrite: bool,
    timeout: i64,
) -> PyResult<String> {
    py.detach(|| {
        runtime().download_and_register_spk(
            body,
            ipl,
            start,
            end,
            naif_id,
            path.as_deref(),
            center,
            overwrite,
            timeout as f64,
        )
    })
    .map(|path| path.to_string_lossy().into_owned())
    .map_err(|error| to_pyerr(py, error))
}
/// Evaluate an explicitly registered library body ipl at TT Julian Day t.
///
/// Return None for an unknown/unservable registration, otherwise a six-tuple
/// (longitude, latitude, distance, longitude_rate, latitude_rate, radial_rate)
/// in degrees, AU, degrees/TT-day and AU/TT-day. Rates require FLG_SPEED.
/// This lower-level helper is not calc(): type 2/3 outputs are ecliptic-of-date;
/// type 21 also supports J2000. EQUATORIAL/XYZ/RADIANS/BARYCTR do not provide
/// their usual public calc representations here; type 21 ignores TOPOCTR.
///
/// Actual target/center-chain coverage is checked without extrapolation or a
/// Keplerian/LEB fallback. A usable target may acquire the main observer kernel
/// (including permitted network I/O). Closed/malformed/out-of-range SPKs and
/// observer/frame/clock errors propagate; forced LEB mode raises RuntimeError
/// if main-kernel access is needed. TT is converted on the native SPK route.
#[pyfunction]
fn calc_spk_body_position(
    py: Python<'_>,
    t: f64,
    ipl: i32,
    iflag: i32,
) -> PyResult<Option<Position>> {
    py.detach(|| runtime().calc_spk_body_position(t, ipl, iflag))
        .map(|value| value.map(|[a, b, c, d, e, f]| (a, b, c, d, e, f)))
        .map_err(|error| match error {
            libephemeris_core::errors::Error::ConfigurationMissing {
                missing_config: "jpl_source_access",
                ..
            } => PyRuntimeError::new_err(error.to_string()),
            _ => to_pyerr(py, error),
        })
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    horizons::register(module)?;
    module.add_function(wrap_pyfunction!(calc_spk_body_position, module)?)?;
    module.add_function(wrap_pyfunction!(download_spk, module)?)?;
    module.add_function(wrap_pyfunction!(register_spk_body, module)?)?;
    module.add_function(wrap_pyfunction!(unregister_spk_body, module)?)?;
    module.add_function(wrap_pyfunction!(get_spk_body_info, module)?)?;
    module.add_function(wrap_pyfunction!(list_spk_bodies, module)?)?;
    module.add_function(wrap_pyfunction!(get_spk_coverage, module)?)?;
    module.add_function(wrap_pyfunction!(download_and_register_spk, module)?)?;
    Ok(())
}

#[path = "spk/horizons.rs"]
mod horizons;
