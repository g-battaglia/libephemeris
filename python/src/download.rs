//! Argument/result conversion for native data inspection and installers.
use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::config::download_surface::DownloadValidator;
use libephemeris_core::errors::Error;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::path::{Path, PathBuf};

fn path_object<'py>(py: Python<'py>, path: &Path) -> PyResult<Bound<'py, PyAny>> {
    py.import("pathlib")?.getattr("Path")?.call1((path,))
}
fn paths_object<'py>(py: Python<'py>, paths: Vec<PathBuf>) -> PyResult<Bound<'py, PyList>> {
    let result = PyList::empty(py);
    for path in paths {
        result.append(path_object(py, &path)?)?;
    }
    Ok(result)
}
fn download_error(py: Python<'_>, error: Error) -> PyErr {
    match error {
        Error::HashMismatch { .. }
        | Error::AssetRejected {
            reason: "format validation failed" | "unknown native format validator",
            ..
        } => PyValueError::new_err(error.to_string()),
        _ => to_pyerr(py, error),
    }
}
/// Return the configured absolute data root as pathlib.Path without creating it.
///
/// This follows resolved data-directory configuration, not necessarily the
/// explicit main-kernel search path set by set_ephe_path(). Configuration or
/// absolute-path resolution errors propagate; no download/file open occurs.
#[pyfunction]
fn get_data_dir(py: Python<'_>) -> PyResult<Bound<'_, PyAny>> {
    let path = py
        .detach(|| runtime().data_directory())
        .map_err(|error| to_pyerr(py, error))?;
    path_object(py, &path)
}
/// Create the data root's leb subdirectory and return it as pathlib.Path.
///
/// This getter has a filesystem side effect, but acquires no astronomical
/// files and changes no active selection. Configuration/filesystem errors raise.
#[pyfunction]
fn get_leb_dir(py: Python<'_>) -> PyResult<Bound<'_, PyAny>> {
    let path = py
        .detach(|| runtime().leb_directory())
        .map_err(|error| to_pyerr(py, error))?;
    path_object(py, &path)
}
/// Return data-root/leb/ephemeris_{tier_name}.leb as pathlib.Path.
///
/// This historical LEB1 naming helper preserves tier_name literally: it does
/// not validate a tier, locate a LEB2 asset or test file existence. It creates
/// the leb directory, but no file, and does not activate the returned path.
/// Configuration/filesystem errors propagate.
#[pyfunction(signature=(tier_name))]
fn get_leb_path_for_tier<'py>(py: Python<'py>, tier_name: &str) -> PyResult<Bound<'py, PyAny>> {
    let path = py
        .detach(|| runtime().leb_path_for_tier(tier_name))
        .map_err(|error| to_pyerr(py, error))?;
    path_object(py, &path)
}
/// Rename temp_path over dest_path, consuming the temporary file on success.
///
/// Both arguments accept str or os.PathLike. On Unix, setting permissions to
/// 0644 is attempted first (failure is ignored); rename failures raise OSError.
/// No parent directory is created and no digest/format check is performed.
/// Atomic replacement requires filesystem support, normally same-filesystem
/// paths. This is a publication primitive, not a downloader; return None.
#[pyfunction(signature=(temp_path,dest_path))]
fn publish_temp_file(py: Python<'_>, temp_path: PathBuf, dest_path: PathBuf) -> PyResult<()> {
    py.detach(|| libephemeris_core::download::publish_temp_file(&temp_path, &dest_path))
        .map_err(Into::into)
}
/// Download url to a str/PathLike dest_path and return True on success.
///
/// Always attempt acquisition, rather than treating an existing file as cached.
/// description labels diagnostics; expected_sha256=None omits digest checking.
/// validator is None, "bsp" or "leb" (a native structural check), not a Python
/// callable. show_progress is accepted but no Python progress callback runs.
/// timeout is positive integer seconds for connect and overall exchange
/// budgets (default 30); there is no internal retry.
///
/// Parent directories and an exclusive temporary sibling are created before
/// policy/timeout checks. A failed attempt cleans its own temporary file and
/// preserves an existing destination; success atomically replaces it. Unknown
/// validators, digest mismatches and format failures raise ValueError; sealed
/// policy, transport and filesystem errors propagate. Native work releases the
/// GIL; Python code is not invoked as a validator or progress handler.
#[pyfunction(signature=(url,dest_path,description="Downloading",expected_sha256=None,show_progress=true,timeout=30,validator=None))]
#[allow(clippy::too_many_arguments)]
fn download_file(
    py: Python<'_>,
    url: &str,
    dest_path: PathBuf,
    description: &str,
    expected_sha256: Option<&str>,
    show_progress: bool,
    timeout: i64,
    validator: Option<&str>,
) -> PyResult<bool> {
    let _ = show_progress;
    let validator = validator
        .map(DownloadValidator::from_name)
        .transpose()
        .map_err(|error| download_error(py, error))?;
    py.detach(|| {
        runtime().download_file(
            url,
            &dest_path,
            description,
            expected_sha256,
            timeout as f64,
            validator,
        )
    })
    .map(|_| true)
    .map_err(|error| download_error(py, error))
}
/// Install the pinned legacy planet_centers.bsp; return its pathlib.Path.
///
/// Keep a digest/structure-verified cached copy unless force=True requests
/// replacement. Missing/replacement acquisition is policy-gated and atomic;
/// directories may be created. Retire default optional-center discovery for
/// future use, while existing immutable center handles remain valid.
/// show_progress/quiet are accepted but do not control Python output/callbacks.
/// Digest/format failures raise ValueError; other acquisition errors propagate.
/// This is the legacy auxiliary file, not all tier-specific/main kernels.
#[pyfunction(signature=(force=false,show_progress=true,quiet=false))]
fn download_planet_centers(
    py: Python<'_>,
    force: bool,
    show_progress: bool,
    quiet: bool,
) -> PyResult<Bound<'_, PyAny>> {
    let _ = (show_progress, quiet);
    let path = py
        .detach(|| runtime().download_planet_centers(force))
        .map_err(|error| download_error(py, error))?;
    path_object(py, &path)
}
/// Install the fixed auxiliary-data set and return a list of pathlib.Path.
///
/// Currently this contains only legacy planet_centers.bsp: it does not mean
/// all tiers, LEB files, minor-body SPKs or ASSIST data. force=False retains
/// verified files; True requests atomic replacement under network policy.
/// show_progress/quiet are accepted but ignored. Directory, policy, transport
/// and validation errors follow download_planet_centers().
#[pyfunction(signature=(force=false,show_progress=true,quiet=false))]
fn download_all(
    py: Python<'_>,
    force: bool,
    show_progress: bool,
    quiet: bool,
) -> PyResult<Bound<'_, PyList>> {
    let _ = (show_progress, quiet);
    let paths = py
        .detach(|| runtime().download_all(force))
        .map_err(|error| download_error(py, error))?;
    paths_object(py, paths)
}
/// Return copied local manifest-asset metadata without downloads or validation.
///
/// Map each filename to a fresh dict with exists (bool), path (pathlib.Path
/// or None), size (bytes or None), description (str) and expected_size_mb
/// (advertised MiB). Missing files have neither path nor size. Presence/size
/// is not proof of a correct digest or scientifically usable coverage. No
/// directory is created; configuration/metadata errors other than missing
/// paths propagate. Editing the result does not modify installation state.
#[pyfunction]
fn check_data_status(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    let status = py
        .detach(|| runtime().check_data_status())
        .map_err(|error| to_pyerr(py, error))?;
    let result = PyDict::new(py);
    for (name, status) in status {
        let row = PyDict::new(py);
        row.set_item("exists", status.exists)?;
        row.set_item(
            "path",
            status
                .path
                .as_deref()
                .map(|path| path_object(py, path))
                .transpose()?,
        )?;
        row.set_item("size", status.size)?;
        row.set_item("description", status.description)?;
        row.set_item("expected_size_mb", status.expected_size_mb)?;
        result.set_item(name, row)?;
    }
    Ok(result)
}
/// Print the resolved data inventory to native stdout; return None.
///
/// as_json=True selects pretty JSON; otherwise verbose>0 adds paths and >1
/// lists concrete DE/LEB/ASSIST/SPK candidates. Presence/size is descriptive,
/// not digest or numerical validation. No ephemeris is downloaded by the
/// inventory itself, but its IERS range query can create/load/acquire the
/// selected IERS cache under automatic-download/network policy. Filesystem,
/// configuration and IERS acquisition errors propagate.
#[pyfunction(signature=(as_json=false,verbose=0))]
fn print_data_status(py: Python<'_>, as_json: bool, verbose: i32) -> PyResult<()> {
    py.detach(|| {
        runtime()
            .data_status_text(as_json, verbose)
            .map(|text| println!("{text}"))
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Eagerly initialize the main kernel, legacy centers and minor-body SPK chunks.
///
/// cache_dir is a str/PathLike SPK output directory; None uses the configured
/// automatic cache, then ~/.libephemeris/spk. Gregorian start_year/end_year
/// labels bound twenty-year chunks (end exclusive, defaults 1550..2650).
/// An empty/reversed year interval schedules no chunks. force=False skips
/// existing chunk paths; True requests new chunks/center replacement. The
/// main image uses normal acquisition/cache behavior, not a forced refresh.
/// show_progress/quiet are accepted but ignored. This synchronous operation
/// may make many serial requests; it has no whole-initialization timeout.
///
/// Return a copied dict: de440 (actual selected-main success, even for another
/// tier), planet_centers (bool); spk_success/failed/skipped/out_of_range (chunk
/// counts); failed_details (list of (body_name, start_year, end_year, message));
/// and spk_blocked (body-name list, omitted when empty). Ordinary main/center
/// failures are logged and SPK failures collected rather than raised. Directory,
/// configuration and invalid-date failures can still raise. Successful files
/// and registrations persist if another step fails; no transactional rollback.
#[pyfunction(signature=(cache_dir=None,force=false,show_progress=true,quiet=false,start_year=1550,end_year=2650))]
fn init_all(
    py: Python<'_>,
    cache_dir: Option<PathBuf>,
    force: bool,
    show_progress: bool,
    quiet: bool,
    start_year: i32,
    end_year: i32,
) -> PyResult<Bound<'_, PyDict>> {
    let _ = (show_progress, quiet);
    let value = py
        .detach(|| runtime().init_all(cache_dir.as_deref(), force, start_year, end_year))
        .map_err(|error| to_pyerr(py, error))?;
    let result = PyDict::new(py);
    result.set_item("de440", value.de440)?;
    result.set_item("planet_centers", value.planet_centers)?;
    result.set_item("spk_success", value.spk_success)?;
    result.set_item("spk_failed", value.spk_failed)?;
    result.set_item("spk_skipped", value.spk_skipped)?;
    result.set_item("spk_out_of_range", value.spk_out_of_range)?;
    result.set_item("failed_details", value.failed_details)?;
    if !value.spk_blocked.is_empty() {
        result.set_item("spk_blocked", value.spk_blocked)?;
    }
    Ok(result)
}
/// Activate tier_name and provision its centers, main kernel and minor-body SPKs.
///
/// tier_name must be "base", "medium" or "extended". Activation clears the
/// explicit JPL filename even for the same tier; explicit LEB selection remains.
/// force=True requests replacement; otherwise verified/registered data can be
/// reused. show_progress/quiet are accepted but ignored. Optional center-file
/// failures are logged; the ensure sweep reports ordinary body fallbacks.
///
/// Return a fresh dict with tier; ephemeris (file/status and a str path only
/// when cached); spks ({body_name: status_or_fallback_reason}); and summary
/// (total/cached/downloaded/fallback/errors counts, including the main file).
/// Unknown tiers raise ValueError and required sealed requests raise
/// NetworkSealedError. Activation and completed installations are not rolled
/// back on later failure. Network/native work releases the GIL.
#[pyfunction(signature=(tier_name,force=false,show_progress=true,quiet=false))]
fn download_for_tier<'py>(
    py: Python<'py>,
    tier_name: &str,
    force: bool,
    show_progress: bool,
    quiet: bool,
) -> PyResult<Bound<'py, PyDict>> {
    let _ = (show_progress, quiet);
    let result = py
        .detach(|| runtime().download_for_tier(tier_name, force))
        .map_err(|error| to_pyerr(py, error))?;
    crate::spk_auto::summary(py, result)
}
/// Install a reviewed tier core and return its pathlib.Path (not necessarily LEB1).
///
/// tier_name is "base", "medium" or "extended". Retain a digest/format-verified
/// file unless force=True; bundled base data can be installed without network,
/// while other acquisition obeys policy. show_progress/quiet are ignored.
/// activate=True selects this explicit LEB path for future use; it does not
/// switch calculation mode or tier. False installs without selecting it.
///
/// Installation can create directories or atomically replace files and retires
/// default reader discovery even if the selected path is unchanged. Retained
/// immutable snapshots stay valid. Unknown tiers/digest/format failures raise
/// ValueError; configuration/policy/transport/filesystem errors propagate.
#[pyfunction(signature=(tier_name,force=false,show_progress=true,quiet=false,activate=true))]
fn download_leb_for_tier<'py>(
    py: Python<'py>,
    tier_name: &str,
    force: bool,
    show_progress: bool,
    quiet: bool,
    activate: bool,
) -> PyResult<Bound<'py, PyAny>> {
    let _ = (show_progress, quiet);
    let path = py
        .detach(|| runtime().download_leb_for_tier(tier_name, force, activate))
        .map_err(|error| download_error(py, error))?;
    path_object(py, &path)
}
/// Install selected reviewed LEB2 groups; return successful pathlib.Path values.
///
/// tier_name is "base", "medium" or "extended". groups=None requests core,
/// asteroids, exotics and apogee; an explicit sequence is processed in supplied
/// order, with unknown/unpublished names skipped. An empty list installs nothing
/// (but still creates the leb directory). force=True requests replacement;
/// otherwise verified cached files are retained. show_progress/quiet are ignored.
///
/// Each group uses atomic verified installation. Per-group acquisition errors
/// are logged and skipped: an empty/partial result is possible without an
/// exception. Tier/directory/configuration errors outside that step can raise.
/// activate=True selects only a successfully installed core, not companions,
/// and changes neither mode nor tier. False leaves selection unchanged.
/// Successful installations retire default readers; existing snapshots keep
/// their images. This is not an all-or-nothing transaction or a data generator.
#[pyfunction(signature=(tier_name,groups=None,force=false,show_progress=true,quiet=false,activate=true))]
fn download_leb2_for_tier<'py>(
    py: Python<'py>,
    tier_name: &str,
    groups: Option<Vec<String>>,
    force: bool,
    show_progress: bool,
    quiet: bool,
    activate: bool,
) -> PyResult<Bound<'py, PyList>> {
    let _ = (show_progress, quiet);
    let paths = py
        .detach(|| runtime().download_leb2_for_tier(tier_name, groups.as_deref(), force, activate))
        .map_err(|error| download_error(py, error))?;
    paths_object(py, paths)
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(get_data_dir, module)?)?;
    module.add_function(wrap_pyfunction!(get_leb_dir, module)?)?;
    module.add_function(wrap_pyfunction!(get_leb_path_for_tier, module)?)?;
    module.add_function(wrap_pyfunction!(publish_temp_file, module)?)?;
    module.add_function(wrap_pyfunction!(download_file, module)?)?;
    module.add_function(wrap_pyfunction!(download_planet_centers, module)?)?;
    module.add_function(wrap_pyfunction!(download_all, module)?)?;
    module.add_function(wrap_pyfunction!(check_data_status, module)?)?;
    module.add_function(wrap_pyfunction!(print_data_status, module)?)?;
    module.add_function(wrap_pyfunction!(init_all, module)?)?;
    module.add_function(wrap_pyfunction!(download_for_tier, module)?)?;
    module.add_function(wrap_pyfunction!(download_leb_for_tier, module)?)?;
    module.add_function(wrap_pyfunction!(download_leb2_for_tier, module)?)?;
    let files = PyDict::new(module.py());
    for entry in libephemeris_core::download::DATA_FILES {
        let row = PyDict::new(module.py());
        if let Some(url) = entry.url() {
            row.set_item("url", url)?;
        }
        if entry.is_bundled_only() {
            row.set_item("bundled_resource", "native-embedded/base_core.leb2")?;
        }
        row.set_item("sha256", entry.sha256())?;
        row.set_item("size_mb", entry.size_mb())?;
        row.set_item("description", entry.description())?;
        if let Some(subdir) = entry.dest_subdir() {
            row.set_item("dest_subdir", subdir)?;
        }
        files.set_item(entry.name(), row)?;
    }
    module.add("_download_DATA_FILES", files)?;
    Ok(())
}
