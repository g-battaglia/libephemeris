//! Automatic SPK configuration, cache inspection and native acquisition helpers.
use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::config::spk_auto::{AutoSpkError, AutoSpkSettings, SpkCacheInfo};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

fn error(py: Python<'_>, value: AutoSpkError) -> PyErr {
    match value {
        AutoSpkError::Invalid(message) => PyValueError::new_err(message),
        AutoSpkError::Native(error) => to_pyerr(py, error),
    }
}
#[derive(FromPyObject)]
enum Target {
    Integer(i64),
    Text(String),
}
impl Target {
    fn text(self) -> String {
        match self {
            Self::Integer(value) => value.to_string(),
            Self::Text(value) => value,
        }
    }
}

/// Mutable automatic-SPK settings, either standalone or shared with the runtime.
///
/// Parameters are a library body ID (ipl), Horizons target string (body_id),
/// optional NAIF target code, start/end Gregorian date labels ('YYYY-MM-DD'),
/// and an optional cache path. Direct construction stores values without eager
/// date/path validation or registration; enabled starts True and spk_path None.
/// Use enable_auto_spk() to install settings, then get_auto_spk_config() to edit
/// the live record. A standalone instance does not enable automatic acquisition.
///
/// Records returned by lookup/list functions share mutable runtime settings.
/// Property assignment does not download, clear prior registration or reset
/// acquisition metadata. Changing ipl does not rekey the runtime's dictionary.
/// Disabling/replacing a record leaves retained handles alive but detached from
/// that dictionary entry. Individual accesses are synchronized, not a transaction
/// across several property updates. There is no close() method.
///
/// Example:
///
/// ```python
/// from libephemeris.spk_auto import enable_auto_spk, get_auto_spk_config
/// enable_auto_spk(17, '1;', naif_id=2000001)
/// config = get_auto_spk_config(17)
/// config.enabled = False  # Edit the installed record; no download.
/// ```
#[pyclass(module = "libephemeris.spk_auto", skip_from_py_object)]
#[derive(Clone)]
struct AutoSpkConfig {
    inner: Arc<Mutex<AutoSpkSettings>>,
}
impl AutoSpkConfig {
    fn snapshot(&self) -> AutoSpkSettings {
        self.inner.lock().unwrap_or_else(|p| p.into_inner()).clone()
    }
}
#[pymethods]
impl AutoSpkConfig {
    /// Create unregistered enabled settings; see class docs for units and ownership.
    #[new]
    #[pyo3(signature=(ipl,body_id,naif_id=None,start="2000-01-01",end="2100-01-01",cache_dir=None))]
    fn new(
        ipl: i32,
        body_id: String,
        naif_id: Option<i32>,
        start: &str,
        end: &str,
        cache_dir: Option<PathBuf>,
    ) -> Self {
        Self {
            inner: Arc::new(Mutex::new(AutoSpkSettings::new(
                ipl,
                body_id,
                naif_id,
                start.to_owned(),
                end.to_owned(),
                cache_dir,
            ))),
        }
    }
    /// Return the effective directory string, creating missing directories.
    ///
    /// Uses nonempty cache_dir, then runtime selection, then the conventional
    /// user cache. A relative explicit path stays relative. No download occurs.
    /// Configuration/path failures raise ConfigurationError, FileNotFoundError
    /// or DataNotFoundError as mapped from the underlying operation.
    fn get_cache_dir(&self, py: Python<'_>) -> PyResult<String> {
        py.detach(|| runtime().auto_spk_settings_dir(&self.snapshot()))
            .map(|path| path.to_string_lossy().into_owned())
            .map_err(|value| to_pyerr(py, value))
    }
    /// Return the opaque target/date-dependent BSP basename without doing I/O.
    ///
    /// The short digest is an identifier, not a file-integrity verification.
    /// This is distinct from the top-level get_cache_path() simple target path.
    fn get_cache_filename(&self, py: Python<'_>) -> String {
        py.detach(|| self.snapshot().cache_filename())
    }
    /// Return effective directory plus opaque filename, creating the directory.
    ///
    /// Does not require a BSP to exist or download/register it. Directory errors
    /// follow get_cache_dir(); returning a path does not certify its coverage.
    fn get_cache_path(&self, py: Python<'_>) -> PyResult<String> {
        py.detach(|| runtime().auto_spk_settings_path(&self.snapshot()))
            .map(|path| path.to_string_lossy().into_owned())
            .map_err(|value| to_pyerr(py, value))
    }
    /// Library body ID; assignment does not change the runtime dictionary key.
    #[getter]
    fn ipl(&self, py: Python<'_>) -> i32 {
        py.detach(|| self.snapshot().ipl)
    }
    /// Store the ID without rekeying settings or updating SPK registrations.
    #[setter]
    fn set_ipl(&self, py: Python<'_>, value: i32) {
        py.detach(|| self.inner.lock().unwrap_or_else(|p| p.into_inner()).ipl = value);
    }
    /// Horizons target identifier/name string; edits do not acquire a new image.
    #[getter]
    fn body_id(&self, py: Python<'_>) -> String {
        py.detach(|| self.snapshot().body_id)
    }
    /// Store a target string without validation or clearing prior spk_path.
    #[setter]
    fn set_body_id(&self, py: Python<'_>, value: String) {
        py.detach(|| self.inner.lock().unwrap_or_else(|p| p.into_inner()).body_id = value);
    }
    /// Preferred integer NAIF code, or None to infer it during acquisition.
    #[getter]
    fn naif_id(&self, py: Python<'_>) -> Option<i32> {
        py.detach(|| self.snapshot().naif_id)
    }
    /// Store the preferred target; this does not update an existing registration.
    #[setter]
    fn set_naif_id(&self, py: Python<'_>, value: Option<i32>) {
        py.detach(|| self.inner.lock().unwrap_or_else(|p| p.into_inner()).naif_id = value);
    }
    /// Requested inclusive Gregorian start-date label; no eager date validation.
    #[getter]
    fn start(&self, py: Python<'_>) -> String {
        py.detach(|| self.snapshot().start)
    }
    /// Store the start label without downloading or resetting acquisition metadata.
    #[setter]
    fn set_start(&self, py: Python<'_>, value: String) {
        py.detach(|| self.inner.lock().unwrap_or_else(|p| p.into_inner()).start = value);
    }
    /// Requested inclusive Gregorian end-date label; no eager date validation.
    #[getter]
    fn end(&self, py: Python<'_>) -> String {
        py.detach(|| self.snapshot().end)
    }
    /// Store the end label without downloading or resetting acquisition metadata.
    #[setter]
    fn set_end(&self, py: Python<'_>, value: String) {
        py.detach(|| self.inner.lock().unwrap_or_else(|p| p.into_inner()).end = value);
    }
    /// Stored cache-directory string or None, not necessarily the resolved directory.
    ///
    /// Reading creates nothing; get_cache_dir() resolves defaults and may create
    /// directories. Assignment accepts a path-like value or None without I/O.
    #[getter(cache_dir)]
    fn cache_directory(&self, py: Python<'_>) -> Option<String> {
        py.detach(|| {
            self.snapshot()
                .cache_dir
                .map(|path| path.to_string_lossy().into_owned())
        })
    }
    /// Store a directory override without moving/deleting previously acquired files.
    #[setter]
    fn set_cache_dir(&self, py: Python<'_>, value: Option<PathBuf>) {
        py.detach(|| {
            self.inner
                .lock()
                .unwrap_or_else(|p| p.into_inner())
                .cache_dir = value
        });
    }
    /// Opportunistic acquisition switch; False does not unregister an acquired SPK.
    #[getter]
    fn enabled(&self, py: Python<'_>) -> bool {
        py.detach(|| self.snapshot().enabled)
    }
    /// Change the opportunistic switch without downloading or removing the record.
    #[setter]
    fn set_enabled(&self, py: Python<'_>, value: bool) {
        py.detach(|| self.inner.lock().unwrap_or_else(|p| p.into_inner()).enabled = value);
    }
    /// Last acquired path string or None; an informational field, not a live file.
    ///
    /// Assignment does not verify, open or register the selected file. Other
    /// property edits do not automatically clear this value.
    #[getter]
    fn spk_path(&self, py: Python<'_>) -> Option<String> {
        py.detach(|| {
            self.snapshot()
                .spk_path
                .map(|path| path.to_string_lossy().into_owned())
        })
    }
    /// Store acquisition metadata only; no file/coverage check or registration.
    #[setter]
    fn set_spk_path(&self, py: Python<'_>, value: Option<PathBuf>) {
        py.detach(|| {
            self.inner
                .lock()
                .unwrap_or_else(|p| p.into_inner())
                .spk_path = value
        });
    }
}

fn cache_info_dict<'py>(py: Python<'py>, value: SpkCacheInfo) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("cache_dir", value.directory.to_string_lossy().as_ref())?;
    result.set_item("num_files", value.files.len())?;
    result.set_item("total_size_mb", value.size_mb())?;
    result.set_item("files", value.files)?;
    Ok(result)
}
/// Resolve/create the automatic-SPK directory and return its absolute path string.
///
/// cache_dir is an optional path-like override; None uses the runtime selection,
/// then ~/.libephemeris/spk (or the relative home fallback). Does not download.
/// Configuration and filesystem failures propagate as public native exceptions.
#[pyfunction(signature=(cache_dir=None))]
fn ensure_cache_dir(py: Python<'_>, cache_dir: Option<PathBuf>) -> PyResult<String> {
    py.detach(|| runtime().auto_spk_cache_dir(cache_dir.as_deref(), true))
        .map(|path| path.to_string_lossy().into_owned())
        .map_err(|value| to_pyerr(py, value))
}
/// Return an absolute simple sanitized-target BSP path, without creating it.
///
/// body_id accepts an integer or Horizons target string. None cache_dir uses
/// the conventional user cache, not the runtime cache override. This helper is
/// distinct from AutoSpkConfig's date-dependent filename. Path-resolution errors
/// propagate; existence, target identity and coverage are not checked.
#[pyfunction(signature=(body_id,cache_dir=None))]
fn get_cache_path(py: Python<'_>, body_id: Target, cache_dir: Option<PathBuf>) -> PyResult<String> {
    let body = body_id.text();
    py.detach(|| runtime().simple_spk_cache_path(&body, cache_dir.as_deref()))
        .map(|path| path.to_string_lossy().into_owned())
        .map_err(|value| to_pyerr(py, value))
}
/// Return a new dict of conventional/explicit cache filesystem statistics.
///
/// Keys: cache_dir (absolute str), num_files (int), total_size_mb (MiB rounded
/// to two decimals), files (list of BSP basenames in filesystem order). None
/// cache_dir uses the conventional root, not the runtime override. A missing
/// directory yields zero files; no directory is created or BSP evaluated.
/// Directory/stat failures propagate as native file/configuration exceptions.
#[pyfunction(signature=(cache_dir=None))]
fn cache_info<'py>(py: Python<'py>, cache_dir: Option<PathBuf>) -> PyResult<Bound<'py, PyDict>> {
    let value = py
        .detach(|| runtime().auto_spk_cache_info(cache_dir.as_deref()))
        .map_err(|value| to_pyerr(py, value))?;
    cache_info_dict(py, value)
}
/// Install/replace automatic acquisition settings for library body ID ipl.
///
/// body_id is a Horizons target string; naif_id is an optional preferred NAIF
/// code. start/end are Gregorian date labels; cache_dir is an optional path.
/// None captures a configured runtime cache override when present. Returns None;
/// no date validation, directory creation or download occurs at this step.
/// Configuration-resolution errors propagate. See AutoSpkConfig for live edits.
#[pyfunction(signature=(ipl,body_id,naif_id=None,start="2000-01-01",end="2100-01-01",cache_dir=None))]
fn enable_auto_spk(
    py: Python<'_>,
    ipl: i32,
    body_id: String,
    naif_id: Option<i32>,
    start: &str,
    end: &str,
    cache_dir: Option<PathBuf>,
) -> PyResult<()> {
    py.detach(|| {
        runtime().enable_auto_spk(AutoSpkSettings::new(
            ipl,
            body_id,
            naif_id,
            start.to_owned(),
            end.to_owned(),
            cache_dir,
        ))
    })
    .map_err(|value| to_pyerr(py, value))
}
/// Remove automatic settings for library ID ipl; absent IDs are a no-op.
///
/// Returns None; acquired registrations/files and retained config handles survive.
#[pyfunction]
fn disable_auto_spk(py: Python<'_>, ipl: i32) {
    py.detach(|| runtime().disable_auto_spk(ipl));
}
/// Return True only if library ID ipl has installed settings with enabled=True.
#[pyfunction]
fn is_auto_spk_enabled(py: Python<'_>, ipl: i32) -> bool {
    py.detach(|| runtime().is_auto_spk_enabled(ipl))
}
/// Return a live shared AutoSpkConfig for library ID ipl, or None if absent.
///
/// Mutating its properties edits the installed settings, not a detached copy.
/// Disabling/replacing the entry does not destroy this retained handle.
#[pyfunction]
fn get_auto_spk_config(py: Python<'_>, ipl: i32) -> Option<AutoSpkConfig> {
    py.detach(|| {
        runtime()
            .auto_spk_config(ipl)
            .map(|inner| AutoSpkConfig { inner })
    })
}
/// Return a new dict of library IDs to live shared AutoSpkConfig handles.
///
/// Editing the dict does not edit the registry; editing a config does edit its
/// shared record. Map keys are not updated when a config's ipl field changes.
#[pyfunction]
fn list_auto_spk_bodies(py: Python<'_>) -> HashMap<i32, AutoSpkConfig> {
    py.detach(|| {
        runtime()
            .auto_spk_configs()
            .into_iter()
            .map(|(key, inner)| (key, AutoSpkConfig { inner }))
            .collect()
    })
}
/// Acquire and register the installed configuration for library ID ipl now.
///
/// Returns its path string and updates spk_path/naif_id. force=True replaces a
/// cached file; otherwise a structurally valid cached file is retained. Explicit
/// acquisition is attempted even if the record's opportunistic switch is False.
/// Creates the cache directory. Raises ValueError if unconfigured or no target
/// can be deduced, NetworkSealedError for denied network access, and native
/// filesystem, transport, SPK or registration exceptions for actual failures.
#[pyfunction(signature=(ipl,force=false))]
fn download_now(py: Python<'_>, ipl: i32, force: bool) -> PyResult<String> {
    py.detach(|| runtime().download_auto_spk_now(ipl, force))
        .map(|path| path.to_string_lossy().into_owned())
        .map_err(|value| error(py, value))
}
/// Opportunistically acquire library ID ipl and return its path string or None.
///
/// None means absent/disabled settings or any acquisition/registration failure;
/// those native failures are suppressed. Use download_now() to preserve errors.
#[pyfunction]
fn try_auto_download(py: Python<'_>, ipl: i32) -> Option<String> {
    py.detach(|| {
        runtime()
            .try_auto_download(ipl)
            .map(|path| path.to_string_lossy().into_owned())
    })
}
/// Delete acquired paths referenced by one config or all configs; return a count.
///
/// ipl=None visits all installed records and ignores individual deletion failures;
/// a specific library ID preserves file-deletion errors. Successful deletions
/// clear their spk_path metadata. No settings/registration is removed, unrelated
/// files are untouched, and images retained by calculations remain owned.
#[pyfunction(signature=(ipl=None))]
fn clear_cache(py: Python<'_>, ipl: Option<i32>) -> PyResult<usize> {
    py.detach(|| runtime().clear_auto_spk_cache(ipl))
        .map_err(|value| to_pyerr(py, value))
}
/// Return cache_info() for the conventional user cache root.
///
/// Returns a fresh dict with directory, file count, rounded MiB total and names.
/// Does not create the directory or include other runtime/custom cache roots;
/// native directory/stat failures propagate as in cache_info().
#[pyfunction]
fn get_cache_info<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let value = py
        .detach(|| runtime().auto_spk_cache_info(None))
        .map_err(|value| to_pyerr(py, value))?;
    cache_info_dict(py, value)
}
/// Return new metadata dicts for BSP files in the selected cache roots.
///
/// An explicit cache_dir selects one root; None visits conventional and runtime
/// roots without duplicating an identical root. Each dict has filename/path,
/// size_bytes, size_mb (MiB rounded to two decimals), jd_start/jd_end (TDB
/// envelope or None), date_start/date_end (civil labels or None), last_accessed
/// (Unix seconds, 0 if unavailable). Coverage is a file envelope, not a continuous
/// target-chain guarantee. Missing/corrupt SPK coverage is None; filesystem or
/// configuration errors propagate. No data is downloaded.
#[pyfunction(signature=(cache_dir=None))]
fn list_cached_spk<'py>(
    py: Python<'py>,
    cache_dir: Option<PathBuf>,
) -> PyResult<Bound<'py, PyList>> {
    let values = py
        .detach(|| runtime().cached_spk_files(cache_dir.as_deref()))
        .map_err(|value| to_pyerr(py, value))?;
    let result = PyList::empty(py);
    for value in values {
        let item = PyDict::new(py);
        let dates = value.dates();
        item.set_item("filename", &value.filename)?;
        item.set_item("path", value.path.to_string_lossy().as_ref())?;
        item.set_item("size_bytes", value.size_bytes)?;
        item.set_item("size_mb", value.size_mb())?;
        item.set_item("jd_start", value.coverage.map(|span| span.0))?;
        item.set_item("jd_end", value.coverage.map(|span| span.1))?;
        item.set_item("date_start", dates.as_ref().map(|span| span.0.as_str()))?;
        item.set_item("date_end", dates.as_ref().map(|span| span.1.as_str()))?;
        item.set_item("last_accessed", value.last_accessed)?;
        result.append(item)?;
    }
    Ok(result)
}
/// Delete BSP files in selected roots and return the successful deletion count.
///
/// cache_dir=None visits conventional and runtime roots; an explicit path selects
/// one. Unrelated files remain. Individual deletion failures are skipped, while
/// directory/configuration failures propagate. Retires missing spk_path metadata
/// without disabling settings or unregistering retained calculation images.
#[pyfunction(signature=(cache_dir=None))]
fn clear_spk_cache(py: Python<'_>, cache_dir: Option<PathBuf>) -> PyResult<usize> {
    py.detach(|| runtime().clear_spk_cache(cache_dir.as_deref()))
        .map_err(|value| to_pyerr(py, value))
}
/// Return selected BSP storage in MiB, rounded to two decimal places.
///
/// None cache_dir visits conventional and runtime roots; an explicit path selects
/// one. Missing directories contribute zero; individual unreadable stat results
/// are skipped. Directory/configuration failures propagate. This is disk size,
/// not memory usage, and the operation does not create directories.
#[pyfunction(signature=(cache_dir=None))]
fn get_cache_size(py: Python<'_>, cache_dir: Option<PathBuf>) -> PyResult<f64> {
    py.detach(|| runtime().spk_cache_size_mb(cache_dir.as_deref()))
        .map_err(|value| to_pyerr(py, value))
}
/// Delete BSP files last accessed more than max_age_days ago; return a count.
///
/// max_age_days must be a positive integer; a day is 86400 seconds. cache_dir
/// selects one root, or None visits conventional/runtime roots. Uses filesystem
/// access times, not ephemeris dates; missing access times/deletion failures are
/// skipped. Raises ValueError for nonpositive age and native directory/config
/// exceptions. Missing acquired paths are retired, not their registrations.
#[pyfunction(signature=(max_age_days,cache_dir=None))]
fn prune_old_cache(
    py: Python<'_>,
    max_age_days: i64,
    cache_dir: Option<PathBuf>,
) -> PyResult<usize> {
    py.detach(|| runtime().prune_spk_cache(max_age_days, cache_dir.as_deref()))
        .map_err(|value| error(py, value))
}
/// Enable Chiron, Pholus, Ceres, Pallas, Juno, Vesta, Eris and Sedna settings.
///
/// start/end are Gregorian date labels; cache_dir is an optional path override.
/// Returns None, replacing these records without downloading or validating dates.
/// Other records are unchanged. Configuration errors propagate.
#[pyfunction(signature=(start="2000-01-01",end="2100-01-01",cache_dir=None))]
fn enable_common_bodies(
    py: Python<'_>,
    start: &str,
    end: &str,
    cache_dir: Option<PathBuf>,
) -> PyResult<()> {
    py.detach(|| runtime().enable_common_spk_bodies(start, end, cache_dir))
        .map_err(|value| to_pyerr(py, value))
}
/// Remove all automatic-SPK settings; return None without deleting cached files.
///
/// Existing registrations and external config/image handles remain alive.
#[pyfunction]
fn disable_all(py: Python<'_>) {
    py.detach(|| runtime().disable_all_auto_spk());
}
/// Test whether a target-prefixed cached BSP passes the native envelope check.
///
/// body_id accepts integer/string Horizons identifiers; jd_start/jd_end are TDB
/// Julian days compared without clock conversion or download padding. cache_dir
/// overrides the runtime/default root. The check permits a 1.5-day start and
/// 0.5-day end boundary margin; it does not validate ordering/finiteness, target
/// identity or chain continuity. NaN comparisons fail; other values follow the
/// comparisons as supplied. Returns False when no matching envelope passes.
/// Does not create the directory/download; path/configuration failures propagate.
#[pyfunction(signature=(body_id,jd_start,jd_end,cache_dir=None))]
fn is_spk_cached(
    py: Python<'_>,
    body_id: Target,
    jd_start: f64,
    jd_end: f64,
    cache_dir: Option<PathBuf>,
) -> PyResult<bool> {
    let body = body_id.text();
    py.detach(|| runtime().is_spk_cached(&body, jd_start, jd_end, cache_dir.as_deref()))
        .map_err(|value| to_pyerr(py, value))
}
/// Reuse a covering cached kernel or acquire a padded interval; return its path.
///
/// body_id accepts integer/string Horizons identifiers. Finite TDB JD bounds
/// require jd_end > jd_start; configured whole-day padding is applied to both
/// ends. Downloads convert JDs to Gregorian date labels, not exact sub-day
/// endpoints. cache_dir selects the runtime/default root and is created.
/// If ipl is supplied, register that library ID using optional preferred naif_id;
/// otherwise no body registration occurs. Cache checks permit native 1.5-day
/// start/0.5-day end boundary margins and do not prove continuous target coverage.
/// Raises ValueError for invalid intervals/target
/// deduction, NetworkSealedError when acquisition is denied, and native file,
/// transport or SPK exceptions. No numerical calculation is performed here.
#[pyfunction(signature=(body_id,jd_start,jd_end,cache_dir=None,ipl=None,naif_id=None))]
fn auto_get_spk(
    py: Python<'_>,
    body_id: Target,
    jd_start: f64,
    jd_end: f64,
    cache_dir: Option<PathBuf>,
    ipl: Option<i32>,
    naif_id: Option<i32>,
) -> PyResult<String> {
    let body = body_id.text();
    py.detach(|| {
        runtime().auto_get_spk(&body, jd_start, jd_end, cache_dir.as_deref(), ipl, naif_id)
    })
    .map(|path| path.to_string_lossy().into_owned())
    .map_err(|value| error(py, value))
}
/// Download a Horizons SPK to output_path and optionally register a library ID.
///
/// body_id is an integer/string target; finite jd_start/jd_end require end >
/// start and are converted to Gregorian date labels without TT/UTC conversion
/// or sub-day endpoint preservation. location is a Horizons center string;
/// '@0' selects the SSB (normalized to '500@0'). Existing output can be replaced.
/// Returns the path string. ipl enables registration; naif_id is its optional
/// preferred NAIF target, unused without ipl. Raises ValueError for invalid
/// bounds/undeducible target, NetworkSealedError for refused acquisition and
/// native filesystem, transport or SPK exceptions. No automatic padding is added.
#[pyfunction(signature=(body_id,jd_start,jd_end,output_path,location="@0",ipl=None,naif_id=None))]
#[allow(clippy::too_many_arguments)]
fn download_spk_from_horizons(
    py: Python<'_>,
    body_id: Target,
    jd_start: f64,
    jd_end: f64,
    output_path: PathBuf,
    location: &str,
    ipl: Option<i32>,
    naif_id: Option<i32>,
) -> PyResult<String> {
    let body = body_id.text();
    py.detach(|| {
        runtime().download_spk_interval(
            &body,
            jd_start,
            jd_end,
            &output_path,
            location,
            ipl,
            naif_id,
        )
    })
    .map(|path| path.to_string_lossy().into_owned())
    .map_err(|value| error(py, value))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<AutoSpkConfig>()?;
    module.add("DEFAULT_START_DATE", "2000-01-01")?;
    module.add("DEFAULT_END_DATE", "2100-01-01")?;
    module.add(
        "DEFAULT_AUTO_SPK_DIR",
        runtime().default_auto_spk_dir().to_string_lossy().as_ref(),
    )?;
    module.add_function(wrap_pyfunction!(ensure_cache_dir, module)?)?;
    module.add_function(wrap_pyfunction!(get_cache_path, module)?)?;
    module.add_function(wrap_pyfunction!(cache_info, module)?)?;
    module.add_function(wrap_pyfunction!(enable_auto_spk, module)?)?;
    module.add_function(wrap_pyfunction!(disable_auto_spk, module)?)?;
    module.add_function(wrap_pyfunction!(is_auto_spk_enabled, module)?)?;
    module.add_function(wrap_pyfunction!(get_auto_spk_config, module)?)?;
    module.add_function(wrap_pyfunction!(list_auto_spk_bodies, module)?)?;
    module.add_function(wrap_pyfunction!(download_now, module)?)?;
    module.add_function(wrap_pyfunction!(try_auto_download, module)?)?;
    module.add_function(wrap_pyfunction!(clear_cache, module)?)?;
    module.add_function(wrap_pyfunction!(get_cache_info, module)?)?;
    module.add_function(wrap_pyfunction!(list_cached_spk, module)?)?;
    module.add_function(wrap_pyfunction!(clear_spk_cache, module)?)?;
    module.add_function(wrap_pyfunction!(get_cache_size, module)?)?;
    module.add_function(wrap_pyfunction!(prune_old_cache, module)?)?;
    module.add_function(wrap_pyfunction!(enable_common_bodies, module)?)?;
    module.add_function(wrap_pyfunction!(disable_all, module)?)?;
    module.add_function(wrap_pyfunction!(is_spk_cached, module)?)?;
    module.add_function(wrap_pyfunction!(auto_get_spk, module)?)?;
    module.add_function(wrap_pyfunction!(download_spk_from_horizons, module)?)?;
    Ok(())
}
