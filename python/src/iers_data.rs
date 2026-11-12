//! Thin IERS table and cache-lifecycle bindings.
use crate::{errors::to_pyerr, positions::runtime};
use libephemeris_core::iers::IersFileKind;
use libephemeris_core::{iers, time};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::path::PathBuf;

/// Mutable value record for one daily UT1-minus-UTC Earth-orientation row.
///
/// mjd uses MJD = JD - 2400000.5; offsets/uncertainties are seconds. Calendar
/// labels and numeric values are stored as supplied, without consistency or
/// finite-value validation. Editing this record does not update loaded tables.
/// It owns no cache/file resource and performs no I/O.
#[pyclass(
    module = "libephemeris.iers_data",
    get_all,
    set_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, PartialEq)]
struct IERSDataPoint {
    /// Row's Modified Julian Day; native queries use the UT1 date coordinate.
    mjd: f64,
    /// Calendar year label; not cross-validated against mjd.
    year: i32,
    /// Calendar month label, conventionally 1 through 12.
    month: i32,
    /// Calendar day-of-month label.
    day: i32,
    /// UT1 minus UTC in seconds, not days.
    ut1_utc: f64,
    /// UT1-minus-UTC uncertainty in seconds, or None when unspecified.
    ut1_utc_err: Option<f64>,
    /// Whether the source row is predicted rather than observed.
    is_prediction: bool,
}
#[pymethods]
impl IERSDataPoint {
    /// Store the required date/UT1-minus-UTC fields; uncertainty defaults to
    /// None and is_prediction to False. No table insertion or validation occurs.
    #[new]
    #[pyo3(signature=(mjd, year, month, day, ut1_utc, ut1_utc_err=None, is_prediction=false))]
    fn new(
        mjd: f64,
        year: i32,
        month: i32,
        day: i32,
        ut1_utc: f64,
        ut1_utc_err: Option<f64>,
        is_prediction: bool,
    ) -> Self {
        Self {
            mjd,
            year,
            month,
            day,
            ut1_utc,
            ut1_utc_err,
            is_prediction,
        }
    }
}

/// Mutable value record for a UTC leap-table entry, with TAI-minus-UTC seconds.
///
/// mjd marks the entry's effective date (JD - 2400000.5); year/month/day are
/// calendar labels. Construction and later assignments do not validate dates,
/// consistency or finite values, and do not change the runtime's leap table.
#[pyclass(
    module = "libephemeris.iers_data",
    get_all,
    set_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, PartialEq)]
struct LeapSecondEntry {
    /// Effective UTC Modified Julian Day (JD - 2400000.5).
    mjd: f64,
    /// TAI minus UTC from this entry onward, in seconds.
    tai_utc: f64,
    /// Calendar year of the effective-date label.
    year: i32,
    /// Calendar month of the effective-date label, conventionally 1 through 12.
    month: i32,
    /// Calendar day of the effective-date label.
    day: i32,
}
#[pymethods]
impl LeapSecondEntry {
    /// Store mjd, tai_utc (seconds), year, month and day; all are required.
    /// This creates only a value record, without validation, I/O or table edits.
    #[new]
    fn new(mjd: f64, tai_utc: f64, year: i32, month: i32, day: i32) -> Self {
        Self {
            mjd,
            tai_utc,
            year,
            month,
            day,
        }
    }
}

/// Mutable value record for observed TT minus UT1 (Delta-T), in seconds.
///
/// mjd is the native UT1 date coordinate (JD - 2400000.5). These seconds differ
/// from deltat()/LEB delta_t() results, which are days. Values/date labels are
/// stored without validation; editing a record never alters the loaded tables.
#[pyclass(
    module = "libephemeris.iers_data",
    get_all,
    set_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, PartialEq)]
struct DeltaTDataPoint {
    /// UT1 Modified Julian Day of the row (JD - 2400000.5).
    mjd: f64,
    /// Calendar year label, not cross-validated against mjd.
    year: i32,
    /// Calendar month label, conventionally 1 through 12.
    month: i32,
    /// Calendar day-of-month label.
    day: i32,
    /// TT minus UT1 in seconds, not days.
    delta_t: f64,
}
#[pymethods]
impl DeltaTDataPoint {
    /// Store mjd, year, month, day and delta_t (seconds); all are required.
    /// No finite/date validation, table insertion or data acquisition occurs.
    #[new]
    fn new(mjd: f64, year: i32, month: i32, day: i32, delta_t: f64) -> Self {
        Self {
            mjd,
            year,
            month,
            day,
            delta_t,
        }
    }
}

macro_rules! download {
    ($name:ident,$kind:ident,$doc:literal) => {
        #[doc = $doc]
        #[pyfunction(signature=(force=false))]
        fn $name(py: Python<'_>, force: bool) -> PyResult<String> {
            py.detach(|| runtime().download_iers_file(IersFileKind::$kind, force))
                .map(|path| path.to_string_lossy().into_owned())
                .map_err(|error| to_pyerr(py, error))
        }
    };
}
download!(
    download_iers_finals,
    Finals,
    r#"Ensure finals2000A.data and return its cache filename as str.

force=False keeps an existing file younger than 30 days; True requests
replacement. Otherwise try the primary then backup URL, honoring network
policy, and publish each successful download atomically. The cache directory
may be created even if transport fails. Transport/policy/filesystem errors
propagate; this explicit download does not consult the auto-download toggle
or reload already parsed rows. Use load_iers_data() to reparse them."#
);
download!(
    download_leap_seconds,
    LeapSeconds,
    r#"Ensure leap_seconds.dat and return its cache filename as str.

force=False keeps a file younger than 30 days; True requests replacement.
Missing/stale files try primary then backup URLs under the network policy.
This may create the cache directory; failed transport does not replace the
installed file. Acquisition errors propagate. This explicit download ignores
the auto-download toggle and does not refresh loaded leap rows; reparse with
load_iers_data(). Offsets in the table are TAI-minus-UTC seconds."#
);
download!(
    download_delta_t_data,
    DeltaT,
    r#"Ensure deltat.data and return its cache filename as str.

force=False keeps an existing file younger than 30 days; True requests
replacement. Otherwise try primary then backup URLs under network policy,
creating the cache directory as needed and publishing atomically. Errors
propagate; auto-download need not be enabled for this explicit operation.
Loaded rows are unchanged until load_iers_data() reparses them. These
observations store TT-minus-UT1 seconds, not the days returned by deltat()."#
);
/// Reparse local finals, leap and observed-Delta-T files into a new generation.
///
/// With automatic download enabled, acquire missing files; force_download=True
/// requests refresh of all three. It does not override that toggle or network
/// policy. Ordinary download failures leave local data usable; configuration
/// and directory/acquisition errors can raise. Completed downloads persist.
///
/// Return True if any table has rows, including the packaged leap fallback;
/// True alone does not establish UT1 or Delta-T coverage. Missing finals retain
/// previous finals rows on explicit reload; missing monthly rows become empty.
/// Active snapshots keep their prior rows. Loading does not enable IERS use
/// in the model or change independent frame clocks.
#[pyfunction(signature=(force_download=false))]
fn load_iers_data(py: Python<'_>, force_download: bool) -> PyResult<bool> {
    py.detach(|| runtime().load_iers_data(force_download))
        .map_err(|error| to_pyerr(py, error))
}
/// Interpolate observed TT-minus-UT1 seconds at finite UT1 Julian Day jd.
///
/// Return None outside the monthly table's inclusive interval or when empty;
/// no finals/model substitution or extrapolation occurs. Lazily load cached
/// rows, possibly creating/downloading missing files under configured policy.
/// Nonfinite jd raises ValueError; acquisition/configuration errors propagate.
#[pyfunction]
fn get_observed_delta_t(py: Python<'_>, jd: f64) -> PyResult<Option<f64>> {
    py.detach(|| {
        let rows = runtime().iers_rows()?;
        iers::observed_delta_t_for_jd(jd, rows.observed_rows())
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Return the observed monthly table's inclusive (start, end) UT1 JD, or None.
///
/// This is JD, not MJD or TT. A cold query lazily loads the shared row cache
/// and may create/acquire missing files; acquisition/configuration errors raise.
/// The range describes available rows, not an extrapolation domain.
#[pyfunction]
fn get_observed_delta_t_data_range(py: Python<'_>) -> PyResult<Option<(f64, f64)>> {
    py.detach(|| {
        runtime()
            .iers_rows()
            .map(|rows| iers::observed_delta_t_range_jd(rows.observed_rows()))
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Return whether finite UT1 Julian Day jd lies in the monthly table interval.
///
/// Empty tables return False. This lazily acquires rows like
/// get_observed_delta_t(), rather than testing a model fallback. Nonfinite
/// jd raises ValueError; configuration/acquisition errors propagate.
#[pyfunction]
fn is_observed_delta_t_available(py: Python<'_>, jd: f64) -> PyResult<bool> {
    py.detach(|| {
        let rows = runtime().iers_rows()?;
        iers::is_observed_delta_t_available(jd, rows.observed_rows())
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Return IERS TT-minus-UT1 seconds at finite UT1 Julian Day jd, or None.
///
/// Covered monthly observations take priority; otherwise finals/leaps give
/// TAI-minus-UTC + 32.184 - (UT1-minus-UTC). None means neither source can
/// answer; no polynomial Delta-T model or extrapolation is substituted.
/// Unlike deltat(), this returns seconds. Rows are acquired lazily under
/// current cache/automatic-download/network settings. Nonfinite epochs raise
/// ValueError; directory/configuration/acquisition errors can propagate.
#[pyfunction]
fn get_delta_t_iers(py: Python<'_>, jd: f64) -> PyResult<Option<f64>> {
    py.detach(|| runtime().iers_rows()?.delta_t_seconds(jd))
        .map_err(|error| to_pyerr(py, error))
}
/// Return UT1-minus-UTC seconds at finite Modified Julian Day mjd, or None.
///
/// mjd is the native UT1 date coordinate (JD - 2400000.5), not a full JD.
/// Daily interpolation uses continuous UT1-minus-TAI across leap seconds and
/// restores TAI-minus-UTC at the query. A single available floor/next-day
/// bracket supplies its value; neither available returns None. There is no
/// remote/model estimate outside those brackets. Lazy row acquisition may
/// create/download cache files. Bad epochs and acquisition errors raise.
#[pyfunction]
fn get_ut1_utc(py: Python<'_>, mjd: f64) -> PyResult<Option<f64>> {
    py.detach(|| {
        let rows = runtime().iers_rows()?;
        iers::ut1_utc_for_mjd(mjd, rows.finals_rows(), rows.leap_seconds())
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Return TAI-minus-UTC seconds at finite UTC Modified Julian Day mjd.
///
/// The pre-1972 rate model is used before the leap era; later epochs use the
/// last applicable leap entry (the final known offset continues afterward,
/// not a prediction of future leap events). Missing leap files use packaged
/// rows. The shared IERS generation is still acquired lazily, potentially
/// creating/downloading files. Nonfinite mjd raises ValueError; acquisition
/// and configuration failures propagate.
#[pyfunction]
fn get_tai_utc(py: Python<'_>, mjd: f64) -> PyResult<f64> {
    py.detach(|| {
        let rows = runtime().iers_rows()?;
        time::tai_utc_for_mjd(mjd, rows.leap_seconds())
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Return the finals rows' (minimum, maximum) UT1 Julian Days, or None if empty.
///
/// This envelope is in full JD, not MJD, and does not prove every daily row
/// exists. A cold query lazily loads rows and can create/download missing
/// cache files; configuration/acquisition failures propagate.
#[pyfunction]
fn get_iers_data_range(py: Python<'_>) -> PyResult<Option<(f64, f64)>> {
    py.detach(|| {
        runtime()
            .iers_rows()
            .map(|rows| iers::iers_data_range_jd(rows.finals_rows()))
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Return whether finals contain a floor/next-day bracket for UT1 JD jd.
///
/// This is a row-availability check, not a monthly/model Delta-T check or
/// simply an envelope test. It lazily loads rows under the current settings.
/// Nonfinite jd raises ValueError; acquisition/configuration errors propagate.
#[pyfunction]
fn is_iers_data_available(py: Python<'_>, jd: f64) -> PyResult<bool> {
    py.detach(|| {
        let rows = runtime().iers_rows()?;
        iers::is_iers_data_available(jd, rows.finals_rows())
    })
    .map_err(|error| to_pyerr(py, error))
}
/// Release the default parsed IERS generation; return None without disk I/O.
///
/// Files and settings are unchanged; the next row query loads again. Active
/// calculations and independently retained snapshots keep their owned rows.
#[pyfunction]
fn clear_iers_cache(py: Python<'_>) {
    py.detach(|| runtime().clear_iers_cache());
}
/// Delete the three selected IERS cache files and clear the default row cache.
///
/// Return the count of successful removals (0 through 3), not attempted
/// removals or bytes. Missing files and individual removal errors are ignored;
/// path/configuration/directory errors can raise. An absent cache directory
/// may be created. Other files, active snapshots and other contexts are not
/// invalidated. Future lazy queries may download these files again.
#[pyfunction]
fn delete_iers_cache_files(py: Python<'_>) -> PyResult<usize> {
    py.detach(|| runtime().delete_iers_cache_files())
        .map_err(|error| to_pyerr(py, error))
}
/// Return a copied dictionary of disk status and loaded IERS row metadata.
///
/// Keys: cache_dir (str); finals/leap_seconds/delta_t suffixed with _exists
/// (bool) and _age_days (one-decimal days or None); data_points,
/// leap_second_entries, delta_t_entries (counts); data_range_start_jd/end_jd
/// and delta_t_range_start_jd/end_jd (UT1 JD or None).
///
/// Counts are sampled before lazy range acquisition: a first call may report
/// zero counts but nonempty ranges. The query may create the directory and
/// acquire missing files under policy. Presence/ages are not digest checks;
/// future modification times have negative ages. Configuration/acquisition
/// errors propagate. Editing the returned dictionary changes no cache state.
#[pyfunction]
fn get_iers_cache_info<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
    let value = py
        .detach(|| runtime().iers_cache_info())
        .map_err(|error| to_pyerr(py, error))?;
    let result = PyDict::new(py);
    result.set_item("cache_dir", value.cache_dir.to_string_lossy().as_ref())?;
    result.set_item("finals_exists", value.finals.exists)?;
    result.set_item("finals_age_days", value.finals.age_days)?;
    result.set_item("leap_seconds_exists", value.leap_seconds.exists)?;
    result.set_item("leap_seconds_age_days", value.leap_seconds.age_days)?;
    result.set_item("delta_t_exists", value.delta_t.exists)?;
    result.set_item("delta_t_age_days", value.delta_t.age_days)?;
    result.set_item("data_points", value.data_points)?;
    result.set_item("leap_second_entries", value.leap_second_entries)?;
    result.set_item("delta_t_entries", value.delta_t_entries)?;
    result.set_item("data_range_start_jd", value.data_range.map(|range| range.0))?;
    result.set_item("data_range_end_jd", value.data_range.map(|range| range.1))?;
    result.set_item(
        "delta_t_range_start_jd",
        value.delta_t_range.map(|range| range.0),
    )?;
    result.set_item(
        "delta_t_range_end_jd",
        value.delta_t_range.map(|range| range.1),
    )?;
    Ok(result)
}
/// Select a str/PathLike IERS cache directory; None restores configured discovery.
///
/// Clears the default parsed generation so the next lookup uses this choice;
/// creates/opens/deletes no files now. Existing row snapshots remain usable.
#[pyfunction(signature=(path))]
fn set_iers_cache_dir(py: Python<'_>, path: Option<PathBuf>) {
    py.detach(|| runtime().set_iers_cache_dir(path));
}
/// Return the explicit IERS cache-directory override as str, or None.
///
/// This does not resolve the environment/default directory or create it.
/// Use get_iers_cache_info() for the effective directory (with lazy I/O).
#[pyfunction]
fn get_iers_cache_dir(py: Python<'_>) -> Option<String> {
    py.detach(|| {
        runtime()
            .iers_cache_dir_override()
            .map(|path| path.to_string_lossy().into_owned())
    })
}
/// Override automatic missing-file acquisition with True/False; None restores
/// captured/default configuration. Return None; no download or cache clear
/// happens now. Enabling this does not bypass a sealed network policy or
/// automatically enable IERS Delta-T in calculations.
#[pyfunction(signature=(enabled))]
fn set_iers_auto_download(py: Python<'_>, enabled: Option<bool>) {
    py.detach(|| runtime().set_iers_auto_download(enabled));
}
/// Return the resolved automatic IERS acquisition toggle, independently of
/// network permission and IERS Delta-T use. No files are read or downloaded;
/// invalid captured configuration can raise an error.
#[pyfunction]
fn get_iers_auto_download(py: Python<'_>) -> PyResult<bool> {
    py.detach(|| runtime().iers_auto_download())
        .map_err(|error| to_pyerr(py, error))
}
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<IERSDataPoint>()?;
    module.add_class::<LeapSecondEntry>()?;
    module.add_class::<DeltaTDataPoint>()?;
    module.add_function(wrap_pyfunction!(download_iers_finals, module)?)?;
    module.add_function(wrap_pyfunction!(download_leap_seconds, module)?)?;
    module.add_function(wrap_pyfunction!(download_delta_t_data, module)?)?;
    module.add_function(wrap_pyfunction!(load_iers_data, module)?)?;
    module.add_function(wrap_pyfunction!(get_observed_delta_t, module)?)?;
    module.add_function(wrap_pyfunction!(get_observed_delta_t_data_range, module)?)?;
    module.add_function(wrap_pyfunction!(is_observed_delta_t_available, module)?)?;
    module.add_function(wrap_pyfunction!(get_delta_t_iers, module)?)?;
    module.add_function(wrap_pyfunction!(get_ut1_utc, module)?)?;
    module.add_function(wrap_pyfunction!(get_tai_utc, module)?)?;
    module.add_function(wrap_pyfunction!(get_iers_data_range, module)?)?;
    module.add_function(wrap_pyfunction!(is_iers_data_available, module)?)?;
    module.add_function(wrap_pyfunction!(clear_iers_cache, module)?)?;
    module.add_function(wrap_pyfunction!(delete_iers_cache_files, module)?)?;
    module.add_function(wrap_pyfunction!(get_iers_cache_info, module)?)?;
    module.add_function(wrap_pyfunction!(set_iers_cache_dir, module)?)?;
    module.add_function(wrap_pyfunction!(get_iers_cache_dir, module)?)?;
    module.add_function(wrap_pyfunction!(set_iers_auto_download, module)?)?;
    module.add_function(wrap_pyfunction!(get_iers_auto_download, module)?)?;
    Ok(())
}
