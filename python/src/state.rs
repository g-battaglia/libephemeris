//! Process-wide state forwarding and native immutable metadata records.

use crate::context::{
    DataLoader, PlanetCenterTarget, PlanetCenters, PlanetaryEphemeris, TimeScale,
};
use libephemeris_core::context::{TIER_TABLE, Tier, spk_date_range_for_tier};
use libephemeris_core::leb::{LebFile, StarEntry};
use pyo3::IntoPyObjectExt;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Mutex;

use crate::positions::runtime;

/// Immutable precision-tier metadata, not a data loader or coverage guarantee.
///
/// Built-in records from list_tiers()/TIERS describe base, medium and extended.
/// Direct construction stores arbitrary strings without validating a tier or
/// changing global settings. All properties are read-only; equality compares
/// the stored fields. No filesystem or network work occurs.
#[pyclass(
    module = "libephemeris.state",
    frozen,
    get_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct PrecisionTier {
    /// Tier label; built-ins are "base", "medium" and "extended".
    name: String,
    /// Main JPL kernel filename, not an opened path.
    ephemeris_file: String,
    /// (start, end) civil-date strings for minor-body SPK provisioning.
    spk_date_range: (String, String),
    /// Human-readable intended range/size description, not a validity test.
    description: String,
}

#[pymethods]
impl PrecisionTier {
    /// Store name, ephemeris_file, (start, end) date strings and description.
    /// All arguments are required; no data/tier validation or activation occurs.
    #[new]
    fn new(
        name: String,
        ephemeris_file: String,
        spk_date_range: (String, String),
        description: String,
    ) -> Self {
        Self {
            name,
            ephemeris_file,
            spk_date_range,
            description,
        }
    }
}

/// An immutable copy of one selected LEB file's stored star record.
///
/// The identifier belongs to the file; it is not interpreted as a HIP number
/// or looked up in the built-in catalog. Values retain their stored units and
/// IEEE values, without physical validation, propagation or frame conversion.
/// Obtain it from LebReader.get_star(); there is no public constructor.
/// All properties are read-only copies. The record retains no file resource
/// and remains usable after reader closure.
#[pyclass(module = "libephemeris.state", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct LebStarRecord {
    inner: StarEntry,
}

#[pymethods]
impl LebStarRecord {
    /// Internal identifier stored in this LEB file, without catalog remapping.
    #[getter]
    fn star_id(&self) -> i32 {
        self.inner.star_id()
    }
    /// J2000 right ascension, degrees.
    #[getter]
    fn ra_j2000_deg(&self) -> f64 {
        self.inner.ra_j2000_deg()
    }
    /// J2000 declination, degrees.
    #[getter]
    fn dec_j2000_deg(&self) -> f64 {
        self.inner.dec_j2000_deg()
    }
    /// RA proper motion including cos(dec), degrees per year.
    #[getter]
    fn pm_ra_deg_per_yr(&self) -> f64 {
        self.inner.pm_ra_deg_per_yr()
    }
    /// Declination proper motion, degrees per year.
    #[getter]
    fn pm_dec_deg_per_yr(&self) -> f64 {
        self.inner.pm_dec_deg_per_yr()
    }
    /// Parallax, arcseconds.
    #[getter]
    fn parallax_arcsec(&self) -> f64 {
        self.inner.parallax_arcsec()
    }
    /// Radial velocity, kilometres per second.
    #[getter]
    fn rv_km_per_s(&self) -> f64 {
        self.inner.rv_km_per_s()
    }
    /// Stored brightness in magnitudes; no photometric correction or catalog substitution.
    #[getter]
    fn magnitude(&self) -> f64 {
        self.inner.magnitude()
    }
}

/// Owned handle to an immutable native LEB1/LEB2 image and shared coefficient cache.
///
/// Obtain from get_leb_reader(); there is no public constructor. Queries use
/// only this selected file, without companion/model fallback, network requests
/// or apparent-place corrections. Indexed backing may read coefficients lazily.
/// close() releases this handle only; other snapshots and in-flight evaluations
/// remain valid. cool() instead clears shared coefficients without closing.
#[pyclass(module = "libephemeris._core", name = "LebReader", frozen)]
pub(crate) struct LebReader {
    file: Mutex<LebFile>,
}

impl LebReader {
    pub(crate) fn from_file(file: LebFile) -> Self {
        Self {
            file: Mutex::new(file),
        }
    }
    pub(crate) fn image(&self) -> LebFile {
        self.file
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .clone()
    }
}

#[pymethods]
impl LebReader {
    /// Original path label as str, still readable after close(); no file is opened.
    #[getter]
    fn path(&self) -> String {
        self.image().path().to_string()
    }
    /// Header-declared (start, end) TT Julian Days, including endpoints.
    /// This is not recomputed from body coverage. Closed handles raise an error.
    #[getter]
    fn jd_range(&self, py: Python<'_>) -> PyResult<(f64, f64)> {
        self.image()
            .jd_range()
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Whether this handle has been closed; other image owners are independent.
    #[getter]
    fn closed(&self) -> bool {
        self.image().is_closed()
    }
    /// Test for library body_id in this file's index, not epoch usability.
    /// Return bool; a closed handle raises a reader error. No fallback or I/O.
    fn has_body(&self, py: Python<'_>, body_id: i32) -> PyResult<bool> {
        self.image()
            .has_body(body_id)
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Return a copied list of library body IDs in stored index order.
    /// Closed handles raise a reader error; mutating the list changes no data.
    fn body_ids(&self, py: Python<'_>) -> PyResult<Vec<i32>> {
        self.image()
            .index()
            .map(|index| index.bodies().iter().map(|body| body.body_id()).collect())
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Evaluate library body_id at finite TT Julian Day jd without extrapolation.
    ///
    /// Return (components, rates), two copied float lists in stored component
    /// order. ICRS channels use Cartesian AU and AU/TT-day (planet center or
    /// system barycenter as stored); ecliptic channels use longitude/latitude
    /// degrees, distance AU and their TT-day rates. No center/frame/apparent
    /// correction is added. Component count belongs to the file, not a fixed
    /// six-slot calc() result. Body absence, coverage, decoding/backing and
    /// closed-handle errors propagate; nonfinite epochs raise ValueError.
    fn eval_body(&self, py: Python<'_>, body_id: i32, jd: f64) -> PyResult<(Vec<f64>, Vec<f64>)> {
        let image = self.image();
        py.detach(|| image.eval_body_any(body_id, jd))
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Probe the stored nutation header; this does not validate coefficients.
    /// Missing, undecodable or empty headers return False. Backing acquisition
    /// failures and closed handles raise their native reader errors.
    fn has_nutation(&self, py: Python<'_>) -> PyResult<bool> {
        let image = self.image();
        py.detach(|| image.has_nutation())
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Evaluate stored (dpsi, deps) in radians at TT Julian Day jd.
    /// The stored interval includes both endpoints; uncovered or nonfinite
    /// epochs, missing/malformed data and closed handles raise reader errors.
    /// This evaluates stored coefficients, without an analytical fallback.
    fn eval_nutation(&self, py: Python<'_>, jd: f64) -> PyResult<(f64, f64)> {
        let image = self.image();
        py.detach(|| image.eval_nutation(jd))
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Interpolate stored TT-minus-UT1 in days at table Julian Day jd.
    /// Finite epochs outside the table clamp to its nearest endpoint. This
    /// uses the stored date coordinate without clock conversion or a model
    /// fallback. Nonfinite epochs, missing/malformed data, backing failures
    /// and closed handles raise reader errors.
    fn delta_t(&self, py: Python<'_>, jd: f64) -> PyResult<f64> {
        let image = self.image();
        py.detach(|| image.delta_t(jd))
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Copy a stored star record by its file-local integer identifier.
    /// Missing catalogs or IDs raise DataNotFoundError; backing and closed
    /// handle errors propagate. Repeated IDs select the last complete record;
    /// incomplete trailing records are ignored by the native decoder. No
    /// built-in catalog, astrometric transformation or physical validation is used.
    fn get_star(&self, py: Python<'_>, star_id: i32) -> PyResult<LebStarRecord> {
        let image = self.image();
        py.detach(|| image.get_star(star_id))
            .map(|inner| LebStarRecord { inner })
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Advise warming byte ranges overlapping the TT JD window; return their count.
    ///
    /// jd_start/jd_end must be finite. Owned/static backing touches pages;
    /// mapped backing receives best-effort OS advice; indexed backing does not
    /// prefetch coefficients. No coefficient evaluation/residency guarantee.
    /// Metadata/backing/closed errors propagate, and bad epochs raise ValueError.
    fn warm(&self, py: Python<'_>, jd_start: f64, jd_end: f64) -> PyResult<usize> {
        let image = self.image();
        py.detach(|| image.warm(jd_start, jd_end))
            .map_err(|error| crate::errors::to_pyerr(py, error))
    }
    /// Clear shared decoded coefficients and advise mapped-page release; return None.
    /// Other clones share that coefficient cache, but active evaluations remain
    /// valid. Owned image bytes/directories remain; no immediate RSS reduction
    /// is guaranteed. Safe after close() or repeated cooling. No network/I/O
    /// acquisition is performed; later reads may decode coefficients again.
    fn cool(&self, py: Python<'_>) {
        py.detach(|| self.image().cool());
    }
    /// Release this handle's reference to its image; repeated closes are a no-op.
    /// Return None without closing other snapshots or invalidating active calls.
    /// Numerical/index queries subsequently fail; path and closed remain readable.
    fn close(&self) {
        self.file
            .lock()
            .unwrap_or_else(|error| error.into_inner())
            .close();
    }
}

/// Return a DataLoader handle bound to the mutable process-default runtime.
///
/// Constructing this handle performs no I/O. Calling it with a filename may
/// open/download an SPK under the then-current paths and policy, and returns
/// an independently retained immutable image rather than a third-party object.
#[pyfunction]
fn get_loader() -> DataLoader {
    DataLoader::default_runtime()
}
/// Return native TimeScale helpers bound to the process-default runtime.
///
/// Handle creation acquires no files. Methods return ordinary float tuples,
/// not third-party time objects; UT1-to-TT uses the independent bundled clock,
/// not the user's model Delta-T override. See TimeScale for split-JD layouts.
#[pyfunction]
fn get_timescale() -> TimeScale {
    TimeScale::default_runtime()
}
/// Eagerly acquire the selected main JPL kernel as a PlanetaryEphemeris handle.
///
/// May open/download data under current paths/policy; a cached immutable owner
/// is reused. The returned image survives later close()/reload/reconfiguration.
/// Forced LEB mode raises RuntimeError before SPK access. File, configuration,
/// validation and NetworkSealedError failures propagate. No position is computed.
#[pyfunction]
fn get_planets(py: Python<'_>) -> PyResult<PlanetaryEphemeris> {
    py.detach(|| runtime().main_image())
        .map(PlanetaryEphemeris::from_owner)
        .map_err(|error| crate::context::planetary_image_error(py, error))
}
/// Return retained optional PlanetCenters images, or None when unavailable.
///
/// Forced LEB returns None. Other modes discover/open and digest-check local
/// tier companions lazily; they do not acquire/download the main kernel.
/// Missing or rejected optional images are skipped (rejections are logged);
/// configuration/retention errors can propagate. Returned handles remain valid
/// after runtime close/reload and do not themselves prove chain coverage.
#[pyfunction]
fn get_planet_centers(py: Python<'_>) -> PyResult<Option<PlanetCenters>> {
    py.detach(|| runtime().planet_center_images())
        .map(|image| image.map(PlanetCenters::from_image))
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Select optional center-offset descriptors for a NAIF target, or return None.
///
/// naif_id is a NAIF planet-center code, not a library planet ID. jd=None
/// retains all matching descriptors; a finite TDB Julian Day selects the last
/// descriptor covering that epoch inclusively. This metadata operation neither
/// interpolates nor invents continuous coverage between segments. It may open
/// local optional center images like get_planet_centers(), never the main
/// kernel or network. When images exist, nonfinite jd raises ValueError;
/// acquisition/configuration errors propagate. The result retains its image.
#[pyfunction(signature=(naif_id,jd=None))]
fn get_planet_center_segment(
    py: Python<'_>,
    naif_id: i32,
    jd: Option<f64>,
) -> PyResult<Option<PlanetCenterTarget>> {
    py.detach(|| {
        runtime().planet_center_images().and_then(|image| {
            image
                .map(|image| image.segment(naif_id, jd))
                .transpose()
                .map(Option::flatten)
        })
    })
    .map(|target| target.map(PlanetCenterTarget::from_target))
    .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Return the selected (longitude_east_deg, latitude_north_deg, height_m), or None.
/// This copies process-default observer settings; it calculates/acquires nothing.
#[pyfunction]
fn get_topo() -> Option<(f64, f64, f64)> {
    runtime()
        .context()
        .observer()
        .map(|observer| (observer.lon_deg(), observer.lat_deg(), observer.height_m()))
}
/// Return the sidereal base mode ID, or (mode, t0, ayan_t0) when full=True.
///
/// t0 is the user defining TT Julian Day unless SIDBIT_USER_UT selects UT1;
/// ayan_t0 is degrees. This reports stored selection, not an ayanamsha value
/// or whether a later calculation requests SIDEREAL. Unset/tropical readback
/// uses the Fagan/Bradley mode with zero epoch/offset. No data acquisition.
#[pyfunction(signature=(full=false))]
fn get_sid_mode(py: Python<'_>, full: bool) -> PyResult<Py<PyAny>> {
    let context = runtime().context();
    if full {
        context.sidereal_mode_full().into_py_any(py)
    } else {
        context.sidereal_mode().into_py_any(py)
    }
}
/// Return a copied {angle_name: longitude_deg} dictionary for the default session.
/// No angle is calculated and mutating this copy does not alter cached values.
#[pyfunction]
fn get_angles_cache() -> HashMap<String, f64> {
    runtime().context().angles_cache().clone()
}
/// Replace cached derived angles from a {str: float} mapping; return None.
/// Values are conventionally longitudes in degrees, copied without physical
/// validation or recomputation. Later edits to the supplied dict have no effect.
#[pyfunction]
fn set_angles_cache(angles: HashMap<String, f64>) {
    runtime().set_angles_cache(&angles);
}
/// Clear cached angle results, retaining observer/sidereal settings and data.
/// Return None; this performs no computation, filesystem or network work.
#[pyfunction]
fn clear_angles_cache() {
    runtime().clear_angles_cache();
}
/// Return whether a missing unrecognized JPL filename selects the tier default.
///
/// This diagnoses configured filename fallback, not every provider fallback
/// or successful data acquisition. Forced LEB returns False. It checks local
/// presence/configuration without opening files, creating directories or
/// downloading. Invalid captured configuration can raise an error.
#[pyfunction]
fn jpl_fallback_active(py: Python<'_>) -> PyResult<bool> {
    py.detach(|| runtime().jpl_fallback_active())
        .map_err(|error| crate::errors::to_pyerr(py, error))
}

/// Set the process-default geodetic observer; return None.
///
/// lon is east-positive degrees in [-180, 360], lat north-positive degrees in
/// [-90, 90], alt finite metres above the ellipsoid (default 0; negative allowed).
/// This selects a site; topocentric output still requires the relevant flag.
/// Invalid/nonfinite coordinates raise CoordinateError without changing the
/// site. No files or network are acquired.
#[pyfunction(signature = (lon, lat, alt=0.0))]
fn set_topo(py: Python<'_>, lon: f64, lat: f64, alt: f64) -> PyResult<()> {
    runtime()
        .set_topo(lon, lat, alt)
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Set the main-kernel search directory; None or an empty path restores discovery.
///
/// An unchanged path preserves resources. A changed path retires default
/// main/reader/center discovery while active snapshots keep their images;
/// existing explicit SPK mappings are retained. A supplied directory is also
/// scanned eagerly for known local minor-body SPKs (discovery failures are
/// ignored here). No directory is created and no download is initiated.
#[pyfunction(signature = (path=None))]
fn set_ephe_path(py: Python<'_>, path: Option<PathBuf>) {
    py.detach(|| runtime().set_ephe_path(path.filter(|path| !path.as_os_str().is_empty())));
}
/// Select an explicit main JPL filename string for future acquisition; return None.
/// No file is opened or validated now. Changed selection retires the current
/// main association, not retained images. Missing/unrecognized names follow
/// the native filename-selection policy when eventually consumed.
#[pyfunction]
fn set_ephemeris_file(filename: &str) {
    runtime().set_ephemeris_file(filename);
}
/// Select the explicit main JPL filename; same behavior as set_ephemeris_file().
/// Takes str, returns None, and defers opening/validation until consumption.
/// It does not register a minor-body SPK or change the precision tier.
#[pyfunction]
fn set_jpl_file(filename: &str) {
    runtime().set_ephemeris_file(filename);
}
/// Store lunar tidal acceleration acc in arcseconds per Julian century squared.
/// TIDAL_AUTOMATIC restores automatic selection. Other floats are stored
/// without validation; consuming calculations retain their own domains.
/// Return None; no data acquisition or Delta-T evaluation occurs.
#[pyfunction]
fn set_tid_acc(acc: f64) {
    runtime().set_tid_acc(acc);
}
/// Return the stored tidal-acceleration override, or TIDAL_DEFAULT if automatic.
/// Units are arcseconds per Julian century squared. This is a setting readback,
/// not a date-dependent tidal/Delta-T calculation; no data acquisition occurs.
#[pyfunction]
fn get_tid_acc() -> f64 {
    runtime().context().effective_tidal_acceleration()
}
/// Return the explicit TT-minus-UT1 override in days, or None for automatic.
/// Unlike IERS table helpers, this is days rather than seconds. No model
/// evaluation or table acquisition occurs; independent frame clocks are separate.
#[pyfunction]
fn get_delta_t_userdef() -> Option<f64> {
    runtime().context().user_delta_t_days()
}
/// Store atmospheric lapse rate lrate in K/m; None restores 0.0065 K/m.
/// Values are stored without finite/physical validation. Return None without
/// I/O; individual refraction calculations still enforce their input domains.
#[pyfunction(signature = (lrate))]
fn set_lapse_rate(lrate: Option<f64>) {
    runtime().set_lapse_rate(lrate);
}
/// Return the effective atmospheric lapse rate in K/m (default 0.0065).
/// This reads the default session's override without calculation or I/O.
#[pyfunction]
fn get_lapse_rate() -> f64 {
    runtime().context().effective_lapse_rate()
}
/// Return the absolute provider data-directory string without creating it.
/// Explicit main-kernel search directory takes priority over configured data
/// root. Path/configuration failures propagate; no library binary path or
/// guarantee of installed ephemeris files is implied.
#[pyfunction]
fn get_library_path(py: Python<'_>) -> PyResult<String> {
    runtime()
        .library_path()
        .map(|path| path.to_string_lossy().into_owned())
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Return (path_str, start_jd, end_jd, source_de_number) for Sun (0) or Moon (1).
///
/// Other ifno values or no selected metadata return ("", 0.0, 0.0, 0).
/// Stored LEB body bounds are TT JD; an already acquired main SPK reports its
/// TDB segment envelope. The DE number may be 0 when provenance is unknown.
/// Coverage is metadata, not target-chain proof. The query can acquire local
/// configured LEB readers, but does not open/download an unacquired main SPK.
/// Configuration/reader errors propagate; no position is computed.
#[pyfunction(signature = (ifno=0))]
fn get_current_file_data(py: Python<'_>, ifno: i32) -> PyResult<(String, f64, f64, i32)> {
    py.detach(|| runtime().current_file_data(ifno))
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Release default runtime resources and reset transient observing state.
///
/// Clears observer/sidereal/angles, SPK/moon registrations, acquired IERS/main/
/// center/LEB owners and the default Horizons-client reference. Resets explicit
/// tier/JPL directory/filename, SPK cache/download/padding, tidal/lapse and
/// IERS-Delta-T selections to their automatic/default settings. Captured
/// configuration may apply again. Preserve calculation mode, explicit LEB
/// selection, model/user Delta-T, strictness, IERS cache/auto-download and
/// network overrides. No disk files are deleted.
///
/// Return None; not a terminal shutdown. Later consumers acquire again, and
/// active calculations/returned handles keep their immutable owned resources.
/// Use reset_session() for an observing reset that retains acquired data.
#[pyfunction]
fn close() {
    runtime().close();
}
/// Clear observer, sidereal selection and cached angles; return None.
///
/// Retain tier, paths, calculation/clock/network settings and acquired data
/// images/registrations. The runtime drops its default Horizons-client
/// reference, but external client handles/in-flight requests survive. No files
/// are deleted or reopened; this is not close() or reload_configuration().
#[pyfunction]
fn reset_session() {
    runtime().reset_session();
}
/// Return a retained LebReader for the first active LEB image, or None.
///
/// Configured LEB files are opened only in auto/LEB modes. In auto, an unavailable
/// or invalid optional container can return None; forced LEB propagates file
/// errors. Resource-policy/configuration errors always propagate. Acquisition
/// opens/checks local data lazily without downloading. The returned handle has
/// an independent close lifecycle and may share coefficients with calculations.
#[pyfunction]
fn get_leb_reader(py: Python<'_>) -> PyResult<Option<LebReader>> {
    py.detach(|| runtime().leb_reader())
        .map(|reader| reader.map(LebReader::from_file))
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Return the resolved mode: "auto", "leb", "skyfield" or "horizons".
/// "skyfield" names the native Rust JPL/SPK route, not a Python dependency.
/// This resolves settings without acquiring data; invalid configuration raises.
#[pyfunction]
fn get_calc_mode(py: Python<'_>) -> PyResult<String> {
    runtime()
        .effective_context()
        .map(|ctx| ctx.mode().as_calc_mode_name().to_string())
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Override automatic minor-body SPK acquisition; None restores configured policy.
/// True permits the ordinary provisioning path, not an immediate download or
/// bypass of sealed-network policy. False disables that automatic leg only.
/// Return None; existing registrations and files are not removed.
#[pyfunction(signature = (enabled))]
fn set_auto_spk_download(enabled: Option<bool>) {
    runtime().set_auto_spk_download(enabled);
}
/// Return the resolved automatic-SPK toggle, independently of network permission.
/// Configuration errors can raise; no files or downloads are acquired by this query.
#[pyfunction]
fn get_auto_spk_download(py: Python<'_>) -> PyResult<bool> {
    runtime()
        .effective_context()
        .map(|ctx| ctx.effective_auto_spk_download())
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Select a str/PathLike SPK cache directory; None restores captured discovery.
/// Return None without creating/opening files or relocating existing kernels.
/// Registrations retain the images and paths with which they were acquired.
#[pyfunction(signature = (path))]
fn set_spk_cache_dir(path: Option<PathBuf>) {
    runtime().set_spk_cache_dir(path);
}
/// Return the resolved configured SPK cache directory as str, or None.
/// None lets individual acquisition helpers choose their documented default
/// root; this does not create the directory. Configuration errors propagate.
#[pyfunction]
fn get_spk_cache_dir(py: Python<'_>) -> PyResult<Option<String>> {
    runtime()
        .effective_context()
        .map(|ctx| {
            ctx.spk_cache_dir()
                .map(|path| path.to_string_lossy().into_owned())
        })
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Set nonnegative integer padding days around automatic SPK coverage requests.
/// Negative values raise ValueError without changing the previous value.
/// Return None; existing kernels are neither extended nor replaced now.
#[pyfunction]
fn set_spk_date_padding(py: Python<'_>, days: i64) -> PyResult<()> {
    runtime()
        .set_spk_date_padding(days)
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Return configured automatic-SPK padding as an integer number of days.
/// This setting is not proof that any downloaded kernel covers that margin.
#[pyfunction]
fn get_spk_date_padding() -> u64 {
    runtime().context().spk_date_padding_days()
}
/// Return whether the resolved model Delta-T chain may use IERS observations.
/// This does not acquire rows, establish their epoch coverage or change the
/// independent frame/observer clocks. Invalid configuration can raise.
#[pyfunction]
fn get_iers_delta_t_enabled(py: Python<'_>) -> PyResult<bool> {
    runtime()
        .effective_context()
        .map(|ctx| ctx.effective_iers_delta_t())
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Return the resolved Delta-T model name ("smh2016" or "espenak_meeus").
/// The default is "smh2016".
/// This is model selection, not an evaluated TT-minus-UT1 value; user overrides
/// and enabled IERS observations may take priority. No tables are acquired;
/// invalid captured configuration propagates as an error.
#[pyfunction]
fn get_delta_t_model(py: Python<'_>) -> PyResult<String> {
    runtime()
        .effective_context()
        .map(|ctx| {
            ctx.delta_t_model_override()
                .unwrap_or(libephemeris_core::time::DeltaTModel::Smh2016)
                .as_name()
                .to_string()
        })
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Override strict minor-body precision; None restores captured/default settings.
/// True refuses the low-precision Keplerian fallback where the provider chain
/// requires an SPK-quality solution; it does not acquire data or expand coverage.
/// Return None without changing existing kernel files/registrations.
#[pyfunction(signature = (enabled))]
fn set_strict_precision(enabled: Option<bool>) {
    runtime().set_strict_precision(enabled);
}
/// Return the resolved strict minor-body precision toggle without acquisition.
/// This reports policy, not whether a particular body/date can be served.
/// Invalid captured configuration can raise an error.
#[pyfunction]
fn get_strict_precision(py: Python<'_>) -> PyResult<bool> {
    runtime()
        .effective_context()
        .map(|ctx| ctx.effective_strict_precision())
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Return the resolved tier name: "base", "medium" or "extended".
/// This reads configuration without acquiring data or asserting actual body
/// coverage. Other invalid captured configuration can still raise an error.
#[pyfunction]
fn get_precision_tier(py: Python<'_>) -> PyResult<String> {
    runtime()
        .effective_context()
        .map(|ctx| ctx.tier().as_name().to_string())
        .map_err(|error| crate::errors::to_pyerr(py, error))
}
/// Return fresh immutable PrecisionTier records in base/medium/extended order.
/// Kernel names and SPK civil-date windows are planning metadata, not installed
/// coverage checks. No files are opened and the active tier is unchanged.
#[pyfunction]
fn list_tiers() -> Vec<PrecisionTier> {
    TIER_TABLE
        .iter()
        .map(|row| PrecisionTier {
            name: row.name.to_string(),
            ephemeris_file: row.ephemeris_file.to_string(),
            spk_date_range: (row.spk_start.to_string(), row.spk_end.to_string()),
            description: row.description.to_string(),
        })
        .collect()
}
/// Return the tier's planned SPK (start, end) civil-date strings.
/// tier_name=None uses the resolved active tier; otherwise supply exactly
/// "base", "medium" or "extended". Unknown names raise ValueError and
/// configuration errors propagate. This performs no file/coverage inspection,
/// clock conversion or download and does not authorize extrapolation.
#[pyfunction(signature = (tier_name=None))]
fn get_spk_date_range_for_tier(
    py: Python<'_>,
    tier_name: Option<&str>,
) -> PyResult<(String, String)> {
    let tier = match tier_name {
        Some(name) => Tier::from_name(name).ok_or_else(|| {
            crate::errors::to_pyerr(
                py,
                libephemeris_core::errors::Error::InvalidTier {
                    value: name.to_string(),
                },
            )
        })?,
        None => runtime()
            .effective_context()
            .map_err(|error| crate::errors::to_pyerr(py, error))?
            .tier(),
    };
    let (start, end) = spk_date_range_for_tier(tier);
    Ok((start.to_string(), end.to_string()))
}
/// Best-effort advice that OS data pages may be reclaimed; return None.
/// This neither closes readers nor clears decoded coefficient/numerical caches,
/// deletes files, resets settings or guarantees an immediate RSS reduction.
/// Immutable owners remain valid. Resolving the data path can raise a
/// configuration/filesystem error; no network request is made.
#[pyfunction]
fn release_data_cache(py: Python<'_>) -> PyResult<()> {
    py.detach(|| runtime().release_data_cache())
        .map_err(|error| crate::errors::to_pyerr(py, error))
}

/// Cool every acquired LEB image without opening or hashing new readers.
/// Returns the number of images visited; OS page release remains advisory.
/// Pass clear_coefficients=False to preserve the bounded warm coefficient cache.
#[pyfunction(signature = (clear_coefficients=true))]
fn cool_leb_readers(py: Python<'_>, clear_coefficients: bool) -> usize {
    py.detach(|| runtime().cool_leb_readers(clear_coefficients))
}
/// Recapture environment, dotenv and TOML settings for future default calls.
/// Explicit API overrides still win. Retire default reader/center discovery,
/// while active snapshots keep their immutable images/settings. Return None;
/// configuration files are read but process environment is never modified and
/// no download starts. Unreadable files are ignored; invalid resolved settings
/// surface at their consumers. This is not reset_session() or a file deletion.
#[pyfunction]
fn reload_configuration() {
    runtime().reload_configuration();
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(get_loader, module)?)?;
    module.add_function(wrap_pyfunction!(get_timescale, module)?)?;
    module.add_function(wrap_pyfunction!(get_planets, module)?)?;
    module.add_function(wrap_pyfunction!(get_planet_centers, module)?)?;
    module.add_function(wrap_pyfunction!(get_planet_center_segment, module)?)?;
    module.add_function(wrap_pyfunction!(get_topo, module)?)?;
    module.add_function(wrap_pyfunction!(get_sid_mode, module)?)?;
    module.add_function(wrap_pyfunction!(get_angles_cache, module)?)?;
    module.add_function(wrap_pyfunction!(set_angles_cache, module)?)?;
    module.add_function(wrap_pyfunction!(clear_angles_cache, module)?)?;
    module.add_function(wrap_pyfunction!(jpl_fallback_active, module)?)?;
    module.add_class::<PrecisionTier>()?;
    module.add_class::<LebReader>()?;
    module.add("_LebStarRecord", module.py().get_type::<LebStarRecord>())?;
    let tiers = PyDict::new(module.py());
    for tier in list_tiers() {
        tiers.set_item(tier.name.clone(), Py::new(module.py(), tier)?)?;
    }
    module.add("TIERS", tiers)?;
    module.add_function(wrap_pyfunction!(set_topo, module)?)?;
    module.add_function(wrap_pyfunction!(set_ephe_path, module)?)?;
    module.add_function(wrap_pyfunction!(set_ephemeris_file, module)?)?;
    module.add_function(wrap_pyfunction!(set_jpl_file, module)?)?;
    module.add_function(wrap_pyfunction!(set_tid_acc, module)?)?;
    module.add_function(wrap_pyfunction!(get_tid_acc, module)?)?;
    module.add_function(wrap_pyfunction!(get_delta_t_userdef, module)?)?;
    module.add_function(wrap_pyfunction!(set_lapse_rate, module)?)?;
    module.add_function(wrap_pyfunction!(get_lapse_rate, module)?)?;
    module.add_function(wrap_pyfunction!(get_library_path, module)?)?;
    module.add_function(wrap_pyfunction!(get_current_file_data, module)?)?;
    module.add_function(wrap_pyfunction!(close, module)?)?;
    module.add_function(wrap_pyfunction!(reset_session, module)?)?;
    module.add_function(wrap_pyfunction!(get_leb_reader, module)?)?;
    module.add_function(wrap_pyfunction!(get_calc_mode, module)?)?;
    module.add_function(wrap_pyfunction!(set_auto_spk_download, module)?)?;
    module.add_function(wrap_pyfunction!(get_auto_spk_download, module)?)?;
    module.add_function(wrap_pyfunction!(set_spk_cache_dir, module)?)?;
    module.add_function(wrap_pyfunction!(get_spk_cache_dir, module)?)?;
    module.add_function(wrap_pyfunction!(set_spk_date_padding, module)?)?;
    module.add_function(wrap_pyfunction!(get_spk_date_padding, module)?)?;
    module.add_function(wrap_pyfunction!(get_iers_delta_t_enabled, module)?)?;
    module.add_function(wrap_pyfunction!(get_delta_t_model, module)?)?;
    module.add_function(wrap_pyfunction!(set_strict_precision, module)?)?;
    module.add_function(wrap_pyfunction!(get_strict_precision, module)?)?;
    module.add_function(wrap_pyfunction!(get_precision_tier, module)?)?;
    module.add_function(wrap_pyfunction!(list_tiers, module)?)?;
    module.add_function(wrap_pyfunction!(get_spk_date_range_for_tier, module)?)?;
    module.add_function(wrap_pyfunction!(release_data_cache, module)?)?;
    module.add_function(wrap_pyfunction!(cool_leb_readers, module)?)?;
    module.add_function(wrap_pyfunction!(reload_configuration, module)?)?;
    Ok(())
}
