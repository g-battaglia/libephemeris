//! Independent native contexts and their owned data/clock handles.

use crate::errors::to_pyerr;
use crate::positions::{PositionResult, position, runtime};
use crate::state::LebReader;
use libephemeris_core::config::Runtime;
use libephemeris_core::positions::MainImageOwner;
use pyo3::IntoPyObjectExt;
use pyo3::exceptions::{PyRuntimeError, PyUserWarning, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyTuple, PyType};
use std::collections::HashMap;
use std::ffi::CString;
use std::path::PathBuf;
use std::sync::Arc;

/// A default-runtime reference or an independently owned context. Handles
/// retain the actual Rust runtime rather than constructing a foreign loader.
enum RuntimeHandle {
    Default,
    Owned(Arc<Runtime>),
}
impl std::ops::Deref for RuntimeHandle {
    type Target = Runtime;
    fn deref(&self) -> &Self::Target {
        match self {
            Self::Default => runtime(),
            Self::Owned(value) => value,
        }
    }
}

/// Open planetary SPK images using a retained runtime's paths and policy.
///
/// Obtain this handle from get_loader() or EphemerisContext.get_loader(); it
/// cannot be constructed directly. Calling loader(filename) returns an owned
/// PlanetaryEphemeris without changing the runtime's selected main image.
/// Acquisition may open or download data when the configured policy allows it.
/// Missing files, invalid SPKs, configuration and network-policy errors propagate.
/// A forced LEB configuration refuses SPK access with Error.
///
/// The handle keeps its runtime alive and observes its later configuration
/// changes. Returned images instead retain the data they actually opened.
/// There is no close() method; dropping the handle releases its runtime reference.
///
/// Example (given a context and an available kernel or permitted acquisition):
///
/// ```python
/// loader = context.get_loader()
/// image = loader("de440.bsp")
/// position, velocity = image.geometric_state_tdb(399, 2451545.0)
/// ```
#[pyclass(module = "libephemeris.context", frozen)]
pub(crate) struct DataLoader {
    runtime: RuntimeHandle,
}
impl DataLoader {
    pub(crate) fn default_runtime() -> Self {
        Self {
            runtime: RuntimeHandle::Default,
        }
    }
}

#[pymethods]
impl DataLoader {
    /// Read-only directory string; querying does not create the directory.
    ///
    /// Default-runtime handles report its data root; independent-context handles
    /// report its absolute planetary search directory. Configuration/path errors
    /// propagate. This is a path label, not a separate Python loader object.
    #[getter]
    fn directory(&self, py: Python<'_>) -> PyResult<String> {
        let directory = match &self.runtime {
            RuntimeHandle::Default => self.runtime.data_directory(),
            RuntimeHandle::Owned(_) => self.runtime.library_path(),
        };
        directory
            .map(|path| path.to_string_lossy().into_owned())
            .map_err(|error| to_pyerr(py, error))
    }
    /// Open filename under this runtime's paths/policy; see DataLoader for errors.
    fn __call__(&self, py: Python<'_>, filename: &str) -> PyResult<PlanetaryEphemeris> {
        py.detach(|| self.runtime.open_planetary_image(filename))
            .map(|owner| PlanetaryEphemeris { owner })
            .map_err(|error| to_pyerr(py, error))
    }
}

/// Native clock helpers returning tuples of ordinary Python floats.
///
/// Obtain from get_timescale() or EphemerisContext.get_timescale(); no public
/// constructor or close() method is provided. UT1 conversion uses the independent
/// frame-clock model, not context-specific user/model Delta-T overrides. These
/// results are numbers, not third-party time objects.
///
/// Example (given an EphemerisContext named context):
///
/// ```python
/// ts = context.get_timescale()
/// assert ts.tt_jd(2451545.125) == (2451545.0, 0.125)
/// ut1, offset_days = ts.ut1_jd(2451545.125)
/// # (ut1, offset_days) represents TT; its first part is not floored.
/// ```
#[pyclass(module = "libephemeris.context", frozen)]
pub(crate) struct TimeScale {
    runtime: RuntimeHandle,
}
impl TimeScale {
    pub(crate) fn default_runtime() -> Self {
        Self {
            runtime: RuntimeHandle::Default,
        }
    }
}

#[pymethods]
impl TimeScale {
    /// Convert finite UT1 JD to two-part TT (original_jd, tt_minus_ut1_days).
    ///
    /// The first part is the entire input, not its integer part. Keep both
    /// parts to retain the small clock offset; their sum is scalar TT JD.
    /// Raises ValueError for nonfinite input or a nonfinite model result.
    fn ut1_jd(&self, py: Python<'_>, jd: f64) -> PyResult<(f64, f64)> {
        py.detach(|| self.runtime.frame_tt_parts(jd))
            .map_err(|error| to_pyerr(py, error))
    }
    /// Split finite TT JD into (floor(jd), jd - floor(jd)); no clock conversion.
    ///
    /// The fractional part is in [0, 1). Raises ValueError for nonfinite input.
    fn tt_jd(&self, py: Python<'_>, jd: f64) -> PyResult<(f64, f64)> {
        libephemeris_core::config::Runtime::tt_parts(jd).map_err(|error| to_pyerr(py, error))
    }
    /// Return inward TT JD bounds corresponding to a closed TDB JD interval.
    ///
    /// start/end must be finite with start <= end. The native TT/TDB model is
    /// inverted and representable endpoints are chosen inside the source domain;
    /// this is not a deep-time accuracy guarantee. Returns (start_tt, end_tt).
    /// Raises InputValidationError for invalid or unrepresentable intervals and
    /// ConvergenceError if the endpoint search cannot settle.
    fn tdb_interval_to_tt(&self, py: Python<'_>, start: f64, end: f64) -> PyResult<(f64, f64)> {
        libephemeris_core::time::tt_interval_from_tdb(start, end)
            .map_err(|error| to_pyerr(py, error))
    }
}

pub(crate) fn planetary_image_error(
    py: Python<'_>,
    error: libephemeris_core::errors::Error,
) -> PyErr {
    if matches!(
        error,
        libephemeris_core::errors::Error::ConfigurationMissing {
            missing_config: "jpl_source_access",
            ..
        }
    ) {
        PyRuntimeError::new_err(error.to_string())
    } else {
        to_pyerr(py, error)
    }
}

/// Immutable owned planetary SPK image, independent of later reconfiguration.
///
/// Obtain from get_planets(), a context's get_planets(), or a DataLoader call;
/// no public constructor or close() method is provided. The image remains usable
/// after runtime resource release while this handle is retained. State queries
/// use only its stored center chains, not alternate images or apparent-position
/// corrections. NAIF target codes are distinct from library planet IDs.
///
/// Returned lists/tuples are Python-owned copies, never NumPy arrays or views
/// into the image; modifying them cannot change future state queries.
#[pyclass(module = "libephemeris.context", frozen)]
pub(crate) struct PlanetaryEphemeris {
    owner: MainImageOwner,
}
type GeometricDescriptor = (i32, i32, i32, f64, f64);
impl PlanetaryEphemeris {
    pub(crate) fn from_owner(owner: MainImageOwner) -> Self {
        Self { owner }
    }
}

#[pymethods]
impl PlanetaryEphemeris {
    /// Read-only exact registry path/key identifying the retained image.
    #[getter]
    fn path(&self) -> String {
        self.owner.registry_key().to_string()
    }
    /// Read-only effective ephemeris filename selected during acquisition.
    #[getter]
    fn filename(&self) -> String {
        self.owner.effective_name().to_string()
    }
    /// Return geometric SSB-relative ICRS/J2000 state at finite TDB JD.
    ///
    /// naif_id is a NAIF target in this image. Returns ([x, y, z], [vx, vy, vz])
    /// as a tuple of new lists, in AU and AU per TDB day. No clock conversion,
    /// apparent correction or extrapolation occurs. Every center link must be
    /// present and cover the epoch. Raises ValueError for nonfinite jd_tdb,
    /// DataNotFoundError for missing targets/centers, EphemerisRangeError for
    /// uncovered dates, and SPKEvaluationError for invalid supported-source data.
    fn geometric_state_tdb(
        &self,
        py: Python<'_>,
        naif_id: i32,
        jd_tdb: f64,
    ) -> PyResult<([f64; 3], [f64; 3])> {
        py.detach(|| self.owner.geometric_state_tdb(naif_id, jd_tdb))
            .map_err(|error| to_pyerr(py, error))
    }
    /// Return new [x, y, z, vx, vy, vz] lists at the supplied TDB JDs.
    ///
    /// naif_id and physical scope/errors follow geometric_state_tdb(). jd_tdb
    /// is a sequence of finite dates; output is a list of rows in input order,
    /// with AU positions and AU per TDB day velocities. Empty input returns [].
    /// Any failed sample raises its error rather than returning a partial batch.
    fn geometric_states_tdb(
        &self,
        py: Python<'_>,
        naif_id: i32,
        jd_tdb: Vec<f64>,
    ) -> PyResult<Vec<[f64; 6]>> {
        py.detach(|| self.owner.geometric_states_tdb(naif_id, &jd_tdb))
            .map_err(|error| to_pyerr(py, error))
    }
    /// Return geometric SSB ICRS/J2000 rows on a finite TT JD grid.
    ///
    /// naif_id selects the retained NAIF target; jd_tt is a sequence of dates.
    /// The native split TT/TDB clock converts epochs and scales velocities to
    /// AU per TT day. Returns new [x, y, z, vx, vy, vz] lists, in input order,
    /// with positions in AU. Empty input returns []; any clock/source error
    /// aborts the batch. Other errors and image scope follow geometric_state_tdb().
    fn geometric_states_tt(
        &self,
        py: Python<'_>,
        naif_id: i32,
        jd_tt: Vec<f64>,
    ) -> PyResult<Vec<[f64; 6]>> {
        py.detach(|| self.owner.geometric_states_tt(naif_id, &jd_tt))
            .map_err(|error| to_pyerr(py, error))
    }
    /// Create a fixed-target GeometricSpkSource retaining these immutable images.
    ///
    /// naif_id must occur in this primary image, or be the SSB origin (0).
    /// supporting is an optional sequence of PlanetaryEphemeris handles supplying
    /// missing center links. Primary descriptors win; later supports take
    /// precedence over earlier supports. No file is opened or downloaded.
    /// The returned source survives deletion of the original Python handles.
    /// Raises DataNotFoundError for missing targets/links and SPKEvaluationError
    /// for conflicting image identities or invalid chain/frame metadata.
    #[pyo3(signature=(naif_id, supporting=None))]
    fn geometric_source(
        &self,
        py: Python<'_>,
        naif_id: i32,
        supporting: Option<Vec<Py<PlanetaryEphemeris>>>,
    ) -> PyResult<GeometricSpkSource> {
        let images: Vec<_> = supporting
            .unwrap_or_default()
            .iter()
            .map(|image| image.borrow(py).owner.clone())
            .collect();
        let borrowed = images.iter().collect::<Vec<_>>();
        py.detach(|| self.owner.geometric_source(naif_id, &borrowed))
            .map(|inner| GeometricSpkSource { inner })
            .map_err(|error| to_pyerr(py, error))
    }
    /// Return file-order (target, center, frame, start_tdb, end_tdb) tuples.
    ///
    /// IDs are NAIF codes; bounds are inclusive TDB JDs. The new list preserves
    /// overlapping/disjoint segments and does not certify a continuous evaluable
    /// center chain. Missing retained image metadata raises SPKNotFoundError.
    fn segment_descriptors_tdb(&self, py: Python<'_>) -> PyResult<Vec<GeometricDescriptor>> {
        py.detach(|| self.owner.segment_descriptors_tdb())
            .map(|rows| {
                rows.into_iter()
                    .map(|row| (row.target, row.center, row.frame, row.start_jd, row.end_jd))
                    .collect()
            })
            .map_err(|error| to_pyerr(py, error))
    }
}

/// Fixed geometric NAIF target with retained primary/supporting SPK images.
///
/// Create with PlanetaryEphemeris.geometric_source(), not a direct constructor.
/// Queries return SSB-relative ICRS/J2000 AU and daily velocities without
/// apparent corrections, source fallback, downloads or extrapolation. Retained
/// images survive runtime close/reconfiguration and original-handle deletion.
/// There is no close() method; data ownership ends with the last retained owner.
#[pyclass(module = "libephemeris.context", frozen, skip_from_py_object)]
pub(crate) struct GeometricSpkSource {
    pub(crate) inner: libephemeris_core::spk::GeometricSpkSource,
}

#[pymethods]
impl GeometricSpkSource {
    /// Read-only fixed NAIF target code, without library planet-ID translation.
    #[getter]
    fn naif_id(&self) -> i32 {
        self.inner.target()
    }
    /// New list of retained image keys in center-chain order; not a mutable view.
    #[getter]
    fn paths(&self) -> Vec<String> {
        self.inner.paths().to_vec()
    }
    /// Return a new list of closed (start, end) TDB JD coverage windows.
    ///
    /// Windows intersect every retained center link and may be disjoint. They
    /// describe metadata, not a guarantee that malformed coefficients evaluate.
    /// The SSB origin has unbounded metadata coverage; queries still require
    /// finite dates.
    fn coverage_tdb(&self) -> Vec<(f64, f64)> {
        self.inner.coverage_tdb().to_vec()
    }
    /// Return ([x, y, z], [vx, vy, vz]) at finite TDB JD, in AU/AU per TDB day.
    ///
    /// The tuple contains two new Python lists in SSB ICRS/J2000 coordinates.
    /// Raises ValueError for nonfinite dates, EphemerisRangeError for uncovered
    /// dates, DataNotFoundError for missing links and SPKEvaluationError for
    /// corrupt/unsupported source descriptors. No partial or alternate state.
    fn state_tdb(&self, py: Python<'_>, jd_tdb: f64) -> PyResult<([f64; 3], [f64; 3])> {
        py.detach(|| self.inner.state_tdb(jd_tdb))
            .map_err(|error| to_pyerr(py, error))
    }
    /// Batch state_tdb() over a sequence of finite TDB JDs in input order.
    ///
    /// Returns new [x, y, z, vx, vy, vz] lists in AU/AU per TDB day, not arrays
    /// or borrowed views. Empty input returns []; an error aborts the whole batch.
    fn states_tdb(&self, py: Python<'_>, jd_tdb: Vec<f64>) -> PyResult<Vec<[f64; 6]>> {
        py.detach(|| self.inner.states_tdb(&jd_tdb))
            .map_err(|error| to_pyerr(py, error))
    }
    /// Return SSB ICRS/J2000 rows at a sequence of finite TT JDs.
    ///
    /// The native split TT/TDB clock and velocity chain rule give AU positions
    /// and AU per TT day velocities. Output is a new list of six-float lists
    /// in input order. Empty input returns []; clock/source failures abort the
    /// batch, with source errors as documented for state_tdb().
    fn states_tt(&self, py: Python<'_>, jd_tt: Vec<f64>) -> PyResult<Vec<[f64; 6]>> {
        py.detach(|| self.inner.states_tt(&jd_tt))
            .map_err(|error| to_pyerr(py, error))
    }
}

/// Immutable scalar descriptor of a physical planet-center offset segment.
///
/// Returned in PlanetCenters.segments or PlanetCenterTarget.segments; no public
/// constructor. Values are copied metadata, not an interpolator or a continuous
/// coverage guarantee. All properties are read-only.
#[pyclass(module = "libephemeris.state", frozen, skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct PlanetCenterSegment {
    metadata: libephemeris_core::config::module_state::CenterSegmentMetadata,
}
#[pymethods]
impl PlanetCenterSegment {
    /// NAIF planet-center target code.
    #[getter]
    fn target(&self) -> i32 {
        self.metadata.target
    }
    /// NAIF origin code, normally the corresponding system barycenter.
    #[getter]
    fn center(&self) -> i32 {
        self.metadata.center
    }
    /// NAIF reference-frame identifier (1 denotes J2000).
    #[getter]
    fn frame(&self) -> i32 {
        self.metadata.frame
    }
    /// Inclusive segment start as a TDB Julian day.
    #[getter]
    fn start_jd(&self) -> f64 {
        self.metadata.start_jd
    }
    /// Inclusive segment end as a TDB Julian day.
    #[getter]
    fn end_jd(&self) -> f64 {
        self.metadata.end_jd
    }
}
/// Metadata handle retaining immutable optional planet-center SPK images.
///
/// Returned by get_planet_centers() when images are available; no public
/// constructor or close() method. Handles survive runtime close/reload.
/// Metadata access does not acquire the planetary main image or calculate states.
#[pyclass(module = "libephemeris.state", frozen)]
pub(crate) struct PlanetCenters {
    image: libephemeris_core::config::module_state::PlanetCenterImages,
}
impl PlanetCenters {
    pub(crate) fn from_image(
        image: libephemeris_core::config::module_state::PlanetCenterImages,
    ) -> Self {
        Self { image }
    }
}
#[pymethods]
impl PlanetCenters {
    /// Read-only property returning a new list of actual retained path labels.
    #[getter]
    fn paths(&self) -> Vec<String> {
        self.image.paths()
    }
    /// New list of immutable descriptors, preserving file/segment order.
    ///
    /// TDB bounds belong to individual descriptors, not a continuous union.
    /// Metadata-conversion errors propagate as native input/data errors.
    #[getter]
    fn segments(&self, py: Python<'_>) -> PyResult<Vec<PlanetCenterSegment>> {
        py.detach(|| self.image.segments())
            .map(|values| {
                values
                    .into_iter()
                    .map(|metadata| PlanetCenterSegment { metadata })
                    .collect()
            })
            .map_err(|error| to_pyerr(py, error))
    }
}
/// One NAIF center target's descriptors and retained immutable image owner.
///
/// Returned by get_planet_center_segment(naif_id, jd=None); no public constructor.
/// An undated query retains all matching descriptors; a TDB-dated query selects
/// the last covering descriptor. Handles survive runtime close/reload and have
/// no close() method. Properties return scalars or independent list copies.
#[pyclass(module = "libephemeris.state", frozen)]
pub(crate) struct PlanetCenterTarget {
    target: libephemeris_core::config::module_state::PlanetCenterTarget,
}
impl PlanetCenterTarget {
    pub(crate) fn from_target(
        target: libephemeris_core::config::module_state::PlanetCenterTarget,
    ) -> Self {
        Self { target }
    }
}
#[pymethods]
impl PlanetCenterTarget {
    /// Read-only selected NAIF planet-center target code.
    #[getter]
    fn target(&self) -> i32 {
        self.target.target
    }
    /// New list of retained image path labels, valid after runtime retirement.
    #[getter]
    fn paths(&self) -> Vec<String> {
        self.target.paths()
    }
    /// New list of copied immutable matching descriptors; coverage may be disjoint.
    #[getter]
    fn segments(&self) -> Vec<PlanetCenterSegment> {
        self.target
            .segments
            .iter()
            .cloned()
            .map(|metadata| PlanetCenterSegment { metadata })
            .collect()
    }
}

/// Independent calculation settings with shared immutable astronomical resources.
///
/// Parameters
/// ----------
/// ephe_path : path-like or None
///     Planetary search-directory override; None inherits the default runtime.
/// ephe_file : str
///     Filename override. The default spelling 'de440.bsp' leaves the inherited
///     selection unchanged; it does not force DE440 over an existing selection.
///
/// Construction copies resolved provider/time settings from the process default
/// but starts fresh observer, sidereal, angle-cache and explicit SPK-registration
/// state. It opens no new image. Later global setters do not configure this
/// context. Use separate contexts for independent requests; a sequence of setters
/// and calculations on one shared context is not an atomic transaction.
///
/// Configuration errors can be raised at construction; data/policy failures are
/// normally reported on acquisition/calculation. Returned resource handles own
/// their images. close() releases caches globally without disabling contexts or
/// invalidating retained handles. This is not a 'with' context manager.
///
/// Example:
///
/// ```python
/// from libephemeris.context import EphemerisContext
/// context = EphemerisContext()
/// context.set_topo(12.5, 41.9, 25.0)
/// assert context.topo == (12.5, 41.9, 25.0)
/// # With configured data: values, flags = context.calc_ut(jd_ut1, 0, 258)
/// ```
#[pyclass(module = "libephemeris.context", frozen)]
pub(crate) struct EphemerisContext {
    inner: Arc<Runtime>,
}

#[pymethods]
impl EphemerisContext {
    /// Copy default provider/time settings without opening data; see class docs.
    #[new]
    #[pyo3(signature = (ephe_path=None, ephe_file="de440.bsp"))]
    fn new(py: Python<'_>, ephe_path: Option<PathBuf>, ephe_file: &str) -> PyResult<Self> {
        Runtime::independent_context(runtime(), ephe_path, ephe_file)
            .map(|inner| Self { inner })
            .map_err(|error| to_pyerr(py, error))
    }
    /// Return a loader retaining this runtime and observing its future path/policy changes.
    fn get_loader(&self) -> DataLoader {
        DataLoader {
            runtime: RuntimeHandle::Owned(Arc::clone(&self.inner)),
        }
    }
    /// Return native clock helpers; their frame-clock model ignores Delta-T overrides.
    fn get_timescale(&self) -> TimeScale {
        TimeScale {
            runtime: RuntimeHandle::Owned(Arc::clone(&self.inner)),
        }
    }
    /// Acquire the selected planetary SPK as an immutable owned handle.
    ///
    /// Later configuration/close calls do not invalidate the returned image.
    /// Raises RuntimeError in forced LEB mode before SPK access; missing files,
    /// malformed data, configuration and NetworkSealedError failures propagate.
    fn get_planets(&self, py: Python<'_>) -> PyResult<PlanetaryEphemeris> {
        py.detach(|| self.inner.main_image())
            .map(PlanetaryEphemeris::from_owner)
            .map_err(|error| planetary_image_error(py, error))
    }
    /// Select a LEB path lazily, or clear the explicit selection with None.
    ///
    /// No file is opened/validated here; existing handles retain their images.
    /// Returns None. Acquisition failures are reported by subsequent consumers.
    #[pyo3(signature = (filepath))]
    fn set_leb_file(&self, filepath: Option<PathBuf>) {
        self.inner.set_leb_file(filepath);
    }
    /// Return an owned reader for the active LEB image, or None when unavailable.
    ///
    /// Auto mode can return None for unavailable/invalid optional data; forced
    /// LEB mode propagates file errors. Resource-policy errors always propagate.
    /// Closing the returned handle does not close other retained reader handles.
    fn get_leb_reader(&self, py: Python<'_>) -> PyResult<Option<LebReader>> {
        py.detach(|| self.inner.leb_reader())
            .map(|file| file.map(LebReader::from_file))
            .map_err(|error| to_pyerr(py, error))
    }
    /// Set this context's geodetic observer; return None on success.
    ///
    /// lon is degrees east in [-180, 360], lat degrees north in [-90, 90],
    /// and alt finite metres above the reference surface (negative allowed).
    /// Raises CoordinateError for invalid/nonfinite components without replacing
    /// the old observer. Ordinary body positions use it with FLG_TOPOCTR;
    /// angle identifiers can use the observer without that flag.
    #[pyo3(signature = (lon, lat, alt=0.0))]
    fn set_topo(&self, py: Python<'_>, lon: f64, lat: f64, alt: f64) -> PyResult<()> {
        self.inner
            .set_topo(lon, lat, alt)
            .map_err(|error| to_pyerr(py, error))
    }
    /// Return (longitude_deg_east, latitude_deg_north, altitude_m), or None if unset.
    fn get_topo(&self) -> Option<(f64, f64, f64)> {
        self.inner
            .context()
            .observer()
            .map(|observer| (observer.lon_deg(), observer.lat_deg(), observer.height_m()))
    }
    /// Read-only observer tuple in degrees east/north and metres, or None if unset.
    #[getter(topo)]
    fn topo_value(&self) -> Option<(f64, f64, f64)> {
        self.get_topo()
    }
    /// Store a SIDM mode plus optional SIDBIT flags and defining epoch/offset.
    ///
    /// t0 is JD TT, or JD UT with SIDBIT_USER_UT; ayan_t0 is degrees. These
    /// values are stored as given, without finite/range validation, and are used
    /// by user-defined modes. Implemented flags are retained; unknown high bits
    /// are stripped with UserWarning. Negative modes remain stored for dispatch.
    /// Returns None. Warning filters may raise after the selection was stored.
    #[pyo3(signature = (mode, t0=0.0, ayan_t0=0.0))]
    fn set_sid_mode(&self, py: Python<'_>, mode: i32, t0: f64, ayan_t0: f64) -> PyResult<()> {
        for warning in self.inner.set_sid_mode(mode, t0, ayan_t0) {
            let message = CString::new(warning.to_string())
                .map_err(|error| PyValueError::new_err(error.to_string()))?;
            PyErr::warn(py, &py.get_type::<PyUserWarning>(), &message, 1)?;
        }
        Ok(())
    }
    /// Return the base mode integer, or (mode, t0, ayan_t0) when full=True.
    ///
    /// The tuple excludes independently stored sidereal_bits. Epoch/offset units
    /// follow set_sid_mode(). Fresh contexts read mode 0 with epoch 2451545.0.
    #[pyo3(signature = (full=false))]
    fn get_sid_mode(&self, py: Python<'_>, full: bool) -> PyResult<Py<PyAny>> {
        let context = self.inner.context();
        if full {
            context.sidereal_mode_full().into_py_any(py)
        } else {
            context.sidereal_mode().into_py_any(py)
        }
    }
    /// Base SIDM ID; assignment preserves epoch/offset but replaces SIDBIT flags.
    ///
    /// Supply any desired flags in the assigned mode. Advisories are discarded
    /// by property assignment; use set_sid_mode() to receive UserWarning.
    #[getter]
    fn sidereal_mode(&self) -> i32 {
        self.inner.context().sidereal_mode()
    }
    /// Set mode/encoded bits with the existing epoch pair, without emitting advisories.
    #[setter]
    fn set_sidereal_mode(&self, mode: i32) {
        let (_, t0, ayan) = self.inner.context().sidereal_mode_full();
        self.inner.set_sid_mode(mode, t0, ayan);
    }
    /// Defining JD, TT unless SIDBIT_USER_UT is set; stored without validation.
    ///
    /// Assignment preserves base mode/offset but clears existing sidereal_bits
    /// and discards advisories. Set the desired bits again after assignment.
    #[getter]
    fn sidereal_t0(&self) -> f64 {
        self.inner.context().sidereal_mode_full().1
    }
    /// Store the epoch unchanged and clear sidereal bits; no advisory is emitted.
    #[setter]
    fn set_sidereal_t0(&self, t0: f64) {
        let (mode, _, ayan) = self.inner.context().sidereal_mode_full();
        self.inner.set_sid_mode(mode, t0, ayan);
    }
    /// Defining ayanamsha in degrees, stored without finite/range validation.
    ///
    /// Assignment preserves base mode/epoch but clears existing sidereal_bits
    /// and discards advisories. It does not calculate an ayanamsha.
    #[getter]
    fn sidereal_ayan_t0(&self) -> f64 {
        self.inner.context().sidereal_mode_full().2
    }
    /// Store the angular offset unchanged and clear sidereal bits without advisories.
    #[setter]
    fn set_sidereal_ayan_t0(&self, ayan: f64) {
        let (mode, t0, _) = self.inner.context().sidereal_mode_full();
        self.inner.set_sid_mode(mode, t0, ayan);
    }
    /// Independent implemented SIDBIT mask; assignment silently discards unknown bits.
    ///
    /// Setting this property preserves base mode/epoch/offset and validates no
    /// numerical values. get_sid_mode(full=True) does not include this mask.
    #[getter]
    fn sidereal_bits(&self) -> i32 {
        self.inner.context().sidereal_bits()
    }
    /// Set only the implemented projection/epoch bits, without emitting warnings.
    #[setter]
    fn set_sidereal_bits(&self, bits: i32) {
        self.inner.set_sidereal_bits(bits);
    }
    /// Return a new dict of cached angle names and scalar values, normally degrees.
    ///
    /// The cache can contain caller-supplied values; reading does not calculate
    /// or validate angles. Editing this copy does not mutate the context.
    fn get_angles_cache(&self) -> HashMap<String, f64> {
        self.inner.context().angles_cache().clone()
    }
    /// Replace the cache with a copy of a string-to-float mapping; return None.
    ///
    /// No angular normalization or physical validation is performed. Later
    /// changes to the supplied Python dict do not affect the retained values.
    fn set_angles_cache(&self, angles: HashMap<String, f64>) {
        self.inner.set_angles_cache(&angles);
    }
    /// Clear this context's angle cache without changing observer/sidereal settings.
    fn clear_angles_cache(&self) {
        self.inner.clear_angles_cache();
    }
    /// Register a library body ID (ipl) against an immutable SPK and NAIF target.
    ///
    /// spk_file accepts a path; unresolved relative paths use this context's
    /// planetary directory. naif_id is an integer NAIF code, not a library ID.
    /// Returns None. Raises SPKNotFoundError for missing files, DataNotFoundError
    /// for absent targets and SPKEvaluationError for corrupt supported SPK data.
    /// Existing calculation snapshots retain their prior image ownership.
    fn register_spk_body(
        &self,
        py: Python<'_>,
        ipl: i32,
        spk_file: PathBuf,
        naif_id: i32,
    ) -> PyResult<()> {
        py.detach(|| self.inner.register_spk_body(ipl, &spk_file, naif_id))
            .map_err(|error| to_pyerr(py, error))
    }
    /// Remove a context-local library body registration; absent IDs are a no-op.
    ///
    /// Returns None and does not invalidate images retained by existing calls.
    fn unregister_spk_body(&self, ipl: i32) {
        self.inner.unregister_spk_body(ipl);
    }
    /// Return a copied (path, naif_id) tuple for library ID ipl, or None.
    ///
    /// This metadata lookup does not remove the registration or reopen its file.
    fn get_spk_body_info(&self, ipl: i32) -> Option<(String, i32)> {
        self.inner.registered_spks().remove(&ipl)
    }
    /// Return a new dict mapping library IDs to (path, NAIF ID) tuples.
    ///
    /// Registrations remain listed after resource release for lazy reopening;
    /// mutating the returned dict does not change this context.
    fn list_spk_bodies(&self) -> HashMap<i32, (String, i32)> {
        self.inner.registered_spks()
    }
    /// Calculate library body ipl at finite UT1 Julian day tjd_ut.
    ///
    /// iflag combines FLG_* constants. Returns ((a, b, c, da, db, dc), flags):
    /// normally longitude/latitude in degrees, distance in AU and rates per TT day.
    /// EQUATORIAL selects RA/Dec, XYZ Cartesian AU, RADIANS angular radians;
    /// J2000, origin, apparent and sidereal choices follow flags. Rates require
    /// SPEED/SPEED3. Returned flags are not a physical data-source identity.
    /// ECL_NUT returns four frame angles in degrees and two zeros. Angle and
    /// Arabic-part IDs instead return degree longitude plus five placeholders,
    /// ignoring ordinary XYZ/RADIANS shaping. Angles consume this context's
    /// observer even without TOPOCTR; parts use its cached chart inputs rather
    /// than recomputing them at the requested epoch.
    /// Raises ValueError for nonfinite time, UnknownBodyError for unknown IDs,
    /// EphemerisRangeError for uncovered dates, and native data/configuration/
    /// network-policy errors. Physical/native advisories emit Python warnings.
    fn calc_ut(
        &self,
        py: Python<'_>,
        tjd_ut: f64,
        ipl: i32,
        iflag: i32,
    ) -> PyResult<PositionResult> {
        py.detach(|| self.inner.calc_ut(tjd_ut, ipl, iflag))
            .map_err(|error| to_pyerr(py, error))
            .and_then(|value| position(py, value))
    }
    /// Calculate library body ipl at finite TT Julian day tjd.
    ///
    /// iflag, owned tuple layout, units, warnings and errors follow calc_ut();
    /// only the input clock differs. TT is not UTC or UT1.
    fn calc(&self, py: Python<'_>, tjd: f64, ipl: i32, iflag: i32) -> PyResult<PositionResult> {
        py.detach(|| self.inner.calc(tjd, ipl, iflag))
            .map_err(|error| to_pyerr(py, error))
            .and_then(|value| position(py, value))
    }
    /// Calculate body ipl relative to body iplctr at UT1 Julian day tjd_ut.
    ///
    /// IMPORTANT: this context method takes UT1; the top-level calc_pctr()
    /// function takes TT. Both IDs are library body IDs, not NAIF codes. This
    /// context's Delta-T model converts UT1 to TT. iflag and the owned
    /// (six-component tuple, flags) layout/units follow calc_ut(). Unlike an
    /// ordinary position request, SPEED3 alone leaves the three rates zero;
    /// use FLG_SPEED to request them. Unsupported flag/body combinations, invalid times and native data/range errors raise
    /// their public exceptions; no substitute center is selected.
    fn calc_pctr(
        &self,
        py: Python<'_>,
        tjd_ut: f64,
        ipl: i32,
        iplctr: i32,
        iflag: i32,
    ) -> PyResult<PositionResult> {
        py.detach(|| self.inner.calc_pctr_ut(tjd_ut, ipl, iplctr, iflag))
            .map_err(|error| to_pyerr(py, error))
            .and_then(|value| position(py, value))
    }
    /// Return (cusps, angles) at UT1 JD and geographic degrees north/east.
    ///
    /// tjd_ut is UT1 JD; lat/lon are latitude/longitude, independent of set_topo.
    /// hsys is an integer character code, e.g. ord('P'), not a string. cusps is
    /// a 12-float tuple (index 0 is house 1). angles is an 8-float tuple ordered
    /// Asc, MC, ARMC, Vertex, equatorial Asc, co-Asc 1, co-Asc 2, polar Asc.
    /// Values are degrees; no implicit polar fallback is used. Raises ValueError
    /// if hsys does not fit a byte. Unknown alphabetic selectors use Placidus
    /// geometry; G/g and nonalphabetic bytes raise a library Error because this
    /// fixed-twelve-cusp entry does not serve Gauquelin sectors. Invalid sites
    /// raise CoordinateError and undefined polar geometry can raise
    /// PolarCircleError; clock/data errors also propagate.
    fn houses<'py>(
        &self,
        py: Python<'py>,
        tjd_ut: f64,
        lat: f64,
        lon: f64,
        hsys: i32,
    ) -> PyResult<Bound<'py, PyTuple>> {
        let system = u8::try_from(hsys)
            .map_err(|_| PyValueError::new_err("house system must fit one byte"))?;
        let result = py
            .detach(|| self.inner.houses(tjd_ut, lat, lon, system))
            .map_err(|error| to_pyerr(py, error))?;
        PyTuple::new(
            py,
            [
                PyTuple::new(py, result.cusps)?,
                PyTuple::new(py, result.angles)?,
            ],
        )
    }
    /// Release cached resources of all shared contexts and the default runtime.
    ///
    /// This classmethod is global even when called through an instance. It
    /// preserves selected settings and registration definitions; contexts remain
    /// reusable and reopen data lazily. Retained image/reader handles and active
    /// calculation snapshots stay valid. It is not the package-level close(),
    /// which also resets transient default-runtime state. Returns None.
    #[classmethod]
    fn close(_cls: &Bound<'_, PyType>, py: Python<'_>) {
        py.detach(|| {
            Runtime::release_shared_resources();
            runtime().release_resources();
        });
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<EphemerisContext>()?;
    module.add_class::<DataLoader>()?;
    module.add_class::<TimeScale>()?;
    module.add_class::<PlanetaryEphemeris>()?;
    module.add_class::<GeometricSpkSource>()?;
    module.add_class::<PlanetCenterSegment>()?;
    module.add_class::<PlanetCenters>()?;
    module.add_class::<PlanetCenterTarget>()?;
    Ok(())
}
