// SPDX-License-Identifier: AGPL-3.0
//! Typed N-body records and direct calls to the native propagation API.

use crate::errors::to_pyerr;
use libephemeris_core::{minor_bodies, nbody, python_support};
use pyo3::exceptions::{PyRuntimeWarning, PyUserWarning, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyDict;
use std::ffi::CString;
use std::path::PathBuf;
use std::sync::Arc;

/// Mutable owned minor-body osculating elements, in the ecliptic-J2000 model.
///
/// Parameters/fields: name is a string; epoch is TT JD; a is AU (perihelion
/// distance q for parabolic rows); e is dimensionless; i, omega, Omega, M0 are
/// degrees; n is stored mean motion in degrees/day, not derived from a.
/// Construction raises ValueError for nonfinite numeric fields but does not
/// reject physically implausible finite elements. Assignment performs ordinary
/// scalar conversion without finiteness/physical checks; validate() rechecks
/// finiteness and returns physical-sanity messages rather than emitting warnings.
/// Propagation functions take an owned value copy, not a mutable live view.
///
/// Import from libephemeris.minor_bodies: the package-root OrbitalElements is a
/// different hypothetical-body record. Example:
///
/// ```python
/// from libephemeris.minor_bodies import OrbitalElements
/// elements = OrbitalElements('example', 2451545.0, 2.5, 0.1,
///                            5.0, 40.0, 80.0, 0.0, 0.25)
/// issues = elements.validate()
/// ```
#[pyclass(module = "libephemeris.minor_bodies", from_py_object)]
#[derive(Clone, Debug)]
#[allow(non_snake_case)]
pub(crate) struct OrbitalElements {
    /// Owned body-name string; assignment does not resolve a catalog identity.
    #[pyo3(get, set)]
    name: String,
    /// Osculating TT Julian day; assignment stores the number without validation.
    #[pyo3(get, set)]
    epoch: f64,
    /// Semi-major axis in AU; q for parabolic rows, negative for hyperbolic rows.
    #[pyo3(get, set)]
    a: f64,
    /// Dimensionless eccentricity; assignment does not enforce an orbit domain.
    #[pyo3(get, set)]
    e: f64,
    /// Inclination in degrees; physical-sanity range is [0, 180].
    #[pyo3(get, set)]
    i: f64,
    /// Argument of perihelion in degrees, distinct from longitude of perihelion.
    #[pyo3(get, set)]
    omega: f64,
    /// Longitude of ascending node in degrees (capital Omega).
    #[pyo3(get, set)]
    Omega: f64,
    /// Mean anomaly at epoch, in degrees.
    #[pyo3(get, set)]
    M0: f64,
    /// Stored mean daily motion in degrees/day; not recomputed on assignment.
    #[pyo3(get, set)]
    n: f64,
    // Retained only for this object obtained from native catalogue acquisition.
    // Ordinary scalar construction supplies an anonymous token, even when all
    // visible fields happen to equal a catalogue record.
    identity: python_support::MinorBodyIdentity,
}

impl OrbitalElements {
    /// Convert this library-owned record to checked native elements.
    pub(crate) fn native(
        &self,
    ) -> libephemeris_core::errors::Result<minor_bodies::OrbitalElements> {
        minor_bodies::OrbitalElements::new(
            self.name.clone(),
            self.epoch,
            self.a,
            self.e,
            self.i,
            self.omega,
            self.Omega,
            self.M0,
            self.n,
        )
        .map(|elements| python_support::retain_minor_body_identity(elements, self.identity))
    }

    /// Convert an owned native element record without changing units or values.
    pub(crate) fn from_native(elements: minor_bodies::OrbitalElements) -> Self {
        Self {
            name: elements.name().into(),
            epoch: elements.epoch(),
            a: elements.a(),
            e: elements.e(),
            i: elements.i(),
            omega: elements.omega(),
            Omega: elements.node(),
            M0: elements.m0(),
            n: elements.n(),
            identity: python_support::minor_body_identity(&elements),
        }
    }
}

#[pymethods]
impl OrbitalElements {
    /// Construct finite owned fields in the units listed by OrbitalElements.
    #[new]
    #[allow(non_snake_case, clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        name: String,
        epoch: f64,
        a: f64,
        e: f64,
        i: f64,
        omega: f64,
        Omega: f64,
        M0: f64,
        n: f64,
    ) -> PyResult<Self> {
        minor_bodies::OrbitalElements::new(name, epoch, a, e, i, omega, Omega, M0, n)
            .map(Self::from_native)
            .map_err(|error| to_pyerr(py, error))
    }

    /// Return a new list of physical-sanity messages without modifying the record.
    ///
    /// Checks axis/eccentricity consistency, inclination and elliptic mean motion;
    /// [] means those checks found no issue, not that integration is guaranteed.
    /// Raises ValueError if assignments introduced nonfinite numeric fields.
    fn validate(&self, py: Python<'_>) -> PyResult<Vec<String>> {
        let elements = self.native().map_err(|error| to_pyerr(py, error))?;
        Ok(py.detach(|| elements.validate()))
    }
}

/// Immutable enum-like selection of one of five REBOUND integration methods.
///
/// value must be exactly 'ias15', 'whfast', 'mercurius', 'trace' or 'leapfrog';
/// unknown/case-mismatched strings raise ValueError. The corresponding class
/// constants are IAS15, WHFAST, MERCURIUS, TRACE and LEAPFROG. This is not a
/// Python enum.Enum and does not provide its iteration/indexing protocols.
/// Propagation integrator arguments require this type, not a plain string.
///
/// Example:
///
/// ```python
/// from libephemeris.rebound_integration import ReboundIntegrator
/// method = ReboundIntegrator('ias15')
/// assert method == ReboundIntegrator.IAS15
/// assert method.value == 'ias15' and method.name == 'IAS15'
/// ```
#[pyclass(
    module = "libephemeris.rebound_integration",
    frozen,
    eq,
    from_py_object
)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct ReboundIntegrator {
    native: nbody::ReboundIntegrator,
}

#[pymethods]
impl ReboundIntegrator {
    /// Parse an exact lowercase method string; raises ValueError for unknown names.
    #[new]
    fn new(py: Python<'_>, value: &str) -> PyResult<Self> {
        nbody::ReboundIntegrator::parse(value)
            .map(|native| Self { native })
            .map_err(|error| to_pyerr(py, error))
    }
    /// Adaptive fifteenth-order Gauss-Radau method; scalar propagation default.
    #[classattr]
    const IAS15: Self = Self {
        native: nbody::ReboundIntegrator::Ias15,
    };
    /// Fixed-step Wisdom-Holman method.
    #[classattr]
    const WHFAST: Self = Self {
        native: nbody::ReboundIntegrator::Whfast,
    };
    /// Hybrid close-encounter method.
    #[classattr]
    const MERCURIUS: Self = Self {
        native: nbody::ReboundIntegrator::Mercurius,
    };
    /// Time-reversible hybrid close-encounter method.
    #[classattr]
    const TRACE: Self = Self {
        native: nbody::ReboundIntegrator::Trace,
    };
    /// Fixed-step leapfrog method.
    #[classattr]
    const LEAPFROG: Self = Self {
        native: nbody::ReboundIntegrator::Leapfrog,
    };

    /// Read-only lowercase constructor value, e.g. 'ias15'.
    #[getter]
    fn value(&self) -> &'static str {
        self.native.as_str()
    }
    /// Read-only uppercase constant name, e.g. 'IAS15'.
    #[getter]
    fn name(&self) -> &'static str {
        match self.native {
            nbody::ReboundIntegrator::Ias15 => "IAS15",
            nbody::ReboundIntegrator::Whfast => "WHFAST",
            nbody::ReboundIntegrator::Mercurius => "MERCURIUS",
            nbody::ReboundIntegrator::Trace => "TRACE",
            nbody::ReboundIntegrator::Leapfrog => "LEAPFROG",
        }
    }
    /// Return the enum-like constant name and lowercase value for inspection.
    fn __repr__(&self) -> String {
        format!("<ReboundIntegrator.{}: '{}'>", self.name(), self.value())
    }
    /// Return 'ReboundIntegrator.NAME'.
    fn __str__(&self) -> String {
        format!("ReboundIntegrator.{}", self.name())
    }
}

/// Mutable ASSIST path selection retaining a reusable native ephemeris session.
///
/// planets_file/asteroids_file are optional path strings; data_dir is an optional
/// discovery-directory string. Missing selections are searched under data_dir,
/// ASSIST_DIR, ~/.libephemeris/assist, then ./data. Construction resolves paths
/// using filesystem metadata but does not map, download or certify their coverage.
/// Explicit missing files remain selected and fail when consumed. Ordinary
/// propagation requires planets; asteroid perturbers are optional.
///
/// Properties expose resolved paths. Assigning any property creates a new lazy
/// session; passing None for a file re-enables discovery of that file. Changing
/// data_dir alone does not replace already resolved file selections. Existing
/// calls/sources retain their previous session. Calls sharing one session
/// serialize native ephemeris access. Selected files must remain immutable while
/// retained; there is no close() method. No process-default settings are changed.
///
/// Example:
///
/// ```python
/// from libephemeris.rebound_integration import AssistEphemConfig
/// config = AssistEphemConfig(data_dir='/path/to/assist')
/// # Opening/mapping is deferred until a propagation/source request.
/// ```
#[pyclass(module = "libephemeris.rebound_integration", from_py_object)]
#[derive(Clone)]
pub(crate) struct AssistEphemConfig {
    session: Arc<nbody::AssistSession>,
}

impl AssistEphemConfig {
    fn replace(&mut self, config: nbody::AssistEphemConfig) {
        self.session = Arc::new(nbody::AssistSession::new(config));
    }
}

#[pymethods]
impl AssistEphemConfig {
    /// Resolve optional path strings without mapping or downloading; see class docs.
    #[new]
    #[pyo3(signature = (planets_file=None, asteroids_file=None, data_dir=None))]
    fn new(
        py: Python<'_>,
        planets_file: Option<String>,
        asteroids_file: Option<String>,
        data_dir: Option<String>,
    ) -> Self {
        let config = nbody::AssistEphemConfig {
            planets_file: planets_file.map(PathBuf::from),
            asteroids_file: asteroids_file.map(PathBuf::from),
            data_dir: data_dir.map(PathBuf::from),
        };
        Self {
            session: Arc::new(py.detach(|| nbody::AssistSession::new(config))),
        }
    }
    /// Resolved planet-force path string or None; assignment replaces the session.
    ///
    /// None requests discovery, not forced absence. Missing paths fail on use.
    #[getter]
    fn planets_file(&self) -> Option<String> {
        self.session
            .config()
            .planets_file
            .as_ref()
            .map(|path| path.to_string_lossy().into_owned())
    }
    /// Select/discover the planet file for future calls; active old owners survive.
    #[setter]
    fn set_planets_file(&mut self, value: Option<String>) {
        let mut config = self.session.config().clone();
        config.planets_file = value.map(PathBuf::from);
        self.replace(config);
    }
    /// Resolved massive-asteroid path string or None; optional for ordinary propagation.
    ///
    /// Assignment creates a new lazy session; None requests discovery rather
    /// than disabling a discoverable file. Explicit-source generation requires it.
    #[getter]
    fn asteroids_file(&self) -> Option<String> {
        self.session
            .config()
            .asteroids_file
            .as_ref()
            .map(|path| path.to_string_lossy().into_owned())
    }
    /// Select/discover the asteroid file without altering an already retained session.
    #[setter]
    fn set_asteroids_file(&mut self, value: Option<String>) {
        let mut config = self.session.config().clone();
        config.asteroids_file = value.map(PathBuf::from);
        self.replace(config);
    }
    /// Optional highest-priority discovery directory, after explicit file paths.
    ///
    /// Assignment replaces the session but preserves its already resolved file
    /// paths. Set a file property to None to discover it under the new directory.
    #[getter]
    fn data_dir(&self) -> Option<String> {
        self.session
            .config()
            .data_dir
            .as_ref()
            .map(|path| path.to_string_lossy().into_owned())
    }
    /// Set discovery directory without discarding already selected planet/asteroid paths.
    #[setter]
    fn set_data_dir(&mut self, value: Option<String>) {
        let mut config = self.session.config().clone();
        config.data_dir = value.map(PathBuf::from);
        self.replace(config);
    }

    /// Integrate a curated library body ID to finite TT Julian day jd_tt.
    ///
    /// Returns an owned PropagationResult: heliocentric ecliptic J2000 AU/AU per
    /// day. TT is used as a TDB approximation (periodic difference about 1.7 ms),
    /// unlike NbodySpkSource's explicit clock conversion. A target present in the
    /// selected asteroid file starts from that image at its catalog epoch; other
    /// curated bodies start from elements. No download or numerical fallback.
    /// Raises UnknownBodyError for absent curated elements, ValueError for
    /// nonfinite input, FileNotFoundError for missing force files, RuntimeError
    /// for native failures, and native SPK/range errors for invalid seed data.
    /// Self-force advisories emit UserWarning; native messages emit RuntimeWarning.
    fn propagate_curated(
        &self,
        py: Python<'_>,
        body: i32,
        jd_tt: f64,
    ) -> PyResult<PropagationResult> {
        let session = self.session.clone();
        let outcome = py
            .detach(|| session.propagate_curated(body, jd_tt))
            .map_err(|error| to_pyerr(py, error))?;
        diagnostics(py, outcome).map(PropagationResult::from)
    }

    /// Create an explicit SPK-seeded NbodySpkSource at finite TDB JD anchor_tdb.
    ///
    /// seed must supply a direct heliocentric J2000 target at the anchor (default
    /// J2000). Both selected force files are required; their common continuous
    /// coverage must contain the anchor after a 32-day guard at each end.
    /// Samples the seed once and retains this configuration's force session,
    /// not the original seed image. No integration to other epochs or download
    /// occurs here. The returned source survives later config assignments.
    /// Raises ValueError for invalid anchor/force-domain selection, FileNotFoundError
    /// for missing files, RuntimeError for native failures and native SPK errors
    /// for unsupported seed centers/frames/coverage. Force files must stay immutable.
    #[pyo3(signature = (seed, anchor_tdb=2451545.0))]
    fn nbody_source(
        &self,
        py: Python<'_>,
        seed: &crate::context::GeometricSpkSource,
        anchor_tdb: f64,
    ) -> PyResult<NbodySpkSource> {
        let session = self.session.clone();
        let source = seed.inner.clone();
        py.detach(|| nbody::NbodySpkSource::new(session, &source, anchor_tdb))
            .map(|inner| NbodySpkSource { inner })
            .map_err(|error| to_pyerr(py, error))
    }
}

/// Explicit SPK-seeded barycentric ICRS/J2000 trajectory for offline sampling.
///
/// Obtain with AssistEphemConfig.nbody_source(); no public constructor/close().
/// Retains a sampled initial state and immutable force-file session. The seed
/// image is needed only at construction. Each batch starts independent monotone
/// IAS15 sweeps forward/backward from the anchor; order is restored afterwards.
/// There is no catalog/model fallback or download. This reduced-force trajectory
/// is not a promise to reproduce a full JPL orbit throughout the force-file span.
/// Returned nested Python lists are owned copies, never NumPy arrays or views.
#[pyclass(
    module = "libephemeris.rebound_integration",
    frozen,
    skip_from_py_object
)]
pub(crate) struct NbodySpkSource {
    inner: nbody::NbodySpkSource,
}

#[pymethods]
impl NbodySpkSource {
    /// Read-only TDB Julian day of the fixed initial SPK sample.
    #[getter]
    fn anchor_tdb(&self) -> f64 {
        self.inner.anchor_tdb()
    }

    /// Return inclusive (start, end) TDB JD bounds with the 32-day force-stage guard.
    ///
    /// This continuous window is shared by the force targets; the seed's narrower
    /// coverage does not limit integration because it was sampled at the anchor.
    fn coverage_tdb(&self) -> (f64, f64) {
        self.inner.coverage_tdb()
    }

    /// Integrate a sequence of finite TDB JDs; return new six-float row lists.
    ///
    /// Rows are [x, y, z, vx, vy, vz] in SSB ICRS/J2000 AU and AU per TDB day.
    /// Input order/duplicates are preserved; [] performs no integration. Raises
    /// ValueError for nonfinite/out-of-domain epochs, RuntimeError for native
    /// integration failures and native file/data errors for force-source failures.
    /// No partial batch is returned. Physical/native diagnostics emit UserWarning
    /// or RuntimeWarning and warning filters may themselves raise.
    fn states_tdb(&self, py: Python<'_>, jd_tdb: Vec<f64>) -> PyResult<Vec<[f64; 6]>> {
        let outcome = py
            .detach(|| self.inner.states_tdb(&jd_tdb))
            .map_err(|error| to_pyerr(py, error))?;
        diagnostics(py, outcome)
    }
    /// Integrate a sequence of finite TT JDs using explicit native TT/TDB conversion.
    ///
    /// Returns new SSB ICRS/J2000 [x, y, z, vx, vy, vz] lists, in AU and AU per
    /// TT day after analytic velocity scaling. Unlike ordinary ASSIST propagation,
    /// TT is not treated as TDB here. Ordering, empty inputs, guarded coverage,
    /// warnings and errors follow states_tdb(); clock failures also abort the batch.
    fn states_tt(&self, py: Python<'_>, jd_tt: Vec<f64>) -> PyResult<Vec<[f64; 6]>> {
        let outcome = py
            .detach(|| self.inner.states_tt(&jd_tt))
            .map_err(|error| to_pyerr(py, error))?;
        diagnostics(py, outcome)
    }
}

/// Mutable owned heliocentric ecliptic-J2000 Cartesian state at TT JD.
///
/// x/y/z are AU; vx/vy/vz AU/day; jd_tt is the scalar TT epoch. Construction
/// raises ValueError for any nonfinite component/time. Later field assignment
/// only converts the scalar: it does not validate, rotate or propagate a state.
/// Derived properties read current fields without revalidating them. Returned
/// results do not borrow an integrator, and editing one cannot change its source.
///
/// Example:
///
/// ```python
/// from libephemeris.rebound_integration import PropagationResult
/// state = PropagationResult(1., 0., 0., 0., 0.01, 0., 2451545.)
/// assert state.distance == 1.0
/// assert len(state.to_tuple()) == 6  # Epoch is a separate property.
/// ```
#[pyclass(
    module = "libephemeris.rebound_integration",
    get_all,
    set_all,
    from_py_object
)]
#[derive(Clone, Copy, Debug)]
pub(crate) struct PropagationResult {
    /// Heliocentric ecliptic-J2000 x in AU; mutable without assignment validation.
    x: f64,
    /// Heliocentric ecliptic-J2000 y in AU; mutable without assignment validation.
    y: f64,
    /// Heliocentric ecliptic-J2000 z in AU; mutable without assignment validation.
    z: f64,
    /// Heliocentric ecliptic-J2000 x velocity in AU/day.
    vx: f64,
    /// Heliocentric ecliptic-J2000 y velocity in AU/day.
    vy: f64,
    /// Heliocentric ecliptic-J2000 z velocity in AU/day.
    vz: f64,
    /// Scalar TT Julian day; editing it does not propagate the Cartesian state.
    jd_tt: f64,
}

impl From<nbody::PropagationResult> for PropagationResult {
    fn from(value: nbody::PropagationResult) -> Self {
        Self {
            x: value.x,
            y: value.y,
            z: value.z,
            vx: value.vx,
            vy: value.vy,
            vz: value.vz,
            jd_tt: value.jd_tt,
        }
    }
}

impl PropagationResult {
    fn native(&self) -> nbody::PropagationResult {
        nbody::PropagationResult {
            x: self.x,
            y: self.y,
            z: self.z,
            vx: self.vx,
            vy: self.vy,
            vz: self.vz,
            jd_tt: self.jd_tt,
        }
    }
}

#[pymethods]
impl PropagationResult {
    /// Construct finite heliocentric ecliptic-J2000 components; see class units.
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        py: Python<'_>,
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
        jd_tt: f64,
    ) -> PyResult<Self> {
        nbody::PropagationResult::from_state(
            nbody::CartesianState {
                position: [x, y, z],
                velocity: [vx, vy, vz],
            },
            jd_tt,
        )
        .map(Self::from)
        .map_err(|error| to_pyerr(py, error))
    }
    /// Read-only longitude in degrees on [0, 360) for a finite nonzero position.
    ///
    /// Uses the current Cartesian fields without numerical revalidation.
    #[getter]
    fn ecliptic_lon(&self) -> f64 {
        self.native().ecliptic_lon()
    }
    /// Read-only latitude in degrees on [-90, 90]; zero position gives 0.
    ///
    /// Uses current fields without numerical revalidation after assignment.
    #[getter]
    fn ecliptic_lat(&self) -> f64 {
        self.native().ecliptic_lat()
    }
    /// Read-only Euclidean heliocentric distance in AU, from current x/y/z.
    #[getter]
    fn distance(&self) -> f64 {
        self.native().distance()
    }
    /// Return an independent (x, y, z, vx, vy, vz) tuple in AU and AU/day.
    ///
    /// jd_tt is deliberately excluded; no propagation or validation is performed.
    #[allow(clippy::wrong_self_convention)]
    fn to_tuple(&self) -> (f64, f64, f64, f64, f64, f64) {
        self.native().to_tuple()
    }
}

fn diagnostics<T>(py: Python<'_>, outcome: nbody::PropagationOutcome<T>) -> PyResult<T> {
    for warning in outcome.warnings {
        let message = CString::new(warning.to_string())
            .map_err(|error| PyValueError::new_err(error.to_string()))?;
        PyErr::warn(py, &py.get_type::<PyUserWarning>(), &message, 1)?;
    }
    for message in outcome.messages {
        let message =
            CString::new(message).map_err(|error| PyValueError::new_err(error.to_string()))?;
        PyErr::warn(py, &py.get_type::<PyRuntimeWarning>(), &message, 1)?;
    }
    Ok(outcome.value)
}

fn session(config: Option<AssistEphemConfig>) -> Arc<nbody::AssistSession> {
    config.map_or_else(nbody::default_assist_session, |config| config.session)
}

/// Return whether this extension was built with the bundled REBOUND backend.
///
/// No simulation, file lookup or external Python REBOUND import is performed.
#[pyfunction]
fn check_rebound_available() -> bool {
    nbody::check_rebound_available()
}
/// Return whether this extension includes ASSIST; does not check its data files.
#[pyfunction]
fn check_assist_available() -> bool {
    nbody::check_assist_available()
}
/// Return the bundled REBOUND version string, or None when its backend is disabled.
#[pyfunction]
fn get_rebound_version() -> Option<&'static str> {
    nbody::get_rebound_version()
}
/// Return the bundled ASSIST version string, or None when its backend is disabled.
#[pyfunction]
fn get_assist_version() -> Option<&'static str> {
    nbody::get_assist_version()
}

/// Propagate copied minor-body elements with a massless target and unit-mass Sun.
///
/// elements is libephemeris.minor_bodies.OrbitalElements. jd_start/jd_end are
/// finite TT JDs; the secular element model initializes at jd_start. Returns an
/// owned heliocentric ecliptic-J2000 PropagationResult at jd_end, in AU/AU per day.
/// This Sun-only calculation has no planetary/asteroid force perturbations.
/// integrator must be a ReboundIntegrator object (default IAS15), not 'ias15'.
/// The preserved compatibility signature displays the default as a string;
/// that display is not an accepted explicit string argument.
/// dt is a finite step in days; zero is allowed only for a zero-length interval.
/// WHFAST/LEAPFROG without dt use 365.25*a**1.5/20 days, requiring positive a.
/// A zero interval returns the initialized state without advancing a solver.
/// Raises ValueError for nonfinite/invalid step or element policy, ImportError
/// if the backend is disabled, RuntimeError for native orbit/integration failures.
/// Native messages emit RuntimeWarning; no alternative orbit replaces failure.
///
/// Example, given a minor-body record named elements:
///
/// ```python
/// from libephemeris.rebound_integration import propagate_orbit_rebound, ReboundIntegrator
/// result = propagate_orbit_rebound(elements, elements.epoch,
///                                  elements.epoch + 10.0,
///                                  integrator=ReboundIntegrator.IAS15)
/// ```
#[pyfunction(signature = (elements, jd_start, jd_end, integrator=ReboundIntegrator::IAS15, dt=None),
    text_signature = "(elements, jd_start, jd_end, integrator='ias15', dt=None)")]
fn propagate_orbit_rebound(
    py: Python<'_>,
    elements: OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    integrator: ReboundIntegrator,
    dt: Option<f64>,
) -> PyResult<PropagationResult> {
    let outcome = py
        .detach(|| {
            nbody::propagate_orbit_rebound_with_warnings(
                &elements.native()?,
                jd_start,
                jd_end,
                integrator.native,
                dt,
            )
        })
        .map_err(|error| to_pyerr(py, error))?;
    diagnostics(py, outcome).map(PropagationResult::from)
}

/// Propagate copied minor-body elements with retained JPL-driven ASSIST forces.
///
/// jd_start/jd_end are finite scalar TT JDs, used as a TDB approximation (about
/// 1.7 ms periodic difference); elements initialize at jd_start. Returns an owned
/// heliocentric ecliptic-J2000 state at jd_end in AU and AU/day. ephem_config
/// selects an AssistEphemConfig session, or None uses discovered default data.
/// Planet data are required; asteroid perturbers are optional. No download or
/// Sun-only fallback occurs. Selected files must remain immutable while retained.
/// include_non_gravitational enables Marsden A1/A2/A3 radial/transverse/normal
/// accelerations in AU/day**2; otherwise those three values are ignored.
/// Zero duration returns the initialized state. Raises ValueError for nonfinite
/// admitted inputs, FileNotFoundError for missing force files, ImportError if
/// disabled and RuntimeError for native construction/coverage/integration failures.
/// Self-force exclusions emit UserWarning; native messages emit RuntimeWarning.
#[pyfunction(signature = (elements, jd_start, jd_end, ephem_config=None, include_non_gravitational=false, A1=0.0, A2=0.0, A3=0.0))]
#[allow(non_snake_case, clippy::too_many_arguments)]
fn propagate_orbit_assist(
    py: Python<'_>,
    elements: OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    ephem_config: Option<AssistEphemConfig>,
    include_non_gravitational: bool,
    A1: f64,
    A2: f64,
    A3: f64,
) -> PyResult<PropagationResult> {
    let outcome = py
        .detach(|| {
            session(ephem_config).propagate(
                &elements.native()?,
                jd_start,
                jd_end,
                include_non_gravitational.then_some(nbody::NonGravitationalParameters {
                    radial: A1,
                    transverse: A2,
                    normal: A3,
                }),
            )
        })
        .map_err(|error| to_pyerr(py, error))?;
    diagnostics(py, outcome).map(PropagationResult::from)
}

/// Return a new list of owned states at equally spaced TT JDs.
///
/// elements is a copied minor-body record. Finite jd_start/jd_end include both
/// endpoints for num_points > 1; one gives the start and nonpositive counts [].
/// Reverse intervals are allowed. States are heliocentric ecliptic J2000 AU and
/// AU/day. use_assist=True tries the selected/discovered ephem_config session
/// (TT-as-TDB); unavailable data/backend, force-model/domain or initialization
/// failures allowed by the native policy fall back to Sun-only REBOUND. Native
/// escape/encounter/interruption failures remain errors. No download occurs.
/// integrator must be a ReboundIntegrator object and applies to the Sun-only
/// route; trajectory fixed-step setup uses the native default, not scalar dt.
/// Its compatibility signature displays 'ias15', but explicitly supplying that
/// string is not supported; use ReboundIntegrator.IAS15.
/// Raises ValueError for invalid/nonfinite admitted inputs or unallocatable
/// sample counts, ImportError for an unavailable required backend, RuntimeError
/// for terminal native failures and native file/data exceptions as applicable.
/// Physical/native diagnostics emit UserWarning/RuntimeWarning. No partial list
/// is returned; warning filters may raise instead of returning successful states.
#[pyfunction(signature = (elements, jd_start, jd_end, num_points=100, use_assist=true, ephem_config=None, integrator=ReboundIntegrator::IAS15),
    text_signature = "(elements, jd_start, jd_end, num_points=100, use_assist=True, ephem_config=None, integrator='ias15')")]
#[allow(clippy::too_many_arguments)]
fn propagate_trajectory(
    py: Python<'_>,
    elements: OrbitalElements,
    jd_start: f64,
    jd_end: f64,
    num_points: i64,
    use_assist: bool,
    ephem_config: Option<AssistEphemConfig>,
    integrator: ReboundIntegrator,
) -> PyResult<Vec<PropagationResult>> {
    let outcome = py
        .detach(|| {
            session(ephem_config).trajectory(
                &elements.native()?,
                jd_start,
                jd_end,
                num_points,
                use_assist,
                integrator.native,
            )
        })
        .map_err(|error| to_pyerr(py, error))?;
    Ok(diagnostics(py, outcome)?
        .into_iter()
        .map(PropagationResult::from)
        .collect())
}

fn comparison_point<'py>(
    py: Python<'py>,
    point: nbody::ComparisonPoint,
) -> PyResult<Bound<'py, PyDict>> {
    let result = PyDict::new(py);
    result.set_item("x", point.position[0])?;
    result.set_item("y", point.position[1])?;
    result.set_item("z", point.position[2])?;
    result.set_item("lon", point.longitude)?;
    result.set_item("lat", point.latitude)?;
    result.set_item("dist", point.distance)?;
    Ok(result)
}

/// Return whether ASSIST is built and default discovered force paths exist.
///
/// Planet data are required; a selected asteroid path must also exist. Does not
/// open/map files, certify their format/coverage, download, or inspect a supplied
/// AssistEphemConfig. False is a capability/discovery result, not a raised error.
#[pyfunction]
fn check_assist_data_available(py: Python<'_>) -> bool {
    py.detach(nbody::check_assist_data_available)
}
/// Release the idle default ASSIST session; return None.
///
/// Future default use rediscovers paths and maps lazily. Active calls and explicit
/// AssistEphemConfig/NbodySpkSource owners keep their files; no data is deleted.
#[pyfunction]
fn reset_assist_data_cache(py: Python<'_>) {
    py.detach(nbody::reset_assist_data_cache);
}
/// Return a new dict of secularly perturbed orbit values at finite TT JD.
///
/// elements is a copied minor-body OrbitalElements record. Keys are a (AU),
/// e (dimensionless), inc/omega/Omega/M (radians). Unlike the input record's
/// degree fields, these angles are radians for the native orbit constructor.
/// No integration, resource acquisition or record mutation occurs. Raises
/// ValueError for nonfinite inputs or nonfinite derived orbit values.
#[pyfunction(signature=(elements, jd_tt))]
fn elements_to_rebound_orbit<'py>(
    py: Python<'py>,
    elements: OrbitalElements,
    jd_tt: f64,
) -> PyResult<Bound<'py, PyDict>> {
    let native = elements.native().map_err(|error| to_pyerr(py, error))?;
    let orbit = py
        .detach(|| nbody::elements_to_rebound_orbit(&native, jd_tt))
        .map_err(|error| to_pyerr(py, error))?;
    let result = PyDict::new(py);
    for (key, value) in [
        ("a", orbit.semi_major_axis),
        ("e", orbit.eccentricity),
        ("inc", orbit.inclination),
        ("omega", orbit.perihelion_argument),
        ("Omega", orbit.ascending_node),
        ("M", orbit.mean_anomaly),
    ] {
        result.set_item(key, value)?;
    }
    Ok(result)
}
/// Explicitly provision selected ASSIST files and return a pathlib.Path directory.
///
/// target_dir is path-like or defaults to ~/.libephemeris/assist. planets requests
/// DE440, asteroids requests massive-asteroid perturbers, planets_de441 additionally
/// requests the deep-time DE441 planet file. Verified existing files are retained
/// unless force=True; downloads validate temporary files before replacement.
/// The destination directory may be created even when network policy then refuses
/// acquisition. show_progress/quiet are accepted but do not alter native output.
/// Does not change runtime configuration or existing explicit ASSIST sessions.
/// Raises ConfigurationError if no default directory can be resolved,
/// NetworkSealedError for refused network access, and native data/transport/file
/// exceptions for failed provisioning. No astronomical calculation is performed.
#[pyfunction(signature=(target_dir=None, planets=true, asteroids=true, planets_de441=false, force=false, show_progress=true, quiet=false))]
#[allow(clippy::too_many_arguments)]
fn download_assist_data<'py>(
    py: Python<'py>,
    target_dir: Option<PathBuf>,
    planets: bool,
    asteroids: bool,
    planets_de441: bool,
    force: bool,
    show_progress: bool,
    quiet: bool,
) -> PyResult<Bound<'py, PyAny>> {
    let _ = (show_progress, quiet);
    let path = py
        .detach(|| {
            crate::positions::runtime().download_assist_data(
                target_dir.as_deref(),
                planets,
                asteroids,
                planets_de441,
                force,
            )
        })
        .map_err(|error| to_pyerr(py, error))?;
    py.import("pathlib")?.getattr("Path")?.call1((path,))
}

/// Compare propagated and secular-Keplerian positions at finite TT JD jd_tt.
///
/// elements is a copied minor-body record; integration begins at elements.epoch.
/// use_assist/ephem_config select ASSIST (TT-as-TDB) or Sun-only IAS15. Permitted
/// ASSIST data/model failures fall back as in propagate_trajectory(); terminal
/// native failures propagate. Returns a new dict with keplerian/nbody dicts
/// (x/y/z/dist in heliocentric AU, lon/lat ecliptic-J2000 degrees),
/// delta_lon_arcsec (wrapped), delta_lat_arcsec, angular_sep_arcsec (small-angle
/// cosine-mean-latitude norm), delta_dist_au, method and propagation_days.
/// Differences are nbody minus keplerian. method reports the requested mode,
/// not proof that ASSIST survived fallback. This is not an accuracy certificate.
/// Errors follow the propagation routes; physical/native diagnostics become
/// UserWarning/RuntimeWarning. No input mutation, borrowed arrays or download.
#[pyfunction(signature = (elements, jd_tt, use_assist=false, ephem_config=None))]
fn compare_with_keplerian<'py>(
    py: Python<'py>,
    elements: OrbitalElements,
    jd_tt: f64,
    use_assist: bool,
    ephem_config: Option<AssistEphemConfig>,
) -> PyResult<Bound<'py, PyDict>> {
    let outcome = py
        .detach(|| session(ephem_config).compare(&elements.native()?, jd_tt, use_assist))
        .map_err(|error| to_pyerr(py, error))?;
    let value = diagnostics(py, outcome)?;
    let result = PyDict::new(py);
    result.set_item("keplerian", comparison_point(py, value.keplerian)?)?;
    result.set_item("nbody", comparison_point(py, value.nbody)?)?;
    result.set_item("delta_lon_arcsec", value.delta_lon_arcsec)?;
    result.set_item("delta_lat_arcsec", value.delta_lat_arcsec)?;
    result.set_item("angular_sep_arcsec", value.angular_sep_arcsec)?;
    result.set_item("delta_dist_au", value.delta_dist_au)?;
    result.set_item("method", value.method)?;
    result.set_item("propagation_days", value.propagation_days)?;
    Ok(result)
}

/// Register the public own records, methods and version/capability queries.
pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add(
        "_MinorOrbitalElements",
        module.py().get_type::<OrbitalElements>(),
    )?;
    module.add_class::<ReboundIntegrator>()?;
    module.add_class::<AssistEphemConfig>()?;
    module.add_class::<NbodySpkSource>()?;
    module.add_class::<PropagationResult>()?;
    module.add_function(wrap_pyfunction!(check_rebound_available, module)?)?;
    module.add_function(wrap_pyfunction!(check_assist_available, module)?)?;
    module.add_function(wrap_pyfunction!(get_rebound_version, module)?)?;
    module.add_function(wrap_pyfunction!(get_assist_version, module)?)?;
    module.add_function(wrap_pyfunction!(check_assist_data_available, module)?)?;
    module.add_function(wrap_pyfunction!(reset_assist_data_cache, module)?)?;
    module.add_function(wrap_pyfunction!(elements_to_rebound_orbit, module)?)?;
    module.add_function(wrap_pyfunction!(download_assist_data, module)?)?;
    module.add_function(wrap_pyfunction!(propagate_orbit_rebound, module)?)?;
    module.add_function(wrap_pyfunction!(propagate_orbit_assist, module)?)?;
    module.add_function(wrap_pyfunction!(propagate_trajectory, module)?)?;
    module.add_function(wrap_pyfunction!(compare_with_keplerian, module)?)?;
    Ok(())
}
