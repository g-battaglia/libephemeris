//! Minor-body catalog, SBDB, resonance and native SPK provisioning.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use crate::rebound_integration::OrbitalElements;
use libephemeris_core::errors::Error;
use libephemeris_core::minor_bodies as bodies;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

/// Immutable Neptune mean-motion resonance descriptor, not a dynamical proof.
///
/// ``p:q`` means body/Neptune mean-motion ratio p/q and corresponding axis
/// ``a_Neptune*(q/p)**(2/3)``. Constructor stores values without validation.
#[derive(Clone, Debug)]
#[pyclass(
    module = "libephemeris.minor_bodies",
    frozen,
    get_all,
    skip_from_py_object
)]
struct ResonanceInfo {
    /// Numerator p of the expected body/Neptune mean-motion ratio p/q.
    p: i32,
    /// Denominator q of the expected mean-motion ratio p/q.
    q: i32,
    /// Descriptive resonance label, such as ``plutino``.
    name: String,
    /// Expected heliocentric semimajor axis in AU.
    a_resonant: f64,
    /// Dimensionless fractional mean-motion-ratio detection threshold.
    tolerance: f64,
}

/// Immutable resonance-test outcome with optional match and advisory text.
///
/// A match is only proximity in mean motion, not confirmation that the
/// resonant argument librates. No-match records have zero deviation and
/// ``resonance=warning_message=None``; constructor itself permits any values.
#[derive(Clone, Debug)]
#[pyclass(
    module = "libephemeris.minor_bodies",
    frozen,
    get_all,
    skip_from_py_object
)]
struct ResonanceResult {
    /// Whether a catalogue ratio passed the requested fractional threshold.
    is_resonant: bool,
    /// Closest matching ``ResonanceInfo``, or ``None`` for no match.
    resonance: Option<ResonanceInfo>,
    /// Fractional ratio deviation; zero when unmatched, not an accuracy estimate.
    deviation: f64,
    /// Secular-theory advisory for a match, otherwise ``None``; not emitted as a warning.
    warning_message: Option<String>,
}

#[pymethods]
impl ResonanceInfo {
    /// Store ratio integers, name, expected axis (AU) and fractional tolerance.
    /// No ratio/physical validation or registry insertion is performed.
    #[new]
    fn new(p: i32, q: i32, name: String, a_resonant: f64, tolerance: f64) -> Self {
        Self {
            p,
            q,
            name,
            a_resonant,
            tolerance,
        }
    }
}

#[pymethods]
impl ResonanceResult {
    /// Store a result; ``resonance`` and ``warning_message`` may explicitly be None.
    /// ``deviation`` is a dimensionless fraction. No warning is emitted or state checked.
    #[new]
    #[pyo3(signature = (is_resonant, resonance, deviation, warning_message))]
    fn new(
        is_resonant: bool,
        resonance: Option<PyRef<'_, ResonanceInfo>>,
        deviation: f64,
        warning_message: Option<String>,
    ) -> Self {
        Self {
            is_resonant,
            resonance: resonance.map(|row| row.clone()),
            deviation,
            warning_message,
        }
    }
}

impl From<bodies::ResonanceResult> for ResonanceResult {
    fn from(value: bodies::ResonanceResult) -> Self {
        Self {
            is_resonant: value.is_resonant,
            resonance: value.resonance.map(|row| ResonanceInfo {
                p: row.p,
                q: row.q,
                name: row.name.into(),
                a_resonant: row.a_resonant,
                tolerance: row.tolerance,
            }),
            deviation: value.deviation,
            warning_message: value.warning_message,
        }
    }
}

/// Resolve an asteroid ``name`` to its catalogue number, or ``None`` on a miss.
///
/// Trims/lowercases names for the local/cache lookup, then queries JPL SBDB
/// using the original name. Successful remote results are cached process-wide.
/// ``timeout=30.0`` seconds governs connection/response acquisition. Blank names,
/// unknown responses and ordinary transport/parse failures return ``None``;
/// sealed-network policy and configuration errors raise. Local hits need no socket.
/// Example:
///
/// ```python
/// from libephemeris.minor_bodies import get_asteroid_number
/// number = get_asteroid_number("Ceres")  # local catalogue number 1
/// ```
#[pyfunction(signature = (name, timeout = 30.0))]
fn get_asteroid_number(py: Python<'_>, name: &str, timeout: f64) -> PyResult<Option<i64>> {
    py.detach(|| runtime().get_asteroid_number(name, timeout))
        .map_err(|error| to_pyerr(py, error))
}

/// Clear process-wide asteroid name mappings and return the number removed.
/// The next nonblank lookup can seed built-in local names again; no files are deleted.
#[pyfunction]
fn clear_asteroid_name_cache(py: Python<'_>) -> usize {
    py.detach(|| runtime().clear_asteroid_name_cache())
}

/// Return numbered asteroid ``(longitude_deg, latitude_deg, range_au)`` at TT JD.
///
/// ``asteroid_number`` is a positive catalogue number, not a library body ID.
/// Output is heliocentric true ecliptic of date. Curated bodies use registered
/// SPKs when available, otherwise the approximate multi-epoch secular/Keplerian
/// model; other numbers try SPK when ``use_spk=True``, then SBDB elements, which
/// may require permitted network access. Curated bodies still probe SPK even
/// with ``use_spk=False``. No speed tuple is returned. Invalid/unserved numbers
/// raise ``ValueError``; non-finite time, configuration and sealed-policy errors
/// propagate. Approximate elements are not equivalent to precision SPK/N-body
/// output; off-epoch, resonant and close-encounter accuracy is limited.
#[pyfunction(signature = (asteroid_number, jd_tt, use_spk = true))]
fn calc_asteroid_by_number(
    py: Python<'_>,
    asteroid_number: i64,
    jd_tt: f64,
    use_spk: bool,
) -> PyResult<(f64, f64, f64)> {
    py.detach(|| runtime().calc_asteroid_by_number(asteroid_number, jd_tt, use_spk))
        .map_err(|error| match error {
            Error::InvalidBody { .. } | Error::UnknownBody { .. } => {
                PyValueError::new_err(error.to_string())
            }
            error => to_pyerr(py, error),
        })
}

/// Fetch/cache JPL SBDB osculating ``minor_bodies.OrbitalElements`` by number.
///
/// Result fields use TT epoch JD, axis AU, dimensionless e, degree angles and
/// degree/day n. ``timeout=30.0`` seconds applies to connection/response budgets.
/// Successful cached rows need no socket. Unknown/nonpositive numbers and
/// ordinary transport/parse failures return ``None``; sealed policy can refuse
/// even a nonpositive number. Invalid consumed timeout and configuration/policy
/// errors raise. This record is NOT ``hypothetical.OrbitalElements``.
#[pyfunction(signature = (asteroid_number, timeout = 30.0))]
fn fetch_orbital_elements_from_sbdb(
    py: Python<'_>,
    asteroid_number: i64,
    timeout: f64,
) -> PyResult<Option<OrbitalElements>> {
    py.detach(|| runtime().fetch_orbital_elements_from_sbdb(asteroid_number, timeout))
        .map(|row| row.map(OrbitalElements::from_native))
        .map_err(|error| to_pyerr(py, error))
}

/// Clear process-wide cached SBDB element rows and return the number removed.
/// Curated element tables and downloaded SPK files/registrations are unchanged.
#[pyfunction]
fn clear_asteroid_elements_cache(py: Python<'_>) -> usize {
    py.detach(|| runtime().clear_asteroid_elements_cache())
}

/// Provision/register an asteroid SPK; return its path string or ``None``.
///
/// ``body_id`` is a library ID. ``jd_start``/``jd_end`` are optional TT JDs;
/// omitted bounds use a wall-clock-derived now +/-3652.5-day estimate (not a
/// precision TT clock). Bounds are clipped to the native Horizons generation
/// envelope and rendered as UTC dates. ``boundary_padding_days=0.0`` expands
/// that clipping envelope, not the requested interval; intended domain is
/// finite/nonnegative. ``force=False`` returns an existing registration without
/// testing requested coverage. Otherwise local cache/download and registration
/// follow configuration (30-second transport budgets).
///
/// Unsupported, unsatisfied, malformed-window and most download/registration
/// failures return ``None``. Sealed-network and configuration errors propagate;
/// a successful registration is shared with future calculations. This helper
/// can write/download files and is not a precision/coverage certification.
#[pyfunction(signature = (body_id, jd_start = None, jd_end = None, force = false, boundary_padding_days = 0.0))]
fn auto_download_asteroid_spk(
    py: Python<'_>,
    body_id: i32,
    jd_start: Option<f64>,
    jd_end: Option<f64>,
    force: bool,
    boundary_padding_days: f64,
) -> PyResult<Option<String>> {
    py.detach(|| {
        runtime().auto_download_asteroid_spk(
            body_id,
            jd_start,
            jd_end,
            force,
            boundary_padding_days,
        )
    })
    .map(|path| path.map(|p| p.to_string_lossy().into_owned()))
    .map_err(|error| to_pyerr(py, error))
}

/// Test explicit runtime SPK registration for library ``body_id``.
/// No coefficients/files are read; true means membership only, not epoch coverage.
#[pyfunction]
fn is_spk_available_for_body(py: Python<'_>, body_id: i32) -> bool {
    py.detach(|| runtime().is_spk_available_for_body(body_id))
}

/// Ensure an asteroid SPK registration, returning a boolean.
///
/// ``body_id`` is a library ID; optional ``jd`` is TT JD for a +/-3652.5-day
/// provisioning window. ``jd=None`` uses the default wall-clock estimate.
/// Existing registration returns true without a coverage probe. Unsupported or
/// failed provisioning returns false; configuration/sealed-network errors raise.
/// May download/write/register a kernel via ``auto_download_asteroid_spk``.
#[pyfunction(signature = (body_id, jd = None))]
fn ensure_major_asteroid_spk(py: Python<'_>, body_id: i32, jd: Option<f64>) -> PyResult<bool> {
    py.detach(|| runtime().ensure_major_asteroid_spk(body_id, jd))
        .map_err(|error| to_pyerr(py, error))
}

/// Test whether library ``body_id`` is in the major-asteroid SPK metadata table.
/// This is not a mass threshold, a live JPL query or an availability/coverage check.
#[pyfunction]
fn is_major_asteroid(body_id: i32) -> bool {
    libephemeris_core::spk::is_major_asteroid(body_id)
}

/// Return ``(asteroid_number, horizons_id, naif_id, display_name)`` or ``None``.
/// ``body_id`` is a library ID. The Horizons identifier is a string, possibly
/// using ``Name;`` syntax. No network, kernel acquisition or coverage test occurs.
#[pyfunction]
fn get_major_asteroid_info(body_id: i32) -> Option<(i32, &'static str, i32, &'static str)> {
    libephemeris_core::spk::get_major_asteroid_info(body_id)
}

/// Return a list of ``(library_body_id, display_name)`` major-asteroid rows.
/// Registry order is retained; the list describes SPK support, not loaded data.
#[pyfunction]
fn list_major_asteroids() -> Vec<(i32, &'static str)> {
    bodies::list_major_asteroids()
}

/// Detect proximity to a listed Neptune mean-motion ratio; return ``ResonanceResult``.
///
/// ``elements`` must be ``libephemeris.minor_bodies.OrbitalElements``, not the
/// hypothetical polynomial record. ``tolerance=0.02`` is a fractional ratio
/// deviation (2%, not degrees). Axes below 20 AU never match; the closest ratio
/// strictly below tolerance wins. A match carries advisory text, not an emitted
/// warning or proof of libration. Finite element conversion errors raise;
/// tolerance has no additional range validation. Example:
///
/// ```python
/// from libephemeris.minor_bodies import (
///     MINOR_BODY_ELEMENTS, detect_mean_motion_resonance)
/// result = detect_mean_motion_resonance(next(iter(MINOR_BODY_ELEMENTS.values())))
/// ```
#[pyfunction(signature = (elements, tolerance = 0.02))]
fn detect_mean_motion_resonance(
    py: Python<'_>,
    elements: OrbitalElements,
    tolerance: f64,
) -> PyResult<ResonanceResult> {
    let elements = elements.native().map_err(|error| to_pyerr(py, error))?;
    Ok(py
        .detach(|| bodies::detect_mean_motion_resonance(&elements, tolerance))
        .into())
}

/// Immutable approximate resonant-TNO libration parameters.
///
/// These project models are not a replacement for integrated resonant dynamics;
/// Gonggong's row is an estimate, not a published fit. Construction is a value
/// operation and does not install a model for a new body.
#[derive(Clone, Debug)]
#[pyclass(
    module = "libephemeris.minor_bodies",
    frozen,
    get_all,
    skip_from_py_object
)]
struct LibrationParameters {
    /// Maximum resonant-argument displacement from center, degrees.
    amplitude: f64,
    /// Libration period in days.
    period: f64,
    /// Resonant-argument center in degrees.
    center: f64,
    /// Sinusoid phase at J2000.0 in RADIANS, unlike amplitude/center.
    phase_j2000: f64,
    /// Numerator p of the resonance mean-motion ratio p/q.
    resonance_p: i32,
    /// Denominator q of the resonance ratio; longitude correction divides by q.
    resonance_q: i32,
}
#[pymethods]
impl LibrationParameters {
    /// Store amplitude/center (deg), period (days), J2000 phase (rad), and p/q.
    /// No finite, positive-period or resonance-ratio validation is performed.
    #[new]
    fn new(
        amplitude: f64,
        period: f64,
        center: f64,
        phase_j2000: f64,
        resonance_p: i32,
        resonance_q: i32,
    ) -> Self {
        Self {
            amplitude,
            period,
            center,
            phase_j2000,
            resonance_p,
            resonance_q,
        }
    }
}
impl From<bodies::LibrationParameters> for LibrationParameters {
    fn from(r: bodies::LibrationParameters) -> Self {
        Self::new(
            r.amplitude,
            r.period,
            r.center,
            r.phase_j2000,
            r.resonance_p,
            r.resonance_q,
        )
    }
}
/// Return whether curated library ``body_id`` is near a listed Neptune ratio.
/// ``tolerance=0.02`` is a fractional mean-motion-ratio threshold, not an angle.
/// Unknown IDs raise ``ValueError``; a true result does not establish libration.
#[pyfunction(signature=(body_id,tolerance=0.02))]
fn is_body_resonant(py: Python<'_>, body_id: i32, tolerance: f64) -> PyResult<bool> {
    py.detach(|| bodies::is_body_resonant(body_id, tolerance))
        .map_err(|e| PyValueError::new_err(e.to_string()))
}
/// Return a default-2%-threshold ``ResonanceResult`` for a curated library ID.
/// Unknown IDs return ``None``; known nonresonant bodies return an unmatched
/// record (not None), with no resonance/advisory and zero deviation.
#[pyfunction]
fn get_resonance_info(py: Python<'_>, body_id: i32) -> Option<ResonanceResult> {
    py.detach(|| bodies::get_resonance_info(body_id))
        .map(ResonanceResult::from)
}
/// Return linear J2000-reference Neptune mean longitude in degrees [0,360).
/// ``jd_tt`` is TT JD. This approximate resonance reference has no apparent
/// corrections or time validation; non-finite inputs may return NaN.
#[pyfunction]
fn calc_neptune_mean_longitude(py: Python<'_>, jd_tt: f64) -> f64 {
    py.detach(|| bodies::calc_neptune_mean_longitude(jd_tt))
}
/// Test whether library ``body_id`` has project libration parameters.
/// True currently for Ixion, Orcus and Gonggong; unknown IDs return false.
#[pyfunction]
fn has_libration_model(body_id: i32) -> bool {
    bodies::has_libration_model(body_id)
}
/// Return approximate ``LibrationParameters`` for library ``body_id``, or None.
/// Angles are degrees except ``phase_j2000`` in radians; period is days. These
/// values describe the project sinusoid, not a precision integrated solution.
#[pyfunction]
fn get_libration_parameters(body_id: i32) -> Option<LibrationParameters> {
    bodies::get_libration_parameters(body_id).map(LibrationParameters::from)
}
/// Return approximate longitude correction in degrees at ``jd_tt`` (TT JD).
/// For a modeled body it is ``(amplitude/q)*sin(2*pi*dt/period + phase_j2000)``.
/// Unknown/unmodeled IDs return 0, including non-finite dates; modeled dates
/// have no finite check and may yield NaN. This is not a resonant event finder.
#[pyfunction]
fn calc_libration_correction(py: Python<'_>, body_id: i32, jd_tt: f64) -> f64 {
    py.detach(|| bodies::calc_libration_correction(body_id, jd_tt))
}
/// Return ``(domega, dOmega, dn)`` in degrees/day for minor-body elements.
///
/// ``elements`` is ``libephemeris.minor_bodies.OrbitalElements`` (axis AU,
/// dimensionless e, degree angles, degree/day n). The reduced giant-planet
/// secular approximation includes a near-resonance heuristic, not full N-body
/// integration. Finite field conversion is checked; physical domains are not.
/// Intended for bound nonencounter orbits: invalid/hyperbolic geometry may yield
/// NaN/infinity rather than a convergence error. No ephemeris acquisition occurs.
#[pyfunction]
fn calc_secular_perturbation_rates(
    py: Python<'_>,
    elements: OrbitalElements,
) -> PyResult<(f64, f64, f64)> {
    let row = elements.native().map_err(|e| to_pyerr(py, e))?;
    Ok(py.detach(|| bodies::secular_perturbation_rates(&row)))
}
/// Return propagated ``(omega, Omega, M, n, e, i)`` for minor-body elements.
///
/// Use ``libephemeris.minor_bodies.OrbitalElements`` and TT ``jd_tt``. Angles
/// are degrees, n degrees/day and e dimensionless; omega is the argument of
/// perihelion, Omega the node longitude. ``include_perturbations=True`` uses
/// approximate Laplace-Lagrange secular evolution; false or |dt|<1 day retains
/// orientation/e/i and advances M uniformly. Perturbed e is clamped to [0,0.999].
/// Not a resonance/encounter model. Finite element conversion errors raise,
/// but time/physical-domain validation is not added; non-finite results may escape.
#[pyfunction(signature=(elements,jd_tt,include_perturbations=true))]
fn apply_secular_perturbations(
    py: Python<'_>,
    elements: OrbitalElements,
    jd_tt: f64,
    include_perturbations: bool,
) -> PyResult<(f64, f64, f64, f64, f64, f64)> {
    let row = elements.native().map_err(|e| to_pyerr(py, e))?;
    Ok(py.detach(|| bodies::apply_secular_perturbations(&row, jd_tt, include_perturbations)))
}
/// Return plutino resonant argument ``3*lambda-2*lambda_Neptune-varpi`` in degrees.
///
/// ``elements`` is ``minor_bodies.OrbitalElements``; ``jd_tt`` TT JD;
/// ``omega_pert`` is the perturbed perihelion ARGUMENT and ``M_pert`` mean
/// anomaly, both degrees. Uses the stored node longitude to form lambda/varpi;
/// result is wrapped to [0,360). Finite element conversion is checked, but
/// scalar date/angle inputs are not; no test of actual libration is performed.
#[pyfunction]
#[allow(non_snake_case)]
fn calc_resonant_argument_plutino(
    py: Python<'_>,
    elements: OrbitalElements,
    jd_tt: f64,
    omega_pert: f64,
    M_pert: f64,
) -> PyResult<f64> {
    let row = elements.native().map_err(|e| to_pyerr(py, e))?;
    Ok(py.detach(|| bodies::calc_resonant_argument_plutino(&row, jd_tt, omega_pert, M_pert)))
}
macro_rules! kepler {
    ($(($name:ident, $doc:literal)),+)=>{$(
        #[doc = $doc]
        ///
        /// ``M`` is the mean-anomaly parameter in radians, ``e`` dimensionless
        /// eccentricity (NOT radians), and ``tol=1e-8`` an absolute successive-
        /// iterate tolerance in radians. Intended inputs are finite, physical e
        /// and positive tol; no validation or convergence exception is provided.
        /// Exhaustion returns the last iterate; invalid arithmetic may return
        /// NaN/infinity. A finite result is not a convergence certificate.
        #[pyfunction(signature=(M,e,tol=1e-8))]
        #[allow(non_snake_case)]
        fn $name(py:Python<'_>,M:f64,e:f64,tol:f64)->f64 { py.detach(|| bodies::$name(M,e,tol)) }
    )+};
}
kepler!(
    (
        solve_kepler_equation_elliptic,
        "Solve M=E-e*sin(E) for eccentric anomaly E in radians (0<=e<1).\n\nThe signed remainder M modulo 2*pi is reduced by 2*pi when above pi; negative remainders are not symmetrically folded. Uses at most 30 Newton steps with the model's high-e starter. Branch wrapping is not a revolution count. Example: E = solve_kepler_equation_elliptic(1.0, 0.1)."
    ),
    (
        solve_kepler_equation_hyperbolic,
        "Solve M=e*sinh(H)-H for hyperbolic anomaly H (radian convention), e>1.\n\nM is not wrapped. Uses at most 50 Newton steps and stops early if the derivative magnitude is below 1e-15, returning the current iterate without raising."
    ),
    (
        solve_kepler_equation,
        "Dispatch by dimensionless eccentricity; the returned anomaly depends on orbit class.\n\nFor abs(e-1)<1e-10 return the Barker-model TRUE anomaly in radians (tol unused); below unity return elliptic eccentric anomaly E, above return hyperbolic anomaly H. See the class-specific solvers for bounded-iteration/nonconvergence behavior."
    )
);
/// Return the parabolic model's true anomaly in radians from scalar parameter M.
///
/// Evaluates ``2*atan(2*sinh(asinh(3*M)/3))`` directly: if D=tan(nu/2),
/// the convention is ``D + D**3/3 = 2*M``. M is dimensionless with the
/// propagation API's radian mean-anomaly convention, not degrees or a time.
/// There is no iteration/tolerance or validation; overflow/non-finite inputs
/// follow IEEE arithmetic. This returns true anomaly, not elliptic E/hyperbolic H.
#[pyfunction]
#[allow(non_snake_case)]
fn solve_barker_equation(py: Python<'_>, M: f64) -> f64 {
    py.detach(|| bodies::solve_barker_equation(M))
}
/// Return heliocentric J2000-ecliptic Cartesian ``(x, y, z)`` in AU.
///
/// ``elements`` must be ``libephemeris.minor_bodies.OrbitalElements``;
/// ``jd_tt`` is TT JD. ``include_perturbations=True`` enables approximate secular
/// and short-period corrections for bound orbits. Optional library ``body_id``
/// selects modeled libration/Yarkovsky terms; None supplies no body-specific
/// correction. These elliptic body-specific terms are applied independently of
/// the perturbation switch. Elliptic, parabolic and hyperbolic branches are supported;
/// for ``abs(e-1)<1e-10`` the parabolic branch uses ``elements.a`` as q in AU.
/// Finite element/time conversion is checked, but not full physical validity;
/// invalid or nonconvergent geometry can produce non-finite coordinates without
/// raising. This is approximate geometry, not an apparent/SPK/N-body position.
/// Example (the two element classes are not interchangeable):
///
/// ```python
/// from libephemeris.minor_bodies import OrbitalElements as MinorElements
/// from libephemeris.minor_bodies import calc_minor_body_position
/// orbit = MinorElements("example", 2451545.0, 2.5, 0.1, 5.0,
///                       40.0, 80.0, 10.0, 0.25)
/// xyz = calc_minor_body_position(orbit, 2451545.0, include_perturbations=False)
/// ```
#[pyfunction(signature=(elements,jd_tt,include_perturbations=true,body_id=None))]
fn calc_minor_body_position(
    py: Python<'_>,
    elements: OrbitalElements,
    jd_tt: f64,
    include_perturbations: bool,
    body_id: Option<i32>,
) -> PyResult<(f64, f64, f64)> {
    let row = elements.native().map_err(|e| to_pyerr(py, e))?;
    py.detach(|| bodies::calc_minor_body_position(&row, jd_tt, include_perturbations, body_id))
        .map(|[x, y, z]| (x, y, z))
        .map_err(|e| to_pyerr(py, e))
}
/// Return heliocentric ``(longitude_deg, latitude_deg, range_au)`` at TT ``jd_tt``.
///
/// ``body_id`` is a curated library ID. ``use_spk=True`` first tries a registered
/// SPK through configured native calculation; absence/range misses fall through
/// to epoch-blended Keplerian/secular geometry, precessed/nutated to true ecliptic
/// of date. False skips the SPK probe. This helper is not the full automatic
/// SPK/ASSIST precision chain. Unknown IDs raise ``ValueError``; non-finite time,
/// configuration, frame and non-range provider failures propagate. Approximate
/// accuracy degrades away from element epochs and near resonances/encounters.
#[pyfunction(signature=(body_id,jd_tt,use_spk=true))]
fn calc_minor_body_heliocentric(
    py: Python<'_>,
    body_id: i32,
    jd_tt: f64,
    use_spk: bool,
) -> PyResult<(f64, f64, f64)> {
    py.detach(|| runtime().calc_minor_body_heliocentric(body_id, jd_tt, use_spk))
        .map_err(|error| match error {
            Error::UnknownBody { .. } => PyValueError::new_err(error.to_string()),
            e => to_pyerr(py, e),
        })
}
/// Test native SPK download eligibility for library ``body_id``.
/// Major-asteroid/mapped targets qualify except explicitly blocked bodies.
/// No live JPL, network-policy, file or coverage test is performed.
#[pyfunction]
fn is_spk_downloadable(body_id: i32) -> bool {
    libephemeris_core::spk::is_spk_downloadable(body_id)
}
/// Compatibility predicate that always returns false, regardless of ``body_id``.
/// Native Horizons requests use ``Name;`` syntax to bypass the historical index
/// restriction; this is not a live query of JPL catalogue membership.
#[pyfunction]
fn is_in_jpl_major_body_index(body_id: i32) -> bool {
    libephemeris_core::spk::is_in_jpl_major_body_index(body_id)
}
/// Log physical-sanity advisories for ``minor_bodies.OrbitalElements``; return None.
///
/// Checks axis sign, eccentricity/orbit sign consistency, inclination domain
/// `0..=180` degrees and positive elliptic mean motion. Advisories are log messages,
/// not returned strings or raised warnings, and elements are unchanged. Non-finite
/// fields fail native conversion before logging. Use ``elements.validate()``
/// when a list of strings is needed; this is not full orbital validation.
#[pyfunction]
fn validate_elements(py: Python<'_>, elements: OrbitalElements) -> PyResult<()> {
    let row = elements.native().map_err(|e| to_pyerr(py, e))?;
    py.detach(|| bodies::validate_elements_and_log(&row));
    Ok(())
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add("_ResonanceInfo", module.py().get_type::<ResonanceInfo>())?;
    module.add(
        "_ResonanceResult",
        module.py().get_type::<ResonanceResult>(),
    )?;
    macro_rules! add { ($($name:ident),+ $(,)?)=>{$(module.add_function(wrap_pyfunction!($name,module)?)?;)+}; }
    add!(
        get_asteroid_number,
        clear_asteroid_name_cache,
        calc_asteroid_by_number,
        fetch_orbital_elements_from_sbdb,
        clear_asteroid_elements_cache,
        auto_download_asteroid_spk,
        is_spk_available_for_body,
        ensure_major_asteroid_spk,
        is_major_asteroid,
        get_major_asteroid_info,
        list_major_asteroids,
        detect_mean_motion_resonance,
        validate_elements,
        is_body_resonant,
        get_resonance_info,
        calc_neptune_mean_longitude,
        calc_resonant_argument_plutino,
        calc_libration_correction,
        has_libration_model,
        get_libration_parameters,
        calc_secular_perturbation_rates,
        apply_secular_perturbations,
        solve_kepler_equation_elliptic,
        solve_kepler_equation_hyperbolic,
        solve_kepler_equation,
        solve_barker_equation,
        calc_minor_body_position,
        calc_minor_body_heliocentric,
        is_spk_downloadable,
        is_in_jpl_major_body_index
    );
    module.add_class::<LibrationParameters>()?;
    let elements = PyDict::new(module.py());
    let libration = PyDict::new(module.py());
    let yarkovsky = PyDict::new(module.py());
    for (id, row) in bodies::ordered_minor_body_elements() {
        elements.set_item(id, OrbitalElements::from_native(row))?;
        if let Some(row) = bodies::get_libration_parameters(id) {
            libration.set_item(id, LibrationParameters::from(row))?;
        }
        if let Some(rate) = bodies::yarkovsky_da_dt_au_per_my(id) {
            yarkovsky.set_item(id, rate)?;
        }
    }
    let multi = PyDict::new(module.py());
    let mut grouped = std::collections::HashMap::<i32, Vec<OrbitalElements>>::new();
    for (body, row) in bodies::ordered_minor_body_multi_elements() {
        grouped
            .entry(body)
            .or_default()
            .push(OrbitalElements::from_native(row));
    }
    for (body, rows) in grouped {
        multi.set_item(body, PyList::new(module.py(), rows)?)?;
    }
    module.add("MINOR_BODY_ELEMENTS_MULTI", multi)?;
    let resonances = bodies::NEPTUNE_RESONANCES
        .iter()
        .map(|&(p, q, name, tolerance)| ResonanceInfo {
            p,
            q,
            name: name.into(),
            a_resonant: bodies::resonant_semi_major_axis(p, q),
            tolerance,
        });
    module.add("NEPTUNE_RESONANCES", PyList::new(module.py(), resonances)?)?;
    module.add("MINOR_BODY_ELEMENTS", elements)?;
    module.add("PLUTINO_LIBRATION_PARAMS", libration)?;
    module.add("YARKOVSKY_DA_DT", yarkovsky)?;
    module.add("_Minor_GM_SUN", bodies::GM_SUN)?;
    module.add("_Minor_MASS_RATIO_JUPITER", bodies::MASS_RATIO_JUPITER)?;
    module.add("_Minor_MASS_RATIO_SATURN", bodies::MASS_RATIO_SATURN)?;
    module.add("_Minor_MASS_RATIO_URANUS", bodies::MASS_RATIO_URANUS)?;
    module.add("_Minor_MASS_RATIO_NEPTUNE", bodies::MASS_RATIO_NEPTUNE)?;
    module.add("_Minor_JUPITER_A", bodies::JUPITER_A)?;
    module.add("_Minor_JUPITER_E", bodies::JUPITER_E)?;
    module.add("_Minor_JUPITER_I", bodies::JUPITER_I)?;
    module.add("_Minor_JUPITER_OMEGA", bodies::JUPITER_OMEGA)?;
    module.add("_Minor_JUPITER_NODE", bodies::JUPITER_NODE)?;
    module.add("_Minor_JUPITER_N", bodies::JUPITER_N)?;
    module.add("_Minor_SATURN_A", bodies::SATURN_A)?;
    module.add("_Minor_SATURN_E", bodies::SATURN_E)?;
    module.add("_Minor_SATURN_I", bodies::SATURN_I)?;
    module.add("_Minor_SATURN_OMEGA", bodies::SATURN_OMEGA)?;
    module.add("_Minor_SATURN_NODE", bodies::SATURN_NODE)?;
    module.add("_Minor_SATURN_N", bodies::SATURN_N)?;
    module.add("_Minor_URANUS_A", bodies::URANUS_A)?;
    module.add("_Minor_URANUS_E", bodies::URANUS_E)?;
    module.add("_Minor_URANUS_I", bodies::URANUS_I)?;
    module.add("_Minor_URANUS_OMEGA", bodies::URANUS_OMEGA)?;
    module.add("_Minor_URANUS_NODE", bodies::URANUS_NODE)?;
    module.add("_Minor_URANUS_N", bodies::URANUS_N)?;
    module.add("_Minor_NEPTUNE_A", bodies::NEPTUNE_A)?;
    module.add("_Minor_NEPTUNE_E", bodies::NEPTUNE_E)?;
    module.add("_Minor_NEPTUNE_I", bodies::NEPTUNE_I)?;
    module.add("_Minor_NEPTUNE_OMEGA", bodies::NEPTUNE_OMEGA)?;
    module.add("_Minor_NEPTUNE_NODE", bodies::NEPTUNE_NODE)?;
    module.add("_Minor_NEPTUNE_N", bodies::NEPTUNE_N)?;
    module.add("_Minor_JUPITER_DE", bodies::JUPITER_DE)?;
    module.add("_Minor_JUPITER_DI", bodies::JUPITER_DI)?;
    module.add("_Minor_JUPITER_DOMEGA", bodies::JUPITER_DOMEGA)?;
    module.add("_Minor_JUPITER_DNODE", bodies::JUPITER_DNODE)?;
    module.add("_Minor_SATURN_DE", bodies::SATURN_DE)?;
    module.add("_Minor_SATURN_DI", bodies::SATURN_DI)?;
    module.add("_Minor_SATURN_DOMEGA", bodies::SATURN_DOMEGA)?;
    module.add("_Minor_SATURN_DNODE", bodies::SATURN_DNODE)?;
    module.add("_Minor_URANUS_DE", bodies::URANUS_DE)?;
    module.add("_Minor_URANUS_DI", bodies::URANUS_DI)?;
    module.add("_Minor_URANUS_DOMEGA", bodies::URANUS_DOMEGA)?;
    module.add("_Minor_URANUS_DNODE", bodies::URANUS_DNODE)?;
    module.add("_Minor_NEPTUNE_DE", bodies::NEPTUNE_DE)?;
    module.add("_Minor_NEPTUNE_DI", bodies::NEPTUNE_DI)?;
    module.add("_Minor_NEPTUNE_DOMEGA", bodies::NEPTUNE_DOMEGA)?;
    module.add("_Minor_NEPTUNE_DNODE", bodies::NEPTUNE_DNODE)?;
    module.add("_Minor_J2000_EPOCH", bodies::J2000_EPOCH)?;
    module.add(
        "_Minor_NEPTUNE_MEAN_LONGITUDE_J2000",
        bodies::NEPTUNE_MEAN_LONGITUDE_J2000,
    )?;
    module.add("_Minor_NEPTUNE_MEAN_MOTION", bodies::NEPTUNE_MEAN_MOTION)?;
    module.add(
        "_Minor_DEFAULT_RESONANCE_TOLERANCE",
        bodies::DEFAULT_RESONANCE_TOLERANCE,
    )?;
    module.add(
        "_Minor_EPOCH_BLEND_HALF_DAYS",
        bodies::EPOCH_BLEND_HALF_DAYS,
    )?;
    module.add("_Minor_EPOCH_MERGE_TOL_DAYS", bodies::EPOCH_MERGE_TOL_DAYS)?;
    module.add("_Minor_SBDB_API_URL", bodies::SBDB_API_URL)?;
    module.add(
        "_Minor_HORIZONS_SPK_DATE_MIN",
        bodies::HORIZONS_SPK_DATE_MIN,
    )?;
    module.add(
        "_Minor_HORIZONS_SPK_DATE_MAX",
        bodies::HORIZONS_SPK_DATE_MAX,
    )?;
    module.add("_Minor_HORIZONS_SPK_JD_MIN", bodies::HORIZONS_SPK_JD_MIN)?;
    module.add("_Minor_HORIZONS_SPK_JD_MAX", bodies::HORIZONS_SPK_JD_MAX)?;
    Ok(())
}
