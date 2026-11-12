//! Satellite SPK configuration and numeric TT positions from the native engine.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::errors::Error;
use libephemeris_core::planetary_moons as moons;
use pyo3::exceptions::{PyFileNotFoundError, PyValueError};
use pyo3::prelude::*;
use std::collections::HashMap;

type MoonPosition = (f64, f64, f64, f64, f64, f64);

/// Register a local satellite SPK for the process-wide moon registry.
///
/// ``spk_file`` is a filename; relative paths search the configured data directory
/// before the current directory. ``moons=None`` discovers satellite targets;
/// otherwise supply library legacy IDs or canonical ``PLMOON_OFFSET + NAIF``
/// IDs. Unknown/unchainable selections are skipped; if none register,
/// ``ValueError`` is raised. Missing files raise ``FileNotFoundError``; SPK
/// decoding/configuration errors propagate. No download is performed. Returns
/// ``None``; registration alone does not guarantee coverage at a given epoch.
#[pyfunction(signature = (spk_file, moons = None))]
fn register_moon_spk(py: Python<'_>, spk_file: &str, moons: Option<Vec<i32>>) -> PyResult<()> {
    py.detach(|| runtime().register_moon_spk(spk_file, moons.as_deref()))
        .map_err(|error| match error {
            Error::SpkFileNotFound { .. } => PyFileNotFoundError::new_err(error.to_string()),
            Error::SpkNoValidMoons { .. } => PyValueError::new_err(error.to_string()),
            error => to_pyerr(py, error),
        })
}

/// Remove registrations matching ``spk_file`` or its nonempty basename.
///
/// Returns ``None`` and does nothing if no kernel matches. Already running
/// calculations retain their immutable images; unrelated providers are unchanged.
#[pyfunction]
fn unregister_moon_spk(py: Python<'_>, spk_file: &str) {
    py.detach(|| runtime().unregister_moon_spk(spk_file));
}

/// Return a new ``{moon_id: spk_path}`` dict for registered satellite aliases.
///
/// Includes canonical/legacy moon IDs, not parent-chain support targets.
/// An empty registry returns ``{}``; this metadata query does not open files.
#[pyfunction]
fn list_registered_moons(py: Python<'_>) -> HashMap<i32, String> {
    py.detach(|| runtime().list_registered_moons())
}

/// Return the display name for a known legacy or canonical moon ID.
///
/// Unknown IDs return ``"Unknown Moon (<id>)"`` rather than raising.
/// Recognition does not imply that a kernel has been registered.
#[pyfunction]
fn get_moon_name(moon_id: i32) -> String {
    moons::get_moon_name(moon_id)
}

/// Test the planetary-moon identifier grammar, not kernel availability.
///
/// Accepts legacy ``MOON_*`` aliases and canonical ``PLMOON_OFFSET + NAIF``
/// IDs with NAIF in 401..999. A true result need not identify a named or
/// registered target. Earth's Moon uses the separate lunar API.
#[pyfunction]
fn is_planetary_moon(ipl: i32) -> bool {
    moons::is_planetary_moon(ipl)
}

/// Calculate a registered moon's configured apparent position at ``t`` (TT JD).
///
/// ``moon_id`` is a legacy/canonical library ID; ``iflag`` is a calculation
/// bitmask. Return ``(lon, lat, range, dlon, dlat, drange)`` in degrees, degrees,
/// AU and corresponding units per TT day, on the ecliptic of date. This direct
/// helper does not apply EQUATORIAL/J2000/XYZ/RADIANS output shaping. Observer
/// precedence is BARYCTR, HELCTR, TOPOCTR, geocenter; TOPOCTR needs ``set_topo``.
/// Flags control light time, aberration, nutation, sidereal reduction and SPEED;
/// absent SPEED, rates are zero. Angular rates use adaptive TT differencing and
/// radial rate differentiates the light-time chain. It is not a planet-relative
/// analytical vector; NOGDEFL has no extra deflection stage to disable.
///
/// Unknown/unregistered targets return ``None`` after finite-time admission,
/// before main-provider acquisition. Otherwise missing main/site data, explicit
/// LEB-only mode, SPK coverage and frame errors raise; no analytical satellite
/// substitutes for a failed kernel. Main data acquisition follows configuration.
#[pyfunction]
fn calc_moon_position(
    py: Python<'_>,
    t: f64,
    moon_id: i32,
    iflag: i32,
) -> PyResult<Option<MoonPosition>> {
    py.detach(|| runtime().calc_moon_position(t, moon_id, iflag))
        .map(|position| position.map(|[a, b, c, d, e, f]| (a, b, c, d, e, f)))
        .map_err(|error| to_pyerr(py, error))
}

/// Return the registered kernel's ``(start_jd, end_jd)`` envelope in TDB.
///
/// ``None`` means missing registration/image or unavailable finite segment bounds.
/// Bounds span the whole kernel, not the intersection needed by this moon's center
/// chain; gaps and target-specific limits can still make calculation fail.
/// These are not TT input bounds. No coefficients are sampled or files opened.
#[pyfunction]
fn get_moon_coverage(py: Python<'_>, moon_id: i32) -> Option<(f64, f64)> {
    py.detach(|| runtime().get_moon_coverage(moon_id))
}

/// Clear all runtime-owned satellite kernels and registrations; return ``None``.
///
/// Other providers are unchanged, and in-flight calculations retain their own
/// shared images until completion. Repeated calls are harmless.
#[pyfunction]
fn close_moon_kernels(py: Python<'_>) {
    py.detach(|| runtime().close_moon_kernels());
}

macro_rules! moon_vector {
    ($name:ident, $target:ident, $doc:literal) => {
        #[doc = $doc]
        #[pyfunction]
        fn $name(py: Python<'_>, jd: f64) -> PyResult<(f64, f64, f64)> {
            py.detach(|| moons::$target(jd))
                .map_err(|error| to_pyerr(py, error))
        }
    };
}
moon_vector!(
    triton_position,
    triton_position,
    "Return Triton's Neptune-center-relative (x, y, z) in AU, J2000 ecliptic axes.\n\njd is a finite TT Julian day. This approximate NEP 097 element model includes nodal precession, not light time or observer corrections; nominal accuracy is roughly 20-50 km, not guaranteed at arbitrary dates. No kernel is required. Non-finite jd raises."
);
moon_vector!(
    charon_position,
    charon_position,
    "Return Charon's Pluto-center-relative (x, y, z) in AU, J2000 ecliptic axes.\n\njd is a finite TT Julian day. The frozen PLU 060 two-body model has no light-time/observer corrections and nominal 1-5 km accuracy, with limited long-term validity. No kernel is required. Non-finite jd raises."
);
moon_vector!(
    uranus_cob_offset,
    uranus_cob_offset,
    "Return the Uranus-system-barycenter to Uranus-center offset in AU.\n\njd is a finite TT Julian day; (x, y, z) is equatorial ICRF/J2000, unlike the ecliptic Triton/Charon vectors. Add it to a system-barycenter vector in the same frame. Uses a negative GM-weighted sum of five approximate circular moon models, not a registered SPK. Non-finite jd raises."
);

/// Return Jupiter-relative ``(io, europa, ganymede, callisto)`` vectors in km.
///
/// Each vector is ``(x, y, z)`` in the model's J2000-ecliptic axes. ``jd`` is
/// used directly as dynamical Julian date (TT/TDB, with no conversion). The
/// approximate E5 series applies no light time or observer correction; nominal
/// accuracy is 100-200 km, not precision SPK output. Non-finite time raises.
/// Example:
///
/// ```python
/// from libephemeris.planetary_moons import galilean_moon_positions
/// io, europa, ganymede, callisto = galilean_moon_positions(2451545.0)
/// ```
#[pyfunction]
fn galilean_moon_positions(py: Python<'_>, jd: f64) -> PyResult<moons::GalileanMoonQuartet> {
    py.detach(|| moons::galilean_moon_positions(jd))
        .map_err(|error| to_pyerr(py, error))
}

/// Return Uranus-relative ``(ariel, umbriel, titania, oberon, miranda)`` in km.
///
/// Each vector is ``(x, y, z)`` in equatorial ICRF/J2000 axes. ``jd`` is a TT
/// Julian day used directly. Circular mean-element geometry omits GUST86 periodic
/// terms and observer/light-time effects; moon errors can be thousands of km.
/// Intended for approximate center-of-body offsets, not precision moon places.
/// Non-finite time raises; no kernel is required.
#[pyfunction]
fn uranian_moon_positions(py: Python<'_>, jd: f64) -> PyResult<moons::UranianMoonQuintet> {
    py.detach(|| moons::uranian_moon_positions(jd))
        .map_err(|error| to_pyerr(py, error))
}

/// Return a Saturn-center-relative TASS1.7 ``(x, y, z)`` vector in AU.
///
/// ``jd`` is a finite TT Julian day; axes are the emitted J2000 ecliptic model
/// frame, with no observer/light-time reduction. ``body`` is a TASS index:
/// 0 Mimas, 1 Enceladus, 2 Tethys, 3 Dione, 4 Rhea, 5 Titan (default),
/// 6 Iapetus, 7 Hyperion. Outside 0..7 raises ``ValueError``; non-finite time
/// raises. Nominal model accuracy is approximately 50-100 km per satellite.
#[pyfunction(signature = (jd, body = 5))]
fn saturn_moon_position(py: Python<'_>, jd: f64, body: i32) -> PyResult<(f64, f64, f64)> {
    py.detach(|| moons::saturn_moon_position(jd, body))
        .map_err(|error| tass_error(py, error))
}

/// Return Saturn-relative TASS1.7 ``(x, y, z, vx, vy, vz)`` in AU and AU/day.
///
/// ``jd`` is TT Julian day; position and velocity use the emitted J2000 ecliptic
/// model axes, without apparent-place corrections. ``body`` is 0 Mimas,
/// 1 Enceladus, 2 Tethys, 3 Dione, 4 Rhea, 5 Titan (default), 6 Iapetus or
/// 7 Hyperion. Invalid indices raise ``ValueError``; non-finite time raises.
/// The approximate model's limits are those of ``saturn_moon_position``.
#[pyfunction(signature = (jd, body = 5))]
fn saturn_moon_position_velocity(py: Python<'_>, jd: f64, body: i32) -> PyResult<MoonPosition> {
    py.detach(|| moons::saturn_moon_position_velocity(jd, body))
        .map_err(|error| tass_error(py, error))
}

/// Translate the analytical TASS index refusal without changing other bodies.
fn tass_error(py: Python<'_>, error: Error) -> PyErr {
    match error.primary_error() {
        Error::InvalidBody { .. } => PyValueError::new_err(error.to_string()),
        _ => to_pyerr(py, error),
    }
}

/// Return system-barycenter to planet-center displacement ``(x, y, z)`` in AU.
///
/// ``t`` is TT Julian day; output axes are equatorial ICRF/J2000. Add the
/// vector to a barycenter position in the same frame. ``planet_name`` accepts
/// Jupiter, Saturn, Uranus, Neptune or Pluto, ignoring ASCII case, surrounding
/// whitespace and ``" barycenter"``. Other names return the zero vector, even
/// for non-finite time; known names reject non-finite time. Uses approximate
/// analytical moon models, not configured SPK apparent positions.
#[pyfunction]
fn get_cob_offset(py: Python<'_>, planet_name: &str, t: f64) -> PyResult<(f64, f64, f64)> {
    py.detach(|| moons::get_cob_offset(planet_name, t))
        .map_err(|error| to_pyerr(py, error))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(triton_position, module)?)?;
    module.add_function(wrap_pyfunction!(charon_position, module)?)?;
    module.add_function(wrap_pyfunction!(uranus_cob_offset, module)?)?;
    module.add_function(wrap_pyfunction!(galilean_moon_positions, module)?)?;
    module.add_function(wrap_pyfunction!(uranian_moon_positions, module)?)?;
    module.add_function(wrap_pyfunction!(saturn_moon_position, module)?)?;
    module.add_function(wrap_pyfunction!(saturn_moon_position_velocity, module)?)?;
    module.add_function(wrap_pyfunction!(get_cob_offset, module)?)?;
    module.add_function(wrap_pyfunction!(register_moon_spk, module)?)?;
    module.add_function(wrap_pyfunction!(unregister_moon_spk, module)?)?;
    module.add_function(wrap_pyfunction!(list_registered_moons, module)?)?;
    module.add_function(wrap_pyfunction!(get_moon_name, module)?)?;
    module.add_function(wrap_pyfunction!(is_planetary_moon, module)?)?;
    module.add_function(wrap_pyfunction!(calc_moon_position, module)?)?;
    module.add_function(wrap_pyfunction!(get_moon_coverage, module)?)?;
    module.add_function(wrap_pyfunction!(close_moon_kernels, module)?)?;
    Ok(())
}
