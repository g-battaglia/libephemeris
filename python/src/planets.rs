// SPDX-License-Identifier: AGPL-3.0
//! Typed orbital, photometric, sidereal and body-metadata calls.

use crate::errors::to_pyerr;
use crate::positions::{PositionResult, position, runtime};
use libephemeris_core::python_support::{Ayanamsha, Session};
use libephemeris_core::{erfa, errors, python_support as api, sidereal};
use pyo3::exceptions::{PyUserWarning, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use std::ffi::CString;

pyo3::create_exception!(
    libephemeris,
    NutationFallbackWarning,
    PyUserWarning,
    "Warning category for nutation-fallback advisories; a UserWarning subclass."
);
pyo3::create_exception!(
    libephemeris,
    HeliocentricNodApsWarning,
    PyUserWarning,
    "Warning category for heliocentric node/apse advisories; a UserWarning subclass."
);

type PosTuple = (f64, f64, f64, f64, f64, f64);
type NodApsResult = (PosTuple, PosTuple, PosTuple, PosTuple);

fn six(values: [f64; 6]) -> PosTuple {
    let [a, b, c, d, e, f] = values;
    (a, b, c, d, e, f)
}

fn session_call<T: Send>(
    py: Python<'_>,
    operation: impl FnOnce(&Session) -> errors::Result<T> + Send,
) -> PyResult<T> {
    py.detach(|| operation(runtime().shared_session()?.as_ref()))
        .map_err(|error| to_pyerr(py, error))
}

fn warn(py: Python<'_>, message: &str) -> PyResult<()> {
    if !message.is_empty() {
        let message =
            CString::new(message).map_err(|error| PyValueError::new_err(error.to_string()))?;
        PyErr::warn(py, &py.get_type::<PyUserWarning>(), &message, 1)?;
    }
    Ok(())
}

/// Calculate a target as seen from another body's center at TT Julian day.
///
/// Parameters
/// ----------
/// tjdet : float
///     TT Julian day. Unlike EphemerisContext.calc_pctr, this function does
///     not take UT1.
/// planet, center : int
///     Target and observer body IDs: Sun, Moon, Mercury through Pluto or Earth.
///     Nodes, stars, asteroids and planetary moons are unsupported endpoints.
/// flags : int, optional
///     FLG_* bitmask, default FLG_SWIEPH | FLG_SPEED. FLG_TOPOCTR does not
///     introduce a surface observer; center remains explicit. FLG_TRUEPOS
///     selects geometric separation; NOABERR/NOGDEFL disable their corrections.
///
/// Returns
/// -------
/// tuple
///     (values, returned_flags). The six values follow calc: spherical
///     longitude/latitude/distance and rates in degrees, AU and per TT day;
///     EQUATORIAL, RADIANS or XYZ select the corresponding representation.
///     FLG_SPEED enables rates; FLG_SPEED3 alone does not enable them here.
///
/// Raises
/// ------
/// ValueError
///     For nonfinite epochs on checked paths.
/// Error
///     For unsupported endpoints, missing data/configuration, uncovered
///     target/support samples or numerical failures. No extrapolation occurs.
///
/// Notes
/// -----
/// Coincident endpoints use the native zero-distance convention. A central
/// date inside a file need not cover all retarded or rate samples.
///
/// Examples
/// --------
/// With data already configured:
/// ```python
/// import libephemeris as eph
/// values, flags = eph.calc_pctr(2451545.0, eph.MARS, eph.EARTH)
/// ```
#[pyfunction(signature = (tjdet, planet, center, flags=258))]
fn calc_pctr(
    py: Python<'_>,
    tjdet: f64,
    planet: i32,
    center: i32,
    flags: i32,
) -> PyResult<PositionResult> {
    session_call(py, |session| {
        session.calc_pctr(tjdet, planet, center, flags)
    })
    .and_then(|value| position(py, value))
}

/// Calculate ascending/descending nodes and perihelion/aphelion at TT JD.
///
/// tjdet is TT Julian day; planet is a body ID. method combines NODBIT_*
/// choices: NODBIT_MEAN selects mean elements where supported, osculating
/// choices use instantaneous geometry, and NODBIT_FOPOINT selects the second
/// focus instead of aphelion. flags is the FLG_* calculation/output bitmask.
///
/// Returns (ascending_node, descending_node, perihelion, aphelion_or_focus),
/// four six-value tuples, without a returned-flags element. Position/rate
/// components follow calc's representation: degrees/AU and rates per TT day,
/// or the requested equatorial, radians or XYZ output. FLG_SPEED enables rates.
/// Undefined interpolated-apse points retain NaNs. Raises Error for unsupported
/// objects or data/coverage/configuration failures, and ValueError for checked
/// nonfinite inputs. See nod_aps_ut for UT1 input.
#[pyfunction(signature = (tjdet, planet, method=1, flags=258))]
fn nod_aps(
    py: Python<'_>,
    tjdet: f64,
    planet: i32,
    method: i32,
    flags: i32,
) -> PyResult<NodApsResult> {
    session_call(py, |session| session.nod_aps(tjdet, planet, method, flags)).map(|result| {
        (
            six(result.ascending),
            six(result.descending),
            six(result.perihelion),
            six(result.aphelion),
        )
    })
}

/// Calculate nodes and apsides from a UT1 Julian day.
///
/// tjdut is UT1, not UTC; planet, method and flags follow nod_aps. Returns
/// four six-component tuples in ascending-node, descending-node, perihelion,
/// aphelion/second-focus order, with no returned-flags element. Angular/length
/// units and rates per TT day follow nod_aps, including undefined-point NaNs.
/// UT1 observation conversion uses the independent bundled timescale; model
/// Delta-T settings still apply to nested model calculations. Data, coverage,
/// configuration and unsupported-object errors propagate as in nod_aps.
#[pyfunction(signature = (tjdut, planet, method=1, flags=258))]
fn nod_aps_ut(
    py: Python<'_>,
    tjdut: f64,
    planet: i32,
    method: i32,
    flags: i32,
) -> PyResult<NodApsResult> {
    session_call(py, |session| {
        session.nod_aps_ut(tjdut, planet, method, flags)
    })
    .map(|result| {
        (
            six(result.ascending),
            six(result.descending),
            six(result.perihelion),
            six(result.aphelion),
        )
    })
}

/// Return a body's 50-component osculating-orbital-element tuple at TT JD.
///
/// Parameters
/// ----------
/// tjdet : float
///     TT Julian day.
/// planet : int
///     Supported planet, lunar, minor-body or hypothetical-body identifier.
/// flags : int
///     Calculation/center flags, including FLG_HELCTR, FLG_BARYCTR and
///     FLG_ORBEL_AA where applicable. This is not a six-coordinate calc result;
///     the element units below do not change to radians with FLG_RADIANS.
///
/// Returns
/// -------
/// tuple of float
///     Indices and units (all angles are degrees):
///
/// ```text
/// 0: semi-major axis (AU); 1: eccentricity; 2: inclination;
/// 3: ascending-node longitude; 4: perihelion argument;
/// 5: perihelion longitude; 6: mean anomaly; 7: true anomaly;
/// 8: eccentric anomaly; 9: mean longitude;
/// 10: sidereal period (years of 365.24219 days);
/// 11: mean motion (degrees/day);
/// 12: tropical period (years of 365.24219 days);
/// 13: signed synodic period (days); 14: perihelion-passage TT JD;
/// 15: perihelion distance (AU); 16: aphelion distance (AU);
/// 17-49: reserved zeros. No returned-flags component is appended.
/// ```
///
/// Notes
/// -----
/// Synodic period can be negative for an interior orbit; degenerate geometry
/// retains native zero, infinity or NaN values. These are instantaneous
/// two-body elements, not a promise of long-term orbit stability.
///
/// Raises
/// ------
/// Error
///     For unsupported objects, missing support/configuration or coverage and
///     calculation failures.
/// ValueError
///     For checked nonfinite numerical inputs.
#[pyfunction]
fn get_orbital_elements<'py>(
    py: Python<'py>,
    tjdet: f64,
    planet: i32,
    flags: i32,
) -> PyResult<Bound<'py, PyTuple>> {
    let values = py
        .detach(|| runtime().get_orbital_elements(tjdet, planet, flags))
        .map_err(|error| to_pyerr(py, error))?;
    PyTuple::new(py, values)
}

/// Return the 50 orbital elements from a UT1 Julian day.
///
/// tjd_ut is UT1, ipl is the body ID and iflag is the calculation bitmask.
/// The tuple's indices, units, reserved zeros and exceptions are documented
/// in get_orbital_elements; there is no appended returned-flags value.
/// UT1 is converted with the independent observation timescale; nested model
/// calculations retain their configured clocks. This is not a UTC interface.
#[pyfunction]
fn get_orbital_elements_ut<'py>(
    py: Python<'py>,
    tjd_ut: f64,
    ipl: i32,
    iflag: i32,
) -> PyResult<Bound<'py, PyTuple>> {
    let values = py
        .detach(|| runtime().get_orbital_elements_ut(tjd_ut, ipl, iflag))
        .map_err(|error| to_pyerr(py, error))?;
    PyTuple::new(py, values)
}

/// Return maximum, minimum and current true orbital distance at TT JD.
///
/// tjdet is TT Julian day; planet is a supported body ID; flags selects
/// calculation/center/output behavior. Returns (maximum, minimum, current)
/// in AU for ordinary spherical output. Extrema use osculating-ellipse
/// geometry, not a dated event search. The current value is calc component
/// two with FLG_TRUEPOS added: FLG_XYZ therefore returns signed z there,
/// not a radius. No returned flags are appended. Unsupported objects,
/// missing data/configuration and coverage/calculation failures raise Error;
/// checked nonfinite input raises ValueError.
#[pyfunction]
fn orbit_max_min_true_distance(
    py: Python<'_>,
    tjdet: f64,
    planet: i32,
    flags: i32,
) -> PyResult<(f64, f64, f64)> {
    py.detach(|| runtime().orbit_max_min_true_distance(tjdet, planet, flags))
        .map_err(|error| to_pyerr(py, error))
}

/// Calculate physical phenomena for a body at TT Julian day.
///
/// tjdet is TT JD; planet is a body ID; flags is the calculation/correction
/// bitmask (default FLG_SWIEPH). FLG_TOPOCTR requires a configured observer.
///
/// Returns a 20-float tuple, without returned flags:
/// 0 phase angle (degrees); 1 illuminated fraction; 2 solar elongation
/// (degrees); 3 apparent angular diameter (degrees); 4 visual magnitude;
/// 5 lunar horizontal parallax (degrees, zero for other bodies);
/// 6-19 unused zeros. Angular units remain degrees, not FLG_RADIANS output.
/// Point/degenerate bodies retain their native NaN/zero sentinels. Magnitudes
/// use the body's physical photometric model, not a fixed-star catalogue field.
/// Missing data/observer configuration, coverage and calculation failures raise
/// Error subclasses; checked nonfinite numerical input raises ValueError.
/// See pheno_ut for UT1 input.
#[pyfunction(signature = (tjdet, planet, flags=2))]
fn pheno<'py>(
    py: Python<'py>,
    tjdet: f64,
    planet: i32,
    flags: i32,
) -> PyResult<Bound<'py, PyTuple>> {
    let values = py
        .detach(|| runtime().pheno(tjdet, planet, flags))
        .map_err(|error| to_pyerr(py, error))?
        .values;
    PyTuple::new(py, values)
}

/// Calculate physical phenomena for a body at UT1 Julian day.
///
/// tjdut is UT1, not UTC; planet is a body ID and flags is the calculation
/// bitmask. Returns the same 20-float tuple as pheno: phase angle, illuminated
/// fraction, solar elongation, diameter, magnitude, lunar parallax, then 14
/// unused zeros. Angles are degrees; there is no returned-flags element.
/// Observer requirements, body-specific sentinels and exceptions follow pheno.
///
/// Examples
/// --------
/// With an ephemeris already configured:
/// ```python
/// import libephemeris as eph
/// attributes = eph.pheno_ut(2451545.0, eph.VENUS)
/// phase_angle, illuminated_fraction = attributes[:2]
/// visual_magnitude = attributes[4]
/// ```
#[pyfunction(signature = (tjdut, planet, flags=2))]
fn pheno_ut<'py>(
    py: Python<'py>,
    tjdut: f64,
    planet: i32,
    flags: i32,
) -> PyResult<Bound<'py, PyTuple>> {
    let values = py
        .detach(|| runtime().pheno_ut(tjdut, planet, flags))
        .map_err(|error| to_pyerr(py, error))?
        .values;
    PyTuple::new(py, values)
}

/// Return the signed target-minus-Sun longitude difference and evening flag.
///
/// tjd_ut is UT1 JD, ipl is the body ID and iflag is passed to both calc_ut
/// calls. Returns (signed_difference, is_evening), with is_evening true
/// exactly when the difference is positive. Ordinary spherical output gives
/// degrees, folded by one turn outside [-180, 180]; this is a longitude
/// difference, not the great-circle solar elongation returned by pheno_ut.
/// Use ordinary spherical degrees: RADIANS/XYZ retain literal component-zero
/// subtraction and the same numeric fold thresholds, not a converted angle.
/// Position acquisition/observer/coverage errors propagate from calc_ut.
#[pyfunction(signature = (tjd_ut, ipl, iflag=0))]
fn get_elongation_from_sun(
    py: Python<'_>,
    tjd_ut: f64,
    ipl: i32,
    iflag: i32,
) -> PyResult<(f64, bool)> {
    session_call(py, |session| {
        session.get_elongation_from_sun(tjd_ut, ipl, iflag)
    })
}

/// Return the signed target-minus-Sun longitude difference at UT1 JD.
///
/// tjd_ut, ipl and iflag follow get_elongation_from_sun. Returns its first
/// component alone, normally degrees; its component-subtraction flag caveat
/// and propagated calculation exceptions apply unchanged.
#[pyfunction(signature = (tjd_ut, ipl, iflag=0))]
fn get_signed_elongation(py: Python<'_>, tjd_ut: f64, ipl: i32, iflag: i32) -> PyResult<f64> {
    session_call(py, |session| {
        session.get_signed_elongation(tjd_ut, ipl, iflag)
    })
}

/// Test the morning-star longitude relation at UT1 Julian day tjd_ut.
///
/// ipl is a body ID; iflag follows get_elongation_from_sun, including its
/// component-subtraction caveat. Returns False for the Sun; other bodies
/// negate the positive evening predicate, so zero or unordered differences
/// return True. This is not an atmospheric visibility prediction. Position
/// errors propagate, except that the Sun returns before calculation.
#[pyfunction(signature = (tjd_ut, ipl, iflag=0))]
fn is_morning_star(py: Python<'_>, tjd_ut: f64, ipl: i32, iflag: i32) -> PyResult<bool> {
    session_call(py, |session| session.is_morning_star(tjd_ut, ipl, iflag))
}

/// Test the evening-star longitude relation at UT1 Julian day tjd_ut.
///
/// ipl is a body ID; iflag follows get_elongation_from_sun. Returns True
/// for a strictly positive target-minus-Sun component-zero difference and
/// False for the Sun. This does not predict visibility; ordinary spherical
/// degree flags are appropriate. Position errors propagate for other bodies.
#[pyfunction(signature = (tjd_ut, ipl, iflag=0))]
fn is_evening_star(py: Python<'_>, tjd_ut: f64, ipl: i32, iflag: i32) -> PyResult<bool> {
    session_call(py, |session| session.is_evening_star(tjd_ut, ipl, iflag))
}

/// Classify a body's longitude relation to the Sun at UT1 Julian day.
///
/// tjd_ut, ipl and iflag follow get_elongation_from_sun. Returns 'eastern'
/// for a positive signed difference, 'western' otherwise, or 'none' for the
/// Sun. This is not a visibility classification. Calculation exceptions and
/// the spherical-degree/component-subtraction caveat follow that function.
#[pyfunction(signature = (tjd_ut, ipl, iflag=0))]
fn get_elongation_type(
    py: Python<'_>,
    tjd_ut: f64,
    ipl: i32,
    iflag: i32,
) -> PyResult<&'static str> {
    session_call(py, |session| {
        session.get_elongation_type(tjd_ut, ipl, iflag)
    })
}

/// Return metadata for the native IAU 2006/2000A nutation model.
///
/// Returns a dict with model (str), terms (int), precision (str), source
/// (str) and skyfield_available (False). The last key is retained for
/// compatibility; no Python astronomical engine is used. precision describes
/// the published model, not a measured error bound at an arbitrary epoch.
#[pyfunction]
fn get_nutation_model(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
    let info = py.detach(erfa::nutation_model_info);
    let result = PyDict::new(py);
    result.set_item("model", info.model)?;
    result.set_item("terms", info.terms)?;
    result.set_item("precision", info.precision)?;
    result.set_item("source", info.source)?;
    // Keep the public key while reporting the native engine truthfully.
    result.set_item("skyfield_available", false)?;
    Ok(result)
}

/// Configure the default runtime's sidereal mode and optional defining epoch.
///
/// Parameters
/// ----------
/// mode : int
///     SIDM_* identifier combined with supported SIDBIT_* projection bits.
/// t0 : float, optional
///     Defining Julian day for SIDM_USER: TT unless SIDBIT_USER_UT is set,
///     in which case it is UT1. Stored as supplied, not a universal reset date.
/// ayan_t0 : float, optional
///     Ayanamsha at t0 in degrees for the user-defined mode.
///
/// Notes
/// -----
/// Use FLG_SIDEREAL on calculations to request sidereal output. Unknown modes
/// or unsupported projections can issue UserWarning and follow the native
/// fallback/ignored-bit policy; this setter does not validate epochs.
///
/// Examples
/// --------
/// ```python
/// import libephemeris as eph
/// eph.set_sid_mode(eph.SIDM_LAHIRI)
/// ayanamsha_deg = eph.get_ayanamsa(2451545.0)
/// ```
#[pyfunction(signature = (mode, t0=0.0, ayan_t0=0.0))]
fn set_sid_mode(py: Python<'_>, mode: i32, t0: f64, ayan_t0: f64) -> PyResult<()> {
    let warnings = py.detach(|| runtime().set_sid_mode(mode, t0, ayan_t0));
    for warning in warnings {
        warn(py, &warning.to_string())?;
    }
    Ok(())
}

fn ayanamsha(
    py: Python<'_>,
    epoch: f64,
    flags: i32,
    is_ut: bool,
    extended: bool,
) -> PyResult<Ayanamsha> {
    let result = session_call(py, |session| {
        session.public_ayanamsha(epoch, flags, is_ut, extended)
    })?;
    warn(py, &result.diagnostic)?;
    Ok(result)
}

/// Return the configured mean ayanamsha in degrees at UT1 JD tjdut.
///
/// The selected mode comes from set_sid_mode. Unlike get_ayanamsa_ex_ut,
/// this scalar excludes the extended entry's nutation term and returns no
/// flags. Dynamic modes may require stellar/data support; acquisition,
/// coverage and configuration failures raise Error subclasses. Unknown-mode
/// fallback advisories can issue UserWarning. Checked nonfinite epochs raise
/// ValueError. Use get_ayanamsa for TT input.
#[pyfunction]
fn get_ayanamsa_ut(py: Python<'_>, tjdut: f64) -> PyResult<f64> {
    ayanamsha(py, tjdut, 0, true, false).map(|result| result.value)
}

/// Return the configured mean ayanamsha in degrees at TT JD tjdet.
///
/// Uses set_sid_mode and returns a float without flags or the extended
/// entry's nutation term. Dynamic defining anchors may require ephemeris
/// data. Data/configuration/coverage errors propagate as Error subclasses;
/// unknown-mode advisories can issue UserWarning. Checked nonfinite epochs
/// raise ValueError. See get_ayanamsa_ut for UT1.
#[pyfunction]
fn get_ayanamsa(py: Python<'_>, tjdet: f64) -> PyResult<f64> {
    ayanamsha(py, tjdet, 0, false, false).map(|result| result.value)
}

/// Return (returned_flags, ayanamsha_degrees) at TT Julian day tjdet.
///
/// flags selects the extended ayanamsha correction choices. Nutation is
/// included unless FLG_NONUT is set; FLG_TRUEPOS/FLG_NOABERR/FLG_NOGDEFL
/// affect applicable dynamic anchors. The result remains degrees even with
/// FLG_RADIANS. Note the flags-first order, unlike calc. Mode selection,
/// propagated data/configuration/coverage errors and UserWarning advisories
/// follow get_ayanamsa.
#[pyfunction(signature = (tjdet, flags=0))]
fn get_ayanamsa_ex(py: Python<'_>, tjdet: f64, flags: i32) -> PyResult<(i32, f64)> {
    ayanamsha(py, tjdet, flags, false, true).map(|result| (result.returned_flags, result.value))
}

/// Return (returned_flags, ayanamsha_degrees) at UT1 Julian day tjdut.
///
/// flags and the degrees-only, flags-first result follow get_ayanamsa_ex:
/// nutation is included unless FLG_NONUT is set. This is UT1, not UTC.
/// The configured sidereal mode, data requirements, propagated exceptions
/// and possible UserWarning advisories are shared with get_ayanamsa_ut.
#[pyfunction(signature = (tjdut, flags=0))]
fn get_ayanamsa_ex_ut(py: Python<'_>, tjdut: f64, flags: i32) -> PyResult<(i32, f64)> {
    ayanamsha(py, tjdut, flags, true, true).map(|result| (result.returned_flags, result.value))
}

/// Look up the predefined ayanamsha name for integer sidmode.
///
/// Returns str, or '' for an unassigned/user-defined mode. SIDBIT_* high
/// bits do not affect this low-byte lookup. No mode is changed or data loaded.
#[pyfunction]
fn get_ayanamsa_name(py: Python<'_>, sidmode: i32) -> &'static str {
    py.detach(|| sidereal::ayanamsha_name(sidmode))
}

/// Return the native display name for an integer library body identifier.
///
/// planet may identify a planet, lunar point or another supported body family.
/// Returns str, or '' for an unknown selector. This is a naming lookup, not
/// a position/coverage check or a network search.
#[pyfunction]
fn get_planet_name(py: Python<'_>, planet: i32) -> String {
    py.detach(|| api::get_planet_name(planet))
}

/// Remove compatibility request bits that are consumed without being echoed.
///
/// flags is an integer bitmask and planet a library body ID. Returns an int:
/// JPLHOR/JPLHOR_APPROX are removed, as is CENTER_BODY for Sun through Mars.
/// Other bits are preserved. This pure helper neither calculates a position
/// nor validates the full request or establishes its physical provider.
#[pyfunction(signature=(flags, planet))]
fn consume_non_echoed_flags(py: Python<'_>, flags: i32, planet: i32) -> i32 {
    py.detach(|| libephemeris_core::positions::consume_non_echoed_flags(flags, planet))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = module.py();
    let nutation_warning = py.get_type::<NutationFallbackWarning>();
    nutation_warning.setattr("__module__", "libephemeris.planets")?;
    module.add("NutationFallbackWarning", nutation_warning)?;
    let nodaps_warning = py.get_type::<HeliocentricNodApsWarning>();
    nodaps_warning.setattr("__module__", "libephemeris.planets")?;
    module.add("HeliocentricNodApsWarning", nodaps_warning)?;
    module.add_function(wrap_pyfunction!(calc_pctr, module)?)?;
    module.add_function(wrap_pyfunction!(nod_aps, module)?)?;
    module.add_function(wrap_pyfunction!(nod_aps_ut, module)?)?;
    module.add_function(wrap_pyfunction!(get_orbital_elements, module)?)?;
    module.add_function(wrap_pyfunction!(get_orbital_elements_ut, module)?)?;
    module.add_function(wrap_pyfunction!(orbit_max_min_true_distance, module)?)?;
    module.add_function(wrap_pyfunction!(pheno, module)?)?;
    module.add_function(wrap_pyfunction!(pheno_ut, module)?)?;
    module.add_function(wrap_pyfunction!(get_elongation_from_sun, module)?)?;
    module.add_function(wrap_pyfunction!(get_signed_elongation, module)?)?;
    module.add_function(wrap_pyfunction!(is_morning_star, module)?)?;
    module.add_function(wrap_pyfunction!(is_evening_star, module)?)?;
    module.add_function(wrap_pyfunction!(get_elongation_type, module)?)?;
    module.add_function(wrap_pyfunction!(get_nutation_model, module)?)?;
    module.add_function(wrap_pyfunction!(set_sid_mode, module)?)?;
    module.add_function(wrap_pyfunction!(get_ayanamsa_ut, module)?)?;
    module.add_function(wrap_pyfunction!(get_ayanamsa, module)?)?;
    module.add_function(wrap_pyfunction!(get_ayanamsa_ex, module)?)?;
    module.add_function(wrap_pyfunction!(get_ayanamsa_ex_ut, module)?)?;
    module.add_function(wrap_pyfunction!(get_ayanamsa_name, module)?)?;
    module.add_function(wrap_pyfunction!(get_planet_name, module)?)?;
    module.add_function(wrap_pyfunction!(consume_non_echoed_flags, module)?)?;
    Ok(())
}
