//! Thin heliacal visibility and event-search bindings.

use crate::errors::to_pyerr;
use crate::positions::runtime;
use libephemeris_core::{heliacal, python_support};
use pyo3::prelude::*;
use pyo3::types::PyTuple;

/// Return True exactly for the Mercury or Venus body identifier.
///
/// body is an integer library ID. This identity test acquires no ephemeris
/// and does not validate any other identifier.
#[pyfunction]
fn is_inner_planet(body: i32) -> bool {
    heliacal::is_inner_planet(body)
}

/// Return whether integer body >= FIXSTAR_OFFSET.
///
/// This tests the identifier range only; it does not establish that a catalog
/// row exists or acquire any ephemeris data.
#[pyfunction]
fn is_fixed_star(body: i32) -> bool {
    heliacal::is_fixed_star(body)
}

/// Search for a heliacal event and its visibility interval (sequence API).
///
/// Parameters
/// ----------
/// tjdut : float
///     Search start, UT1 Julian day, not TT or UTC.
/// geopos : sequence of float
///     (longitude east degrees, latitude north degrees, height metres).
///     Missing words default to zero; extra words are ignored.
/// atmo : sequence of float
///     (pressure mbar, temperature Celsius, humidity percent, met_range).
///     Missing words default to (1013.25, 15, 40, 0). Nonpositive search
///     pressure selects 1013.25. Range 0 selects the model default, 0<range<1
///     supplies total extinction mag/airmass, >=1 means visibility kilometres.
/// observer : sequence of float
///     (age years, Snellen ratio, binocular, magnification, aperture mm,
///     transmission fraction). Missing defaults: (36, 1, 0, 0, 0, 1).
///     The last four words are used only with HELFLAG_OPTICAL_PARAMS.
/// objname : str
///     Planet name or catalog-star name. This named entry does not parse
///     decimal body-ID strings. Sun is not a valid search target.
/// eventtype : int
///     1 morning first, 2 evening last, 3 evening first, 4 morning last.
///     Moon accepts 3/4 only; 3/4 otherwise require Mercury or Venus. Stars
///     and outer planets accept 1/2.
/// flags : int
///     Heliacal HELFLAG_* options (including OPTICAL_PARAMS, HIGH_PRECISION,
///     NO_DETAILS), not arbitrary position-output flags. Default 2.
///
/// Returns
/// -------
/// tuple of float
///     (visibility_start, optimum, visibility_end), all UT1 JDs. With
///     HELFLAG_NO_DETAILS: (detected_event_jd, 0, 0). A completed 800-day miss
///     is (0, 0, 0), not an exception. Detailed refinement may fail separately.
///
/// Raises
/// ------
/// ValueError, libephemeris.Error
///     Invalid event/name/site or consumed nonfinite optical/atmospheric
///     input, data/coverage failure, or unusable detailed visibility interval.
///
/// Notes
/// -----
/// Root libephemeris.heliacal_ut uses this sequence interface, unlike
/// libephemeris.eclipse.heliacal_ut's compact scalar (time, event_type) model.
/// Native apparent geometry uses independent observation/site clocks;
/// photometry retains the configured model clock. Stored mode retains its
/// selected reader. Bounded scans/refinement do not guarantee convergence or
/// actual observational visibility. No observer setting is changed.
#[pyfunction(signature=(tjdut,geopos,atmo,observer,objname,eventtype,flags=2))]
#[allow(clippy::too_many_arguments)]
fn heliacal_ut(
    py: Python<'_>,
    tjdut: f64,
    geopos: Vec<f64>,
    atmo: Vec<f64>,
    observer: Vec<f64>,
    objname: &str,
    eventtype: i32,
    flags: i32,
) -> PyResult<(f64, f64, f64)> {
    py.detach(|| {
        let session = runtime().shared_session()?;
        python_support::heliacal_ut_words(
            &session, tjdut, &geopos, &atmo, &observer, objname, eventtype, flags,
        )
    })
    .map(|[a, b, c]| (a, b, c))
    .map_err(|error| to_pyerr(py, error))
}

/// Return fifty fixed-instant heliacal phenomena and local-window values.
///
/// Parameters
/// ----------
/// tjdut : float
///     UT1 Julian day, not UTC/TT.
/// geopos, atmo, observer : sequence of float
///     Layouts/defaults as heliacal_ut: east longitude, north latitude, metres;
///     mbar, Celsius, humidity percent (default 40), range km/direct coefficient;
///     age years, Snellen ratio and flag-enabled optics (aperture mm,
///     transmission fraction). Missing words use those defaults.
/// objname : str
///     Planet or catalog-star name (not a decimal ID string).
/// eventtype : int
///     Codes 1..4 select the associated rise/set visibility window: 1/4 rise,
///     2/3 set. Unlike the outer event search, unsupported body/window pairings
///     may still return instantaneous geometry.
/// flags : int
///     Heliacal options, default 2; not arbitrary coordinate-output flags.
///
/// Returns
/// -------
/// tuple of float
///     Fifty slots, zero-based; angular quantities below are degrees:
///     0 geometric topocentric altitude; 1 refracted apparent altitude;
///     2 geocentric altitude; 3 object azimuth; 4/5 Sun altitude/azimuth;
///     6 topocentric vertical arc; 7 geocentric vertical arc; 8 signed
///     Sun-minus-object azimuth [-180,180]; 9 arc of light; 10 extinction
///     mag/airmass; 11 required vertical arc; 12/13/14 visibility start/optimum/
///     end UT1 JDs; 15 lunar Yallop optimum UT1 JD; 16 crescent width degrees;
///     17 Yallop q; 18 class code 1..6 as float; 19 altitude parallax proxy;
///     20 visual magnitude (catalog for stars); 21/22 object/Sun disc-centre
///     rise or set UT1 JDs; 23 crossing lag days; 24 visibility duration days;
///     25 crescent length degrees; 26 solar elongation degrees; 27 illumination
///     PERCENT, not fraction; 28..49 reserved zeros. Azimuths start north and
///     increase eastward, unlike azalt's south-west convention.
///
/// Notes
/// -----
/// Moon-only slots 16..18/25 are zero otherwise. Missing window/crossing times
/// and unavailable durations retain TJD_INVALID, not zero or a valid date.
/// Missing lag is zero; unsupported 3/4 windows have duration zero. Stars have
/// zero parallax and 100% illumination. Selected numerical photometry failures
/// can use angular separation with zero magnitude/phase; data/range failures
/// propagate. This is not an 800-day event search or convergence guarantee.
/// Invalid event/name/site and consumed input/data/window failures raise their
/// native Python exceptions. Clock/provider distinctions match heliacal_ut.
#[pyfunction(signature=(tjdut,geopos,atmo,observer,objname,eventtype,flags=2))]
#[allow(clippy::too_many_arguments)]
fn heliacal_pheno_ut<'py>(
    py: Python<'py>,
    tjdut: f64,
    geopos: Vec<f64>,
    atmo: Vec<f64>,
    observer: Vec<f64>,
    objname: &str,
    eventtype: i32,
    flags: i32,
) -> PyResult<Bound<'py, PyTuple>> {
    let result = py
        .detach(|| {
            let session = runtime().shared_session()?;
            python_support::heliacal_pheno_ut_words(
                &session, tjdut, &geopos, &atmo, &observer, objname, eventtype, flags,
            )
        })
        .map_err(|error| to_pyerr(py, error))?;
    PyTuple::new(py, result)
}

/// Return the VISLIMIT vision flag and ten limiting-magnitude detail slots.
///
/// Parameters
/// ----------
/// tjdut : float
///     Fixed UT1 Julian day, not TT/UTC.
/// geopos : sequence of float
///     (east longitude degrees, north latitude degrees, height metres),
///     missing words default zero.
/// atmo : sequence of float
///     (pressure mbar, temperature Celsius, humidity PERCENT, met_range).
///     Missing defaults (1013.25, 15, 50, 0), unlike search's 40% humidity.
///     Explicit zero pressure stays zero (altitude-form molecular extinction).
///     Range 0 is default, 0<range<1 direct total mag/airmass, >=1 km.
/// observer : sequence of float
///     Age years/Snellen ratio default 36/1. With HELFLAG_OPTICAL_PARAMS,
///     additional words are binocular, magnification, aperture mm and
///     transmission fraction, default 0/0/0/1; otherwise ignored.
/// objname : str
///     Planet name, supported numeric-ID string, or catalog-star name.
/// flags : int
///     Heliacal/visibility bitmask. DARK/NOMOON suppress Sun/Moon background;
///     PHOTOPIC/SCOTOPIC can override the detected vision regime.
///
/// Returns
/// -------
/// tuple
///     (vision_flag_as_float, details10): 0 limiting catalog V magnitude,
///     1 object altitude, 2 object azimuth, 3 Sun altitude, 4 Sun azimuth,
///     5 Moon altitude, 6 Moon azimuth, 7 object visual magnitude,
///     8/9 reserved zero. Angles are degrees; azimuths start north and increase
///     eastward (including the converted stored/star horizon results).
///     Below horizon: HELFLAG_BELOW_HORIZON and (-100., 0., ..., 0.).
///
/// Notes
/// -----
/// Limiting magnitude already includes extinction; compare directly with the
/// catalog magnitude. This is not eclipse.vis_limit_mag's eight-slot compact
/// model. Data/provider, name/site and consumed nonfinite-input failures
/// propagate; no visibility/convergence guarantee is supplied.
#[pyfunction(signature=(tjdut,geopos,atmo,observer,objname,flags=2))]
fn vis_limit_mag<'py>(
    py: Python<'py>,
    tjdut: f64,
    geopos: Vec<f64>,
    atmo: Vec<f64>,
    observer: Vec<f64>,
    objname: &str,
    flags: i32,
) -> PyResult<(f64, Bound<'py, PyTuple>)> {
    let (flag, values) = py
        .detach(|| {
            let session = runtime().shared_session()?;
            python_support::vis_limit_mag_words(
                &session, tjdut, &geopos, &atmo, &observer, objname, flags,
            )
        })
        .map_err(|error| to_pyerr(py, error))?;
    Ok((flag, PyTuple::new(py, values)?))
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(is_inner_planet, module)?)?;
    module.add_function(wrap_pyfunction!(is_fixed_star, module)?)?;
    module.add_function(wrap_pyfunction!(heliacal_ut, module)?)?;
    module.add_function(wrap_pyfunction!(heliacal_pheno_ut, module)?)?;
    module.add_function(wrap_pyfunction!(vis_limit_mag, module)?)?;
    Ok(())
}
