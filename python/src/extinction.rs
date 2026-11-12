//! Thin atmospheric extinction, twilight and naked-eye visibility bindings.

use crate::errors::to_pyerr;
use libephemeris_core::errors::Error;
use libephemeris_core::extinction as native;
use pyo3::prelude::*;

/// Mutable atmospheric extinction components, in magnitudes per airmass.
///
/// calc_extinction_coefficient returns a component sum in k_total. Manual
/// construction and field assignment do not validate or recompute that sum.
#[pyclass(
    module = "libephemeris.extinction",
    get_all,
    set_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct ExtinctionCoefficients {
    /// Molecular (Rayleigh) extinction, magnitudes per airmass.
    k_rayleigh: f64,
    /// Aerosol extinction, magnitudes per airmass.
    k_aerosol: f64,
    /// Ozone absorption, magnitudes per airmass.
    k_ozone: f64,
    /// Water-vapour absorption, magnitudes per airmass.
    k_water: f64,
    /// Total extinction, magnitudes per airmass; independently mutable.
    k_total: f64,
}
#[pymethods]
impl ExtinctionCoefficients {
    /// Store five coefficients in mag/airmass without validation or summation.
    #[new]
    fn new(k_rayleigh: f64, k_aerosol: f64, k_ozone: f64, k_water: f64, k_total: f64) -> Self {
        Self {
            k_rayleigh,
            k_aerosol,
            k_ozone,
            k_water,
            k_total,
        }
    }
}
impl From<native::ExtinctionCoefficients> for ExtinctionCoefficients {
    fn from(value: native::ExtinctionCoefficients) -> Self {
        Self {
            k_rayleigh: value.k_rayleigh(),
            k_aerosol: value.k_aerosol(),
            k_ozone: value.k_ozone(),
            k_water: value.k_water(),
            k_total: value.k_total(),
        }
    }
}

/// Mutable twilight-model result with surface brightness and visual limit.
///
/// Higher surface brightness in mag/arcsec² means a darker sky. Manual field
/// changes do not update the other values or validate physical consistency.
#[pyclass(
    module = "libephemeris.extinction",
    get_all,
    set_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct TwilightSkyBrightness {
    /// Sky surface brightness, magnitudes per square arcsecond.
    surface_brightness: f64,
    /// Model phase: 'day', 'civil', 'nautical', 'astronomical' or 'night'.
    twilight_phase: String,
    /// Input Sun altitude above the horizon, degrees.
    sun_altitude_deg: f64,
    /// Dimensionless Sun-azimuth factor blended toward one as twilight fades.
    azimuth_factor: f64,
    /// Approximate naked-eye limiting magnitude, not surface brightness.
    limiting_magnitude: f64,
    /// Sky luminance in nanolamberts (not lamberts).
    nanolamberts: f64,
}
#[pymethods]
impl TwilightSkyBrightness {
    /// Store the result fields verbatim; units are documented on each property.
    ///
    /// No phase-token validation or brightness-unit conversion is performed.
    #[new]
    fn new(
        surface_brightness: f64,
        twilight_phase: String,
        sun_altitude_deg: f64,
        azimuth_factor: f64,
        limiting_magnitude: f64,
        nanolamberts: f64,
    ) -> Self {
        Self {
            surface_brightness,
            twilight_phase,
            sun_altitude_deg,
            azimuth_factor,
            limiting_magnitude,
            nanolamberts,
        }
    }
}
impl From<native::TwilightSkyBrightness> for TwilightSkyBrightness {
    fn from(value: native::TwilightSkyBrightness) -> Self {
        Self {
            surface_brightness: value.surface_brightness,
            twilight_phase: value.twilight_phase.as_str().to_string(),
            sun_altitude_deg: value.sun_altitude_deg,
            azimuth_factor: value.azimuth_factor,
            limiting_magnitude: value.limiting_magnitude,
            nanolamberts: value.nanolamberts,
        }
    }
}

/// Mutable result of the compact point-source visibility model.
///
/// This is not the full heliacal retinal/optical model. Manual assignment does
/// not re-evaluate the decision or enforce consistency between fields.
#[pyclass(
    module = "libephemeris.extinction",
    get_all,
    set_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct VisibilityResult {
    /// Whether effective object magnitude is strictly below the limiting magnitude.
    is_visible: bool,
    /// Effective object magnitude, including extinction when requested/applied.
    object_magnitude: f64,
    /// Faintest model magnitude; equality is not classified as visible.
    limiting_magnitude: f64,
    /// Input sky surface brightness, magnitudes per square arcsecond.
    sky_brightness: f64,
    /// Dimensionless object-to-one-arcsec²-sky-patch flux ratio.
    contrast: f64,
    /// Required dimensionless flux ratio, not a percentage or Blackwell contrast.
    threshold_contrast: f64,
    /// Limiting minus effective object magnitude; positive means visible.
    visibility_margin: f64,
    /// Canonical eye state: 'photopic', 'mesopic' or 'dark'.
    eye_adaptation: String,
    /// Skill ID: 1 novice, 2 average, 3 experienced, 4 expert; unknown IDs use average.
    observer_skill: i32,
}
#[pymethods]
impl VisibilityResult {
    /// Store the supplied decision and diagnostics without recomputing them.
    ///
    /// Magnitude, sky-brightness and ratio units are documented on the fields;
    /// eye_adaptation and observer_skill are stored without validation.
    #[new]
    #[allow(clippy::too_many_arguments)]
    fn new(
        is_visible: bool,
        object_magnitude: f64,
        limiting_magnitude: f64,
        sky_brightness: f64,
        contrast: f64,
        threshold_contrast: f64,
        visibility_margin: f64,
        eye_adaptation: String,
        observer_skill: i32,
    ) -> Self {
        Self {
            is_visible,
            object_magnitude,
            limiting_magnitude,
            sky_brightness,
            contrast,
            threshold_contrast,
            visibility_margin,
            eye_adaptation,
            observer_skill,
        }
    }
}
impl From<native::VisibilityResult> for VisibilityResult {
    fn from(value: native::VisibilityResult) -> Self {
        Self {
            is_visible: value.is_visible,
            object_magnitude: value.object_magnitude,
            limiting_magnitude: value.limiting_magnitude,
            sky_brightness: value.sky_brightness,
            contrast: value.contrast,
            threshold_contrast: value.threshold_contrast,
            visibility_margin: value.visibility_margin,
            eye_adaptation: value.eye_adaptation.as_str().to_string(),
            observer_skill: value.observer_skill,
        }
    }
}

fn eye(
    state: Option<&str>,
    entry: &'static str,
    py: Python<'_>,
) -> PyResult<Option<native::EyeAdaptation>> {
    state
        .map(|state| {
            native::normalize_eye_adaptation(state).ok_or_else(|| {
                to_pyerr(
                    py,
                    Error::InvalidGeneratorInput {
                        entry,
                        detail: "unknown eye adaptation state",
                    },
                )
            })
        })
        .transpose()
}

/// Calculate dimensionless relative optical airmass from altitude in degrees.
///
/// Parameters
/// ----------
/// altitude_deg : float
///     Object altitude above the horizon; physical domain 0 to 90 degrees.
/// method : str
///     Case-sensitive 'kasten_young' (default), 'secant' or 'rozenberg'. The
///     plane-parallel secant model is best used above roughly 15 degrees.
///
/// Returns
/// -------
/// float
///     Relative path length; 40 at or below the horizon. Secant also returns
///     40 at altitudes up to 5 degrees; the other formulae have their own clamps.
///
/// Raises
/// ------
/// libephemeris.InputValidationError
///     Unknown method. Numeric inputs are not physically range-validated.
///
/// Examples
/// --------
/// >>> from libephemeris.extinction import calc_airmass
/// >>> calc_airmass(45.0) > 1.0
/// True
#[pyfunction(signature=(altitude_deg,method="kasten_young"))]
fn calc_airmass(py: Python<'_>, altitude_deg: f64, method: &str) -> PyResult<f64> {
    let method = native::parse_airmass_method(method).ok_or_else(|| {
        to_pyerr(
            py,
            Error::InvalidGeneratorInput {
                entry: "calc_airmass",
                detail: "unknown airmass method",
            },
        )
    })?;
    let result = py.detach(|| native::calc_airmass(altitude_deg, method));
    Ok(result)
}

/// Return extinction components and their sum in magnitudes per airmass.
///
/// Parameters
/// ----------
/// pressure_mbar, temperature_c : float
///     Local pressure in mbar (hPa) and air temperature in Celsius.
/// humidity_percent : float
///     Relative humidity, normally 0 to 100 (not a fraction).
/// altitude_m : float
///     Observer height in metres.
/// wavelength_nm : float
///     Positive wavelength in nanometres; default 550 is V band.
/// visibility_km : float or None
///     Positive meteorological visibility in km overrides the aerosol term
///     using Koschmieder. None/nonpositive selects humidity/height modelling.
///
/// Returns
/// -------
/// ExtinctionCoefficients
///     Rayleigh, aerosol, ozone, water-vapour and total coefficients.
///
/// Notes
/// -----
/// A non-default pressure sets the Rayleigh air column; at exactly 1013.25
/// mbar the altitude sets it instead. Temperature is accepted but unused by
/// these component formulae. Numeric ranges/finiteness are not validated;
/// nonphysical inputs may produce nonfinite or nonphysical results.
#[pyfunction(signature=(pressure_mbar=1013.25,temperature_c=15.0,humidity_percent=50.0,altitude_m=0.0,wavelength_nm=550.0,visibility_km=None))]
fn calc_extinction_coefficient(
    py: Python<'_>,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    altitude_m: f64,
    wavelength_nm: f64,
    visibility_km: Option<f64>,
) -> PyResult<ExtinctionCoefficients> {
    let result = py.detach(|| {
        native::calc_extinction_coefficient(
            pressure_mbar,
            temperature_c,
            humidity_percent,
            altitude_m,
            wavelength_nm,
            visibility_km,
        )
    });
    Ok(result.into())
}

/// Return atmospheric dimming in magnitudes at altitude_deg above the horizon.
///
/// Uses the total extinction coefficient times Kasten-Young airmass. Altitude
/// is degrees, pressure_mbar is mbar, temperature_c is Celsius,
/// humidity_percent is percent (0..100), observer_altitude_m is metres and
/// wavelength_nm is a positive wavelength in nm. visibility_km is an optional
/// positive aerosol visibility override, as in calc_extinction_coefficient.
/// Returns 99.0 at or below zero altitude, not infinity or a missing value.
/// Numeric physical ranges are not validated; this is a scalar model, not an
/// ephemeris or refraction calculation.
#[pyfunction(signature=(altitude_deg,pressure_mbar=1013.25,temperature_c=15.0,humidity_percent=50.0,observer_altitude_m=0.0,wavelength_nm=550.0,visibility_km=None))]
#[allow(clippy::too_many_arguments)]
fn calc_extinction_magnitude(
    py: Python<'_>,
    altitude_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    observer_altitude_m: f64,
    wavelength_nm: f64,
    visibility_km: Option<f64>,
) -> PyResult<f64> {
    let result = py.detach(|| {
        native::calc_extinction_magnitude(
            altitude_deg,
            pressure_mbar,
            temperature_c,
            humidity_percent,
            observer_altitude_m,
            wavelength_nm,
            visibility_km,
        )
    });
    Ok(result)
}

/// Return simple plane-parallel extinction k/sin(altitude_deg), in magnitudes.
///
/// altitude_deg is degrees above the horizon; k is magnitudes per airmass.
/// Intended above roughly 15 degrees. Returns 99.0 at/below the horizon;
/// positive altitudes with sin(altitude) <= 0.001 return k*40. Numeric values
/// are not range-validated; this is not the Kasten-Young model.
#[pyfunction(signature=(altitude_deg,k=0.28))]
fn calc_simple_extinction(py: Python<'_>, altitude_deg: f64, k: f64) -> PyResult<f64> {
    let result = py.detach(|| native::calc_simple_extinction(altitude_deg, k));
    Ok(result)
}

/// Add model extinction to an extra-atmospheric catalog magnitude.
///
/// Returns catalog_magnitude + calc_extinction_magnitude(...). altitude_deg
/// is degrees above the horizon; pressure_mbar is mbar, temperature_c Celsius,
/// humidity_percent percent, observer_altitude_m metres and wavelength_nm nm
/// (positive; default V band). At/below the horizon 99 magnitudes are added.
/// Numeric physical ranges are not validated; no position is acquired.
#[pyfunction(signature=(catalog_magnitude,altitude_deg,pressure_mbar=1013.25,temperature_c=15.0,humidity_percent=50.0,observer_altitude_m=0.0,wavelength_nm=550.0))]
#[allow(clippy::too_many_arguments)]
fn apparent_magnitude_with_extinction(
    py: Python<'_>,
    catalog_magnitude: f64,
    altitude_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    observer_altitude_m: f64,
    wavelength_nm: f64,
) -> PyResult<f64> {
    let result = py.detach(|| {
        native::apparent_magnitude_with_extinction(
            catalog_magnitude,
            altitude_deg,
            pressure_mbar,
            temperature_c,
            humidity_percent,
            observer_altitude_m,
            wavelength_nm,
        )
    });
    Ok(result)
}

/// Return the total V-band extinction coefficient in magnitudes per airmass.
///
/// zenith_angle_deg is accepted but unused: this returns a coefficient, not
/// the extinction along that line of sight. pressure_mbar is mbar,
/// temperature_c Celsius (unused by these components), humidity_percent
/// percent and altitude_m metres. Uses 550 nm and no visibility override;
/// physical numeric ranges are not validated.
#[pyfunction(signature=(zenith_angle_deg,pressure_mbar=1013.25,temperature_c=15.0,humidity_percent=50.0,altitude_m=0.0))]
fn get_extinction_for_heliacal(
    py: Python<'_>,
    zenith_angle_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    altitude_m: f64,
) -> PyResult<f64> {
    let result = py.detach(|| {
        native::get_extinction_for_heliacal(
            zenith_angle_deg,
            pressure_mbar,
            temperature_c,
            humidity_percent,
            altitude_m,
        )
    });
    Ok(result)
}

/// Return molecular scattering extinction in magnitudes per airmass.
///
/// wavelength_nm is positive nm; pressure_mbar is mbar and altitude_m metres.
/// A pressure different from 1013.25 scales the air column directly; otherwise
/// height supplies an exponential scale. temperature_c is Celsius but unused.
/// The wavelength dependence is proportional to wavelength**-4. No physical
/// range validation is performed; zero wavelength can produce infinity.
#[pyfunction(signature=(wavelength_nm=550.0,pressure_mbar=1013.25,temperature_c=15.0,altitude_m=0.0))]
fn calc_rayleigh_coefficient(
    py: Python<'_>,
    wavelength_nm: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    altitude_m: f64,
) -> PyResult<f64> {
    let result = py.detach(|| {
        native::calc_rayleigh_coefficient(wavelength_nm, pressure_mbar, temperature_c, altitude_m)
    });
    Ok(result)
}

/// Return aerosol extinction in magnitudes per airmass.
///
/// humidity_percent is percent, altitude_m metres and wavelength_nm positive
/// nm. Positive visibility_km uses max(3.912/visibility_km - 0.1066, 0) and
/// bypasses the other operands. None/nonpositive visibility selects the
/// humidity/height/Ångström wavelength model. Numeric ranges are not validated.
#[pyfunction(signature=(humidity_percent=50.0,altitude_m=0.0,wavelength_nm=550.0,visibility_km=None))]
fn calc_aerosol_coefficient(
    py: Python<'_>,
    humidity_percent: f64,
    altitude_m: f64,
    wavelength_nm: f64,
    visibility_km: Option<f64>,
) -> PyResult<f64> {
    let result = py.detach(|| {
        native::calc_aerosol_coefficient(humidity_percent, altitude_m, wavelength_nm, visibility_km)
    });
    Ok(result)
}

/// Return approximate ozone absorption in magnitudes per airmass.
///
/// wavelength_nm is a positive wavelength in nm. A step model changes at
/// 400, 500 and 700 nm; it is not a resolved absorption spectrum. No numeric
/// range validation is performed.
#[pyfunction(signature=(wavelength_nm=550.0))]
fn calc_ozone_coefficient(py: Python<'_>, wavelength_nm: f64) -> PyResult<f64> {
    let result = py.detach(|| native::calc_ozone_coefficient(wavelength_nm));
    Ok(result)
}

/// Return water-vapour absorption in magnitudes per airmass.
///
/// humidity_percent is percent (not a fraction); wavelength_nm is positive
/// nm. Piecewise coefficients at 600/800 nm scale linearly with humidity.
/// temperature_c is accepted in Celsius but unused. Physical ranges are not
/// validated; this is an approximate band model, not a resolved spectrum.
#[pyfunction(signature=(humidity_percent=50.0,temperature_c=15.0,wavelength_nm=550.0))]
fn calc_water_vapor_coefficient(
    py: Python<'_>,
    humidity_percent: f64,
    temperature_c: f64,
    wavelength_nm: f64,
) -> PyResult<f64> {
    let result = py.detach(|| {
        native::calc_water_vapor_coefficient(humidity_percent, temperature_c, wavelength_nm)
    });
    Ok(result)
}

/// Classify Sun altitude in degrees as a twilight phase string.
///
/// Returns 'day' above 0, 'civil' above -6, 'nautical' above -12,
/// 'astronomical' above -18, otherwise 'night'. Exact boundaries select the
/// darker phase. This performs no clock conversion or ephemeris calculation.
#[pyfunction(signature=(sun_altitude_deg))]
fn get_twilight_phase(py: Python<'_>, sun_altitude_deg: f64) -> PyResult<&'static str> {
    let result = py.detach(|| native::get_twilight_phase(sun_altitude_deg));
    Ok(result.as_str())
}

/// Estimate twilight sky brightness in a supplied horizontal direction.
///
/// Parameters
/// ----------
/// sun_altitude_deg, target_altitude_deg : float
///     Sun and target altitude above the horizon, degrees.
/// sun_azimuth_deg, target_azimuth_deg : float
///     Azimuth degrees in the same convention, normally 0..360. Only their
///     angular difference is used; defaults point opposite the Sun.
/// pressure_mbar, temperature_c, humidity_percent : float
///     Pressure mbar, temperature Celsius and humidity percent, not fraction.
/// altitude_m, wavelength_nm : float
///     Observer height metres and positive wavelength nm (default 550).
///
/// Returns
/// -------
/// TwilightSkyBrightness
///     Surface brightness (mag/arcsec², clamped 0..22.5), phase, Sun altitude,
///     azimuth factor, approximate limiting magnitude (-2..7), nanolamberts.
///
/// Notes
/// -----
/// Empirical twilight-anchor interpolation with azimuth, altitude and aerosol
/// adjustments; not the full heliacal optical model. No ephemeris/UTC clock
/// is needed. Numeric ranges are not validated and visibility is not guaranteed.
#[pyfunction(signature=(sun_altitude_deg,target_altitude_deg=90.0,sun_azimuth_deg=0.0,target_azimuth_deg=180.0,pressure_mbar=1013.25,temperature_c=15.0,humidity_percent=50.0,altitude_m=0.0,wavelength_nm=550.0))]
#[allow(clippy::too_many_arguments)]
fn calc_twilight_sky_brightness(
    py: Python<'_>,
    sun_altitude_deg: f64,
    target_altitude_deg: f64,
    sun_azimuth_deg: f64,
    target_azimuth_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    altitude_m: f64,
    wavelength_nm: f64,
) -> PyResult<TwilightSkyBrightness> {
    let result = py.detach(|| {
        native::calc_twilight_sky_brightness(
            sun_altitude_deg,
            target_altitude_deg,
            sun_azimuth_deg,
            target_azimuth_deg,
            pressure_mbar,
            temperature_c,
            humidity_percent,
            altitude_m,
            wavelength_nm,
        )
    });
    Ok(result.into())
}

/// Return approximate zenith twilight surface brightness, in mag/arcsec².
///
/// sun_altitude_deg is degrees above the horizon. Evaluates
/// calc_twilight_sky_brightness at the zenith, opposite the Sun, with standard
/// atmosphere and V band. Higher returned magnitude means a darker sky.
#[pyfunction(signature=(sun_altitude_deg))]
fn calc_twilight_brightness_simple(py: Python<'_>, sun_altitude_deg: f64) -> PyResult<f64> {
    let result = py.detach(|| native::calc_twilight_brightness_simple(sun_altitude_deg));
    Ok(result)
}

/// Estimate a naked-eye limiting magnitude during twilight, clamped -2..6.5.
///
/// Sun/target altitudes and azimuths are degrees in a common horizontal
/// convention. pressure_mbar is mbar, humidity_percent percent and altitude_m
/// observer height in metres. Temperature is fixed at 15 Celsius and wavelength
/// at 550 nm. Subtracts k_total*(airmass-1) from the twilight brightness fit;
/// not the full heliacal detection model and not a below-horizon visibility
/// test. Numeric physical ranges are not validated.
#[pyfunction(signature=(sun_altitude_deg,target_altitude_deg=45.0,sun_azimuth_deg=0.0,target_azimuth_deg=180.0,pressure_mbar=1013.25,humidity_percent=50.0,altitude_m=0.0))]
#[allow(clippy::too_many_arguments)]
fn calc_limiting_magnitude_twilight(
    py: Python<'_>,
    sun_altitude_deg: f64,
    target_altitude_deg: f64,
    sun_azimuth_deg: f64,
    target_azimuth_deg: f64,
    pressure_mbar: f64,
    humidity_percent: f64,
    altitude_m: f64,
) -> PyResult<f64> {
    let result = py.detach(|| {
        native::calc_limiting_magnitude_twilight(
            sun_altitude_deg,
            target_altitude_deg,
            sun_azimuth_deg,
            target_azimuth_deg,
            pressure_mbar,
            humidity_percent,
            altitude_m,
        )
    });
    Ok(result)
}

/// Return 'photopic', 'mesopic' or 'dark' from sky surface brightness.
///
/// sky_brightness is mag/arcsec² (not nanolamberts). Values below 16 select
/// photopic, below 20 mesopic, and at/above 20 dark-adapted vision. This is a
/// brightness classification, not a time-dependent adaptation simulation.
#[pyfunction(signature=(sky_brightness))]
fn calc_eye_adaptation_state(py: Python<'_>, sky_brightness: f64) -> PyResult<&'static str> {
    let result = py.detach(|| native::calc_eye_adaptation_state(sky_brightness));
    Ok(result.as_str())
}

/// Return the model's dimensionless point-source detection flux ratio.
///
/// Parameters
/// ----------
/// sky_brightness : float
///     Surface brightness in mag/arcsec²; higher means darker.
/// eye_adaptation : str or None
///     None auto-selects; exact tokens 'photopic', 'mesopic', 'dark' and
///     'scotopic' (alias of dark) override it.
/// observer_skill : int
///     1 novice, 2 average, 3 experienced, 4 expert; other IDs use average.
///
/// Returns
/// -------
/// float
///     Required object-to-one-arcsec²-sky-patch flux ratio, clamped 1e-8..1e10.
///     Not a percentage or classical Blackwell contrast fraction.
///
/// Raises
/// ------
/// libephemeris.InputValidationError
///     Unknown eye_adaptation token. Numeric ranges are not validated.
#[pyfunction(signature=(sky_brightness,eye_adaptation=None,observer_skill=2))]
fn calc_contrast_threshold(
    py: Python<'_>,
    sky_brightness: f64,
    eye_adaptation: Option<&str>,
    observer_skill: i32,
) -> PyResult<f64> {
    let eye_adaptation = eye(eye_adaptation, "calc_contrast_threshold", py)?;
    let result = py
        .detach(|| native::calc_contrast_threshold(sky_brightness, eye_adaptation, observer_skill));
    Ok(result)
}

/// Estimate point-source visibility against a supplied sky surface brightness.
///
/// Parameters
/// ----------
/// object_magnitude, sky_brightness : float
///     Object magnitude and sky surface brightness in mag/arcsec² respectively.
/// eye_adaptation : str or None
///     None auto-selects; 'photopic', 'mesopic', 'dark' or 'scotopic'.
/// observer_skill : int
///     1 novice, 2 average, 3 experienced, 4 expert; unknown IDs use average.
/// object_altitude_deg : float
///     Object altitude above horizon in degrees.
/// apply_extinction : bool
///     Add V-band extinction only when True and altitude is strictly positive.
/// pressure_mbar, humidity_percent, observer_altitude_m : float
///     Pressure mbar, humidity percent and height metres; extinction uses 15 C.
///
/// Returns
/// -------
/// VisibilityResult
///     Effective magnitude, limiting magnitude, sky brightness, flux ratios,
///     margin, adaptation/skill and strict magnitude < limit decision.
///
/// Raises
/// ------
/// libephemeris.InputValidationError
///     Unknown adaptation token.
///
/// Notes
/// -----
/// Compact model treating sky flux as a one-arcsec² patch, not the full
/// heliacal retinal/optics model. Does not reject a below-horizon object;
/// altitude only controls the optional extinction correction. Other numeric
/// ranges are not validated; this does not guarantee observational visibility.
#[pyfunction(signature=(object_magnitude,sky_brightness,eye_adaptation=None,observer_skill=2,object_altitude_deg=90.0,apply_extinction=false,pressure_mbar=1013.25,humidity_percent=50.0,observer_altitude_m=0.0))]
#[allow(clippy::too_many_arguments)]
fn calc_visibility_threshold(
    py: Python<'_>,
    object_magnitude: f64,
    sky_brightness: f64,
    eye_adaptation: Option<&str>,
    observer_skill: i32,
    object_altitude_deg: f64,
    apply_extinction: bool,
    pressure_mbar: f64,
    humidity_percent: f64,
    observer_altitude_m: f64,
) -> PyResult<VisibilityResult> {
    let eye_adaptation = eye(eye_adaptation, "calc_visibility_threshold", py)?;
    let result = py.detach(|| {
        native::calc_visibility_threshold(
            object_magnitude,
            sky_brightness,
            eye_adaptation,
            observer_skill,
            object_altitude_deg,
            apply_extinction,
            pressure_mbar,
            humidity_percent,
            observer_altitude_m,
        )
    });
    Ok(result.into())
}

/// Return only the boolean decision from calc_visibility_threshold.
///
/// object_magnitude is a point-source magnitude; sky_brightness is
/// mag/arcsec², observer_skill is 1..4 (unknown IDs use average), and
/// object_altitude_deg is degrees. Adaptation is automatic. If apply_extinction
/// is True and altitude > 0, uses 1013.25 mbar, 15 C, 50% humidity, sea level
/// and V band. Equality at the limit is False. There is no horizon rejection;
/// the result is an approximate brightness decision, not an event forecast.
#[pyfunction(signature=(object_magnitude,sky_brightness,observer_skill=2,object_altitude_deg=90.0,apply_extinction=false))]
fn is_object_visible(
    py: Python<'_>,
    object_magnitude: f64,
    sky_brightness: f64,
    observer_skill: i32,
    object_altitude_deg: f64,
    apply_extinction: bool,
) -> PyResult<bool> {
    let result = py.detach(|| {
        native::is_object_visible(
            object_magnitude,
            sky_brightness,
            observer_skill,
            object_altitude_deg,
            apply_extinction,
        )
    });
    Ok(result)
}

/// Return the compact model's limiting point-source magnitude for a sky.
///
/// sky_brightness is mag/arcsec²; observer_skill is 1 novice, 2 average,
/// 3 experienced or 4 expert (unknown IDs use average). eye_adaptation is None
/// for automatic selection or 'photopic', 'mesopic', 'dark'/'scotopic'. Unknown
/// tokens raise InputValidationError. Uses no extinction correction and does
/// not acquire ephemeris data. Numeric physical ranges are not validated.
///
/// Examples
/// --------
/// >>> from libephemeris.extinction import calc_limiting_magnitude_for_sky
/// >>> limit = calc_limiting_magnitude_for_sky(21.5, observer_skill=3)
/// >>> isinstance(limit, float)
/// True
#[pyfunction(signature=(sky_brightness,observer_skill=2,eye_adaptation=None))]
fn calc_limiting_magnitude_for_sky(
    py: Python<'_>,
    sky_brightness: f64,
    observer_skill: i32,
    eye_adaptation: Option<&str>,
) -> PyResult<f64> {
    let eye_adaptation = eye(eye_adaptation, "calc_limiting_magnitude_for_sky", py)?;
    let result = py.detach(|| {
        native::calc_limiting_magnitude_for_sky(sky_brightness, observer_skill, eye_adaptation)
    });
    Ok(result)
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_class::<ExtinctionCoefficients>()?;
    module.add_class::<TwilightSkyBrightness>()?;
    module.add_class::<VisibilityResult>()?;
    module.add_function(wrap_pyfunction!(calc_airmass, module)?)?;
    module.add_function(wrap_pyfunction!(calc_extinction_coefficient, module)?)?;
    module.add_function(wrap_pyfunction!(calc_extinction_magnitude, module)?)?;
    module.add_function(wrap_pyfunction!(calc_simple_extinction, module)?)?;
    module.add_function(wrap_pyfunction!(
        apparent_magnitude_with_extinction,
        module
    )?)?;
    module.add_function(wrap_pyfunction!(get_extinction_for_heliacal, module)?)?;
    module.add_function(wrap_pyfunction!(calc_rayleigh_coefficient, module)?)?;
    module.add_function(wrap_pyfunction!(calc_aerosol_coefficient, module)?)?;
    module.add_function(wrap_pyfunction!(calc_ozone_coefficient, module)?)?;
    module.add_function(wrap_pyfunction!(calc_water_vapor_coefficient, module)?)?;
    module.add_function(wrap_pyfunction!(get_twilight_phase, module)?)?;
    module.add_function(wrap_pyfunction!(calc_twilight_sky_brightness, module)?)?;
    module.add_function(wrap_pyfunction!(calc_twilight_brightness_simple, module)?)?;
    module.add_function(wrap_pyfunction!(calc_limiting_magnitude_twilight, module)?)?;
    module.add_function(wrap_pyfunction!(calc_eye_adaptation_state, module)?)?;
    module.add_function(wrap_pyfunction!(calc_contrast_threshold, module)?)?;
    module.add_function(wrap_pyfunction!(calc_visibility_threshold, module)?)?;
    module.add_function(wrap_pyfunction!(is_object_visible, module)?)?;
    module.add_function(wrap_pyfunction!(calc_limiting_magnitude_for_sky, module)?)?;
    Ok(())
}
