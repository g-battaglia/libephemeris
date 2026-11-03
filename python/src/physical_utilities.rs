//! Direct bindings for scalar atmospheric and house-construction utilities.

use crate::errors::to_pyerr;
use libephemeris_core::{houses, refraction, schaefer};
use pyo3::exceptions::{PyOverflowError, PyValueError, PyZeroDivisionError};
use pyo3::prelude::*;
use pyo3::types::PyDict;

fn schaefer_error(error: schaefer::SchaeferError) -> PyErr {
    match error {
        schaefer::SchaeferError::MathDomain => PyValueError::new_err(error.to_string()),
        schaefer::SchaeferError::MathRange | schaefer::SchaeferError::PowerRange => {
            PyOverflowError::new_err(error.to_string())
        }
        schaefer::SchaeferError::DivisionByZero => PyZeroDivisionError::new_err(error.to_string()),
    }
}

fn submodule<'py>(parent: &Bound<'py, PyModule>, name: &str) -> PyResult<Bound<'py, PyModule>> {
    let full_name = format!("libephemeris._core._{name}");
    let module = PyModule::new(parent.py(), &full_name)?;
    parent
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item(full_name, &module)?;
    parent.add(format!("_{name}"), &module)?;
    Ok(module)
}

#[allow(non_snake_case)]
mod schaefer_bindings {
    use super::*;
    /// Return molecular V-band extinction in magnitudes per airmass.
    ///
    /// pressure_mbar is local pressure in mbar/hPa; altitude_m is site height
    /// in metres. Uses ``0.1451 * (pressure/1013.25) * exp(-height/8500)``.
    /// Nonnegative pressure is physical but not enforced. Finite exponential
    /// overflow raises OverflowError; other IEEE NaN/infinity arithmetic may
    /// propagate. This column model reads no local atmosphere or ephemeris.
    #[pyfunction(signature=(pressure_mbar, altitude_m))]
    fn calc_rayleigh_extinction(
        py: Python<'_>,
        pressure_mbar: f64,
        altitude_m: f64,
    ) -> PyResult<f64> {
        py.detach(|| schaefer::calc_rayleigh_extinction(pressure_mbar, altitude_m))
            .map_err(schaefer_error)
    }
    /// Return aerosol extinction in magnitudes per airmass.
    ///
    /// humidity_fraction normally lies in `0..=1`. met_range_km defaults to 0:
    /// values >= 1 are meteorological range in kilometres; values strictly
    /// between 0 and 1 are an already local extinction coefficient, returned
    /// unchanged; other values select the humidity polynomial. altitude_m
    /// (default 0 metres) scales range/humidity estimates, NOT direct inputs.
    /// No physical bounds are imposed. Finite power/exponential overflow
    /// raises OverflowError; ordinary NaN/infinity arithmetic can propagate.
    #[pyfunction(signature=(humidity_fraction, met_range_km=0.0, altitude_m=0.0))]
    fn calc_aerosol_extinction(
        py: Python<'_>,
        humidity_fraction: f64,
        met_range_km: f64,
        altitude_m: f64,
    ) -> PyResult<f64> {
        py.detach(|| schaefer::calc_aerosol_extinction(humidity_fraction, met_range_km, altitude_m))
            .map_err(schaefer_error)
    }
    /// Return the fixed representative V-band ozone extinction, 0.016 mag/airmass.
    ///
    /// Takes no arguments/data and models no seasonal or ozone-column changes.
    #[pyfunction(signature=())]
    fn calc_ozone_extinction(py: Python<'_>) -> PyResult<f64> {
        Ok(py.detach(schaefer::calc_ozone_extinction))
    }
    /// Return water-vapour extinction ``0.014 * humidity_fraction`` in mag/airmass.
    ///
    /// Humidity normally lies in `0..=1` but is not clamped or validated.
    /// Nonfinite values follow IEEE arithmetic; no atmosphere/data is read.
    #[pyfunction(signature=(humidity_fraction))]
    fn calc_water_vapor_extinction(py: Python<'_>, humidity_fraction: f64) -> PyResult<f64> {
        Ok(py.detach(|| schaefer::calc_water_vapor_extinction(humidity_fraction)))
    }
    /// Sum Rayleigh, aerosol, ozone and water extinction in mag/airmass.
    ///
    /// pressure_mbar defaults to 1013.25 mbar/hPa, temperature_c to 15 Celsius,
    /// humidity_fraction to 0.5, met_range_km to 0, altitude_m to 0 metres.
    /// Temperature is accepted but UNUSED. met_range_km has the dual range/
    /// coefficient meaning of calc_aerosol_extinction. Inputs are not clamped
    /// to physical ranges. Component overflow raises OverflowError; ordinary
    /// nonfinite arithmetic may propagate. This is a representative V-band
    /// atmosphere, not a measured vertical profile or a full visibility model.
    #[pyfunction(signature=(pressure_mbar=1013.25, temperature_c=15.0, humidity_fraction=0.5, met_range_km=0.0, altitude_m=0.0))]
    fn calc_total_extinction(
        py: Python<'_>,
        pressure_mbar: f64,
        temperature_c: f64,
        humidity_fraction: f64,
        met_range_km: f64,
        altitude_m: f64,
    ) -> PyResult<f64> {
        py.detach(|| {
            schaefer::calc_total_extinction(
                pressure_mbar,
                temperature_c,
                humidity_fraction,
                met_range_km,
                altitude_m,
            )
        })
        .map_err(schaefer_error)
    }
    /// Return Kasten-Young relative airmass for object altitude_deg in degrees.
    ///
    /// The dimensionless result is clamped to `1..=40`: altitude <= 0 returns
    /// 40 and altitude >= 90 returns 1, including infinities. NaN propagates.
    /// No refraction, atmosphere acquisition or physical-input error is added.
    #[pyfunction(signature=(altitude_deg))]
    fn calc_airmass(py: Python<'_>, altitude_deg: f64) -> PyResult<f64> {
        Ok(py.detach(|| schaefer::calc_airmass(altitude_deg)))
    }
    /// Return full-path object dimming ``k_total * airmass`` in magnitudes.
    ///
    /// altitude_deg is object altitude in degrees; k_total is mag/airmass.
    /// Uses calc_airmass, including its horizon clamps and NaN propagation.
    /// No physical validation or additional exception is imposed.
    #[pyfunction(signature=(altitude_deg, k_total))]
    fn calc_extinction_magnitude(py: Python<'_>, altitude_deg: f64, k_total: f64) -> PyResult<f64> {
        Ok(py.detach(|| schaefer::calc_extinction_magnitude(altitude_deg, k_total)))
    }
    /// Return heuristic twilight sky brightness in nanoLamberts (nL).
    ///
    /// sun_altitude_deg is solar altitude in degrees: >= 0 returns 1e10 nL,
    /// <= -18 returns zero; intermediate civil/nautical/astronomical branches
    /// use project curves. NaN propagates, and branch joins need not be
    /// continuous. No azimuth, aerosol profile, clock or event is modeled.
    #[pyfunction(signature=(sun_altitude_deg))]
    fn calc_twilight_brightness(py: Python<'_>, sun_altitude_deg: f64) -> PyResult<f64> {
        Ok(py.detach(|| schaefer::calc_twilight_brightness(sun_altitude_deg)))
    }
    /// Return the simplified moonlight sky-brightness estimate in nL.
    ///
    /// All three inputs are DEGREES: Moon altitude, phase angle (0 = full,
    /// 180 = new, not an illuminated fraction) and object-Moon separation.
    /// Moon at/below the horizon returns zero before reading other inputs;
    /// separation below 1 degree is raised to 1. Infinite phase in the active
    /// branch raises ValueError; NaNs can propagate. This project heuristic
    /// is not the separate VISLIMIT scattered-moonlight model.
    #[pyfunction(signature=(moon_altitude_deg, moon_phase, angular_distance_deg))]
    fn calc_moon_brightness(
        py: Python<'_>,
        moon_altitude_deg: f64,
        moon_phase: f64,
        angular_distance_deg: f64,
    ) -> PyResult<f64> {
        py.detach(|| {
            schaefer::calc_moon_brightness(moon_altitude_deg, moon_phase, angular_distance_deg)
        })
        .map_err(schaefer_error)
    }
    /// Return constant natural-night airglow brightness, 145 nanoLamberts.
    ///
    /// No arguments or errors; seasonal, geographic and temporal variability
    /// are outside this scalar convenience model.
    #[pyfunction(signature=())]
    fn calc_airglow_brightness(py: Python<'_>) -> PyResult<f64> {
        Ok(py.detach(schaefer::calc_airglow_brightness))
    }
    /// Return heuristic zodiacal sky brightness in nanoLamberts.
    ///
    /// ecliptic_latitude_deg and sun_elongation_deg are degrees. Uses an
    /// exponential latitude factor and piecewise sine elongation factor,
    /// zero below 20 degrees. Inputs are not physically bounded: elongations
    /// beyond 180 can yield negative brightness. Infinite elongation in a
    /// sine branch raises ValueError; NaN arithmetic propagates. No dust
    /// distribution, observed sky or ephemeris is acquired.
    #[pyfunction(signature=(ecliptic_latitude_deg, sun_elongation_deg))]
    fn calc_zodiacal_brightness(
        py: Python<'_>,
        ecliptic_latitude_deg: f64,
        sun_elongation_deg: f64,
    ) -> PyResult<f64> {
        py.detach(|| schaefer::calc_zodiacal_brightness(ecliptic_latitude_deg, sun_elongation_deg))
            .map_err(schaefer_error)
    }
    /// Sum twilight, moonlight, airglow and ramped zodiacal brightness in nL.
    ///
    /// All inputs are degrees. sun_altitude_deg is required; defaults are
    /// moon_altitude_deg=-90, moon_phase_deg=180 (new Moon), moon_distance_deg
    /// =90 (object-Moon angular separation), ecliptic_latitude_deg=0 and
    /// sun_elongation_deg=90. Zodiacal weight ramps from zero at Sun -12 to
    /// one at -18 degrees, but the component is evaluated even at zero weight.
    /// Consumed infinite trigonometric arguments raise ValueError; component
    /// NaNs can propagate. These are convenience curves, not a measured sky
    /// or an ephemeris-driven visibility/event calculation.
    #[pyfunction(signature=(sun_altitude_deg, moon_altitude_deg=-90.0, moon_phase_deg=180.0, moon_distance_deg=90.0, ecliptic_latitude_deg=0.0, sun_elongation_deg=90.0),text_signature="(sun_altitude_deg, moon_altitude_deg=-90.0, moon_phase_deg=180.0, moon_distance_deg=90.0, ecliptic_latitude_deg=0.0, sun_elongation_deg=90.0)")]
    fn calc_total_sky_brightness(
        py: Python<'_>,
        sun_altitude_deg: f64,
        moon_altitude_deg: f64,
        moon_phase_deg: f64,
        moon_distance_deg: f64,
        ecliptic_latitude_deg: f64,
        sun_elongation_deg: f64,
    ) -> PyResult<f64> {
        py.detach(|| {
            schaefer::calc_total_sky_brightness(
                sun_altitude_deg,
                moon_altitude_deg,
                moon_phase_deg,
                moon_distance_deg,
                ecliptic_latitude_deg,
                sun_elongation_deg,
            )
        })
        .map_err(schaefer_error)
    }
    /// Convert nL brightness to surface brightness in magnitudes/arcsecond^2.
    ///
    /// Returns ``26.33 - 2.5*log10(brightness_nl)`` for positive brightness;
    /// nonpositive input returns the 22-mag/arcsecond^2 dark-sky sentinel.
    /// NaN propagates and positive infinity returns negative infinity. This
    /// scalar conversion performs no spectral correction or domain refusal.
    #[pyfunction(signature=(brightness_nl))]
    fn brightness_to_mag_arcsec2(py: Python<'_>, brightness_nl: f64) -> PyResult<f64> {
        Ok(py.detach(|| schaefer::brightness_to_mag_arcsec2(brightness_nl)))
    }
    /// Convert surface magnitude per square arcsecond to nanoLamberts.
    ///
    /// mag_arcsec2 maps to ``10**((26.33 - mag_arcsec2)/2.5)``. This inverts
    /// the positive-brightness branch, not its nonpositive sentinel. Finite
    /// power overflow raises OverflowError; NaNs/infinities follow arithmetic
    /// and very faint finite values may underflow to zero.
    #[pyfunction(signature=(mag_arcsec2))]
    fn mag_arcsec2_to_brightness(py: Python<'_>, mag_arcsec2: f64) -> PyResult<f64> {
        py.detach(|| schaefer::mag_arcsec2_to_brightness(mag_arcsec2))
            .map_err(schaefer_error)
    }
    /// Estimate the faintest visible INTRINSIC naked-eye magnitude.
    ///
    /// sky_brightness_nl is nL, extinction_coeff mag/airmass, object_altitude_deg
    /// degrees, observer_age years (default 25), snellen_ratio dimensionless
    /// visual acuity (default 1). The linear convenience model subtracts sky,
    /// age, acuity and ``k*(airmass-1)`` losses from 6.5; extinction is already
    /// included in the intrinsic threshold. It is not the full retinal model.
    /// Inputs are not physically clamped; nonpositive Snellen takes no acuity
    /// loss. A nonpositive logarithm argument raises ValueError (including
    /// infinite positive Snellen); ordinary IEEE infinities/NaNs can propagate.
    #[pyfunction(signature=(sky_brightness_nl, extinction_coeff, object_altitude_deg, observer_age=25.0, snellen_ratio=1.0))]
    fn calc_limiting_magnitude(
        py: Python<'_>,
        sky_brightness_nl: f64,
        extinction_coeff: f64,
        object_altitude_deg: f64,
        observer_age: f64,
        snellen_ratio: f64,
    ) -> PyResult<f64> {
        py.detach(|| {
            schaefer::calc_limiting_magnitude(
                sky_brightness_nl,
                extinction_coeff,
                object_altitude_deg,
                observer_age,
                snellen_ratio,
            )
        })
        .map_err(schaefer_error)
    }
    /// Estimate instantaneous naked-eye visibility from supplied scalar conditions.
    ///
    /// Parameters
    /// ----------
    /// object_magnitude : float
    ///     Intrinsic magnitude, before atmospheric dimming.
    /// object_altitude_deg, sun_altitude_deg : float
    ///     Object/Sun altitudes in degrees.
    /// pressure_mbar, temperature_c, humidity_fraction, met_range_km, altitude_m
    ///     Defaults 1013.25 mbar/hPa, 15 Celsius (unused), 0.5 relative humidity,
    ///     0 dual-use aerosol range/coefficient and 0 metres site height.
    /// observer_age, snellen_ratio
    ///     Defaults 25 years and 1 (dimensionless acuity).
    /// moon_altitude_deg, moon_phase_deg, sun_elongation_deg, include_moon
    ///     Defaults -90 degrees, 180 degrees (new Moon), 90 degrees and True.
    ///     Disabling moonlight does not disable the other sky components.
    ///
    /// Returns
    /// -------
    /// (visible, intrinsic_limit_mag, apparent_object_mag)
    ///     bool and two floats. At object altitude <= 0 returns (False,0,m);
    ///     otherwise Sun altitude >= 0 returns (False,-5,m). Active night
    ///     evaluation compares INTRINSIC m to the extinction-adjusted limit
    ///     and applies magnitude-dependent 7/10/15-degree elongation gates.
    ///     Apparent magnitude m+k*airmass is informational, not compared again.
    ///
    /// Notes
    /// -----
    /// Uses fixed object-Moon distance 90 degrees and ecliptic latitude 0.
    /// Physical ranges are not enforced. Consumed component domain failures
    /// raise ValueError and finite power/exponential overflow OverflowError;
    /// early horizon/daylight returns skip those components. No epoch or
    /// visibility event is solved, and nonfinite arithmetic may propagate.
    #[allow(clippy::too_many_arguments)]
    #[pyfunction(signature=(object_magnitude, object_altitude_deg, sun_altitude_deg, pressure_mbar=1013.25, temperature_c=15.0, humidity_fraction=0.5, met_range_km=0.0, altitude_m=0.0, observer_age=25.0, snellen_ratio=1.0, moon_altitude_deg=-90.0, moon_phase_deg=180.0, sun_elongation_deg=90.0, include_moon=true),text_signature="(object_magnitude, object_altitude_deg, sun_altitude_deg, pressure_mbar=1013.25, temperature_c=15.0, humidity_fraction=0.5, met_range_km=0.0, altitude_m=0.0, observer_age=25.0, snellen_ratio=1.0, moon_altitude_deg=-90.0, moon_phase_deg=180.0, sun_elongation_deg=90.0, include_moon=True)")]
    fn is_object_visible(
        py: Python<'_>,
        object_magnitude: f64,
        object_altitude_deg: f64,
        sun_altitude_deg: f64,
        pressure_mbar: f64,
        temperature_c: f64,
        humidity_fraction: f64,
        met_range_km: f64,
        altitude_m: f64,
        observer_age: f64,
        snellen_ratio: f64,
        moon_altitude_deg: f64,
        moon_phase_deg: f64,
        sun_elongation_deg: f64,
        include_moon: bool,
    ) -> PyResult<(bool, f64, f64)> {
        py.detach(|| {
            schaefer::is_object_visible(
                object_magnitude,
                object_altitude_deg,
                sun_altitude_deg,
                pressure_mbar,
                temperature_c,
                humidity_fraction,
                met_range_km,
                altitude_m,
                observer_age,
                snellen_ratio,
                moon_altitude_deg,
                moon_phase_deg,
                sun_elongation_deg,
                include_moon,
            )
        })
        .map_err(schaefer_error)
    }
    /// Return a project magnitude-band arcus-visionis default in degrees.
    ///
    /// object_magnitude is intrinsic magnitude. Table bins cover [-5,6);
    /// values >= 6 return 17 degrees, otherwise the out-of-table result is 4
    /// degrees (also for NaN). No event, atmospheric model or physical guard
    /// is evaluated; the historical table name is not an accuracy claim.
    #[pyfunction(signature=(object_magnitude))]
    fn get_arcus_visionis(py: Python<'_>, object_magnitude: f64) -> PyResult<f64> {
        Ok(py.detach(|| schaefer::get_arcus_visionis(object_magnitude)))
    }
    /// Return heuristic optimal solar altitude in degrees for an intrinsic magnitude.
    ///
    /// is_morning is required but UNUSED: both directions use the same
    /// magnitude bins, returning -8/-9/-10/-11/-12 degrees. NaN takes -12.
    /// There is no epoch, event solver, atmosphere acquisition or input error.
    #[pyfunction(signature=(object_magnitude, is_morning))]
    fn get_optimal_sun_altitude(
        py: Python<'_>,
        object_magnitude: f64,
        is_morning: bool,
    ) -> PyResult<f64> {
        Ok(py.detach(|| schaefer::get_optimal_sun_altitude(object_magnitude, is_morning)))
    }
    /// Return the heuristic minimum heliacal object altitude in degrees.
    ///
    /// object_magnitude is intrinsic, sun_altitude_deg degrees and
    /// extinction_coeff mag/airmass. Scales the magnitude-band arcus by the
    /// extinction and solar-altitude factors, with an ordered one-degree
    /// floor. No physical validation or event search occurs; arithmetic NaNs
    /// select the floor and infinities may propagate without an exception.
    #[pyfunction(signature=(object_magnitude, sun_altitude_deg, extinction_coeff))]
    fn calc_heliacal_altitude_threshold(
        py: Python<'_>,
        object_magnitude: f64,
        sun_altitude_deg: f64,
        extinction_coeff: f64,
    ) -> PyResult<f64> {
        Ok(py.detach(|| {
            schaefer::calc_heliacal_altitude_threshold(
                object_magnitude,
                sun_altitude_deg,
                extinction_coeff,
            )
        }))
    }
    /// Return a thirteen-key atmospheric/visibility diagnostic dictionary.
    ///
    /// Sun/object altitudes are degrees and object_magnitude is intrinsic.
    /// Defaults: pressure_mbar=1013.25 mbar/hPa, humidity_fraction=0.5 and
    /// altitude_m=0 metres. Uses fixed 15 C, range 0, no Moon, ecliptic latitude
    /// 0, elongation 90 degrees, observer age 25 and Snellen ratio 1.
    ///
    /// Keys 'extinction_total', 'extinction_rayleigh', 'extinction_aerosol'
    /// are mag/airmass; 'airmass' dimensionless; 'extinction_magnitude',
    /// 'limiting_magnitude', 'apparent_magnitude', 'magnitude_margin' magnitudes;
    /// 'sky_brightness_nl' nL; 'sky_brightness_mag_arcsec2' mag/arcsecond^2;
    /// 'arcus_visionis', 'optimal_sun_altitude' degrees; 'is_visible' bool.
    /// Margin is intrinsic limit minus intrinsic object magnitude. Unlike
    /// is_object_visible, the dict's comparison has NO horizon/daylight/
    /// elongation gate. Component ValueError/OverflowError failures propagate;
    /// physical ranges are not guarded and IEEE NaNs may remain in the dict.
    #[pyfunction(signature=(sun_altitude_deg, object_altitude_deg, object_magnitude, pressure_mbar=1013.25, humidity_fraction=0.5, altitude_m=0.0))]
    fn get_visibility_conditions<'py>(
        py: Python<'py>,
        sun_altitude_deg: f64,
        object_altitude_deg: f64,
        object_magnitude: f64,
        pressure_mbar: f64,
        humidity_fraction: f64,
        altitude_m: f64,
    ) -> PyResult<Bound<'py, PyDict>> {
        let value = py
            .detach(|| {
                schaefer::get_visibility_conditions(
                    sun_altitude_deg,
                    object_altitude_deg,
                    object_magnitude,
                    pressure_mbar,
                    humidity_fraction,
                    altitude_m,
                )
            })
            .map_err(schaefer_error)?;
        let output = PyDict::new(py);
        output.set_item("extinction_total", value.extinction_total)?;
        output.set_item("extinction_rayleigh", value.extinction_rayleigh)?;
        output.set_item("extinction_aerosol", value.extinction_aerosol)?;
        output.set_item("airmass", value.airmass)?;
        output.set_item("extinction_magnitude", value.extinction_magnitude)?;
        output.set_item("sky_brightness_nl", value.sky_brightness_nl)?;
        output.set_item(
            "sky_brightness_mag_arcsec2",
            value.sky_brightness_mag_arcsec2,
        )?;
        output.set_item("limiting_magnitude", value.limiting_magnitude)?;
        output.set_item("apparent_magnitude", value.apparent_magnitude)?;
        output.set_item("magnitude_margin", value.magnitude_margin)?;
        output.set_item("is_visible", value.is_visible)?;
        output.set_item("arcus_visionis", value.arcus_visionis)?;
        output.set_item("optimal_sun_altitude", value.optimal_sun_altitude)?;
        Ok(output)
    }
    pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
        module.add_function(wrap_pyfunction!(calc_rayleigh_extinction, module)?)?;
        module.add_function(wrap_pyfunction!(calc_aerosol_extinction, module)?)?;
        module.add_function(wrap_pyfunction!(calc_ozone_extinction, module)?)?;
        module.add_function(wrap_pyfunction!(calc_water_vapor_extinction, module)?)?;
        module.add_function(wrap_pyfunction!(calc_total_extinction, module)?)?;
        module.add_function(wrap_pyfunction!(calc_airmass, module)?)?;
        module.add_function(wrap_pyfunction!(calc_extinction_magnitude, module)?)?;
        module.add_function(wrap_pyfunction!(calc_twilight_brightness, module)?)?;
        module.add_function(wrap_pyfunction!(calc_moon_brightness, module)?)?;
        module.add_function(wrap_pyfunction!(calc_airglow_brightness, module)?)?;
        module.add_function(wrap_pyfunction!(calc_zodiacal_brightness, module)?)?;
        module.add_function(wrap_pyfunction!(calc_total_sky_brightness, module)?)?;
        module.add_function(wrap_pyfunction!(brightness_to_mag_arcsec2, module)?)?;
        module.add_function(wrap_pyfunction!(mag_arcsec2_to_brightness, module)?)?;
        module.add_function(wrap_pyfunction!(calc_limiting_magnitude, module)?)?;
        module.add_function(wrap_pyfunction!(is_object_visible, module)?)?;
        module.add_function(wrap_pyfunction!(get_arcus_visionis, module)?)?;
        module.add_function(wrap_pyfunction!(get_optimal_sun_altitude, module)?)?;
        module.add_function(wrap_pyfunction!(calc_heliacal_altitude_threshold, module)?)?;
        module.add_function(wrap_pyfunction!(get_visibility_conditions, module)?)?;
        module.add(
            "RAYLEIGH_COEFF_SEA_LEVEL",
            schaefer::RAYLEIGH_COEFF_SEA_LEVEL,
        )?;
        module.add("RAYLEIGH_SCALE_HEIGHT", schaefer::RAYLEIGH_SCALE_HEIGHT)?;
        module.add("AEROSOL_COEFF_CLEAR", schaefer::AEROSOL_COEFF_CLEAR)?;
        module.add("AEROSOL_COEFF_AVERAGE", schaefer::AEROSOL_COEFF_AVERAGE)?;
        module.add("AEROSOL_COEFF_POOR", schaefer::AEROSOL_COEFF_POOR)?;
        module.add("OZONE_COEFF", schaefer::OZONE_COEFF)?;
        module.add("WATER_VAPOR_COEFF", schaefer::WATER_VAPOR_COEFF)?;
        module.add("AIRGLOW_BRIGHTNESS_NL", schaefer::AIRGLOW_BRIGHTNESS_NL)?;
        module.add("AIRGLOW_MAG_ARCSEC2", schaefer::AIRGLOW_MAG_ARCSEC2)?;
        module.add("ZODIACAL_MIN_NL", schaefer::ZODIACAL_MIN_NL)?;
        module.add("ZODIACAL_MAX_NL", schaefer::ZODIACAL_MAX_NL)?;
        module.add(
            "FULL_MOON_SKY_BRIGHTNESS_NL",
            schaefer::FULL_MOON_SKY_BRIGHTNESS_NL,
        )?;
        module.add(
            "TWILIGHT_TRANSITION_ZENITH",
            schaefer::TWILIGHT_TRANSITION_ZENITH,
        )?;
        module.add("DARK_ADAPTED_PUPIL_MM", schaefer::DARK_ADAPTED_PUPIL_MM)?;
        module.add("LIGHT_ADAPTED_PUPIL_MM", schaefer::LIGHT_ADAPTED_PUPIL_MM)?;
        module.add(
            "IDEAL_NAKED_EYE_LIMIT_MAG",
            schaefer::IDEAL_NAKED_EYE_LIMIT_MAG,
        )?;
        module.add(
            "AGE_DEGRADATION_PER_DECADE",
            schaefer::AGE_DEGRADATION_PER_DECADE,
        )?;
        module.add("AGE_BASELINE", schaefer::AGE_BASELINE)?;
        module.add(
            "PTOLEMY_ARCUS_VISIONIS",
            schaefer::PTOLEMY_ARCUS_VISIONIS.to_vec(),
        )?;
        module.add("MIN_ELONGATION_BRIGHT", schaefer::MIN_ELONGATION_BRIGHT)?;
        module.add("MIN_ELONGATION_NORMAL", schaefer::MIN_ELONGATION_NORMAL)?;
        module.add("MIN_ELONGATION_FAINT", schaefer::MIN_ELONGATION_FAINT)?;
        Ok(())
    }
}

#[allow(non_snake_case)]
mod refraction_bindings {
    use super::*;
    /// Return APPARENT altitude, not a refraction increment, from true altitude.
    ///
    /// true_alt is geometric degrees; pressure defaults to 1013.25 mbar and
    /// temperature_C to 15 Celsius. Uses Saemundsson/upper-altitude closed
    /// forms. Nonfinite true_alt, values <= -5 degrees and candidates not
    /// strictly above zero return the input unchanged. Atmosphere values are
    /// not physically validated; IEEE values may propagate without errors.
    /// This independent approximation need not invert the reverse helper.
    #[pyfunction(signature=(true_alt, pressure=1013.25, temperature_C=15.0))]
    fn calc_refrac_compat_true_to_app(
        py: Python<'_>,
        true_alt: f64,
        pressure: f64,
        temperature_C: f64,
    ) -> PyResult<f64> {
        Ok(py.detach(|| {
            refraction::calc_refrac_compat_true_to_app(true_alt, pressure, temperature_C)
        }))
    }
    /// Return TRUE altitude, not a refraction increment, from apparent altitude.
    ///
    /// apparent_alt is degrees; pressure defaults to 1013.25 mbar and
    /// temperature_C to 15 Celsius. Uses Bennett's closed form, not numerical
    /// inversion of the forward fit. Nonfinite altitude, the exact -4.4 pole,
    /// guarded near-zenith geometry and nonpositive candidate altitudes return
    /// the input unchanged. No physical atmosphere guard or domain exception
    /// is added; IEEE arithmetic may propagate. Exact round trips are not promised.
    #[pyfunction(signature=(apparent_alt, pressure=1013.25, temperature_C=15.0))]
    fn calc_refrac_compat_app_to_true(
        py: Python<'_>,
        apparent_alt: f64,
        pressure: f64,
        temperature_C: f64,
    ) -> PyResult<f64> {
        Ok(py.detach(|| {
            refraction::calc_refrac_compat_app_to_true(apparent_alt, pressure, temperature_C)
        }))
    }
    /// Return ICAO ray-traced refraction CORRECTION in degrees for true altitude.
    ///
    /// true_alt is geometric altitude in degrees; add the result to obtain
    /// apparent altitude. pressure defaults to 1013.25 mbar, temperature_C to
    /// 15 Celsius, obs_alt to 0 metres and lapse_rate to 0.0065 K/m (positive
    /// temperature decrease with height). Uses a sampled atmosphere and
    /// bounded quadrature/fixed-point iteration, with below-horizon extension.
    /// Nonpositive pressure or nonphysical absolute temperature gives zero;
    /// no implicit pressure estimate is made. Other unphysical inputs follow
    /// IEEE/branch semantics rather than new errors; no ephemeris is read.
    #[pyfunction(signature=(true_alt, pressure=1013.25, temperature_C=15.0, obs_alt=0.0, lapse_rate=0.0065))]
    fn calc_refraction_true_to_app(
        py: Python<'_>,
        true_alt: f64,
        pressure: f64,
        temperature_C: f64,
        obs_alt: f64,
        lapse_rate: f64,
    ) -> PyResult<f64> {
        Ok(py.detach(|| {
            refraction::calc_refraction_true_to_app(
                true_alt,
                pressure,
                temperature_C,
                obs_alt,
                lapse_rate,
            )
        }))
    }
    /// Return ICAO ray-traced refraction CORRECTION for apparent altitude.
    ///
    /// apparent_alt is degrees; subtract the returned nonnegative degree
    /// correction for estimated true altitude. Defaults/units: pressure 1013.25
    /// mbar, temperature_C 15 Celsius, obs_alt 0 metres, lapse_rate 0.0065 K/m.
    /// Newton inversion uses at most eight steps and can stop at its derivative
    /// floor; convergence failure is not raised. Nonpositive pressure returns
    /// zero, not an estimated pressure. Other scalar inputs follow native
    /// IEEE/branch rules without physical validation; no clock/data is read.
    #[pyfunction(signature=(apparent_alt, pressure=1013.25, temperature_C=15.0, obs_alt=0.0, lapse_rate=0.0065))]
    fn calc_refraction_app_to_true(
        py: Python<'_>,
        apparent_alt: f64,
        pressure: f64,
        temperature_C: f64,
        obs_alt: f64,
        lapse_rate: f64,
    ) -> PyResult<f64> {
        Ok(py.detach(|| {
            refraction::calc_refraction_app_to_true(
                apparent_alt,
                pressure,
                temperature_C,
                obs_alt,
                lapse_rate,
            )
        }))
    }
    /// Return the same degree correction as calc_refraction_app_to_true.
    ///
    /// apparent_alt is degrees; pressure=1013.25 mbar, temperature_C=15 C,
    /// obs_alt=0 metres and lapse_rate=0.0065 K/m. Subtract the result from
    /// apparent altitude. This retained spelling uses the identical bounded
    /// ICAO inversion and IEEE/nonphysical-input behavior, not a different
    /// reference model. No convergence exception or data acquisition is added.
    #[pyfunction(signature=(apparent_alt, pressure=1013.25, temperature_C=15.0, obs_alt=0.0, lapse_rate=0.0065))]
    fn calc_refraction_ref_app_to_true(
        py: Python<'_>,
        apparent_alt: f64,
        pressure: f64,
        temperature_C: f64,
        obs_alt: f64,
        lapse_rate: f64,
    ) -> PyResult<f64> {
        Ok(py.detach(|| {
            refraction::calc_refraction_ref_app_to_true(
                apparent_alt,
                pressure,
                temperature_C,
                obs_alt,
                lapse_rate,
            )
        }))
    }
    /// Return the same degree correction as calc_refraction_true_to_app.
    ///
    /// true_alt is degrees; pressure=1013.25 mbar, temperature_C=15 C,
    /// obs_alt=0 metres and lapse_rate=0.0065 K/m. Add the correction to obtain
    /// apparent altitude. This retained spelling uses identical bounded ICAO
    /// ray tracing and IEEE/nonphysical-input behavior, not a different
    /// reference model or a converted-altitude return. No data is acquired.
    #[pyfunction(signature=(true_alt, pressure=1013.25, temperature_C=15.0, obs_alt=0.0, lapse_rate=0.0065))]
    fn calc_refraction_ref_true_to_app(
        py: Python<'_>,
        true_alt: f64,
        pressure: f64,
        temperature_C: f64,
        obs_alt: f64,
        lapse_rate: f64,
    ) -> PyResult<f64> {
        Ok(py.detach(|| {
            refraction::calc_refraction_ref_true_to_app(
                true_alt,
                pressure,
                temperature_C,
                obs_alt,
                lapse_rate,
            )
        }))
    }
    /// Return Thom refraction-corrected horizon dip in degrees (normally negative).
    ///
    /// obs_alt is height in metres; lapse_rate defaults to 0.0065 K/m,
    /// atpress to 1013.25 mbar and attemp to 10 Celsius. Height <= 0 or
    /// nonpositive model absolute temperature returns zero; the refractive
    /// square-root argument is floored at zero. Other inputs are not physically
    /// validated and IEEE values may propagate. No clock/data is acquired.
    #[pyfunction(signature=(obs_alt, lapse_rate=0.0065, atpress=1013.25, attemp=10.0))]
    fn calc_dip(
        py: Python<'_>,
        obs_alt: f64,
        lapse_rate: f64,
        atpress: f64,
        attemp: f64,
    ) -> PyResult<f64> {
        Ok(py.detach(|| refraction::calc_dip(obs_alt, lapse_rate, atpress, attemp)))
    }
    pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
        module.add_function(wrap_pyfunction!(calc_refrac_compat_true_to_app, module)?)?;
        module.add_function(wrap_pyfunction!(calc_refrac_compat_app_to_true, module)?)?;
        module.add_function(wrap_pyfunction!(calc_refraction_true_to_app, module)?)?;
        module.add_function(wrap_pyfunction!(calc_refraction_app_to_true, module)?)?;
        module.add_function(wrap_pyfunction!(calc_refraction_ref_app_to_true, module)?)?;
        module.add_function(wrap_pyfunction!(calc_refraction_ref_true_to_app, module)?)?;
        module.add_function(wrap_pyfunction!(calc_dip, module)?)?;
        Ok(())
    }
}

#[allow(non_snake_case)]
mod house_constructions_bindings {
    use super::*;
    /// Construct Savard-A cusps from caller-supplied, consistent degree angles.
    ///
    /// armc is active MC right ascension, lat north geographic latitude,
    /// eps true obliquity, asc Ascendant longitude and mc MC longitude.
    /// Returns a Python LIST of 13 degree values: index 0 is a LEADING ZERO,
    /// indices 1..12 are houses 1..12 (unlike the high-level houses tuples).
    /// Takes the supplied cardinals as-is; no clock, ephemeris, polar fallback
    /// or physical/finite-input validation is added. IEEE values can propagate.
    #[pyfunction(signature=(armc, lat, eps, asc, mc))]
    fn houses_savard_a(
        py: Python<'_>,
        armc: f64,
        lat: f64,
        eps: f64,
        asc: f64,
        mc: f64,
    ) -> PyResult<Vec<f64>> {
        let cusps = py.detach(|| houses::houses_savard_a(armc, lat, eps, asc, mc));
        Ok(std::iter::once(0.0).chain(cusps).collect())
    }
    /// Construct Krusinski-Pisa cusps from caller-supplied degree geometry.
    ///
    /// armc is active MC right ascension, lat north geographic latitude,
    /// eps true obliquity and asc Ascendant longitude. mc is accepted but
    /// UNUSED: cusp 10 comes from the zenith construction. Returns a Python
    /// LIST [0.0, cusp1, ..., cusp12] in degrees, not an unpadded cusp tuple.
    /// A degenerate Ascendant-zenith circle raises the native calculation
    /// exception. No clock/data, physical-range checks or fallback are added;
    /// other invalid geometry can retain IEEE NaN/infinity behavior.
    #[pyfunction(signature=(armc, lat, eps, asc, mc))]
    fn houses_krusinski(
        py: Python<'_>,
        armc: f64,
        lat: f64,
        eps: f64,
        asc: f64,
        mc: f64,
    ) -> PyResult<Vec<f64>> {
        let cusps = py
            .detach(|| houses::houses_krusinski(armc, lat, eps, asc, mc))
            .map_err(|error| to_pyerr(py, error))?;
        Ok(std::iter::once(0.0).chain(cusps).collect())
    }
    /// Construct Sunshine/Makransky cusps from explicit geometry in degrees.
    ///
    /// armc is active MC right ascension, lat north geographic latitude,
    /// eps true obliquity, asc/mc consistent ecliptic cardinals, sun_dec solar
    /// declination. Returns a Python LIST [0.0, cusp1, ..., cusp12] in degrees.
    /// A circumpolar Sun (abs(tan(sun_dec)*tan(lat)) >= 1) raises
    /// PolarCircleError; no alternative house system or solar declination is
    /// acquired. This direct kernel takes supplied cardinals as-is and adds
    /// no general finite/range admission; IEEE values may propagate.
    #[pyfunction(signature=(armc, lat, eps, asc, mc, sun_dec))]
    fn houses_sunshine_makransky(
        py: Python<'_>,
        armc: f64,
        lat: f64,
        eps: f64,
        asc: f64,
        mc: f64,
        sun_dec: f64,
    ) -> PyResult<Vec<f64>> {
        let cusps = py
            .detach(|| houses::houses_sunshine_makransky(armc, lat, eps, asc, mc, sun_dec))
            .map_err(|error| to_pyerr(py, error))?;
        Ok(std::iter::once(0.0).chain(cusps).collect())
    }
    /// Return one APC cusp longitude in DEGREES from RADIAN geometry.
    ///
    /// index is a one-based house number 1..12; lat_rad is geographic north
    /// latitude, eps_rad true obliquity and armc_rad MC right ascension, all
    /// radians. An index outside 1..12 raises KeyError. Other geometric inputs
    /// are not finite/range-checked: degenerate geometry retains its native
    /// floating-point behavior. No date, ephemeris or observer is consulted.
    #[pyfunction(signature=(index, lat_rad, eps_rad, armc_rad))]
    fn apc_cusp(
        py: Python<'_>,
        index: i64,
        lat_rad: f64,
        eps_rad: f64,
        armc_rad: f64,
    ) -> PyResult<f64> {
        py.detach(|| houses::apc_cusp_checked(index, lat_rad, eps_rad, armc_rad))
            .map_err(|error| to_pyerr(py, error))
    }
    pub(super) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
        module.add_function(wrap_pyfunction!(houses_savard_a, module)?)?;
        module.add_function(wrap_pyfunction!(houses_krusinski, module)?)?;
        module.add_function(wrap_pyfunction!(houses_sunshine_makransky, module)?)?;
        module.add_function(wrap_pyfunction!(apc_cusp, module)?)?;
        Ok(())
    }
}

pub(crate) fn register(module: &Bound<'_, PyModule>) -> PyResult<()> {
    schaefer_bindings::register(&submodule(module, "schaefer")?)?;
    refraction_bindings::register(&submodule(module, "refraction")?)?;
    house_constructions_bindings::register(&submodule(module, "house_constructions")?)?;
    Ok(())
}
