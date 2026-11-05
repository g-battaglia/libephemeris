//! Typed public records and methods over the core VISLIMIT optical model.
use crate::errors::to_pyerr;
use libephemeris_core::heliacal::{HeliacalObserver, VislimAtmosphere, VislimModel};
use libephemeris_core::schaefer;
use pyo3::prelude::*;
use pyo3::types::PyDict;

/// Immutable atmosphere inputs for SchaeferModel, not an ephemeris context.
///
/// Constructor values are stored verbatim. Model construction checks consumed
/// values for finiteness, but does not certify physical ranges.
#[pyclass(
    module = "libephemeris.heliacal",
    frozen,
    get_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, Copy, PartialEq)]
pub(crate) struct AtmosphericConditions {
    /// Pressure in mbar (hPa); zero selects the altitude-form molecular term.
    pressure: f64,
    /// Air temperature in degrees Celsius.
    temperature: f64,
    /// Relative humidity in percent, normally 0..100, not a fraction.
    humidity: f64,
    /// Zero: default; 0 < value < 1: total mag/airmass; >=1: visibility km.
    met_range: f64,
    /// Observer height in metres above sea level.
    altitude: f64,
}
#[pymethods]
impl AtmosphericConditions {
    /// Store pressure (mbar), temperature (Celsius), humidity (percent),
    /// met_range (km or direct extinction coefficient) and altitude (metres).
    ///
    /// No defaults or validation are applied by this record constructor.
    #[new]
    fn new(pressure: f64, temperature: f64, humidity: f64, met_range: f64, altitude: f64) -> Self {
        Self {
            pressure,
            temperature,
            humidity,
            met_range,
            altitude,
        }
    }
}
impl Default for AtmosphericConditions {
    fn default() -> Self {
        Self::new(1013.25, 15.0, 50.0, 0.0, 0.0)
    }
}
impl AtmosphericConditions {
    fn core(&self, latitude: f64, jd: f64) -> libephemeris_core::errors::Result<VislimAtmosphere> {
        VislimAtmosphere::new(
            self.pressure,
            self.temperature,
            self.humidity,
            self.met_range,
            self.altitude,
            latitude,
            jd,
        )
    }
}

/// Immutable observer and instrument inputs for the VISLIMIT optical model.
///
/// Values are stored without validation; SchaeferModel checks their finiteness.
/// The model applies instrument corrections only for magnification >1 and
/// positive aperture. Binocular viewing modifies that instrument correction.
#[pyclass(
    module = "libephemeris.heliacal",
    frozen,
    get_all,
    eq,
    skip_from_py_object
)]
#[derive(Clone, Copy, PartialEq)]
pub(crate) struct ObserverParams {
    /// Observer age in years; 36 is the model's naked-eye reference.
    age: f64,
    /// Dimensionless Snellen acuity ratio; 1 means normal acuity.
    snellen: f64,
    /// Whether instrument viewing uses both eyes.
    binocular: bool,
    /// Dimensionless telescope magnification; zero denotes naked-eye viewing.
    telescope_mag: f64,
    /// Instrument objective diameter in millimetres.
    aperture: f64,
    /// Optical transmission fraction, normally 0..1 (not percent).
    transmission: f64,
}
#[pymethods]
impl ObserverParams {
    /// Store age (years), Snellen ratio, binocular bool, magnification,
    /// aperture diameter (mm) and transmission fraction without validation.
    #[new]
    fn new(
        age: f64,
        snellen: f64,
        binocular: bool,
        telescope_mag: f64,
        aperture: f64,
        transmission: f64,
    ) -> Self {
        Self {
            age,
            snellen,
            binocular,
            telescope_mag,
            aperture,
            transmission,
        }
    }
}
impl Default for ObserverParams {
    fn default() -> Self {
        Self::new(36.0, 1.0, false, 0.0, 0.0, 1.0)
    }
}
impl ObserverParams {
    fn core(&self) -> libephemeris_core::errors::Result<HeliacalObserver> {
        HeliacalObserver::new(
            self.age,
            self.snellen,
            self.binocular,
            self.telescope_mag,
            self.aperture,
            self.transmission,
        )
    }
}

/// Reusable V-band VISLIMIT atmospheric, sky-radiance and detection model.
///
/// Geometry is supplied by the caller in degrees; this object does not acquire
/// ephemeris data. Limiting magnitudes include extinction and compare directly
/// with extra-atmospheric catalog magnitudes. This approximate visual model is
/// not a guarantee that an object or heliacal event will be observed.
///
/// Examples
/// --------
/// >>> from libephemeris.heliacal import SchaeferModel
/// >>> model = SchaeferModel()
/// >>> limit = model.limiting_magnitude(sun_alt=-12.0, obj_alt=15.0)
/// >>> isinstance(limit, float)
/// True
#[pyclass(module = "libephemeris.heliacal", skip_from_py_object)]
#[derive(Clone)]
pub(crate) struct SchaeferModel {
    native: VislimModel,
    atmo: AtmosphericConditions,
    observer: ObserverParams,
    latitude: f64,
}
impl SchaeferModel {
    fn build(
        py: Python<'_>,
        atmo: AtmosphericConditions,
        observer: ObserverParams,
        latitude: f64,
        jd: f64,
        pressure_scales_extinction: bool,
    ) -> PyResult<Self> {
        let native = py
            .detach(|| {
                Ok::<_, libephemeris_core::errors::Error>(VislimModel::new(
                    &atmo.core(latitude, jd)?,
                    observer.core()?,
                    pressure_scales_extinction,
                ))
            })
            .map_err(|error| to_pyerr(py, error))?;
        Ok(Self {
            native,
            atmo,
            observer,
            latitude,
        })
    }
}
#[pymethods]
impl SchaeferModel {
    /// Construct the optical model with immutable atmosphere/observer inputs.
    ///
    /// Parameters
    /// ----------
    /// atmo : AtmosphericConditions or None
    ///     None uses 1013.25 mbar, 15 C, 50% humidity, default range, sea level.
    /// observer : ObserverParams or None
    ///     None uses age 36, Snellen 1, monocular, no instrument, transmission 1.
    /// latitude : float
    ///     Geographic latitude in degrees north, normally -90..90.
    /// jd : float
    ///     UT1 Julian day used for season/solar-cycle state, default J2000.
    /// pressure_scales_extinction : bool
    ///     True uses positive pressure for molecular extinction; False uses
    ///     the standard-atmosphere altitude form regardless of stored pressure.
    ///
    /// Raises
    /// ------
    /// ValueError
    ///     A consumed atmosphere, observer, latitude or epoch value is nonfinite.
    ///     Finite physical ranges are not checked.
    #[new]
    #[pyo3(signature=(atmo=None, observer=None, latitude=0.0, jd=2451545.0, pressure_scales_extinction=false))]
    fn new(
        py: Python<'_>,
        atmo: Option<PyRef<'_, AtmosphericConditions>>,
        observer: Option<PyRef<'_, ObserverParams>>,
        latitude: f64,
        jd: f64,
        pressure_scales_extinction: bool,
    ) -> PyResult<Self> {
        Self::build(
            py,
            atmo.as_deref().copied().unwrap_or_default(),
            observer.as_deref().copied().unwrap_or_default(),
            latitude,
            jd,
            pressure_scales_extinction,
        )
    }
    /// Original immutable atmosphere record; seasonal updates do not alter it.
    #[getter]
    fn atmo(&self) -> AtmosphericConditions {
        self.atmo
    }
    /// Immutable observer/optics parameters used by this model.
    #[getter]
    fn observer(&self) -> ObserverParams {
        self.observer
    }
    /// Geographic latitude, degrees north.
    #[getter]
    fn latitude(&self) -> f64 {
        self.latitude
    }
    /// Current season-source UT1 Julian day, including update_season changes.
    #[getter]
    fn jd(&self) -> f64 {
        self.native.jd_ut()
    }
    /// Current molecular extinction coefficient, magnitudes per airmass.
    #[getter]
    fn k_rayleigh(&self) -> f64 {
        self.native.coeffs().k_rayleigh()
    }
    /// Current aerosol extinction coefficient, magnitudes per airmass.
    #[getter]
    fn k_aerosol(&self) -> f64 {
        self.native.coeffs().k_aerosol()
    }
    /// Current seasonal ozone coefficient, magnitudes per airmass.
    #[getter]
    fn k_ozone(&self) -> f64 {
        self.native.coeffs().k_ozone()
    }
    /// Current water-vapour extinction coefficient, magnitudes per airmass.
    #[getter]
    fn k_water(&self) -> f64 {
        self.native.coeffs().k_water()
    }
    /// Sum of current extinction components, magnitudes per airmass.
    #[getter]
    fn k_total(&self) -> f64 {
        self.native.coeffs().k_total()
    }
    /// Refresh seasonal coefficients and sky time at jd, a UT1 Julian day.
    ///
    /// Returns None. Repeating the current epoch is a no-op; a nonfinite epoch
    /// leaves the previous model unchanged. Original atmo/observer stay fixed.
    fn update_season(&mut self, py: Python<'_>, jd: f64) {
        py.detach(|| self.native.update_season(jd));
    }
    /// Return Rozenberg relative sky airmass for altitude_deg in degrees.
    ///
    /// Returns dimensionless path length, or 40 at/below the horizon. Intended
    /// for physical horizontal altitudes; no range/finiteness validation occurs.
    fn airmass(&self, py: Python<'_>, altitude_deg: f64) -> f64 {
        py.detach(|| libephemeris_core::heliacal::vl_airmass(altitude_deg))
    }
    /// Return dimming in magnitudes at altitude_deg above the horizon.
    ///
    /// Sums component-specific gas, aerosol and ozone path lengths, not merely
    /// k_total*airmass. Intended above the horizon; no horizon sentinel or
    /// physical range validation is applied by this method.
    fn extinction(&self, py: Python<'_>, altitude_deg: f64) -> f64 {
        py.detach(|| self.native.coeffs().extinction_at_altitude(altitude_deg))
    }
    /// Return total sky brightness in nanolamberts, despite the '_bl' name.
    ///
    /// Parameters
    /// ----------
    /// sun_alt, moon_alt, obj_alt : float
    ///     Sun, Moon and object altitude above the horizon, degrees.
    /// moon_phase : float
    ///     Illuminated fraction, 0 new to 1 full (not phase angle or percent).
    /// sun_obj_angle, moon_obj_angle : float
    ///     Angular separations from the object, degrees, normally 0..180.
    ///
    /// Notes
    /// -----
    /// Combines dark sky, daylight or twilight and risen-Moon scattered light.
    /// Angles are caller-supplied, with no ephemeris or validation. Result is
    /// not lamberts or mag/arcsec²; nonphysical inputs may yield invalid values.
    #[allow(clippy::too_many_arguments)]
    fn sky_brightness_bl(
        &self,
        py: Python<'_>,
        sun_alt: f64,
        moon_alt: f64,
        moon_phase: f64,
        obj_alt: f64,
        sun_obj_angle: f64,
        moon_obj_angle: f64,
    ) -> f64 {
        py.detach(|| {
            self.native.sky_brightness_nl(
                sun_alt,
                moon_alt,
                moon_phase,
                obj_alt,
                sun_obj_angle,
                moon_obj_angle,
            )
        })
    }
    /// Return the limiting extra-atmospheric catalog V magnitude.
    ///
    /// sun_alt, moon_alt and obj_alt are horizontal altitudes in degrees;
    /// sun_obj_angle/moon_obj_angle are separations in degrees. moon_phase is
    /// illuminated fraction (0 new, 1 full). Default Moon altitude -90 suppresses
    /// moonlight. Extinction, acuity, age or instrument corrections are already
    /// included: do not add extinction again when comparing catalog magnitude.
    /// Geometry is not validated and this method does not reject daylight or
    /// below-horizon objects; use is_visible for those gates.
    #[pyo3(signature=(sun_alt, moon_alt=-90.0, moon_phase=0.0, obj_alt=45.0, sun_obj_angle=180.0, moon_obj_angle=180.0),text_signature="($self, sun_alt, moon_alt=-90.0, moon_phase=0.0, obj_alt=45.0, sun_obj_angle=180.0, moon_obj_angle=180.0)")]
    #[allow(clippy::too_many_arguments)]
    fn limiting_magnitude(
        &self,
        py: Python<'_>,
        sun_alt: f64,
        moon_alt: f64,
        moon_phase: f64,
        obj_alt: f64,
        sun_obj_angle: f64,
        moon_obj_angle: f64,
    ) -> f64 {
        py.detach(|| {
            self.native.limiting_magnitude(
                sun_alt,
                moon_alt,
                moon_phase,
                obj_alt,
                sun_obj_angle,
                moon_obj_angle,
            )
        })
    }
    /// Estimate the minimum vertical arc required for catalog V magnitude.
    ///
    /// Returns degrees by scanning solar depressions -1..-22 degrees and
    /// bounded altitude bisection with the object above the Sun. Returns 0.0
    /// if no sampled geometry is visible; zero is not proof of zero required
    /// separation. body_type is accepted but ignored. No convergence or
    /// observational-visibility guarantee is provided.
    #[pyo3(signature=(body_mag, body_type="planet"))]
    fn arcus_visionis_required(&self, py: Python<'_>, body_mag: f64, body_type: &str) -> f64 {
        let _ = body_type;
        py.detach(|| self.native.arcus_visionis_required(body_mag))
    }
    /// Test supplied twilight geometry against the catalog-magnitude limit.
    ///
    /// body_alt, sun_alt, elongation, moon_alt and moon_obj_angle are degrees;
    /// body_mag is catalog V magnitude and moon_phase illuminated fraction
    /// (0 new, 1 full). margin is magnitudes (default 0.5; larger is stricter).
    /// Returns False for body altitude <0, Sun altitude >0 or elongation <7;
    /// otherwise tests body_mag <= limiting_magnitude - margin. No ephemeris
    /// is acquired and numeric physical ranges are not validated.
    #[allow(clippy::too_many_arguments)]
    #[pyo3(signature=(body_alt, body_mag, sun_alt, elongation, moon_alt=-90.0, moon_phase=0.0, moon_obj_angle=180.0, margin=0.5),text_signature="($self, body_alt, body_mag, sun_alt, elongation, moon_alt=-90.0, moon_phase=0.0, moon_obj_angle=180.0, margin=0.5)")]
    fn is_visible(
        &self,
        py: Python<'_>,
        body_alt: f64,
        body_mag: f64,
        sun_alt: f64,
        elongation: f64,
        moon_alt: f64,
        moon_phase: f64,
        moon_obj_angle: f64,
        margin: f64,
    ) -> bool {
        py.detach(|| {
            self.native.is_visible(
                body_alt,
                body_mag,
                sun_alt,
                elongation,
                moon_alt,
                moon_phase,
                moon_obj_angle,
                margin,
            )
        })
    }
    /// Return max(sun_alt + required_arcus_visionis(body_mag), 0), in degrees.
    ///
    /// body_mag is catalog V magnitude; sun_alt is solar altitude in degrees.
    /// This convenience altitude estimate does not search for an event. The
    /// underlying required-arc calculation uses zero for no sampled solution.
    #[pyo3(signature=(body_mag, sun_alt=-8.0),text_signature="($self, body_mag, sun_alt=-8.0)")]
    fn heliacal_visibility_angle(&self, py: Python<'_>, body_mag: f64, sun_alt: f64) -> f64 {
        py.detach(|| self.native.heliacal_visibility_angle(body_mag, sun_alt))
    }
}

/// Reference V-band atmosphere, eye and sky constants for model users.
///
/// These public reference values are not a complete specification of the
/// VISLIMIT equations; the model does not consume every class attribute.
/// K_RAYLEIGH_SEA_LEVEL/K_OZONE/K_AEROSOL_BASE are mag/airmass;
/// SCALE_HEIGHT_PRESSURE/SCALE_HEIGHT_AEROSOL are km; AIRGLOW/ZODIACAL_LIGHT
/// are nanolamberts; DARK_SKY_MAG is mag/arcsec²; FULL_MOON_MAG_DIFF is a
/// magnitude difference; PERFECT_SKY_LIM_MAG is stellar magnitude;
/// THRESHOLD_CONTRAST is a fraction; DARK_PUPIL_DIAMETER/BRIGHT_PUPIL_DIAMETER
/// are mm; ARCUS_VISIONIS maps body/category names to reference degrees.
#[pyclass(module = "libephemeris.heliacal")]
pub(crate) struct SchaeferConstants;
#[pymethods]
impl SchaeferConstants {
    /// Construct the stateless constants namespace; no ephemeris is acquired.
    #[new]
    fn new() -> Self {
        Self
    }
    /// Reference sea-level molecular extinction, magnitudes per airmass.
    #[classattr]
    const K_RAYLEIGH_SEA_LEVEL: f64 = schaefer::RAYLEIGH_COEFF_SEA_LEVEL;
    /// Reference ozone extinction, magnitudes per airmass.
    #[classattr]
    const K_OZONE: f64 = schaefer::OZONE_COEFF;
    /// Reference aerosol extinction, magnitudes per airmass.
    #[classattr]
    const K_AEROSOL_BASE: f64 = 0.08;
    /// Reference atmospheric pressure scale height, kilometres.
    #[classattr]
    const SCALE_HEIGHT_PRESSURE: f64 = 8.5;
    /// Reference aerosol scale height, kilometres.
    #[classattr]
    const SCALE_HEIGHT_AEROSOL: f64 = 1.5;
    /// Reference airglow sky brightness, nanolamberts.
    #[classattr]
    const AIRGLOW: f64 = schaefer::AIRGLOW_BRIGHTNESS_NL;
    /// Reference zodiacal-light brightness, nanolamberts.
    #[classattr]
    const ZODIACAL_LIGHT: f64 = 100.0;
    /// Reference dark-sky V surface brightness, mag/arcsec².
    #[classattr]
    const DARK_SKY_MAG: f64 = 21.6;
    /// Reference full-Moon versus dark-sky brightness difference, magnitudes.
    #[classattr]
    const FULL_MOON_MAG_DIFF: f64 = 14.0;
    /// Ideal reference naked-eye limiting stellar magnitude.
    #[classattr]
    const PERFECT_SKY_LIM_MAG: f64 = schaefer::IDEAL_NAKED_EYE_LIMIT_MAG;
    /// Reference dimensionless threshold contrast fraction (not percent).
    #[classattr]
    const THRESHOLD_CONTRAST: f64 = 0.0094;
    /// Reference dark-adapted pupil diameter, millimetres.
    #[classattr]
    const DARK_PUPIL_DIAMETER: f64 = schaefer::DARK_ADAPTED_PUPIL_MM;
    /// Reference bright-adapted pupil diameter, millimetres.
    #[classattr]
    const BRIGHT_PUPIL_DIAMETER: f64 = 2.0;
    /// Reference body/category-to-vertical-arc dictionary, values in degrees.
    ///
    /// These nominal values do not replace the model's magnitude-based solve.
    #[classattr]
    #[allow(non_snake_case)]
    fn ARCUS_VISIONIS(py: Python<'_>) -> PyResult<Bound<'_, PyDict>> {
        let result = PyDict::new(py);
        for (key, value) in [
            ("venus", 5.0),
            ("mercury", 10.0),
            ("mars", 11.0),
            ("jupiter", 9.0),
            ("saturn", 11.0),
            ("star_bright", 7.0),
            ("star_medium", 10.0),
            ("star_faint", 13.0),
        ] {
            result.set_item(key, value)?;
        }
        Ok(result)
    }
}

/// Construct SchaeferModel from scalar atmospheric and observer parameters.
///
/// Parameters
/// ----------
/// pressure, temperature, humidity : float
///     Pressure mbar, temperature Celsius, relative humidity percent (not fraction).
/// met_range : float
///     Zero for default, 0<value<1 for total mag/airmass, >=1 for visibility km.
/// altitude, latitude, jd : float
///     Observer height metres, geographic latitude degrees north, UT1 Julian day.
/// observer_age, snellen : float
///     Age years and dimensionless visual acuity ratio (1 normal).
/// binocular, telescope_mag, aperture, transmission : bool, float, float, float
///     Two-eye instrument flag, magnification, objective diameter mm and
///     transmission fraction (normally 0..1). No instrument is the default.
/// pressure_scales_extinction : bool
///     False uses altitude for molecular extinction; True uses positive pressure.
///
/// Returns
/// -------
/// SchaeferModel
///     Reusable optical model; no ephemeris download or acquisition is performed.
///
/// Raises
/// ------
/// ValueError
///     Nonfinite consumed inputs. Finite physical ranges are not validated.
#[pyfunction(signature=(pressure=1013.25, temperature=15.0, humidity=50.0, met_range=0.0, altitude=0.0, observer_age=36.0, snellen=1.0, binocular=false, telescope_mag=0.0, aperture=0.0, transmission=1.0, latitude=0.0, jd=2451545.0, pressure_scales_extinction=false))]
#[allow(clippy::too_many_arguments)]
fn create_schaefer_model(
    py: Python<'_>,
    pressure: f64,
    temperature: f64,
    humidity: f64,
    met_range: f64,
    altitude: f64,
    observer_age: f64,
    snellen: f64,
    binocular: bool,
    telescope_mag: f64,
    aperture: f64,
    transmission: f64,
    latitude: f64,
    jd: f64,
    pressure_scales_extinction: bool,
) -> PyResult<SchaeferModel> {
    SchaeferModel::build(
        py,
        AtmosphericConditions::new(pressure, temperature, humidity, met_range, altitude),
        ObserverParams::new(
            observer_age,
            snellen,
            binocular,
            telescope_mag,
            aperture,
            transmission,
        ),
        latitude,
        jd,
        pressure_scales_extinction,
    )
}

pub(crate) fn register(parent: &Bound<'_, PyModule>) -> PyResult<()> {
    let module = PyModule::new(parent.py(), "libephemeris._core._heliacal_model")?;
    module.add_class::<AtmosphericConditions>()?;
    module.add_class::<ObserverParams>()?;
    module.add_class::<SchaeferConstants>()?;
    module.add_class::<SchaeferModel>()?;
    module.add_function(wrap_pyfunction!(create_schaefer_model, &module)?)?;
    parent
        .py()
        .import("sys")?
        .getattr("modules")?
        .set_item("libephemeris._core._heliacal_model", &module)?;
    parent.add("_heliacal_model", &module)?;
    Ok(())
}
