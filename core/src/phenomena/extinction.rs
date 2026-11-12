//! Atmospheric extinction, optical airmass and limiting magnitude.
//!
//! The models combine Kasten & Young (1989) airmass, Hayes & Latham Rayleigh
//! scattering, Koschmieder visibility, Ångström aerosol extinction,
//! ozone/water-vapour terms, twilight interpolation and Pogson flux ratios.
//! No epoch, frame transformation or ephemeris acquisition occurs here.
//!
//! Angles are degrees, pressure is mbar, temperature is Celsius, humidity
//! is percent, height is metres, wavelength is nanometres, and extinction
//! is magnitudes or magnitudes per airmass as identified by each entry.
//! Objects at or below the horizon use the documented extinction sentinel
//! and practical airmass cap.
//!
//! [`ObserverConditions`] validates finite atmospheric inputs.
//! [`ExtinctionCoefficients`] owns the component sum. Scalar kernels use
//! ordinary binary64 arithmetic; zero or negative wavelengths lie outside
//! their physical domain. The heliacal search composes these optical models
//! with provider-owned positions in its separate module.

use crate::errors::{Error, Result};

// ---------------------------------------------------------------------------
// Bands, defaults and scale heights.
// ---------------------------------------------------------------------------

/// U-band reference wavelength in nm.
pub const WAVELENGTH_U: f64 = 365.0;
/// B-band reference wavelength in nm.
pub const WAVELENGTH_B: f64 = 440.0;
/// V-band reference wavelength in nm.
pub const WAVELENGTH_V: f64 = 550.0;
/// R-band reference wavelength in nm.
pub const WAVELENGTH_R: f64 = 640.0;
/// I-band reference wavelength in nm.
pub const WAVELENGTH_I: f64 = 790.0;

/// Default sea-level pressure in mbar.
pub const DEFAULT_PRESSURE_MBAR: f64 = 1013.25;
/// Default air temperature in C.
pub const DEFAULT_TEMPERATURE_C: f64 = 15.0;
/// Default relative humidity in percent.
pub const DEFAULT_HUMIDITY_PERCENT: f64 = 50.0;
/// Default observer altitude in metres.
pub const DEFAULT_ALTITUDE_M: f64 = 0.0;

/// Atmospheric density scale height in metres.
pub const SCALE_HEIGHT_M: f64 = 8500.0;

/// Extinction sentinel for objects at or below the horizon, in magnitudes
pub const BELOW_HORIZON_EXTINCTION: f64 = 99.0;

/// Practical airmass cap for horizon objects.
pub const HORIZON_AIRMASS_CAP: f64 = 40.0;

/// Default simple-model extinction coefficient in mag/airmass.
pub const SIMPLE_EXTINCTION_DEFAULT_K: f64 = 0.28;

// Rayleigh.

/// V-band sea-level Rayleigh coefficient in mag/airmass, Hayes & Latham
/// (1975) via the release (`k_ref = 0.1451`).
pub const RAYLEIGH_REF_COEFF: f64 = 0.1451;
/// Rayleigh reference wavelength in nm.
pub const RAYLEIGH_REF_WAVELENGTH_NM: f64 = 550.0;

// Aerosol.

/// Aerosol scale height in metres.
pub const AEROSOL_SCALE_HEIGHT_M: f64 = 1500.0;
/// Clear-air aerosol baseline in mag/airmass.
pub const AEROSOL_BASE_COEFF: f64 = 0.08;
/// Koschmieder constant for the visibility arm.
pub const KOSCHMIEDER_CONSTANT: f64 = 3.912;
/// Molecular part removed from the Koschmieder estimate, mag/airmass
pub const KOSCHMIEDER_MOLECULAR_PART: f64 = 0.1066;

// Ozone pieces.

/// Ozone coefficient below 400 nm.
pub const OZONE_COEFF_UV: f64 = 0.10;
/// Ozone coefficient for 400-500 nm.
pub const OZONE_COEFF_BLUE: f64 = 0.04;
/// Ozone coefficient for 500-700 nm.
pub const OZONE_COEFF_VISUAL: f64 = 0.016;
/// Ozone coefficient at and above 700 nm.
pub const OZONE_COEFF_RED: f64 = 0.008;

// ---------------------------------------------------------------------------
// Twilight anchors.
// ---------------------------------------------------------------------------

/// Sun altitude of the day/civil boundary in degrees.
pub const TWILIGHT_CIVIL_START: f64 = 0.0;
/// Sun altitude of the civil/nautical boundary in degrees.
pub const TWILIGHT_CIVIL_END: f64 = -6.0;
/// Sun altitude of the nautical/astronomical boundary in degrees.
pub const TWILIGHT_NAUTICAL_END: f64 = -12.0;
/// Sun altitude of the astronomical/night boundary in degrees.
pub const TWILIGHT_ASTRONOMICAL_END: f64 = -18.0;

/// Typical dark-sky V-band surface brightness, mag/arcsec^2.
pub const DARK_SKY_BRIGHTNESS_V: f64 = 21.7;
/// Zenith dark-sky brightness, mag/arcsec^2.
pub const ZENITH_DARK_SKY: f64 = 21.9;

/// Twilight zenith anchor at Sun altitude 0, mag/arcsec^2.
pub const SKY_BRIGHTNESS_SUN_0: f64 = 3.0;
/// Twilight zenith anchor at Sun altitude -6, mag/arcsec^2.
pub const SKY_BRIGHTNESS_SUN_MINUS6: f64 = 8.0;
/// Twilight zenith anchor at Sun altitude -12, mag/arcsec^2.
pub const SKY_BRIGHTNESS_SUN_MINUS12: f64 = 17.0;
/// Twilight zenith anchor at Sun altitude -18, mag/arcsec^2.
pub const SKY_BRIGHTNESS_SUN_MINUS18: f64 = 21.7;

/// Photometric zero point of the mag/arcsec^2 to nanoLambert conversion
pub const NANOLAMBERT_ZERO_POINT: f64 = 26.33;

// ---------------------------------------------------------------------------
// Observer skill.
// ---------------------------------------------------------------------------

/// Inexperienced observer skill id.
pub const OBSERVER_SKILL_INEXPERIENCED: i32 = 1;
/// Average observer skill id.
pub const OBSERVER_SKILL_AVERAGE: i32 = 2;
/// Experienced observer skill id.
pub const OBSERVER_SKILL_EXPERIENCED: i32 = 3;
/// Expert observer skill id.
pub const OBSERVER_SKILL_EXPERT: i32 = 4;

// ---------------------------------------------------------------------------
// Tokens.
// ---------------------------------------------------------------------------

/// Relative-optical-path airmass model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AirmassMethod {
    /// Plane-parallel `sec(z)`, accurate only above roughly 15 deg altitude.
    Secant,
    /// Kasten & Young (1989) empirical form.
    KastenYoung,
    /// Rozenberg (1966) low-altitude form.
    Rozenberg,
}

impl AirmassMethod {
    /// Release method token spelled by this variant.
    pub fn as_str(self) -> &'static str {
        match self {
            AirmassMethod::Secant => "secant",
            AirmassMethod::KastenYoung => "kasten_young",
            AirmassMethod::Rozenberg => "rozenberg",
        }
    }
}

/// Exact, case-sensitive method token lookup: only the three registered tokens match.
/// `None` is the release `InputValidationError` arm; the Python exception
/// itself is raised on the binding side, like every release subclass
/// mapping in [`crate::errors`].
pub fn parse_airmass_method(method: &str) -> Option<AirmassMethod> {
    if method == "secant" {
        Some(AirmassMethod::Secant)
    } else if method == "kasten_young" {
        Some(AirmassMethod::KastenYoung)
    } else if method == "rozenberg" {
        Some(AirmassMethod::Rozenberg)
    } else {
        None
    }
}

/// Eye light-adaptation regime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EyeAdaptation {
    /// Cone vision, bright conditions.
    Photopic,
    /// Mixed rod/cone vision, twilight.
    Mesopic,
    /// Dark-adapted rod vision.
    Dark,
}

impl EyeAdaptation {
    /// Canonical release token spelled by this variant.
    pub fn as_str(self) -> &'static str {
        match self {
            EyeAdaptation::Photopic => "photopic",
            EyeAdaptation::Mesopic => "mesopic",
            EyeAdaptation::Dark => "dark",
        }
    }
}

/// Canonical eye-adaptation token lookup: exactly `"photopic"`, `"mesopic"`,
/// `"dark"`, plus the `"scotopic"` alias for dark. Matching is
/// case-sensitive; `None` is the release `InputValidationError` arm
/// (binding-side). `None` as a *missing* argument remains the
/// auto-detection sentinel handled by the callers, not by this function.
pub fn normalize_eye_adaptation(state: &str) -> Option<EyeAdaptation> {
    if state == "scotopic" {
        Some(EyeAdaptation::Dark)
    } else if state == "photopic" {
        Some(EyeAdaptation::Photopic)
    } else if state == "mesopic" {
        Some(EyeAdaptation::Mesopic)
    } else if state == "dark" {
        Some(EyeAdaptation::Dark)
    } else {
        None
    }
}

/// Twilight phase from the Sun altitude.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TwilightPhase {
    /// Sun above the horizon.
    Day,
    /// Sun 0 to -6 deg.
    Civil,
    /// Sun -6 to -12 deg.
    Nautical,
    /// Sun -12 to -18 deg.
    Astronomical,
    /// Sun at or below -18 deg.
    Night,
}

impl TwilightPhase {
    /// Release phase word spelled by this variant.
    pub fn as_str(self) -> &'static str {
        match self {
            TwilightPhase::Day => "day",
            TwilightPhase::Civil => "civil",
            TwilightPhase::Nautical => "nautical",
            TwilightPhase::Astronomical => "astronomical",
            TwilightPhase::Night => "night",
        }
    }
}

// ---------------------------------------------------------------------------
// Checked values.
// ---------------------------------------------------------------------------

/// Atmospheric extinction coefficients by component, mag/airmass. Fields are private: `k_total` is fixed at
/// construction as the component sum.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExtinctionCoefficients {
    k_rayleigh: f64,
    k_aerosol: f64,
    k_ozone: f64,
    k_water: f64,
    k_total: f64,
}

impl ExtinctionCoefficients {
    /// Checked constructor: records the four components and fixes
    /// `k_total` as their sum, the same association as the release.
    pub fn new(k_rayleigh: f64, k_aerosol: f64, k_ozone: f64, k_water: f64) -> Self {
        Self {
            k_rayleigh,
            k_aerosol,
            k_ozone,
            k_water,
            k_total: k_rayleigh + k_aerosol + k_ozone + k_water,
        }
    }

    /// Rayleigh (molecular) component, mag/airmass.
    pub fn k_rayleigh(&self) -> f64 {
        self.k_rayleigh
    }
    /// Aerosol (Mie) component, mag/airmass.
    pub fn k_aerosol(&self) -> f64 {
        self.k_aerosol
    }
    /// Ozone component, mag/airmass.
    pub fn k_ozone(&self) -> f64 {
        self.k_ozone
    }
    /// Water-vapour component, mag/airmass.
    pub fn k_water(&self) -> f64 {
        self.k_water
    }
    /// Total extinction coefficient, mag/airmass.
    pub fn k_total(&self) -> f64 {
        self.k_total
    }
}

/// Composite atmospheric inputs for the extinction composition (the
/// keyword bundle of `calc_extinction_coefficient` /
/// `calc_extinction_magnitude`). Fields are private with accessors; the
/// constructor refuses non-finite words fail-closed. Finite out-of-range
/// words (negative pressure, humidity above 100) are admitted exactly
/// like the release, which validates nothing beyond the token enums.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ObserverConditions {
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    altitude_m: f64,
    wavelength_nm: f64,
    visibility_km: Option<f64>,
}

impl ObserverConditions {
    /// Checked constructor. Units: pressure mbar, temperature C, humidity
    /// percent, observer altitude metres, wavelength nm, meteorological
    /// visibility km (`None` selects the humidity/altitude aerosol
    /// estimate). Non-finite numerics are refused; a zero wavelength is
    /// admitted and reports `+inf` coefficients (float division) where
    /// the release raises `ZeroDivisionError` — outside the swept domain.
    pub fn new(
        pressure_mbar: f64,
        temperature_c: f64,
        humidity_percent: f64,
        altitude_m: f64,
        wavelength_nm: f64,
        visibility_km: Option<f64>,
    ) -> Result<Self> {
        for value in [
            pressure_mbar,
            temperature_c,
            humidity_percent,
            altitude_m,
            wavelength_nm,
        ] {
            if !value.is_finite() {
                return Err(Error::NonFiniteArgument { value });
            }
        }
        if let Some(vis) = visibility_km {
            if !vis.is_finite() {
                return Err(Error::NonFiniteArgument { value: vis });
            }
        }
        Ok(Self {
            pressure_mbar,
            temperature_c,
            humidity_percent,
            altitude_m,
            wavelength_nm,
            visibility_km,
        })
    }

    /// Atmospheric pressure in mbar.
    pub fn pressure_mbar(&self) -> f64 {
        self.pressure_mbar
    }
    /// Air temperature in degrees Celsius.
    pub fn temperature_c(&self) -> f64 {
        self.temperature_c
    }
    /// Relative humidity in percent.
    pub fn humidity_percent(&self) -> f64 {
        self.humidity_percent
    }
    /// Observer altitude in metres.
    pub fn altitude_m(&self) -> f64 {
        self.altitude_m
    }
    /// Observation wavelength in nm.
    pub fn wavelength_nm(&self) -> f64 {
        self.wavelength_nm
    }
    /// Meteorological visibility in km, if supplied.
    pub fn visibility_km(&self) -> Option<f64> {
        self.visibility_km
    }

    /// Total extinction coefficient and its components under these
    /// conditions.
    pub fn coefficients(&self) -> ExtinctionCoefficients {
        calc_extinction_coefficient(
            self.pressure_mbar,
            self.temperature_c,
            self.humidity_percent,
            self.altitude_m,
            self.wavelength_nm,
            self.visibility_km,
        )
    }
}

/// Twilight sky-brightness outcome.
/// A plain datum: no validation, hence public fields.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TwilightSkyBrightness {
    /// Sky surface brightness, mag/arcsec^2 (higher = darker).
    pub surface_brightness: f64,
    /// Current twilight phase.
    pub twilight_phase: TwilightPhase,
    /// Input Sun altitude in degrees.
    pub sun_altitude_deg: f64,
    /// Angular-distance factor from the Sun (0 = toward the Sun, 1 = away).
    pub azimuth_factor: f64,
    /// Approximate naked-eye limiting magnitude.
    pub limiting_magnitude: f64,
    /// Sky brightness in nanoLamberts.
    pub nanolamberts: f64,
}

/// Visibility-threshold outcome. A plain
/// datum: no validation, hence public fields.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VisibilityResult {
    /// Whether the object clears the detection threshold.
    pub is_visible: bool,
    /// Effective (possibly extinction-corrected) object magnitude.
    pub object_magnitude: f64,
    /// Limiting magnitude at this sky brightness.
    pub limiting_magnitude: f64,
    /// Input sky surface brightness, mag/arcsec^2.
    pub sky_brightness: f64,
    /// Achieved object-to-sky flux ratio.
    pub contrast: f64,
    /// Minimum flux ratio required for detection.
    pub threshold_contrast: f64,
    /// `limiting_magnitude - object_magnitude` (positive = visible).
    pub visibility_margin: f64,
    /// Canonical adaptation state used.
    pub eye_adaptation: EyeAdaptation,
    /// Observer skill id used.
    pub observer_skill: i32,
}

// ---------------------------------------------------------------------------
// Airmass.
// ---------------------------------------------------------------------------

/// Relative optical path length through the atmosphere, dimensionless
/// Objects at or below the horizon report the
/// practical `40.0` cap on every arm. Pure kernel: non-finite altitudes
/// propagate per IEEE arithmetic (see the module-level NaN note); the
/// method is a type, so the release unknown-method refusal has no numeric
/// arm here — it lives in [`parse_airmass_method`].
pub fn calc_airmass(altitude_deg: f64, method: AirmassMethod) -> f64 {
    if altitude_deg <= 0.0 {
        return HORIZON_AIRMASS_CAP;
    }
    let zenith_deg = 90.0 - altitude_deg;
    match method {
        AirmassMethod::Secant => {
            if zenith_deg >= 85.0 {
                return HORIZON_AIRMASS_CAP;
            }
            1.0 / zenith_deg.to_radians().cos()
        }
        AirmassMethod::KastenYoung => {
            let sin_h = altitude_deg.to_radians().sin();
            let h = altitude_deg;
            let denominator = sin_h + 0.50572 * (h + 6.07995).powf(-1.6364);
            if denominator <= 0.0 {
                // Out-of-range caller altitudes can make this empirical
                // denominator non-positive; retain the horizon floor.
                return HORIZON_AIRMASS_CAP;
            }
            (1.0 / denominator).clamp(1.0, HORIZON_AIRMASS_CAP)
        }
        AirmassMethod::Rozenberg => {
            let cos_z = zenith_deg.to_radians().cos();
            (1.0 / (cos_z + 0.025 * (-11.0 * cos_z).exp())).max(1.0)
        }
    }
}

// ---------------------------------------------------------------------------
// Extinction components.
// ---------------------------------------------------------------------------

/// Rayleigh scattering coefficient in mag/airmass: `0.1451 * column * (550/lambda)^4`.
/// The air column uses exactly one of its two descriptions — a
/// non-default site pressure already encodes the reduced column, else
/// the observer altitude drives the barometric factor — selected by the
/// same `pressure != DEFAULT` float comparison as the release. `temperature_c` is
/// part of the release signature and unused there; it is unused here.
/// Wavelengths at or below zero are outside the swept domain (zero
/// divides: `+inf` here, `ZeroDivisionError` there).
pub fn calc_rayleigh_coefficient(
    wavelength_nm: f64,
    pressure_mbar: f64,
    _temperature_c: f64,
    altitude_m: f64,
) -> f64 {
    let column_factor = if pressure_mbar != DEFAULT_PRESSURE_MBAR {
        pressure_mbar / DEFAULT_PRESSURE_MBAR
    } else {
        (-altitude_m / SCALE_HEIGHT_M).exp()
    };
    let wavelength_factor = (RAYLEIGH_REF_WAVELENGTH_NM / wavelength_nm).powf(4.0);
    RAYLEIGH_REF_COEFF * column_factor * wavelength_factor
}

/// Aerosol (Mie) coefficient in mag/airmass. A positive supplied visibility drives the
/// Koschmieder estimate (`3.912/V - 0.1066`, floored at zero); otherwise
/// the humidity/altitude/wavelength estimate applies (project-calibrated
/// hygroscopic-growth fallback under the Angstrom 1929 power law).
pub fn calc_aerosol_coefficient(
    humidity_percent: f64,
    altitude_m: f64,
    wavelength_nm: f64,
    visibility_km: Option<f64>,
) -> f64 {
    if let Some(vis) = visibility_km {
        if vis > 0.0 {
            return (KOSCHMIEDER_CONSTANT / vis - KOSCHMIEDER_MOLECULAR_PART).max(0.0);
        }
    }
    let humidity_fraction = humidity_percent / 100.0;
    let humidity_factor = 1.0 + 2.0 * humidity_fraction.powf(2.0);
    let altitude_factor = (-altitude_m / AEROSOL_SCALE_HEIGHT_M).exp();
    let wavelength_factor = (WAVELENGTH_V / wavelength_nm).powf(1.3);
    AEROSOL_BASE_COEFF * humidity_factor * altitude_factor * wavelength_factor
}

/// Ozone absorption coefficient in mag/airmass: project-calibrated step approximation whose
/// ~0.016 visual-band figure matches the small constant used by Schaefer
/// (1990).
pub fn calc_ozone_coefficient(wavelength_nm: f64) -> f64 {
    if wavelength_nm < 400.0 {
        OZONE_COEFF_UV
    } else if wavelength_nm < 500.0 {
        OZONE_COEFF_BLUE
    } else if wavelength_nm < 700.0 {
        OZONE_COEFF_VISUAL
    } else {
        OZONE_COEFF_RED
    }
}

/// Water-vapour absorption coefficient in mag/airmass: project-calibrated per-band
/// coefficients scaled linearly by humidity. `temperature_c` is part of
/// the release signature and unused there; it is unused here.
pub fn calc_water_vapor_coefficient(
    humidity_percent: f64,
    _temperature_c: f64,
    wavelength_nm: f64,
) -> f64 {
    let humidity_fraction = humidity_percent / 100.0;
    if wavelength_nm < 600.0 {
        0.002 * humidity_fraction
    } else if wavelength_nm < 800.0 {
        0.005 * humidity_fraction
    } else {
        0.02 * humidity_fraction
    }
}

/// Total extinction coefficient and its components in mag/airmass
/// Units and domains are those of [`ObserverConditions::new`]; this composes the same six scalars.
#[allow(clippy::too_many_arguments)]
pub fn calc_extinction_coefficient(
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    altitude_m: f64,
    wavelength_nm: f64,
    visibility_km: Option<f64>,
) -> ExtinctionCoefficients {
    let k_rayleigh =
        calc_rayleigh_coefficient(wavelength_nm, pressure_mbar, temperature_c, altitude_m);
    let k_aerosol =
        calc_aerosol_coefficient(humidity_percent, altitude_m, wavelength_nm, visibility_km);
    let k_ozone = calc_ozone_coefficient(wavelength_nm);
    let k_water = calc_water_vapor_coefficient(humidity_percent, temperature_c, wavelength_nm);
    ExtinctionCoefficients::new(k_rayleigh, k_aerosol, k_ozone, k_water)
}

/// Atmospheric extinction for an object altitude in magnitudes: the total coefficient times the
/// Kasten-Young airmass. Altitudes at or below zero report `99.0`.
#[allow(clippy::too_many_arguments)]
pub fn calc_extinction_magnitude(
    altitude_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    observer_altitude_m: f64,
    wavelength_nm: f64,
    visibility_km: Option<f64>,
) -> f64 {
    if altitude_deg <= 0.0 {
        return BELOW_HORIZON_EXTINCTION;
    }
    let airmass = calc_airmass(altitude_deg, AirmassMethod::KastenYoung);
    let coeff = calc_extinction_coefficient(
        pressure_mbar,
        temperature_c,
        humidity_percent,
        observer_altitude_m,
        wavelength_nm,
        visibility_km,
    );
    coeff.k_total() * airmass
}

/// Extinction under the simple plane-parallel model in magnitudes
/// `k / sin(altitude)`, intended
/// above roughly 15 deg. Altitudes at or below zero report `99.0`;
/// altitudes with `sin <= 0.001` report the practical `k * 40.0` cap.
pub fn calc_simple_extinction(altitude_deg: f64, k: f64) -> f64 {
    if altitude_deg <= 0.0 {
        return BELOW_HORIZON_EXTINCTION;
    }
    let sin_alt = altitude_deg.to_radians().sin();
    if sin_alt <= 0.001 {
        return k * HORIZON_AIRMASS_CAP;
    }
    k / sin_alt
}

/// Apparent magnitude including extinction: catalog (extra-atmospheric)
/// magnitude plus the extinction dimming.
#[allow(clippy::too_many_arguments)]
pub fn apparent_magnitude_with_extinction(
    catalog_magnitude: f64,
    altitude_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    observer_altitude_m: f64,
    wavelength_nm: f64,
) -> f64 {
    let extinction = calc_extinction_magnitude(
        altitude_deg,
        pressure_mbar,
        temperature_c,
        humidity_percent,
        observer_altitude_m,
        wavelength_nm,
        None,
    );
    catalog_magnitude + extinction
}

/// Extinction coefficient for heliacal visibility calculations,
/// mag/airmass: the total
/// V-band coefficient. The zenith angle is accepted for the public interface
/// and unused, exactly like the release.
pub fn get_extinction_for_heliacal(
    _zenith_angle_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    altitude_m: f64,
) -> f64 {
    calc_extinction_coefficient(
        pressure_mbar,
        temperature_c,
        humidity_percent,
        altitude_m,
        WAVELENGTH_V,
        None,
    )
    .k_total()
}

// ---------------------------------------------------------------------------
// Twilight sky brightness.
// ---------------------------------------------------------------------------

/// Twilight phase for a Sun altitude in degrees.
/// Boundaries are strict `>` comparisons against the `TWILIGHT_*` consts, so an exact boundary value falls into the darker phase.
pub fn get_twilight_phase(sun_altitude_deg: f64) -> TwilightPhase {
    if sun_altitude_deg > TWILIGHT_CIVIL_START {
        TwilightPhase::Day
    } else if sun_altitude_deg > TWILIGHT_CIVIL_END {
        TwilightPhase::Civil
    } else if sun_altitude_deg > TWILIGHT_NAUTICAL_END {
        TwilightPhase::Nautical
    } else if sun_altitude_deg > TWILIGHT_ASTRONOMICAL_END {
        TwilightPhase::Astronomical
    } else {
        TwilightPhase::Night
    }
}

/// Azimuthal brightness factor from the angular distance to the Sun
/// The gradient
/// weight follows the Sun altitude: 0.7 above -6 deg, 0.4 above -12 deg,
/// 0.2 above -18 deg, 0.0 at night.
/// The unweighted factor runs 0 (toward the Sun) to 1 (away); the altitude weight blends it toward 1.
fn azimuth_factor(sun_azimuth_deg: f64, target_azimuth_deg: f64, sun_altitude_deg: f64) -> f64 {
    let mut delta_az = (sun_azimuth_deg - target_azimuth_deg).abs();
    if delta_az > 180.0 {
        delta_az = 360.0 - delta_az;
    }
    let az_factor = delta_az / 180.0;
    let weight = if sun_altitude_deg > -6.0 {
        0.7
    } else if sun_altitude_deg > -12.0 {
        0.4
    } else if sun_altitude_deg > -18.0 {
        0.2
    } else {
        0.0
    };
    az_factor * weight + (1.0 - weight)
}

/// Target-altitude brightness factor: the Van Rhijn secant brightening toward the
/// horizon, 0.5 below the horizon, 1.5 past 85 deg zenith angle.
fn altitude_factor(target_altitude_deg: f64) -> f64 {
    if target_altitude_deg <= 0.0 {
        return 0.5;
    }
    let zenith_angle = 90.0 - target_altitude_deg;
    if zenith_angle >= 85.0 {
        return 1.5;
    }
    let cos_z = zenith_angle.to_radians().cos();
    let sec_z = if cos_z > 0.1 {
        (1.0 / cos_z).min(10.0)
    } else {
        10.0
    };
    (0.2 * (sec_z - 1.0) + 1.0).min(2.0)
}

/// Sky surface brightness during twilight. The base brightness interpolates the
/// zenith anchors across the twilight phases (project-calibrated
/// representative values, trend-consistent with Patat 2003 and Rozenberg
/// 1966); the azimuthal, altitude and aerosol adjustments, the `22.5`
/// clamp, the piecewise-linear limiting-magnitude fit (a project fit of
/// limiting magnitude against mag/arcsec^2, not Schaefer's closed form)
/// with its `[-2, 7]` clamp, and the Garstang nanoLambert conversion all
/// follow the release spelling.
/// Altitudes and azimuths are in degrees, pressure in mbar, humidity in percent, altitude in metres, wavelength in nm; the release defaults (zenith target, opposite-Sun azimuths, standard atmosphere, V-band) are explicit arguments here, and the daytime arm darkens linearly from 3.0 with a 0.0 floor.
#[allow(clippy::too_many_arguments)]
pub fn calc_twilight_sky_brightness(
    sun_altitude_deg: f64,
    target_altitude_deg: f64,
    sun_azimuth_deg: f64,
    target_azimuth_deg: f64,
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_percent: f64,
    altitude_m: f64,
    wavelength_nm: f64,
) -> TwilightSkyBrightness {
    let phase = get_twilight_phase(sun_altitude_deg);
    let base_brightness = if sun_altitude_deg >= 0.0 {
        (3.0 - sun_altitude_deg * 0.05).max(0.0)
    } else if sun_altitude_deg >= TWILIGHT_CIVIL_END {
        let fraction = -sun_altitude_deg / 6.0;
        SKY_BRIGHTNESS_SUN_0 + fraction * (SKY_BRIGHTNESS_SUN_MINUS6 - SKY_BRIGHTNESS_SUN_0)
    } else if sun_altitude_deg >= TWILIGHT_NAUTICAL_END {
        let fraction = (-sun_altitude_deg - 6.0) / 6.0;
        SKY_BRIGHTNESS_SUN_MINUS6
            + fraction * (SKY_BRIGHTNESS_SUN_MINUS12 - SKY_BRIGHTNESS_SUN_MINUS6)
    } else if sun_altitude_deg >= TWILIGHT_ASTRONOMICAL_END {
        let fraction = (-sun_altitude_deg - 12.0) / 6.0;
        SKY_BRIGHTNESS_SUN_MINUS12
            + fraction * (SKY_BRIGHTNESS_SUN_MINUS18 - SKY_BRIGHTNESS_SUN_MINUS12)
    } else {
        DARK_SKY_BRIGHTNESS_V
    };

    let az_factor = azimuth_factor(sun_azimuth_deg, target_azimuth_deg, sun_altitude_deg);
    let alt_factor = altitude_factor(target_altitude_deg);

    let az_adjustment = if sun_altitude_deg >= TWILIGHT_CIVIL_END {
        (1.0 - az_factor) * 2.0
    } else if sun_altitude_deg >= TWILIGHT_NAUTICAL_END {
        (1.0 - az_factor) * 1.0
    } else if sun_altitude_deg >= TWILIGHT_ASTRONOMICAL_END {
        (1.0 - az_factor) * 0.3
    } else {
        0.0
    };

    let alt_adjustment = if alt_factor > 1.0 {
        -(alt_factor - 1.0) * 1.5
    } else {
        0.0
    };

    let coeff = calc_extinction_coefficient(
        pressure_mbar,
        temperature_c,
        humidity_percent,
        altitude_m,
        wavelength_nm,
        None,
    );
    let aerosol_adjustment = -(coeff.k_aerosol() - 0.1) * 0.5;

    let surface_brightness =
        (base_brightness - az_adjustment + alt_adjustment + aerosol_adjustment).clamp(0.0, 22.5);

    let limiting_mag = if surface_brightness < 10.0 {
        -2.0 + surface_brightness * 0.5
    } else if surface_brightness < 18.0 {
        3.0 + (surface_brightness - 10.0) * 0.35
    } else {
        5.8 + (surface_brightness - 18.0) * 0.2
    };
    let limiting_magnitude = limiting_mag.clamp(-2.0, 7.0);

    let nanolamberts = 10.0_f64.powf((NANOLAMBERT_ZERO_POINT - surface_brightness) / 2.5);

    TwilightSkyBrightness {
        surface_brightness,
        twilight_phase: phase,
        sun_altitude_deg,
        azimuth_factor: az_factor,
        limiting_magnitude,
        nanolamberts,
    }
}

/// Approximate zenith sky brightness during twilight, mag/arcsec^2
/// the full model at the
/// zenith opposite the Sun.
pub fn calc_twilight_brightness_simple(sun_altitude_deg: f64) -> f64 {
    calc_twilight_sky_brightness(
        sun_altitude_deg,
        90.0,
        0.0,
        180.0,
        DEFAULT_PRESSURE_MBAR,
        DEFAULT_TEMPERATURE_C,
        DEFAULT_HUMIDITY_PERCENT,
        DEFAULT_ALTITUDE_M,
        WAVELENGTH_V,
    )
    .surface_brightness
}

/// Experience factor for an observer skill id: unknown ids default to
/// `1.0`, the same default as the release.
pub fn experience_factor(observer_skill: i32) -> f64 {
    match observer_skill {
        OBSERVER_SKILL_INEXPERIENCED => 1.3,
        OBSERVER_SKILL_AVERAGE => 1.0,
        OBSERVER_SKILL_EXPERIENCED => 0.85,
        OBSERVER_SKILL_EXPERT => 0.7,
        _ => 1.0,
    }
}

/// Ordered contents of the release `EXPERIENCE_FACTORS` table: one
/// `(skill id, factor)` row per known observer skill, in the release
/// dict insertion order (inexperienced, average, experienced, expert).
/// Unknown skill ids are not table rows; they resolve to the `1.0`
/// default through [`experience_factor`], exactly like the release
/// `EXPERIENCE_FACTORS.get(observer_skill, 1.0)` lookup.
pub const EXPERIENCE_TABLE: [(i32, f64); 4] = [
    (OBSERVER_SKILL_INEXPERIENCED, 1.3),
    (OBSERVER_SKILL_AVERAGE, 1.0),
    (OBSERVER_SKILL_EXPERIENCED, 0.85),
    (OBSERVER_SKILL_EXPERT, 0.7),
];

/// Number of rows in [`EXPERIENCE_TABLE`].
pub fn experience_table_len() -> usize {
    EXPERIENCE_TABLE.len()
}

/// Row of [`EXPERIENCE_TABLE`] by position (`None` past the end).
pub fn experience_at(index: usize) -> Option<(i32, f64)> {
    EXPERIENCE_TABLE.get(index).copied()
}

/// Naked-eye limiting magnitude during twilight: the sky-brightness fit minus the
/// extinction penalty `k_total * (X - 1)` for the viewing direction,
/// clamped to `[-2, 6.5]`.
/// Temperature and wavelength are fixed to `DEFAULT_TEMPERATURE_C` and `WAVELENGTH_V`; remaining altitudes are in degrees, pressure in mbar, humidity in percent, altitude in metres.
#[allow(clippy::too_many_arguments)]
pub fn calc_limiting_magnitude_twilight(
    sun_altitude_deg: f64,
    target_altitude_deg: f64,
    sun_azimuth_deg: f64,
    target_azimuth_deg: f64,
    pressure_mbar: f64,
    humidity_percent: f64,
    altitude_m: f64,
) -> f64 {
    let result = calc_twilight_sky_brightness(
        sun_altitude_deg,
        target_altitude_deg,
        sun_azimuth_deg,
        target_azimuth_deg,
        pressure_mbar,
        DEFAULT_TEMPERATURE_C,
        humidity_percent,
        altitude_m,
        WAVELENGTH_V,
    );
    let airmass = calc_airmass(target_altitude_deg, AirmassMethod::KastenYoung);
    let coeff = calc_extinction_coefficient(
        pressure_mbar,
        DEFAULT_TEMPERATURE_C,
        humidity_percent,
        altitude_m,
        WAVELENGTH_V,
        None,
    );
    let extinction_penalty = coeff.k_total() * (airmass - 1.0);
    let adjusted_limit = result.limiting_magnitude - extinction_penalty;
    adjusted_limit.clamp(-2.0, 6.5)
}

// ---------------------------------------------------------------------------
// Visibility threshold.
// ---------------------------------------------------------------------------

/// Eye adaptation state for a sky brightness, mag/arcsec^2: photopic below 16, mesopic below 20,
/// dark (the scotopic regime) at and above 20.
pub fn calc_eye_adaptation_state(sky_brightness: f64) -> EyeAdaptation {
    if sky_brightness < 16.0 {
        EyeAdaptation::Photopic
    } else if sky_brightness < 20.0 {
        EyeAdaptation::Mesopic
    } else {
        EyeAdaptation::Dark
    }
}

/// Point-source detection threshold as a star-to-sky flux ratio
/// the model limiting magnitude for
/// the adaptation regime and observer skill, inverted through
/// `10^((B - m_lim)/2.5)` and clamped to `[1e-8, 1e10]`. Brighter skies
/// and less-adapted eyes give larger thresholds. This is the inverse of
/// the visibility limiting-magnitude logic, not the classical Blackwell
/// contrast fraction. `None` auto-detects the regime from the sky
/// brightness; explicit states are already canonical here (the string
/// refusal arm lives in [`normalize_eye_adaptation`]).
/// Branches: dark `6.0 + (B - 21.5) * 0.3`, mesopic `3.4 + (B - 16.0) * 0.5`, photopic three-segment descent; skill adjustment `(1 - factor) * 3.0` with the limiting magnitude clamped to `[-5.0, 7.5]`.
pub fn calc_contrast_threshold(
    sky_brightness: f64,
    eye_adaptation: Option<EyeAdaptation>,
    observer_skill: i32,
) -> f64 {
    calc_contrast_threshold_with_factor(
        sky_brightness,
        eye_adaptation,
        experience_factor(observer_skill),
    )
}

/// Detection contrast using an already consumed observer experience factor.
/// Sky brightness is magnitudes per square arcsecond; the dimensionless
/// `factor` contributes `(1 - factor) * 3` magnitudes before the same limiting
/// magnitude and contrast clamps. `None` selects eye adaptation from the sky.
/// This entry performs no skill-table lookup or caller-object conversion.
pub fn calc_contrast_threshold_with_factor(
    sky_brightness: f64,
    eye_adaptation: Option<EyeAdaptation>,
    factor: f64,
) -> f64 {
    let exponent = contrast_threshold_exponent_with_factor(sky_brightness, eye_adaptation, factor);
    10.0_f64.powf(exponent).clamp(1e-8, 1e10)
}

/// Base-ten contrast exponent from sky surface brightness, canonical eye
/// adaptation and an already consumed dimensionless experience factor.
/// This retains the limiting magnitude assembly and its [-5, 7.5] clamp,
/// returning `(sky - limiting) / 2.5` before any power or contrast clamp.
/// It performs no table lookup or power evaluation; callers can preserve
/// their arithmetic backend's actual overflow or operator exceptions.
pub fn contrast_threshold_exponent_with_factor(
    sky_brightness: f64,
    eye_adaptation: Option<EyeAdaptation>,
    factor: f64,
) -> f64 {
    let state = eye_adaptation.unwrap_or_else(|| calc_eye_adaptation_state(sky_brightness));
    let limiting_mag = match state {
        EyeAdaptation::Dark => 6.0 + (sky_brightness - 21.5) * 0.3,
        EyeAdaptation::Mesopic => 3.4 + (sky_brightness - 16.0) * 0.5,
        EyeAdaptation::Photopic => {
            if sky_brightness >= 10.0 {
                1.0 + (sky_brightness - 10.0) * 0.4
            } else if sky_brightness >= 5.0 {
                -2.0 + (sky_brightness - 5.0) * 0.6
            } else {
                -4.5 + sky_brightness * 0.5
            }
        }
    };
    let experience_adj = (1.0 - factor) * 3.0;
    let limiting_mag = (limiting_mag + experience_adj).clamp(-5.0, 7.5);
    let delta_m = sky_brightness - limiting_mag;
    delta_m / 2.5
}

/// Visibility-threshold decision:
/// a compact Schaefer-inspired convenience model, not Schaefer's full
/// retinal/telescope model. Pogson fluxes carry a shared unit
/// normalization (it cancels in the ratio); the sky value is treated as
/// the flux of a one-arcsec^2 patch, a stable API approximation. The
/// object is visible under the strict `magnitude < limiting` comparison.
#[allow(clippy::too_many_arguments)]
pub fn calc_visibility_threshold(
    object_magnitude: f64,
    sky_brightness: f64,
    eye_adaptation: Option<EyeAdaptation>,
    observer_skill: i32,
    object_altitude_deg: f64,
    apply_extinction: bool,
    pressure_mbar: f64,
    humidity_percent: f64,
    observer_altitude_m: f64,
) -> VisibilityResult {
    calc_visibility_threshold_with_factor(
        object_magnitude,
        sky_brightness,
        eye_adaptation,
        observer_skill,
        object_altitude_deg,
        apply_extinction,
        pressure_mbar,
        humidity_percent,
        observer_altitude_m,
        experience_factor(observer_skill),
    )
}

/// Visibility with an already consumed experience factor. All geometry,
/// atmosphere, magnitude and strict visibility semantics match
/// [`calc_visibility_threshold`]. `observer_skill` is retained result metadata;
/// only the supplied dimensionless `factor` affects the contrast threshold.
/// The caller owns any mutable table lookup and its failure chronology.
#[allow(clippy::too_many_arguments)]
pub fn calc_visibility_threshold_with_factor(
    object_magnitude: f64,
    sky_brightness: f64,
    eye_adaptation: Option<EyeAdaptation>,
    observer_skill: i32,
    object_altitude_deg: f64,
    apply_extinction: bool,
    pressure_mbar: f64,
    humidity_percent: f64,
    observer_altitude_m: f64,
    factor: f64,
) -> VisibilityResult {
    let mut effective_magnitude = object_magnitude;
    if apply_extinction && object_altitude_deg > 0.0 {
        let extinction = calc_extinction_magnitude(
            object_altitude_deg,
            pressure_mbar,
            DEFAULT_TEMPERATURE_C,
            humidity_percent,
            observer_altitude_m,
            WAVELENGTH_V,
            None,
        );
        effective_magnitude = object_magnitude + extinction;
    }

    let state = eye_adaptation.unwrap_or_else(|| calc_eye_adaptation_state(sky_brightness));
    let threshold = calc_contrast_threshold_with_factor(sky_brightness, Some(state), factor);

    let sky_flux = 10.0_f64.powf(-0.4 * sky_brightness);
    let object_flux = 10.0_f64.powf(-0.4 * effective_magnitude);
    let contrast = if sky_flux > 0.0 {
        object_flux / sky_flux
    } else {
        f64::INFINITY
    };

    let limiting_magnitude = if threshold > 0.0 {
        sky_brightness - 2.5 * threshold.log10()
    } else {
        sky_brightness + 10.0
    };

    let is_visible = effective_magnitude < limiting_magnitude;
    let visibility_margin = limiting_magnitude - effective_magnitude;

    VisibilityResult {
        is_visible,
        object_magnitude: effective_magnitude,
        limiting_magnitude,
        sky_brightness,
        contrast,
        threshold_contrast: threshold,
        visibility_margin,
        eye_adaptation: state,
        observer_skill,
    }
}

/// Convenience visibility check: the
/// threshold decision with default atmosphere, reduced to its boolean.
pub fn is_object_visible(
    object_magnitude: f64,
    sky_brightness: f64,
    observer_skill: i32,
    object_altitude_deg: f64,
    apply_extinction: bool,
) -> bool {
    calc_visibility_threshold(
        object_magnitude,
        sky_brightness,
        None,
        observer_skill,
        object_altitude_deg,
        apply_extinction,
        DEFAULT_PRESSURE_MBAR,
        DEFAULT_HUMIDITY_PERCENT,
        DEFAULT_ALTITUDE_M,
    )
    .is_visible
}

/// Limiting magnitude for a sky brightness: the threshold decision for a
/// reference magnitude-0 object, reduced to its limiting magnitude.
pub fn calc_limiting_magnitude_for_sky(
    sky_brightness: f64,
    observer_skill: i32,
    eye_adaptation: Option<EyeAdaptation>,
) -> f64 {
    calc_visibility_threshold(
        0.0,
        sky_brightness,
        eye_adaptation,
        observer_skill,
        90.0,
        false,
        DEFAULT_PRESSURE_MBAR,
        DEFAULT_HUMIDITY_PERCENT,
        DEFAULT_ALTITUDE_M,
    )
    .limiting_magnitude
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_coeff() -> ExtinctionCoefficients {
        calc_extinction_coefficient(
            DEFAULT_PRESSURE_MBAR,
            DEFAULT_TEMPERATURE_C,
            DEFAULT_HUMIDITY_PERCENT,
            DEFAULT_ALTITUDE_M,
            WAVELENGTH_V,
            None,
        )
    }

    #[test]
    fn airmass_horizon_cap_on_every_arm() {
        for method in [
            AirmassMethod::Secant,
            AirmassMethod::KastenYoung,
            AirmassMethod::Rozenberg,
        ] {
            assert_eq!(calc_airmass(0.0, method), HORIZON_AIRMASS_CAP);
            assert_eq!(calc_airmass(-5.0, method), HORIZON_AIRMASS_CAP);
        }
        // Secant low-altitude guard: zenith >= 85 deg caps as well.
        assert_eq!(
            calc_airmass(5.0, AirmassMethod::Secant),
            HORIZON_AIRMASS_CAP
        );
        assert_eq!(
            calc_airmass(4.9, AirmassMethod::Secant),
            HORIZON_AIRMASS_CAP
        );
    }

    #[test]
    fn airmass_zenith_is_unity() {
        // The Kasten-Young/Rozenberg empirical denominators sit below 1 at
        // the zenith, so the release 1.0 floor reports exactly 1.0.
        assert_eq!(calc_airmass(90.0, AirmassMethod::KastenYoung), 1.0);
        assert_eq!(calc_airmass(90.0, AirmassMethod::Rozenberg), 1.0);
        assert_eq!(calc_airmass(90.0, AirmassMethod::Secant), 1.0);
    }

    #[test]
    fn airmass_kasten_young_spelling() {
        // Operation-order pin: recompute with the release association.
        let alt = 30.0_f64;
        let denom = alt.to_radians().sin() + 0.50572 * (alt + 6.07995).powf(-1.6364);
        let expected = (1.0 / denom).clamp(1.0, 40.0);
        assert_eq!(
            calc_airmass(alt, AirmassMethod::KastenYoung).to_bits(),
            expected.to_bits()
        );
        // Nominal value near 2 at 30 deg altitude.
        let x = calc_airmass(30.0, AirmassMethod::KastenYoung);
        assert!(x > 1.9 && x < 2.1, "airmass {x}");
    }

    #[test]
    fn method_tokens_are_exact_and_case_sensitive() {
        assert_eq!(parse_airmass_method("secant"), Some(AirmassMethod::Secant));
        assert_eq!(
            parse_airmass_method("kasten_young"),
            Some(AirmassMethod::KastenYoung)
        );
        assert_eq!(
            parse_airmass_method("rozenberg"),
            Some(AirmassMethod::Rozenberg)
        );
        assert_eq!(parse_airmass_method("Kasten_Young"), None);
        assert_eq!(parse_airmass_method(""), None);
        assert_eq!(parse_airmass_method("sec"), None);
    }

    #[test]
    fn rayleigh_defaults_and_branches() {
        // Standard conditions reproduce the Hayes & Latham V-band figure.
        assert_eq!(
            calc_rayleigh_coefficient(550.0, 1013.25, 15.0, 0.0).to_bits(),
            0.1451_f64.to_bits()
        );
        // Blue light scatters more than red.
        let blue = calc_rayleigh_coefficient(440.0, 1013.25, 15.0, 0.0);
        let red = calc_rayleigh_coefficient(640.0, 1013.25, 15.0, 0.0);
        assert!(blue > 0.3 && blue < 0.4, "blue {blue}");
        assert!(red > 0.05 && red < 0.1, "red {red}");
        // A supplied site pressure encodes the column: observer altitude
        // is then ignored, release branch.
        let low = calc_rayleigh_coefficient(550.0, 900.0, 15.0, 0.0);
        let high = calc_rayleigh_coefficient(550.0, 900.0, 15.0, 5000.0);
        assert_eq!(low.to_bits(), high.to_bits());
        assert!((low - 0.1451 * 900.0 / 1013.25).abs() < 1e-15);
        // Default pressure derives the column from altitude instead.
        let elevated = calc_rayleigh_coefficient(550.0, 1013.25, 15.0, 8500.0);
        assert!((elevated - 0.1451 / std::f64::consts::E).abs() < 1e-12);
    }

    #[test]
    fn aerosol_visibility_and_fallback_arms() {
        // Koschmieder arm: 3.912/V minus the molecular part, floored.
        let k = calc_aerosol_coefficient(50.0, 0.0, 550.0, Some(20.0));
        assert!((k - (3.912 / 20.0 - 0.1066)).abs() < 1e-15, "k {k}");
        assert_eq!(calc_aerosol_coefficient(50.0, 0.0, 550.0, Some(200.0)), 0.0);
        // Non-positive visibility falls back to the humidity estimate.
        let fallback = calc_aerosol_coefficient(50.0, 0.0, 550.0, Some(0.0));
        let direct = calc_aerosol_coefficient(50.0, 0.0, 550.0, None);
        assert_eq!(fallback.to_bits(), direct.to_bits());
        // Humidity estimate at defaults: 0.08 * 1.5 = 0.12.
        assert!((direct - 0.12).abs() < 1e-12, "direct {direct}");
        // Dry air is clearer; altitude thins the aerosol column.
        let dry = calc_aerosol_coefficient(0.0, 0.0, 550.0, None);
        assert!(dry < direct);
        let high = calc_aerosol_coefficient(50.0, 1500.0, 550.0, None);
        assert!((high - direct / std::f64::consts::E).abs() < 1e-12);
    }

    #[test]
    fn ozone_and_water_steps() {
        assert_eq!(calc_ozone_coefficient(300.0), OZONE_COEFF_UV);
        assert_eq!(calc_ozone_coefficient(399.9), OZONE_COEFF_UV);
        assert_eq!(calc_ozone_coefficient(400.0), OZONE_COEFF_BLUE);
        assert_eq!(calc_ozone_coefficient(499.9), OZONE_COEFF_BLUE);
        assert_eq!(calc_ozone_coefficient(500.0), OZONE_COEFF_VISUAL);
        assert_eq!(calc_ozone_coefficient(550.0), 0.016);
        assert_eq!(calc_ozone_coefficient(699.9), OZONE_COEFF_VISUAL);
        assert_eq!(calc_ozone_coefficient(700.0), OZONE_COEFF_RED);
        // Water vapour scales linearly with humidity within each band.
        assert_eq!(calc_water_vapor_coefficient(0.0, 15.0, 550.0), 0.0);
        assert!((calc_water_vapor_coefficient(50.0, 15.0, 550.0) - 0.001).abs() < 1e-15);
        assert!((calc_water_vapor_coefficient(100.0, 15.0, 700.0) - 0.005).abs() < 1e-15);
        assert!((calc_water_vapor_coefficient(50.0, 15.0, 900.0) - 0.01).abs() < 1e-15);
    }

    #[test]
    fn coefficient_total_is_the_component_sum() {
        let coeff = default_coeff();
        assert_eq!(
            coeff.k_total().to_bits(),
            (coeff.k_rayleigh() + coeff.k_aerosol() + coeff.k_ozone() + coeff.k_water()).to_bits()
        );
        // V-band defaults: ~0.1451 + 0.12 + 0.016 + 0.001.
        assert!(
            (coeff.k_total() - 0.2821).abs() < 1e-4,
            "total {}",
            coeff.k_total()
        );
        let via_conditions = ObserverConditions::new(
            DEFAULT_PRESSURE_MBAR,
            DEFAULT_TEMPERATURE_C,
            DEFAULT_HUMIDITY_PERCENT,
            DEFAULT_ALTITUDE_M,
            WAVELENGTH_V,
            None,
        )
        .unwrap()
        .coefficients();
        assert_eq!(via_conditions, coeff);
    }

    #[test]
    fn extinction_magnitude_guards_and_composition() {
        assert_eq!(
            calc_extinction_magnitude(0.0, 1013.25, 15.0, 50.0, 0.0, 550.0, None),
            BELOW_HORIZON_EXTINCTION
        );
        assert_eq!(
            calc_extinction_magnitude(-3.0, 1013.25, 15.0, 50.0, 0.0, 550.0, None),
            BELOW_HORIZON_EXTINCTION
        );
        // Composition identity: k_total times the Kasten-Young airmass.
        let coeff = default_coeff();
        for alt in [0.5, 5.0, 30.0, 90.0] {
            let mag = calc_extinction_magnitude(alt, 1013.25, 15.0, 50.0, 0.0, 550.0, None);
            let expected = coeff.k_total() * calc_airmass(alt, AirmassMethod::KastenYoung);
            assert_eq!(mag.to_bits(), expected.to_bits(), "alt {alt}");
        }
        // Zenith extinction is about one total coefficient.
        let zenith = calc_extinction_magnitude(90.0, 1013.25, 15.0, 50.0, 0.0, 550.0, None);
        assert!((zenith - coeff.k_total()).abs() < 1e-12);
    }

    #[test]
    fn simple_extinction_branches() {
        assert_eq!(calc_simple_extinction(0.0, 0.28), BELOW_HORIZON_EXTINCTION);
        assert_eq!(calc_simple_extinction(-1.0, 0.28), BELOW_HORIZON_EXTINCTION);
        // Near-horizon sine floor reports the practical cap.
        assert_eq!(calc_simple_extinction(0.05, 0.28), 0.28 * 40.0);
        // Nominal arm: k / sin(alt).
        let m = calc_simple_extinction(30.0, 0.28);
        assert!((m - 0.56).abs() < 1e-12, "m {m}");
    }

    #[test]
    fn apparent_magnitude_adds_the_dimming() {
        let app = apparent_magnitude_with_extinction(-4.0, 30.0, 1013.25, 15.0, 50.0, 0.0, 550.0);
        let dim = calc_extinction_magnitude(30.0, 1013.25, 15.0, 50.0, 0.0, 550.0, None);
        assert_eq!(app.to_bits(), (-4.0 + dim).to_bits());
        assert!(app > -4.0);
        // Below the horizon the sentinel dominates the catalog value.
        let sunk = apparent_magnitude_with_extinction(6.0, -2.0, 1013.25, 15.0, 50.0, 0.0, 550.0);
        assert_eq!(sunk, 6.0 + BELOW_HORIZON_EXTINCTION);
    }

    #[test]
    fn heliacal_coefficient_ignores_zenith_angle() {
        let a = get_extinction_for_heliacal(60.0, 1013.25, 15.0, 50.0, 0.0);
        let b = get_extinction_for_heliacal(10.0, 1013.25, 15.0, 50.0, 0.0);
        assert_eq!(a.to_bits(), b.to_bits());
        assert_eq!(a.to_bits(), default_coeff().k_total().to_bits());
    }

    #[test]
    fn conditions_constructor_refuses_non_finite() {
        for bad in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                ObserverConditions::new(bad, 15.0, 50.0, 0.0, 550.0, None),
                Err(Error::NonFiniteArgument { .. })
            ));
            assert!(matches!(
                ObserverConditions::new(1013.25, 15.0, 50.0, 0.0, 550.0, Some(bad)),
                Err(Error::NonFiniteArgument { .. })
            ));
        }
        // Finite out-of-range words are admitted like the release.
        assert!(ObserverConditions::new(1013.25, 15.0, 150.0, 0.0, 550.0, None).is_ok());
        assert!(ObserverConditions::new(-5.0, 15.0, 50.0, 0.0, 550.0, None).is_ok());
    }

    #[test]
    fn twilight_phase_boundaries() {
        assert_eq!(get_twilight_phase(5.0).as_str(), "day");
        assert_eq!(get_twilight_phase(0.0).as_str(), "civil");
        assert_eq!(get_twilight_phase(-3.0).as_str(), "civil");
        assert_eq!(get_twilight_phase(-6.0).as_str(), "nautical");
        assert_eq!(get_twilight_phase(-12.0).as_str(), "astronomical");
        assert_eq!(get_twilight_phase(-18.0).as_str(), "night");
        assert_eq!(get_twilight_phase(-25.0).as_str(), "night");
    }

    #[test]
    fn twilight_anchors_and_clamps() {
        // Anchor interpolation: Sun at horizon, zenith opposite the Sun.
        let civil =
            calc_twilight_sky_brightness(0.0, 90.0, 0.0, 180.0, 1013.25, 15.0, 50.0, 0.0, 550.0);
        assert_eq!(civil.twilight_phase, TwilightPhase::Civil);
        assert!(
            (civil.surface_brightness - 2.99).abs() < 1e-12,
            "civil {}",
            civil.surface_brightness
        );
        assert_eq!(civil.sun_altitude_deg, 0.0);
        assert_eq!(civil.azimuth_factor, 1.0);
        // End of astronomical twilight settles at the dark anchor minus
        // the default aerosol adjustment.
        let dark = calc_twilight_brightness_simple(-18.0);
        assert!((dark - 21.69).abs() < 1e-9, "dark {dark}");
        // Full night never exceeds the 22.5 clamp; day never goes negative.
        let night = calc_twilight_brightness_simple(-30.0);
        assert!((21.0..=22.5).contains(&night), "night {night}");
        let day = calc_twilight_brightness_simple(60.0);
        assert!((0.0..3.0).contains(&day), "day {day}");
        // Looking toward the Sun is brighter than looking away.
        let toward =
            calc_twilight_sky_brightness(-3.0, 45.0, 270.0, 270.0, 1013.25, 15.0, 50.0, 0.0, 550.0);
        let away =
            calc_twilight_sky_brightness(-3.0, 45.0, 270.0, 90.0, 1013.25, 15.0, 50.0, 0.0, 550.0);
        assert!(toward.surface_brightness < away.surface_brightness);
        assert!(toward.azimuth_factor < away.azimuth_factor);
        let unknown_direction = calc_twilight_sky_brightness(
            -3.0,
            45.0,
            f64::NAN,
            90.0,
            1013.25,
            15.0,
            50.0,
            0.0,
            550.0,
        );
        assert!(unknown_direction.surface_brightness.is_nan());
        assert!(unknown_direction.limiting_magnitude.is_nan());
        assert!(unknown_direction.nanolamberts.is_nan());
    }

    #[test]
    fn twilight_limiting_range() {
        // Civil twilight admits only bright objects; full night reaches
        // the naked-eye band (extinction penalty applied at 45 deg).
        let civil = calc_limiting_magnitude_twilight(-3.0, 45.0, 0.0, 180.0, 1013.25, 50.0, 0.0);
        assert!(civil > 0.0 && civil < 4.5, "civil {civil}");
        let night = calc_limiting_magnitude_twilight(-25.0, 45.0, 0.0, 180.0, 1013.25, 50.0, 0.0);
        assert!(night > 5.0 && night <= 6.5, "night {night}");
        // Composition identity with the sky model and the penalty.
        let sky =
            calc_twilight_sky_brightness(-9.0, 45.0, 0.0, 180.0, 1013.25, 15.0, 50.0, 0.0, 550.0);
        let penalty =
            default_coeff().k_total() * (calc_airmass(45.0, AirmassMethod::KastenYoung) - 1.0);
        let expected = (sky.limiting_magnitude - penalty).clamp(-2.0, 6.5);
        assert_eq!(
            calc_limiting_magnitude_twilight(-9.0, 45.0, 0.0, 180.0, 1013.25, 50.0, 0.0).to_bits(),
            expected.to_bits()
        );
    }

    #[test]
    fn eye_adaptation_tokens_and_states() {
        assert_eq!(
            normalize_eye_adaptation("photopic"),
            Some(EyeAdaptation::Photopic)
        );
        assert_eq!(
            normalize_eye_adaptation("mesopic"),
            Some(EyeAdaptation::Mesopic)
        );
        assert_eq!(normalize_eye_adaptation("dark"), Some(EyeAdaptation::Dark));
        assert_eq!(
            normalize_eye_adaptation("scotopic"),
            Some(EyeAdaptation::Dark)
        );
        assert_eq!(normalize_eye_adaptation("SCOTOPIC"), None);
        assert_eq!(normalize_eye_adaptation(""), None);
        assert_eq!(normalize_eye_adaptation("scotopic "), None);
        assert_eq!(calc_eye_adaptation_state(15.9), EyeAdaptation::Photopic);
        assert_eq!(calc_eye_adaptation_state(16.0), EyeAdaptation::Mesopic);
        assert_eq!(calc_eye_adaptation_state(19.9), EyeAdaptation::Mesopic);
        assert_eq!(calc_eye_adaptation_state(20.0), EyeAdaptation::Dark);
    }

    #[test]
    fn contrast_threshold_reference_points() {
        // Dark sky, average observer: 10^((21.5-6.0)/2.5), the release
        // docstring value.
        let dark = calc_contrast_threshold(21.5, None, OBSERVER_SKILL_AVERAGE);
        assert_eq!(dark.round(), 1_584_893.0);
        // Nautical twilight, experienced observer: the second release
        // docstring value (auto state is photopic at B=15).
        let twilight = calc_contrast_threshold(15.0, None, OBSERVER_SKILL_EXPERIENCED);
        assert_eq!(twilight.round(), 41_687.0);
        // Explicit states select their limiting-magnitude branch.
        let mesopic = calc_contrast_threshold(18.0, Some(EyeAdaptation::Mesopic), 2);
        let expected = 10.0_f64.powf((18.0 - (3.4 + 2.0 * 0.5)) / 2.5);
        assert_eq!(mesopic.to_bits(), expected.to_bits());
        // Skill ordering: experts detect fainter (smaller threshold).
        let expert = calc_contrast_threshold(21.5, None, OBSERVER_SKILL_EXPERT);
        let novice = calc_contrast_threshold(21.5, None, OBSERVER_SKILL_INEXPERIENCED);
        assert!(expert < dark && dark < novice);
        // Unknown skill ids fall back to the average factor.
        assert_eq!(
            calc_contrast_threshold(21.5, None, 99).to_bits(),
            dark.to_bits()
        );
        // Thresholds stay inside the release clamp on extreme skies.
        assert!(calc_contrast_threshold(40.0, None, 2) <= 1e10);
        assert!(calc_contrast_threshold(-40.0, None, 2) >= 1e-8);
    }

    #[test]
    fn experience_table_matches_factor_lookup() {
        // Release `EXPERIENCE_FACTORS` insertion order with exact words.
        assert_eq!(experience_table_len(), 4);
        let expected = [
            (OBSERVER_SKILL_INEXPERIENCED, 1.3_f64),
            (OBSERVER_SKILL_AVERAGE, 1.0_f64),
            (OBSERVER_SKILL_EXPERIENCED, 0.85_f64),
            (OBSERVER_SKILL_EXPERT, 0.7_f64),
        ];
        for (i, (skill, factor)) in expected.iter().enumerate() {
            let (got_skill, got_factor) = experience_at(i).expect("row present");
            assert_eq!(got_skill, *skill);
            assert_eq!(got_factor.to_bits(), factor.to_bits());
            assert_eq!(experience_factor(*skill).to_bits(), factor.to_bits());
        }
        assert!(experience_at(4).is_none());
        // Unknown ids are not rows; the lookup defaults to 1.0.
        assert_eq!(experience_factor(99).to_bits(), 1.0_f64.to_bits());
    }

    #[test]
    fn visibility_docstring_outcomes() {
        // Venus in civil twilight: visible; a 6th-mag star is not.
        assert!(
            calc_visibility_threshold(-4.0, 8.0, None, 2, 90.0, false, 1013.25, 50.0, 0.0)
                .is_visible
        );
        assert!(
            !calc_visibility_threshold(6.0, 8.0, None, 2, 90.0, false, 1013.25, 50.0, 0.0)
                .is_visible
        );
        // 4th magnitude under a dark sky: visible.
        let dark = calc_visibility_threshold(4.0, 21.5, None, 2, 90.0, false, 1013.25, 50.0, 0.0);
        assert!(dark.is_visible);
        assert!(dark.visibility_margin > 0.0);
        assert_eq!(dark.eye_adaptation, EyeAdaptation::Dark);
        // An expert reaches a 6th-magnitude star the average eye misses
        // (the comparison is strict: mag == limiting is not visible).
        assert!(
            calc_visibility_threshold(6.0, 21.5, None, 4, 90.0, false, 1013.25, 50.0, 0.0)
                .is_visible
        );
        assert!(
            !calc_visibility_threshold(6.0, 21.5, None, 2, 90.0, false, 1013.25, 50.0, 0.0)
                .is_visible
        );
        // Margin sign agrees with the boolean, exactly.
        let edge = calc_visibility_threshold(6.0, 21.5, None, 2, 90.0, false, 1013.25, 50.0, 0.0);
        assert_eq!(edge.is_visible, edge.visibility_margin > 0.0);
        // Extinction dims low objects out of visibility.
        let high = calc_visibility_threshold(2.0, 21.5, None, 2, 60.0, true, 1013.25, 50.0, 0.0);
        let low = calc_visibility_threshold(2.0, 21.5, None, 2, 2.0, true, 1013.25, 50.0, 0.0);
        assert!(high.is_visible && !low.is_visible);
        assert!(low.object_magnitude > high.object_magnitude);
    }

    #[test]
    fn visibility_helpers_agree() {
        assert!(is_object_visible(-1.5, 21.5, 2, 90.0, false));
        assert!(!is_object_visible(6.0, 8.0, 2, 90.0, false));
        // Limiting magnitude for a dark sky, average observer.
        let lim = calc_limiting_magnitude_for_sky(21.5, 2, None);
        assert!((lim - 6.0).abs() < 1e-9, "lim {lim}");
        let expert = calc_limiting_magnitude_for_sky(21.5, 4, None);
        assert!(expert > lim);
        // Twilight skies limit visibility.
        let twilight = calc_limiting_magnitude_for_sky(15.0, 2, None);
        assert!((3.0..4.0).contains(&twilight), "twilight {twilight}");
        // Explicit adaptation selects the branch: mesopic at B=18.
        let mesopic = calc_limiting_magnitude_for_sky(18.0, 2, Some(EyeAdaptation::Mesopic));
        assert!((mesopic - 4.4).abs() < 1e-9, "mesopic {mesopic}");
    }

    #[test]
    fn unknown_altitude_does_not_become_one_airmass() {
        assert!(calc_airmass(f64::NAN, AirmassMethod::KastenYoung).is_nan());
        assert_eq!(calc_airmass(90.0, AirmassMethod::KastenYoung), 1.0);
        assert_eq!(
            calc_airmass(0.0, AirmassMethod::KastenYoung),
            HORIZON_AIRMASS_CAP
        );
    }
}
