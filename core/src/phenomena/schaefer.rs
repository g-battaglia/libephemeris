//! Standalone atmospheric and naked-eye visibility convenience relations.
//!
//! Inputs are caller-supplied scalar words; no ephemeris, clock, observer
//! registry or I/O is consulted. Extinction is measured in magnitudes per
//! airmass, sky brightness in nanoLamberts (nL), angular inputs in degrees,
//! altitude in metres and pressure in mbar/hPa. Magnitudes are logarithmic.
//! The Rayleigh and airmass relations have published observational origins;
//! the aerosol humidity estimate, twilight/moon curves, linear limiting
//! magnitude and magnitude-band visibility tables are project heuristics.
//! They are not the full retinal contrast model or the separate VISLIMIT
//! scattered-moonlight model.
//!
//! Physical input ranges are documented but are not imposed as new guards.
//! Ordered comparisons, early returns and clamps retain their literal
//! semantics for negative values, NaNs and infinities. Finite exponential or
//! power overflow and infinite trigonometric arguments produce typed errors;
//! ordinary floating-point arithmetic may still return infinity or NaN.
//! Functions take all arguments explicitly. The accompanying source mapping
//! records omitted-argument defaults for a compatibility entry to materialize.

use std::fmt;

/// V-band molecular extinction at sea level, in mag/airmass (Hayes & Latham).
pub const RAYLEIGH_COEFF_SEA_LEVEL: f64 = 0.1451;
/// Molecular column-density scale height in metres.
pub const RAYLEIGH_SCALE_HEIGHT: f64 = 8500.0;
/// Project label for very clear aerosol extinction, in mag/airmass.
pub const AEROSOL_COEFF_CLEAR: f64 = 0.05;
/// Project label for average aerosol extinction, in mag/airmass.
pub const AEROSOL_COEFF_AVERAGE: f64 = 0.12;
/// Project label for poor aerosol extinction, in mag/airmass.
pub const AEROSOL_COEFF_POOR: f64 = 0.25;
/// Representative visual ozone absorption in mag/airmass.
pub const OZONE_COEFF: f64 = 0.016;
/// Water-vapour coefficient multiplier in mag/airmass per humidity fraction.
pub const WATER_VAPOR_COEFF: f64 = 0.014;
/// Natural-night airglow brightness in nanoLamberts.
pub const AIRGLOW_BRIGHTNESS_NL: f64 = 145.0;
/// Approximate airglow surface-brightness label in mag/arcsec².
pub const AIRGLOW_MAG_ARCSEC2: f64 = 22.0;
/// Representative zodiacal pole brightness label in nL; not a formula clamp.
pub const ZODIACAL_MIN_NL: f64 = 40.0;
/// Representative zodiacal plane brightness label in nL; not a formula clamp.
pub const ZODIACAL_MAX_NL: f64 = 500.0;
/// Representative full-Moon sky-brightness label in nL; not the Moon formula.
pub const FULL_MOON_SKY_BRIGHTNESS_NL: f64 = 20000.0;
/// Project dark-sky transition label, in solar-depression degrees.
pub const TWILIGHT_TRANSITION_ZENITH: f64 = 16.0;
/// Representative dark-adapted pupil diameter in millimetres.
pub const DARK_ADAPTED_PUPIL_MM: f64 = 7.0;
/// Representative light-adapted pupil diameter in millimetres.
pub const LIGHT_ADAPTED_PUPIL_MM: f64 = 2.5;
/// Ideal young-observer intrinsic naked-eye magnitude threshold.
pub const IDEAL_NAKED_EYE_LIMIT_MAG: f64 = 6.5;
/// Project limiting-magnitude loss per decade past the age baseline.
pub const AGE_DEGRADATION_PER_DECADE: f64 = 0.1;
/// Observer-age baseline in years.
pub const AGE_BASELINE: f64 = 20.0;
/// Project magnitude-band arcus defaults `(inclusive_min, exclusive_max, deg)`.
/// The legacy name does not imply an ancient measured or published table.
pub const PTOLEMY_ARCUS_VISIONIS: [(f64, f64, f64); 6] = [
    (-5.0, -2.0, 5.0),
    (-2.0, 0.0, 7.0),
    (0.0, 1.5, 9.0),
    (1.5, 3.0, 11.0),
    (3.0, 4.5, 13.0),
    (4.5, 6.0, 15.0),
];
/// Project minimum elongation for intrinsic magnitude below −2, in degrees.
pub const MIN_ELONGATION_BRIGHT: f64 = 7.0;
/// Project minimum elongation for intrinsic magnitude from −2 to below 2.
pub const MIN_ELONGATION_NORMAL: f64 = 10.0;
/// Project minimum elongation for intrinsic magnitude at least 2, in degrees.
pub const MIN_ELONGATION_FAINT: f64 = 15.0;

/// Observable mathematical exception categories of these scalar relations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SchaeferError {
    /// Undefined real mathematical argument; compatible with `ValueError`.
    MathDomain,
    /// Finite exponential overflow; compatible with `OverflowError`.
    MathRange,
    /// Finite power overflow, with the host ERANGE description in its message.
    /// Compatible with errno-based `OverflowError` from scalar `**` on macOS/Unix.
    PowerRange,
    /// A selected range formula consumed a zero scalar denominator.
    /// Ordinary range selection skips zero; explicit acquired comparisons may
    /// select this branch. The binding restores actual `ZeroDivisionError` args.
    DivisionByZero,
}

impl fmt::Display for SchaeferError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MathDomain => f.write_str("math domain error"),
            Self::MathRange => f.write_str("math range error"),
            Self::DivisionByZero => f.write_str("division by zero"),
            Self::PowerRange => {
                // Float power raises an errno-based exception. ERANGE is 34
                // on the supported macOS/Unix hosts; its description comes
                // from the host OS, including its actual locale.
                let error = std::io::Error::from_raw_os_error(34).to_string();
                let description = error.strip_suffix(" (os error 34)").unwrap_or(&error);
                write!(f, "(34, '{description}')")
            }
        }
    }
}

impl std::error::Error for SchaeferError {}

/// Result carrying the precise mathematical refusal category.
pub type SchaeferResult<T> = Result<T, SchaeferError>;

/// Exponential with finite-argument overflow distinguished from `exp(+∞)`.
fn checked_exp(value: f64) -> SchaeferResult<f64> {
    let result = value.exp();
    if value.is_finite() && result.is_infinite() {
        Err(SchaeferError::MathRange)
    } else {
        Ok(result)
    }
}

/// Real power with finite-input overflow distinguished from infinite inputs.
/// Bases are positive except for the humidity square, whose exponent is the
/// integer two. A negative base with a noninteger exponent is not required.
fn checked_pow(base: f64, exponent: f64) -> SchaeferResult<f64> {
    let result = base.powf(exponent);
    if base.is_finite() && exponent.is_finite() && result.is_infinite() {
        Err(SchaeferError::PowerRange)
    } else {
        Ok(result)
    }
}

/// Real sine with an infinite-argument domain refusal; NaN remains NaN.
fn checked_sin(value: f64) -> SchaeferResult<f64> {
    if value.is_infinite() {
        Err(SchaeferError::MathDomain)
    } else {
        Ok(value.sin())
    }
}

/// Real cosine with an infinite-argument domain refusal; NaN remains NaN.
fn checked_cos(value: f64) -> SchaeferResult<f64> {
    if value.is_infinite() {
        Err(SchaeferError::MathDomain)
    } else {
        Ok(value.cos())
    }
}

/// Decimal logarithm with its nonpositive-argument mathematical refusal.
fn checked_log10(value: f64) -> SchaeferResult<f64> {
    if value <= 0.0 {
        Err(SchaeferError::MathDomain)
    } else {
        Ok(value.log10())
    }
}

/// Ordered two-argument maximum, retaining the first operand on NaN/ties.
fn ordered_max(first: f64, second: f64) -> f64 {
    if second > first { second } else { first }
}

/// Ordered two-argument minimum, retaining the first operand on NaN/ties.
fn ordered_min(first: f64, second: f64) -> f64 {
    if second < first { second } else { first }
}

/// Molecular V-band extinction in mag/airmass.
///
/// For pressure `p` in mbar/hPa and site altitude `h` in metres, computes
/// `0.1451 * (p / 1013.25) * exp(-h / 8500)` in that order. Physical pressure
/// is nonnegative; no pressure or altitude bounds are enforced. Negative
/// pressure can produce negative extinction. NaNs and infinite inputs follow
/// arithmetic semantics; finite exponential overflow returns `MathRange`.
/// The column-density model does not account for local vertical profiles.
pub fn calc_rayleigh_extinction(pressure_mbar: f64, altitude_m: f64) -> SchaeferResult<f64> {
    let pressure_factor = pressure_mbar / 1013.25;
    let altitude_factor = checked_exp(-altitude_m / RAYLEIGH_SCALE_HEIGHT)?;
    Ok(RAYLEIGH_COEFF_SEA_LEVEL * pressure_factor * altitude_factor)
}

/// Aerosol extinction in mag/airmass with a dual-use visibility parameter.
///
/// Humidity `q` is normally a fraction in `[0,1]`; `met_range_km` defaults to
/// zero and altitude `h` in metres defaults to zero. At range `V >= 1`, use
/// `max(0, 3.912/V - 0.106)`; at `0 < V < 1`, return `V` as an already local
/// coefficient without altitude scaling. Otherwise use the project estimate
/// `0.04 + 0.18*q + 0.12*q**2`. Derived coefficients multiply `exp(-h/1500)`.
/// Inputs are not clamped to physical domains. Finite square overflow returns
/// `PowerRange`; finite exponential overflow returns `MathRange`, in acquisition
/// order. This is a representative atmosphere, not a fog/profile solver.
pub fn calc_aerosol_extinction(
    humidity_fraction: f64,
    met_range_km: f64,
    altitude_m: f64,
) -> SchaeferResult<f64> {
    let selection = if met_range_km >= 1.0 {
        AerosolSelection::Range
    } else if 0.0 < met_range_km && met_range_km < 1.0 {
        AerosolSelection::Direct
    } else {
        AerosolSelection::Humidity
    };
    calc_aerosol_extinction_selected(humidity_fraction, met_range_km, altitude_m, selection)
}

/// Already acquired branch of the dual-use aerosol argument.
/// A caller may retain comparison results independently of its scalar value;
/// numerical evaluation does not repeat those comparisons.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AerosolSelection {
    /// Meteorological range formula followed by altitude scaling.
    Range,
    /// An already local extinction coefficient; all other inputs are unused.
    Direct,
    /// Humidity-derived coefficient followed by altitude scaling.
    Humidity,
}

/// Aerosol extinction from an acquired range/coefficient/humidity selection.
///
/// Units, formulas, order and approximation limits match
/// [`calc_aerosol_extinction`]. Selection is consumed, including when
/// an external comparison selected a branch different from the scalar's usual
/// interval. Direct returns the coefficient without reading humidity/altitude.
/// Humidity power errors precede altitude exponential errors. A zero denominator
/// in the selected range branch returns [`SchaeferError::DivisionByZero`].
pub fn calc_aerosol_extinction_selected(
    humidity_fraction: f64,
    met_range_km: f64,
    altitude_m: f64,
    selection: AerosolSelection,
) -> SchaeferResult<f64> {
    let k_aerosol = match selection {
        AerosolSelection::Range => {
            if met_range_km == 0.0 {
                return Err(SchaeferError::DivisionByZero);
            }
            ordered_max(0.0, 3.912 / met_range_km - 0.106)
        }
        AerosolSelection::Direct => return Ok(met_range_km),
        AerosolSelection::Humidity => {
            0.04 + 0.18 * humidity_fraction + 0.12 * checked_pow(humidity_fraction, 2.0)?
        }
    };
    let altitude_factor = checked_exp(-altitude_m / 1500.0)?;
    Ok(k_aerosol * altitude_factor)
}

/// Constant representative V-band ozone extinction, `0.016` mag/airmass.
/// Has no inputs or errors; it does not model season or ozone column changes.
pub fn calc_ozone_extinction() -> f64 {
    OZONE_COEFF
}

/// Water-vapour extinction, `0.014 * humidity_fraction`, in mag/airmass.
/// Humidity is normally `[0,1]` but is not clamped. Arithmetic NaN/infinity
/// propagate without a mathematical exception. This is a minor V-band estimate.
pub fn calc_water_vapor_extinction(humidity_fraction: f64) -> f64 {
    WATER_VAPOR_COEFF * humidity_fraction
}

/// Sum Rayleigh, aerosol, ozone and water-vapour extinction in that order.
///
/// Inputs/defaults are pressure `1013.25` mbar/hPa, temperature `15` °C,
/// humidity fraction `0.5`, dual-use range `0`, altitude `0` metres.
/// Temperature is accepted but unused, including non-finite temperatures.
/// Returns mag/airmass; Rayleigh errors precede aerosol errors. Domains,
/// NaNs and mathematical overflow follow the component functions; no
/// independent validation or atmosphere-profile correction is added.
pub fn calc_total_extinction(
    pressure_mbar: f64,
    temperature_c: f64,
    humidity_fraction: f64,
    met_range_km: f64,
    altitude_m: f64,
) -> SchaeferResult<f64> {
    let _ = temperature_c;
    let k_rayleigh = calc_rayleigh_extinction(pressure_mbar, altitude_m)?;
    let k_aerosol = calc_aerosol_extinction(humidity_fraction, met_range_km, altitude_m)?;
    let k_ozone = calc_ozone_extinction();
    let k_water = calc_water_vapor_extinction(humidity_fraction);
    Ok(calc_total_extinction_from_components(
        k_rayleigh, k_aerosol, k_ozone, k_water,
    ))
}

/// Left-associated extinction sum from four already acquired components.
/// Components and result are mag/airmass. Rayleigh, aerosol, ozone and water
/// retain this arithmetic order, including IEEE NaNs/infinities and signed zero.
/// This function reads no atmosphere, repeats no selection and adds no guards.
pub fn calc_total_extinction_from_components(
    rayleigh: f64,
    aerosol: f64,
    ozone: f64,
    water: f64,
) -> f64 {
    rayleigh + aerosol + ozone + water
}

/// Kasten & Young (1989) relative atmospheric path length.
///
/// Altitude `a <= 0` degrees returns `40`; `a >= 90` returns `1`, including
/// infinities. Otherwise set zenith angle `z = 90-a` and compute
/// `X = 1 / (cos(radians(z)) + 0.50572*(96.07995-z)**(-1.6364))`, then ordered
/// `min(max(X,1),40)`. The open `(0,90)` domain has a positive denominator;
/// its mathematical-domain/division fallback is unreachable. NaN remains
/// NaN through the ordered clamp. No refraction or layered atmosphere is used.
pub fn calc_airmass(altitude_deg: f64) -> f64 {
    if altitude_deg <= 0.0 {
        return 40.0;
    }
    if altitude_deg >= 90.0 {
        return 1.0;
    }
    let z = 90.0 - altitude_deg;
    let cos_z = z.to_radians().cos();
    let denominator = cos_z + 0.50572 * (96.07995 - z).powf(-1.6364);
    let airmass = 1.0 / denominator;
    ordered_min(ordered_max(airmass, 1.0), 40.0)
}

/// Full-path magnitude dimming `k_total * calc_airmass(altitude_deg)`.
/// Altitude is degrees and `k_total` is mag/airmass. Returns magnitudes;
/// domains and NaNs follow airmass/arithmetic, without new errors or guards.
pub fn calc_extinction_magnitude(altitude_deg: f64, k_total: f64) -> f64 {
    let airmass = calc_airmass(altitude_deg);
    k_total * airmass
}

/// Project three-phase twilight brightness in nL from solar altitude degrees.
///
/// Day (`a >= 0`) returns `1e10`; night (`a <= -18`) returns zero. For
/// depression `h=-a`, civil `h<=6` uses `1e8*exp(-1.5*h)`; nautical `h<=12`
/// log-linearly bridges the civil anchor at six to `1e3` at twelve; otherwise
/// use `1e3*exp(-0.4*(h-12))`. The finite twilight domain cannot overflow;
/// NaN propagates. The zero/night and day transitions remain discontinuous.
/// These fitted curves do not model azimuth, aerosols or scattering geometry.
pub fn calc_twilight_brightness(sun_altitude_deg: f64) -> f64 {
    if sun_altitude_deg >= 0.0 {
        return 1e10;
    }
    if sun_altitude_deg <= -18.0 {
        return 0.0;
    }
    let h = -sun_altitude_deg;
    let b_civil_at_6 = 1e8 * (-6.0_f64 * 1.5).exp();
    let b_astro_at_12 = 1e3;
    if h <= 6.0 {
        1e8 * (-h * 1.5).exp()
    } else if h <= 12.0 {
        let frac = (h - 6.0) / 6.0;
        b_civil_at_6 * (b_astro_at_12 / b_civil_at_6).powf(frac)
    } else {
        let h_eff = h - 12.0;
        b_astro_at_12 * (-h_eff * 0.4).exp()
    }
}

/// Simplified project moonlight brightness in nL; all angular inputs degrees.
///
/// Moon altitude at/below zero returns zero before inspecting other inputs.
/// Distance below one is raised to one. Phase is `0` full, `180` new; use
/// `f=(1+cos(radians(phase)))/2`, returning zero when `f<0.01`. Otherwise
/// multiply `3000*f*exp(-0.4*(X-1))*d`, with Moon airmass `X` and distance
/// factor `d=10/(rho+0.1)` for `rho<10`, `1+(45-rho)/35` for `rho<45`, else
/// one. Infinite phase returns `MathDomain` unless an early return applies;
/// NaNs retain branch/arithmetic semantics. This is a heuristic, not the
/// Krisciunas–Schaefer scattered-light equations used by VISLIMIT.
pub fn calc_moon_brightness(
    moon_altitude_deg: f64,
    moon_phase: f64,
    angular_distance_deg: f64,
) -> SchaeferResult<f64> {
    if moon_altitude_deg <= 0.0 {
        return Ok(0.0);
    }
    let angular_distance_deg = if angular_distance_deg < 1.0 {
        1.0
    } else {
        angular_distance_deg
    };
    let phase_rad = moon_phase.to_radians();
    let illuminated_fraction = (1.0 + checked_cos(phase_rad)?) / 2.0;
    if illuminated_fraction < 0.01 {
        return Ok(0.0);
    }
    let moon_airmass = calc_airmass(moon_altitude_deg);
    let alt_factor = (-0.4 * (moon_airmass - 1.0)).exp();
    let dist_factor = if angular_distance_deg < 10.0 {
        10.0 / (angular_distance_deg + 0.1)
    } else if angular_distance_deg < 45.0 {
        1.0 + (45.0 - angular_distance_deg) / 35.0
    } else {
        1.0
    };
    Ok(3000.0 * illuminated_fraction * alt_factor * dist_factor)
}

/// Constant natural-night airglow of `145` nL, with no inputs or errors.
/// Seasonal, latitude and time variability are outside this convenience model.
pub fn calc_airglow_brightness() -> f64 {
    AIRGLOW_BRIGHTNESS_NL
}

/// Project zodiacal sky brightness, `150*exp(-abs(beta_rad))*e`, in nL.
///
/// Ecliptic latitude `beta` and solar elongation `epsilon` are degrees.
/// `e=0` below `20`, `sin(epsilon)` below `90`, otherwise `sin(180-epsilon)`.
/// No latitude or elongation bounds are imposed; beyond 180 the result can
/// be negative. Infinite elongation in a sine branch returns `MathDomain`;
/// NaNs propagate through the applicable calculations. Latitude is evaluated
/// even when the elongation factor is zero. No dust distribution is solved.
pub fn calc_zodiacal_brightness(
    ecliptic_latitude_deg: f64,
    sun_elongation_deg: f64,
) -> SchaeferResult<f64> {
    let lat_rad = ecliptic_latitude_deg.abs().to_radians();
    let lat_factor = (-lat_rad).exp();
    let elong_factor = if sun_elongation_deg < 20.0 {
        0.0
    } else if sun_elongation_deg < 90.0 {
        checked_sin(sun_elongation_deg.to_radians())?
    } else {
        checked_sin((180.0 - sun_elongation_deg).to_radians())?
    };
    Ok(150.0 * lat_factor * elong_factor)
}

/// Sum twilight, moonlight, airglow and ramped zodiacal brightness in nL.
///
/// All inputs are degrees. Only Sun altitude is required by the convenience
/// entry; defaults are Moon altitude `-90`, phase `180`, distance `90`,
/// ecliptic latitude `0`, solar elongation `90`. Zodiacal weight is one at
/// Sun `<=-18`, `(-12-Sun)/6` below `-12`, otherwise zero. Its component is
/// still evaluated at zero weight, retaining domain errors/NaNs. Component
/// errors propagate in evaluation order. No geometry or ephemeris is acquired.
#[allow(clippy::too_many_arguments)]
pub fn calc_total_sky_brightness(
    sun_altitude_deg: f64,
    moon_altitude_deg: f64,
    moon_phase_deg: f64,
    moon_distance_deg: f64,
    ecliptic_latitude_deg: f64,
    sun_elongation_deg: f64,
) -> SchaeferResult<f64> {
    let b_twilight = calc_twilight_brightness(sun_altitude_deg);
    let b_moon = calc_moon_brightness(moon_altitude_deg, moon_phase_deg, moon_distance_deg)?;
    let b_airglow = calc_airglow_brightness();
    let zodiacal_fraction = if sun_altitude_deg <= -18.0 {
        1.0
    } else if sun_altitude_deg < -12.0 {
        (-12.0 - sun_altitude_deg) / 6.0
    } else {
        0.0
    };
    let b_zodiacal =
        zodiacal_fraction * calc_zodiacal_brightness(ecliptic_latitude_deg, sun_elongation_deg)?;
    Ok(b_twilight + b_moon + b_airglow + b_zodiacal)
}

/// Surface brightness `26.33 - 2.5*log10(B)` in mag/arcsec² for nL brightness.
/// Nonpositive `B` returns `22`, preserving a dark-sky sentinel instead of
/// taking its logarithm. Positive infinity returns negative infinity; NaN
/// remains NaN. The guarded logarithm has no domain errors. The zero point
/// follows the Garstang relation; this does not perform spectral conversion.
pub fn brightness_to_mag_arcsec2(brightness_nl: f64) -> f64 {
    if brightness_nl <= 0.0 {
        return 22.0;
    }
    26.33 - 2.5 * brightness_nl.log10()
}

/// nL brightness `10**((26.33-m)/2.5)` from surface magnitude `m` in mag/arcsec².
/// This inverts the positive-brightness branch, not its nonpositive sentinel.
/// No bounds are enforced. Finite power overflow returns `PowerRange`;
/// infinities/NaNs propagate, and very faint finite magnitudes may underflow.
pub fn mag_arcsec2_to_brightness(mag_arcsec2: f64) -> SchaeferResult<f64> {
    checked_pow(10.0, (26.33 - mag_arcsec2) / 2.5)
}

/// Intrinsic naked-eye magnitude threshold from a linear convenience model.
///
/// Sky `B` is nL, extinction `k` mag/airmass, object altitude degrees, age
/// years (default `25`) and Snellen ratio dimensionless (default `1`). Start
/// at `6.5`; subtract `2.5*log10(B/250)` only for `B>250`, age loss
/// `0.1*(age-20)/10` only past twenty, `2.5*log10(1/Snellen)` only for positive
/// nonunit Snellen, and `k*(airmass-1)`, in that order. Nonpositive Snellen
/// takes zero loss; no physical validation is imposed. An overflowing
/// reciprocal can become infinity without an exception. Infinite positive
/// Snellen instead produces a zero reciprocal and `MathDomain`. This is not the
/// full Schaefer retinal contrast-threshold model. Extinction is already
/// incorporated into this intrinsic threshold.
pub fn calc_limiting_magnitude(
    sky_brightness_nl: f64,
    extinction_coeff: f64,
    object_altitude_deg: f64,
    observer_age: f64,
    snellen_ratio: f64,
) -> SchaeferResult<f64> {
    let sky_factor = if sky_brightness_nl > 250.0 {
        2.5 * (sky_brightness_nl / 250.0).log10()
    } else {
        0.0
    };
    let age_factor = if observer_age > AGE_BASELINE {
        let decades_past_baseline = (observer_age - AGE_BASELINE) / 10.0;
        AGE_DEGRADATION_PER_DECADE * decades_past_baseline
    } else {
        0.0
    };
    let snellen_factor = if snellen_ratio > 0.0 && snellen_ratio != 1.0 {
        2.5 * checked_log10(1.0 / snellen_ratio)?
    } else {
        0.0
    };
    let airmass = calc_airmass(object_altitude_deg);
    let extinction_loss = extinction_coeff * (airmass - 1.0);
    Ok(IDEAL_NAKED_EYE_LIMIT_MAG - sky_factor - age_factor - snellen_factor - extinction_loss)
}

/// Visibility tuple `(visible, intrinsic_limit_mag, apparent_object_mag)`.
///
/// Object/Sun/Moon altitude, phase and solar elongation use degrees; object
/// magnitude is intrinsic. Defaults: pressure `1013.25` mbar, temperature
/// `15` °C (unused), humidity `0.5`, range `0`, site altitude `0` metres,
/// age `25` years, Snellen `1`, Moon altitude `-90`, phase `180`, elongation
/// `90`, include Moon `true`. Object at/below horizon returns `(false,0,m)`
/// first; Sun at/above horizon returns `(false,-5,m)` next. Otherwise sum
/// extinction/sky at fixed Moon distance `90` and ecliptic latitude `0`,
/// compare intrinsic magnitude with its already-extincted threshold, and
/// apply solar-elongation minima `7/10/15` for magnitude `<-2/<2/otherwise`.
/// The apparent magnitude `m+k*airmass` is informational. Component overflow
/// and domain errors propagate unless an early return skips acquisition.
/// This convenience geometry does not predict rise/set or visibility events.
#[allow(clippy::too_many_arguments)]
pub fn is_object_visible(
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
) -> SchaeferResult<(bool, f64, f64)> {
    if object_altitude_deg <= 0.0 {
        return Ok((false, 0.0, object_magnitude));
    }
    if sun_altitude_deg >= 0.0 {
        return Ok((false, -5.0, object_magnitude));
    }
    let k_total = calc_total_extinction(
        pressure_mbar,
        temperature_c,
        humidity_fraction,
        met_range_km,
        altitude_m,
    )?;
    let sky_brightness = if include_moon {
        calc_total_sky_brightness(
            sun_altitude_deg,
            moon_altitude_deg,
            moon_phase_deg,
            90.0,
            0.0,
            sun_elongation_deg,
        )?
    } else {
        calc_total_sky_brightness(
            sun_altitude_deg,
            -90.0,
            180.0,
            90.0,
            0.0,
            sun_elongation_deg,
        )?
    };
    let limiting_mag = calc_limiting_magnitude(
        sky_brightness,
        k_total,
        object_altitude_deg,
        observer_age,
        snellen_ratio,
    )?;
    Ok(is_object_visible_from_limit(
        object_magnitude,
        object_altitude_deg,
        sun_elongation_deg,
        k_total,
        limiting_mag,
    ))
}

/// Visibility and apparent magnitude from an acquired extinction and limit.
///
/// Object altitude and solar elongation are degrees, extinction mag/airmass
/// and both magnitudes intrinsic/apparent Pogson magnitudes. The same airmass,
/// informational apparent magnitude, intrinsic threshold comparison and
/// seven/ten/fifteen-degree elongation gates of [`is_object_visible`] apply in
/// that order. This tail repeats no atmosphere or sky acquisition. Early
/// horizon/daylight sentinels belong to the full entry, before this tail.
pub fn is_object_visible_from_limit(
    object_magnitude: f64,
    object_altitude_deg: f64,
    sun_elongation_deg: f64,
    k_total: f64,
    limiting_mag: f64,
) -> (bool, f64, f64) {
    let airmass = calc_airmass(object_altitude_deg);
    let apparent_mag = object_magnitude + k_total * airmass;
    let mut visible = object_magnitude <= limiting_mag;
    let min_elong = if object_magnitude < -2.0 {
        MIN_ELONGATION_BRIGHT
    } else if object_magnitude < 2.0 {
        MIN_ELONGATION_NORMAL
    } else {
        MIN_ELONGATION_FAINT
    };
    if sun_elongation_deg < min_elong {
        visible = false;
    }
    (visible, limiting_mag, apparent_mag)
}

/// Project arcus-visionis default in degrees for intrinsic magnitude.
/// Uses inclusive-lower/exclusive-upper table bins, then `17` for `m>=6`,
/// otherwise `4`. NaN follows the final `4` branch; infinities retain their
/// ordered comparisons. No errors or observational-history claim are implied.
pub fn get_arcus_visionis(object_magnitude: f64) -> f64 {
    for (mag_min, mag_max, arcus) in PTOLEMY_ARCUS_VISIONIS {
        if mag_min <= object_magnitude && object_magnitude < mag_max {
            return arcus;
        }
    }
    if object_magnitude >= 6.0 { 17.0 } else { 4.0 }
}

/// Project optimal solar altitude in degrees from intrinsic magnitude bins.
/// Returns `-8/-9/-10/-11/-12` for `m<-2/<0/<2/<4/otherwise`. The required
/// `is_morning` argument is intentionally unused: both event directions use
/// the same table. NaN takes `-12`; there are no errors, epochs or event solves.
pub fn get_optimal_sun_altitude(object_magnitude: f64, is_morning: bool) -> f64 {
    let _ = is_morning;
    if object_magnitude < -2.0 {
        -8.0
    } else if object_magnitude < 0.0 {
        -9.0
    } else if object_magnitude < 2.0 {
        -10.0
    } else if object_magnitude < 4.0 {
        -11.0
    } else {
        -12.0
    }
}

/// Project minimum heliacal object altitude in degrees.
/// Multiply the magnitude-band arcus by extinction factor `1+0.5*(k-0.25)`
/// (`k` mag/airmass) and Sun-altitude factor `0.8` below `-12`, `1` below
/// `-9`, otherwise `1.2+0.1*(Sun+9)` with Sun in degrees. Ordered
/// `max(1,result)` floors at one and returns one for NaN. Negative/unphysical
/// inputs are not independently rejected; arithmetic infinity may propagate.
/// This is a heuristic threshold, with no errors or heliacal event solver.
pub fn calc_heliacal_altitude_threshold(
    object_magnitude: f64,
    sun_altitude_deg: f64,
    extinction_coeff: f64,
) -> f64 {
    let base_av = get_arcus_visionis(object_magnitude);
    calc_heliacal_altitude_threshold_with_arcus(base_av, sun_altitude_deg, extinction_coeff)
}

/// Minimum object altitude in degrees from an already acquired arcus angle.
/// The input arcus is degrees and extinction mag/airmass. Extinction factor,
/// solar-altitude bins, multiplication order and ordered one-degree floor match
/// [`calc_heliacal_altitude_threshold`]. No table lookup, physical validation or
/// astronomical provider is repeated; IEEE behavior remains unchanged.
pub fn calc_heliacal_altitude_threshold_with_arcus(
    base_av: f64,
    sun_altitude_deg: f64,
    extinction_coeff: f64,
) -> f64 {
    let extinction_factor = 1.0 + 0.5 * (extinction_coeff - 0.25);
    let sun_factor = if sun_altitude_deg < -12.0 {
        0.8
    } else if sun_altitude_deg < -9.0 {
        1.0
    } else {
        1.2 + (sun_altitude_deg + 9.0) * 0.1
    };
    let min_altitude = base_av * extinction_factor * sun_factor;
    ordered_max(1.0, min_altitude)
}

/// The thirteen visibility-dictionary fields in their public insertion order.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct VisibilityConditions {
    /// Summed extinction in mag/airmass.
    pub extinction_total: f64,
    /// Molecular extinction in mag/airmass.
    pub extinction_rayleigh: f64,
    /// Aerosol extinction in mag/airmass.
    pub extinction_aerosol: f64,
    /// Relative atmospheric path length.
    pub airmass: f64,
    /// Full-path object dimming in magnitudes.
    pub extinction_magnitude: f64,
    /// Summed default sky brightness in nL.
    pub sky_brightness_nl: f64,
    /// Sky surface brightness in mag/arcsec².
    pub sky_brightness_mag_arcsec2: f64,
    /// Faintest visible intrinsic magnitude.
    pub limiting_magnitude: f64,
    /// Informational object magnitude after full-path extinction.
    pub apparent_magnitude: f64,
    /// Intrinsic threshold minus intrinsic object magnitude.
    pub magnitude_margin: f64,
    /// Intrinsic comparison only, without horizon/day/elongation gates.
    pub is_visible: bool,
    /// Project magnitude-band arcus default in degrees.
    pub arcus_visionis: f64,
    /// Project optimal solar altitude in degrees, with morning selection.
    pub optimal_sun_altitude: f64,
}

/// Detailed atmospheric/visibility values preserving all thirteen dict keys.
///
/// Required inputs: Sun/object altitude in degrees and intrinsic magnitude.
/// Defaults: pressure `1013.25` mbar, humidity `0.5`, site altitude `0` metres.
/// Use temperature `15`, range `0`, default Moon/sky geometry and age/Snellen
/// `25/1`. Acquire total/Rayleigh/aerosol, then sky brightness/conversion,
/// airmass, full extinction, limiting and apparent magnitudes. The margin
/// and visibility compare intrinsic magnitude; unlike `is_object_visible`,
/// this diagnostic dictionary has no horizon, daytime or elongation gates.
/// Component errors and non-finite arithmetic retain their original order.
/// No caller state, astronomical acquisition or event search is consulted.
pub fn get_visibility_conditions(
    sun_altitude_deg: f64,
    object_altitude_deg: f64,
    object_magnitude: f64,
    pressure_mbar: f64,
    humidity_fraction: f64,
    altitude_m: f64,
) -> SchaeferResult<VisibilityConditions> {
    let k_total = calc_total_extinction(pressure_mbar, 15.0, humidity_fraction, 0.0, altitude_m)?;
    let k_rayleigh = calc_rayleigh_extinction(pressure_mbar, altitude_m)?;
    let k_aerosol = calc_aerosol_extinction(humidity_fraction, 0.0, altitude_m)?;
    let sky_brightness =
        calc_total_sky_brightness(sun_altitude_deg, -90.0, 180.0, 90.0, 0.0, 90.0)?;
    let sky_mag = brightness_to_mag_arcsec2(sky_brightness);
    let airmass = calc_airmass(object_altitude_deg);
    let extinction_mag = k_total * airmass;
    let limiting_mag =
        calc_limiting_magnitude(sky_brightness, k_total, object_altitude_deg, 25.0, 1.0)?;
    let apparent_mag = object_magnitude + extinction_mag;
    Ok(VisibilityConditions {
        extinction_total: k_total,
        extinction_rayleigh: k_rayleigh,
        extinction_aerosol: k_aerosol,
        airmass,
        extinction_magnitude: extinction_mag,
        sky_brightness_nl: sky_brightness,
        sky_brightness_mag_arcsec2: sky_mag,
        limiting_magnitude: limiting_mag,
        apparent_magnitude: apparent_mag,
        magnitude_margin: limiting_mag - object_magnitude,
        is_visible: object_magnitude <= limiting_mag,
        arcus_visionis: get_arcus_visionis(object_magnitude),
        optimal_sun_altitude: get_optimal_sun_altitude(object_magnitude, true),
    })
}

#[cfg(test)]
#[path = "tests/schaefer.rs"]
mod tests;
