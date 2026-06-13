# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Atmospheric extinction model for libephemeris.

This module provides functions to calculate atmospheric extinction,
which increases the apparent magnitude of celestial objects as they
approach the horizon. This is essential for heliacal visibility calculations.

The extinction model is based on the work of:
- Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
- Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"

Key concepts:
- Extinction is measured in magnitudes per airmass (mag/airmass)
- Airmass = sec(z) where z is zenith angle (for plane-parallel atmosphere)
- Total extinction = extinction_coefficient * airmass

The extinction coefficient k is composed of:
- Rayleigh scattering (molecular): wavelength-dependent, ~0.14 mag/airmass at V
- Aerosol scattering: varies with atmospheric conditions, typically 0.05-0.25
- Ozone absorption: small contribution, ~0.016 mag/airmass

References:
    - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
    - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    - Young, A.T. (1974) "Observational Technique and Data Reduction"
    - Kasten, F. & Young, A.T. (1989) "Revised optical air mass tables"
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


# Standard wavelengths in nanometers for common photometric bands
WAVELENGTH_U = 365.0  # U band (ultraviolet)
WAVELENGTH_B = 440.0  # B band (blue)
WAVELENGTH_V = 550.0  # V band (visual, standard)
WAVELENGTH_R = 640.0  # R band (red)
WAVELENGTH_I = 790.0  # I band (infrared)

# Default atmospheric parameters
DEFAULT_PRESSURE_MBAR = 1013.25  # Standard sea level pressure
DEFAULT_TEMPERATURE_C = 15.0  # Standard temperature
DEFAULT_HUMIDITY_PERCENT = 50.0  # Moderate humidity
DEFAULT_ALTITUDE_M = 0.0  # Sea level

# Atmospheric scale height in meters (for pressure/density variation)
SCALE_HEIGHT_M = 8500.0


@dataclass
class ExtinctionCoefficients:
    """
    Atmospheric extinction coefficients broken down by component.

    Attributes:
        k_rayleigh: Rayleigh (molecular) scattering coefficient [mag/airmass]
        k_aerosol: Aerosol (Mie) scattering coefficient [mag/airmass]
        k_ozone: Ozone absorption coefficient [mag/airmass]
        k_water: Water vapor absorption coefficient [mag/airmass]
        k_total: Total extinction coefficient [mag/airmass]
    """

    k_rayleigh: float
    k_aerosol: float
    k_ozone: float
    k_water: float
    k_total: float


def calc_airmass(altitude_deg: float, method: str = "kasten_young") -> float:
    """
    Calculate the airmass (relative path length through atmosphere).

    Airmass indicates how much atmosphere light must traverse compared
    to looking straight up (zenith, airmass = 1.0). At the horizon,
    airmass is approximately 38-40.

    Args:
        altitude_deg: Altitude of object above horizon in degrees.
                     Negative values indicate object below horizon.
        method: Airmass calculation method:
            - "secant": Simple sec(z) formula, accurate above ~15 degrees
            - "kasten_young": Kasten & Young (1989) formula, accurate to horizon
            - "rozenberg": Rozenberg (1966) formula for very low altitudes

    Returns:
        Airmass value (dimensionless). Returns 40.0 for objects at or below
        horizon as practical maximum.

    Algorithm:
        The simple plane-parallel approximation uses X = sec(z) where z is
        the zenith angle. This breaks down near the horizon due to Earth's
        curvature and atmospheric refraction.

        The Kasten & Young (1989) formula provides accurate airmass even
        near the horizon:
            X = 1 / [sin(h) + 0.50572 * (h + 6.07995)^(-1.6364)]
        where h is altitude in degrees.

    Example:
        >>> calc_airmass(90.0)  # Zenith
        1.0
        >>> calc_airmass(30.0)  # 30 degrees altitude
        2.0
        >>> calc_airmass(0.0)   # Horizon
        37.9...

    References:
        - Kasten, F. & Young, A.T. (1989) Applied Optics 28, 4735-4738
        - Rozenberg, G.V. (1966) Twilight: A Study in Atmospheric Optics
    """
    # Handle objects at or below horizon
    if altitude_deg <= 0:
        return 40.0  # Practical maximum airmass

    # Zenith angle in degrees and radians
    zenith_deg = 90.0 - altitude_deg
    zenith_rad = math.radians(zenith_deg)
    altitude_rad = math.radians(altitude_deg)

    if method == "secant":
        # Simple plane-parallel approximation
        # Only accurate above ~15 degrees altitude
        if zenith_deg >= 85:
            return 40.0  # Avoid division by very small cosine
        return 1.0 / math.cos(zenith_rad)

    elif method == "kasten_young":
        # Kasten & Young (1989) empirical formula
        # Accurate from zenith to horizon
        sin_h = math.sin(altitude_rad)
        # Altitude in degrees for polynomial term
        h = altitude_deg
        denominator = sin_h + 0.50572 * ((h + 6.07995) ** (-1.6364))
        if denominator <= 0:  # pragma: no cover - unreachable: a non-positive
            # denominator requires a negative power base, which raises before
            # this guard can return; kept as a defensive floor.
            return 40.0
        return min(1.0 / denominator, 40.0)

    elif method == "rozenberg":
        # Rozenberg (1966) formula, good for very low altitudes
        cos_z = math.cos(zenith_rad)
        return 1.0 / (cos_z + 0.025 * math.exp(-11.0 * cos_z))

    else:
        # Default to Kasten & Young
        return calc_airmass(altitude_deg, method="kasten_young")


def calc_rayleigh_coefficient(
    wavelength_nm: float = WAVELENGTH_V,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    altitude_m: float = DEFAULT_ALTITUDE_M,
) -> float:
    """
    Calculate the Rayleigh scattering extinction coefficient.

    Rayleigh scattering is caused by molecules in the atmosphere
    (primarily N2 and O2). It is strongly wavelength-dependent,
    following a lambda^(-4) law, which is why the sky is blue.

    Args:
        wavelength_nm: Wavelength of observation in nanometers.
                       Default is V-band (550 nm).
        pressure_mbar: Atmospheric pressure in millibars (hPa).
        temperature_c: Temperature in degrees Celsius.
        altitude_m: Observer altitude in meters above sea level.

    Returns:
        Rayleigh scattering coefficient in magnitudes per airmass.

    Algorithm:
        k_R = 0.1451 * (P/1013.25) * (lambda_0/lambda)^4
        where lambda_0 = 550 nm (V-band reference)

        The coefficient 0.1451 is the standard Rayleigh coefficient
        at sea level for V-band observations.

    Example:
        >>> calc_rayleigh_coefficient(550.0, 1013.25)  # Standard conditions
        0.1451
        >>> calc_rayleigh_coefficient(440.0, 1013.25)  # B-band (bluer)
        0.359...  # Higher extinction for blue light

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Young, A.T. (1974) "Observational Technique and Data Reduction"
    """
    # Reference values
    k_ref = 0.1451  # Rayleigh coefficient at sea level, V-band
    wavelength_ref = 550.0  # V-band reference wavelength

    # Pressure correction (higher pressure = more scattering)
    pressure_factor = pressure_mbar / DEFAULT_PRESSURE_MBAR

    # Altitude correction (less atmosphere at higher altitudes)
    altitude_factor = math.exp(-altitude_m / SCALE_HEIGHT_M)

    # Wavelength dependence (Rayleigh scattering ~ lambda^-4)
    wavelength_factor = (wavelength_ref / wavelength_nm) ** 4

    return k_ref * pressure_factor * altitude_factor * wavelength_factor


def calc_aerosol_coefficient(
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
    visibility_km: Optional[float] = None,
) -> float:
    """
    Calculate the aerosol (Mie) scattering extinction coefficient.

    Aerosol scattering is caused by dust, pollen, smoke, and other
    particulates in the atmosphere. It is highly variable and depends
    on local conditions.

    Args:
        humidity_percent: Relative humidity (0-100).
                         Higher humidity increases aerosol scattering.
        altitude_m: Observer altitude in meters above sea level.
                   Aerosol concentration decreases with altitude.
        wavelength_nm: Wavelength of observation in nanometers.
        visibility_km: Meteorological visibility in km. If provided,
                       this is used to estimate aerosol content directly.

    Returns:
        Aerosol scattering coefficient in magnitudes per airmass.
        Typical values range from 0.05 (excellent) to 0.30 (hazy).

    Algorithm:
        If visibility is given:
            k_A = 3.912 / V - 0.1066 (empirical formula)
        Otherwise, estimate from humidity:
            k_A = 0.08 * (1 + 2 * (humidity/100)^2) * exp(-altitude/H_A)
        where H_A ~ 1500m is the aerosol scale height.

    Example:
        >>> calc_aerosol_coefficient(30.0, 0.0)   # Low humidity, sea level
        0.10...
        >>> calc_aerosol_coefficient(80.0, 0.0)   # High humidity
        0.18...
        >>> calc_aerosol_coefficient(50.0, 2000)  # Moderate humidity, 2km altitude
        0.07...

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    """
    if visibility_km is not None and visibility_km > 0:
        # Direct calculation from meteorological visibility
        # Using Koschmieder formula: V = 3.912 / beta
        # where beta is extinction coefficient per km
        # Convert to mag/airmass
        k_aerosol = 3.912 / visibility_km - 0.1066
        return max(0.0, k_aerosol)

    # Aerosol scale height (aerosols concentrated in lower atmosphere)
    aerosol_scale_height = 1500.0  # meters

    # Base aerosol coefficient (clean atmosphere)
    k_base = 0.08

    # Humidity factor (hygroscopic growth of aerosols)
    humidity_fraction = humidity_percent / 100.0
    humidity_factor = 1.0 + 2.0 * humidity_fraction**2

    # Altitude factor (aerosols decrease exponentially with altitude)
    altitude_factor = math.exp(-altitude_m / aerosol_scale_height)

    # Wavelength dependence (weaker than Rayleigh: ~lambda^-1.3)
    wavelength_factor = (WAVELENGTH_V / wavelength_nm) ** 1.3

    return k_base * humidity_factor * altitude_factor * wavelength_factor


def calc_ozone_coefficient(wavelength_nm: float = WAVELENGTH_V) -> float:
    """
    Calculate the ozone absorption coefficient.

    Ozone (O3) in the stratosphere absorbs UV and some visible light.
    The effect is small in the visual band but significant in UV/blue.

    Args:
        wavelength_nm: Wavelength of observation in nanometers.

    Returns:
        Ozone absorption coefficient in magnitudes per airmass.
        Typically 0.01-0.04 for visual wavelengths.

    Algorithm:
        Ozone absorption has a complex wavelength dependence.
        For visual wavelengths, we use an empirical approximation.

    Example:
        >>> calc_ozone_coefficient(550.0)  # V-band
        0.016
    """
    # Ozone absorption peaks in UV (Hartley band) and has Chappuis bands
    # in visible (around 600nm). For V-band, the contribution is small.
    if wavelength_nm < 400:
        # UV region - strong ozone absorption
        return 0.10
    elif wavelength_nm < 500:
        # Blue region
        return 0.04
    elif wavelength_nm < 700:
        # Visual region (includes Chappuis bands)
        return 0.016
    else:
        # Red/IR region
        return 0.008


def calc_water_vapor_coefficient(
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    wavelength_nm: float = WAVELENGTH_V,
) -> float:
    """
    Calculate the water vapor absorption coefficient.

    Water vapor absorption is significant primarily in the infrared
    and has relatively minor effect in the visual band.

    Args:
        humidity_percent: Relative humidity (0-100).
        temperature_c: Temperature in degrees Celsius.
        wavelength_nm: Wavelength of observation in nanometers.

    Returns:
        Water vapor absorption coefficient in magnitudes per airmass.
        Typically very small (<0.01) for visual wavelengths.

    Example:
        >>> calc_water_vapor_coefficient(50.0, 15.0, 550.0)
        0.003...
    """
    # Water vapor effect is minimal in visual bands
    # but included for completeness
    humidity_fraction = humidity_percent / 100.0

    if wavelength_nm < 600:
        # Blue-green: minimal water vapor effect
        return 0.002 * humidity_fraction
    elif wavelength_nm < 800:
        # Red: small water vapor bands
        return 0.005 * humidity_fraction
    else:
        # Near-IR: significant water vapor bands
        return 0.02 * humidity_fraction


def calc_extinction_coefficient(
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
    visibility_km: Optional[float] = None,
) -> ExtinctionCoefficients:
    """
    Calculate the total atmospheric extinction coefficient and its components.

    This is the main function for determining how much light is lost
    as it passes through the atmosphere. The extinction coefficient k
    is used with airmass X to calculate total extinction:
        delta_m = k * X
    where delta_m is the increase in apparent magnitude.

    Args:
        pressure_mbar: Atmospheric pressure in millibars (hPa).
                      Default is 1013.25 (standard sea level).
        temperature_c: Temperature in degrees Celsius.
                      Default is 15.0.
        humidity_percent: Relative humidity (0-100).
                         Default is 50.0.
        altitude_m: Observer altitude in meters above sea level.
                   Default is 0.0 (sea level).
        wavelength_nm: Wavelength of observation in nanometers.
                      Default is 550 (V-band).
        visibility_km: Meteorological visibility in km.
                      If provided, used to estimate aerosol content.

    Returns:
        ExtinctionCoefficients dataclass with individual components
        and total extinction coefficient in mag/airmass.

    Example:
        >>> coeff = calc_extinction_coefficient()
        >>> print(f"Total extinction: {coeff.k_total:.3f} mag/airmass")
        Total extinction: 0.28... mag/airmass

        >>> coeff = calc_extinction_coefficient(humidity_percent=90)
        >>> print(f"Humid extinction: {coeff.k_total:.3f} mag/airmass")
        Humid extinction: 0.36... mag/airmass

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    """
    k_rayleigh = calc_rayleigh_coefficient(
        wavelength_nm=wavelength_nm,
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        altitude_m=altitude_m,
    )

    k_aerosol = calc_aerosol_coefficient(
        humidity_percent=humidity_percent,
        altitude_m=altitude_m,
        wavelength_nm=wavelength_nm,
        visibility_km=visibility_km,
    )

    k_ozone = calc_ozone_coefficient(wavelength_nm=wavelength_nm)

    k_water = calc_water_vapor_coefficient(
        humidity_percent=humidity_percent,
        temperature_c=temperature_c,
        wavelength_nm=wavelength_nm,
    )

    k_total = k_rayleigh + k_aerosol + k_ozone + k_water

    return ExtinctionCoefficients(
        k_rayleigh=k_rayleigh,
        k_aerosol=k_aerosol,
        k_ozone=k_ozone,
        k_water=k_water,
        k_total=k_total,
    )


def calc_extinction_magnitude(
    altitude_deg: float,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    observer_altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
    visibility_km: Optional[float] = None,
) -> float:
    """
    Calculate the total atmospheric extinction in magnitudes for an object.

    This function combines the extinction coefficient with airmass to give
    the actual extinction experienced by an object at a given altitude.

    For objects near the horizon, the extinction can be several magnitudes,
    making faint objects invisible.

    Args:
        altitude_deg: Altitude of object above horizon in degrees.
        pressure_mbar: Atmospheric pressure in millibars (hPa).
        temperature_c: Temperature in degrees Celsius.
        humidity_percent: Relative humidity (0-100).
        observer_altitude_m: Observer altitude in meters above sea level.
        wavelength_nm: Wavelength of observation in nanometers.
        visibility_km: Meteorological visibility in km.

    Returns:
        Total extinction in magnitudes. This value should be added
        to the object's catalog magnitude to get apparent magnitude.
        Returns a large value (99.0) for objects below the horizon.

    Algorithm:
        extinction = k * X
        where k is the extinction coefficient and X is airmass.

        Using the basic formula (approximately):
            extinction ~ 0.28 * sec(z)
        where z is zenith angle and 0.28 is typical k for V-band.

    Example:
        >>> calc_extinction_magnitude(90.0)  # Zenith
        0.28...  # Minimal extinction
        >>> calc_extinction_magnitude(30.0)  # 30 degrees altitude
        0.56...  # 2x airmass = 2x extinction
        >>> calc_extinction_magnitude(5.0)   # Near horizon
        3.2...   # Significant extinction
        >>> calc_extinction_magnitude(0.0)   # Horizon
        11.2...  # Very high extinction

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    """
    # Objects below horizon have essentially infinite extinction
    if altitude_deg < 0:
        return 99.0

    # Calculate airmass
    airmass = calc_airmass(altitude_deg, method="kasten_young")

    # Calculate extinction coefficient
    coeff = calc_extinction_coefficient(
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        humidity_percent=humidity_percent,
        altitude_m=observer_altitude_m,
        wavelength_nm=wavelength_nm,
        visibility_km=visibility_km,
    )

    # Total extinction = k * X
    return coeff.k_total * airmass


def calc_simple_extinction(altitude_deg: float, k: float = 0.28) -> float:
    """
    Calculate atmospheric extinction using the simple Schaefer/Green formula.

    This is the basic formula mentioned in the task:
        extinction = k * sec(z) = k / cos(z) = k / sin(h)
    where z is zenith angle, h is altitude, and k ~ 0.28 for typical
    visual observations.

    Args:
        altitude_deg: Altitude of object above horizon in degrees.
        k: Extinction coefficient in mag/airmass. Default 0.28 is
           typical for V-band at a good observing site.

    Returns:
        Extinction in magnitudes.

    Note:
        This simple formula is accurate for altitudes above ~15 degrees.
        For lower altitudes, use calc_extinction_magnitude() which
        applies proper airmass corrections for atmospheric curvature.

    Example:
        >>> calc_simple_extinction(90.0, 0.28)  # Zenith
        0.28
        >>> calc_simple_extinction(30.0, 0.28)  # 30 degrees (sec(60) = 2)
        0.56
        >>> calc_simple_extinction(19.47, 0.28)  # airmass = 3
        0.84

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction"
    """
    if altitude_deg <= 0:
        return 99.0

    # Simple secant formula: sec(z) = 1/sin(h)
    sin_alt = math.sin(math.radians(altitude_deg))
    if sin_alt <= 0.001:  # Very near horizon
        return k * 40.0  # Practical maximum

    return k / sin_alt


def apparent_magnitude_with_extinction(
    catalog_magnitude: float,
    altitude_deg: float,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    observer_altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
) -> float:
    """
    Calculate the apparent magnitude of an object including extinction.

    This applies atmospheric extinction to convert from catalog (extra-
    atmospheric) magnitude to what an observer would actually see.

    Args:
        catalog_magnitude: The object's magnitude outside the atmosphere.
        altitude_deg: Altitude of object above horizon in degrees.
        pressure_mbar: Atmospheric pressure in millibars.
        temperature_c: Temperature in degrees Celsius.
        humidity_percent: Relative humidity (0-100).
        observer_altitude_m: Observer altitude in meters.
        wavelength_nm: Wavelength of observation in nanometers.

    Returns:
        Apparent magnitude as seen by observer (fainter = higher value).
        Objects below horizon return a very large magnitude (99+).

    Example:
        >>> # Venus (mag -4.0) at 30 degrees altitude
        >>> apparent_magnitude_with_extinction(-4.0, 30.0)
        -3.4...  # Slightly dimmed by extinction

        >>> # Sirius (mag -1.46) near horizon (5 degrees)
        >>> apparent_magnitude_with_extinction(-1.46, 5.0)
        1.7...  # Significantly dimmed

        >>> # A 6th magnitude star at the zenith
        >>> apparent_magnitude_with_extinction(6.0, 90.0)
        6.28...  # Minimal extinction
    """
    extinction = calc_extinction_magnitude(
        altitude_deg=altitude_deg,
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        humidity_percent=humidity_percent,
        observer_altitude_m=observer_altitude_m,
        wavelength_nm=wavelength_nm,
    )

    return catalog_magnitude + extinction


def get_extinction_for_heliacal(
    zenith_angle_deg: float,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
) -> float:
    """
    Calculate extinction coefficient for heliacal visibility calculations.

    This is a convenience function that returns the extinction coefficient
    in the format expected by heliacal event calculations.

    Args:
        zenith_angle_deg: Zenith angle of the object in degrees (0 = zenith, 90 = horizon).
        pressure_mbar: Atmospheric pressure in millibars.
        temperature_c: Temperature in degrees Celsius.
        humidity_percent: Relative humidity (0-100).
        altitude_m: Observer altitude in meters.

    Returns:
        Extinction coefficient (k value) suitable for heliacal calculations.
        This is the total extinction coefficient in mag/airmass.

    Example:
        >>> k = get_extinction_for_heliacal(60.0)  # 30 degrees altitude
        >>> print(f"Extinction coefficient: {k:.3f}")
        Extinction coefficient: 0.28...
    """
    coeff = calc_extinction_coefficient(
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        humidity_percent=humidity_percent,
        altitude_m=altitude_m,
        wavelength_nm=WAVELENGTH_V,
    )
    return coeff.k_total


# =============================================================================
# TWILIGHT SKY BRIGHTNESS MODEL
# =============================================================================
# The twilight sky brightness model calculates the surface brightness of the sky
# during the three phases of twilight:
#   - Civil twilight: Sun altitude 0° to -6°
#   - Nautical twilight: Sun altitude -6° to -12°
#   - Astronomical twilight: Sun altitude -12° to -18°
#
# The model is based on the work of:
#   - Patat, F. (2003) "UBVRI twilight sky brightness at ESO-Paranal"
#   - Krisciunas, K. & Schaefer, B.E. (1991) "A model of the brightness of moonlight"
#   - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
#   - Rozenberg, G.V. (1966) "Twilight: A Study in Atmospheric Optics"
# =============================================================================

# Twilight phase boundaries (Sun altitude in degrees)
TWILIGHT_CIVIL_START = 0.0  # Sun at horizon
TWILIGHT_CIVIL_END = -6.0  # End of civil twilight
TWILIGHT_NAUTICAL_END = -12.0  # End of nautical twilight
TWILIGHT_ASTRONOMICAL_END = -18.0  # End of astronomical twilight (full darkness)

# Sky brightness constants in mag/arcsec^2 (V-band)
# Dark sky brightness is typically 21.5-22.0 mag/arcsec^2 at excellent sites
DARK_SKY_BRIGHTNESS_V = 21.7  # Typical dark sky, V-band
ZENITH_DARK_SKY = 21.9  # Zenith dark sky brightness

# Twilight sky brightness at key Sun altitudes (approximate, V-band)
# These are representative values for zenith in mag/arcsec^2
SKY_BRIGHTNESS_SUN_0 = 3.0  # Sun at horizon (very bright)
SKY_BRIGHTNESS_SUN_MINUS6 = 8.0  # End of civil twilight
SKY_BRIGHTNESS_SUN_MINUS12 = 17.0  # End of nautical twilight
SKY_BRIGHTNESS_SUN_MINUS18 = 21.7  # End of astronomical twilight


@dataclass
class TwilightSkyBrightness:
    """
    Result of twilight sky brightness calculation.

    Attributes:
        surface_brightness: Sky surface brightness in mag/arcsec^2 (V-band).
                           Higher values = darker sky.
        twilight_phase: Current twilight phase ("day", "civil", "nautical",
                        "astronomical", or "night").
        sun_altitude_deg: Sun altitude in degrees (negative = below horizon).
        azimuth_factor: Angular distance factor from Sun (0-1, 1 = away from Sun).
        limiting_magnitude: Approximate naked-eye limiting magnitude.
        nanolamberts: Sky brightness in nanoLamberts (alternative unit).
    """

    surface_brightness: float  # mag/arcsec^2
    twilight_phase: str
    sun_altitude_deg: float
    azimuth_factor: float
    limiting_magnitude: float
    nanolamberts: float


def get_twilight_phase(sun_altitude_deg: float) -> str:
    """
    Determine the current twilight phase based on Sun altitude.

    Args:
        sun_altitude_deg: Sun altitude in degrees (negative = below horizon).

    Returns:
        Twilight phase as a string:
            - "day": Sun above horizon (altitude > 0)
            - "civil": Civil twilight (0 >= altitude > -6)
            - "nautical": Nautical twilight (-6 >= altitude > -12)
            - "astronomical": Astronomical twilight (-12 >= altitude > -18)
            - "night": Full darkness (altitude <= -18)

    Example:
        >>> get_twilight_phase(-3.0)
        'civil'
        >>> get_twilight_phase(-9.0)
        'nautical'
        >>> get_twilight_phase(-15.0)
        'astronomical'
        >>> get_twilight_phase(-20.0)
        'night'
    """
    if sun_altitude_deg > TWILIGHT_CIVIL_START:
        return "day"
    elif sun_altitude_deg > TWILIGHT_CIVIL_END:
        return "civil"
    elif sun_altitude_deg > TWILIGHT_NAUTICAL_END:
        return "nautical"
    elif sun_altitude_deg > TWILIGHT_ASTRONOMICAL_END:
        return "astronomical"
    else:
        return "night"


def _calc_azimuth_factor(
    sun_azimuth_deg: float,
    target_azimuth_deg: float,
    sun_altitude_deg: float,
) -> float:
    """
    Calculate the azimuthal brightness factor based on angular distance from Sun.

    The sky is brighter in the direction of the Sun during twilight.
    This function returns a factor between 0 and 1, where:
    - 0 = looking toward the Sun (brightest)
    - 1 = looking away from the Sun (darkest)

    Args:
        sun_azimuth_deg: Sun azimuth in degrees (0-360).
        target_azimuth_deg: Target azimuth in degrees (0-360).
        sun_altitude_deg: Sun altitude in degrees.

    Returns:
        Azimuth factor between 0 and 1.
    """
    # Calculate azimuth difference (0-180 degrees)
    delta_az = abs(sun_azimuth_deg - target_azimuth_deg)
    if delta_az > 180.0:
        delta_az = 360.0 - delta_az

    # Normalize to 0-1 range (0 = toward Sun, 1 = opposite Sun)
    az_factor = delta_az / 180.0

    # The effect is stronger when Sun is closer to horizon
    # At lower Sun altitudes, the brightness gradient is more pronounced
    if sun_altitude_deg > -6.0:
        # Civil twilight: strong azimuthal gradient
        weight = 0.7
    elif sun_altitude_deg > -12.0:
        # Nautical twilight: moderate gradient
        weight = 0.4
    elif sun_altitude_deg > -18.0:
        # Astronomical twilight: weak gradient
        weight = 0.2
    else:
        # Night: negligible gradient
        weight = 0.0

    # Return weighted azimuth factor
    return az_factor * weight + (1.0 - weight)


def _calc_altitude_factor(target_altitude_deg: float) -> float:
    """
    Calculate brightness factor based on target altitude.

    Sky brightness varies with zenith angle. The sky is generally
    brighter near the horizon due to longer light path through atmosphere.

    Args:
        target_altitude_deg: Target altitude in degrees.

    Returns:
        Altitude factor (multiplicative adjustment to brightness).
    """
    if target_altitude_deg <= 0:
        return 0.5  # Below horizon - very bright sky glow at horizon

    # Zenith angle in degrees
    zenith_angle = 90.0 - target_altitude_deg

    # Van Rhijn function approximation for sky brightness vs zenith angle
    # B(z) = B_0 * (1 + k * sec(z))
    # At zenith (z=0): factor = 1
    # At horizon (z=90): factor is higher (brighter near horizon)
    if zenith_angle >= 85:
        return 1.5  # Near horizon - significantly brighter

    zenith_rad = math.radians(zenith_angle)
    cos_z = math.cos(zenith_rad)

    # Simple model: brightness increases as secant of zenith angle
    # but limited to prevent infinities near horizon
    sec_z = min(1.0 / cos_z, 10.0) if cos_z > 0.1 else 10.0

    # Normalize so zenith = 1.0
    factor = 0.2 * (sec_z - 1.0) + 1.0

    return min(factor, 2.0)


def calc_twilight_sky_brightness(
    sun_altitude_deg: float,
    target_altitude_deg: float = 90.0,
    sun_azimuth_deg: float = 0.0,
    target_azimuth_deg: float = 180.0,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    temperature_c: float = DEFAULT_TEMPERATURE_C,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
) -> TwilightSkyBrightness:
    """
    Calculate sky surface brightness during twilight.

    This function models the sky brightness as a function of Sun altitude,
    azimuth relative to the target, and atmospheric conditions. The model
    covers all three twilight phases:

    - Civil twilight (Sun 0° to -6°): Brightest stars visible, horizon visible
    - Nautical twilight (Sun -6° to -12°): Horizon line barely visible
    - Astronomical twilight (Sun -12° to -18°): Sky still not fully dark

    Args:
        sun_altitude_deg: Sun altitude in degrees (negative = below horizon).
                         Range: typically -18 to 0 for twilight.
        target_altitude_deg: Altitude of viewing direction in degrees.
                            Default 90.0 (zenith).
        sun_azimuth_deg: Sun azimuth in degrees (0-360, N=0, E=90).
        target_azimuth_deg: Target azimuth in degrees.
        pressure_mbar: Atmospheric pressure in millibars.
        temperature_c: Temperature in degrees Celsius.
        humidity_percent: Relative humidity (0-100).
        altitude_m: Observer altitude in meters above sea level.
        wavelength_nm: Wavelength of observation in nanometers.

    Returns:
        TwilightSkyBrightness dataclass containing:
            - surface_brightness: Sky brightness in mag/arcsec^2
            - twilight_phase: Current twilight phase
            - sun_altitude_deg: Input Sun altitude
            - azimuth_factor: Angular distance factor from Sun
            - limiting_magnitude: Approximate naked-eye limiting magnitude
            - nanolamberts: Sky brightness in nanoLamberts

    Algorithm:
        The model uses a multi-component approach:

        1. Base brightness from Sun altitude:
           The primary driver of twilight sky brightness is the Sun's
           altitude below the horizon. We use empirical relations from
           Patat (2003) and Rozenberg (1966).

        2. Azimuthal variation:
           The sky is brighter toward the Sun. This effect is most
           pronounced during civil twilight and diminishes as the
           Sun goes deeper below the horizon.

        3. Altitude/zenith angle variation:
           The sky is generally brighter near the horizon due to
           longer atmospheric path length (Van Rhijn effect).

        4. Atmospheric conditions:
           Aerosols and humidity affect scattering and can brighten
           the twilight sky, especially near the horizon.

    Example:
        >>> # Civil twilight, looking at zenith
        >>> result = calc_twilight_sky_brightness(-3.0, 90.0)
        >>> print(f"Phase: {result.twilight_phase}")
        Phase: civil
        >>> print(f"Brightness: {result.surface_brightness:.1f} mag/arcsec^2")
        Brightness: 5.5 mag/arcsec^2

        >>> # Nautical twilight, looking away from Sun
        >>> result = calc_twilight_sky_brightness(-9.0, 45.0, 270.0, 90.0)
        >>> print(f"Phase: {result.twilight_phase}")
        Phase: nautical
        >>> print(f"Brightness: {result.surface_brightness:.1f} mag/arcsec^2")
        Brightness: 14.2 mag/arcsec^2

        >>> # Astronomical twilight, looking at zenith
        >>> result = calc_twilight_sky_brightness(-15.0)
        >>> print(f"Phase: {result.twilight_phase}")
        Phase: astronomical
        >>> print(f"Limit mag: {result.limiting_magnitude:.1f}")
        Limit mag: 5.5

    References:
        - Patat, F. (2003) "UBVRI twilight sky brightness at ESO-Paranal"
          A&A 400, 1183-1198
        - Krisciunas, K. & Schaefer, B.E. (1991) PASP 103, 1033
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Rozenberg, G.V. (1966) "Twilight: A Study in Atmospheric Optics"
    """
    # Determine twilight phase
    phase = get_twilight_phase(sun_altitude_deg)

    # Calculate base sky brightness from Sun altitude
    if sun_altitude_deg >= 0:
        # Daytime - very bright sky
        base_brightness = 3.0 - sun_altitude_deg * 0.05  # Brighter as Sun goes up
        base_brightness = max(0.0, base_brightness)
    elif sun_altitude_deg >= TWILIGHT_CIVIL_END:
        # Civil twilight (0 to -6 degrees)
        # Interpolate between Sun at horizon (3.0) and end of civil (-6 -> 8.0)
        fraction = -sun_altitude_deg / 6.0
        base_brightness = SKY_BRIGHTNESS_SUN_0 + fraction * (
            SKY_BRIGHTNESS_SUN_MINUS6 - SKY_BRIGHTNESS_SUN_0
        )
    elif sun_altitude_deg >= TWILIGHT_NAUTICAL_END:
        # Nautical twilight (-6 to -12 degrees)
        fraction = (-sun_altitude_deg - 6.0) / 6.0
        base_brightness = SKY_BRIGHTNESS_SUN_MINUS6 + fraction * (
            SKY_BRIGHTNESS_SUN_MINUS12 - SKY_BRIGHTNESS_SUN_MINUS6
        )
    elif sun_altitude_deg >= TWILIGHT_ASTRONOMICAL_END:
        # Astronomical twilight (-12 to -18 degrees)
        fraction = (-sun_altitude_deg - 12.0) / 6.0
        base_brightness = SKY_BRIGHTNESS_SUN_MINUS12 + fraction * (
            SKY_BRIGHTNESS_SUN_MINUS18 - SKY_BRIGHTNESS_SUN_MINUS12
        )
    else:
        # Full night
        base_brightness = DARK_SKY_BRIGHTNESS_V

    # Calculate azimuth factor (looking toward vs away from Sun)
    az_factor = _calc_azimuth_factor(
        sun_azimuth_deg, target_azimuth_deg, sun_altitude_deg
    )

    # Calculate altitude factor (zenith vs horizon)
    alt_factor = _calc_altitude_factor(target_altitude_deg)

    # Adjust brightness for azimuth
    # Looking toward Sun = brighter (lower mag/arcsec^2)
    # Looking away = darker (higher mag/arcsec^2)
    # Effect is up to ~2 magnitudes during civil twilight
    if sun_altitude_deg >= TWILIGHT_CIVIL_END:
        az_adjustment = (1.0 - az_factor) * 2.0  # Up to 2 mag brighter toward Sun
    elif sun_altitude_deg >= TWILIGHT_NAUTICAL_END:
        az_adjustment = (1.0 - az_factor) * 1.0  # Up to 1 mag brighter
    elif sun_altitude_deg >= TWILIGHT_ASTRONOMICAL_END:
        az_adjustment = (1.0 - az_factor) * 0.3  # Small effect
    else:
        az_adjustment = 0.0  # No effect at night

    # Adjust for altitude factor
    # Looking at horizon = brighter (lower mag)
    # Looking at zenith = darker (higher mag)
    if alt_factor > 1.0:
        alt_adjustment = -(alt_factor - 1.0) * 1.5  # Brighter near horizon
    else:
        alt_adjustment = 0.0

    # Atmospheric conditions affect brightness
    # Higher aerosol content = more scattering = brighter twilight sky
    coeff = calc_extinction_coefficient(
        pressure_mbar=pressure_mbar,
        temperature_c=temperature_c,
        humidity_percent=humidity_percent,
        altitude_m=altitude_m,
        wavelength_nm=wavelength_nm,
    )

    # Aerosol effect: more aerosols = brighter twilight
    # Typical k_aerosol is 0.1-0.2; higher values brighten the sky
    aerosol_adjustment = -(coeff.k_aerosol - 0.1) * 0.5  # Subtle effect

    # Calculate final surface brightness
    surface_brightness = (
        base_brightness - az_adjustment + alt_adjustment + aerosol_adjustment
    )

    # Clamp to reasonable range
    surface_brightness = max(0.0, min(surface_brightness, 22.5))

    # Calculate limiting magnitude from sky brightness
    # Empirical relation: limiting mag depends on sky brightness and observer
    # Based on Schaefer (1990): m_lim = 7.93 - 5*log10(1 + 10^(4.316 - B/5))
    # Simplified approximation:
    if surface_brightness < 10:
        limiting_mag = -2.0 + surface_brightness * 0.5  # Very bright sky
    elif surface_brightness < 18:
        limiting_mag = 3.0 + (surface_brightness - 10) * 0.35  # Twilight
    else:
        limiting_mag = 5.8 + (surface_brightness - 18) * 0.2  # Approaching dark sky

    limiting_mag = max(-2.0, min(limiting_mag, 7.0))

    # Convert mag/arcsec^2 to nanoLamberts
    # Standard photometric relation: nL = 10^((26.33 - B)/2.5)
    # (B = 26.33 mag/arcsec^2 corresponds to 1 nL; the previous
    # 10^(35.96 - 0.4*B) form was off by ~25 orders of magnitude)
    nanolamberts = 10 ** ((26.33 - surface_brightness) / 2.5)

    return TwilightSkyBrightness(
        surface_brightness=surface_brightness,
        twilight_phase=phase,
        sun_altitude_deg=sun_altitude_deg,
        azimuth_factor=az_factor,
        limiting_magnitude=limiting_mag,
        nanolamberts=nanolamberts,
    )


def calc_twilight_brightness_simple(
    sun_altitude_deg: float,
) -> float:
    """
    Calculate approximate sky brightness during twilight (simplified model).

    This is a simplified version of calc_twilight_sky_brightness that
    returns only the zenith sky brightness based on Sun altitude.

    Args:
        sun_altitude_deg: Sun altitude in degrees (negative = below horizon).

    Returns:
        Sky surface brightness in mag/arcsec^2 (V-band).
        Higher values indicate a darker sky.

    Example:
        >>> calc_twilight_brightness_simple(0.0)   # Sun at horizon
        3.0
        >>> calc_twilight_brightness_simple(-6.0)  # End of civil twilight
        8.0
        >>> calc_twilight_brightness_simple(-12.0) # End of nautical twilight
        17.0
        >>> calc_twilight_brightness_simple(-18.0) # End of astronomical twilight
        21.7

    References:
        - Patat, F. (2003) "UBVRI twilight sky brightness at ESO-Paranal"
    """
    result = calc_twilight_sky_brightness(
        sun_altitude_deg=sun_altitude_deg,
        target_altitude_deg=90.0,  # Zenith
        sun_azimuth_deg=0.0,
        target_azimuth_deg=180.0,  # Opposite to Sun
    )
    return result.surface_brightness


# =============================================================================
# OBSERVER EXPERIENCE CONSTANTS
# =============================================================================
# Constants defining observer experience levels for visibility calculations.
# These affect the threshold contrast that an observer can detect.

OBSERVER_SKILL_INEXPERIENCED = 1
OBSERVER_SKILL_AVERAGE = 2
OBSERVER_SKILL_EXPERIENCED = 3
OBSERVER_SKILL_EXPERT = 4

# Experience factor applied to threshold - lower = can detect fainter objects
EXPERIENCE_FACTORS = {
    OBSERVER_SKILL_INEXPERIENCED: 1.3,  # 30% worse than average
    OBSERVER_SKILL_AVERAGE: 1.0,  # Baseline
    OBSERVER_SKILL_EXPERIENCED: 0.85,  # 15% better than average
    OBSERVER_SKILL_EXPERT: 0.7,  # 30% better than average
}


def calc_limiting_magnitude_twilight(
    sun_altitude_deg: float,
    target_altitude_deg: float = 45.0,
    sun_azimuth_deg: float = 0.0,
    target_azimuth_deg: float = 180.0,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
) -> float:
    """
    Calculate the naked-eye limiting magnitude during twilight.

    This function estimates the faintest star visible to the naked eye
    given the current twilight conditions. The limiting magnitude depends
    primarily on sky brightness, which varies with Sun altitude and
    viewing direction.

    Args:
        sun_altitude_deg: Sun altitude in degrees (negative = below horizon).
        target_altitude_deg: Altitude of viewing direction in degrees.
        sun_azimuth_deg: Sun azimuth in degrees.
        target_azimuth_deg: Target azimuth in degrees.
        pressure_mbar: Atmospheric pressure in millibars.
        humidity_percent: Relative humidity (0-100).
        altitude_m: Observer altitude in meters.

    Returns:
        Limiting visual magnitude. Stars fainter than this value
        will not be visible to the naked eye.

    Example:
        >>> # Civil twilight - only bright objects visible
        >>> calc_limiting_magnitude_twilight(-3.0)
        3.0...

        >>> # Nautical twilight - more stars visible
        >>> calc_limiting_magnitude_twilight(-9.0)
        4.5...

        >>> # Astronomical twilight - approaching dark sky limit
        >>> calc_limiting_magnitude_twilight(-15.0)
        5.5...

        >>> # Full night - typical naked-eye limit
        >>> calc_limiting_magnitude_twilight(-25.0)
        6.0...

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Crumey, A. (2014) "Human contrast threshold and astronomical visibility"
    """
    result = calc_twilight_sky_brightness(
        sun_altitude_deg=sun_altitude_deg,
        target_altitude_deg=target_altitude_deg,
        sun_azimuth_deg=sun_azimuth_deg,
        target_azimuth_deg=target_azimuth_deg,
        pressure_mbar=pressure_mbar,
        humidity_percent=humidity_percent,
        altitude_m=altitude_m,
    )

    # Apply atmospheric extinction to the limiting magnitude
    # Objects near horizon appear fainter due to extinction
    airmass = calc_airmass(target_altitude_deg)
    coeff = calc_extinction_coefficient(
        pressure_mbar=pressure_mbar,
        humidity_percent=humidity_percent,
        altitude_m=altitude_m,
    )

    # Reduce limiting magnitude for objects seen through more atmosphere
    # (fainter objects cannot be seen through thick atmosphere)
    extinction_penalty = coeff.k_total * (airmass - 1.0)

    adjusted_limit = result.limiting_magnitude - extinction_penalty

    return max(-2.0, min(adjusted_limit, 6.5))


# =============================================================================
# SCHAEFER VISIBILITY THRESHOLD MODEL
# =============================================================================
# The Schaefer visibility model determines whether a celestial object of a
# given apparent magnitude can be detected against a sky of given brightness.
#
# The model accounts for:
# - Object brightness (apparent magnitude)
# - Sky surface brightness (background)
# - Observer's eye adaptation state
# - Observer's experience and skill level
#
# Based on the work of:
#   - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
#   - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
#   - Blackwell, H.R. (1946) "Contrast thresholds of the human eye"
#   - Crumey, A. (2014) "Human contrast threshold and astronomical visibility"
# =============================================================================


@dataclass
class VisibilityResult:
    """
    Result of visibility threshold calculation.

    Attributes:
        is_visible: Whether the object is visible to the observer.
        object_magnitude: The apparent magnitude of the object.
        limiting_magnitude: The limiting magnitude at this sky brightness.
        sky_brightness: Sky surface brightness in mag/arcsec^2.
        contrast: The contrast ratio between object and sky.
        threshold_contrast: The minimum contrast required for detection.
        visibility_margin: Magnitude margin (positive = visible, negative = not visible).
        eye_adaptation: The eye adaptation state used ("dark", "mesopic", or "photopic").
        observer_skill: The observer skill level used (1-4).
    """

    is_visible: bool
    object_magnitude: float
    limiting_magnitude: float
    sky_brightness: float
    contrast: float
    threshold_contrast: float
    visibility_margin: float
    eye_adaptation: str
    observer_skill: int


def calc_eye_adaptation_state(sky_brightness: float) -> str:
    """
    Determine the eye's adaptation state based on sky brightness.

    The human eye adapts to ambient light levels, transitioning between
    three states: photopic (cone vision, bright conditions), mesopic
    (mixed rod/cone vision, twilight), and scotopic (rod vision, dark).

    Args:
        sky_brightness: Sky surface brightness in mag/arcsec^2.
                       Lower values = brighter sky.

    Returns:
        Adaptation state as a string:
            - "photopic": Bright conditions (sky < 16 mag/arcsec^2)
            - "mesopic": Twilight conditions (16 <= sky < 20)
            - "dark": Dark-adapted (sky >= 20)

    Algorithm:
        The transitions are based on approximate luminance levels:
        - Photopic: > 3 cd/m^2 (roughly < 16 mag/arcsec^2)
        - Mesopic: 0.001 to 3 cd/m^2 (16-20 mag/arcsec^2)
        - Scotopic: < 0.001 cd/m^2 (> 20 mag/arcsec^2)

    Example:
        >>> calc_eye_adaptation_state(15.0)  # Civil twilight
        'photopic'
        >>> calc_eye_adaptation_state(18.0)  # Late twilight
        'mesopic'
        >>> calc_eye_adaptation_state(21.5)  # Dark sky
        'dark'

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
    """
    if sky_brightness < 16.0:
        return "photopic"
    elif sky_brightness < 20.0:
        return "mesopic"
    else:
        return "dark"


def calc_contrast_threshold(
    sky_brightness: float,
    eye_adaptation: Optional[str] = None,
    observer_skill: int = OBSERVER_SKILL_AVERAGE,
) -> float:
    """
    Calculate the minimum contrast required for object detection.

    The contrast threshold represents the minimum brightness difference
    (expressed as a ratio) that an observer can detect between an object
    and the background sky. This threshold depends on the observer's
    eye adaptation state, background brightness, and experience.

    Args:
        sky_brightness: Sky surface brightness in mag/arcsec^2.
        eye_adaptation: Eye adaptation state ("photopic", "mesopic", or "dark").
                       If None, automatically determined from sky brightness.
        observer_skill: Observer skill level (1-4):
            - 1 (OBSERVER_SKILL_INEXPERIENCED): Inexperienced observer
            - 2 (OBSERVER_SKILL_AVERAGE): Average observer
            - 3 (OBSERVER_SKILL_EXPERIENCED): Experienced observer
            - 4 (OBSERVER_SKILL_EXPERT): Expert observer

    Returns:
        Contrast threshold (dimensionless ratio). Lower values mean the
        observer can detect fainter contrasts.

    Algorithm:
        The contrast threshold follows the Blackwell-Schaefer model:

        For dark-adapted vision (scotopic):
            C_threshold = 0.017 * (1 + sqrt(B_sky / 10^5))

        For mesopic vision:
            C_threshold = 0.025 * (1 + B_sky^0.4 / 100)

        For photopic vision:
            C_threshold = 0.05 * (1 + B_sky^0.3 / 10)

        where B_sky is the sky brightness in nanoLamberts.

        The threshold is then adjusted for observer experience.

    Example:
        >>> # Dark sky, average observer
        >>> calc_contrast_threshold(21.5)
        0.018...
        >>> # Twilight, experienced observer
        >>> calc_contrast_threshold(15.0, observer_skill=OBSERVER_SKILL_EXPERIENCED)
        0.04...

    References:
        - Blackwell, H.R. (1946) "Contrast thresholds of the human eye"
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Crumey, A. (2014) "Human contrast threshold and astronomical visibility"
    """
    # Determine eye adaptation state if not provided
    if eye_adaptation is None:
        eye_adaptation = calc_eye_adaptation_state(sky_brightness)

    # The Schaefer visibility model relates limiting magnitude to sky brightness.
    # We compute the contrast threshold C such that:
    #   limiting_magnitude = sky_brightness - 2.5 * log10(C)
    #
    # Rearranging from known limiting magnitudes at various sky brightnesses:
    #   C = 10^((sky_brightness - limiting_magnitude) / 2.5)
    #
    # For typical conditions:
    #   - Dark sky (B=21.5): m_lim ~6.0, so delta_m = 15.5
    #   - Mesopic (B=18): m_lim ~5.0, so delta_m = 13.0
    #   - Nautical twilight (B=15): m_lim ~4.0, so delta_m = 11.0
    #   - Civil twilight (B=8): m_lim ~0.5, so delta_m = 7.5
    #   - Day sky (B=3): m_lim ~-4.0, so delta_m = 7.0

    # Calculate limiting magnitude based on sky brightness and adaptation
    if eye_adaptation == "dark":
        # Dark-adapted vision: best sensitivity
        # At B=21.5, limit ~6.0; at B=20, limit ~5.5
        base_limit = 6.0
        # Adjust for actual sky brightness relative to reference (21.5)
        brightness_factor = (sky_brightness - 21.5) * 0.3
        limiting_mag = base_limit + brightness_factor
    elif eye_adaptation == "mesopic":
        # Mesopic (twilight) vision
        # At B=18, limit ~4.5-5.0; at B=16, limit ~3.5-4.0
        # Linear interpolation in this range
        limiting_mag = 2.0 + (sky_brightness - 16.0) * 0.5
    else:  # photopic
        # Photopic (daylight) vision: worst sensitivity for point sources
        # At B=8 (civil twilight), limit ~0-1.0
        # At B=3 (daylight), only planets like Venus visible (limit ~ -4)
        if sky_brightness >= 10:
            limiting_mag = 1.0 + (sky_brightness - 10.0) * 0.4
        elif sky_brightness >= 5:
            limiting_mag = -2.0 + (sky_brightness - 5.0) * 0.6
        else:
            limiting_mag = -4.5 + (sky_brightness) * 0.5

    # Apply observer experience adjustment
    # Better observers can see about 0.5-1.0 magnitude fainter
    experience_factor = EXPERIENCE_FACTORS.get(observer_skill, 1.0)
    # Factor 1.0 = average, <1.0 = better
    # Each 0.1 decrease in factor = ~0.3 mag improvement
    experience_adj = (1.0 - experience_factor) * 3.0  # +0.9 mag for experts
    limiting_mag = limiting_mag + experience_adj

    # Clamp limiting magnitude to reasonable range
    limiting_mag = max(-5.0, min(limiting_mag, 7.5))

    # Calculate contrast threshold from limiting magnitude
    # C = 10^((sky_brightness - limiting_mag) / 2.5)
    delta_m = sky_brightness - limiting_mag
    threshold = 10 ** (delta_m / 2.5)

    # Clamp to reasonable range
    return max(1e-8, min(threshold, 1e10))


def calc_visibility_threshold(
    object_magnitude: float,
    sky_brightness: float,
    eye_adaptation: Optional[str] = None,
    observer_skill: int = OBSERVER_SKILL_AVERAGE,
    object_altitude_deg: float = 90.0,
    apply_extinction: bool = False,
    pressure_mbar: float = DEFAULT_PRESSURE_MBAR,
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    observer_altitude_m: float = DEFAULT_ALTITUDE_M,
) -> VisibilityResult:
    """
    Determine if an object is visible using the Schaefer visibility model.

    This function implements the Schaefer visibility threshold model, which
    determines whether a celestial object of given apparent magnitude can
    be detected against a sky of given surface brightness, accounting for
    the observer's eye adaptation and experience level.

    The model calculates:
    1. The limiting magnitude for the given sky brightness
    2. The contrast between the object and sky
    3. Whether this contrast exceeds the detection threshold

    Args:
        object_magnitude: Apparent visual magnitude of the object.
                         Brighter objects have lower (more negative) values.
        sky_brightness: Sky surface brightness in mag/arcsec^2.
                       Typical values: 3 (civil twilight) to 22 (dark sky).
        eye_adaptation: Eye adaptation state ("photopic", "mesopic", or "dark").
                       If None, automatically determined from sky brightness.
        observer_skill: Observer skill level (1-4):
            - 1 (OBSERVER_SKILL_INEXPERIENCED): Inexperienced observer
            - 2 (OBSERVER_SKILL_AVERAGE): Average observer
            - 3 (OBSERVER_SKILL_EXPERIENCED): Experienced observer
            - 4 (OBSERVER_SKILL_EXPERT): Expert observer
        object_altitude_deg: Altitude of object above horizon in degrees.
                            Only used if apply_extinction is True.
        apply_extinction: If True, add atmospheric extinction to object magnitude.
        pressure_mbar: Atmospheric pressure in millibars (for extinction).
        humidity_percent: Relative humidity (for extinction).
        observer_altitude_m: Observer altitude in meters (for extinction).

    Returns:
        VisibilityResult dataclass containing:
            - is_visible: Whether the object is visible
            - object_magnitude: The (possibly extinction-corrected) magnitude
            - limiting_magnitude: The limiting magnitude at this sky brightness
            - sky_brightness: The input sky brightness
            - contrast: The contrast ratio achieved
            - threshold_contrast: The minimum contrast required
            - visibility_margin: Magnitude margin (positive = visible)
            - eye_adaptation: The adaptation state used
            - observer_skill: The skill level used

    Algorithm:
        1. Calculate the contrast threshold for the given conditions
        2. Convert sky brightness to flux (nanoLamberts)
        3. Calculate object flux from magnitude
        4. Compute contrast = (object_flux - sky_flux) / sky_flux
        5. Compare contrast to threshold

        The limiting magnitude is derived from the contrast threshold:
            m_lim = sky_brightness - 2.5 * log10(C_threshold)

        An object is visible if:
            object_magnitude < limiting_magnitude

    Example:
        >>> # Venus (mag -4) during civil twilight (bright sky)
        >>> result = calc_visibility_threshold(-4.0, 8.0)
        >>> result.is_visible
        True

        >>> # 6th magnitude star during civil twilight
        >>> result = calc_visibility_threshold(6.0, 8.0)
        >>> result.is_visible
        False

        >>> # 4th magnitude star under dark sky
        >>> result = calc_visibility_threshold(4.0, 21.5)
        >>> result.is_visible
        True

        >>> # Expert observer can see fainter objects
        >>> result = calc_visibility_threshold(
        ...     6.0, 21.5, observer_skill=OBSERVER_SKILL_EXPERT
        ... )
        >>> result.is_visible
        True

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
        - Blackwell, H.R. (1946) "Contrast thresholds of the human eye"
        - Crumey, A. (2014) "Human contrast threshold and astronomical visibility"
    """
    # Apply atmospheric extinction if requested
    effective_magnitude = object_magnitude
    if apply_extinction and object_altitude_deg > 0:
        extinction = calc_extinction_magnitude(
            altitude_deg=object_altitude_deg,
            pressure_mbar=pressure_mbar,
            humidity_percent=humidity_percent,
            observer_altitude_m=observer_altitude_m,
        )
        effective_magnitude = object_magnitude + extinction

    # Determine eye adaptation state
    if eye_adaptation is None:
        eye_adaptation = calc_eye_adaptation_state(sky_brightness)

    # Calculate contrast threshold
    threshold = calc_contrast_threshold(
        sky_brightness=sky_brightness,
        eye_adaptation=eye_adaptation,
        observer_skill=observer_skill,
    )

    # Convert magnitudes to flux for contrast calculation
    # Using: flux ~ 10^(-0.4 * magnitude)
    # The zero-point is arbitrary since we only need ratios
    sky_flux = 10 ** (-0.4 * sky_brightness)
    object_flux = 10 ** (-0.4 * effective_magnitude)

    # Calculate contrast
    # For a point source against extended background:
    # C = (object_flux) / (sky_flux_per_unit_area)
    # Since we're comparing magnitudes, this simplifies to flux ratio
    if sky_flux > 0:
        contrast = object_flux / sky_flux
    else:
        contrast = float("inf")

    # Calculate limiting magnitude from threshold
    # An object is at the limit when contrast = threshold
    # So: 10^(-0.4 * m_lim) / 10^(-0.4 * sky_b) = threshold
    # => m_lim = sky_b - 2.5 * log10(threshold)
    if threshold > 0:
        limiting_mag = sky_brightness - 2.5 * math.log10(threshold)
    else:
        limiting_mag = sky_brightness + 10.0  # Very faint limit

    # An object is visible if its magnitude is brighter than the limit
    # (lower magnitude = brighter object)
    is_visible = effective_magnitude < limiting_mag

    # Calculate visibility margin (positive means visible)
    visibility_margin = limiting_mag - effective_magnitude

    return VisibilityResult(
        is_visible=is_visible,
        object_magnitude=effective_magnitude,
        limiting_magnitude=limiting_mag,
        sky_brightness=sky_brightness,
        contrast=contrast,
        threshold_contrast=threshold,
        visibility_margin=visibility_margin,
        eye_adaptation=eye_adaptation,
        observer_skill=observer_skill,
    )


def is_object_visible(
    object_magnitude: float,
    sky_brightness: float,
    observer_skill: int = OBSERVER_SKILL_AVERAGE,
    object_altitude_deg: float = 90.0,
    apply_extinction: bool = False,
) -> bool:
    """
    Simple check if an object is visible (convenience function).

    This is a simplified interface to calc_visibility_threshold() that
    returns only the visibility boolean.

    Args:
        object_magnitude: Apparent visual magnitude of the object.
        sky_brightness: Sky surface brightness in mag/arcsec^2.
        observer_skill: Observer skill level (1-4).
        object_altitude_deg: Altitude of object above horizon in degrees.
        apply_extinction: If True, add atmospheric extinction.

    Returns:
        True if the object is visible, False otherwise.

    Example:
        >>> is_object_visible(-1.5, 21.5)  # Sirius under dark sky
        True
        >>> is_object_visible(6.0, 8.0)  # 6th mag star in twilight
        False
    """
    result = calc_visibility_threshold(
        object_magnitude=object_magnitude,
        sky_brightness=sky_brightness,
        observer_skill=observer_skill,
        object_altitude_deg=object_altitude_deg,
        apply_extinction=apply_extinction,
    )
    return result.is_visible


def calc_limiting_magnitude_for_sky(
    sky_brightness: float,
    observer_skill: int = OBSERVER_SKILL_AVERAGE,
    eye_adaptation: Optional[str] = None,
) -> float:
    """
    Calculate the limiting magnitude for a given sky brightness.

    This function returns the faintest star magnitude that can be detected
    by an observer under the given sky conditions.

    Args:
        sky_brightness: Sky surface brightness in mag/arcsec^2.
        observer_skill: Observer skill level (1-4).
        eye_adaptation: Eye adaptation state. If None, auto-determined.

    Returns:
        Limiting visual magnitude. Stars fainter than this cannot be seen.

    Example:
        >>> # Dark sky, average observer
        >>> calc_limiting_magnitude_for_sky(21.5)
        6.0...
        >>> # Expert observer sees fainter
        >>> calc_limiting_magnitude_for_sky(21.5, OBSERVER_SKILL_EXPERT)
        6.5...
        >>> # Twilight sky limits visibility
        >>> calc_limiting_magnitude_for_sky(15.0)
        3.5...

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
    """
    # Use a reference object at mag 0 to get the limiting magnitude
    result = calc_visibility_threshold(
        object_magnitude=0.0,
        sky_brightness=sky_brightness,
        observer_skill=observer_skill,
        eye_adaptation=eye_adaptation,
    )
    return result.limiting_magnitude
