# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Schaefer (1990) Atmospheric Model for Heliacal Visibility.

This module implements the Schaefer atmospheric model for determining
heliacal visibility of celestial objects. The model is based on:

- Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
  PASP 102, 212-229
- Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
  Vistas in Astronomy 36, 311-361

The model calculates:
1. Atmospheric extinction (Rayleigh scattering, aerosol scattering, ozone absorption)
2. Sky brightness from twilight, moonlight, zodiacal light, and airglow
3. Visibility thresholds based on contrast sensitivity

This implementation follows the standard approach for
compatibility with archaeoastronomical applications.

Constants:
    Most constants are from Schaefer (1990).
"""

from __future__ import annotations

import math
from typing import Tuple

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

# Rayleigh scattering constants (at sea level, wavelength ~550nm)
# From Schaefer (1990) Table 1
RAYLEIGH_COEFF_SEA_LEVEL = 0.1451  # mag/airmass at 550nm
RAYLEIGH_SCALE_HEIGHT = 8500.0  # meters (atmospheric scale height)

# Aerosol scattering constants
# k_aerosol = 0.12 * V^(-1.3) where V is visual range in km
# For V=25km (clear conditions): k_aerosol ~ 0.05
# For V=10km (average): k_aerosol ~ 0.12
# For V=5km (poor): k_aerosol ~ 0.25
AEROSOL_COEFF_CLEAR = 0.05  # Very clear (V > 40km)
AEROSOL_COEFF_AVERAGE = 0.12  # Average (V ~ 20km)
AEROSOL_COEFF_POOR = 0.25  # Poor (V ~ 5km)

# Ozone absorption coefficient
OZONE_COEFF = 0.016  # mag/airmass (essentially constant)

# Water vapor absorption (simplified)
WATER_VAPOR_COEFF = 0.014  # mag/airmass at ~50% humidity

# =============================================================================
# SKY BRIGHTNESS CONSTANTS
# =============================================================================

# Sky brightness in nanoLamberts (nL) / mag/arcsec^2
# Reference: dark sky = 21.5 mag/arcsec^2 = ~250 nL

# Airglow brightness (minimum sky brightness, no Moon, astronomical twilight)
AIRGLOW_BRIGHTNESS_NL = 145.0  # nanoLamberts
AIRGLOW_MAG_ARCSEC2 = 22.0  # mag/arcsec^2 (approximate)

# Zodiacal light brightness (varies with ecliptic coordinates)
ZODIACAL_MIN_NL = 40.0  # At ecliptic poles
ZODIACAL_MAX_NL = 500.0  # Near ecliptic plane, elongation ~90

# Full Moon sky brightness increase
FULL_MOON_SKY_BRIGHTNESS_NL = 20000.0  # Can increase sky to ~18 mag/arcsec^2

# Twilight sky brightness parameters
# Based on Rozenberg and Schaefer models
TWILIGHT_TRANSITION_ZENITH = 16.0  # degrees below horizon for dark sky

# =============================================================================
# OBSERVER PARAMETERS
# =============================================================================

# Eye adaptation parameters
# Dark adapted eye: pupil ~7mm diameter
# Light adapted: pupil ~2mm diameter
DARK_ADAPTED_PUPIL_MM = 7.0
LIGHT_ADAPTED_PUPIL_MM = 2.5

# Naked eye limiting magnitude under ideal conditions
# Young observer, dark sky, high altitude
IDEAL_NAKED_EYE_LIMIT_MAG = 6.5

# Age degradation factor (Schaefer 1993)
# Vision degrades ~0.1 mag per decade after age 20
AGE_DEGRADATION_PER_DECADE = 0.1
AGE_BASELINE = 20.0

# Snellen ratio factor
# log10(limiting mag change) ~ -2.5 * log10(Snellen ratio)
# Perfect vision (Snellen 1.0) = no change

# =============================================================================
# PTOLEMAIC ARCUS VISIONIS
# =============================================================================

# Minimum arcus visionis values from Ptolemy's Almagest
# These are the minimum altitude of object above horizon when
# Sun is at dawn/dusk, for first/last visibility
# Values depend on object magnitude

# Ptolemaic arcus visionis table (Schaefer's interpretation)
# Format: (magnitude_min, magnitude_max, arcus_visionis_degrees)
PTOLEMY_ARCUS_VISIONIS = [
    (-5.0, -2.0, 5.0),  # Very bright (Venus at brightest)
    (-2.0, 0.0, 7.0),  # Bright (Jupiter, Sirius)
    (0.0, 1.5, 9.0),  # Moderately bright
    (1.5, 3.0, 11.0),  # Average
    (3.0, 4.5, 13.0),  # Faint
    (4.5, 6.0, 15.0),  # Very faint
]

# Minimum elongation from Sun for visibility
# Depends on object brightness
MIN_ELONGATION_BRIGHT = 7.0  # Venus at brightest
MIN_ELONGATION_NORMAL = 10.0  # Typical planets
MIN_ELONGATION_FAINT = 15.0  # Faint stars


# =============================================================================
# EXTINCTION CALCULATIONS
# =============================================================================


def calc_rayleigh_extinction(pressure_mbar: float, altitude_m: float) -> float:
    """
    Calculate Rayleigh scattering extinction coefficient.

    Rayleigh scattering is caused by air molecules and depends on
    atmospheric pressure (density).

    Args:
        pressure_mbar: Atmospheric pressure in mbar/hPa
        altitude_m: Observer altitude in meters above sea level

    Returns:
        Rayleigh extinction coefficient in mag/airmass

    Reference:
        Schaefer (1990) Eq. 1
    """
    # Pressure correction (Rayleigh scattering proportional to density)
    pressure_factor = pressure_mbar / 1013.25

    # Altitude correction using barometric formula
    altitude_factor = math.exp(-altitude_m / RAYLEIGH_SCALE_HEIGHT)

    return RAYLEIGH_COEFF_SEA_LEVEL * pressure_factor * altitude_factor


def calc_aerosol_extinction(
    humidity_fraction: float,
    met_range_km: float = 0.0,
    altitude_m: float = 0.0,
) -> float:
    """
    Calculate aerosol scattering extinction coefficient.

    Aerosol scattering depends on atmospheric particles (dust, water droplets).
    Higher humidity and lower visibility increase aerosol extinction.

    Args:
        humidity_fraction: Relative humidity (0.0-1.0)
        met_range_km: Dual-use parameter:
            - ``0`` → compute the coefficient from ``humidity_fraction``;
            - ``>= 1.0`` → interpreted as the meteorological visibility range in
              km (Koschmieder: ``k = 3.912 / V - 0.106``);
            - ``0 < value < 1.0`` → interpreted DIRECTLY as the aerosol
              coefficient in mag/airmass, NOT as a sub-km visibility range. To
              model dense fog (sub-km visibility) pass the resulting coefficient
              explicitly rather than the range.
        altitude_m: Observer altitude in meters

    Returns:
        Aerosol extinction coefficient in mag/airmass

    Reference:
        Schaefer (1990) Section 3.2
    """
    # If meteorological range is given directly
    if met_range_km >= 1.0:
        # k_aerosol = 3.912 / V - k_rayleigh_approx
        # Simplified: k_aerosol = 3.912 / V - 0.106 (for visual wavelengths)
        k_aerosol = max(0.0, 3.912 / met_range_km - 0.106)
    elif 0 < met_range_km < 1.0:
        # Direct aerosol coefficient given at the observer site: it is already
        # the final mag/airmass value, so return it without altitude scaling
        # (which applies only to coefficients derived from visibility/humidity).
        return met_range_km
    else:
        # Estimate from humidity
        # Higher humidity leads to more aerosol scattering
        # Base aerosol ~0.05 at 0% humidity, increases with humidity
        k_aerosol = 0.04 + 0.18 * humidity_fraction + 0.12 * humidity_fraction**2

    # Altitude correction - aerosols decrease faster than air
    # Aerosol scale height ~1.5-2 km
    aerosol_scale_height = 1500.0  # meters
    altitude_factor = math.exp(-altitude_m / aerosol_scale_height)

    return k_aerosol * altitude_factor


def calc_ozone_extinction() -> float:
    """
    Calculate ozone absorption coefficient.

    Ozone absorption is relatively constant and small in visual wavelengths.

    Returns:
        Ozone extinction coefficient in mag/airmass
    """
    return OZONE_COEFF


def calc_water_vapor_extinction(humidity_fraction: float) -> float:
    """
    Calculate water vapor absorption.

    Water vapor has minor absorption in visual wavelengths.

    Args:
        humidity_fraction: Relative humidity (0.0-1.0)

    Returns:
        Water vapor extinction coefficient in mag/airmass
    """
    return WATER_VAPOR_COEFF * humidity_fraction


def calc_total_extinction(
    pressure_mbar: float = 1013.25,
    temperature_c: float = 15.0,
    humidity_fraction: float = 0.5,
    met_range_km: float = 0.0,
    altitude_m: float = 0.0,
) -> float:
    """
    Calculate total atmospheric extinction coefficient.

    Combines all extinction sources: Rayleigh, aerosol, ozone, water vapor.

    Args:
        pressure_mbar: Atmospheric pressure in mbar/hPa
        temperature_c: Temperature in Celsius (not currently used, for future)
        humidity_fraction: Relative humidity (0.0-1.0)
        met_range_km: Meteorological visibility range in km
        altitude_m: Observer altitude in meters

    Returns:
        Total extinction coefficient k in mag/airmass

    Example:
        >>> k = calc_total_extinction(1013.25, 15.0, 0.5, 0.0, 0.0)
        >>> print(f"Extinction: {k:.3f} mag/airmass")
    """
    k_rayleigh = calc_rayleigh_extinction(pressure_mbar, altitude_m)
    k_aerosol = calc_aerosol_extinction(humidity_fraction, met_range_km, altitude_m)
    k_ozone = calc_ozone_extinction()
    k_water = calc_water_vapor_extinction(humidity_fraction)

    return k_rayleigh + k_aerosol + k_ozone + k_water


def calc_airmass(altitude_deg: float) -> float:
    """
    Calculate airmass using Kasten & Young (1989) formula.

    Airmass is the path length through the atmosphere relative
    to the zenith path. At zenith, airmass = 1. At horizon, airmass ~ 38.

    Args:
        altitude_deg: Object altitude above horizon in degrees

    Returns:
        Airmass value (minimum 1.0, maximum ~40)

    Reference:
        Kasten, F. & Young, A.T. (1989) Applied Optics 28, 4735
    """
    if altitude_deg <= 0:
        return 40.0  # Maximum airmass at/below horizon

    if altitude_deg >= 90:
        return 1.0  # Zenith

    z = 90.0 - altitude_deg  # Zenith angle in degrees
    z_rad = math.radians(z)

    # Kasten & Young (1989) formula - good to horizon
    try:
        cos_z = math.cos(z_rad)
        denominator = cos_z + 0.50572 * (96.07995 - z) ** (-1.6364)
        airmass = 1.0 / denominator
    except (ValueError, ZeroDivisionError):  # pragma: no cover - altitude in (0,90) keeps the denominator strictly positive
        airmass = 40.0

    return min(max(airmass, 1.0), 40.0)


def calc_extinction_magnitude(
    altitude_deg: float,
    k_total: float,
) -> float:
    """
    Calculate magnitude loss due to atmospheric extinction.

    Args:
        altitude_deg: Object altitude in degrees
        k_total: Total extinction coefficient in mag/airmass

    Returns:
        Magnitude increase (dimming) due to extinction
    """
    airmass = calc_airmass(altitude_deg)
    return k_total * airmass


# =============================================================================
# SKY BRIGHTNESS CALCULATIONS
# =============================================================================


def calc_twilight_brightness(sun_altitude_deg: float) -> float:
    """
    Calculate sky brightness contribution from twilight.

    Based on Schaefer's model for twilight sky brightness.

    Args:
        sun_altitude_deg: Sun altitude in degrees (negative when below horizon)

    Returns:
        Sky brightness in nanoLamberts

    Reference:
        Schaefer (1993) twilight model; Rozenberg (1966)
    """
    if sun_altitude_deg >= 0:
        # Daytime - very bright sky
        return 1e10  # Essentially infinite

    if sun_altitude_deg <= -18:
        # Astronomical night - no twilight contribution
        return 0.0

    # Twilight zone (-18 to 0)
    # Sky brightness decreases exponentially with Sun depression
    # Using Schaefer's empirical formula

    h = -sun_altitude_deg  # Sun depression below horizon (positive)

    # Three-phase twilight model:
    # Civil twilight (0 to -6): rapid brightness change
    # Nautical twilight (-6 to -12): moderate change
    # Astronomical twilight (-12 to -18): slow change

    if h <= 6:
        # Civil twilight - exponential decrease
        # At Sun = 0: very bright (~1e8 nL)
        # At Sun = -6: bright twilight (~5e4 nL)
        b = 1e8 * math.exp(-h * 1.5)
    elif h <= 12:
        # Nautical twilight
        # At Sun = -6: ~5e4 nL
        # At Sun = -12: ~1e3 nL
        h_eff = h - 6
        b = 5e4 * math.exp(-h_eff * 0.6)
    else:
        # Astronomical twilight
        # At Sun = -12: ~1e3 nL
        # At Sun = -18: ~100 nL (approaching airglow)
        h_eff = h - 12
        b = 1e3 * math.exp(-h_eff * 0.4)

    return b


def calc_moon_brightness(
    moon_altitude_deg: float,
    moon_phase: float,
    angular_distance_deg: float,
) -> float:
    """
    Calculate sky brightness contribution from moonlight.

    Based on Krisciunas & Schaefer (1991) model.

    Args:
        moon_altitude_deg: Moon altitude in degrees
        moon_phase: Moon phase angle in degrees (0 = full, 180 = new)
        angular_distance_deg: Angular distance from Moon in degrees

    Returns:
        Sky brightness contribution in nanoLamberts

    Reference:
        Krisciunas, K. & Schaefer, B.E. (1991) PASP 103, 1033
    """
    if moon_altitude_deg <= 0:
        return 0.0  # Moon below horizon

    if angular_distance_deg < 1:
        angular_distance_deg = 1.0  # Avoid division by zero

    # Moon illuminance at full moon: ~0.25 lux
    # Moon brightness decreases as phase increases (towards new moon)

    # Illuminated fraction approximation
    phase_rad = math.radians(moon_phase)
    illuminated_fraction = (1 + math.cos(phase_rad)) / 2.0

    # Very dim at new moon
    if illuminated_fraction < 0.01:
        return 0.0

    # Moon's contribution to sky brightness
    # Depends on:
    # 1. Moon altitude (airmass of moonlight)
    # 2. Moon phase (total illumination)
    # 3. Angular distance from Moon (scattered light decreases with distance)

    # Altitude factor
    moon_airmass = calc_airmass(moon_altitude_deg)
    alt_factor = math.exp(-0.4 * (moon_airmass - 1.0))  # Extinction

    # Angular distance factor - scattered light decreases with distance
    # At 10 from Moon: high scattered light
    # At 90 from Moon: much less scattered light
    if angular_distance_deg < 10:
        dist_factor = 10.0 / (angular_distance_deg + 0.1)
    elif angular_distance_deg < 45:
        dist_factor = 1.0 + (45 - angular_distance_deg) / 35.0
    else:
        dist_factor = 1.0

    # Base brightness at full moon, Moon at zenith, 45 from Moon
    base_brightness = 3000.0  # nanoLamberts

    brightness = base_brightness * illuminated_fraction * alt_factor * dist_factor

    return brightness


def calc_airglow_brightness() -> float:
    """
    Calculate airglow sky brightness.

    Airglow is the minimum natural sky brightness from
    atmospheric chemiluminescence.

    Returns:
        Airglow brightness in nanoLamberts
    """
    return AIRGLOW_BRIGHTNESS_NL


def calc_zodiacal_brightness(
    ecliptic_latitude_deg: float,
    sun_elongation_deg: float,
) -> float:
    """
    Calculate zodiacal light sky brightness.

    Zodiacal light is sunlight scattered by interplanetary dust.
    Brightest near the ecliptic, at elongations 30-90 from Sun.

    Args:
        ecliptic_latitude_deg: Ecliptic latitude of viewing direction
        sun_elongation_deg: Angular distance from Sun

    Returns:
        Zodiacal light brightness in nanoLamberts
    """
    # Zodiacal light is brightest:
    # - Near the ecliptic (low ecliptic latitude)
    # - At moderate elongations from Sun (30-90)

    # Ecliptic latitude factor
    lat_rad = math.radians(abs(ecliptic_latitude_deg))
    lat_factor = math.exp(-lat_rad)  # Decreases away from ecliptic

    # Elongation factor - peaks around 45
    if sun_elongation_deg < 20:
        # Too close to Sun - overwhelmed by twilight
        elong_factor = 0.0
    elif sun_elongation_deg < 90:
        # Peak brightness region
        elong_factor = math.sin(math.radians(sun_elongation_deg))
    else:
        # Decreases toward anti-Sun point
        elong_factor = math.sin(math.radians(180 - sun_elongation_deg))

    # Base zodiacal brightness near ecliptic at 45 elongation
    base_brightness = 150.0  # nanoLamberts

    return base_brightness * lat_factor * elong_factor


def calc_total_sky_brightness(
    sun_altitude_deg: float,
    moon_altitude_deg: float = -90.0,
    moon_phase_deg: float = 180.0,
    moon_distance_deg: float = 90.0,
    ecliptic_latitude_deg: float = 0.0,
    sun_elongation_deg: float = 90.0,
) -> float:
    """
    Calculate total sky brightness from all sources.

    Args:
        sun_altitude_deg: Sun altitude in degrees
        moon_altitude_deg: Moon altitude in degrees
        moon_phase_deg: Moon phase angle (0=full, 180=new)
        moon_distance_deg: Angular distance from Moon
        ecliptic_latitude_deg: Ecliptic latitude of viewing direction
        sun_elongation_deg: Elongation from Sun

    Returns:
        Total sky brightness in nanoLamberts
    """
    # Twilight contribution (dominates during twilight)
    b_twilight = calc_twilight_brightness(sun_altitude_deg)

    # Moon contribution (if above horizon and not new)
    b_moon = calc_moon_brightness(moon_altitude_deg, moon_phase_deg, moon_distance_deg)

    # Airglow (always present, minimum brightness)
    b_airglow = calc_airglow_brightness()

    # Zodiacal light (only visible during astronomical twilight/night)
    if sun_altitude_deg < -12:
        b_zodiacal = calc_zodiacal_brightness(ecliptic_latitude_deg, sun_elongation_deg)
    else:
        b_zodiacal = 0.0

    # Total brightness (sum of all sources)
    return b_twilight + b_moon + b_airglow + b_zodiacal


def brightness_to_mag_arcsec2(brightness_nl: float) -> float:
    """
    Convert sky brightness from nanoLamberts to mag/arcsec^2.

    Args:
        brightness_nl: Brightness in nanoLamberts

    Returns:
        Surface brightness in mag/arcsec^2 (higher = darker)

    Reference:
        1 nL = 27.78 mag/arcsec^2 at 0 nL limit
        Dark sky (~22 mag/arcsec^2) = ~250 nL
    """
    if brightness_nl <= 0:
        return 22.0  # Dark sky limit

    # Garstang (1989) formula
    # m = 26.33 - 2.5 * log10(B) where B is in nanoLamberts
    mag_arcsec2 = 26.33 - 2.5 * math.log10(brightness_nl)

    return mag_arcsec2


def mag_arcsec2_to_brightness(mag_arcsec2: float) -> float:
    """
    Convert sky brightness from mag/arcsec^2 to nanoLamberts.

    Args:
        mag_arcsec2: Surface brightness in mag/arcsec^2

    Returns:
        Brightness in nanoLamberts
    """
    return 10 ** ((26.33 - mag_arcsec2) / 2.5)


# =============================================================================
# LIMITING MAGNITUDE CALCULATIONS
# =============================================================================


def calc_limiting_magnitude(
    sky_brightness_nl: float,
    extinction_coeff: float,
    object_altitude_deg: float,
    observer_age: float = 25.0,
    snellen_ratio: float = 1.0,
) -> float:
    """
    Calculate limiting visual magnitude for naked eye observation.

    Based on Schaefer (1990) visibility model.

    Args:
        sky_brightness_nl: Sky brightness in nanoLamberts
        extinction_coeff: Atmospheric extinction coefficient (mag/airmass)
        object_altitude_deg: Object altitude in degrees
        observer_age: Observer age in years
        snellen_ratio: Observer's Snellen ratio (1.0 = normal vision)

    Returns:
        Limiting magnitude (objects brighter than this are visible)

    Reference:
        Schaefer (1990) Eq. 18-22
    """
    # Base limiting magnitude for dark, clear sky
    # Young observer with good vision
    base_limit = IDEAL_NAKED_EYE_LIMIT_MAG

    # Sky brightness reduction
    # Brighter sky = lower limiting magnitude
    # Dark sky (~250 nL, 22 mag/arcsec^2) = base limit
    # Each factor of 10 in brightness reduces limit by ~2.5 mag
    if sky_brightness_nl > 250:
        sky_factor = 2.5 * math.log10(sky_brightness_nl / 250.0)
    else:
        sky_factor = 0.0

    # Age factor - vision degrades with age
    if observer_age > AGE_BASELINE:
        decades_past_baseline = (observer_age - AGE_BASELINE) / 10.0
        age_factor = AGE_DEGRADATION_PER_DECADE * decades_past_baseline
    else:
        age_factor = 0.0

    # Snellen ratio factor
    if snellen_ratio > 0 and snellen_ratio != 1.0:
        snellen_factor = 2.5 * math.log10(1.0 / snellen_ratio)
    else:
        snellen_factor = 0.0

    # Airmass extinction (objects appear fainter near horizon)
    # This is applied to object magnitude, but effectively reduces limiting mag
    airmass = calc_airmass(object_altitude_deg)
    extinction_loss = extinction_coeff * (
        airmass - 1.0
    )  # Extra extinction above zenith

    # Total limiting magnitude
    limiting_mag = (
        base_limit - sky_factor - age_factor - snellen_factor - extinction_loss
    )

    return limiting_mag


def is_object_visible(
    object_magnitude: float,
    object_altitude_deg: float,
    sun_altitude_deg: float,
    pressure_mbar: float = 1013.25,
    temperature_c: float = 15.0,
    humidity_fraction: float = 0.5,
    met_range_km: float = 0.0,
    altitude_m: float = 0.0,
    observer_age: float = 25.0,
    snellen_ratio: float = 1.0,
    moon_altitude_deg: float = -90.0,
    moon_phase_deg: float = 180.0,
    sun_elongation_deg: float = 90.0,
    include_moon: bool = True,
) -> Tuple[bool, float, float]:
    """
    Determine if a celestial object is visible.

    Full Schaefer model implementation.

    Args:
        object_magnitude: Visual magnitude of the object
        object_altitude_deg: Object altitude in degrees
        sun_altitude_deg: Sun altitude in degrees
        pressure_mbar: Atmospheric pressure in mbar
        temperature_c: Temperature in Celsius
        humidity_fraction: Relative humidity (0.0-1.0)
        met_range_km: Meteorological visibility range
        altitude_m: Observer altitude in meters
        observer_age: Observer age in years
        snellen_ratio: Observer's Snellen ratio
        moon_altitude_deg: Moon altitude in degrees
        moon_phase_deg: Moon phase angle
        sun_elongation_deg: Object's elongation from Sun
        include_moon: Whether to include Moon brightness

    Returns:
        Tuple of (is_visible, limiting_magnitude, apparent_magnitude)
    """
    # Object must be above horizon
    if object_altitude_deg <= 0:
        return False, 0.0, object_magnitude

    # Sun must be below horizon for heliacal visibility
    if sun_altitude_deg >= 0:
        return False, -5.0, object_magnitude

    # Calculate atmospheric extinction
    k_total = calc_total_extinction(
        pressure_mbar, temperature_c, humidity_fraction, met_range_km, altitude_m
    )

    # Calculate sky brightness
    if include_moon:
        sky_brightness = calc_total_sky_brightness(
            sun_altitude_deg,
            moon_altitude_deg,
            moon_phase_deg,
            90.0,  # Assumed distance from Moon
            0.0,  # Assumed ecliptic latitude
            sun_elongation_deg,
        )
    else:
        sky_brightness = calc_total_sky_brightness(
            sun_altitude_deg,
            -90.0,  # Moon below horizon
            180.0,
            90.0,
            0.0,
            sun_elongation_deg,
        )

    # Calculate limiting magnitude
    limiting_mag = calc_limiting_magnitude(
        sky_brightness,
        k_total,
        object_altitude_deg,
        observer_age,
        snellen_ratio,
    )

    # Apparent (observed) magnitude after full-path extinction — informational.
    airmass = calc_airmass(object_altitude_deg)
    apparent_mag = object_magnitude + k_total * airmass

    # Visibility: calc_limiting_magnitude already folds the altitude-dependent
    # extinction into the limit (it returns the faintest INTRINSIC magnitude
    # visible at this altitude), so compare the object's intrinsic magnitude
    # directly. Comparing apparent_mag here too would double-count extinction —
    # penalising the object by k*(2*airmass-1) instead of the correct k*airmass.
    is_visible = object_magnitude <= limiting_mag

    # Additional check: minimum elongation from Sun
    # Bright objects can be seen closer to Sun than faint ones
    if object_magnitude < -2:
        min_elong = MIN_ELONGATION_BRIGHT
    elif object_magnitude < 2:
        min_elong = MIN_ELONGATION_NORMAL
    else:
        min_elong = MIN_ELONGATION_FAINT

    if sun_elongation_deg < min_elong:
        is_visible = False

    return is_visible, limiting_mag, apparent_mag


# =============================================================================
# ARCUS VISIONIS (VISIBILITY ARC)
# =============================================================================


def get_arcus_visionis(object_magnitude: float) -> float:
    """
    Get the minimum arcus visionis for an object of given magnitude.

    Arcus visionis is the minimum altitude of the object above horizon
    when the Sun is at the horizon, for the object to be visible.

    Args:
        object_magnitude: Visual magnitude of the object

    Returns:
        Minimum arcus visionis in degrees

    Reference:
        Based on Ptolemy's Almagest with Schaefer's modern interpretation
    """
    for mag_min, mag_max, arcus in PTOLEMY_ARCUS_VISIONIS:
        if mag_min <= object_magnitude < mag_max:
            return arcus

    # For very faint objects
    if object_magnitude >= 6.0:
        return 17.0

    # For extremely bright objects
    return 4.0


def get_optimal_sun_altitude(object_magnitude: float, is_morning: bool) -> float:
    """
    Get the optimal Sun altitude for heliacal observation.

    The optimal time for heliacal visibility is when sky is dark enough
    but object is still visible before/after sunrise/sunset.

    Args:
        object_magnitude: Visual magnitude of the object
        is_morning: True for heliacal rising, False for setting

    Returns:
        Optimal Sun altitude in degrees (negative = below horizon)
    """
    # Brighter objects can be seen with Sun closer to horizon
    if object_magnitude < -2:
        optimal_depression = -8.0  # Venus at brightest
    elif object_magnitude < 0:
        optimal_depression = -9.0  # Bright planets/stars
    elif object_magnitude < 2:
        optimal_depression = -10.0  # Moderate stars
    elif object_magnitude < 4:
        optimal_depression = -11.0  # Faint stars
    else:
        optimal_depression = -12.0  # Very faint

    return optimal_depression


def calc_heliacal_altitude_threshold(
    object_magnitude: float,
    sun_altitude_deg: float,
    extinction_coeff: float,
) -> float:
    """
    Calculate the minimum object altitude for heliacal visibility.

    This combines arcus visionis with extinction effects.

    Args:
        object_magnitude: Visual magnitude of the object
        sun_altitude_deg: Sun altitude in degrees
        extinction_coeff: Total extinction coefficient

    Returns:
        Minimum object altitude in degrees
    """
    # Base arcus visionis
    base_av = get_arcus_visionis(object_magnitude)

    # Extinction increases the required altitude
    # Higher extinction = need object higher to compensate
    extinction_factor = 1.0 + 0.5 * (extinction_coeff - 0.25)

    # Sun altitude affects threshold
    # When Sun is deeper below horizon, lower altitude is acceptable
    if sun_altitude_deg < -12:
        sun_factor = 0.8
    elif sun_altitude_deg < -9:
        sun_factor = 1.0
    else:
        sun_factor = 1.2 + (sun_altitude_deg + 9) * 0.1

    min_altitude = base_av * extinction_factor * sun_factor

    # Minimum is always positive (above horizon)
    return max(1.0, min_altitude)


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================


def get_visibility_conditions(
    sun_altitude_deg: float,
    object_altitude_deg: float,
    object_magnitude: float,
    pressure_mbar: float = 1013.25,
    humidity_fraction: float = 0.5,
    altitude_m: float = 0.0,
) -> dict:
    """
    Get a dictionary of visibility-related conditions.

    Useful for debugging and detailed analysis.

    Args:
        sun_altitude_deg: Sun altitude in degrees
        object_altitude_deg: Object altitude in degrees
        object_magnitude: Object's visual magnitude
        pressure_mbar: Atmospheric pressure
        humidity_fraction: Relative humidity
        altitude_m: Observer altitude

    Returns:
        Dictionary with visibility parameters
    """
    k_total = calc_total_extinction(
        pressure_mbar, 15.0, humidity_fraction, 0.0, altitude_m
    )
    k_rayleigh = calc_rayleigh_extinction(pressure_mbar, altitude_m)
    k_aerosol = calc_aerosol_extinction(humidity_fraction, 0.0, altitude_m)

    sky_brightness = calc_total_sky_brightness(sun_altitude_deg)
    sky_mag = brightness_to_mag_arcsec2(sky_brightness)

    airmass = calc_airmass(object_altitude_deg)
    extinction_mag = k_total * airmass

    limiting_mag = calc_limiting_magnitude(sky_brightness, k_total, object_altitude_deg)

    apparent_mag = object_magnitude + extinction_mag

    return {
        "extinction_total": k_total,
        "extinction_rayleigh": k_rayleigh,
        "extinction_aerosol": k_aerosol,
        "airmass": airmass,
        "extinction_magnitude": extinction_mag,
        "sky_brightness_nl": sky_brightness,
        "sky_brightness_mag_arcsec2": sky_mag,
        "limiting_magnitude": limiting_mag,
        "apparent_magnitude": apparent_mag,
        # limiting_mag is an INTRINSIC-magnitude threshold (extinction already
        # folded in by calc_limiting_magnitude), so compare the object's
        # intrinsic magnitude — comparing apparent_mag would double-count it.
        "magnitude_margin": limiting_mag - object_magnitude,
        "is_visible": object_magnitude <= limiting_mag,
        "arcus_visionis": get_arcus_visionis(object_magnitude),
        "optimal_sun_altitude": get_optimal_sun_altitude(object_magnitude, True),
    }
