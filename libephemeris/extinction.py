# SPDX-License-Identifier: AGPL-3.0-only
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
    - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes", PASP 102,
      212-229, DOI 10.1086/132629
    - Green, D.W.E. (1992) "Magnitude Corrections for Atmospheric Extinction",
      International Comet Quarterly 14, 55-59
    - Hayes, D.S. & Latham, D.W. (1975) "A rediscussion of the atmospheric
      extinction and the absolute spectral-energy distribution of Vega",
      ApJ 197, 593-601 (V-band Rayleigh coefficient)
    - Angstrom, A. (1929) "On the atmospheric transmission of sun radiation
      and on dust in the air", Geografiska Annaler 11, 156-166 (aerosol
      wavelength power law)
    - Koschmieder, H. (1924) "Theorie der horizontalen Sichtweite", Beitr.
      Phys. freien Atmos. 12, 33-53 (visibility-range relation)
    - Young, A.T. (1974) "Observational Technique and Data Reduction"
    - Kasten, F. & Young, A.T. (1989) "Revised optical air mass tables and
      approximation formula", Applied Optics 28, 4735-4738,
      DOI 10.1364/AO.28.004735

Provenance:
    Rayleigh, aerosol, ozone, water-vapour, and airmass relations are empirical
    source models and are identified beside their formulas. Standard-pressure,
    temperature, humidity, wavelength, clamping, and below-horizon behavior are
    explicitly labelled defaults or numerical guards, not source observations.
    The module performs unit conversions and composition independently and does
    not recover constants from compatibility output.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple


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


# The airmass models this module implements. A method outside this set is a
# caller mistake, not a request for the default: silently computing
# Kasten-Young for it returns a valid-looking number from a model the caller
# did not ask for, and nothing in the result says so.
AIRMASS_METHODS: Tuple[str, ...] = ("secant", "kasten_young", "rozenberg")

# The adaptation states this module returns and accepts. "dark" is the token;
# "scotopic" is the name of the same physiological regime and is accepted as
# an input alias, because the documentation used to name it as though it were
# a returned value.
EYE_ADAPTATION_STATES: Tuple[str, ...] = ("photopic", "mesopic", "dark")
_EYE_ADAPTATION_ALIASES: dict = {"scotopic": "dark"}


def _validate_airmass_method(method: str) -> str:
    """Return ``method`` if this module implements it, else raise."""
    from .exceptions import InputValidationError

    if method in AIRMASS_METHODS:
        return method
    raise InputValidationError(
        f"calc_airmass: unknown method {method!r}. "
        f"Valid methods: {', '.join(AIRMASS_METHODS)}"
    )


def normalize_eye_adaptation(state: str, function: str) -> str:
    """Return the canonical adaptation token for ``state``, else raise.

    ``"scotopic"`` maps to ``"dark"``: they name the same regime, and the
    documentation used to present the former as a returned value.
    """
    from .exceptions import InputValidationError

    if isinstance(state, str):
        canonical = _EYE_ADAPTATION_ALIASES.get(state, state)
        if canonical in EYE_ADAPTATION_STATES:
            return canonical
    raise InputValidationError(
        f"{function}: unknown eye adaptation state {state!r}. "
        f"Valid states: {', '.join(EYE_ADAPTATION_STATES)} "
        f"(or 'scotopic' for 'dark')"
    )


def calc_airmass(altitude_deg: float, method: str = "kasten_young") -> float:
    """Calculate relative optical path length through the atmosphere.

    Args:
        altitude_deg: Object altitude above the horizon in degrees.
        method: One of ``"secant"``, ``"kasten_young"``, or
            ``"rozenberg"``. The default Kasten-Young model remains useful
            down to the horizon. Any other value is rejected.

    Returns:
        Dimensionless airmass. Objects at or below the horizon use the
        practical cap of 40.0.

    Raises:
        InputValidationError: If ``method`` is not one of the implemented
            models.

    Notes:
        The plane-parallel model uses ``X = sec(z)``. The Kasten-Young model
        corrects the low-altitude behavior for atmospheric curvature.

    References:
        Kasten and Young (1989), *Applied Optics* 28, 4735-4738; Rozenberg
        (1966), *Twilight: A Study in Atmospheric Optics*.
    """
    method = _validate_airmass_method(method)

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
        # Physical airmass floor is 1.0 (at the zenith); the empirical
        # denominator slightly exceeds 1 there and would return <1 otherwise.
        return min(max(1.0 / denominator, 1.0), 40.0)

    # Rozenberg (1966) formula, good for very low altitudes. The method was
    # validated on entry, so this is the last of the three.
    cos_z = math.cos(zenith_rad)
    return max(1.0 / (cos_z + 0.025 * math.exp(-11.0 * cos_z)), 1.0)


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
        - Hayes, D.S. & Latham, D.W. (1975), ApJ 197, 593-601 (V-band
          sea-level Rayleigh coefficient 0.1451 mag/airmass)
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes", PASP 102,
          212-229
        - Young, A.T. (1974) "Observational Technique and Data Reduction"
    """
    # Reference values. The V-band sea-level Rayleigh coefficient 0.1451
    # mag/airmass comes from Hayes & Latham (1975), ApJ 197, 593.
    k_ref = 0.1451  # Rayleigh coefficient at sea level, V-band
    wavelength_ref = 550.0  # V-band reference wavelength

    # Air-column reduction: pressure and altitude both describe the *same*
    # physical quantity (the mass of air above the observer). Applying both
    # would double-count it, so we use exactly one:
    #   - If the caller supplied a site pressure, it already encodes the
    #     reduced column (P ~= 1013.25 * exp(-h / H)), so use it directly.
    #   - Otherwise (pressure left at the sea-level default) derive the
    #     reduction from the observer altitude via the barometric formula.
    if pressure_mbar != DEFAULT_PRESSURE_MBAR:
        column_factor = pressure_mbar / DEFAULT_PRESSURE_MBAR
    else:
        column_factor = math.exp(-altitude_m / SCALE_HEIGHT_M)

    # Wavelength dependence (Rayleigh scattering ~ lambda^-4)
    wavelength_factor = (wavelength_ref / wavelength_nm) ** 4

    return k_ref * column_factor * wavelength_factor


def calc_aerosol_coefficient(
    humidity_percent: float = DEFAULT_HUMIDITY_PERCENT,
    altitude_m: float = DEFAULT_ALTITUDE_M,
    wavelength_nm: float = WAVELENGTH_V,
    visibility_km: Optional[float] = None,
) -> float:
    """Calculate the aerosol (Mie) extinction coefficient.

    Args:
        humidity_percent: Relative humidity from 0 to 100 percent.
        altitude_m: Observer altitude above sea level in metres.
        wavelength_nm: Observation wavelength in nanometres.
        visibility_km: Optional meteorological visibility in kilometres.

    Returns:
        Aerosol extinction in magnitudes per airmass. Visibility, when
        supplied, drives a Koschmieder estimate; otherwise humidity, altitude,
        and wavelength provide the estimate.

    References:
        Koschmieder, H. (1924), Beitr. Phys. freien Atmos. 12, 33-53, and
        Green, D.W.E. (1992), International Comet Quarterly 14, 55-59, for the
        visibility-range aerosol relation; Angstrom, A. (1929), Geografiska
        Annaler 11, 156-166, for the aerosol wavelength power law. The
        humidity/altitude fallback below is a project-calibrated estimate.
    """
    if visibility_km is not None and visibility_km > 0:
        # Direct calculation from meteorological visibility using the
        # Koschmieder (1924) relation V = 3.912 / beta (beta = extinction per
        # unit path). The constant is -ln(0.02) for the canonical 2% human
        # contrast threshold; the standard ICQ extinction procedure of
        # Green (1992, International Comet Quarterly 14, 55) uses the same
        # 3.912 value. The molecular part (~0.1066
        # mag/airmass at V) is removed to leave the aerosol coefficient.
        k_aerosol = 3.912 / visibility_km - 0.1066
        return max(0.0, k_aerosol)

    # Aerosol scale height: aerosols are concentrated in the lower troposphere
    # (order 1-2 km; e.g. Allen, "Astrophysical Quantities").
    aerosol_scale_height = 1500.0  # meters

    # Base aerosol coefficient for a clean atmosphere at V. Aerosol extinction
    # follows the Angstrom (1929) turbidity power law; 0.08 mag/airmass is a
    # project-chosen clear-air baseline.
    k_base = 0.08

    # Humidity factor (project-calibrated hygroscopic-growth term).
    humidity_fraction = humidity_percent / 100.0
    humidity_factor = 1.0 + 2.0 * humidity_fraction**2

    # Altitude factor (aerosols decrease exponentially with altitude)
    altitude_factor = math.exp(-altitude_m / aerosol_scale_height)

    # Wavelength dependence: Angstrom (1929) power law, exponent ~1.3 (weaker
    # than Rayleigh's lambda^-4).
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
        Ozone absorption has a complex wavelength dependence (strong Hartley/
        Huggins bands in the UV, weaker Chappuis bands near 600 nm). The
        piecewise values returned here are a project-calibrated step
        approximation; the ~0.016 mag/airmass V-band figure matches the small,
        roughly constant visual-band value used by Schaefer (1990), PASP 102,
        212.

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

    Note:
        The piecewise per-band coefficients (linearly scaled by humidity) are
        a project-calibrated approximation. Water-vapour absorption is
        negligible in the blue-green and grows into the red/near-IR bands.

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
    """Calculate total atmospheric extinction and its components.

    Args:
        pressure_mbar: Atmospheric pressure in hPa.
        temperature_c: Air temperature in degrees Celsius.
        humidity_percent: Relative humidity from 0 to 100 percent.
        altitude_m: Observer altitude above sea level in metres.
        wavelength_nm: Observation wavelength in nanometres.
        visibility_km: Optional meteorological visibility in kilometres.

    Returns:
        An ``ExtinctionCoefficients`` value containing the Rayleigh, aerosol,
        ozone, water-vapour, and total coefficients in magnitudes per airmass.

    Example:
        >>> coeff = calc_extinction_coefficient()
        >>> coeff.k_total > 0
        True
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
    """Calculate atmospheric extinction for an object altitude.

    Args:
        altitude_deg: Object altitude above the horizon in degrees.
        pressure_mbar: Atmospheric pressure in hPa.
        temperature_c: Air temperature in degrees Celsius.
        humidity_percent: Relative humidity from 0 to 100 percent.
        observer_altitude_m: Observer altitude above sea level in metres.
        wavelength_nm: Observation wavelength in nanometres.
        visibility_km: Optional meteorological visibility in kilometres.

    Returns:
        Extinction in magnitudes. Objects at or below the horizon return 99.0.

    Notes:
        The result is the total extinction coefficient multiplied by the
        selected airmass model.
    """
    # Objects at or below the horizon have essentially infinite extinction.
    # Use <= 0 (not < 0) to match calc_simple_extinction and avoid an ~88-mag
    # discontinuity between altitude 0.0 (airmass capped at 40) and -0.0.
    if altitude_deg <= 0:
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
    """Calculate extinction with the simple plane-parallel model.

    Args:
        altitude_deg: Object altitude above the horizon in degrees.
        k: Extinction coefficient in magnitudes per airmass.

    Returns:
        Extinction in magnitudes, or 99.0 at and below the horizon.

    Notes:
        The formula is ``k / sin(altitude)`` and is intended for altitudes
        above roughly 15 degrees. Use ``calc_extinction_magnitude()`` near
        the horizon.
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
#   - Patat, F. (2003) "UBVRI twilight sky brightness at ESO-Paranal",
#     A&A 400, 1183-1198
#   - Krisciunas, K. & Schaefer, B.E. (1991) "A model of the brightness of
#     moonlight", PASP 103, 1033-1039, DOI 10.1086/132921
#   - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes", PASP 102, 212-229
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

# Twilight zenith sky brightness anchors at key Sun altitudes (V-band,
# mag/arcsec^2). These are project-calibrated representative values,
# consistent in trend with the measured twilight profiles of Patat, F. (2003),
# A&A 400, 1183-1198, and Rozenberg (1966); linearly interpolated between
# anchors by calc_twilight_sky_brightness. Not transcribed measured constants.
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
        # Daytime - very bright sky. Project-calibrated linear darkening of the
        # zenith surface brightness with increasing Sun altitude (not a
        # published relation).
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

    # Calculate limiting magnitude from sky brightness.
    #
    # This is a project-calibrated PIECEWISE-LINEAR fit of naked-eye limiting
    # magnitude as a function of the zenith surface brightness (mag/arcsec^2),
    # NOT an evaluation of Schaefer's closed-form limiting-magnitude relation
    # (Schaefer (1990), PASP 102, 212, whose sky-brightness argument is a
    # luminance, not a mag/arcsec^2 value). The three segments are chosen to
    # reproduce the standard limits: ~ -2 mag in a bright twilight sky rising
    # to ~6-6.5 mag under a dark sky.
    if surface_brightness < 10:
        limiting_mag = -2.0 + surface_brightness * 0.5  # Very bright sky
    elif surface_brightness < 18:
        limiting_mag = 3.0 + (surface_brightness - 10) * 0.35  # Twilight
    else:
        limiting_mag = 5.8 + (surface_brightness - 18) * 0.2  # Approaching dark sky

    limiting_mag = max(-2.0, min(limiting_mag, 7.0))

    # Convert mag/arcsec^2 to nanoLamberts: nL = 10^((26.33 - B)/2.5), the
    # inverse of the surface-brightness relation in the form used by Garstang,
    # R.H. (1989), PASP 101, 306-329. The 26.33 term is the photometric zero
    # point (B = 26.33 mag/arcsec^2 corresponds to 1 nL; the previous
    # 10^(35.96 - 0.4*B) form was off by ~25 orders of magnitude).
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
        - Patat, F. (2003) "UBVRI twilight sky brightness at ESO-Paranal",
          A&A 400, 1183-1198
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

# Experience factor applied to threshold - lower = can detect fainter objects.
# Project convention (not a published calibration).
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
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes", PASP 102,
          212-229
        - Crumey, A. (2014) "Human contrast threshold and astronomical
          visibility", MNRAS 442, 2600-2619, DOI 10.1093/mnras/stu992
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
    three regimes: photopic (cone vision, bright conditions), mesopic
    (mixed rod/cone vision, twilight), and scotopic (rod vision, dark).

    The returned token for the scotopic regime is ``"dark"``. The whole
    module speaks that vocabulary — it is what the eye_adaptation parameters
    accept and what VisibilityResult reports — while the prose here used to
    name "scotopic" as though it were a returned value, so code comparing
    against that name never matched. ``"scotopic"`` is accepted wherever a
    state is taken as input.

    Args:
        sky_brightness: Sky surface brightness in mag/arcsec^2.
                       Lower values = brighter sky.

    Returns:
        One of ``EYE_ADAPTATION_STATES``:
            - "photopic": Bright conditions (sky < 16 mag/arcsec^2)
            - "mesopic": Twilight conditions (16 <= sky < 20)
            - "dark": Dark-adapted, the scotopic regime (sky >= 20)

    Algorithm:
        The transitions are based on approximate luminance levels:
        - Photopic: > 3 cd/m^2 (roughly < 16 mag/arcsec^2) -> "photopic"
        - Mesopic: 0.001 to 3 cd/m^2 (16-20 mag/arcsec^2) -> "mesopic"
        - Scotopic: < 0.001 cd/m^2 (> 20 mag/arcsec^2) -> "dark"

    Example:
        >>> calc_eye_adaptation_state(15.0)  # Civil twilight
        'photopic'
        >>> calc_eye_adaptation_state(18.0)  # Late twilight
        'mesopic'
        >>> calc_eye_adaptation_state(21.5)  # Dark sky (scotopic regime)
        'dark'

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
    """
    if sky_brightness < 16.0:
        return "photopic"
    if sky_brightness < 20.0:
        return "mesopic"
    return "dark"


def calc_contrast_threshold(
    sky_brightness: float,
    eye_adaptation: Optional[str] = None,
    observer_skill: int = OBSERVER_SKILL_AVERAGE,
) -> float:
    """
    Calculate the point-source detection threshold as a star-to-sky flux ratio.

    Inverts the model limiting magnitude for the given eye-adaptation state
    and observer skill into a dimensionless flux-ratio threshold::

        threshold = 10 ** ((sky_brightness - limiting_mag) / 2.5)

    i.e. the ratio between the sky surface brightness (per arcsec^2, expressed
    as a magnitude) and the faintest detectable point source. Brighter skies
    and less-adapted eyes give LARGER thresholds (detection is harder). This
    is the inverse of :func:`calc_visibility_threshold`'s limiting-magnitude
    logic, not the classical Blackwell contrast *fraction* (which is of order
    0.01-0.1); typical return values here span ~10^2 (daylight) to ~10^6
    (dark sky).

    Args:
        sky_brightness: Sky surface brightness in mag/arcsec^2.
        eye_adaptation: Eye adaptation state, one of ``EYE_ADAPTATION_STATES``
                       ("photopic", "mesopic", "dark"); "scotopic" is accepted
                       as an alias for "dark". If None, automatically
                       determined from sky brightness. Any other value is
                       rejected.
        observer_skill: Observer skill level (1-4):
            - 1 (OBSERVER_SKILL_INEXPERIENCED): Inexperienced observer
            - 2 (OBSERVER_SKILL_AVERAGE): Average observer
            - 3 (OBSERVER_SKILL_EXPERIENCED): Experienced observer
            - 4 (OBSERVER_SKILL_EXPERT): Expert observer

    Returns:
        Flux-ratio detection threshold (dimensionless, clamped to
        [1e-8, 1e10]). Lower values mean fainter sources are detectable
        relative to the sky.

    Algorithm:
        1. Pick the limiting magnitude for the adaptation regime (piecewise
           linear in sky brightness: dark-adapted anchored at m_lim ~6.0 for
           B=21.5; mesopic at m_lim = 3.4 + (B-16)*0.5, continuous with the
           photopic branch at B=16; photopic descending to m_lim ~ -4 in
           daylight).
        2. Add the observer-experience adjustment (up to +0.9 mag for expert).
        3. Convert to a flux ratio: 10 ** ((B - m_lim) / 2.5).

    Example:
        >>> # Dark sky, average observer: 10**((21.5-6.0)/2.5) ~ 1.6e6
        >>> round(calc_contrast_threshold(21.5))
        1584893
        >>> # Nautical twilight, experienced observer
        >>> round(calc_contrast_threshold(15.0, observer_skill=3))
        41687

    References:
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes" (the
          limiting-magnitude framing this function inverts)
        - Blackwell, H.R. (1946); Crumey, A. (2014) for the classical
          contrast-fraction formulation (NOT what this function returns)
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

    # Calculate limiting magnitude based on sky brightness and adaptation.
    # The state is validated rather than allowed to fall through: an
    # unrecognised one used to land on the photopic branch, which is the
    # same silent substitution as an unknown airmass method.
    eye_adaptation = normalize_eye_adaptation(
        eye_adaptation, "calc_contrast_threshold"
    )
    if eye_adaptation == "dark":
        # Dark-adapted vision: best sensitivity
        # At B=21.5, limit ~6.0; at B=20, limit ~5.5
        base_limit = 6.0
        # Adjust for actual sky brightness relative to reference (21.5)
        brightness_factor = (sky_brightness - 21.5) * 0.3
        limiting_mag = base_limit + brightness_factor
    elif eye_adaptation == "mesopic":
        # Mesopic (twilight) vision
        # At B=16, limit ~3.4; at B=18, limit ~4.4; at B=20, limit ~5.4
        # The intercept (3.4) makes this branch continuous with the photopic
        # branch at the 16 mag/arcsec^2 boundary: the photopic branch reaches
        # 1.0 + (16.0 - 10.0) * 0.4 = 3.4 there. Without this offset the
        # mesopic branch restarted near 2.0, so a *darker* sky produced a
        # *smaller* limiting magnitude (physically impossible). The intercept
        # keeps limiting magnitude monotonically non-decreasing as the sky
        # darkens across the photopic/mesopic transition, and stays below the
        # dark branch's value at B=20 (5.55) so monotonicity holds there too.
        limiting_mag = 3.4 + (sky_brightness - 16.0) * 0.5
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

    # Apply observer experience adjustment. The EXPERIENCE_FACTORS and the
    # factor-to-magnitude scaling (3.0) are a project convention, not a
    # published calibration: better observers can see ~0.5-1.0 mag fainter.
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
    Determine whether an object clears the project visibility threshold.

    This is a deliberately compact, Schaefer-inspired convenience model, not
    a verbatim reproduction of Schaefer's complete retinal/telescope model.
    :func:`calc_contrast_threshold` first constructs a piecewise limiting
    magnitude from sky brightness, adaptation, and observer skill; this
    function applies exactly the inverse magnitude/flux relation. The regime
    boundaries, limiting-magnitude anchors, and skill offsets are declared
    project heuristics in that helper and are not presented as measured
    universal constants.

    The model calculates:
    1. The limiting magnitude for the given sky brightness
    2. The contrast between the object and sky
    3. Whether this contrast exceeds the detection threshold

    Args:
        object_magnitude: Apparent visual magnitude of the object.
                         Brighter objects have lower (more negative) values.
        sky_brightness: Sky surface brightness in mag/arcsec^2.
                       Typical values: 3 (civil twilight) to 22 (dark sky).
        eye_adaptation: Eye adaptation state, one of ``EYE_ADAPTATION_STATES``
                       ("photopic", "mesopic", "dark"); "scotopic" is accepted
                       as an alias for "dark". If None, automatically
                       determined from sky brightness. Any other value is
                       rejected.
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
        1. Calculate the model's star-to-sky flux-ratio threshold.
        2. Apply optional atmospheric extinction to the point-source
           magnitude.
        3. Use Pogson's relation ``F/F0 = 10**(-0.4*m)`` for both magnitude
           channels.
        4. Compare ``F_object/F_sky`` with the threshold and recover the
           algebraically equivalent limiting magnitude.

        Because ``sky_brightness`` is a surface brightness (mag/arcsec^2)
        while ``object_magnitude`` is an integrated point-source magnitude,
        the ratio assumes a one-square-arcsecond effective background patch.
        It is therefore a stable API approximation, not a dimensionally
        complete model of a telescope point-spread function or the eye's
        integration area.

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
        - Schaefer, B.E. (1990), "Telescopic Limiting Magnitudes", PASP
          102, 212-229, DOI 10.1086/132629
        - Blackwell, H.R. (1946), "Contrast Thresholds of the Human Eye",
          JOSA 36, 624-643
        - Crumey, A. (2014), "Human contrast threshold and astronomical
          visibility", MNRAS 442, 2600-2619, DOI 10.1093/mnras/stu992
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

    # Determine eye adaptation state. Normalise here as well as in the
    # threshold call, so the state reported back in VisibilityResult is the
    # canonical token rather than whatever alias the caller passed.
    if eye_adaptation is None:
        eye_adaptation = calc_eye_adaptation_state(sky_brightness)
    else:
        eye_adaptation = normalize_eye_adaptation(
            eye_adaptation, "calc_visibility_threshold"
        )

    # Calculate contrast threshold
    threshold = calc_contrast_threshold(
        sky_brightness=sky_brightness,
        eye_adaptation=eye_adaptation,
        observer_skill=observer_skill,
    )

    # Pogson magnitudes obey F = F0 * 10**(-0.4*m).  We set the shared
    # normalization F0 to one because it cancels exactly in
    # F_object/F_sky; no physical photometric zero point is being estimated.
    # Treating the sky value as the flux of a one-arcsec^2 patch is the
    # project approximation documented above.
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
