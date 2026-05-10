"""
Utility functions for libephemeris.

Provides helper functions compatible with the reference API including
angular calculations and other mathematical utilities.
"""

from __future__ import annotations

import math
import erfa
from typing import Optional, Sequence, Tuple, Union

# Azalt calculation method flags (compatible with the reference API)
SE_ECL2HOR: int = 0  # Ecliptic coordinates to horizontal
SE_EQU2HOR: int = 1  # Equatorial coordinates to horizontal

# Azalt_rev calculation method flags (compatible with the reference API)
SE_HOR2ECL: int = 0  # Horizontal to ecliptic coordinates
SE_HOR2EQU: int = 1  # Horizontal to equatorial coordinates

# Refraction calculation flags (compatible with the reference API)
SE_TRUE_TO_APP: int = 0  # True altitude to apparent altitude
SE_APP_TO_TRUE: int = 1  # Apparent altitude to true altitude

# reference API-compatible aliases (without SE_ prefix)
ECL2HOR: int = SE_ECL2HOR
EQU2HOR: int = SE_EQU2HOR
HOR2ECL: int = SE_HOR2ECL
HOR2EQU: int = SE_HOR2EQU
TRUE_TO_APP: int = SE_TRUE_TO_APP
APP_TO_TRUE: int = SE_APP_TO_TRUE


def cotrans_sp(
    coord: "Sequence[float]",
    eps_or_speed: "Union[float, Sequence[float]]",
    eps: "Optional[float]" = None,
) -> "Union[Tuple[float, float, float, float, float, float], Tuple[Tuple[float, float, float], Tuple[float, float, float]]]":
    """
    Transform coordinates and velocities between ecliptic and equatorial systems.

    Supports two calling conventions:
      1. PySwissEph-compatible: cotrans_sp(coord_6, eps) -> flat 6-tuple
      2. Split form: cotrans_sp(coord_3, speed_3, eps) -> (coord_3, speed_3)

    The direction of transformation depends on the sign of obliquity:
    - Negative obliquity: ecliptic (lon, lat) -> equatorial (RA, Dec)
    - Positive obliquity: equatorial (RA, Dec) -> ecliptic (lon, lat)

    Args:
        coord: Sequence of 6 floats (lon,lat,dist,lon_sp,lat_sp,dist_sp),
               or 3 floats (lon,lat,dist) when speed is given separately.
        eps_or_speed: Obliquity (float) for 2-arg form, or speed tuple for 3-arg form.
        eps: Obliquity (float) when using the 3-arg form.

    Returns:
        2-arg form: flat tuple of 6 floats
        3-arg form: (coord_3, speed_3) tuple pair
    """
    if eps is not None:
        # 3-arg form: cotrans_sp(coord_3, speed_3, eps)
        speed_seq = eps_or_speed
        obliquity = float(eps)
        lon = float(coord[0])
        lat = float(coord[1])
        dist = float(coord[2])
        lon_speed = float(speed_seq[0])
        lat_speed = float(speed_seq[1])
        dist_speed = float(speed_seq[2])
        _split_return = True
    else:
        # 2-arg form: cotrans_sp(coord_6, eps)
        obliquity = float(eps_or_speed)
        lon = float(coord[0])
        lat = float(coord[1])
        dist = float(coord[2])
        lon_speed = float(coord[3])
        lat_speed = float(coord[4])
        dist_speed = float(coord[5])
        _split_return = False
    # Convert to radians
    # Negate obliquity to match the pyswisseph API convention
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)
    eps_rad = math.radians(-obliquity)

    # Precompute trig values
    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Calculate the new latitude (Dec for ecl->eq, beta for eq->ecl)
    # sin(new_lat) = sin(lat) * cos(eps) + cos(lat) * sin(eps) * sin(lon)
    sin_new_lat = sin_lat * cos_eps + cos_lat * sin_eps * sin_lon
    # Clamp to avoid numerical issues with asin
    sin_new_lat = max(-1.0, min(1.0, sin_new_lat))
    new_lat_rad = math.asin(sin_new_lat)
    cos_new_lat = math.cos(new_lat_rad)

    # Calculate the new longitude (RA for ecl->eq, lambda for eq->ecl)
    # tan(new_lon) = (sin(lon) * cos(eps) - tan(lat) * sin(eps)) / cos(lon)
    tan_lat = math.tan(lat_rad)
    y = sin_lon * cos_eps - tan_lat * sin_eps
    x = cos_lon

    new_lon_rad = math.atan2(y, x)

    # Convert back to degrees
    new_lon = math.degrees(new_lon_rad)
    new_lat = math.degrees(new_lat_rad)

    # Normalize longitude to [0, 360)
    new_lon = new_lon % 360.0

    # --- Velocity transformation ---
    # Convert speed to radians/day for the calculation
    lon_speed_rad = math.radians(lon_speed)
    lat_speed_rad = math.radians(lat_speed)

    # Derivative of new latitude:
    # d/dt[sin(new_lat)] = cos(new_lat) * d(new_lat)/dt
    # d/dt[sin(lat)*cos(eps) + cos(lat)*sin(eps)*sin(lon)]
    #   = cos(lat)*cos(eps)*d(lat)/dt - sin(lat)*sin(eps)*sin(lon)*d(lat)/dt
    #     + cos(lat)*sin(eps)*cos(lon)*d(lon)/dt
    #
    # new_lat_speed = (1/cos(new_lat)) * [
    #   (cos(lat)*cos(eps) - sin(lat)*sin(eps)*sin(lon)) * lat_speed
    #   + cos(lat)*sin(eps)*cos(lon) * lon_speed
    # ]
    if abs(cos_new_lat) > 1e-10:
        new_lat_speed_rad = (
            (cos_lat * cos_eps - sin_lat * sin_eps * sin_lon) * lat_speed_rad
            + cos_lat * sin_eps * cos_lon * lon_speed_rad
        ) / cos_new_lat
    else:
        # At poles, latitude speed is undefined; use 0
        new_lat_speed_rad = 0.0

    # Derivative of new longitude:
    # new_lon = atan2(y, x) where y = sin(lon)*cos(eps) - tan(lat)*sin(eps), x = cos(lon)
    # d(new_lon)/dt = (x * dy/dt - y * dx/dt) / (x^2 + y^2)
    #
    # dx/dt = -sin(lon) * d(lon)/dt
    # dy/dt = cos(lon)*cos(eps)*d(lon)/dt - sec^2(lat)*sin(eps)*d(lat)/dt
    #       = cos(lon)*cos(eps)*d(lon)/dt - sin(eps)/(cos^2(lat))*d(lat)/dt
    dx_dt = -sin_lon * lon_speed_rad
    cos_lat_sq = cos_lat * cos_lat
    if abs(cos_lat_sq) > 1e-10:
        dy_dt = (
            cos_lon * cos_eps * lon_speed_rad - (sin_eps / cos_lat_sq) * lat_speed_rad
        )
    else:
        # At poles of input coordinates
        dy_dt = cos_lon * cos_eps * lon_speed_rad

    denom = x * x + y * y
    if abs(denom) > 1e-10:
        new_lon_speed_rad = (x * dy_dt - y * dx_dt) / denom
    else:
        new_lon_speed_rad = 0.0

    # Convert speeds back to degrees/day
    new_lon_speed = math.degrees(new_lon_speed_rad)
    new_lat_speed = math.degrees(new_lat_speed_rad)

    if _split_return:
        return (new_lon, new_lat, dist), (new_lon_speed, new_lat_speed, dist_speed)
    return (new_lon, new_lat, dist, new_lon_speed, new_lat_speed, dist_speed)


def azalt(
    tjdut: float,
    flag: int,
    geopos: Tuple[float, float, float],
    atpress: float,
    attemp: float,
    xin: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """
    Convert equatorial or ecliptic coordinates to horizontal (azimuth/altitude).

    This function transforms celestial coordinates to horizontal coordinates
    (azimuth and altitude) for a given observer location and time.
    It accounts for atmospheric refraction.

    Compatible with the reference swe.azalt() API.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        flag: Coordinate type flag:
            - SE_ECL2HOR (0): Input is ecliptic (longitude, latitude, distance)
            - SE_EQU2HOR (1): Input is equatorial (RA, Dec, distance)
        geopos: Tuple of (longitude, latitude, altitude):
            - longitude: Geographic longitude of observer in degrees (East positive)
            - latitude: Geographic latitude of observer in degrees (North positive)
            - altitude: Observer altitude above sea level in meters
        atpress: Atmospheric pressure in mbar (hPa). Use 0 for no refraction.
        attemp: Atmospheric temperature in Celsius
        xin: Tuple of (longitude/RA, latitude/Dec, distance) in degrees

    Returns:
        Tuple of (azimuth, true_altitude, apparent_altitude) where:
            - azimuth: Degrees from South, westward (0=South, 90=West, 180=North, 270=East)
            - true_altitude: Geometric altitude without refraction (degrees)
            - apparent_altitude: Altitude with atmospheric refraction applied (degrees)

    Note:
        - Reference API convention: Azimuth is measured from South, westward.
          This differs from the common convention (from North, eastward).
          To convert: az_from_north = (azimuth + 180) % 360
        - Refraction is calculated using the Bennett formula, which is accurate
          to about 0.07 arcmin for altitudes above 15°.
        - When pressure=0, no refraction is applied (apparent_alt = true_alt)
        - For objects below the horizon (negative altitude), refraction is
          extrapolated but becomes less accurate.

    Example:
        >>> from libephemeris import azalt, SE_EQU2HOR, julday
        >>> jd = julday(2024, 6, 15, 12.0)
        >>> # RA=90°, Dec=23.5°, observer at lat=41.9°N, lon=12.5°E
        >>> geopos = (12.5, 41.9, 0)  # (lon, lat, altitude)
        >>> az, alt_true, alt_app = azalt(jd, SE_EQU2HOR, geopos, 1013.25, 15, (90, 23.5, 1))
    """
    # Extract geopos components
    lon = geopos[0]
    lat = geopos[1]
    altitude = geopos[2]
    pressure = atpress
    temperature = attemp
    coord = xin

    from .time_utils import _sidtime_internal
    from .state import get_timescale

    # Get true obliquity (mean obliquity + nutation in obliquity)
    # This is needed for accurate ecliptic to equatorial conversion
    ts = get_timescale()
    t = ts.ut1_jd(tjdut)
    jd_tt = t.tt

    # Mean obliquity IAU 2006 (Hilton et al. 2006, via pyerfa)
    eps0_rad = erfa.obl06(2451545.0, jd_tt - 2451545.0)
    eps0 = math.degrees(eps0_rad)

    # Nutation IAU 2006/2000A via pyerfa (~0.01-0.05 mas precision)
    dpsi_rad, deps_rad = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
    deps_deg = math.degrees(deps_rad)
    dpsi_deg = math.degrees(dpsi_rad)  # Nutation in longitude
    eps = eps0 + deps_deg  # True obliquity

    # Convert input coordinates to equatorial (RA, Dec) if ecliptic
    if flag == SE_ECL2HOR:
        # Input is ecliptic: convert to equatorial
        ecl_lon = coord[0]
        ecl_lat = coord[1]
        dist = coord[2]

        # Convert ecliptic to equatorial
        # cotrans convention: negative obliquity = ecliptic→equatorial
        eq_coord = cotrans((ecl_lon, ecl_lat, dist), -eps)
        ra = eq_coord[0]
        dec = eq_coord[1]
    else:
        # Input is already equatorial (RA, Dec)
        ra = coord[0]
        dec = coord[1]

    # Calculate Local Sidereal Time
    # Use apparent sidereal time (with nutation) for accuracy
    lst_hours = _sidtime_internal(tjdut, lon, eps, dpsi_deg)
    lst_deg = lst_hours * 15.0  # Convert hours to degrees

    # Calculate Hour Angle
    # H = LST - RA
    ha = (lst_deg - ra) % 360.0

    # Convert to radians for calculation
    ha_rad = math.radians(ha)
    dec_rad = math.radians(dec)
    lat_rad = math.radians(lat)

    # Calculate altitude (geometric/true)
    # sin(alt) = sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(H)
    sin_alt = math.sin(lat_rad) * math.sin(dec_rad) + math.cos(lat_rad) * math.cos(
        dec_rad
    ) * math.cos(ha_rad)
    # Clamp to avoid numerical issues with asin
    sin_alt = max(-1.0, min(1.0, sin_alt))
    alt_true = math.degrees(math.asin(sin_alt))

    # Calculate azimuth
    # tan(A) = sin(H) / (cos(H) * sin(lat) - tan(dec) * cos(lat))
    # Using atan2 for proper quadrant handling
    sin_ha = math.sin(ha_rad)
    cos_ha = math.cos(ha_rad)
    cos_dec = math.cos(dec_rad)

    # Numerator and denominator for azimuth calculation
    # Standard convention: azimuth from South, westward
    y = sin_ha * cos_dec
    x = cos_ha * cos_dec * math.sin(lat_rad) - math.sin(dec_rad) * math.cos(lat_rad)

    az_rad = math.atan2(y, x)
    azimuth = math.degrees(az_rad)

    # Normalize azimuth to 0-360 range
    azimuth = azimuth % 360.0
    if azimuth < 0:
        azimuth += 360.0

    # pressure=0 means "no refraction" (pyswisseph convention).
    # Only apply refraction when pressure is explicitly positive.

    # Calculate atmospheric refraction via ICAO ray-tracing
    if pressure > 0 and alt_true > -2.0:
        alt_apparent, _ = refrac_extended(
            alt_true,
            altitude,
            pressure,
            temperature,
            flag=SE_TRUE_TO_APP,
        )
    else:
        # No refraction correction
        alt_apparent = alt_true

    return (azimuth, alt_true, alt_apparent)


def azalt_rev(
    tjdut: float,
    flag: int,
    geopos: Tuple[float, float, float],
    azimuth: float,
    true_altitude: float,
) -> Tuple[float, float]:
    """
    Convert horizontal coordinates (azimuth/altitude) to equatorial or ecliptic.

    This function is the inverse of azalt(): it transforms horizontal coordinates
    (azimuth and true altitude) to celestial coordinates (equatorial or ecliptic)
    for a given observer location and time.

    Compatible with the reference swe.azalt_rev() API.

    Note: This function is not precisely the reverse of azalt(). If only an
    apparent altitude is available, the true altitude must first be computed
    using a refraction correction.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        flag: Output coordinate type flag:
            - SE_HOR2EQU (1): Output is equatorial (RA, Dec)
            - SE_HOR2ECL (0): Output is ecliptic (longitude, latitude)
        geopos: Tuple of (longitude, latitude, altitude):
            - longitude: Geographic longitude of observer in degrees (East positive)
            - latitude: Geographic latitude of observer in degrees (North positive)
            - altitude: Observer altitude above sea level in meters
        azimuth: Azimuth in degrees from South, westward (0 = South, 90 = West, etc.)
        true_altitude: True (geometric) altitude above horizon in degrees

    Returns:
        Tuple of (x1, x2) where:
            - If flag == SE_HOR2EQU: (Right Ascension, Declination) in degrees
            - If flag == SE_HOR2ECL: (Ecliptic longitude, Ecliptic latitude) in degrees

    Example:
        >>> from libephemeris import azalt_rev, SE_HOR2EQU, julday
        >>> jd = julday(2024, 6, 15, 12.0)
        >>> # Object at azimuth 90° (West), altitude 45°, observer at Rome
        >>> geopos = (12.5, 41.9, 0)  # (lon, lat, altitude)
        >>> ra, dec = azalt_rev(jd, SE_HOR2EQU, geopos, 90.0, 45.0)
    """
    # Extract geopos components
    lon = geopos[0]
    lat = geopos[1]
    altitude = geopos[2]

    from .time_utils import _sidtime_internal
    from .state import get_timescale

    # Get true obliquity (mean obliquity + nutation in obliquity)
    # This is needed for accurate equatorial to ecliptic conversion
    ts = get_timescale()
    t = ts.ut1_jd(tjdut)
    jd_tt = t.tt

    # Mean obliquity IAU 2006 (Hilton et al. 2006, via pyerfa)
    eps0_rad = erfa.obl06(2451545.0, jd_tt - 2451545.0)
    eps0 = math.degrees(eps0_rad)

    # Nutation IAU 2006/2000A via pyerfa (~0.01-0.05 mas precision)
    dpsi_rad, deps_rad = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
    deps_deg = math.degrees(deps_rad)
    dpsi_deg = math.degrees(dpsi_rad)  # Nutation in longitude
    eps = eps0 + deps_deg  # True obliquity

    # Convert to radians for calculation
    az_rad = math.radians(azimuth)
    alt_rad = math.radians(true_altitude)
    lat_rad = math.radians(lat)

    # Calculate declination from horizontal coordinates
    # For azimuth measured from South (westward), the formula is:
    # sin(dec) = sin(lat) * sin(alt) - cos(lat) * cos(alt) * cos(az)
    # Note: The minus sign is because azimuth is from South, not North
    sin_dec = math.sin(lat_rad) * math.sin(alt_rad) - math.cos(lat_rad) * math.cos(
        alt_rad
    ) * math.cos(az_rad)
    # Clamp to avoid numerical issues with asin
    sin_dec = max(-1.0, min(1.0, sin_dec))
    dec_rad = math.asin(sin_dec)
    dec = math.degrees(dec_rad)

    # Calculate hour angle from horizontal coordinates
    # From the forward azalt transform, we have:
    # y = sin(H) * cos(dec)  -> used in atan2 to get azimuth
    # x = cos(H) * cos(dec) * sin(lat) - sin(dec) * cos(lat)
    #
    # Therefore:
    # sin(H) = sin(az) * cos(alt) / cos(dec)
    # cos(H) = (cos(az) * cos(alt) + sin(dec) * cos(lat)) / (cos(dec) * sin(lat))
    cos_dec = math.cos(dec_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)

    if abs(cos_dec) < 1e-10 or abs(sin_lat) < 1e-10:
        # At poles or object at celestial pole, hour angle is undefined
        # Use 0 as default
        ha = 0.0
    else:
        sin_H = math.sin(az_rad) * math.cos(alt_rad) / cos_dec
        cos_H = (math.cos(az_rad) * math.cos(alt_rad) + sin_dec * cos_lat) / (
            cos_dec * sin_lat
        )

        ha_rad = math.atan2(sin_H, cos_H)
        ha = math.degrees(ha_rad)

    # Calculate Local Sidereal Time
    # Use apparent sidereal time (with nutation) for accuracy
    lst_hours = _sidtime_internal(tjdut, lon, eps, dpsi_deg)
    lst_deg = lst_hours * 15.0  # Convert hours to degrees

    # Calculate Right Ascension: RA = LST - H
    ra = (lst_deg - ha) % 360.0

    if flag == SE_HOR2EQU:
        # Return equatorial coordinates (RA, Dec)
        return (ra, dec)
    else:
        # SE_HOR2ECL: Convert equatorial to ecliptic
        # cotrans with positive obliquity converts equatorial to ecliptic
        ecl_coord = cotrans((ra, dec, 1.0), eps)
        return (ecl_coord[0], ecl_coord[1])


def refrac(
    alt: float,
    atpress: float = 1013.25,
    attemp: float = 15.0,
    flag: int = SE_TRUE_TO_APP,
) -> float:
    """
    Calculate true altitude from apparent altitude, or vice-versa.

    Atmospheric refraction makes celestial objects appear higher than their
    true (geometric) position. The effect is strongest near the horizon
    (about 34 arcminutes at 0 degrees) and negligible at high altitudes.

    This function converts between true (geometric) altitude and apparent
    (observed) altitude by adding or removing the refraction correction.

    Compatible with the reference swe.refrac() API.

    Args:
        alt: Altitude in degrees. For SE_TRUE_TO_APP, this is the true
                  (geometric) altitude. For SE_APP_TO_TRUE, this is the apparent
                  (observed) altitude.
        atpress: Atmospheric pressure in mbar (hPa). Default is 1013.25 (sea level).
                  Use 0 to disable refraction correction (returns input altitude).
        attemp: Atmospheric temperature in Celsius. Default is 15.0.
        flag: Direction of conversion:
            - SE_TRUE_TO_APP (0): Convert true altitude to apparent altitude
              (add refraction - object appears higher)
            - SE_APP_TO_TRUE (1): Convert apparent altitude to true altitude
              (subtract refraction - object's true position)

    Returns:
        The converted altitude in degrees:
        - For SE_TRUE_TO_APP: apparent altitude (= true altitude + refraction)
        - For SE_APP_TO_TRUE: true altitude (= apparent altitude - refraction)

    Notes:
        - At the horizon (0 degrees), refraction is approximately 34 arcminutes
          (0.567 degrees) under standard atmospheric conditions.
        - The formula used is the Saemund-Bennett formula, which is accurate
          to about 0.07 arcminutes for altitudes above 15 degrees.
        - For altitudes below -2 degrees, refraction is extrapolated linearly.
        - Pressure and temperature adjustments follow the standard formula:
          correction = (pressure / 1010) * (283 / (273 + temperature))

    Examples:
        >>> # True altitude at horizon -> apparent altitude is higher
        >>> refrac(0.0, 1013.25, 15.0, SE_TRUE_TO_APP)
        0.476...  # apparent altitude
        >>> # Apparent altitude at horizon -> true altitude (0° apparent = ~-0.5° true)
        >>> refrac(0.5, 1013.25, 15.0, SE_APP_TO_TRUE)
        0.0...  # approximately 0° true altitude
        >>> # No refraction when pressure is 0
        >>> refrac(10.0, 0, 15.0, SE_TRUE_TO_APP)
        10.0  # returns input altitude unchanged
    """
    from .refraction import calc_refraction_true_to_app, calc_refraction_app_to_true

    # No refraction if pressure is zero or negative
    if atpress <= 0:
        return alt

    if flag == SE_TRUE_TO_APP:
        refr = calc_refraction_true_to_app(alt, atpress, attemp)
        apparent = alt + refr
        # Match reference API behaviour: if refraction can't bring
        # the object above 0°, return input unchanged.
        if apparent < 0:
            return alt
        return apparent

    else:
        # SE_APP_TO_TRUE
        # Compute the refraction at the geometric horizon to establish a
        # threshold: apparent altitudes below this value are returned
        # unchanged (the object is geometrically below the horizon).
        horizon_refr = calc_refraction_true_to_app(0.0, atpress, attemp)
        if alt < horizon_refr:
            return alt

        refr = calc_refraction_app_to_true(alt, atpress, attemp)
        true_alt = alt - refr
        return true_alt


def refrac_extended(
    alt: float,
    geoalt: float,
    atpress: float = 1013.25,
    attemp: float = 15.0,
    lapserate: Optional[float] = None,
    flag: int = SE_TRUE_TO_APP,
) -> Tuple[float, Tuple[float, float, float, float]]:
    """
    Calculate true altitude from apparent altitude, or vice-versa (extended).

    This is an extended version of refrac() that includes:
    - Observer altitude above sea level
    - Atmospheric lapse rate (temperature variation with altitude)
    - Dip of the horizon calculation

    Compatible with the reference swe.refrac_extended() API.

    Args:
        alt: Altitude of object above geometric horizon in degrees.
                  For SE_TRUE_TO_APP, this is the true (geometric) altitude.
                  For SE_APP_TO_TRUE, this is the apparent (observed) altitude.
        geoalt: Altitude of observer above sea level in meters.
        atpress: Atmospheric pressure in mbar (hPa). Default is 1013.25 (sea level).
                  Use 0 to disable refraction correction.
        attemp: Atmospheric temperature at observer in degrees Celsius.
                     Default is 15.0.
        lapserate: Temperature lapse rate dT/dh in degrees Kelvin per meter.
                    Default is None (uses global value from set_lapse_rate(),
                    or 0.0065 K/m if not set).
                    Typical values range from 0.0034 to 0.010 K/m.
        flag: Direction of conversion:
            - SE_TRUE_TO_APP (0): Convert true altitude to apparent altitude
            - SE_APP_TO_TRUE (1): Convert apparent altitude to true altitude

    Returns:
        A tuple of (converted_altitude, details) where:
        - converted_altitude: The converted altitude in degrees
        - details: A tuple of 4 floats:
            - [0]: True altitude (input or computed)
            - [1]: Apparent altitude (input or computed)
            - [2]: Refraction amount in degrees
            - [3]: Dip of the horizon in degrees (negative, as horizon dips below
                   geometric horizontal for elevated observers)

    Notes:
        - The dip of the horizon accounts for both geometric effects and
          atmospheric refraction. An elevated observer sees the horizon
          below the geometric horizontal plane.
        - The lapse rate affects the dip calculation through atmospheric
          refraction near the horizon. Higher lapse rates result in less
          refraction of the horizon, yielding more negative dip values.
        - Standard atmospheric lapse rate is 0.0065 K/m (6.5°C per 1000m).

    Examples:
        >>> # At sea level, no dip of horizon
        >>> alt, (true, app, ref, dip) = refrac_extended(0.0, 0.0)
        >>> round(ref, 4)
        0.4721
        >>> round(dip, 4)
        0.0

        >>> # At 1000m elevation, horizon dips about 0.88 degrees
        >>> alt, (true, app, ref, dip) = refrac_extended(0.0, 1000.0)
        >>> round(dip, 2)
        -0.88
    """
    from .state import get_lapse_rate
    from .refraction import (
        calc_refraction_true_to_app,
        calc_refraction_app_to_true,
        calc_dip,
    )

    # Use global lapse rate if none provided
    if lapserate is None:
        lapserate = get_lapse_rate()

    # Compute refraction via ICAO ray-tracing (observer altitude aware)
    if flag == SE_TRUE_TO_APP:
        true_alt = alt
        refraction = calc_refraction_true_to_app(
            alt, atpress, attemp, geoalt, lapserate
        )
        apparent_alt = true_alt + refraction
    else:
        apparent_alt = alt
        refraction = calc_refraction_app_to_true(
            alt, atpress, attemp, geoalt, lapserate
        )
        true_alt = apparent_alt - refraction

    # Dip of the horizon for elevated observers
    dip = calc_dip(geoalt, lapserate)

    # Return the converted altitude and detail tuple
    if flag == SE_TRUE_TO_APP:
        return (apparent_alt, (true_alt, apparent_alt, refraction, dip))
    else:
        return (true_alt, (true_alt, apparent_alt, refraction, dip))


def cotrans(
    coord: Tuple[float, float, float], eps: float
) -> Tuple[float, float, float]:
    """
    Transform coordinates between ecliptic and equatorial systems.

    Compatible with the reference swe.cotrans() API.

    The direction of transformation depends on the sign of obliquity:
    - Negative obliquity: ecliptic (lon, lat) → equatorial (RA, Dec)
    - Positive obliquity: equatorial (RA, Dec) → ecliptic (lon, lat)

    Args:
        coord: Tuple of (longitude/RA, latitude/Dec, distance) in degrees
        eps: Obliquity of the ecliptic in degrees.
                   Negative for ecliptic→equatorial, positive for equatorial→ecliptic.

    Returns:
        Tuple of (transformed_lon/RA, transformed_lat/Dec, distance)
        Distance is unchanged by the transformation.

    Examples:
        >>> # Ecliptic to equatorial (negative obliquity)
        >>> cotrans((90.0, 0.0, 1.0), -23.4)
        (90.0, 23.4, 1.0)
        >>> # Equatorial to ecliptic (positive obliquity)
        >>> cotrans((90.0, 23.4, 1.0), 23.4)
        (90.0, 0.0, 1.0)
    """
    lon = coord[0]
    lat = coord[1]
    dist = coord[2]

    # Convert to radians
    # Negate obliquity to match the pyswisseph API convention
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)
    eps_rad = math.radians(-eps)

    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Calculate the new latitude (Dec for ecl→eq, β for eq→ecl)
    # sin(new_lat) = sin(lat) * cos(eps) + cos(lat) * sin(eps) * sin(lon)
    sin_new_lat = sin_lat * cos_eps + cos_lat * sin_eps * sin_lon
    # Clamp to avoid numerical issues with asin
    sin_new_lat = max(-1.0, min(1.0, sin_new_lat))
    new_lat_rad = math.asin(sin_new_lat)

    # Calculate the new longitude (RA for ecl→eq, λ for eq→ecl)
    # tan(new_lon) = (sin(lon) * cos(eps) - tan(lat) * sin(eps)) / cos(lon)
    tan_lat = math.tan(lat_rad)
    y = sin_lon * cos_eps - tan_lat * sin_eps
    x = cos_lon

    new_lon_rad = math.atan2(y, x)

    # Convert back to degrees
    new_lon = math.degrees(new_lon_rad)
    new_lat = math.degrees(new_lat_rad)

    # Normalize longitude to [0, 360)
    new_lon = new_lon % 360.0

    return (new_lon, new_lat, dist)


def degnorm(x: float) -> float:
    """
    Normalize an angle to the range [0, 360).

    Compatible with the reference swe.degnorm() API.
    Equivalent to `angle % 360` but correctly handles negative numbers.

    Args:
        x: Angle in degrees (any value)

    Returns:
        Normalized angle in range [0, 360)

    Examples:
        >>> degnorm(45)
        45.0
        >>> degnorm(-45)
        315.0
        >>> degnorm(370)
        10.0
        >>> degnorm(-370)
        350.0
        >>> degnorm(360)
        0.0
        >>> degnorm(720)
        0.0
    """
    return x % 360.0


TWO_PI = 2.0 * math.pi


def radnorm(x: float) -> float:
    """
    Normalize an angle to the range [0, 2*pi).

    Compatible with the reference swe.radnorm() API.
    Equivalent to `angle % (2*pi)` but correctly handles negative numbers.

    Args:
        x: Angle in radians (any value)

    Returns:
        Normalized angle in range [0, 2*pi)

    Examples:
        >>> import math
        >>> radnorm(math.pi / 4)  # 45 degrees
        0.7853981633974483
        >>> radnorm(-math.pi / 4)  # -45 degrees -> 315 degrees
        5.497787143782138
        >>> radnorm(3 * math.pi)  # 540 degrees -> 180 degrees
        3.141592653589793
        >>> radnorm(-3 * math.pi)  # -540 degrees -> 180 degrees
        3.141592653589793
        >>> radnorm(2 * math.pi)  # 360 degrees -> 0
        0.0
        >>> radnorm(4 * math.pi)  # 720 degrees -> 0
        0.0
    """
    return x % TWO_PI


def difdeg2n(p1: float, p2: float) -> float:
    """
    Calculate distance in degrees p1 - p2 normalized to [-180;180].

    Compatible with the reference swe.difdeg2n() API.
    Computes the signed angular difference, handling 360° wrapping.

    Args:
        p1: First angle in degrees
        p2: Second angle in degrees

    Returns:
        Normalized difference in range [-180, 180]

    Examples:
        >>> difdeg2n(10, 20)
        -10.0
        >>> difdeg2n(350, 10)
        -20.0
        >>> difdeg2n(10, 350)
        20.0
        >>> difdeg2n(180, 0)
        -180.0
    """
    diff = (p1 - p2) % 360.0
    if diff >= 180.0:
        diff -= 360.0
    return diff


def difdegn(p1: float, p2: float) -> float:
    """
    Calculate distance in degrees p1 - p2 normalized to [0, 360).

    Compatible with the reference swe.difdegn() API.
    Computes the difference between two angles, always returning a positive
    value in the range [0, 360). Unlike difdeg2n() which returns [-180, 180],
    this function always returns a positive value.

    Args:
        p1: First angle in degrees
        p2: Second angle in degrees

    Returns:
        Normalized difference in range [0, 360)

    Examples:
        >>> difdegn(10, 20)
        350.0
        >>> difdegn(20, 10)
        10.0
        >>> difdegn(350, 10)
        340.0
        >>> difdegn(10, 350)
        20.0
        >>> difdegn(0, 0)
        0.0
    """
    return (p1 - p2) % 360.0


def difrad2n(p1: float, p2: float) -> float:
    """
    Calculate distance in radians p1 - p2 normalized to [-π, π].

    This is the radians equivalent of difdeg2n().
    Computes the signed angular difference, handling 2π wrapping.

    Args:
        p1: First angle in radians
        p2: Second angle in radians

    Returns:
        Normalized difference in range [-π, π]

    Examples:
        >>> import math
        >>> difrad2n(0.1, 0.2)
        -0.1
        >>> difrad2n(6.0, 0.2)  # Near 2π wraparound
        -0.48318530717958...
        >>> difrad2n(0.2, 6.0)  # Opposite direction
        0.48318530717958...
        >>> difrad2n(math.pi, 0)
        -3.141592653589793
    """
    diff = (p1 - p2) % TWO_PI
    if diff >= math.pi:
        diff -= TWO_PI
    return diff


# Constants for centiseconds calculations
# 1 centisecond = 1/100 arcsecond
# 360 degrees = 360 * 3600 * 100 = 129,600,000 centiseconds
CS360 = 360 * 3600 * 100  # 129600000 centiseconds in a full circle
CS180 = 180 * 3600 * 100  # 64800000 centiseconds in a half circle


def difcs2n(p1: int, p2: int) -> int:
    """
    Calculate distance in centiseconds p1 - p2 normalized to [-180°, +180°].

    This function computes the signed angular difference between two angles
    expressed in centiseconds (1/100 of an arcsecond), handling 360° wrapping.
    The result is normalized to the equivalent of [-180°, +180°] in centiseconds.

    Compatible with the reference swe.difcs2n() API.

    Args:
        p1: First angle in centiseconds
        p2: Second angle in centiseconds

    Returns:
        Normalized difference in range [-64800000, 64800000) centiseconds
        (equivalent to [-180°, +180°))

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 360° = 129,600,000 centiseconds
        - 180° = 64,800,000 centiseconds
        - When the difference is exactly 180°, returns -180° (matching the reference API)

    Examples:
        >>> difcs2n(360000, 720000)  # 1° - 2° = -1° = -360000 cs
        -360000
        >>> difcs2n(129240000, 360000)  # 359° - 1° = -2° = -720000 cs (shorter path)
        -720000
        >>> difcs2n(360000, 129240000)  # 1° - 359° = 2° = 720000 cs (shorter path)
        720000
        >>> difcs2n(64800000, 0)  # 180° - 0° = -180° (reference API convention)
        -64800000
    """
    diff = (p1 - p2) % CS360
    if diff >= CS180:
        diff -= CS360
    return diff


def difcsn(p1: int, p2: int) -> int:
    """
    Calculate distance in centiseconds p1 - p2 normalized to [0, 360°).

    This function computes the angular difference between two angles
    expressed in centiseconds (1/100 of an arcsecond), handling 360° wrapping.
    The result is normalized to the equivalent of [0, 360°) in centiseconds,
    always returning a non-negative value.

    Compatible with the reference swe.difcsn() API.

    Args:
        p1: First angle in centiseconds
        p2: Second angle in centiseconds

    Returns:
        Normalized difference in range [0, 129600000) centiseconds
        (equivalent to [0°, 360°))

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 360° = 129,600,000 centiseconds
        - Unlike difcs2n() which returns [-180°, +180°], this function
          always returns a positive value in [0°, 360°)

    Examples:
        >>> difcsn(720000, 360000)  # 2° - 1° = 1° = 360000 cs
        360000
        >>> difcsn(360000, 720000)  # 1° - 2° = 359° = 129240000 cs (positive)
        129240000
        >>> difcsn(360000, 129240000)  # 1° - 359° = 2° = 720000 cs
        720000
        >>> difcsn(0, 0)  # 0° - 0° = 0°
        0
    """
    return (p1 - p2) % CS360


def csnorm(cs: int) -> int:
    """
    Normalize a value in centiseconds to the range [0, 360°).

    This function normalizes an angle expressed in centiseconds (1/100 of an
    arcsecond) to the equivalent of [0°, 360°) in centiseconds, i.e., the range
    [0, 129600000). Correctly handles negative numbers.

    Compatible with the reference swe.csnorm() API.

    Args:
        cs: Angle in centiseconds (any value)

    Returns:
        Normalized angle in range [0, 129600000) centiseconds
        (equivalent to [0°, 360°))

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 360° = 129,600,000 centiseconds
        - This is the centisecond equivalent of degnorm() for degrees
          and radnorm() for radians

    Examples:
        >>> csnorm(360000)  # 1° stays as 1°
        360000
        >>> csnorm(-360000)  # -1° -> 359° = 129240000 cs
        129240000
        >>> csnorm(129960000)  # 361° -> 1° = 360000 cs
        360000
        >>> csnorm(-129960000)  # -361° -> 359° = 129240000 cs
        129240000
        >>> csnorm(129600000)  # 360° -> 0°
        0
        >>> csnorm(0)
        0
    """
    return cs % CS360


def csroundsec(cs: int) -> int:
    """
    Round a value in centiseconds to the nearest arcsecond.

    This function rounds a centisecond value to the nearest arcsecond,
    returning the result still in centiseconds (i.e., rounded to the
    nearest multiple of 100).

    Compatible with the reference swe.csroundsec() API.

    Args:
        cs: Angle in centiseconds (1/100 arcsecond)

    Returns:
        Centisecond value rounded to the nearest arcsecond (multiple of 100)

    Notes:
        - 1 centisecond = 1/100 arcsecond
        - Uses truncation-toward-zero integer division with +50 offset
        - Special handling at degree boundaries (30° for positive, 90° for negative)

    Examples:
        >>> csroundsec(150)  # 1.50 arcseconds -> 2 arcseconds = 200 cs
        200
        >>> csroundsec(149)  # 1.49 arcseconds -> 1 arcsecond = 100 cs
        100
        >>> csroundsec(50)   # 0.50 arcseconds -> 1 arcsecond = 100 cs
        100
        >>> csroundsec(100)  # 1.00 arcsecond -> 1 arcsecond = 100 cs
        100
        >>> csroundsec(0)    # 0 centiseconds -> 0 arcseconds
        0
        >>> csroundsec(-150) # -1.50 arcseconds -> -1 arcsecond = -100 cs
        -100
    """
    # Truncation-toward-zero integer division (not floor division):
    # For cs + 50, we need truncation toward zero, not floor division
    cs_plus_50 = cs + 50
    if cs_plus_50 >= 0:
        result = (cs_plus_50 // 100) * 100
    else:
        # Truncation toward zero for negative values
        result = -((-cs_plus_50) // 100) * 100

    # Boundary correction: prevent rounding up across zodiacal sign boundaries
    # (30° = 10,800,000 cs). When a value just below a sign boundary would round
    # up to the boundary, keep it in the current sign instead. This is
    # mathematically correct for angular display — rounding 29°59'59.7" should
    # yield 30°00'00" only if the value is >= the midpoint of the last bin.
    if cs > 0 and result % 10800000 == 0 and result != 0 and cs < result:
        return result - 100

    # Boundary correction for negative values at quadrant boundaries
    # (90° = 32,400,000 cs). Same principle: prevent rounding across -90°/-180° etc.
    if cs < 0:
        if result != 0 and result % 32400000 == 0 and cs <= result - 100:
            return result - 100
        # Negative values in (-100, 0) round to 0; values <= -100 round to -100
        # when the standard formula yields 0 (truncation toward zero artifact)
        if cs <= -100 and result == 0:
            return -100

    return result


def cs2degstr(cs: int) -> str:
    """
    Convert a value in centiseconds to a formatted degrees string.

    This function converts an angular measurement in centiseconds (1/100 of an
    arcsecond) to a human-readable string in the format "D°M'S.ss\"" where D is
    degrees, M is arcminutes, S is arcseconds, and ss is centiseconds (hundredths
    of arcsecond).

    Compatible with the reference swe.cs2degstr() API.

    Args:
        cs: Angle in centiseconds (any integer value)

    Returns:
        Formatted string representing the angle in degrees, minutes, seconds.
        Format: "DDD°MM'SS.ss\"" (e.g., "123°45'06.78\"")

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 360° = 129,600,000 centiseconds
        - Negative values produce negative degree strings (e.g., "-1°00'00.00\"")
        - The seconds field includes two decimal places for centiseconds

    Examples:
        >>> cs2degstr(0)
        '  0° 0\\' 0.00"'
        >>> cs2degstr(360000)  # 1 degree
        '  1° 0\\' 0.00"'
        >>> cs2degstr(3723456)  # 10°20'34.56"
        ' 10°20\\'34.56"'
        >>> cs2degstr(-360000)  # -1 degree
        ' -1° 0\\' 0.00"'
    """
    # Handle sign
    if cs < 0:
        sign = -1
        cs = -cs
    else:
        sign = 1

    # Extract degrees, minutes, seconds, and centiseconds
    # 1 degree = 3600 * 100 = 360000 centiseconds
    # 1 minute = 60 * 100 = 6000 centiseconds
    # 1 second = 100 centiseconds
    degrees = cs // 360000
    remainder = cs % 360000
    minutes = remainder // 6000
    remainder = remainder % 6000
    seconds = remainder // 100
    centisecs = remainder % 100

    # Apply sign to degrees
    if sign < 0:
        degrees = -degrees

    # Format the string matching reference API format
    # Format: "%3d°%2d'%2d.%02d"" with proper spacing
    return f"{degrees:3d}°{minutes:2d}'{seconds:2d}.{centisecs:02d}\""


def cs2lonlatstr(cs: int, plus: "str | bytes", minus: "str | bytes") -> str:
    """
    Convert a value in centiseconds to a formatted longitude/latitude string.

    This function converts an angular measurement in centiseconds (1/100 of an
    arcsecond) to a human-readable string with a directional character.

    Compatible with the reference swe.cs2lonlatstr() API.

    Args:
        cs: Angle in centiseconds (any integer value)
        plus: Character to use for positive values (e.g., b"N" or b"E")
        minus: Character to use for negative values (e.g., b"S" or b"W")

    Returns:
        Formatted string representing the angle with directional character.
        Format: "{deg}{char}{min:02d}" when seconds == 0,
        or "{deg}{char}{min:02d}'{sec:02d}" when seconds > 0.

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - Positive values use plus, negative values use minus
        - Seconds are rounded to whole numbers (no centiseconds shown)
        - If seconds round to 0, the seconds portion is omitted

    Examples:
        >>> cs2lonlatstr(360000, b"N", b"S")
        '1N00'
        >>> cs2lonlatstr(-360000, b"N", b"S")
        '1S00'
        >>> cs2lonlatstr(12345678, b"E", b"W")
        "34E17'37"
    """
    # Decode bytes to str if needed
    if isinstance(plus, bytes):
        plus = plus.decode("ascii")
    if isinstance(minus, bytes):
        minus = minus.decode("ascii")

    # Determine direction character based on sign
    if cs >= 0:
        direction = plus
    else:
        direction = minus
        cs = -cs

    # Extract degrees, minutes, and seconds
    # 1 degree = 3600 * 100 = 360000 centiseconds
    # 1 minute = 60 * 100 = 6000 centiseconds
    # 1 second = 100 centiseconds
    degrees = cs // 360000
    remainder = cs % 360000
    minutes = remainder // 6000
    remainder = remainder % 6000
    # Round centiseconds to nearest arcsecond
    seconds = (remainder + 50) // 100

    # Handle carry-over from rounding
    if seconds >= 60:
        seconds = 0
        minutes += 1
        if minutes >= 60:
            minutes = 0
            degrees += 1

    # Format matching pyswisseph: "{deg}{char}{min:02d}" or
    # "{deg}{char}{min:02d}'{sec:02d}"
    if seconds == 0:
        return f"{degrees}{direction}{minutes:02d}"
    else:
        return f"{degrees}{direction}{minutes:02d}'{seconds:02d}"


def cs2timestr(cs: int, sep: "str | bytes" = ":", suppresszero: bool = False) -> str:
    """
    Convert a value in centiseconds to a formatted time string.

    This function converts a time measurement in centiseconds (1/100 of a second)
    to a human-readable string in the format "HH:MM:SS" where HH is hours, MM is
    minutes, and SS is seconds.

    Compatible with the reference swe.cs2timestr() API.

    Args:
        cs: Time in centiseconds (any integer value)

    Returns:
        Formatted string representing the time in hours, minutes, seconds.
        Format: "HH:MM:SS" (e.g., "12:34:56")

    Notes:
        - 1 centisecond = 1/100 second
        - 1 second = 100 centiseconds
        - 1 minute = 6000 centiseconds
        - 1 hour = 360000 centiseconds
        - Negative values produce negative hour strings (e.g., "-1:00:00")
        - Seconds are rounded to whole numbers (centiseconds are rounded)

    Examples:
        >>> cs2timestr(0)
        ' 0:00:00'
        >>> cs2timestr(360000)  # 1 hour
        ' 1:00:00'
        >>> cs2timestr(4526050)  # 12:34:21 (with rounding from .50)
        '12:34:21'
        >>> cs2timestr(-360000)  # -1 hour
        '-1:00:00'
    """
    # Accept bytes separator (pyswisseph uses b':')
    if isinstance(sep, bytes):
        sep = sep.decode("ascii")

    # Handle sign
    if cs < 0:
        sign = -1
        cs = -cs
    else:
        sign = 1

    # Extract hours, minutes, and seconds
    # 1 hour = 60 * 60 * 100 = 360000 centiseconds
    # 1 minute = 60 * 100 = 6000 centiseconds
    # 1 second = 100 centiseconds
    hours = cs // 360000
    remainder = cs % 360000
    minutes = remainder // 6000
    remainder = remainder % 6000
    # Round centiseconds to nearest second
    seconds = (remainder + 50) // 100

    # Handle carry-over from rounding
    if seconds >= 60:
        seconds = 0
        minutes += 1
        if minutes >= 60:
            minutes = 0
            hours += 1

    # Apply sign to hours
    if sign < 0:
        hours = -hours

    # Format the string matching reference API format
    # Format: "%02d<sep>%02d<sep>%02d"
    # Hours are mod 24 to match pyswisseph behavior
    hours = hours % 24
    if suppresszero:
        if seconds == 0:
            return f"{hours:02d}{sep}{minutes:02d}"
    return f"{hours:02d}{sep}{minutes:02d}{sep}{seconds:02d}"


def deg_midp(x1: float, x2: float) -> float:
    """
    Calculate the midpoint between two angles in degrees.

    Handles wraparound at 360° correctly by finding the midpoint along the
    shorter arc between the two angles. For example, the midpoint between
    350° and 10° is 0° (or equivalently 360°), not 180°.

    Args:
        x1: First angle in degrees (any value, will be normalized)
        x2: Second angle in degrees (any value, will be normalized)

    Returns:
        Midpoint angle in range [0, 360)

    Examples:
        >>> deg_midp(0, 90)
        45.0
        >>> deg_midp(350, 10)
        0.0
        >>> deg_midp(10, 350)
        0.0
        >>> deg_midp(180, 0)
        270.0
        >>> deg_midp(170, 190)
        180.0
        >>> deg_midp(-10, 10)
        0.0
    """
    # Normalize both angles to [0, 360)
    x1 = x1 % 360.0
    x2 = x2 % 360.0

    # Calculate the difference
    diff = x2 - x1

    # When both arcs are equally long (diff exactly ±180) we follow the
    # pyswisseph convention: always take the positive (clockwise) half,
    # i.e. treat -180 the same as +180.
    if diff > 180.0:
        diff -= 360.0
    elif diff < -180.0:
        diff += 360.0
    elif diff == -180.0:
        diff = 180.0

    # Calculate midpoint along the chosen arc
    midp = x1 + diff / 2.0

    # Normalize result to [0, 360)
    return midp % 360.0


def rad_midp(x: float, y: float) -> float:
    """
    Calculate the midpoint between two angles in radians.

    Handles wraparound at 2*pi correctly by finding the midpoint along the
    shorter arc between the two angles. For example, the midpoint between
    5.5 radians and 0.5 radians is 0 (or equivalently 2*pi), not pi.

    Args:
        x: First angle in radians (any value, will be normalized)
        y: Second angle in radians (any value, will be normalized)

    Returns:
        Midpoint angle in range [0, 2*pi)

    Examples:
        >>> import math
        >>> rad_midp(0, math.pi / 2)  # 0 to 90 degrees -> 45 degrees
        0.7853981633974483
        >>> rad_midp(5.5, 0.5)  # Near 2*pi wraparound
        0.0  # approximately
        >>> rad_midp(math.pi, 0)  # 180 to 0 -> 270 degrees
        4.71238898038469
    """
    # Normalize both angles to [0, 2*pi)
    x = x % TWO_PI
    y = y % TWO_PI

    # Calculate the difference
    diff = y - x

    # When both arcs are equally long (diff exactly ±π) we follow the
    # pyswisseph convention: always take the positive (clockwise) half.
    if diff > math.pi:
        diff -= TWO_PI
    elif diff < -math.pi:
        diff += TWO_PI
    elif diff == -math.pi:
        diff = math.pi

    # Calculate midpoint along the chosen arc
    midp = x + diff / 2.0

    # Normalize result to [0, 2*pi)
    return midp % TWO_PI


def d2l(d: float) -> int:
    """
    Convert a double (float) to a long integer with rounding.

    This function rounds a floating-point number to the nearest integer using
    "round half away from zero" semantics (also known as commercial rounding).
    This is standard practice in astronomical computation software.

    Compatible with the pyswisseph swe.d2l() API.

    Args:
        d: A floating-point number to convert.

    Returns:
        The nearest integer to d. For values exactly halfway between two
        integers, rounds away from zero (e.g., 0.5 -> 1, -0.5 -> -1).

    Notes:
        - This differs from Python's built-in round() function, which uses
          "round half to even" (banker's rounding) for Python 3.
        - Used internally for coordinate conversions and also exposed
          publicly for consistency with the pyswisseph API.

    Examples:
        >>> d2l(1.4)
        1
        >>> d2l(1.5)
        2
        >>> d2l(1.6)
        2
        >>> d2l(-1.4)
        -1
        >>> d2l(-1.5)
        -2
        >>> d2l(-1.6)
        -2
        >>> d2l(0.5)
        1
        >>> d2l(-0.5)
        -1
        >>> d2l(2.5)
        3
        >>> d2l(-2.5)
        -3
    """
    if d >= 0:
        return int(d + 0.5)
    else:
        return int(d - 0.5)


# Split degree flags (imported from constants for convenience)
SPLIT_DEG_ROUND_SEC: int = 1  # Round to seconds
SPLIT_DEG_ROUND_MIN: int = 2  # Round to minutes
SPLIT_DEG_ROUND_DEG: int = 4  # Round to degrees
SPLIT_DEG_ZODIACAL: int = 8  # Return zodiac sign number (0-11)
SPLIT_DEG_NAKSHATRA: int = 1024  # Return nakshatra number (0-26)
SPLIT_DEG_KEEP_SIGN: int = 16  # Don't round to next zodiac sign/nakshatra
SPLIT_DEG_KEEP_DEG: int = 32  # Don't round to next degree


def _rounding_offset(roundflag: int) -> float:
    """Compute the rounding offset based on which rounding flag is active.

    Only one rounding level applies at a time, checked in order of
    coarseness: degrees first, then minutes, then seconds.
    """
    if roundflag & SPLIT_DEG_ROUND_DEG:
        return 0.5
    if roundflag & SPLIT_DEG_ROUND_MIN:
        return 0.5 / 60.0
    if roundflag & SPLIT_DEG_ROUND_SEC:
        return 0.5 / 3600.0
    return 0.0


def _decompose_to_dms(ddeg: float, has_rounding: bool) -> Tuple[int, int, int, float]:
    """Break a non-negative decimal degree value into deg, min, sec, secfr.

    Uses successive subtraction: each component is extracted as an integer
    then removed from the running remainder before extracting the next one.

    When *has_rounding* is ``True`` the sub-second fraction ``secfr`` is
    set to the integer seconds value (as a float) rather than the true
    fractional part.
    """
    ideg = int(ddeg)
    ddeg = max(ddeg - ideg, 0.0)
    imin = int(ddeg * 60.0)
    ddeg = max(ddeg - imin / 60.0, 0.0)
    isec = int(ddeg * 3600.0)
    if has_rounding:
        secfr = float(isec)
    else:
        secfr = max(ddeg * 3600.0 - isec, 0.0)
    return ideg, imin, isec, secfr


def _split_deg_nakshatra(
    ddeg: float, roundflag: int
) -> Tuple[int, int, int, float, int]:
    """Nakshatra-mode split for non-negative degree values.

    The ecliptic is divided into 27 equal segments (nakshatras) of
    13°20' each.  Returns the position within the current nakshatra.
    """
    nakshatra_span = 360.0 / 27.0  # 13.33333...°

    # Reduce to within-nakshatra position for KEEP-flag boundary checks.
    # Plain fmod(ddeg, span) fails at exact nakshatra multiples (40°, 120°,
    # 360°) because span = 360/27 is irrational in IEEE 754.
    # fmod(ddeg*27, 360)/27 operates on integer multiples at those
    # boundaries, giving exact zeros.
    pos_in_nak = math.fmod(ddeg * 27.0, 360.0) / 27.0

    # Determine and conditionally suppress rounding offset
    offset = _rounding_offset(roundflag)
    if offset > 0.0:
        if roundflag & SPLIT_DEG_KEEP_DEG:
            # Suppress when offset would bump the integer-degree part
            if int(pos_in_nak + offset) - int(pos_in_nak) > 0:
                offset = 0.0
        elif roundflag & SPLIT_DEG_KEEP_SIGN:
            # Suppress when offset would cross into next nakshatra
            if pos_in_nak + offset >= nakshatra_span:
                offset = 0.0

    # Apply offset to the original degree and re-reduce using *27/27.
    # This preserves FP precision better than adding offset to the
    # already-reduced pos_in_nak (avoids 1-ULP loss at boundaries
    # like 90° + ROUND_MIN offset).
    ddeg += offset
    nak_idx = int(ddeg / nakshatra_span)
    if nak_idx >= 27:
        nak_idx = 0
    pos_in_nak = math.fmod(ddeg * 27.0, 360.0) / 27.0

    has_rounding = bool(
        roundflag & (SPLIT_DEG_ROUND_DEG | SPLIT_DEG_ROUND_MIN | SPLIT_DEG_ROUND_SEC)
    )
    ideg, imin, isec, secfr = _decompose_to_dms(pos_in_nak, has_rounding)
    return (ideg, imin, isec, secfr, nak_idx)


def split_deg(degree: float, roundflag: int = 0) -> Tuple[int, int, int, float, int]:
    """
    Split a degree value into sign/nakshatra, degrees, minutes, seconds, and fraction.

    This function decomposes an ecliptic longitude (or any angle in degrees) into
    its constituent parts: zodiac sign (or nakshatra), degrees within that sign,
    minutes, seconds, and fraction of second.

    Compatible with swe_split_deg().

    Args:
        degree: Position in decimal degrees (can be negative or > 360)
        roundflag: Bit flags combination indicating how to round and format:
            - 0: Don't round, return all components with full precision
            - SPLIT_DEG_ROUND_SEC (1): Round to nearest second
            - SPLIT_DEG_ROUND_MIN (2): Round to nearest minute
            - SPLIT_DEG_ROUND_DEG (4): Round to nearest degree
            - SPLIT_DEG_ZODIACAL (8): Return zodiac sign number (0-11, each 30 deg)
            - SPLIT_DEG_NAKSHATRA (1024): Return nakshatra number (0-26, each 13deg20')
            - SPLIT_DEG_KEEP_SIGN (16): Don't round to next zodiac sign/nakshatra
            - SPLIT_DEG_KEEP_DEG (32): Don't round to next degree

    Returns:
        Tuple of (deg, min, sec, secfr, sign) where:
        - deg: Degrees within sign (0-29 with ZODIACAL, 0-13 with NAKSHATRA,
               or total degrees without either flag)
        - min: Arc minutes (0-59)
        - sec: Arc seconds (0-59)
        - secfr: Fraction of arc second (0.0-0.999...), or the rounded seconds
                 value when rounding flags are used
        - sign: Zodiac sign (0-11), nakshatra (0-26), or +1/-1 for positive/negative

    Notes:
        - Without ZODIACAL or NAKSHATRA flag, sign returns +1 (positive) or -1 (negative)
        - With ZODIACAL flag: uses absolute value, sign is 0-11
          (0=Aries, 1=Taurus, ..., 11=Pisces)
        - With NAKSHATRA flag: uses absolute value, sign is 0-26
          (each nakshatra spans 13 deg 20 min = 13.333... degrees)
        - Rounding flags affect how values are truncated/rounded
        - KEEP_SIGN prevents rounding from advancing to the next sign
        - KEEP_DEG prevents rounding from advancing to the next degree

    Examples:
        >>> split_deg(45.5, 0)  # No flags: 45deg 30min
        (45, 30, 0, 0.0, 1)
        >>> split_deg(45.5, SPLIT_DEG_ZODIACAL)  # With zodiac: 15deg 30min Taurus
        (15, 30, 0, 0.0, 1)
        >>> split_deg(-30.5, 0)  # Negative: sign = -1
        (30, 30, 0, 0.0, -1)
        >>> split_deg(-30.5, SPLIT_DEG_ZODIACAL)  # Negative with zodiac: uses absolute
        (0, 30, 0, 0.0, 1)
    """
    # --- Negative handling and nakshatra dispatch ---
    # Negative inputs never enter nakshatra mode; they are made positive
    # and flagged with sign_out = -1.  Non-negative inputs with the
    # NAKSHATRA flag take a dedicated path.
    if degree < 0:
        sign_out = -1
        ddeg = -degree
    elif roundflag & SPLIT_DEG_NAKSHATRA:
        return _split_deg_nakshatra(degree, roundflag)
    else:
        sign_out = 1
        ddeg = degree

    # --- Rounding offset ---
    # The offset is applied to the *full* degree value before any
    # zodiacal division.  KEEP_DEG suppresses it when the integer
    # degree part would change.  KEEP_SIGN suppresses it when a
    # 30° sign boundary would be crossed.
    has_rounding = bool(
        roundflag & (SPLIT_DEG_ROUND_DEG | SPLIT_DEG_ROUND_MIN | SPLIT_DEG_ROUND_SEC)
    )
    offset = _rounding_offset(roundflag)

    if offset > 0.0:
        if roundflag & SPLIT_DEG_KEEP_DEG:
            # Would the integer part of ddeg change?
            if int(ddeg + offset) - int(ddeg) > 0:
                offset = 0.0
        elif roundflag & SPLIT_DEG_KEEP_SIGN:
            # Would we cross the next 30° sign boundary?
            if math.fmod(ddeg, 30.0) + offset >= 30.0:
                offset = 0.0

    ddeg += offset

    # --- Zodiacal sign extraction (after rounding) ---
    if roundflag & SPLIT_DEG_ZODIACAL:
        sign_out = int(ddeg / 30.0)
        if sign_out == 12:  # exactly 360° maps back to sign 0
            sign_out = 0
        ddeg = math.fmod(ddeg, 30.0)

    # --- Decompose into deg / min / sec / secfr ---
    ideg, imin, isec, secfr = _decompose_to_dms(ddeg, has_rounding)

    return (ideg, imin, isec, secfr, sign_out)


def swe_calc_angles(jd_ut: float, lat: float, lon: float):
    """
    Pre-calculate and cache astrological angles and planet positions
    for use with Arabic parts.

    Args:
        jd_ut: Julian Day (UT)
        lat: Latitude (degrees)
        lon: Longitude (degrees)

    Returns:
        Dictionary with calculated positions
    """
    from .state import set_angles_cache, set_topo
    from .angles import calc_angles
    from .planets import swe_calc_ut
    from .constants import SE_SUN, SE_MOON, SE_MERCURY, SE_VENUS

    # Set observer location
    set_topo(lon, lat, 0)

    # Calculate angles
    angles_dict = calc_angles(jd_ut, lat, lon)

    # Calculate and add planet positions for Arabic parts
    sun_pos, _ = swe_calc_ut(jd_ut, SE_SUN, 0)
    moon_pos, _ = swe_calc_ut(jd_ut, SE_MOON, 0)
    mercury_pos, _ = swe_calc_ut(jd_ut, SE_MERCURY, 0)
    venus_pos, _ = swe_calc_ut(jd_ut, SE_VENUS, 0)

    angles_dict["Sun"] = sun_pos[0]
    angles_dict["Moon"] = moon_pos[0]
    angles_dict["Mercury"] = mercury_pos[0]
    angles_dict["Venus"] = venus_pos[0]

    # Cache for Arabic parts
    set_angles_cache(angles_dict)

    return angles_dict


def angular_separation(
    lon1: float, lat1: float, lon2: float, lat2: float
) -> float:
    """Compute the great-circle angular separation between two points.

    Uses the Vincenty formula, which is numerically stable at all
    separations (including < 1 arcsecond where the spherical law of
    cosines loses precision).

    Args:
        lon1: Longitude of the first point in degrees.
        lat1: Latitude of the first point in degrees.
        lon2: Longitude of the second point in degrees.
        lat2: Latitude of the second point in degrees.

    Returns:
        Angular separation in degrees, in the range [0, 180].
    """
    lat1_rad = math.radians(lat1)
    lat2_rad = math.radians(lat2)
    dlon = math.radians(lon2 - lon1)
    sin_lat1 = math.sin(lat1_rad)
    cos_lat1 = math.cos(lat1_rad)
    sin_lat2 = math.sin(lat2_rad)
    cos_lat2 = math.cos(lat2_rad)
    sin_dlon = math.sin(dlon)
    cos_dlon = math.cos(dlon)

    num = math.sqrt(
        (cos_lat2 * sin_dlon) ** 2
        + (cos_lat1 * sin_lat2 - sin_lat1 * cos_lat2 * cos_dlon) ** 2
    )
    den = sin_lat1 * sin_lat2 + cos_lat1 * cos_lat2 * cos_dlon

    return math.degrees(math.atan2(num, den))
