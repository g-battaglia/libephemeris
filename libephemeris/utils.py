# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Utility functions for libephemeris.

Provides helper functions compatible with the reference API including
angular calculations and other mathematical utilities.

Provenance:
    Coordinate rotations, longitude normalization, angular differences, and
    vector/spherical conversions are direct mathematical identities with units
    documented per function. Nutation/obliquity and time inputs delegate to the
    registered ERFA/IERS/Vondrak pipeline; refraction delegates to
    ``refraction.py``. API-compatible names and edge-case policy are project
    presentation choices, not sources of numerical formulas.
"""

from __future__ import annotations

import math
import erfa
from typing import Optional, Sequence, Tuple, Union, cast, overload

# Azalt calculation method flags (compatible with the reference API)
ECL2HOR: int = 0  # Ecliptic coordinates to horizontal
EQU2HOR: int = 1  # Equatorial coordinates to horizontal

# Azalt_rev calculation method flags (compatible with the reference API)
HOR2ECL: int = 0  # Horizontal to ecliptic coordinates
HOR2EQU: int = 1  # Horizontal to equatorial coordinates

# Refraction calculation flags (compatible with the reference API)
TRUE_TO_APP: int = 0  # True altitude to apparent altitude
APP_TO_TRUE: int = 1  # Apparent altitude to true altitude


@overload
def cotrans_sp(
    coord: "Sequence[float]", eps_or_speed: float
) -> "Tuple[float, float, float, float, float, float]": ...


@overload
def cotrans_sp(
    coord: "Sequence[float]", eps_or_speed: "Sequence[float]", eps: float
) -> "Tuple[Tuple[float, float, float], Tuple[float, float, float]]": ...


def cotrans_sp(
    coord: "Sequence[float]",
    eps_or_speed: "Union[float, Sequence[float]]",
    eps: "Optional[float]" = None,
) -> "Union[Tuple[float, float, float, float, float, float], Tuple[Tuple[float, float, float], Tuple[float, float, float]]]":
    """
    Transform coordinates and velocities between ecliptic and equatorial systems.

    Supports two calling conventions:
      1. Reference-API-compatible: cotrans_sp(coord_6, eps) -> flat 6-tuple
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
        speed_seq = cast(Sequence[float], eps_or_speed)
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
        obliquity = float(cast(float, eps_or_speed))
        lon = float(coord[0])
        lat = float(coord[1])
        dist = float(coord[2])
        lon_speed = float(coord[3])
        lat_speed = float(coord[4])
        dist_speed = float(coord[5])
        _split_return = False
    # Convert to radians
    # Negate obliquity to match the reference API convention
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

    # Calculate the new longitude (RA for ecl->eq, lambda for eq->ecl) from the
    # full 3D rotation (see cotrans): using cos_lat-scaled components rather
    # than dividing by cos_lat keeps it correct for |lat| > 90°.
    y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps
    x = cos_lat * cos_lon

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

    # Derivative of new longitude, using the SAME full-3D components as the
    # position above (x = cos_lat*cos_lon, y = cos_lat*sin_lon*cos_eps -
    # sin_lat*sin_eps): d(new_lon)/dt = (x*dy/dt - y*dx/dt) / (x^2 + y^2).
    # These derivatives match the cos_lat-scaled x,y (no 1/cos_lat term), so
    # the speed stays correct for |lat| > 90° and is unchanged for |lat| <= 90°.
    dx_dt = -sin_lat * cos_lon * lat_speed_rad - cos_lat * sin_lon * lon_speed_rad
    dy_dt = (
        -cos_eps * sin_lat * sin_lon * lat_speed_rad
        + cos_eps * cos_lat * cos_lon * lon_speed_rad
        - cos_lat * sin_eps * lat_speed_rad
    )

    # denom = x^2 + y^2 = cos^2(new_lat); it is >= 0 and vanishes only at the
    # exact pole. The analytic longitude rate (x*dy - y*dx)/denom is stable for
    # every denom > 0, so no artificial guard band is used: an earlier
    # ``denom > 1e-10`` guard zeroed the rate within ~0.0006 deg of the pole,
    # whereas the compatibility contract keeps the (large, finite) analytic
    # value there.
    # Only the genuine singularity (denom exactly 0.0) falls back to 0.0.
    denom = x * x + y * y
    if denom > 0.0:
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

    Compatible with the reference azalt() API.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        flag: Coordinate type flag:
            - ECL2HOR (0): Input is ecliptic (longitude, latitude, distance)
            - EQU2HOR (1): Input is equatorial (RA, Dec, distance)
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
        - Refraction follows the clean-room closed-form compatibility model in
          ``libephemeris/refraction.py``.
        - When pressure=0, no refraction is applied (apparent_alt = true_alt)
        - For objects below the horizon (negative altitude), refraction is
          extrapolated but becomes less accurate.

    Example:
        >>> from libephemeris import azalt, EQU2HOR, julday
        >>> jd = julday(2024, 6, 15, 12.0)
        >>> # RA=90°, Dec=23.5°, observer at lat=41.9°N, lon=12.5°E
        >>> geopos = (12.5, 41.9, 0)  # (lon, lat, altitude)
        >>> az, alt_true, alt_app = azalt(jd, EQU2HOR, geopos, 1013.25, 15, (90, 23.5, 1))
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

    # Mean obliquity via the Vondrák 2011 long-term model — the same of-date
    # mean obliquity calc_ut/fast_calc use to express ecliptic-of-date
    # coordinates, so this horizon conversion shares their ecliptic frame. The
    # IAU 2006 obliquity polynomial (erfa.obl06) diverges by several arcsec at
    # deep-BCE dates and would put the ecliptic here in a different frame.
    from .precession_vondrak import vondrak_mean_obliquity_deg

    eps0 = vondrak_mean_obliquity_deg(jd_tt)

    # Nutation IAU 2006/2000A via pyerfa (~0.01-0.05 mas precision)
    dpsi_rad, deps_rad = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
    deps_deg = math.degrees(deps_rad)
    dpsi_deg = math.degrees(dpsi_rad)  # Nutation in longitude
    eps = eps0 + deps_deg  # True obliquity

    # Convert input coordinates to equatorial (RA, Dec) if ecliptic
    if flag == ECL2HOR:
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
    if azimuth < 0:  # pragma: no cover - modulo by positive 360.0 is always >= 0
        azimuth += 360.0

    # Reference API convention for azalt(): pressure 0 does NOT disable
    # refraction here (that is refrac()'s convention) — it means "estimate
    # the pressure from the observer's altitude" using the standard
    # atmosphere, then apply refraction.
    if pressure <= 0:
        # Standard-atmosphere barometric formula. The base goes negative above
        # ~44.3 km, where a fractional power would yield a complex number; clamp
        # it to 0 so extreme altitudes give pressure ~0 (refraction ~0) instead
        # of crashing.
        base = 1.0 - 0.0065 * altitude / 288.0
        pressure = 1013.25 * (max(base, 0.0) ** 5.255)

    # Apparent altitude via refrac_extended: the reference's azalt apparent
    # altitude is identical to its refrac_extended output (including the
    # below-horizon dip clamp), so route through the same function and do not
    # re-clamp here.
    alt_apparent, _ = refrac_extended(
        alt_true,
        altitude,
        pressure,
        temperature,
        flag=TRUE_TO_APP,
    )

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

    Compatible with the reference azalt_rev() API.

    Note: This function is not precisely the reverse of azalt(). If only an
    apparent altitude is available, the true altitude must first be computed
    using a refraction correction.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        flag: Output coordinate type flag:
            - HOR2EQU (1): Output is equatorial (RA, Dec)
            - HOR2ECL (0): Output is ecliptic (longitude, latitude)
        geopos: Tuple of (longitude, latitude, altitude):
            - longitude: Geographic longitude of observer in degrees (East positive)
            - latitude: Geographic latitude of observer in degrees (North positive)
            - altitude: Observer altitude above sea level in meters
        azimuth: Azimuth in degrees from South, westward (0 = South, 90 = West, etc.)
        true_altitude: True (geometric) altitude above horizon in degrees

    Returns:
        Tuple of (x1, x2) where:
            - If flag == HOR2EQU: (Right Ascension, Declination) in degrees
            - If flag == HOR2ECL: (Ecliptic longitude, Ecliptic latitude) in degrees

    Example:
        >>> from libephemeris import azalt_rev, HOR2EQU, julday
        >>> jd = julday(2024, 6, 15, 12.0)
        >>> # Object at azimuth 90° (West), altitude 45°, observer at Rome
        >>> geopos = (12.5, 41.9, 0)  # (lon, lat, altitude)
        >>> ra, dec = azalt_rev(jd, HOR2EQU, geopos, 90.0, 45.0)
    """
    # Extract geopos components
    lon = geopos[0]
    lat = geopos[1]
    geopos[2]

    from .time_utils import _sidtime_internal
    from .state import get_timescale

    # Get true obliquity (mean obliquity + nutation in obliquity)
    # This is needed for accurate equatorial to ecliptic conversion
    ts = get_timescale()
    t = ts.ut1_jd(tjdut)
    jd_tt = t.tt

    # Mean obliquity via the Vondrák 2011 long-term model — the same of-date
    # mean obliquity calc_ut/fast_calc use to express ecliptic-of-date
    # coordinates, so this horizon conversion shares their ecliptic frame. The
    # IAU 2006 obliquity polynomial (erfa.obl06) diverges by several arcsec at
    # deep-BCE dates and would put the ecliptic here in a different frame.
    from .precession_vondrak import vondrak_mean_obliquity_deg

    eps0 = vondrak_mean_obliquity_deg(jd_tt)

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

    if abs(cos_dec) < 1e-10 or abs(cos_lat) < 1e-10:
        # Object at a celestial pole, or observer at a geographic pole:
        # the hour angle is genuinely undefined there.
        ha = 0.0
    else:
        sin_H = math.sin(az_rad) * math.cos(alt_rad) / cos_dec
        # cos(H) from the altitude identity
        #   sin(alt) = sin(lat)·sin(dec) + cos(lat)·cos(dec)·cos(H)
        # which avoids the sin(lat) division of the azimuth identity —
        # geographic latitude 0 (a perfectly normal input) is no longer a
        # degenerate case returning RA = LST for every azimuth/altitude.
        cos_H = (math.sin(alt_rad) - sin_lat * sin_dec) / (cos_dec * cos_lat)

        ha_rad = math.atan2(sin_H, cos_H)
        ha = math.degrees(ha_rad)

    # Calculate Local Sidereal Time
    # Use apparent sidereal time (with nutation) for accuracy
    lst_hours = _sidtime_internal(tjdut, lon, eps, dpsi_deg)
    lst_deg = lst_hours * 15.0  # Convert hours to degrees

    # Calculate Right Ascension: RA = LST - H
    ra = (lst_deg - ha) % 360.0

    if flag == HOR2EQU:
        # Return equatorial coordinates (RA, Dec)
        return (ra, dec)
    else:
        # HOR2ECL: Convert equatorial to ecliptic
        # cotrans with positive obliquity converts equatorial to ecliptic
        ecl_coord = cotrans((ra, dec, 1.0), eps)
        return (ecl_coord[0], ecl_coord[1])


def refrac(
    alt: float,
    atpress: float = 1013.25,
    attemp: float = 15.0,
    flag: int = TRUE_TO_APP,
) -> float:
    """Convert between true and apparent altitude.

    Args:
        alt: Input altitude in degrees.
        atpress: Atmospheric pressure in hPa. Zero disables refraction;
            negative values are accepted as compatibility extrapolations.
        attemp: Atmospheric temperature in degrees Celsius.
        flag: ``TRUE_TO_APP`` adds refraction; ``APP_TO_TRUE`` removes it.

    Returns:
        Apparent altitude for ``TRUE_TO_APP``, or true altitude for
        ``APP_TO_TRUE``, in degrees.

    Notes:
        Plain refraction follows a clean-room closed-form behavioral model of
        the compatibility API. The forward and reverse directions are not
        exact numerical inverses and preserve the API's horizon and zenith
        clamps. At true altitude zero under the defaults, apparent altitude is
        about 28.6 arcminutes higher.
    """
    from .refraction import (
        calc_refrac_compat_app_to_true,
        calc_refrac_compat_true_to_app,
    )

    if flag == TRUE_TO_APP:
        return calc_refrac_compat_true_to_app(alt, atpress, attemp)
    return calc_refrac_compat_app_to_true(alt, atpress, attemp)


def refrac_extended(
    alt: float,
    geoalt: float,
    atpress: float = 1013.25,
    attemp: float = 15.0,
    lapserate: Optional[float] = None,
    flag: int = TRUE_TO_APP,
) -> Tuple[float, Tuple[float, float, float, float]]:
    """
    Calculate true altitude from apparent altitude, or vice-versa (extended).

    This is an extended version of refrac() that includes:
    - Observer altitude above sea level
    - Atmospheric lapse rate (temperature variation with altitude)
    - Dip of the horizon calculation

    Compatible with the reference refrac_extended() API.

    Args:
        alt: Altitude of object above geometric horizon in degrees.
                  For TRUE_TO_APP, this is the true (geometric) altitude.
                  For APP_TO_TRUE, this is the apparent (observed) altitude.
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
            - TRUE_TO_APP (0): Convert true altitude to apparent altitude
            - APP_TO_TRUE (1): Convert apparent altitude to true altitude

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
        calc_dip,
        calc_refraction_app_to_true,
        calc_refraction_true_to_app,
    )

    # Use global lapse rate if none provided
    if lapserate is None:
        lapserate = get_lapse_rate()

    # Dip of the horizon for elevated observers (pressure/temperature scale the
    # refraction part of the dip). Below the dipped horizon the input altitude
    # is returned untouched.
    dip = calc_dip(geoalt, lapserate, atpress, attemp)

    # Independently implemented ICAO Standard Atmosphere ray tracing.
    if flag == TRUE_TO_APP:
        true_alt = alt
        refraction = calc_refraction_true_to_app(
            alt, atpress, attemp, geoalt, lapserate
        )
        apparent_alt = true_alt + refraction
        # Below the (dipped) horizon: report the input unchanged.
        if apparent_alt < dip:
            return (true_alt, (true_alt, true_alt, 0.0, dip))
        return (apparent_alt, (true_alt, apparent_alt, refraction, dip))
    else:
        apparent_alt = alt
        # Apparent altitude below the (dipped) horizon: no refraction removed.
        if apparent_alt < dip:
            return (apparent_alt, (apparent_alt, apparent_alt, 0.0, dip))
        refraction = calc_refraction_app_to_true(
            alt, atpress, attemp, geoalt, lapserate
        )
        true_alt = apparent_alt - refraction
        return (true_alt, (true_alt, apparent_alt, refraction, dip))


def cotrans(
    coord: Tuple[float, float, float], eps: float
) -> Tuple[float, float, float]:
    """
    Transform coordinates between ecliptic and equatorial systems.

    Compatible with the reference cotrans() API.

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
    # Negate obliquity to match the reference API convention
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

    # Calculate the new longitude (RA for ecl→eq, λ for eq→ecl) from the full
    # 3D rotation about the x-axis. Using the true cartesian components
    # y' = cos_lat·sin_lon·cos_eps - sin_lat·sin_eps, x' = cos_lat·cos_lon
    # (rather than dividing them by cos_lat, i.e. tan_lat) keeps the result
    # correct when |lat| > 90°, where cos_lat < 0 would otherwise flip both
    # signs and rotate the longitude by 180°. Identical to the naive-trig form
    # for every |lat| <= 90°.
    y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps
    x = cos_lat * cos_lon

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

    Compatible with the reference degnorm() API.
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
    # Near-zero inputs normalize to exactly 0.0 rather than wrapping to
    # ~360: a value like -5e-14 is the angle zero with roundoff on it, and
    # the modulo alone would report it as 359.99999999999994 — the same
    # direction, but a full turn away for any consumer that compares or
    # bins the result. The width of the snap band is this library's own
    # choice: 1e-13 deg sits above the roundoff floor of the wrap near zero
    # (the ulp of 360.0 is ~5.7e-14) and far below any angular resolution
    # the library reports (1e-13 deg = 3.6e-10 arcsec). The boundary is
    # observable — |x| < 1e-13 snaps, |x| >= 1e-13 wraps normally, so an
    # input like -1e-12 comes back as 359.999999999999 — and once chosen it
    # is held fixed as documented behaviour of this function. The snap runs
    # before the wrap so a tiny negative never reaches the modulo.
    if abs(x) < 1e-13:
        return 0.0
    result = x % 360.0
    # Python's modulo can return exactly 360.0 for tiny negative inputs
    # (e.g. (-1e-17) % 360.0 == 360.0); snap that artifact to 0.0 to keep
    # the documented [0, 360) contract (legitimate near-360 values pass).
    if result >= 360.0:
        return 0.0
    return result


TWO_PI = 2.0 * math.pi


def radnorm(x: float) -> float:
    """
    Normalize an angle to the range [0, 2*pi).

    Compatible with the reference radnorm() API.
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
    # Same snap as degnorm, and the same figure for the band: 1e-13 rad.
    # The choice is made per unit, not by converting the degree band (that
    # would be 1.7e-15 rad, within a couple of ulps of the wrap's own
    # roundoff). The band only has to clear the roundoff floor of the
    # modulo near zero (the ulp of 2*pi is ~8.9e-16, two orders below) and
    # stay under any angular resolution the library reports (1e-13 rad ~
    # 2e-8 arcsec); one round figure that does both in either unit keeps
    # the two helpers symmetric. As in degnorm the boundary is observable
    # and is held fixed as documented behaviour of this function.
    if abs(x) < 1e-13:
        return 0.0
    result = x % TWO_PI
    # Same snap as degnorm: keep the [0, 2*pi) contract for tiny negative
    # inputs whose modulo lands exactly on the upper bound.
    if result >= TWO_PI:
        return 0.0
    return result


def difdeg2n(p1: float, p2: float) -> float:
    """
    Calculate distance in degrees p1 - p2 normalized to [-180;180].

    Compatible with the reference difdeg2n() API.
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

    Compatible with the reference difdegn() API.
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

    Compatible with the reference difcs2n() API.

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

    Compatible with the reference difcsn() API.

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

    Compatible with the reference csnorm() API.

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


#: Hundredths of an arcsecond in one whole arcsecond, and half of that: the
#: bias that turns a truncation into a rounding.
_CS_PER_ARCSEC: int = 100
_CS_HALF_ARCSEC: int = _CS_PER_ARCSEC // 2

#: Hundredths of an arcsecond in one 30-degree sign.
_CS_PER_SIGN: int = CS360 // 12


def csroundsec(cs: int) -> int:
    """
    Round an angle in centiseconds to the nearest whole arcsecond.

    The argument is an angle counted in hundredths of an arcsecond; the answer
    is the same angle to the nearest whole arcsecond, still counted in
    hundredths, so it is always a multiple of 100.

    Non-negative values round half up. Negative values take the same
    half-arcsecond bias and are then truncated toward zero, so the negative
    half of the domain rounds toward positive infinity: the function is
    neither odd nor idempotent there. One interval is restored explicitly, so
    that an angle at or below one whole arcsecond below zero does not collapse
    to zero.

    A rounded value is not allowed to arrive exactly on a multiple of 30
    degrees from strictly inside the sector below it: it is held one whole
    arcsecond short of the boundary, so that a reading stays inside the sign
    it belongs to. The hold is a convention of this formatting surface and
    carries no astronomical meaning.

    Compatible with the reference csroundsec() API.

    Args:
        cs: Angle in centiseconds (any integer, of either direction).

    Returns:
        The angle rounded to a whole arcsecond, in centiseconds; always a
        multiple of 100, and never more than 149 centiseconds from ``cs``.

    Notes:
        - 1 centisecond = 1/100 arcsecond = 1/360000 degree
        - 30° = 10,800,000 centiseconds; 360° = 129,600,000 centiseconds

    Examples:
        >>> csroundsec(149)
        100
        >>> csroundsec(150)
        200
        >>> csroundsec(-150)
        -100
        >>> csroundsec(-120)  # the restored interval
        -100
        >>> csroundsec(10799960)  # held short of the 30° boundary
        10799900
        >>> csroundsec(129600000)  # already on a boundary: unchanged
        129600000
    """
    biased = cs + _CS_HALF_ARCSEC
    if biased >= 0:
        rounded = biased // _CS_PER_ARCSEC * _CS_PER_ARCSEC
    else:
        # Truncation toward zero: this is what makes the negative half of the
        # domain asymmetric with the positive one.
        rounded = -(-biased // _CS_PER_ARCSEC) * _CS_PER_ARCSEC
    if rounded == 0:
        # The bias swallows the last arcsecond below zero; give it back.
        return -_CS_PER_ARCSEC if cs <= -_CS_PER_ARCSEC else 0
    if rounded % _CS_PER_SIGN == 0 and cs < rounded:
        # The value came from strictly inside the sector below the boundary,
        # so it stops one whole arcsecond short of it.
        return rounded - _CS_PER_ARCSEC
    return rounded


def cs2degstr(cs: int) -> str:
    """
    Convert a value in centiseconds to a formatted degrees string.

    This function converts an angular measurement in centiseconds (1/100 of an
    arcsecond) to a human-readable string in the format "D°M'S.ss\"" where D is
    degrees, M is arcminutes, S is arcseconds, and ss is centiseconds (hundredths
    of arcsecond).

    Compatible with the reference cs2degstr() API.

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

    Compatible with the reference cs2lonlatstr() API.

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

    # Format matching the reference API: "{deg}{char}{min:02d}" or
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

    Compatible with the reference cs2timestr() API.

    Args:
        cs: Time in centiseconds (any integer value)

    Returns:
        Formatted string representing the time in hours, minutes, seconds.
        Format: "HH:MM:SS" with zero-padded fields (e.g., "12:34:56")

    Notes:
        - 1 centisecond = 1/100 second
        - 1 second = 100 centiseconds
        - 1 minute = 6000 centiseconds
        - 1 hour = 360000 centiseconds
        - Hours are zero-padded and wrap modulo 24, so negative values map
          into the 0-23 range (e.g., -360000 cs -> "23:00:00")
        - Seconds are rounded to whole numbers (centiseconds are rounded)

    Examples:
        >>> cs2timestr(0)
        '00:00:00'
        >>> cs2timestr(360000)  # 1 hour
        '01:00:00'
        >>> cs2timestr(4526050)  # 12:34:21 (with rounding from .50)
        '12:34:21'
        >>> cs2timestr(-360000)  # -1 hour wraps mod 24
        '23:00:00'
    """
    # Black-box compatibility: the public API accepts a bytes separator.
    if isinstance(sep, bytes):
        sep = sep.decode("ascii")

    # Wrap the whole magnitude into one day [0, 24h) before extracting fields.
    # Applying the sign to the hours field alone left sub-hour negatives (e.g.
    # -6000 cs = -1 min) indistinguishable from their positive counterpart and
    # broke the documented wrap ("-360000 cs -> 23:00:00"). Python's modulo maps
    # negatives into the positive range, so -6000 -> 23:59:00 as intended.
    cs = cs % (24 * 360000)

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

    # Format the string matching reference API format
    # Format: "%02d<sep>%02d<sep>%02d"
    # A rounding carry can push 23:59:59.5x up to 24:00:00; wrap back to 00.
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
    # Choose the arc from the RAW arguments. Normalizing first can move the
    # difference a few ULP off exactly -180, so the equal-arcs tie handler
    # below never fires and the result flips by half a turn
    # (deg_midp(-1e-12, 179.999999999999) must be ~90, not ~270). Snapping
    # first is wrong for the mirror reason — it turns a slightly shorter
    # negative arc INTO an exact tie — so neither normalization belongs
    # before this point; only the final midpoint is normalized.
    diff = x2 - x1
    x1 = x1 % 360.0

    # Reduce the RAW difference into (-180, 180] with a modulo, not with a
    # single conditional shift: an argument more than one turn from the other
    # (deg_midp(720, 0)) needs more than one wrap. Reducing the difference
    # rather than the arguments keeps the exact-opposition tie detectable —
    # normalizing the arguments first can move it a few ULP off 180 and flip
    # the result by half a turn.
    #
    # A raw difference of exactly -180 maps to +180 here, which is also the
    # convention for the equal-arcs tie: with both arcs the same length, take
    # the positive (clockwise) half.
    diff = diff % 360.0
    if diff > 180.0:
        diff -= 360.0

    # Calculate midpoint along the chosen arc
    midp = x1 + diff / 2.0

    # Normalize through degnorm so the midpoint inherits the same near-zero
    # snap and half-open-range guard as every other public angle (a midpoint
    # a hair below zero, e.g. deg_midp(0.0, -1e-13), otherwise reported
    # ~359.99999999999994 instead of 0.0).
    return degnorm(midp)


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
    # Arc chosen from the raw arguments, as in deg_midp (see there).
    _raw_diff = y - x
    x = x % TWO_PI

    # Same reduction of the raw difference as deg_midp (see there).
    diff = _raw_diff % TWO_PI
    if diff > math.pi:
        diff -= TWO_PI

    # Calculate midpoint along the chosen arc
    midp = x + diff / 2.0

    # Normalize through radnorm so the midpoint inherits the near-zero snap
    # and the half-open range (a bare modulo returned exactly 2*pi, outside
    # this function's own documented [0, 2*pi)).
    return radnorm(midp)


def d2l(d: float) -> int:
    """
    Convert a double (float) to a long integer with rounding.

    This function rounds a floating-point number to the nearest integer using
    "round half away from zero" semantics (also known as commercial rounding).
    This is standard practice in astronomical computation software.

    Compatible with the reference d2l() API.

    Args:
        d: A floating-point number to convert.

    Returns:
        The nearest integer to d. For values exactly halfway between two
        integers, rounds away from zero (e.g., 0.5 -> 1, -0.5 -> -1).

    Notes:
        - This differs from Python's built-in round() function, which uses
          "round half to even" (banker's rounding) for Python 3.
        - Used internally for coordinate conversions and also exposed
          publicly for consistency with the reference API.

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
SPLIT_DEG_ZODIACAL: int = 8  # Return raw 30-degree zodiac segment
SPLIT_DEG_NAKSHATRA: int = 1024  # Return raw nakshatra segment number
SPLIT_DEG_KEEP_SIGN: int = 16  # Rounding stays inside the current sign/nakshatra
SPLIT_DEG_KEEP_DEG: int = 32  # Rounding keeps the whole-degree part unchanged


#: Every rounding bit, coarsest first, with half of the unit it asks for in
#: degrees. When several bits are set the coarsest one wins and the finer ones
#: have no effect, so the first match in this order is the answer.
_HALF_ROUNDING_UNIT: Tuple[Tuple[int, float], ...] = (
    (SPLIT_DEG_ROUND_DEG, 0.5),
    (SPLIT_DEG_ROUND_MIN, 0.5 / 60.0),
    (SPLIT_DEG_ROUND_SEC, 0.5 / 3600.0),
)

#: The three rounding bits together: their presence, not their effect, decides
#: what the sub-arcsecond field carries.
_ROUNDING_BITS: int = SPLIT_DEG_ROUND_SEC | SPLIT_DEG_ROUND_MIN | SPLIT_DEG_ROUND_DEG

#: The twelve equal arcs of the ecliptic: their width in degrees, and how many
#: of them make one turn.
_SIGN_SPAN_DEG: float = 30.0
_SIGNS_PER_TURN: int = 12


def _rounding_offset(roundflag: int) -> float:
    """Half of the rounding unit the flag word asks for, in degrees.

    Adding this to a magnitude and then truncating is the "round half up"
    rule: a value exactly half a unit above a whole multiple goes up.

    Args:
        roundflag: The SPLIT_DEG_* flag word.

    Returns:
        0.5 degree, half an arcminute or half an arcsecond, whichever the
        coarsest rounding bit in the word asks for, and 0.0 when the word
        carries no rounding bit at all.
    """
    for bit, half_unit in _HALF_ROUNDING_UNIT:
        if roundflag & bit:
            return half_unit
    return 0.0


def _rounding_is_held_back(
    position: float,
    offset: float,
    roundflag: int,
    degree_span: float,
    division_span: float,
) -> bool:
    """Say whether KEEP_DEG or KEEP_SIGN cancels the rounding for this value.

    Both flags act by suppressing the rounding entirely rather than by
    clamping the answer afterwards. When both are set only the whole-degree
    test governs.

    Args:
        position: The position about to be reported, in whatever unit the
            caller reduces in (degrees, or arcseconds for the lunar segments).
        offset: Half of the rounding unit, in that same unit.
        roundflag: The SPLIT_DEG_* flag word.
        degree_span: One whole degree, in that same unit.
        division_span: The width of the division the position must not leave,
            in that same unit.

    Returns:
        True when the rounding must be dropped for this value.
    """
    if roundflag & SPLIT_DEG_KEEP_DEG:
        # With both flags set, the whole-degree test governs alone.
        return int((position + offset) / degree_span) != int(position / degree_span)
    if roundflag & SPLIT_DEG_KEEP_SIGN:
        return int((position + offset) / division_span) != int(position / division_span)
    return False


def _decompose_to_dms(ddeg: float, has_rounding: bool) -> Tuple[int, int, int, float]:
    """Break a non-negative angle in degrees into deg, min, sec, secfr.

    The plain and zodiacal splits reduce in the degree domain: each part is
    read off the running remainder and taken out of it before the next one is
    read, all in binary64. An angle that is an exact arcsecond count in
    decimal can therefore still leave a sub-arcsecond residue, and how much
    depends on the magnitude of the angle, because taking out a large
    whole-degree part costs mantissa bits.

    Args:
        ddeg: Non-negative angle in degrees. Whole degrees are not reduced
            modulo a turn.
        has_rounding: True when the flag word carries a rounding bit, in
            which case everything below the requested unit is carry and the
            sub-arcsecond field repeats the arcseconds instead.

    Returns:
        Whole degrees (unbounded), arcminutes and arcseconds in 0..59, and
        the sub-arcsecond field: the true fraction in [0.0, 1.0), or the
        arcseconds field as a float when a rounding bit is present.
    """
    ideg = int(ddeg)
    remainder = max(ddeg - ideg, 0.0)
    imin = int(remainder * 60.0)
    remainder = max(remainder - imin / 60.0, 0.0)
    arcsec = remainder * 3600.0
    isec = int(arcsec)
    secfr = float(isec) if has_rounding else max(arcsec - isec, 0.0)
    return ideg, imin, isec, secfr


def _decompose_arcsec_to_dms(
    pos_arcsec: float, has_rounding: bool
) -> Tuple[int, int, int, float]:
    """Break a non-negative arc-second position into deg, min, sec, secfr.

    The nakshatra split reduces in the arc-second domain (where the 13°20'
    segment is exactly 48000") so that positions on an exact degree, minute
    or segment boundary decompose without the sub-arcsecond truncation the
    degree-domain path would show. ``has_rounding`` mirrors
    :func:`_decompose_to_dms`.
    """
    total = max(pos_arcsec, 0.0)
    ideg = int(total / 3600.0)
    total = max(total - ideg * 3600.0, 0.0)
    imin = int(total / 60.0)
    total = max(total - imin * 60.0, 0.0)
    isec = int(total)
    if has_rounding:
        secfr = float(isec)
    else:
        secfr = max(total - isec, 0.0)
    return ideg, imin, isec, secfr


_NAK_SPAN_ARCSEC: float = 48000.0  # 360° * 3600 / 27, exact in the arcsec domain


#: How many lunar segments make one turn.
_NAK_PER_TURN: int = 27

#: One degree counted in arcseconds, the unit the lunar segments are reduced
#: in.
_ARCSEC_PER_DEG: float = 3600.0


def _split_deg_nakshatra(
    ddeg: float, roundflag: int
) -> Tuple[int, int, int, float, int]:
    """Split a non-negative angle against the twenty-seven lunar segments.

    The ecliptic is divided into 27 exact equal segments of 360°/27 = 13°20',
    and the reduction runs in the arcsecond domain, where that span is the
    exact integer 48000. An angle that is a whole number of arcseconds
    therefore decomposes cleanly here, where the degree-domain split of the
    same angle would show a sub-arcsecond residue.

    Args:
        ddeg: Non-negative angle in degrees, not reduced modulo a turn.
        roundflag: The SPLIT_DEG_* flag word.

    Returns:
        The position inside the segment (degrees 0..13, arcminutes,
        arcseconds, sub-arcsecond field) and the index of the segment. The
        index is not reduced modulo a turn, except that an index of exactly
        27 — a value carried through one whole turn — is reported as 0.
    """
    arcsec = ddeg * _ARCSEC_PER_DEG
    index = int(arcsec / _NAK_SPAN_ARCSEC)
    position = arcsec - index * _NAK_SPAN_ARCSEC
    offset = _rounding_offset(roundflag) * _ARCSEC_PER_DEG
    # The two KEEP tests read the position inside the segment, because a
    # segment boundary falls in the middle of a degree: a rounding that leaves
    # the whole-degree part alone can still carry into the next segment.
    if offset and _rounding_is_held_back(
        position, offset, roundflag, _ARCSEC_PER_DEG, _NAK_SPAN_ARCSEC
    ):
        offset = 0.0
    position += offset
    if position >= _NAK_SPAN_ARCSEC:
        position -= _NAK_SPAN_ARCSEC
        index += 1
    if index == _NAK_PER_TURN:
        index = 0
    ideg, imin, isec, secfr = _decompose_arcsec_to_dms(
        position, bool(roundflag & _ROUNDING_BITS)
    )
    return ideg, imin, isec, secfr, index


def split_deg(degree: float, roundflag: int = 0) -> Tuple[int, int, int, float, int]:
    """
    Split an angle in decimal degrees into sexagesimal parts.

    The angle is reported as whole degrees, arcminutes, arcseconds and a
    sub-arcsecond field, plus one further integer that says either which side
    of zero the angle was on or which division of the ecliptic it fell in.
    The flag word chooses how finely the value is rounded, against which
    division of the circle it is reported, and whether the rounding may carry
    across a whole degree or across the boundary of that division.

    The angle is not normalised: nothing is reduced modulo a turn, so 1080°
    reports 1080 whole degrees. Without a division flag the magnitude is split
    and the direction is reported separately, ``-0.0`` counting as
    non-negative.

    Rounding is half up on the magnitude, so a negative input rounds away from
    zero. When several rounding bits are set the coarsest wins. Everything
    below the requested unit is carry, not result: a caller asking for
    arcminutes reads the first two fields and ignores the rest.

    ZODIACAL reports the index of the 30° arc the angle falls in, counting
    from the vernal point, and the first three fields give the position
    inside that arc. NAKSHATRA reports the same ecliptic divided instead
    into 27 exact equal segments of 13°20'. Neither index is reduced modulo a
    turn, except that a raw index of exactly one whole turn — 12 signs, or 27
    segments — is reported as 0.

    A negative input with NAKSHATRA alone follows the ordinary signed split
    and reports -1; with ZODIACAL also set it takes the zodiacal split of its
    magnitude. With both bits set the non-negative NAKSHATRA split has
    precedence, while a negative input takes the zodiacal one.

    Compatible with the reference split_deg() API.

    Args:
        degree: Angle in decimal degrees, of either direction.
        roundflag: Bit mask of the SPLIT_DEG_* constants. Bits outside that
            set are ignored rather than refused.

    Returns:
        A tuple ``(degrees, minutes, seconds, fraction, sign)`` of four
        native ints and one native float. ``minutes`` and ``seconds`` are in
        0..59; ``fraction`` is the true fraction of an arcsecond in
        [0.0, 1.0), or the arcseconds field repeated as a float whenever the
        flag word carries a rounding bit; ``sign`` is +1 or -1 without a
        division flag, and the index of the division with one.

    Raises:
        ValueError: If ``degree`` is a NaN.
        OverflowError: If ``degree`` is an infinity.

    Examples:
        >>> split_deg(45.5)
        (45, 30, 0, 0.0, 1)
        >>> split_deg(-30.5)
        (30, 30, 0, 0.0, -1)
        >>> split_deg(45.5, SPLIT_DEG_ZODIACAL)  # 15°30' of the second sign
        (15, 30, 0, 0.0, 1)
        >>> split_deg(45.0, SPLIT_DEG_NAKSHATRA)  # 5° of the fourth segment
        (5, 0, 0, 0.0, 3)
        >>> split_deg(10.584, SPLIT_DEG_ROUND_SEC)  # 10°35'02"
        (10, 35, 2, 2.0, 1)
    """
    sign_out = 1 if degree >= 0.0 else -1
    magnitude = abs(degree)
    if roundflag & SPLIT_DEG_NAKSHATRA and sign_out > 0:
        return _split_deg_nakshatra(magnitude, roundflag)
    offset = _rounding_offset(roundflag)
    if offset and _rounding_is_held_back(
        magnitude, offset, roundflag, 1.0, _SIGN_SPAN_DEG
    ):
        offset = 0.0
    rounded = magnitude + offset
    has_rounding = bool(roundflag & _ROUNDING_BITS)
    if roundflag & SPLIT_DEG_ZODIACAL:
        index = int(rounded / _SIGN_SPAN_DEG)
        # The position is taken out with the raw index, before a whole turn is
        # folded back onto the first sign.
        position = rounded - index * _SIGN_SPAN_DEG
        ideg, imin, isec, secfr = _decompose_to_dms(position, has_rounding)
        if index == _SIGNS_PER_TURN:
            index = 0
        return ideg, imin, isec, secfr, index
    ideg, imin, isec, secfr = _decompose_to_dms(rounded, has_rounding)
    return ideg, imin, isec, secfr, sign_out


def calc_angles(jd_ut: float, lat: float, lon: float):
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
    from .planets import calc_ut
    from .constants import SUN, MOON, MERCURY, VENUS

    # Set observer location
    set_topo(lon, lat, 0)

    # Calculate angles
    angles_dict = calc_angles(jd_ut, lat, lon)

    # Calculate and add planet positions for Arabic parts
    sun_pos, _ = calc_ut(jd_ut, SUN, 0)
    moon_pos, _ = calc_ut(jd_ut, MOON, 0)
    mercury_pos, _ = calc_ut(jd_ut, MERCURY, 0)
    venus_pos, _ = calc_ut(jd_ut, VENUS, 0)

    angles_dict["Sun"] = sun_pos[0]
    angles_dict["Moon"] = moon_pos[0]
    angles_dict["Mercury"] = mercury_pos[0]
    angles_dict["Venus"] = venus_pos[0]

    # Cache for Arabic parts
    set_angles_cache(angles_dict)

    return angles_dict


def angular_separation(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
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
