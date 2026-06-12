"""
Astrological house system calculations for libephemeris.

Implements 25 distinct house systems across 26 character codes (A as alias of E):
- Placidus (P): Most common, time-based, fails at polar latitudes
- Koch (K): Birthplace system, similar to Placidus
- Porphyry (O): Space-based trisection
- Regiomontanus (R): Medieval rational system
- Campanus (C): Prime vertical system
- Equal (A/E): Equal 30° divisions from Ascendant
- Whole Sign (W): Whole zodiac signs from Ascendant sign
- Meridian (X): Equatorial meridian divisions
- Azimuthal/Horizontal (H): Based on horizon
- Polich-Page/Topocentric (T)
- Alcabitius (B): Ancient Arabic system
- Morinus (M): Equatorial divisions
- Krusinski-Pisa-Goelzer (U)
- Gauquelin (G): 36 sectors
- Vehlow (V): Equal from midpoint
- APC (Y): Astronomical Planetary Cusps
- Carter Poli-Equatorial (F)
- Pullen SD (L): Sinusoidal Delta / Neo-Porphyry
- Pullen SR (Q): Sinusoidal Ratio
- Sripati (S): Divide quadrants equally
- Natural Gradient (N): Equal from 1=Aries
- Equal from MC (D)
- Sunshine/Treindl (I)
- Sunshine/Makransky (i)
- Savard-A (J)

Main Functions:
- houses(): Calculate house cusps and angles (ASCMC)
- houses_ex(): Extended version with sidereal support
- house_pos(): Find which house a point is in
- house_name(): Get house system name

Precision Notes:
Core Calculations (Phase 1 improved):
- GAST (sidereal time): ~0.001 sec = ~0.015 arcsec (Skyfield IAU SOFA)
- Obliquity: ~0.01 arcsec (Laskar 1986 + IAU 2000B nutation)
- Vertex: Exact at equator, ~0.001° elsewhere (rigorous limiting formula)
- Iterative systems (Placidus/Koch): ~0.00036 arcsec convergence (1e-7°)
- Non-iterative systems: Limited by obliquity/GAST precision (~0.01 arcsec)

Expected Accuracy by System:
- Simple (Equal, Whole Sign, Vehlow): Exact (no astronomical calculations)
- Geometric (Porphyry, Morinus, Meridian): ~0.01-0.1 arcsec
- Iterative (Placidus, Koch): ~0.001-0.01° (a few arcseconds)
- Complex (Campanus, Regiomontanus, Topocentric): ~0.01°
- Horizontal: ~0.01° (with convergence fallback to Porphyry)

Comparison with pyswisseph:
- Typical agreement: 0.001-0.1° depending on system and location
- Test suite validates against 130+ reference cases with tolerances 0.1-1.0°

Polar Latitudes:
- Placidus, Koch undefined > ~66° latitude (circumpolar ecliptic points)
- Automatic fallback to Porphyry when iteration fails
- Equal/Whole Sign work at all latitudes

Algorithm Sources:
- Placidus: Time divisions of diurnal/nocturnal arcs
- Regiomontanus: Equator trisection projected to ecliptic
- Campanus: Prime vertical trisection
- Equal: Simple 30° additions
- Algorithms from Meeus "Astronomical Algorithms"

References:
- Meeus "Astronomical Algorithms" 2nd Ed., Ch. 13 (coordinate systems)
- Reference documentation (house systems)
- Hand "Astrological Houses" (comprehensive house treatise)
- IERS Conventions 2003 (nutation models)
"""

from __future__ import annotations

import math
from typing import Any, List, Optional, Tuple, Union, overload
from .constants import *
from .constants import (
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_EQUATORIAL,
    SUN,
    FLG_JPLEPH,
    FLG_SWIEPH,
    FLG_TOPOCTR,
)
from .state import get_timescale
from .planets import get_ayanamsa_ut, calc_ut
from .cache import get_true_obliquity
from .exceptions import PolarCircleError, validate_coordinates
from .utils import difdeg2n


def _init_cardinal_cusps(asc: float, mc: float) -> list:
    """Create 13-element cusp list with cardinal points and their opposites."""
    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[7] = (asc + 180) % 360.0
    cusps[4] = (mc + 180) % 360.0
    return cusps


def _set_opposite_cusps(cusps: list) -> None:
    """Set cusps 5,6,8,9 as 180-degree opposites of 11,12,2,3."""
    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0


def _is_polar_circle(lat: float, eps: float) -> bool:
    """
    Check if latitude is within the polar circle for house calculations.

    At polar latitudes (approximately >66.5°), some house systems like Placidus
    and Koch cannot be calculated because the ecliptic does not properly
    intersect the horizon. This occurs when abs(lat) + eps > 90°.

    Args:
        lat: Geographic latitude in degrees
        eps: True obliquity of the ecliptic in degrees

    Returns:
        True if within polar circle (house calculation will fail)
    """
    # The polar circle is where the sun can be circumpolar
    # This happens when |lat| + obliquity > 90°
    # For typical obliquity of ~23.44°, this is lat > ~66.56°
    return abs(lat) + eps > 90.0


def get_polar_latitude_threshold(obliquity: float = 23.44) -> float:
    """
    Calculate the polar latitude threshold for a given obliquity.

    At latitudes beyond this threshold, some house systems (Placidus, Koch,
    Gauquelin) cannot be calculated because ecliptic points can be circumpolar
    (never rising or setting).

    The threshold is calculated as: 90° - obliquity

    For the current epoch (J2000), obliquity is approximately 23.44°,
    giving a threshold of approximately 66.56°.

    Args:
        obliquity: True obliquity of the ecliptic in degrees.
                   Default is 23.44° (approximate J2000 value).

    Returns:
        The polar latitude threshold in degrees. Latitudes with
        abs(lat) > threshold will trigger polar circle errors.

    Example:
        >>> get_polar_latitude_threshold()
        66.56
        >>> get_polar_latitude_threshold(23.5)
        66.5
    """
    return 90.0 - obliquity


def _get_polar_circle_info(lat: float, eps: float, house_system: str) -> dict:
    """
    Get detailed information about a polar circle condition.

    Args:
        lat: Geographic latitude in degrees
        eps: True obliquity of the ecliptic in degrees
        house_system: House system character (e.g., 'P', 'K', 'G')

    Returns:
        Dictionary with polar circle information:
        - is_polar: Whether the location is within the polar circle
        - latitude: The input latitude
        - threshold: The polar circle threshold for the given obliquity
        - obliquity: The obliquity used
        - excess: How many degrees beyond the threshold (if polar)
        - house_system: The house system character
        - hemisphere: 'N' for north, 'S' for south
    """
    threshold = get_polar_latitude_threshold(eps)
    is_polar = abs(lat) > threshold
    excess = abs(lat) - threshold if is_polar else 0.0
    hemisphere = "N" if lat >= 0 else "S"

    return {
        "is_polar": is_polar,
        "latitude": lat,
        "threshold": threshold,
        "obliquity": eps,
        "excess": excess,
        "house_system": house_system,
        "hemisphere": hemisphere,
    }


# Default threshold for extreme latitude (80°)
EXTREME_LATITUDE_THRESHOLD = 80.0


def _is_extreme_latitude(
    lat: float, threshold: float = EXTREME_LATITUDE_THRESHOLD
) -> bool:
    """
    Check if latitude is at an extreme value where house calculations may be numerically unstable.

    At latitudes beyond 80°, many house systems (especially Campanus, Regiomontanus,
    and Topocentric) may produce results with reduced accuracy or numerical instability,
    even if they don't technically fail like Placidus/Koch at polar latitudes.

    This is distinct from _is_polar_circle() which checks the mathematical limit
    where Placidus/Koch become undefined (~66.5°).

    Args:
        lat: Geographic latitude in degrees
        threshold: Latitude threshold for "extreme" (default 80°)

    Returns:
        True if latitude is at an extreme value (abs(lat) >= threshold)

    Example:
        >>> _is_extreme_latitude(85.0)
        True
        >>> _is_extreme_latitude(70.0)
        False
    """
    return abs(lat) >= threshold


def get_extreme_latitude_info(lat: float, obliquity: float = 23.44) -> dict:
    """
    Get detailed information about a location's latitude classification for house calculations.

    This function provides comprehensive information about how latitude affects
    house calculations, including whether the location is at extreme latitude,
    within the polar circle, and which house systems are affected.

    Args:
        lat: Geographic latitude in degrees
        obliquity: True obliquity of the ecliptic in degrees.
                   Default is 23.44° (approximate J2000 value).

    Returns:
        Dictionary with latitude information:
        - latitude: The input latitude
        - is_extreme: True if abs(lat) >= 80° (may have numerical instability)
        - is_polar_circle: True if abs(lat) > 90° - obliquity (Placidus/Koch/Gauquelin fail)
        - polar_threshold: The polar circle threshold latitude
        - extreme_threshold: The extreme latitude threshold (80°)
        - hemisphere: 'N' for north, 'S' for south
        - affected_systems: List of house system codes that fail at this latitude
        - unstable_systems: List of house system codes that may be numerically unstable
        - stable_systems: List of house system codes that work reliably at all latitudes

    Example:
        >>> info = get_extreme_latitude_info(85.0)
        >>> info['is_extreme']
        True
        >>> info['is_polar_circle']
        True
        >>> info['affected_systems']
        ['P', 'K', 'G']
        >>> info['stable_systems']
        ['E', 'W', 'O', 'M', 'X', 'V', 'N']
    """
    polar_threshold = get_polar_latitude_threshold(obliquity)
    is_polar = abs(lat) > polar_threshold
    is_extreme = _is_extreme_latitude(lat)
    hemisphere = "N" if lat >= 0 else "S"

    # Systems that completely fail at polar latitudes
    affected_systems = ["P", "K", "G"] if is_polar else []

    # Systems that may be numerically unstable at extreme latitudes (>80°)
    # but don't technically fail
    unstable_systems = []
    if is_extreme:
        unstable_systems = ["C", "R", "T", "B", "H", "U", "Y", "F"]

    # Systems that work reliably at all latitudes
    stable_systems = ["E", "W", "O", "M", "X", "V", "N"]

    return {
        "latitude": lat,
        "is_extreme": is_extreme,
        "is_polar_circle": is_polar,
        "polar_threshold": polar_threshold,
        "extreme_threshold": EXTREME_LATITUDE_THRESHOLD,
        "hemisphere": hemisphere,
        "affected_systems": affected_systems,
        "unstable_systems": unstable_systems,
        "stable_systems": stable_systems,
    }


def _validate_cusps(cusps: list | tuple) -> tuple[bool, str | None]:
    """
    Validate house cusps for numerical sanity.

    Checks that all cusp values are within valid range [0, 360) and are finite.
    This helps detect numerical instability at extreme latitudes.

    Args:
        cusps: List or tuple of house cusp longitudes

    Returns:
        Tuple of (is_valid, error_message or None)

    Example:
        >>> _validate_cusps([0.0, 30.0, 60.0, 90.0])
        (True, None)
        >>> _validate_cusps([float('nan'), 30.0])
        (False, 'House cusp contains NaN value')
    """
    for i, cusp in enumerate(cusps):
        # Check for NaN
        if cusp != cusp:  # NaN check (NaN != NaN is True)
            return False, f"House cusp {i + 1} contains NaN value"

        # Check for infinity
        if cusp == float("inf") or cusp == float("-inf"):
            return False, f"House cusp {i + 1} contains infinite value"

        # Check range (should be in [0, 360))
        if cusp < 0 or cusp >= 360:
            return False, f"House cusp {i + 1} out of range [0, 360): {cusp}"

    return True, None


def _raise_polar_circle_error(
    lat: float, eps: float, house_system: str, func_name: str
) -> None:
    """
    Raise a PolarCircleError with detailed information.

    Args:
        lat: Geographic latitude in degrees
        eps: True obliquity of the ecliptic in degrees
        house_system: House system character (e.g., 'P', 'K', 'G')
        func_name: Name of the function raising the error

    Raises:
        PolarCircleError: Always raised with detailed polar circle information
    """
    info = _get_polar_circle_info(lat, eps, house_system)
    threshold = info["threshold"]
    hemisphere = "Northern" if info["hemisphere"] == "N" else "Southern"

    # Map house system codes to names
    system_names = {
        "P": "Placidus",
        "K": "Koch",
        "G": "Gauquelin",
    }
    system_name = system_names.get(house_system, house_system)

    message = (
        f"{func_name}: {system_name} house system cannot be calculated at "
        f"latitude {abs(lat):.2f}°{info['hemisphere']} (within {hemisphere} polar circle). "
        f"Polar threshold for obliquity {eps:.2f}° is ±{threshold:.2f}°. "
        f"Consider using Porphyry ('O'), Equal ('E'), or Whole Sign ('W') house systems "
        f"which work at all latitudes, or use houses_with_fallback() for automatic fallback."
    )

    raise PolarCircleError(
        message=message,
        latitude=lat,
        threshold=threshold,
        obliquity=eps,
        house_system=house_system,
    )


def angular_diff(a: float, b: float) -> float:
    """
    Calculate signed angular difference (a - b) handling 360° wrapping.

    Used by horizontal house system to find ecliptic longitude for given azimuth.

    Args:
        a: First angle in degrees
        b: Second angle in degrees

    Returns:
        Signed difference in range [-180, 180]
    """
    diff = (a - b) % 360.0
    if diff > 180.0:
        diff -= 360.0
    return diff


def _calc_vertex(armc_deg: float, eps: float, lat: float, mc: float) -> float:
    """
    Calculate the Vertex (intersection of Prime Vertical and Ecliptic in Western hemisphere).

    The Vertex is an auxiliary angle representing where the Prime Vertical (great circle
    through zenith perpendicular to meridian) intersects the ecliptic in the western sky.
    Often used in astrology for fateful encounters or significant relationships.

    At the equator (lat=0), the formula has a 1/tan(lat) singularity. We clamp
    latitude to a tiny positive value so the formula evaluates to the correct
    limiting value, matching Swiss Ephemeris behavior.

    Args:
        armc_deg: Right Ascension of Midheaven (sidereal time) in degrees
        eps: True obliquity of ecliptic in degrees
        lat: Geographic latitude in degrees
        mc: Midheaven longitude in degrees (for hemisphere verification)

    Returns:
        Vertex longitude in degrees (western hemisphere)

    Precision: ~0.001° for non-equatorial latitudes
    """
    eps_rad = math.radians(eps)

    # At equator (lat=0), the Vertex formula has a 1/tan(lat) singularity.
    # Clamp to a tiny positive latitude so the formula evaluates to the
    # correct limiting value (matches Swiss Ephemeris behavior).
    if abs(lat) < 1e-10:
        lat = 1e-10

    # Standard formula: Vertex is where Prime Vertical intersects ecliptic in West
    armc_rad = math.radians(armc_deg)
    lat_rad = math.radians(lat)

    num = -math.cos(armc_rad)
    den = math.sin(armc_rad) * math.cos(eps_rad) - math.sin(eps_rad) / math.tan(lat_rad)

    vtx_rad = math.atan2(num, den)
    vtx = math.degrees(vtx_rad) % 360.0

    # Ensure Vertex is in Western Hemisphere relative to MC
    # Vertex should be West of MC (i.e., behind it in diurnal motion)
    diff = (vtx - mc) % 360.0
    if diff < 180.0:
        vtx = (vtx + 180.0) % 360.0

    return vtx


def _ra_to_ecliptic_longitude(
    ra_deg: float, pole_height_deg: float, sin_obliquity: float, cos_obliquity: float
) -> float:
    """
    Convert right ascension to ecliptic longitude via spherical trigonometry.

    Given a right ascension and pole height, computes the corresponding ecliptic
    longitude using the standard spherical trigonometric relation:

        tan(λ) = sin(α) / (cos(ε)·cos(α) - sin(ε)·tan(φ))

    where α is the right ascension, ε is the obliquity, and φ is the pole height.
    Uses atan2 for unambiguous quadrant determination.

    Reference: Smart, "Textbook on Spherical Astronomy", Ch. 3;
               Meeus, "Astronomical Algorithms", Ch. 13

    Args:
        ra_deg: Right ascension in degrees
        pole_height_deg: Pole height (geographic latitude for ascendant) in degrees
        sin_obliquity: Pre-computed sin(obliquity)
        cos_obliquity: Pre-computed cos(obliquity)

    Returns:
        Ecliptic longitude in degrees [0, 360)
    """
    _NEAR_ZERO = 1e-10

    # Polar degenerate cases
    if abs(90.0 - pole_height_deg) < _NEAR_ZERO:
        return 180.0
    if abs(90.0 + pole_height_deg) < _NEAR_ZERO:
        return 0.0

    ra_rad = math.radians(ra_deg % 360.0)
    tan_pole = math.tan(math.radians(pole_height_deg))

    sin_ra = math.sin(ra_rad)
    cos_ra = math.cos(ra_rad)

    # Standard spherical trigonometry formula
    numerator = sin_ra
    denominator = cos_obliquity * cos_ra - sin_obliquity * tan_pole

    # atan2 handles all quadrants correctly
    longitude = math.degrees(math.atan2(numerator, denominator)) % 360.0

    # Snap cardinal points to exact values (avoid floating-point drift)
    for cardinal in (90.0, 180.0, 270.0):
        if abs(longitude - cardinal) < _NEAR_ZERO:
            return cardinal
    if abs(longitude - 360.0) < _NEAR_ZERO or abs(longitude) < _NEAR_ZERO:
        return 0.0

    return longitude


def _calc_ascendant(
    armc_deg: float, eps: float, lat: float, pole_height: float
) -> float:
    """
    Calculate ecliptic longitude from ARMC and pole height.

    Wrapper around _ra_to_ecliptic_longitude that accepts obliquity in degrees
    and computes the trigonometric values internally.

    Based on standard spherical trigonometry:
        tan(λ) = sin(α) / (cos(ε)·cos(α) - sin(ε)·tan(φ))

    Reference: Smart, "Textbook on Spherical Astronomy";
               Meeus, "Astronomical Algorithms"

    Args:
        armc_deg: Right Ascension of MC in degrees
        eps: True obliquity of ecliptic in degrees
        lat: Geographic latitude (unused, kept for call-site compatibility)
        pole_height: Pole height (latitude parameter) in degrees

    Returns:
        Ecliptic longitude in degrees (0-360)
    """
    sin_obliquity = math.sin(math.radians(eps))
    cos_obliquity = math.cos(math.radians(eps))
    return _ra_to_ecliptic_longitude(
        armc_deg, pole_height, sin_obliquity, cos_obliquity
    )


def houses(
    tjdut: float, lat: float, lon: float, hsys: int = ord("P"), iflag: int = 0
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Calculate astrological house cusps and angles for a given time and location.

    Reference API compatible function. Computes house divisions according to
    the specified house system and returns both house cusps and major angles (ASCMC).

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        lat: Geographic latitude in degrees (positive North, negative South)
        lon: Geographic longitude in degrees (positive East, negative West)
        hsys: House system identifier (e.g., ord('P') for Placidus, ord('K') for Koch)
        iflag: Calculation flags (optional, default 0). Use FLG_MOSEPH for extended
               date range (-3000 to +3000 CE). The ephemeris flag is propagated to
               internal calculations (e.g., Sun position for Sunshine houses).

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes (houses 1-12) in degrees
            - ascmc: Tuple of 8 angles: [Asc, MC, ARMC, Vertex, EquAsc, CoAsc, CoAscKoch, PolarAsc]

    House Systems Supported:
        'P' = Placidus, 'K' = Koch, 'R' = Regiomontanus, 'C' = Campanus,
        'E'/'A' = Equal, 'W' = Whole Sign, 'O' = Porphyry, 'B' = Alcabitius,
        'T' = Topocentric, 'M' = Morinus, 'X' = Meridian, 'H' = Horizontal,
        'V' = Vehlow, 'G' = Gauquelin, 'U' = Krusinski, 'F' = Carter,
        'Y' = APC, 'N' = Natural Gradient

    Example:
        >>> cusps, ascmc = houses(2451545.0, 51.5, -0.12, ord('P'))  # London, Placidus
        >>> asc, mc = ascmc[0], ascmc[1]
        >>> house_1_start = cusps[0]  # First house cusp
    """
    # Validate latitude and longitude ranges
    validate_coordinates(lat, lon, "houses")

    # 1. Calculate Sidereal Time (ARMC)
    # ARMC = GMST + lon
    ts = get_timescale()
    t = ts.ut1_jd(tjdut)

    # Use Skyfield's GAST (Greenwich Apparent Sidereal Time) for house calculations
    # GAST includes nutation in right ascension, which is critical for accurate house cusps.
    # Unlike GMST (mean sidereal time), GAST accounts for the true position of the
    # equinox affected by nutation, providing ~0.015 arcsec precision in final cusps.
    # Skyfield GAST precision: ~0.001 seconds of time = ~0.015 arcsec in RA
    # Reference: Skyfield documentation, IAU SOFA standards (iau2000b nutation model)
    gast = float(t.gast)  # in hours (convert numpy.float64 to Python float)
    armc_deg = (gast * 15.0 + lon) % 360.0

    # True Obliquity of Ecliptic - uses cached nutation calculation
    # This is a hot path optimization: obliquity calculation with nutation
    # is expensive (~0.03ms) and called frequently. Caching provides ~50x speedup.
    eps = get_true_obliquity(t.tt)

    # 2. Calculate Ascendant and MC
    # MC is intersection of Meridian and Ecliptic.
    # tan(MC) = tan(ARMC) / cos(eps)
    # Quadrant check needed.

    # Determine if we need to flip MC (and thus ARMC) for specific systems
    # Regiomontanus (R), Campanus (C), Polich/Page (T) flip MC if below horizon.

    hsys_char: str
    if isinstance(hsys, int):
        hsys_char = chr(hsys)
    elif isinstance(hsys, bytes):
        hsys_char = hsys.decode("utf-8")
    else:
        hsys_char = str(hsys)

    # Determine if we need to flip MC (and thus ARMC) for specific systems
    # Regiomontanus (R), Campanus (C), Polich/Page (T) flip MC if below horizon.

    armc_active = armc_deg

    if hsys_char in ["R", "C", "T", "J"]:
        # Check altitude of MC calculated from original ARMC
        mc_dec_rad = math.atan(
            math.sin(math.radians(armc_deg)) * math.tan(math.radians(eps))
        )
        lat_rad = math.radians(lat)
        sin_alt = math.sin(lat_rad) * math.sin(mc_dec_rad) + math.cos(
            lat_rad
        ) * math.cos(mc_dec_rad)

        if sin_alt < 0:
            armc_active = (armc_deg + 180.0) % 360.0

    # 2. Calculate Ascendant and MC

    # MC uses armc_active (flipped if needed)
    mc_rad = math.atan2(
        math.tan(math.radians(armc_active)), math.cos(math.radians(eps))
    )
    mc = math.degrees(mc_rad)
    # Adjust quadrant to match armc_active
    if mc < 0:
        mc += 360.0

    if 90.0 < armc_active <= 270.0:
        if mc < 90.0 or mc > 270.0:
            mc += 180.0
    elif armc_active > 270.0:
        if mc < 270.0:
            mc += 180.0
    elif armc_active <= 90.0:
        if mc > 90.0:
            mc += 180.0

    mc = mc % 360.0

    # Ascendant uses armc_deg (Original)
    num = math.cos(math.radians(armc_deg))
    den = -(
        math.sin(math.radians(armc_deg)) * math.cos(math.radians(eps))
        + math.tan(math.radians(lat)) * math.sin(math.radians(eps))
    )
    asc_rad = math.atan2(num, den)
    asc = math.degrees(asc_rad)
    asc = asc % 360.0

    # Ensure Ascendant is on the Eastern Horizon (Azimuth in [0, 180])
    # We check Azimuth relative to the TRUE ARMC (armc_deg)

    asc_r = math.radians(asc)
    eps_r = math.radians(eps)

    # RA
    y = math.cos(eps_r) * math.sin(asc_r)
    x = math.cos(asc_r)
    ra_r = math.atan2(y, x)
    ra = math.degrees(ra_r) % 360.0

    # Dec
    dec_r = math.asin(math.sin(eps_r) * math.sin(asc_r))

    # Hour Angle using TRUE ARMC
    h_deg = (armc_deg - ra + 360.0) % 360.0
    h_r = math.radians(h_deg)

    # Azimuth
    # tan(Az) = sin(H) / (sin(lat)cos(H) - cos(lat)tan(Dec))
    lat_r = math.radians(lat)

    num_az = math.sin(h_r)
    den_az = math.sin(lat_r) * math.cos(h_r) - math.cos(lat_r) * math.tan(dec_r)
    az_r = math.atan2(num_az, den_az)
    az = math.degrees(az_r)
    az = (az + 180.0) % 360.0

    # Check if H is West (0-180). If so, Asc is Setting (Descendant).
    # We want Rising.
    if 0.0 < h_deg < 180.0:
        asc = (asc + 180.0) % 360.0

    # Ensure ASC is in [0, 360) range (handle floating-point near-360 values)
    if asc >= 360.0 or (360.0 - asc) < 1e-10:
        asc = 0.0

    # Vertex uses armc_deg (Original)
    # Hemisphere check relative to TRUE ARMC (West of True ARMC)
    # The horizontal system computes its equator-degenerate angles on
    # the negative-latitude side (reference behavior — its vertex and
    # coasc2 at lat 0 equal the lat -> 0- limits, other systems the
    # lat -> 0+ limits).
    _vtx_lat = -1e-9 if (hsys_char == "H" and abs(lat) < 1e-10) else lat
    vertex = _calc_vertex(armc_deg, eps, _vtx_lat, armc_deg)

    # Equatorial Ascendant (East Point)
    # This is the intersection of the ecliptic with the celestial equator in the east
    # It's the ecliptic longitude where RA = ARMC + 90°
    equ_asc_ra = (armc_deg + 90.0) % 360.0
    # Convert RA to ecliptic longitude
    # tan(Lon) = tan(RA) / cos(eps)
    # y = sin(RA)
    # x = cos(RA) * cos(eps)
    equ_asc_ra_r = math.radians(equ_asc_ra)
    eps_r = math.radians(eps)
    y = math.sin(equ_asc_ra_r)
    x = math.cos(equ_asc_ra_r) * math.cos(eps_r)
    equ_asc = math.degrees(math.atan2(y, x)) % 360.0

    # Co-Ascendant W. Koch (coasc1)
    # Co-Ascendant (Koch) formula:
    # coasc1 = Asc(ARMC - 90°, latitude) + 180°
    # This is the Ascendant calculated 90° westward on the equator, then opposite point
    coasc_armc = (armc_deg - 90.0) % 360.0
    co_asc_koch = _calc_ascendant(coasc_armc, eps, lat, lat)

    # Add 180° to get opposite point
    co_asc_koch = (co_asc_koch + 180.0) % 360.0

    # Co-Ascendant M. Munkasey (coasc2)
    # Co-Ascendant (Munkasey) formula:
    # If lat >= 0: coasc2 = Asc(ARMC + 90°, 90° - lat)
    # If lat < 0:  coasc2 = Asc(ARMC + 90°, -90° - lat)
    # At equator (lat=0), coasc2_lat becomes 90° which is undefined.
    # The one-sided limits are 180 (lat -> 0+) and 0 (lat -> 0-).
    # The reference outputs 180.0 at exactly lat 0 for every system
    # except the horizontal system 'H', whose internal latitude
    # handling lands on the negative-side limit 0.0 (verified per
    # system against pyswisseph at armc 0/90/180/270).
    if abs(lat) < 1e-10:
        co_asc = 0.0 if hsys_char == "H" else 180.0
    else:
        coasc2_armc = (armc_deg + 90.0) % 360.0
        if lat >= 0:
            coasc2_lat = 90.0 - lat
        else:
            coasc2_lat = -90.0 - lat
        co_asc = _calc_ascendant(coasc2_armc, eps, coasc2_lat, coasc2_lat)

    # Polar Ascendant M. Munkasey (polasc)
    # Polar Ascendant formula:
    # polasc = Asc(ARMC - 90°, latitude)
    # Note: This is the same as coasc1 but WITHOUT the +180°
    polar_asc = _calc_ascendant(coasc_armc, eps, lat, lat)

    # Build ASCMC array with 8 elements (reference API compatible)
    ascmc = [0.0] * 8
    ascmc[0] = asc
    ascmc[1] = mc
    ascmc[2] = armc_deg
    ascmc[3] = vertex
    ascmc[4] = equ_asc
    ascmc[5] = co_asc_koch  # coasc1 (W. Koch) at index 5
    ascmc[6] = co_asc  # coasc2 (M. Munkasey) at index 6
    ascmc[7] = polar_asc

    # 3. House Cusps
    # Use armc_active for house calculations
    # If MC was flipped, we might need to flip latitude for intermediate cusps (Regiomontanus, etc.)
    # Verified for Regiomontanus: using -lat with flipped MC matches reference.

    calc_lat = lat
    if armc_active != armc_deg:
        # MC was flipped. Flip latitude for intermediate cusp calculations.
        calc_lat = -lat

    # Check for polar circle condition for Placidus/Koch/Gauquelin
    # These systems cannot be calculated when abs(lat) + eps > 90°
    # Raise detailed PolarCircleError with useful information
    if hsys_char in ["P", "K", "G"] and _is_polar_circle(lat, eps):
        _raise_polar_circle_error(lat, eps, hsys_char, "houses")

    # Calculate Sun's declination for Sunshine houses ('I' or 'i')
    # This is needed before the house dispatch since only houses has jd_ut
    sun_dec = 0.0
    if hsys_char in ("I", "i"):
        try:
            # Extract ephemeris flags: FLG_JPLEPH=1, FLG_SWIEPH=2
            eph_flags = iflag & (FLG_JPLEPH | FLG_SWIEPH)
            sun_pos, _ = calc_ut(tjdut, SUN, FLG_EQUATORIAL | eph_flags)
            sun_dec = sun_pos[1]  # Declination is second element in equatorial coords
        except (IndexError, TypeError, ValueError):
            # Fallback to 0 declination (same as equinox behavior)
            sun_dec = 0.0

    cusps = [0.0] * 13

    if hsys_char == "P":  # Placidus
        cusps = _houses_placidus(
            armc_active, lat, eps, asc, mc
        )  # Placidus fails anyway
    elif hsys_char == "K":  # Koch
        cusps = _houses_koch(armc_active, lat, eps, asc, mc)  # Koch fails anyway
    elif hsys_char == "R":  # Regiomontanus
        cusps = _houses_regiomontanus(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "C":  # Campanus
        cusps = _houses_campanus(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "E":  # Equal (Ascendant)
        cusps = _houses_equal(asc)
    elif hsys_char == "A":  # Equal (MC)
        cusps = _houses_equal(asc)  # Equal MC uses the same algorithm as Equal Asc
    elif hsys_char == "W":  # Whole Sign
        cusps = _houses_whole_sign(asc)
    elif hsys_char == "O":  # Porphyry
        cusps = _houses_porphyry(asc, mc)
    elif hsys_char == "B":  # Alcabitius
        cusps = _houses_alcabitius(armc_active, lat, eps, asc, mc)
    elif hsys_char == "T":  # Polich/Page (Topocentric)
        cusps = _houses_polich_page(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "M":  # Morinus
        cusps = _houses_morinus(armc_active, lat, eps, asc, mc)
    elif hsys_char == "X":  # Meridian (Axial)
        cusps = _houses_meridian(armc_active, lat, eps, asc, mc)
    elif hsys_char == "V":  # Vehlow
        cusps = _houses_vehlow(asc)
    elif hsys_char == "H":  # Horizontal (Azimuthal)
        cusps = _houses_horizontal(armc_active, lat, eps, asc, mc)
    elif hsys_char == "Y":  # APC Houses
        cusps = _houses_apc(armc_active, lat, eps, asc, mc)
        # APC at polar latitudes needs cusps and MC flipped if MC is below horizon
        # (reference behavior - different from R/C/T which flip armc_active)
        mc_dec_rad = math.atan(
            math.sin(math.radians(armc_deg)) * math.tan(math.radians(eps))
        )
        lat_rad = math.radians(lat)
        sin_alt = math.sin(lat_rad) * math.sin(mc_dec_rad) + math.cos(
            lat_rad
        ) * math.cos(mc_dec_rad)
        if sin_alt < 0:
            # Flip MC in ascmc
            ascmc[1] = (ascmc[1] + 180.0) % 360.0
            # Flip all cusps by 180°
            for i in range(1, 13):
                cusps[i] = (cusps[i] + 180.0) % 360.0
    elif hsys_char == "F":  # Carter (Poli-Equatorial)
        cusps = _houses_carter(armc_active, lat, eps, asc, mc)
    elif hsys_char == "U":  # Krusinski
        cusps = _houses_krusinski(armc_active, lat, eps, asc, mc)
    elif hsys_char == "N":  # Natural Gradient
        cusps = _houses_natural_gradient(armc_active, lat, eps, asc, mc)
    elif hsys_char == "G":  # Gauquelin
        cusps = _houses_gauquelin(armc_active, lat, eps, asc, mc)
    elif hsys_char == "S":  # Sripati
        cusps = _houses_sripati(asc, mc)
    elif hsys_char == "L":  # Pullen SD (Sinusoidal Delta / Neo-Porphyry)
        cusps = _houses_pullen_sd(asc, mc)
    elif hsys_char == "Q":  # Pullen SR (Sinusoidal Ratio)
        cusps = _houses_pullen_sr(asc, mc)
    elif hsys_char == "D":  # Equal from MC
        cusps = _houses_equal_mc(asc, mc)
    elif hsys_char == "I":  # Sunshine (Treindl)
        cusps = _houses_sunshine(armc_active, lat, eps, asc, mc, sun_dec)
        # Sunshine may flip MC internally (mc_under_horizon handling).
        # Sync ascmc[1] back from cusps[10] so the returned MC matches.
        ascmc[1] = cusps[10]
    elif hsys_char == "i":  # Sunshine (Makransky)
        cusps = _houses_sunshine_makransky(armc_active, lat, eps, asc, mc, sun_dec)
        ascmc[1] = cusps[10]
    elif hsys_char == "J":  # Savard-A
        cusps = _houses_savard_a(armc_active, calc_lat, eps, asc, mc)
    else:
        # Default to Placidus
        cusps = _houses_placidus(armc_active, lat, eps, asc, mc)

    # Return cusps array (reference API compatible: no padding at index 0)
    # For Gauquelin ('G'), return 36 sectors; otherwise return 12 houses
    if hsys_char == "G":
        return tuple(cusps[1:37]), tuple(ascmc)
    return tuple(cusps[1:13]), tuple(ascmc)


def houses_with_fallback(
    jd_ut: float,
    lat: float,
    lon: float,
    hsys: int,
    fallback_hsys: int = ord("O"),
    validate_cusps: bool = True,
) -> tuple[tuple[float, ...], tuple[float, ...], bool, str | None]:
    """
    Calculate house cusps with automatic fallback for polar latitudes.

    This convenience function attempts to calculate houses using the requested
    system, and automatically falls back to a polar-safe system (default: Porphyry)
    if the location is within the polar circle.

    For extreme latitudes (>80°), this function also:
    - Validates cusp values for numerical sanity (NaN, infinity, range)
    - Includes warnings about potential numerical instability
    - Falls back to a stable system if cusps are invalid

    This is useful for applications that need to handle polar latitudes gracefully
    without explicit error handling for each calculation.

    Args:
        jd_ut: Julian Day in Universal Time
        lat: Geographic latitude in degrees (positive North, negative South)
        lon: Geographic longitude in degrees (positive East, negative West)
        hsys: Primary house system identifier (e.g., ord('P') for Placidus)
        fallback_hsys: Fallback house system for polar latitudes.
                       Default: ord('O') (Porphyry).
                       Other good choices: ord('E') (Equal), ord('W') (Whole Sign)
        validate_cusps: If True (default), validate cusp values for numerical sanity
                        and fall back if invalid cusps are detected.

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes (houses 1-12) in degrees
            - ascmc: Tuple of 8 angles
            - used_fallback: True if fallback system was used
            - warning_message: Informational message if fallback was used or
                               extreme latitude warning, else None

    Example:
        >>> import libephemeris as ephem
        >>> jd = 2451545.0
        >>> # This would fail with Placidus at 70°N, but uses fallback
        >>> cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
        ...     jd, 70.0, 0.0, ord('P')
        ... )
        >>> if used_fallback:
        ...     print(f"Used fallback: {warning}")
        >>>
        >>> # At extreme latitudes (>80°), you may get a warning even for non-failing systems
        >>> cusps, ascmc, used_fallback, warning = ephem.houses_with_fallback(
        ...     jd, 85.0, 0.0, ord('C')  # Campanus at 85°N
        ... )
        >>> if warning:
        ...     print(f"Warning: {warning}")

    See Also:
        houses: Standard function that raises PolarCircleError at polar latitudes
        get_polar_latitude_threshold: Returns the threshold for a given obliquity
        get_extreme_latitude_info: Returns detailed latitude classification information
    """
    # House system names for error messages
    system_names = {
        "P": "Placidus",
        "K": "Koch",
        "G": "Gauquelin",
        "O": "Porphyry",
        "E": "Equal",
        "W": "Whole Sign",
        "R": "Regiomontanus",
        "C": "Campanus",
        "T": "Topocentric",
        "B": "Alcabitius",
        "M": "Morinus",
        "X": "Meridian",
        "H": "Horizontal",
        "V": "Vehlow",
        "U": "Krusinski",
        "F": "Carter",
        "Y": "APC",
        "N": "Natural Gradient",
    }

    hsys_char = chr(hsys) if isinstance(hsys, int) else hsys
    fallback_char = (
        chr(fallback_hsys) if isinstance(fallback_hsys, int) else fallback_hsys
    )
    primary_name = system_names.get(hsys_char, hsys_char)
    fallback_name = system_names.get(fallback_char, fallback_char)

    # Check for extreme latitude before calculation
    is_extreme = _is_extreme_latitude(lat)

    try:
        cusps, ascmc = houses(jd_ut, lat, lon, hsys)

        # Validate cusps if requested
        if validate_cusps:
            is_valid, validation_error = _validate_cusps(cusps)
            if not is_valid:
                # Invalid cusps detected - fall back to stable system
                cusps, ascmc = houses(jd_ut, lat, lon, fallback_hsys)
                warning = (
                    f"{primary_name} house system produced invalid cusps at latitude "
                    f"{abs(lat):.2f}° ({validation_error}). Using {fallback_name} as fallback."
                )
                return cusps, ascmc, True, warning

        # Generate warning for extreme latitudes even if calculation succeeded
        if is_extreme:
            info = get_extreme_latitude_info(lat)
            if hsys_char in info["unstable_systems"]:
                warning = (
                    f"{primary_name} house system may have reduced accuracy at "
                    f"extreme latitude {abs(lat):.2f}°{info['hemisphere']}. "
                    f"Consider using a stable system like Porphyry, Equal, or Whole Sign."
                )
                return cusps, ascmc, False, warning

        return cusps, ascmc, False, None
    except PolarCircleError as e:
        # Use fallback house system
        cusps, ascmc = houses(jd_ut, lat, lon, fallback_hsys)

        warning = (
            f"{primary_name} house system unavailable at latitude {abs(lat):.2f}° "
            f"(polar circle threshold: {e.threshold:.2f}°). "
            f"Using {fallback_name} as fallback."
        )
        return cusps, ascmc, True, warning


def houses_armc_with_fallback(
    armc: float,
    lat: float,
    eps: float,
    hsys: int,
    fallback_hsys: int = ord("O"),
    validate_cusps: bool = True,
) -> tuple[tuple[float, ...], tuple[float, ...], bool, str | None]:
    """
    Calculate house cusps from ARMC with automatic fallback for polar latitudes.

    Similar to houses_with_fallback, but calculates from ARMC instead of
    Julian Day. Includes the same extreme latitude handling and cusp validation.

    Args:
        armc: Right Ascension of Medium Coeli in degrees (0-360)
        lat: Geographic latitude in degrees (positive North, negative South)
        eps: True obliquity of the ecliptic in degrees
        hsys: Primary house system identifier (e.g., ord('P') for Placidus)
        fallback_hsys: Fallback house system for polar latitudes.
                       Default: ord('O') (Porphyry).
        validate_cusps: If True (default), validate cusp values for numerical sanity
                        and fall back if invalid cusps are detected.

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes (houses 1-12) in degrees
            - ascmc: Tuple of 8 angles
            - used_fallback: True if fallback system was used
            - warning_message: Informational message if fallback was used or
                               extreme latitude warning, else None

    See Also:
        houses_armc: Standard function that raises PolarCircleError at polar latitudes
        get_extreme_latitude_info: Returns detailed latitude classification information
    """
    # House system names for error messages
    system_names = {
        "P": "Placidus",
        "K": "Koch",
        "G": "Gauquelin",
        "O": "Porphyry",
        "E": "Equal",
        "W": "Whole Sign",
        "R": "Regiomontanus",
        "C": "Campanus",
        "T": "Topocentric",
        "B": "Alcabitius",
        "M": "Morinus",
        "X": "Meridian",
        "H": "Horizontal",
        "V": "Vehlow",
        "U": "Krusinski",
        "F": "Carter",
        "Y": "APC",
        "N": "Natural Gradient",
    }

    hsys_char = chr(hsys) if isinstance(hsys, int) else hsys
    fallback_char = (
        chr(fallback_hsys) if isinstance(fallback_hsys, int) else fallback_hsys
    )
    primary_name = system_names.get(hsys_char, hsys_char)
    fallback_name = system_names.get(fallback_char, fallback_char)

    # Check for extreme latitude before calculation
    is_extreme = _is_extreme_latitude(lat)

    try:
        cusps, ascmc = houses_armc(armc, lat, eps, hsys)

        # Validate cusps if requested
        if validate_cusps:
            is_valid, validation_error = _validate_cusps(cusps)
            if not is_valid:
                # Invalid cusps detected - fall back to stable system
                cusps, ascmc = houses_armc(armc, lat, eps, fallback_hsys)
                warning = (
                    f"{primary_name} house system produced invalid cusps at latitude "
                    f"{abs(lat):.2f}° ({validation_error}). Using {fallback_name} as fallback."
                )
                return cusps, ascmc, True, warning

        # Generate warning for extreme latitudes even if calculation succeeded
        if is_extreme:
            info = get_extreme_latitude_info(lat, eps)
            if hsys_char in info["unstable_systems"]:
                warning = (
                    f"{primary_name} house system may have reduced accuracy at "
                    f"extreme latitude {abs(lat):.2f}°{info['hemisphere']}. "
                    f"Consider using a stable system like Porphyry, Equal, or Whole Sign."
                )
                return cusps, ascmc, False, warning

        return cusps, ascmc, False, None
    except PolarCircleError as e:
        # Use fallback house system
        cusps, ascmc = houses_armc(armc, lat, eps, fallback_hsys)

        warning = (
            f"{primary_name} house system unavailable at latitude {abs(lat):.2f}° "
            f"(polar circle threshold: {e.threshold:.2f}°). "
            f"Using {fallback_name} as fallback."
        )
        return cusps, ascmc, True, warning


def houses_armc(
    armc: float, lat: float, eps: float, hsys: int = ord("P"), ascmc9: float = 0.0
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Calculate house cusps and angles from ARMC (Right Ascension of Medium Coeli).

    This function calculates house cusps directly from the ARMC value instead of
    from a Julian Day. This is useful when you have a pre-calculated ARMC or when
    working with house systems that depend only on ARMC, latitude, and obliquity.

    Reference API compatible function (houses_armc equivalent).

    Args:
        armc: Right Ascension of Medium Coeli in degrees (0-360)
        lat: Geographic latitude in degrees (positive North, negative South)
        eps: True obliquity of the ecliptic in degrees
        hsys: House system identifier (e.g., ord('P') for Placidus, ord('K') for Koch)

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes (houses 1-12) in degrees
            - ascmc: Tuple of 8 angles: [Asc, MC, ARMC, Vertex, EquAsc, CoAsc, CoAscKoch, PolarAsc]

    House Systems Supported:
        'P' = Placidus, 'K' = Koch, 'R' = Regiomontanus, 'C' = Campanus,
        'E'/'A' = Equal, 'W' = Whole Sign, 'O' = Porphyry, 'B' = Alcabitius,
        'T' = Topocentric, 'M' = Morinus, 'X' = Meridian, 'H' = Horizontal,
        'V' = Vehlow, 'G' = Gauquelin, 'U' = Krusinski, 'F' = Carter,
        'Y' = APC, 'N' = Natural Gradient

    Example:
        >>> # Calculate obliquity for J2000.0
        >>> eps = 23.4393  # approximate true obliquity
        >>> armc = 292.957  # ARMC in degrees
        >>> cusps, ascmc = houses_armc(armc, 41.9, eps, ord('P'))
        >>> asc, mc = ascmc[0], ascmc[1]
    """
    # Validate latitude range (must be in [-90, 90])
    from .exceptions import validate_latitude

    validate_latitude(lat, "houses_armc")

    # Normalize ARMC to 0-360
    armc_deg = armc % 360.0

    # Convert house system identifier to character
    hsys_char: str
    if isinstance(hsys, int):
        hsys_char = chr(hsys)
    elif isinstance(hsys, bytes):
        hsys_char = hsys.decode("utf-8")
    else:
        hsys_char = str(hsys)

    # Determine if we need to flip MC (and thus ARMC) for specific systems
    # Regiomontanus (R), Campanus (C), Polich/Page (T) flip MC if below horizon.
    armc_active = armc_deg

    if hsys_char in ["R", "C", "T", "J"]:
        # Check altitude of MC calculated from original ARMC
        mc_dec_rad = math.atan(
            math.sin(math.radians(armc_deg)) * math.tan(math.radians(eps))
        )
        lat_rad = math.radians(lat)
        sin_alt = math.sin(lat_rad) * math.sin(mc_dec_rad) + math.cos(
            lat_rad
        ) * math.cos(mc_dec_rad)

        if sin_alt < 0:
            armc_active = (armc_deg + 180.0) % 360.0

    # Calculate MC from armc_active (flipped if needed)
    mc_rad = math.atan2(
        math.tan(math.radians(armc_active)), math.cos(math.radians(eps))
    )
    mc = math.degrees(mc_rad)
    # Adjust quadrant to match armc_active
    if mc < 0:
        mc += 360.0

    # Quadrant correction: MC should be in the same half of the ecliptic as ARMC
    # ARMC in (90, 270] -> MC should be in (90, 270]
    # ARMC in (270, 360] or [0, 90] -> MC should be in (270, 360] or [0, 90]
    if 90.0 < armc_active <= 270.0:
        if mc <= 90.0 or mc > 270.0:
            mc += 180.0
    elif armc_active > 270.0:
        if mc <= 270.0:
            mc += 180.0
    elif armc_active <= 90.0:
        if mc > 90.0:
            mc += 180.0

    mc = mc % 360.0

    # Ascendant uses armc_deg (Original)
    num = math.cos(math.radians(armc_deg))
    den = -(
        math.sin(math.radians(armc_deg)) * math.cos(math.radians(eps))
        + math.tan(math.radians(lat)) * math.sin(math.radians(eps))
    )
    asc_rad = math.atan2(num, den)
    asc = math.degrees(asc_rad)
    asc = asc % 360.0

    # Ensure Ascendant is on the Eastern Horizon (Azimuth in [0, 180])
    # We check Azimuth relative to the TRUE ARMC (armc_deg)
    asc_r = math.radians(asc)
    eps_r = math.radians(eps)

    # RA
    y = math.cos(eps_r) * math.sin(asc_r)
    x = math.cos(asc_r)
    ra_r = math.atan2(y, x)
    ra = math.degrees(ra_r) % 360.0

    # Dec
    dec_r = math.asin(math.sin(eps_r) * math.sin(asc_r))

    # Hour Angle using TRUE ARMC
    h_deg = (armc_deg - ra + 360.0) % 360.0
    h_r = math.radians(h_deg)

    # Azimuth
    # tan(Az) = sin(H) / (sin(lat)cos(H) - cos(lat)tan(Dec))
    lat_r = math.radians(lat)

    num_az = math.sin(h_r)
    den_az = math.sin(lat_r) * math.cos(h_r) - math.cos(lat_r) * math.tan(dec_r)
    az_r = math.atan2(num_az, den_az)
    az = math.degrees(az_r)
    az = (az + 180.0) % 360.0

    # Check if H is West (0-180). If so, Asc is Setting (Descendant).
    # We want Rising.
    if 0.0 < h_deg < 180.0:
        asc = (asc + 180.0) % 360.0

    # Ensure ASC is in [0, 360) range (handle floating-point near-360 values)
    if asc >= 360.0 or (360.0 - asc) < 1e-10:
        asc = 0.0

    # Vertex uses armc_deg (Original)
    # Hemisphere check relative to TRUE ARMC (West of True ARMC)
    # The horizontal system computes its equator-degenerate angles on
    # the negative-latitude side (reference behavior — its vertex and
    # coasc2 at lat 0 equal the lat -> 0- limits, other systems the
    # lat -> 0+ limits).
    _vtx_lat = -1e-9 if (hsys_char == "H" and abs(lat) < 1e-10) else lat
    vertex = _calc_vertex(armc_deg, eps, _vtx_lat, armc_deg)

    # Equatorial Ascendant (East Point)
    # This is the intersection of the ecliptic with the celestial equator in the east
    # It's the ecliptic longitude where RA = ARMC + 90°
    equ_asc_ra = (armc_deg + 90.0) % 360.0
    # Convert RA to ecliptic longitude
    equ_asc_ra_r = math.radians(equ_asc_ra)
    eps_r = math.radians(eps)
    y = math.sin(equ_asc_ra_r)
    x = math.cos(equ_asc_ra_r) * math.cos(eps_r)
    equ_asc = math.degrees(math.atan2(y, x)) % 360.0

    # Co-Ascendant W. Koch (coasc1)
    # coasc1 = Asc(ARMC - 90°, latitude) + 180°
    coasc_armc = (armc_deg - 90.0) % 360.0
    co_asc_koch = _calc_ascendant(coasc_armc, eps, lat, lat)

    # Add 180° to get opposite point
    co_asc_koch = (co_asc_koch + 180.0) % 360.0

    # Co-Ascendant M. Munkasey (coasc2)
    # If lat >= 0: coasc2 = Asc(ARMC + 90°, 90° - lat)
    # If lat < 0:  coasc2 = Asc(ARMC + 90°, -90° - lat)
    # At equator (lat=0), coasc2_lat becomes 90° which is undefined.
    # The one-sided limits are 180 (lat -> 0+) and 0 (lat -> 0-).
    # The reference outputs 180.0 at exactly lat 0 for every system
    # except the horizontal system 'H', whose internal latitude
    # handling lands on the negative-side limit 0.0 (verified per
    # system against pyswisseph at armc 0/90/180/270).
    if abs(lat) < 1e-10:
        co_asc = 0.0 if hsys_char == "H" else 180.0
    else:
        coasc2_armc = (armc_deg + 90.0) % 360.0
        if lat >= 0:
            coasc2_lat = 90.0 - lat
        else:
            coasc2_lat = -90.0 - lat
        co_asc = _calc_ascendant(coasc2_armc, eps, coasc2_lat, coasc2_lat)

    # Polar Ascendant M. Munkasey (polasc)
    # polasc = Asc(ARMC - 90°, latitude)
    # Note: This is the same as coasc1 but WITHOUT the +180°
    polar_asc = _calc_ascendant(coasc_armc, eps, lat, lat)

    # Build ASCMC array with 8 elements (reference API compatible)
    ascmc = [0.0] * 8
    ascmc[0] = asc
    ascmc[1] = mc
    ascmc[2] = armc_deg
    ascmc[3] = vertex
    ascmc[4] = equ_asc
    ascmc[5] = co_asc_koch  # coasc1 (W. Koch) at index 5
    ascmc[6] = co_asc  # coasc2 (M. Munkasey) at index 6
    ascmc[7] = polar_asc

    # House Cusps
    # Use armc_active for house calculations
    # If MC was flipped, we might need to flip latitude for intermediate cusps
    calc_lat = lat
    if armc_active != armc_deg:
        # MC was flipped. Flip latitude for intermediate cusp calculations.
        calc_lat = -lat

    # Check for polar circle condition for Placidus/Koch/Gauquelin
    # These systems cannot be calculated when abs(lat) + eps > 90°
    # Raise detailed PolarCircleError with useful information
    if hsys_char in ["P", "K", "G"] and _is_polar_circle(lat, eps):
        _raise_polar_circle_error(lat, eps, hsys_char, "houses_armc")

    cusps = [0.0] * 13

    if hsys_char == "P":  # Placidus
        cusps = _houses_placidus(armc_active, lat, eps, asc, mc)
    elif hsys_char == "K":  # Koch
        cusps = _houses_koch(armc_active, lat, eps, asc, mc)
    elif hsys_char == "R":  # Regiomontanus
        cusps = _houses_regiomontanus(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "C":  # Campanus
        cusps = _houses_campanus(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "E":  # Equal (Ascendant)
        cusps = _houses_equal(asc)
    elif hsys_char == "A":  # Equal (MC)
        cusps = _houses_equal(asc)  # Equal MC uses the same algorithm as Equal Asc
    elif hsys_char == "W":  # Whole Sign
        cusps = _houses_whole_sign(asc)
    elif hsys_char == "O":  # Porphyry
        cusps = _houses_porphyry(asc, mc)
    elif hsys_char == "B":  # Alcabitius
        cusps = _houses_alcabitius(armc_active, lat, eps, asc, mc)
    elif hsys_char == "T":  # Polich/Page (Topocentric)
        cusps = _houses_polich_page(armc_active, calc_lat, eps, asc, mc)
    elif hsys_char == "M":  # Morinus
        cusps = _houses_morinus(armc_active, lat, eps, asc, mc)
    elif hsys_char == "X":  # Meridian (Axial)
        cusps = _houses_meridian(armc_active, lat, eps, asc, mc)
    elif hsys_char == "V":  # Vehlow
        cusps = _houses_vehlow(asc)
    elif hsys_char == "H":  # Horizontal (Azimuthal)
        cusps = _houses_horizontal(armc_active, lat, eps, asc, mc)
    elif hsys_char == "Y":  # APC Houses
        cusps = _houses_apc(armc_active, lat, eps, asc, mc)
        # APC at polar latitudes needs MC flipped in ascmc if MC is below horizon
        # (reference behavior - different from R/C/T which flip armc_active)
        mc_dec_rad = math.atan(
            math.sin(math.radians(armc_deg)) * math.tan(math.radians(eps))
        )
        lat_rad = math.radians(lat)
        sin_alt = math.sin(lat_rad) * math.sin(mc_dec_rad) + math.cos(
            lat_rad
        ) * math.cos(mc_dec_rad)
        if sin_alt < 0:
            # Flip MC in ascmc
            ascmc[1] = (ascmc[1] + 180.0) % 360.0
            # Flip all cusps by 180°
            for i in range(1, 13):
                cusps[i] = (cusps[i] + 180.0) % 360.0
    elif hsys_char == "F":  # Carter (Poli-Equatorial)
        cusps = _houses_carter(armc_active, lat, eps, asc, mc)
    elif hsys_char == "U":  # Krusinski
        cusps = _houses_krusinski(armc_active, lat, eps, asc, mc)
    elif hsys_char == "N":  # Natural Gradient
        cusps = _houses_natural_gradient(armc_active, lat, eps, asc, mc)
    elif hsys_char == "G":  # Gauquelin
        cusps = _houses_gauquelin(armc_active, lat, eps, asc, mc)
    elif hsys_char == "S":  # Sripati
        cusps = _houses_sripati(asc, mc)
    elif hsys_char == "L":  # Pullen SD (Sinusoidal Delta / Neo-Porphyry)
        cusps = _houses_pullen_sd(asc, mc)
    elif hsys_char == "Q":  # Pullen SR (Sinusoidal Ratio)
        cusps = _houses_pullen_sr(asc, mc)
    elif hsys_char == "D":  # Equal from MC
        cusps = _houses_equal_mc(asc, mc)
    elif hsys_char == "I":  # Sunshine (Treindl)
        # ascmc9 carries the Sun's declination (reference convention)
        cusps = _houses_sunshine(armc_active, lat, eps, asc, mc, ascmc9)
        ascmc[1] = cusps[10]
    elif hsys_char == "i":  # Sunshine (Makransky)
        cusps = _houses_sunshine_makransky(armc_active, lat, eps, asc, mc, ascmc9)
        ascmc[1] = cusps[10]
    elif hsys_char == "J":  # Savard-A
        cusps = _houses_savard_a(armc_active, calc_lat, eps, asc, mc)
    else:
        # Default to Placidus
        cusps = _houses_placidus(armc_active, lat, eps, asc, mc)

    # Return cusps array (reference API compatible: no padding at index 0)
    # For Gauquelin ('G'), return 36 sectors; otherwise return 12 houses
    if hsys_char == "G":
        return tuple(cusps[1:37]), tuple(ascmc)
    return tuple(cusps[1:13]), tuple(ascmc)


def houses_armc_ex2(
    armc: float,
    lat: float,
    eps: float,
    hsys: int = ord("P"),
    ascmc9: float = 0.0,
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    """
    Extended house calculation from ARMC returning cusps, angles, and their velocities.

    This function combines houses_armc() with velocity calculations similar to
    houses_ex2(). It calculates house cusps directly from the ARMC value and
    also returns the velocities (derivatives) of house cusps and angles.

    Velocities are always calculated, matching pyswisseph behavior.
    The ``ascmc9`` parameter is accepted for API compatibility (used by the
    Sunshine house system) but is otherwise unused.

    Velocities are calculated using centered finite differences, with ARMC
    shifted by ±1 second (Koch/Placidus) or ±1 minute (other systems).

    Note: The sidereal flag (FLG_SIDEREAL) has no effect since sidereal mode
    requires a Julian Day for ayanamsa calculation, which is not available when
    using ARMC directly.

    Args:
        armc: Right Ascension of Medium Coeli in degrees (0-360)
        lat: Geographic latitude in degrees (positive North, negative South)
        eps: True obliquity of the ecliptic in degrees
        hsys: House system identifier (e.g., ord('P') for Placidus, ord('K') for Koch)
        ascmc9: Optional parameter for Sunshine house system (default 0.0)

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes in degrees
            - ascmc: Tuple of 8 angles (Asc, MC, ARMC, Vertex, EquAsc, CoAsc, CoAscKoch, PolarAsc)
            - cusps_speed: Tuple of 12 house cusp velocities in degrees/day
            - ascmc_speed: Tuple of 8 angle velocities in degrees/day

    Example:
        >>> eps = 23.4393  # true obliquity
        >>> armc = 292.957  # ARMC in degrees
        >>> cusps, ascmc, cusps_speed, ascmc_speed = houses_armc_ex2(
        ...     armc, 41.9, eps, ord('P')
        ... )
        >>> # cusps_speed[0] is the velocity of the 1st house cusp (same as ASC)
    """
    # Calculate positions at current ARMC
    cusps, ascmc = houses_armc(armc, lat, eps, hsys)

    # Always calculate velocities (matching pyswisseph behavior).
    # Compute d(cusp)/d(ARMC) via centered finite differences, then
    # scale by the sidereal rotation rate to obtain deg/day.
    #
    # The Earth completes 360.98564736629° of ARMC per mean solar day
    # (one extra rotation relative to the Sun).  Previous code used
    # 360°/day (solar rate), which under-estimated all speeds by the
    # ratio 360/360.9856 ≈ 0.27 %.
    _SIDEREAL_RATE = 360.98564736629  # ARMC degrees per mean solar day

    # ARMC step for the finite difference.
    # Koch and Placidus use a 1-second step (nested trig amplifies
    # truncation error at the 1-minute step).
    # All other systems use a 1-minute step.
    if hsys in (ord("K"), ord("P")):
        d_armc = _SIDEREAL_RATE / 86400.0  # sidereal degrees per 1 second
    else:
        d_armc = _SIDEREAL_RATE / 1440.0  # sidereal degrees per 1 minute

    # Calculate positions at ARMC ± d_armc
    cusps_before, ascmc_before = houses_armc(armc - d_armc, lat, eps, hsys)
    cusps_after, ascmc_after = houses_armc(armc + d_armc, lat, eps, hsys)

    def angular_diff_local(pos2: float, pos1: float) -> float:
        """Calculate angular difference handling 360° wraparound."""
        diff = pos2 - pos1
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        return diff

    # Velocities in deg/day = d(cusp)/d(ARMC) * (ARMC deg/day)
    # Since d_armc already equals _SIDEREAL_RATE * dt, dividing by
    # (2*d_armc) and then multiplying by _SIDEREAL_RATE is equivalent to
    # dividing by (2*dt).  Using d_armc directly avoids a separate dt
    # variable: speed = Δcusp / (2*d_armc) * _SIDEREAL_RATE.
    cusps_speed = tuple(
        angular_diff_local(cusps_after[i], cusps_before[i])
        / (2 * d_armc)
        * _SIDEREAL_RATE
        for i in range(len(cusps))
    )

    ascmc_speed = tuple(
        angular_diff_local(ascmc_after[i], ascmc_before[i])
        / (2 * d_armc)
        * _SIDEREAL_RATE
        for i in range(len(ascmc))
    )

    # ── System-specific cusp speed overrides ──────────────────────
    if hsys == ord("W"):
        # Whole Sign: cusps are at fixed sign boundaries (0°, 30°, …).
        # Most cusps have zero speed (they jump discontinuously).
        # Cusps 1,7 (ASC/DESC) get ASC speed; cusps 4,10 (IC/MC) get
        # MC speed — matching pyswisseph behaviour.
        v_asc = ascmc_speed[0]
        v_mc = ascmc_speed[1]
        cs = [0.0] * len(cusps)
        cs[0] = v_asc  # cusp 1  = ASC
        cs[3] = v_mc  # cusp 4  = IC
        cs[6] = v_asc  # cusp 7  = DESC
        cs[9] = v_mc  # cusp 10 = MC
        cusps_speed = tuple(cs)
    # Other systems (Placidus, Koch, Porphyry, etc.): use numerical
    # differentiation directly. Note: Koch and Placidus use a smaller
    # step size (1 second vs 1 minute) to reduce truncation error.

    return cusps, ascmc, cusps_speed, ascmc_speed


def houses_ex(
    tjdut: float, lat: float, lon: float, hsys: int = ord("P"), flags: int = 0
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Extended house calculation with sidereal zodiac support.

    Similar to houses() but applies ayanamsa correction when FLG_SIDEREAL
    flag is set, converting tropical positions to sidereal.

    For Ascendant-based house systems (Whole Sign, Equal, Vehlow), the cusps
    are recalculated using the sidereal Ascendant to ensure proper sign alignment.

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        lat: Geographic latitude in degrees
        lon: Geographic longitude in degrees
        hsys: House system identifier (int or bytes)
        flags: Calculation flags bitmask (e.g., FLG_SIDEREAL, FLG_MOSEPH).
               Use FLG_MOSEPH for extended date range (-3000 to +3000 CE).
               The ephemeris flag is propagated to all internal calculations.

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes in degrees
            - ascmc: Tuple of 8 angles (Asc, MC corrected if sidereal)

    Example:
        >>> from libephemeris import set_sid_mode, SIDM_LAHIRI
        >>> set_sid_mode(SIDM_LAHIRI)
        >>> cusps, ascmc = houses_ex(2451545.0, 51.5, -0.12, ord('P'), FLG_SIDEREAL)
    """
    # Propagate flags to houses() so ephemeris flags (FLG_MOSEPH etc.) are used
    cusps, ascmc = houses(tjdut, lat, lon, hsys, flags)

    if flags & FLG_SIDEREAL:
        ayanamsa = get_ayanamsa_ut(tjdut)

        # Compute sidereal angles
        # All ecliptic longitudes in ascmc get ayanamsa correction EXCEPT
        # ascmc[2] (ARMC) which is an equatorial/geometric quantity.
        ascmc_list = list(ascmc)
        sid_asc = (ascmc_list[0] - ayanamsa) % 360.0
        sid_mc = (ascmc_list[1] - ayanamsa) % 360.0
        ascmc_list[0] = sid_asc
        ascmc_list[1] = sid_mc
        # ascmc[2] = ARMC — geometric, NOT corrected
        ascmc_list[3] = (ascmc_list[3] - ayanamsa) % 360.0  # Vertex
        ascmc_list[4] = (ascmc_list[4] - ayanamsa) % 360.0  # EquAsc
        ascmc_list[5] = (ascmc_list[5] - ayanamsa) % 360.0  # CoAsc Koch
        ascmc_list[6] = (ascmc_list[6] - ayanamsa) % 360.0  # CoAsc Munkasey
        ascmc_list[7] = (ascmc_list[7] - ayanamsa) % 360.0  # PolarAsc
        ascmc = tuple(ascmc_list)

        # Normalize hsys to a character for comparison
        if isinstance(hsys, int):
            hsys_char = chr(hsys)
        elif isinstance(hsys, bytes):
            hsys_char = hsys.decode("utf-8")
        else:
            hsys_char = str(hsys)

        # For Ascendant-based house systems, recalculate using sidereal Ascendant
        # Whole Sign (W), Equal (A/E), Vehlow (V)
        if hsys_char == "W":
            # Whole Sign: each house is 0° of consecutive signs from Asc sign
            # cusps[0] = House 1, cusps[1] = House 2, etc.
            start = math.floor(sid_asc / 30.0) * 30.0
            cusps = tuple([(start + i * 30.0) % 360.0 for i in range(12)])
        elif hsys_char in ("A", "E"):
            # Equal: 30° divisions starting from sidereal Asc
            cusps = tuple([(sid_asc + i * 30.0) % 360.0 for i in range(12)])
        elif hsys_char == "V":
            # Vehlow Equal: sidereal Asc at middle of 1st house
            start = (sid_asc - 15.0) % 360.0
            cusps = tuple([(start + i * 30.0) % 360.0 for i in range(12)])
        elif hsys_char in ("I", "i"):
            # Sunshine: Recalculate using sidereal Asc/MC
            # 'I' = Treindl algorithm, 'i' = Makransky algorithm
            try:
                eph_flags = flags & (FLG_JPLEPH | FLG_SWIEPH)
                sun_pos, _ = calc_ut(tjdut, SUN, FLG_EQUATORIAL | eph_flags)
                sun_dec = sun_pos[1]
            except (IndexError, TypeError, ValueError):
                sun_dec = 0.0
            ts = get_timescale()
            t = ts.ut1_jd(tjdut)
            eps = get_true_obliquity(t.tt)
            armc = ascmc[2]
            if hsys_char == "I":
                sunshine_cusps = _houses_sunshine(
                    armc, lat, eps, sid_asc, sid_mc, sun_dec
                )
            else:
                sunshine_cusps = _houses_sunshine_makransky(
                    armc, lat, eps, sid_asc, sid_mc, sun_dec
                )
            # Angular cusps (1, 4, 7, 10) are already correct with sid_asc/sid_mc
            # Intermediate cusps (2, 3, 5, 6, 8, 9, 11, 12) are calculated
            # geometrically and need to be converted to sidereal
            intermediate_houses = {2, 3, 5, 6, 8, 9, 11, 12}  # 1-indexed
            cusps = tuple(
                (sunshine_cusps[i] - ayanamsa) % 360.0
                if i in intermediate_houses
                else sunshine_cusps[i]
                for i in range(1, 13)
            )
        else:
            # For other systems, just subtract ayanamsa from tropical cusps
            cusps = tuple([(c - ayanamsa) % 360.0 for c in cusps])

    return cusps, ascmc


def houses_ex2(
    tjdut: float, lat: float, lon: float, hsys: int = ord("P"), flags: int = 0
) -> tuple[tuple[float, ...], tuple[float, ...], tuple[float, ...], tuple[float, ...]]:
    """
    Extended house calculation returning cusps, angles, and their velocities.

    This function is an extended version of houses_ex() that also returns
    the velocities (derivatives) of house cusps and angles.

    Velocities are only calculated when the FLG_SPEED flag is set in the flags
    parameter. When FLG_SPEED is not set, zero velocities are returned for
    efficiency. This is useful for progressed chart applications where the rate
    of change of house cusps is needed.

    Velocities are computed via the ARMC-based derivative path
    (houses_armc_ex2), which varies ARMC with fixed obliquity and
    scales by the sidereal rotation rate (~360.986°/day).

    Args:
        tjdut: Julian Day in Universal Time (UT1)
        lat: Geographic latitude in degrees
        lon: Geographic longitude in degrees
        hsys: House system identifier (int or bytes)
        flags: Calculation flags bitmask. Use FLG_SPEED to compute velocities.
               FLG_SIDEREAL can also be used for sidereal calculations.

    Returns:
        Tuple containing:
            - cusps: Tuple of 12 house cusp longitudes in degrees
            - ascmc: Tuple of 8 angles (Asc, MC, etc.)
            - cusps_speed: Tuple of 12 house cusp velocities in degrees/day (0.0 if FLG_SPEED not set)
            - ascmc_speed: Tuple of 8 angle velocities in degrees/day (0.0 if FLG_SPEED not set)

    Example:
        >>> cusps, ascmc, cusps_speed, ascmc_speed = houses_ex2(
        ...     2451545.0, 41.9, 12.5, ord('P'), FLG_SPEED
        ... )
        >>> # cusps_speed[0] is the velocity of the 1st house cusp (ASC)
    """
    # Calculate positions at current time
    cusps, ascmc = houses_ex(tjdut, lat, lon, hsys, flags)

    # Always calculate velocities (matching pyswisseph behavior).
    # Delegate speed computation to the ARMC-based path.
    # This varies ARMC (with fixed obliquity) and scales by the
    # sidereal rotation rate, matching the internal approach used by
    # pyswisseph's houses_ex2.  Direct JD-based finite differences
    # mix ARMC, obliquity, and nutation changes, producing systematic
    # ~0.003 deg/day offsets on angular cusps.
    #
    # We extract ARMC and true obliquity from the ascmc tuple returned
    # by houses_ex (index 2 = ARMC) and compute obliquity via the
    # same cached path used by houses().
    armc_val = ascmc[2]  # ARMC stored by houses
    ts = get_timescale()
    t = ts.ut1_jd(tjdut)
    eps = get_true_obliquity(t.tt)

    _, _, cusps_speed, ascmc_speed = houses_armc_ex2(armc_val, lat, eps, hsys)

    return cusps, ascmc, cusps_speed, ascmc_speed


def _houses_with_context(
    tjdut: float, lat: float, lon: float, hsys: int, ctx
) -> tuple[tuple[float, ...], tuple[float, ...]]:
    """
    Calculate houses using an explicit EphemerisContext (thread-safe).

    Thread-safe wrapper around houses that uses context state.

    Args:
        tjdut: Julian Day in Universal Time
        lat: Geographic latitude in degrees
        lon: Geographic longitude in degrees
        hsys: House system identifier
        ctx: EphemerisContext instance

    Returns:
        Same as houses: (cusps, ascmc)

    Thread Safety:
        This function acquires state._CONTEXT_SWAP_LOCK to ensure that the
        save-set-restore cycle is atomic across threads.
    """
    from . import state

    with state._CONTEXT_SWAP_LOCK:
        # Save current global state
        old_sid_mode = state._SIDEREAL_MODE
        old_sid_t0 = state._SIDEREAL_T0
        old_sid_ayan_t0 = state._SIDEREAL_AYAN_T0

        try:
            # Temporarily set global state from context
            state._SIDEREAL_MODE = ctx.sidereal_mode
            state._SIDEREAL_T0 = ctx.sidereal_t0
            state._SIDEREAL_AYAN_T0 = ctx.sidereal_ayan_t0

            # Use existing house calculation logic
            return houses(tjdut, lat, lon, hsys)
        finally:
            # Restore global state
            state._SIDEREAL_MODE = old_sid_mode
            state._SIDEREAL_T0 = old_sid_t0
            state._SIDEREAL_AYAN_T0 = old_sid_ayan_t0


def house_name(hsys: int) -> str:
    """
    Get the name of a house system.
    """
    hsys_char: str
    if isinstance(hsys, int):
        hsys_char = chr(hsys)
    elif isinstance(hsys, bytes):
        hsys_char = hsys.decode("utf-8")
    else:
        hsys_char = str(hsys)

    names = {
        "P": "Placidus",
        "K": "Koch",
        "O": "Porphyry",
        "R": "Regiomontanus",
        "C": "Campanus",
        "E": "equal",
        "A": "equal",
        "W": "equal/ whole sign",
        "M": "Morinus",
        "B": "Alcabitius",
        "T": "Polich/Page",
        "U": "Krusinski-Pisa-Goelzer",
        "G": "Gauquelin sectors",
        "V": "equal/Vehlow",
        "X": "axial rotation system/Meridian houses",
        "H": "horizon/azimut",
        "F": "Carter poli-equ.",
        "S": "Sripati",
        "L": "Pullen SD",
        "Q": "Pullen SR",
        "N": "equal/1=Aries",
        "Y": "APC houses",
        "D": "equal (MC)",
        "I": "Sunshine",
        "i": "Sunshine/alt.",
        "J": "Savard-A",
    }
    return names.get(hsys_char, "Unknown")


def _houses_placidus(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Placidus house system (time-based divisions of diurnal/nocturnal arcs).

    Most popular house system in modern Western astrology. Divides the time a point
    takes to travel from horizon to meridian (and meridian to horizon) into thirds.

    Algorithm:
        1. Trisect semi-diurnal arc (rising to culmination) for houses 11, 12
        2. Trisect semi-nocturnal arc (setting to anti-culmination) for houses 2, 3
        3. Use iterative solution to find ecliptic longitude at each time division
        4. Calculate opposite houses by adding 180°

    Mathematical Formulas:
        Ascensional Difference: AD = arcsin(tan(φ) · tan(δ))
        Semi-diurnal Arc: SA = 90° + AD (above horizon)
        Semi-nocturnal Arc: SA = 90° - AD (below horizon)

        House 11: H₁₁ = SA/3 = (90° + AD)/3; RA₁₁ = ARMC + H₁₁
        House 12: H₁₂ = 2·SA/3 = 2·(90° + AD)/3; RA₁₂ = ARMC + H₁₂
        House 2:  H₂ = 2·(90° - AD)/3; RA₂ = ARMC + 180° - H₂
        House 3:  H₃ = (90° - AD)/3; RA₃ = ARMC + 180° - H₃

        Iterative convergence (threshold 1e-7°):
            1. tan(δ) = sin(RA) · tan(ε)
            2. AD = arcsin(tan(φ) · tan(δ))
            3. RA_new = ARMC ± H(AD)
            4. Repeat until |RA_new - RA| < 1e-7°

        RA to ecliptic: λ = atan2(sin(RA), cos(RA) · cos(ε))

    Note: Polar latitude handling
        Placidus is undefined at latitudes > ~66° where some ecliptic points never
        rise/set (circumpolar behavior). At polar latitudes, houses() raises
        PolarCircleError with detailed information about the threshold and
        recommended alternatives. Use houses_with_fallback() for automatic
        fallback to Porphyry. The polar threshold is approximately 90° - obliquity.
        Internal fallback to Porphyry is retained as a safety net.

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees (house 1)
        mc: Midheaven longitude in degrees (house 10)

    Returns:
        List of 13 house cusp longitudes (index 0 is 0.0, indices 1-12 are cusps)
    """
    # Actually, standard Placidus:
    # 11, 12 are in SE quadrant (MC to Asc).
    # 2, 3 are in NE quadrant (Asc to IC).

    cusps = _init_cardinal_cusps(asc, mc)

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    # Helper to find cusp
    # factor: 1.0/3.0 for 11 and 3, 2.0/3.0 for 12 and 2
    # quadrant: 1 for 11/12, 2 for 2/3?
    # Actually, we solve for RA.

    def iterate_placidus(offset_deg, is_below_horizon):
        # Initial guess: RAMC + offset
        # For 11: offset = 30
        # For 12: offset = 60
        # For 2: offset = 120
        # For 3: offset = 150

        ra = (armc + offset_deg) % 360.0
        new_ra = ra  # Initialize for type safety
        prev_diff = float("inf")  # Track convergence/divergence
        converged = False

        # Convergence threshold: 1e-7° = 0.00036 arcsec
        CONVERGENCE_THRESHOLD = 1e-7
        # Maximum iterations - if not converged by then, likely diverging
        MAX_ITERATIONS = 50

        # Iterate to convergence (typically 15-20 iterations needed for 1e-7° threshold)
        for iteration in range(MAX_ITERATIONS):
            # Calculate declination for point at this RA on ecliptic
            # Using spherical astronomy formula: tan(dec) = sin(RA) * tan(eps)

            sin_ra = math.sin(math.radians(ra))
            tan_dec = sin_ra * math.tan(rad_eps)
            math.atan(tan_dec)

            # Calculate semi-arc (or part of it)
            # tan(lat) * tan(dec)
            # Within one float-noise step of the polar circle the
            # candidate RA can sweep ecliptic degrees whose declination
            # makes |tan(lat)*tan(dec)| exceed 1 mid-iteration even
            # though the converged cusp itself is well-defined.  Clamp
            # the ascensional-difference argument (AD saturates at
            # +/-90 deg) and keep iterating — the latitude pre-check
            # has already raised for genuinely polar latitudes.
            # Bailing out here sent these cases to the fallback cusps,
            # up to 20 deg away from the reference at lat +/-66.55.
            prod = math.tan(rad_lat) * tan_dec
            if prod > 1.0:
                prod = 1.0
            elif prod < -1.0:
                prod = -1.0

            # AD (Ascensional Difference) = asin(prod)
            # SA (Semi-Arc) = 90 + AD (if decl north and lat north)
            # H = RAMC - RA
            # Placidus condition:
            # H = f * SA?
            # Or sin(H) = f * sin(SA)? No.
            # Placidus: H is proportional to SA.
            # H = (offset / 90) * SA?
            # For 11 (30 deg from MC): H = SA/3.
            # For 12 (60 deg from MC): H = 2*SA/3.

            # Wait, offset 30 means 30 degrees of RA? No.
            # It implies the trisection.
            # Factor f.

            if offset_deg == 30 or offset_deg == 210:
                pass
            if offset_deg == 60 or offset_deg == 240:
                pass
            if offset_deg == 120 or offset_deg == 300:
                pass  # From IC?
            if offset_deg == 150 or offset_deg == 330:
                pass

            # If below horizon (houses 2, 3), semi-arc is nocturnal.
            # SA_noct = 180 - SA_diurnal = 90 - AD.
            # H is measured from IC (RAMC + 180).

            ad = math.asin(prod)

            if is_below_horizon:
                # Houses 2, 3
                # Measured from IC (RAMC + 180)
                # H = f * (90 - AD)
                # RA = IC - H = RAMC + 180 - H?
                # Or RA = IC + H?
                # House 2 is East of IC. RA > IC.
                # RA = RAMC + 180 + f * (90 - AD) ?
                # No, House 2 is "following" IC.
                # Let's stick to standard formula:
                # R = RAMC + 180 + f * (90 + AD)?
                # Note: AD has sign of dec.

                # Let's use the formula:
                # R = RAMC + 180 + f * (90 - AD) ?
                # If lat > 0, dec > 0, AD > 0.
                # Nocturnal arc = 180 - (90 + AD) = 90 - AD.
                # House 2 is 2/3 of way from Asc to IC? No.
                # House 2 is 1/3 of way from IC to Asc? No.
                # House 2 start is 2/3 SA_noct from IC?
                # House 3 start is 1/3 SA_noct from IC?

                # Correct mapping:
                # 11: RAMC + SA/3
                # 12: RAMC + 2*SA/3
                # 2: RAMC + 180 - 2*SA_noct/3 ? No.
                # 2 is after Asc (House 1).
                # Asc is at RAMC + 90 + AD? No.

                # Let's use the standard iterative formula directly:
                # RA_new = RAMC + const + AD? No.

                pass

            # Simplified Placidus Iteration:
            # R(n+1) = RAMC + asin( tan(lat)*tan(dec(Rn)) * factor ) + base_offset
            # Where factor depends on house.
            # House 11: factor = 1/3? No.
            # House 11 condition: sin(RA - RAMC) = tan(dec)*tan(lat)/3 ? No.
            # The condition is on time.
            # H = RA - RAMC.
            # H = SA / 3.
            # SA = 90 + AD.
            # H = 30 + AD/3 ? No.

            # Correct Placidus Formula (from Munkasey):
            # House 11: tan(H) = f * tan(H_Asc)? No.

            # Let's use the "Pole" method which is equivalent.
            # tan(Pole) = sin(HouseAngle) * tan(Lat) ? No.

            # Let's go back to basics:
            # House 11: 1/3 of semi-arc from MC.
            # H = (90 + AD) / 3.
            # RA = RAMC - H (since 11 is East of MC).
            # RA = RAMC - (90 + asin(tan(lat)tan(dec))) / 3.

            # House 12: 2/3 of semi-arc.
            # RA = RAMC - 2*(90 + asin(tan(lat)tan(dec))) / 3.

            # House 2: 2/3 of nocturnal semi-arc from IC (West of IC? No, East).
            # House 2 is below horizon, West of Asc? No, East.
            # 1 -> 2 -> 3 -> 4(IC).
            # House 2 is 1/3 of way from Asc to IC? No.
            # House 2 cusp is 2/3 of semi-arc from IC?
            # SA_noct = 90 - AD.
            # H_from_IC = 2 * SA_noct / 3.
            # RA = (RAMC + 180) - H_from_IC.
            # RA = RAMC + 180 - 2*(90 - AD)/3.

            # House 3: 1/3 of semi-arc from IC.
            # RA = RAMC + 180 - 1*(90 - AD)/3.

            # Let's verify signs.
            # 11 is SE. MC is S. 11 is East of S. RA < RAMC?
            # No, stars rise in East, RA increases Eastwards?
            # RA increases Eastwards.
            # MC is on Meridian.
            # Point East of Meridian has RA > RAMC?
            # H = LST - RA.
            # East of Meridian -> H < 0.
            # So RA > LST.
            # So House 11 (SE) has RA > RAMC?
            # No, House 11 is "before" MC in diurnal motion.
            # Sun is in 11 before it culminates.
            # So RA_Sun > RAMC? No.
            # H = LST - RA.
            # If Sun is East, H is negative (e.g. -2h).
            # RA = LST - H = LST + 2h.
            # So RA > RAMC.

            # So for House 11:
            # H_from_MC = SA / 3.
            # RA = RAMC + H_from_MC = RAMC + (90 + AD)/3.

            # House 12:
            # RA = RAMC + 2*(90 + AD)/3.

            # House 2:
            # Below horizon.
            # H_from_IC = 2 * SA_noct / 3.
            # RA = (RAMC + 180) + H_from_IC.
            # RA = RAMC + 180 + 2*(90 - AD)/3.

            # House 3:
            # RA = RAMC + 180 + 1*(90 - AD)/3.

            # Let's implement this.

            ad = math.asin(prod)
            ad_deg = math.degrees(ad)

            if offset_deg == 30:  # House 11
                h_deg = (90.0 + ad_deg) / 3.0
                new_ra = (armc - h_deg) % 360.0  # Wait, 11 is East of MC.
                # If RA > RAMC, H < 0.
                # H = RAMC - RA.
                # If H is "distance from MC", then RA = RAMC - H (if East).
                # Wait.
                # Sun rises. H = -SA. RA = RAMC + SA.
                # Sun culminates. H = 0. RA = RAMC.
                # House 11 is between Rise and Culminate.
                # So RA should be between RAMC+SA and RAMC.
                # So RA > RAMC.
                # So RA = RAMC + part_of_SA.
                new_ra = (
                    armc + h_deg
                ) % 360.0  # Wait, H is usually defined positive West.
                # If H is positive West, then East is negative H.
                # H = RAMC - RA.
                # RA = RAMC - H.
                # If we want East, H must be negative.
                # H = - (SA/3).
                # RA = RAMC - (-SA/3) = RAMC + SA/3.
                # Correct.

                # But wait, standard Placidus 11th cusp is usually *South-East*.
                # It is 30 degrees "above" Ascendant? No.
                # It is 30 degrees "before" MC?
                # House 10 starts at MC. House 11 starts at Cusp 11.
                # Order: 10, 11, 12, 1.
                # 10 is MC. 1 is Asc.
                # So 11 is between MC and Asc.
                # So RA is between RAMC and RAMC+SA.
                # So RA = RAMC + SA/3?
                # No, 10 is MC. 11 is next.
                # So 11 is "later" in RA?
                # Houses increase in counter-clockwise direction on Ecliptic.
                # 10 (MC) -> 11 -> 12 -> 1 (Asc).
                # MC RA approx 270 (if Aries rising). Asc RA approx 0.
                # So RA increases from 10 to 1.
                # So RA_11 > RA_10.
                # So RA_11 = RAMC + something.
                # Correct.

                new_ra = (armc + h_deg) % 360.0

            elif offset_deg == 60:  # House 12
                h_deg = 2.0 * (90.0 + ad_deg) / 3.0
                new_ra = (armc + h_deg) % 360.0

            elif offset_deg == 120:  # House 2
                # House 2 is after Asc (1).
                # 1 -> 2 -> 3 -> 4 (IC).
                # Asc RA approx 0. IC RA approx 90.
                # So RA increases.
                # RA_2 = RAMC + 180 - something?
                # IC is RAMC + 180.
                # House 2 is "before" IC.
                # So RA_2 < RA_IC.
                # So RA_2 = RAMC + 180 - H_from_IC.
                # H_from_IC = 2 * SA_noct / 3.
                # SA_noct = 90 - AD_deg.
                h_deg = 2.0 * (90.0 - ad_deg) / 3.0
                new_ra = (armc + 180.0 - h_deg) % 360.0

            elif offset_deg == 150:  # House 3
                h_deg = 1.0 * (90.0 - ad_deg) / 3.0
                new_ra = (armc + 180.0 - h_deg) % 360.0

            else:
                # Should never reach here - all valid offsets are handled above
                new_ra = ra  # Fallback to prevent unbound variable error

            # Update RA
            diff = abs(new_ra - ra)
            if diff > 180:
                diff = 360 - diff
            ra = new_ra

            # Check for convergence
            if diff < CONVERGENCE_THRESHOLD:
                converged = True
                break

            # Check for divergence: if diff is increasing, iteration is unstable
            # Allow some oscillation (diff can increase slightly) but detect runaway
            if iteration > 5 and diff > prev_diff * 1.5:
                # Diverging - this happens near polar latitudes
                return None

            prev_diff = diff

        # If we exhausted iterations without converging, return None for fallback
        if not converged:
            return None

        # Converged RA. Find Ecliptic Longitude.
        # tan(lon) = tan(ra) / cos(eps)
        # Use atan2 to handle quadrants correctly
        # sin(lon) ~ sin(ra)
        # cos(lon) ~ cos(ra) * cos(eps)
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))

        return lon % 360.0

    def _ad_deg_at(ra_deg: float) -> float:
        """Clamped ascensional difference (deg) of the ecliptic point at RA."""
        tan_dec = math.sin(math.radians(ra_deg)) * math.tan(rad_eps)
        prod = math.tan(rad_lat) * tan_dec
        if prod > 1.0:
            prod = 1.0
        elif prod < -1.0:
            prod = -1.0
        return math.degrees(math.asin(prod))

    def bisect_placidus(offset_deg: float, is_below_horizon: bool):
        """Solve the Placidus hour-angle condition by bisection.

        The fixed-point iteration above oscillates within float-noise of
        the polar circle (lat +/-66.55 at eps 23.4393), where the
        reference still converges.  The cusp condition is a continuous
        root problem with a guaranteed sign change on its bracket, so
        bisection always lands it:

          house 11: x = (90 + AD(armc + x)) / 3          on x in [0, 60]
          house 12: x = 2 (90 + AD(armc + x)) / 3        on x in [0, 120]
          house 2:  x = 2 (90 - AD(armc + 180 - x)) / 3  on x in [0, 120]
          house 3:  x = (90 - AD(armc + 180 - x)) / 3    on x in [0, 60]

        (x is the hour-angle fraction from the relevant meridian; AD is
        clamped so the bracket endpoints always satisfy g(lo) <= 0 and
        g(hi) >= 0.)
        """
        if offset_deg == 30:
            frac, hi, from_ic = 1.0 / 3.0, 60.0, False
        elif offset_deg == 60:
            frac, hi, from_ic = 2.0 / 3.0, 120.0, False
        elif offset_deg == 120:
            frac, hi, from_ic = 2.0 / 3.0, 120.0, True
        else:  # 150
            frac, hi, from_ic = 1.0 / 3.0, 60.0, True

        def g(x: float) -> float:
            if from_ic:
                ad = _ad_deg_at(armc + 180.0 - x)
                return x - frac * (90.0 - ad)
            ad = _ad_deg_at(armc + x)
            return x - frac * (90.0 + ad)

        lo_x, hi_x = 0.0, hi
        g_lo = g(lo_x)
        if g_lo > 0.0:
            return None  # bracket failed (cannot happen with clamped AD)
        for _ in range(80):
            mid = 0.5 * (lo_x + hi_x)
            if g(mid) <= 0.0:
                lo_x = mid
            else:
                hi_x = mid
        x = 0.5 * (lo_x + hi_x)
        ra = (armc + 180.0 - x) % 360.0 if from_ic else (armc + x) % 360.0
        y = math.sin(math.radians(ra))
        xx = math.cos(math.radians(ra)) * math.cos(rad_eps)
        return math.degrees(math.atan2(y, xx)) % 360.0

    # Calculate cusps: fast fixed point first, bisection where it fails
    c11 = iterate_placidus(30, False)
    c12 = iterate_placidus(60, False)
    c2 = iterate_placidus(120, True)
    c3 = iterate_placidus(150, True)
    if c11 is None:
        c11 = bisect_placidus(30, False)
    if c12 is None:
        c12 = bisect_placidus(60, False)
    if c2 is None:
        c2 = bisect_placidus(120, True)
    if c3 is None:
        c3 = bisect_placidus(150, True)

    if c11 is None or c12 is None or c2 is None or c3 is None:
        # Safety net only: the bisection bracket is mathematically
        # guaranteed with clamped AD.
        return _houses_porphyry(asc, mc)

    cusps[11] = c11
    cusps[12] = c12
    cusps[2] = c2
    cusps[3] = c3

    # Opposites
    cusps[5] = (c11 + 180) % 360.0
    cusps[6] = (c12 + 180) % 360.0
    cusps[8] = (c2 + 180) % 360.0
    cusps[9] = (c3 + 180) % 360.0

    return cusps


def _houses_koch(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Koch (Birthplace/GOH) house system (code 'K').

    The Koch system (Geburtsort-Häusersystem = Birthplace House System) was
    developed by Walter Koch. It divides houses based on the ascensional
    difference of the MC's declination.

    Algorithm derived from:
    Holden, 'The Elements of House Division' (1977), Koch section.

    Derivation from semi-arc concept:
        1. Compute sin(dec_MC) = sin(MC) * sin(eps) (MC declination)
        2. Compute elevation angle c = atan(tan(lat) / cos(dec_MC))
        3. Koch arc = asin(sin(c) * sin(dec_MC)) / 3
           (one-third of the ascensional difference)
        4. Cusps are the ecliptic rising points at ARMC offsets of
           30, 60, 120, 150 degrees, adjusted by multiples of the Koch arc.

    Reference: Meeus, 'Astronomical Algorithms' 2nd ed., Ch. 13
    (spherical trigonometry for rising-point formula).

    Note: Polar latitude handling
        Koch is undefined at latitudes > ~66 deg where some ecliptic points never
        rise/set (circumpolar behavior). At polar latitudes, houses() raises
        PolarCircleError with detailed information about the threshold and
        recommended alternatives. Use houses_with_fallback() for automatic
        fallback to Porphyry. The polar threshold is approximately 90 deg - obliquity.
        Internal fallback to Porphyry is retained as a safety net.

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = _init_cardinal_cusps(asc, mc)

    # Check for polar circle - Koch undefined at |lat| >= 90 - eps
    if abs(lat) >= 90 - eps:
        return _houses_porphyry(asc, mc)

    # Precompute trigonometric values for the obliquity and latitude
    sin_obliquity = math.sin(math.radians(eps))
    tan_lat = math.tan(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))

    # Declination of the MC: sin(dec_MC) = sin(MC) * sin(ε) / cos(φ)
    # Ref: Meeus, "Astronomical Algorithms" 2nd ed., Ch. 13
    mc_dec_sin = max(
        -1.0, min(1.0, math.sin(math.radians(mc)) * sin_obliquity / cos_lat)
    )
    mc_dec_cos = math.sqrt(1.0 - mc_dec_sin * mc_dec_sin)  # always >= 0

    # Auxiliary angle c: elevation of the MC above the prime vertical
    # atan2 handles the degenerate case mc_dec_cos == 0 natively
    zenith_angle = math.degrees(math.atan2(tan_lat, mc_dec_cos))

    # Ascensional difference divided by 3 (used to offset each intermediate cusp)
    # asc_diff_third = asin(sin(zenith_angle) * sin(dec_MC)) / 3
    asc_diff_arg = max(
        -1.0, min(1.0, math.sin(math.radians(zenith_angle)) * mc_dec_sin)
    )
    asc_diff_third = math.degrees(math.asin(asc_diff_arg)) / 3.0

    # Calculate cusps using _calc_ascendant (ecliptic rising-point formula)
    cusps[11] = _calc_ascendant(armc + 30 - 2 * asc_diff_third, eps, lat, lat)
    cusps[12] = _calc_ascendant(armc + 60 - asc_diff_third, eps, lat, lat)
    cusps[2] = _calc_ascendant(armc + 120 + asc_diff_third, eps, lat, lat)
    cusps[3] = _calc_ascendant(armc + 150 + 2 * asc_diff_third, eps, lat, lat)

    # Opposite houses
    _set_opposite_cusps(cusps)

    return cusps


def _houses_regiomontanus(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Regiomontanus (Medieval rational) house system.

    Divides the celestial equator into 12 equal 30° arcs, then projects these
    divisions onto the ecliptic using great circles through the celestial poles.

    Algorithm:
        1. Divide equator into 30° segments from MC
        2. For each segment, calculate pole: tan(Pole) = tan(lat) * sin(H)
        3. Project to ecliptic using spherical trigonometry
        4. Calculate cusp longitude from pole and RAMC offset

    Mathematical Formulas:
        Hour angle offset: H = 30°, 60°, 120°, 150° for houses 11, 12, 2, 3

        Pole of projection:
            tan(P) = tan(φ) · sin(H)

        Right Ascension offset:
            R = ARMC + H - 90°

        Ecliptic longitude:
            λ = atan2(cos(R), -(sin(R)·cos(ε) + tan(P)·sin(ε)))

        Opposite houses: λ_{i+6} = (λᵢ + 180°) mod 360°

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """

    cusps = _init_cardinal_cusps(asc, mc)

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    def calc_cusp(offset_deg):
        h_rad = math.radians(offset_deg)
        tan_pole = math.tan(rad_lat) * math.sin(h_rad)

        # R = RAMC + offset - 90
        r_deg = (armc + offset_deg - 90.0) % 360.0
        r_rad = math.radians(r_deg)

        # Flip signs for East intersection (Ascendant formula)
        num = math.cos(r_rad)
        den = -(math.sin(r_rad) * math.cos(rad_eps) + tan_pole * math.sin(rad_eps))

        lon = math.degrees(math.atan2(num, den))
        return lon % 360.0

    cusps[11] = calc_cusp(30)
    cusps[12] = calc_cusp(60)
    cusps[2] = calc_cusp(120)
    cusps[3] = calc_cusp(150)

    _set_opposite_cusps(cusps)

    return cusps


def _houses_campanus(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Campanus (Prime Vertical) house system.

    Divides the prime vertical (great circle through Zenith, East, Nadir, West)
    into 12 equal 30° arcs, then projects onto the ecliptic via great circles
    through the North and South points of the horizon.

    Derivation from spherical trigonometry
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The North point of the horizon (N) is the **pole** of the prime vertical
    (PV), since it is exactly 90° from every point on the PV.

    For a division point D at angular offset *h* from the East point along the
    PV (positive towards zenith, negative towards nadir), the house circle is
    the great circle through N, D, and S.  Because N is the pole of the PV,
    the distance ND = 90° and the house circle is perpendicular to the PV at D.

    The pole of this house circle lies on the PV at 90° from D.  Converting
    the pole's horizontal coordinates (azimuth 270°, altitude 90° − h) to
    equatorial gives:

        pole_height  = arcsin(sin φ · cos h)
        H_pole       = atan2(sin h,  cos h · cos φ)

    The equator crossing of the house circle is at:

        RA_crossing  = ARMC + 90° − H_pole

    Ecliptic longitude follows from _ra_to_ecliptic_longitude().

    Ref: Smart, "Textbook on Spherical Astronomy", Ch. 3 (pole/great-circle
         geometry); Meeus, "Astronomical Algorithms", Ch. 13 (coordinate
         transforms).

    Args:
        armc: Right Ascension of the Midheaven (RAMC) in degrees.
        lat: Geographic latitude in degrees.
        eps: True obliquity of the ecliptic in degrees.
        asc: Ascendant longitude in degrees.
        mc: Midheaven longitude in degrees.

    Returns:
        List of 13 house cusp longitudes (index 0 unused).
    """
    cusps = _init_cardinal_cusps(asc, mc)

    sin_lat = math.sin(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))
    sin_eps = math.sin(math.radians(eps))
    cos_eps = math.cos(math.radians(eps))

    def _cusp_from_pv_offset(h_deg: float) -> float:
        """Ecliptic longitude for a Campanus house circle at PV offset *h_deg*."""
        h = math.radians(h_deg)
        sin_h = math.sin(h)
        cos_h = math.cos(h)

        # Pole height of the house circle (declination of the circle's pole)
        pole_height = math.degrees(math.asin(
            max(-1.0, min(1.0, sin_lat * cos_h))
        ))

        # Hour angle of the pole → RA of equator crossing
        h_pole = math.degrees(math.atan2(sin_h, cos_h * cos_lat))
        ra_crossing = armc + 90.0 - h_pole

        return _ra_to_ecliptic_longitude(ra_crossing, pole_height, sin_eps, cos_eps)

    # Cusp mapping: angular offset from East point along prime vertical
    #   +30° / +60° = above horizon (cusps 12, 11 — between ASC and MC)
    #   −30° / −60° = below horizon (cusps  2,  3 — between ASC and IC)
    cusps[12] = _cusp_from_pv_offset(30.0)
    cusps[11] = _cusp_from_pv_offset(60.0)
    cusps[2] = _cusp_from_pv_offset(-30.0)
    cusps[3] = _cusp_from_pv_offset(-60.0)

    _set_opposite_cusps(cusps)

    return cusps


def _houses_savard_a(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Savard-A house system (code 'J', John Savard's "Albategnius" houses).

    Like Campanus, the intermediate cusps are the ecliptic crossings of
    position circles through the North and South points of the horizon
    ("house circles").  Campanus anchors those circles at fixed 30°/60°
    arcs along the prime vertical; Savard-A instead anchors them where the
    prime vertical crosses the declination parallels at one third and two
    thirds of the geographic latitude.

    Derivation from spherical trigonometry
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    A point of the prime vertical (PV) at arc *x* from the East point has
    declination δ with

        sin δ = sin φ · sin x

    (right spherical triangle East-point / equator / PV).  Requiring
    δ = k·φ/3 therefore puts the anchor points at the PV arcs

        x_k = arcsin( sin(k·φ/3) / sin φ ),   k = 1, 2

    with the equator limits x_k → arcsin(k/3) as φ → 0.

    A house circle anchored at PV arc *x* is handled exactly like a
    Campanus circle (see _houses_campanus): its anchor's hour-angle offset
    from the East point and its pole height follow from the same right
    triangle,

        H_x = arctan( tan x / cos φ )          (hour-angle offset)
        p_x = arcsin( sin φ · cos x )          (pole height)

    and the cusp is the ecliptic crossing computed by _calc_ascendant at
    RA = ARMC + 90° ∓ H_x with pole height p_x (east side −, west side +).

    Ref: Smart, "Textbook on Spherical Astronomy", Ch. 3 (right spherical
    triangles, pole/great-circle geometry).

    Args:
        armc: Right Ascension of the Midheaven (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = _init_cardinal_cusps(asc, mc)

    _NEAR_ZERO = 1e-10
    sin_lat = math.sin(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))

    # PV arcs of the anchor points: declination φ/3 (inner -> cusps 12/2)
    # and 2φ/3 (outer -> cusps 11/3).
    if abs(lat) < _NEAR_ZERO:
        # Equator limit: sin(kφ/3)/sin(φ) -> k/3 as φ -> 0
        pv_arc_inner = math.degrees(math.asin(1.0 / 3.0))
        pv_arc_outer = math.degrees(math.asin(2.0 / 3.0))
    else:
        pv_arc_inner = math.degrees(
            math.asin(math.sin(math.radians(lat / 3.0)) / sin_lat)
        )
        pv_arc_outer = math.degrees(
            math.asin(math.sin(math.radians(2.0 * lat / 3.0)) / sin_lat)
        )

    # Hour-angle offsets of the anchor points from the East point
    if abs(cos_lat) < _NEAR_ZERO:
        # Observer at a pole: the prime vertical degenerates onto the
        # meridian, so both offsets collapse to it.
        if lat > 0:
            ha_outer = ha_inner = 90.0
        else:
            ha_outer = ha_inner = 270.0
    else:
        ha_outer = math.degrees(
            math.atan(math.tan(math.radians(pv_arc_outer)) / cos_lat)
        )
        ha_inner = math.degrees(
            math.atan(math.tan(math.radians(pv_arc_inner)) / cos_lat)
        )

    # Pole heights of the two house circles
    pole_outer = math.degrees(
        math.asin(sin_lat * math.cos(math.radians(pv_arc_outer)))
    )
    pole_inner = math.degrees(
        math.asin(sin_lat * math.cos(math.radians(pv_arc_inner)))
    )

    # Ecliptic crossings: east side (cusps 11/12), west side (cusps 2/3)
    cusps[11] = _calc_ascendant(armc + 90.0 - ha_outer, eps, lat, pole_outer)
    cusps[12] = _calc_ascendant(armc + 90.0 - ha_inner, eps, lat, pole_inner)
    cusps[2] = _calc_ascendant(armc + 90.0 + ha_inner, eps, lat, pole_inner)
    cusps[3] = _calc_ascendant(armc + 90.0 + ha_outer, eps, lat, pole_outer)

    _set_opposite_cusps(cusps)

    return cusps


def _houses_equal(asc: float) -> List[float]:
    """
    Equal house system (30° divisions from Ascendant).

    Simplest house system: each house is exactly 30° of ecliptic longitude.
    House 1 starts at Ascendant, each subsequent house adds 30°.

    Mathematical Formula:
        λᵢ = (λ_Asc + (i-1) × 30°) mod 360°
        where i = 1, 2, ..., 12

    Properties:
        - Works at all latitudes including polar regions
        - MC may not coincide with 10th house cusp
        - Computationally exact (no astronomical calculations)

    Args:
        asc: Ascendant longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    for i in range(1, 13):
        cusps[i] = (asc + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_whole_sign(asc: float) -> List[float]:
    """
    Whole Sign house system (ancient Hellenistic method).

    Each house occupies one complete zodiac sign. House 1 starts at 0° of the
    sign containing the Ascendant. Used extensively in ancient astrology.

    Algorithm:
        1. Find zodiac sign of Ascendant (floor(asc / 30) * 30)
        2. Each house = one complete sign (30° intervals from sign 0°)

    Mathematical Formula:
        start = floor(λ_Asc / 30°) × 30°
        λᵢ = (start + (i-1) × 30°) mod 360°
        where i = 1, 2, ..., 12

    Example:
        If Asc = 15° Taurus (45°), then:
        - House 1 starts at 30° (0° Taurus)
        - House 2 starts at 60° (0° Gemini)
        - etc.

    Properties:
        - Works at all latitudes including polar regions
        - Computationally exact (no astronomical calculations)
        - House cusps always at 0° of each sign

    Args:
        asc: Ascendant longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    # Start of sign containing Asc.  The 1e-9 deg nudge resolves the
    # boundary case: at cardinal ARMC the ascendant lands within one
    # float ULP of an exact sign cusp (e.g. 180 - 1.4e-13 at ARMC 90,
    # lat -66) and the reference assigns the NEW sign there.
    start = math.floor((asc + 1e-9) / 30.0) * 30.0
    for i in range(1, 13):
        cusps[i] = (start + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_porphyry(asc: float, mc: float) -> List[float]:
    """
    Porphyry house system (space-based trisection).

    Divides each quadrant (Asc-MC, MC-Desc, Desc-IC, IC-Asc) into three equal
    30° sections along the ecliptic. Simple and well-defined at all latitudes.

    Algorithm:
        1. Calculate arc from Asc to MC, divide by 3
        2. Calculate arc from MC to Desc (MC+180), divide by 3
        3. Opposite houses are 180° from each other

    Mathematical Formulas:
        Quadrant MC→Asc (houses 11, 12):
            step = (λ_Asc - λ_MC) mod 360° / 3
            λ₁₁ = λ_MC + step
            λ₁₂ = λ_MC + 2·step

        Quadrant Asc→IC (houses 2, 3):
            step = (λ_IC - λ_Asc) mod 360° / 3
            λ₂ = λ_Asc + step
            λ₃ = λ_Asc + 2·step

        Opposite houses:
            λ₅ = (λ₁₁ + 180°) mod 360°
            λ₆ = (λ₁₂ + 180°) mod 360°
            λ₈ = (λ₂ + 180°) mod 360°
            λ₉ = (λ₃ + 180°) mod 360°

    Properties:
        - Works at all latitudes including polar regions
        - Used as fallback when Placidus/Koch fail
        - Computationally simple (no trigonometry)

    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = _init_cardinal_cusps(asc, mc)

    # Trisect the ecliptic arc between angles
    # Arc 10-1
    diff = (asc - mc) % 360.0
    step = diff / 3.0
    cusps[11] = (mc + step) % 360.0
    cusps[12] = (mc + 2 * step) % 360.0

    # Arc 1-4
    ic = cusps[4]
    diff = (ic - asc) % 360.0
    step = diff / 3.0
    cusps[2] = (asc + step) % 360.0
    cusps[3] = (asc + 2 * step) % 360.0

    # Opposites
    _set_opposite_cusps(cusps)

    return cusps


def _houses_sripati(asc: float, mc: float) -> List[float]:
    """
    Sripati house system (traditional Indian method).

    Creates house cusps at the midpoints of Porphyry house cusps. Each Sripati
    cusp is the midpoint between the previous and current Porphyry cusp.

    Algorithm:
        1. Calculate Porphyry house cusps
        2. For each house i, the Sripati cusp is:
           H'[i] = (H[i-1] + H[i]) / 2  (midpoint on the ecliptic circle)

    Mathematical Formula:
        λ'ᵢ = λᵢ₋₁ + ((λᵢ - λᵢ₋₁) mod 360°) / 2

    This effectively shifts the house boundaries so that the "core" of each
    Porphyry house becomes the cusp, making the Porphyry cusp the center of
    each Sripati house.

    Properties:
        - Works at all latitudes (inherits Porphyry's polar stability)
        - Used in traditional Indian astrology (Jyotish)
        - Computationally simple (no trigonometry beyond Porphyry)

    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are house cusps)
    """
    # Get Porphyry cusps as base
    porphyry = _houses_porphyry(asc, mc)

    cusps = [0.0] * 13
    for i in range(1, 13):
        prev_i = 12 if i == 1 else i - 1
        # Calculate midpoint between previous and current Porphyry cusp
        # Handle wrap-around at 360 degrees
        diff = (porphyry[i] - porphyry[prev_i]) % 360.0
        cusps[i] = (porphyry[prev_i] + diff / 2.0) % 360.0

    return cusps


def _houses_pullen_sd(asc: float, mc: float) -> List[float]:
    """
    Pullen SD (Sinusoidal Delta) house system, also known as Neo-Porphyry.

    Invented by Walter Pullen in 1994. Like Porphyry, based on ecliptic quadrant
    divisions, but fits house widths to a sine wave pattern rather than equal
    trisection.

    Algorithm:
        - Ideal house size = 30°
        - For each quadrant, compute deviation d = quadrant_size - 90°
        - Middle house of quadrant (2nd, 5th, 8th, 11th) gets d/2 added to 30°
        - Side houses each get d/4 added to 30°
        - If middle house size would be negative, set to 0

    House width pattern for quadrant:
        x+n, x, x+n  where x = 30 and n = d/4
        Middle house = x + 2n = 30 + d/2

    Properties:
        - Works at all latitudes (based on ecliptic only)
        - Houses in larger quadrants are larger
        - Middle house absorbs most of the size variation

    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are house cusps)
    """
    cusps = _init_cardinal_cusps(asc, mc)

    # Calculate quadrant sizes
    # Quadrant 1: MC to Asc (houses 11, 12)
    q1 = (asc - mc) % 360.0
    # Quadrant 2: Asc to IC (houses 2, 3)
    ic = cusps[4]
    q2 = (ic - asc) % 360.0
    # Quadrant 3: IC to Desc (houses 5, 6)
    desc = cusps[7]
    q3 = (desc - ic) % 360.0
    # Quadrant 4: Desc to MC (houses 8, 9)
    q4 = (mc - desc) % 360.0

    def calc_house_sizes(quadrant_size: float) -> tuple:
        """Calculate house sizes for a quadrant using sinusoidal delta."""
        d = quadrant_size - 90.0  # deviation from ideal 90°
        # Middle house gets d/2 added to 30°
        middle = 30.0 + d / 2.0
        # If middle would be negative, set to 0 and redistribute
        if middle < 0:
            middle = 0.0
            # Remaining space divided between side houses
            side = quadrant_size / 2.0
        else:
            # Side houses each get d/4 added to 30°
            side = 30.0 + d / 4.0
        return (side, middle, side)

    # Quadrant 1: MC to Asc (houses 11, 12, then Asc)
    sizes1 = calc_house_sizes(q1)
    cusps[11] = (mc + sizes1[0]) % 360.0
    cusps[12] = (cusps[11] + sizes1[1]) % 360.0

    # Quadrant 2: Asc to IC (houses 2, 3, then IC)
    sizes2 = calc_house_sizes(q2)
    cusps[2] = (asc + sizes2[0]) % 360.0
    cusps[3] = (cusps[2] + sizes2[1]) % 360.0

    # Quadrant 3: IC to Desc (houses 5, 6, then Desc)
    sizes3 = calc_house_sizes(q3)
    cusps[5] = (ic + sizes3[0]) % 360.0
    cusps[6] = (cusps[5] + sizes3[1]) % 360.0

    # Quadrant 4: Desc to MC (houses 8, 9, then MC)
    sizes4 = calc_house_sizes(q4)
    cusps[8] = (desc + sizes4[0]) % 360.0
    cusps[9] = (cusps[8] + sizes4[1]) % 360.0

    return cusps


def _houses_pullen_sr(asc: float, mc: float) -> List[float]:
    """
    Pullen SR (Sinusoidal Ratio) house system.

    Proposed by Walter Pullen in 2016 as an improvement over Pullen SD. Uses
    ratio multipliers instead of additive offsets, avoiding negative house
    sizes for small quadrants.

    Algorithm:
        For quadrant size q:
        - Small quadrant houses: rx, x, rx  (sizes sum to q)
        - Large quadrant houses: r³x, r⁴x, r³x  (sizes sum to 180-q)

        Solving: rx + x + rx = q  →  x(2r + 1) = q
                 r³x + r⁴x + r³x = 180 - q  →  x·r³(2 + r) = 180 - q

        Dividing: r³(2 + r) / (2r + 1) = (180 - q) / q
        This is solved numerically for r.

    Properties:
        - Works at all latitudes
        - No negative house sizes even for extreme quadrants
        - More elegant mathematical formulation

    Args:
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are house cusps)
    """
    cusps = _init_cardinal_cusps(asc, mc)

    # Calculate quadrant sizes
    q1 = (asc - mc) % 360.0  # MC to Asc
    ic = cusps[4]
    q2 = (ic - asc) % 360.0  # Asc to IC
    desc = cusps[7]
    q3 = (desc - ic) % 360.0  # IC to Desc
    q4 = (mc - desc) % 360.0  # Desc to MC

    def calc_house_sizes_ratio(quadrant_size: float) -> tuple:
        """Calculate house sizes using sinusoidal ratio method."""
        q = quadrant_size
        q_opp = 180.0 - q  # opposite quadrant

        if abs(q - 90.0) < 0.0001:
            # Quadrant is exactly 90°, r = 1, equal houses of 30°
            return (30.0, 30.0, 30.0)

        # Target ratio: r³(2 + r) / (2r + 1) = q_opp / q
        target = q_opp / q if q > 0.0001 else 1000.0

        # Solve for r using Newton-Raphson
        # f(r) = r³(2 + r) / (2r + 1) - target = 0
        r = 1.0  # initial guess
        for _ in range(20):  # max iterations
            numerator = r * r * r * (2.0 + r)
            denominator = 2.0 * r + 1.0
            f = numerator / denominator - target

            # Derivative: d/dr[r³(2+r)/(2r+1)]
            # Using quotient rule
            dn = 3.0 * r * r * (2.0 + r) + r * r * r  # derivative of numerator
            dd = 2.0  # derivative of denominator
            df = (dn * denominator - numerator * dd) / (denominator * denominator)

            if abs(df) < 1e-12:
                break
            r_new = r - f / df
            if r_new < 0.01:
                r_new = 0.01  # keep r positive
            if abs(r_new - r) < 1e-10:
                break
            r = r_new

        # Calculate x from: x(2r + 1) = q
        x = q / (2.0 * r + 1.0) if (2.0 * r + 1.0) > 0.0001 else q / 3.0

        # House sizes: rx, x, rx
        side = r * x
        middle = x

        return (side, middle, side)

    # Quadrant 1: MC to Asc (houses 11, 12)
    sizes1 = calc_house_sizes_ratio(q1)
    cusps[11] = (mc + sizes1[0]) % 360.0
    cusps[12] = (cusps[11] + sizes1[1]) % 360.0

    # Quadrant 2: Asc to IC (houses 2, 3)
    sizes2 = calc_house_sizes_ratio(q2)
    cusps[2] = (asc + sizes2[0]) % 360.0
    cusps[3] = (cusps[2] + sizes2[1]) % 360.0

    # Quadrant 3: IC to Desc (houses 5, 6)
    sizes3 = calc_house_sizes_ratio(q3)
    cusps[5] = (ic + sizes3[0]) % 360.0
    cusps[6] = (cusps[5] + sizes3[1]) % 360.0

    # Quadrant 4: Desc to MC (houses 8, 9)
    sizes4 = calc_house_sizes_ratio(q4)
    cusps[8] = (desc + sizes4[0]) % 360.0
    cusps[9] = (cusps[8] + sizes4[1]) % 360.0

    return cusps


def _houses_alcabitius(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Alcabitius (Alchabitius) house system (ancient Arabic method).

    Medieval Arabic system that divides the diurnal and nocturnal arcs differently
    than Placidus, using a simpler geometric approach.

    Algorithm:
        1. Calculate RA of Ascendant and MC
        2. Divide RA intervals between angles
        3. Convert RA divisions back to ecliptic longitude

    Mathematical Formulas:
        Right Ascension of Ascendant:
            RA_Asc = atan2(sin(λ_Asc)·cos(ε), cos(λ_Asc))

        Quadrant MC→Asc (houses 11, 12):
            arc = (RA_Asc - ARMC) mod 360°
            step = arc / 3
            RA₁₁ = ARMC + step
            RA₁₂ = ARMC + 2·step

        Quadrant Asc→IC (houses 2, 3):
            RA_IC = (ARMC + 180°) mod 360°
            arc = (RA_IC - RA_Asc) mod 360°
            step = arc / 3
            RA₂ = RA_Asc + step
            RA₃ = RA_Asc + 2·step

        RA to ecliptic conversion:
            λ = atan2(sin(RA), cos(RA)·cos(ε))

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    # Alcabitius
    # Time trisection of Ascendant's diurnal arc, projected by Hour Circles.
    # RA_11 = RAMC + SA/3.
    # SA = 90 + AD_asc.

    cusps = _init_cardinal_cusps(asc, mc)

    math.radians(lat)
    rad_eps = math.radians(eps)

    # RA of Ascendant
    y = math.sin(math.radians(asc)) * math.cos(rad_eps)
    x = math.cos(math.radians(asc))
    ra_asc = math.degrees(math.atan2(y, x)) % 360.0

    # Arc from MC to Asc
    arc = (ra_asc - armc) % 360.0
    step = arc / 3.0

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[11] = get_lon_from_ra(armc + step)
    cusps[12] = get_lon_from_ra(armc + 2 * step)

    # Sector 2: Asc to IC
    ra_ic = (armc + 180.0) % 360.0
    arc2 = (ra_ic - ra_asc) % 360.0
    step2 = arc2 / 3.0

    cusps[2] = get_lon_from_ra(ra_asc + step2)
    cusps[3] = get_lon_from_ra(ra_asc + 2 * step2)

    _set_opposite_cusps(cusps)

    return cusps


def _houses_polich_page(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Polich-Page (Topocentric) house system.

    Developed in 1960s to account for observer's actual position on Earth's surface
    rather than at Earth's center. Uses modified pole calculations.

    Algorithm:
        Similar to Regiomontanus but with modified pole factors that account
        for the topocentric (observer-centered) perspective.

    Mathematical Formulas:
        Pole factors:
            For houses 11, 3: factor = 1/3
            For houses 12, 2: factor = 2/3

        Pole calculation (modified from Regiomontanus):
            tan(P) = tan(φ) × factor

        Hour angle offsets:
            H = 30°, 60°, 120°, 150° for houses 11, 12, 2, 3

        Right Ascension offset:
            R = ARMC + H - 90°

        Ecliptic projection:
            λ = atan2(cos(R), -(sin(R)·cos(ε) + tan(P)·sin(ε)))

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    # Polich/Page (Topocentric)
    # Uses Pole method with tan(Pole) = tan(lat) * factor.
    # factor = 1/3 for 11/3, 2/3 for 12/2.

    cusps = _init_cardinal_cusps(asc, mc)

    rad_lat = math.radians(lat)
    rad_eps = math.radians(eps)

    def calc_cusp(offset_deg, factor):
        tan_pole = math.tan(rad_lat) * factor

        # R = RAMC + offset - 90
        r_deg = (armc + offset_deg - 90.0) % 360.0
        r_rad = math.radians(r_deg)

        num = math.cos(r_rad)
        den = -(math.sin(r_rad) * math.cos(rad_eps) + tan_pole * math.sin(rad_eps))

        lon = math.degrees(math.atan2(num, den))
        return lon % 360.0

    cusps[11] = calc_cusp(30, 1.0 / 3.0)
    cusps[12] = calc_cusp(60, 2.0 / 3.0)
    cusps[2] = calc_cusp(120, 2.0 / 3.0)
    cusps[3] = calc_cusp(150, 1.0 / 3.0)

    _set_opposite_cusps(cusps)

    return cusps


def _houses_morinus(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Morinus house system (equatorial divisions).

    Divides the celestial equator into 12 equal 30° sections starting from 0° Aries,
    then projects to ecliptic. Independent of observer location.

    Algorithm:
        1. Divide equator into 30° RA sections from 0h RA
        2. Convert each RA to ecliptic longitude using obliquity

    Mathematical Formulas:
        For each house at RA = ARMC + (i-10)×30° where i = 10, 11, 12, 1, 2, 3:

        Projection via ecliptic poles:
            tan(λ) = tan(RA) × cos(ε)

        Using atan2 for proper quadrant:
            λ = atan2(sin(RA)·cos(ε), cos(RA))

    Properties:
        - Location-independent (ignores latitude)
        - Morinus 1st cusp differs from standard Ascendant
        - Morinus 10th cusp differs from standard MC

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    # Morinus
    # Projects Equator points (RAMC + 30) to Ecliptic via Ecliptic Poles.
    # Lon = Lon of point on Equator.
    # tan(lon) = tan(ra) * cos(eps).

    cusps = [0.0] * 13
    # Morinus Ascendant? Standard houses returns standard Asc.
    # But cusps[1] should be Morinus Ascendant (RAMC+90 projected).
    # Let's follow standard behavior: cusps array contains system cusps.

    rad_eps = math.radians(eps)

    def get_lon(ra):
        # tan(lon) = tan(ra) * cos(eps)
        y = math.sin(math.radians(ra)) * math.cos(rad_eps)
        x = math.cos(math.radians(ra))
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[10] = get_lon(armc)
    cusps[11] = get_lon(armc + 30)
    cusps[12] = get_lon(armc + 60)
    cusps[1] = get_lon(armc + 90)
    cusps[2] = get_lon(armc + 120)
    cusps[3] = get_lon(armc + 150)

    cusps[4] = (cusps[10] + 180) % 360.0
    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[7] = (cusps[1] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_meridian(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Meridian (Zariel/Axial Rotation) house system.

    Based on meridian passages, divides RA from MC in equal 30° intervals.
    Related to Morinus but starts from MC instead of 0° Aries.

    Algorithm:
        Projects equator points (ARMC + n×30°) to ecliptic via celestial poles.

    Mathematical Formulas:
        For each house at RA = ARMC + (i-10)×30° where i = 10, 11, 12, 1, 2, 3:

        Projection via celestial poles:
            tan(λ) = tan(RA) / cos(ε)

        Using atan2 for proper quadrant:
            λ = atan2(sin(RA), cos(RA)·cos(ε))

    Difference from Morinus:
        Morinus:  λ = atan2(sin(RA)·cos(ε), cos(RA))
        Meridian: λ = atan2(sin(RA), cos(RA)·cos(ε))

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    # Meridian (Axial)
    # Projects Equator points to Ecliptic via Celestial Poles.
    # tan(lon) = tan(ra) / cos(eps).

    cusps = [0.0] * 13

    rad_eps = math.radians(eps)

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    cusps[10] = get_lon_from_ra(armc)
    cusps[11] = get_lon_from_ra(armc + 30)
    cusps[12] = get_lon_from_ra(armc + 60)
    cusps[1] = get_lon_from_ra(armc + 90)
    cusps[2] = get_lon_from_ra(armc + 120)
    cusps[3] = get_lon_from_ra(armc + 150)

    cusps[4] = (cusps[10] + 180) % 360.0
    cusps[5] = (cusps[11] + 180) % 360.0
    cusps[6] = (cusps[12] + 180) % 360.0
    cusps[7] = (cusps[1] + 180) % 360.0
    cusps[8] = (cusps[2] + 180) % 360.0
    cusps[9] = (cusps[3] + 180) % 360.0

    return cusps


def _houses_vehlow(asc: float) -> List[float]:
    """
    Vehlow house system (Equal with Asc in middle of House 1).

    Variant of Equal houses where the Ascendant falls at 15° into House 1
    rather than at the cusp. Each house is still 30°.

    Mathematical Formula:
        start = (λ_Asc - 15°) mod 360°
        λᵢ = (start + (i-1) × 30°) mod 360°
        where i = 1, 2, ..., 12

    Properties:
        - Ascendant falls exactly at 15° into House 1
        - Works at all latitudes including polar regions
        - Computationally exact (no astronomical calculations)

    Args:
        asc: Ascendant longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    start = (asc - 15.0) % 360.0
    for i in range(1, 13):
        cusps[i] = (start + (i - 1) * 30.0) % 360.0
    return cusps


def _houses_carter(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Carter Poli-Equatorial house system.

    Equal 30° divisions on the celestial equator starting from RA of Ascendant,
    projected to ecliptic via hour circles.

    Algorithm:
        1. Calculate RA of Ascendant
        2. Add 30° RA increments for each house
        3. Convert each RA to ecliptic longitude

    Mathematical Formulas:
        Right Ascension of Ascendant:
            RA_Asc = atan2(sin(λ_Asc)·cos(ε), cos(λ_Asc))

        For each house i (1-12):
            RA = RA_Asc + (i-1) × 30°

        RA to ecliptic conversion via celestial poles:
            λ = atan2(sin(RA), cos(RA)·cos(ε))

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """

    cusps = [0.0] * 13

    rad_eps = math.radians(eps)

    # Get RA of Ascendant
    y = math.sin(math.radians(asc)) * math.cos(rad_eps)
    x = math.cos(math.radians(asc))
    ra_asc = math.degrees(math.atan2(y, x)) % 360.0

    def get_lon_from_ra(ra):
        y = math.sin(math.radians(ra))
        x = math.cos(math.radians(ra)) * math.cos(rad_eps)
        lon = math.degrees(math.atan2(y, x))
        return lon % 360.0

    # Equal 30-degree divisions from RA of Asc
    for i in range(1, 13):
        ra = (ra_asc + (i - 1) * 30.0) % 360.0
        cusps[i] = get_lon_from_ra(ra)

    return cusps


def _gauquelin_cusp_for_sector(
    sector_offset: int,
    base_ramc: float,
    ramc_sign: float,
    ascensional_diff: float,
    tan_obliquity: float,
    sin_obliquity: float,
    tan_lat: float,
    eps: float,
    lat: float,
    niter_max: int,
    convergence_threshold: float,
    near_zero: float,
) -> float:
    """
    Compute one intermediate Gauquelin sector cusp via iterative pole-height refinement.

    The algorithm projects a horizon-fraction (sector_offset/9 of the ascensional
    difference arc) onto the ecliptic by iterating the pole height of the house
    circle until the cusp longitude converges.

    This is the spherical-trigonometry core shared by both quadrant loops of the
    36-sector system (Ref: Gauquelin, "L'influence des astres", 1955; geometric
    derivation follows Meeus, "Astronomical Algorithms", Ch. 13).

    Args:
        sector_offset: Integer fraction index (1-8) of the quadrant arc.
        base_ramc: Base RAMC value for the quadrant (armc or armc+180).
        ramc_sign: +1.0 or -1.0 controlling the direction of the RAMC offset.
        ascensional_diff: Ascensional difference angle (degrees) for the latitude.
        tan_obliquity: tan(ε) — tangent of the ecliptic obliquity.
        sin_obliquity: sin(ε) — sine of the ecliptic obliquity.
        tan_lat: tan(φ) — tangent of the geographic latitude.
        eps: True obliquity of ecliptic (degrees).
        lat: Geographic latitude (degrees).
        niter_max: Maximum number of refinement iterations.
        convergence_threshold: Angular convergence tolerance (degrees).
        near_zero: Threshold for treating a tangent as effectively zero.

    Returns:
        Ecliptic longitude of the sector cusp (degrees, 0-360).
    """
    # Initial pole height: latitude of the great circle pole for this sector fraction
    arc_fraction_deg = ascensional_diff * sector_offset / 9.0
    initial_pole_height = math.degrees(
        math.atan(math.sin(math.radians(arc_fraction_deg)) / tan_obliquity)
    )

    # Intermediate RAMC for this sector
    intermediate_ramc = (base_ramc + ramc_sign * (90.0 / 9.0) * sector_offset) % 360.0

    # Seed estimate: ecliptic point rising at the intermediate RAMC with the initial pole
    cusp = _calc_ascendant(intermediate_ramc, eps, lat, initial_pole_height)

    # Declination of the seed cusp point
    cusp_dec = math.asin(sin_obliquity * math.sin(math.radians(cusp)))
    tan_cusp_dec = math.tan(cusp_dec)

    # Degenerate case: cusp on the equator, no pole-height correction possible
    if abs(tan_cusp_dec) < near_zero:
        return intermediate_ramc

    # First pole-height refinement from seed
    pole_arg = max(-1.0, min(1.0, tan_lat * tan_cusp_dec))
    pole_height = math.degrees(
        math.atan(
            math.sin(
                math.radians(math.degrees(math.asin(pole_arg)) * sector_offset / 9.0)
            )
            / tan_cusp_dec
        )
    )
    cusp = _calc_ascendant(intermediate_ramc, eps, lat, pole_height)

    # Iterative refinement until convergence
    prev_cusp = 0.0
    for iteration in range(1, niter_max + 1):
        updated_dec = math.asin(sin_obliquity * math.sin(math.radians(cusp)))
        tan_cusp_dec_updated = math.tan(updated_dec)

        if abs(tan_cusp_dec_updated) < near_zero:
            cusp = intermediate_ramc
            break

        pole_arg = max(-1.0, min(1.0, tan_lat * tan_cusp_dec_updated))
        pole_height = math.degrees(
            math.atan(
                math.sin(
                    math.radians(
                        math.degrees(math.asin(pole_arg)) * sector_offset / 9.0
                    )
                )
                / tan_cusp_dec_updated
            )
        )
        cusp = _calc_ascendant(intermediate_ramc, eps, lat, pole_height)

        if iteration > 1:
            angular_diff = abs(cusp - prev_cusp)
            if angular_diff > 180.0:
                angular_diff = 360.0 - angular_diff
            if angular_diff < convergence_threshold:
                break
        prev_cusp = cusp

    return cusp


def _houses_gauquelin(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Gauquelin 36-Sector house system.

    Divides the celestial sphere into 36 sectors of 10° each using a
    Placidus-inspired construction: the diurnal and nocturnal semi-arcs are
    divided into ninths, and the resulting great circles are projected onto
    the ecliptic via iterative pole-height refinement.

    Reference: Michel Gauquelin, "L'influence des astres" (1955).
    Mathematical derivation: Meeus, "Astronomical Algorithms" 2nd ed., Ch. 13.

    Algorithm:
        1. Compute the ascensional difference:
               a = arcsin(tan(φ) · tan(ε))
           where φ = geographic latitude, ε = ecliptic obliquity.
        2. Divide each quadrant arc into ninths.
        3. For each of the 8 intermediate sectors per quadrant, call
           _gauquelin_cusp_for_sector() to find the ecliptic intersection via
           iterative pole-height refinement.
        4. Opposite sectors are set by adding 180°.
        5. Cardinals: sector 1 = Asc, 10 = MC, 19 = Desc, 28 = IC.

    Polar Limitation:
        Fails within polar circle (|φ| >= 90° - ε); falls back to equal
        10°-divisions from the Ascendant.

    Args:
        armc: Right Ascension of the Midheaven (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 37 sector longitudes (index 0 unused, sectors 1-36)
    """
    NEAR_ZERO = 1e-10
    CONVERGENCE_THRESHOLD = 1e-8
    niter_max = 100

    cusps = [0.0] * 37  # index 0 unused; sectors 1-36

    # Polar circle fallback: equal 10° divisions from the Ascendant
    if abs(lat) >= 90.0 - eps:
        for i in range(1, 37):
            cusps[i] = (asc - (i - 1) * 10.0) % 360.0
        return cusps

    # Precompute obliquity and latitude trig quantities
    sin_obliquity = math.sin(math.radians(eps))
    tan_obliquity = math.tan(math.radians(eps))
    tan_lat = math.tan(math.radians(lat))

    # Ascensional difference: a = arcsin(tan(φ) · tan(ε))
    asc_diff_arg = max(-1.0, min(1.0, tan_lat * tan_obliquity))
    ascensional_diff = math.degrees(math.asin(asc_diff_arg))

    # Shared keyword arguments for the cusp helper
    _cusp_kwargs: dict[str, Any] = dict(
        ascensional_diff=ascensional_diff,
        tan_obliquity=tan_obliquity,
        sin_obliquity=sin_obliquity,
        tan_lat=tan_lat,
        eps=eps,
        lat=lat,
        niter_max=niter_max,
        convergence_threshold=CONVERGENCE_THRESHOLD,
        near_zero=NEAR_ZERO,
    )

    # --- Fourth/Second quadrant: sectors 2-9 (RAMC offsets +10° … +80°) ---
    for ih in range(2, 10):
        sector_offset = 10 - ih  # counts down 8..1 as ih goes 2..9
        cusps[ih] = _gauquelin_cusp_for_sector(
            sector_offset=sector_offset,
            base_ramc=armc,
            ramc_sign=+1.0,
            **_cusp_kwargs,
        )
        cusps[ih + 18] = (cusps[ih] + 180.0) % 360.0  # opposite sector

    # --- First/Third quadrant: sectors 29-36 (RAMC offsets -10° … -80°) ---
    for ih in range(29, 37):
        sector_offset = ih - 28  # counts up 1..8 as ih goes 29..36
        cusps[ih] = _gauquelin_cusp_for_sector(
            sector_offset=sector_offset,
            base_ramc=(armc + 180.0) % 360.0,
            ramc_sign=-1.0,
            **_cusp_kwargs,
        )
        cusps[ih - 18] = (cusps[ih] + 180.0) % 360.0  # opposite sector

    # Cardinal sectors
    cusps[1] = asc
    cusps[10] = mc
    cusps[19] = (asc + 180.0) % 360.0
    cusps[28] = (mc + 180.0) % 360.0

    return cusps


def _ecliptic_to_ra_simple(lon: float, eps: float) -> float:
    """Convert ecliptic longitude to right ascension."""
    rad_lon = math.radians(lon)
    rad_eps = math.radians(eps)

    y = math.sin(rad_lon) * math.cos(rad_eps)
    x = math.cos(rad_lon)
    ra = math.degrees(math.atan2(y, x))

    return ra % 360.0


def _ra_to_ecliptic_simple(ra: float, eps: float) -> float:
    """Convert right ascension to ecliptic longitude."""
    rad_ra = math.radians(ra)
    rad_eps = math.radians(eps)

    y = math.sin(rad_ra)
    x = math.cos(rad_ra) * math.cos(rad_eps)
    lon = math.degrees(math.atan2(y, x))

    return lon % 360.0


def _rotate_spherical_x_axis(x: List[float], eps: float) -> List[float]:
    """
    Standard rotation of spherical coordinates around the x-axis.

    Rotates spherical coordinates [lon, lat, r] by angle eps.
    This is the standard x-axis rotation matrix applied to spherical
    coordinates, used for ecliptic/equatorial conversion and similar
    coordinate frame rotations.

    Reference: Montenbruck & Pfleger, 'Astronomy on the Personal Computer',
    par. 1.3; Meeus, 'Astronomical Algorithms' 2nd ed., ch. 13.

    Args:
        x: [longitude, latitude, radius] in degrees
        eps: Rotation angle in degrees

    Returns:
        Transformed [longitude, latitude, radius]
    """
    eps_rad = math.radians(eps)
    lon_rad = math.radians(x[0])
    lat_rad = math.radians(x[1])
    r = x[2]

    # Convert to Cartesian
    cos_lat = math.cos(lat_rad)
    x_cart = r * cos_lat * math.cos(lon_rad)
    y_cart = r * cos_lat * math.sin(lon_rad)
    z_cart = r * math.sin(lat_rad)

    # Rotate around x-axis by eps
    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)

    x_new = x_cart
    y_new = y_cart * cos_eps + z_cart * sin_eps
    z_new = -y_cart * sin_eps + z_cart * cos_eps

    # Convert back to spherical
    r_new = math.sqrt(x_new**2 + y_new**2 + z_new**2)

    if r_new > 0:
        lon_new = math.degrees(math.atan2(y_new, x_new))
        lat_new = math.degrees(math.asin(z_new / r_new))
    else:
        lon_new = 0.0
        lat_new = 0.0

    return [(lon_new % 360.0), lat_new, r_new]


def _houses_krusinski(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Krusinski-Pisa house system (code 'U').

    The Krusinski-Pisa system uses the great circle passing through the
    Ascendant and Zenith as the fundamental plane. This circle is divided
    into 12 equal 30 deg arcs, and each division point is projected back
    onto the ecliptic to define the house cusps.

    Algorithm derived from the geometric definition by:
    Bogdan Krusinski, presentation at XIV Szkola Astrologii Humanistycznej (1995);
    Milan Pisa, Konstelace No. 22, Czech Astrological Society (1997);
    Goelzer, 'Der Ich-Kosmos', Goetheanum (1995).

    Coordinate transformation sequence (standard spherical rotations):
        Forward transform (ecliptic -> house circle):
            1. Ecliptic -> equatorial (x-axis rotation by -eps)
            2. Equatorial -> hour angle frame (RA rotation by -(ARMC - 90))
            3. Hour angle -> horizontal (x-axis rotation by -(90 - lat))
            4. Horizontal -> house circle (x-axis rotation by -90)

        Inverse transform (house circle -> ecliptic):
            Reverse of the above sequence for each 30 deg division point.

    Reference: Montenbruck & Pfleger, 'Astronomy on the Personal Computer',
    par. 1.3 (coordinate rotation matrices).

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are cusps)
    """
    cusps = [0.0] * 13

    # Within the polar circle the ecliptic can intersect the horizon
    # "backwards", leaving the computed Ascendant west of the MC; swap to the
    # opposite point so the Asc-Zenith circle is anchored at the eastern
    # horizon intersection.
    if difdeg2n(asc, mc) < 0:
        asc = (asc + 180.0) % 360.0

    # --- Forward transform: locate the Ascendant in the house-circle frame ---
    # Spherical triple [longitude, latitude, radius] of the Ascendant
    asc_vec = [asc, 0.0, 1.0]

    # Ecliptic -> equatorial (rotate by obliquity)
    asc_vec = _rotate_spherical_x_axis(asc_vec, -eps)

    # Equatorial -> hour angle frame (subtract local sidereal time offset)
    asc_vec[0] = (asc_vec[0] - (armc - 90.0)) % 360.0

    # Hour angle -> horizontal (rotate by co-latitude)
    asc_vec = _rotate_spherical_x_axis(asc_vec, -(90.0 - lat))

    # Save horizon longitude of Asc to restore during inverse transform
    asc_horizon_lon = asc_vec[0]

    # Zero the horizon longitude (align Asc with reference meridian)
    asc_vec[0] = 0.0

    # Horizontal -> house circle (rotate to Asc-Zenith great circle plane)
    asc_vec = _rotate_spherical_x_axis(asc_vec, -90.0)

    # --- Inverse transform: house circle -> ecliptic ---
    # Divide the great circle into 12 equal arcs of 30 deg each
    sin_obliquity = math.sin(math.radians(eps))
    cos_obliquity = math.cos(math.radians(eps))

    for i in range(6):
        # Division point on the house circle (0, 30, 60, 90, 120, 150 deg)
        point = [30.0 * i, 0.0, 1.0]

        # House circle -> horizontal (reverse the final forward rotation)
        point = _rotate_spherical_x_axis(point, 90.0)

        # Restore horizon longitude
        point[0] = (point[0] + asc_horizon_lon) % 360.0

        # Horizontal -> equatorial (rotate by co-latitude, reverse direction)
        point = _rotate_spherical_x_axis(point, 90.0 - lat)

        # Restore RA from hour angle frame
        point[0] = (point[0] + (armc - 90.0)) % 360.0

        # Equatorial -> ecliptic (RA to ecliptic longitude, zero pole height)
        lon = _ra_to_ecliptic_longitude(point[0], 0.0, sin_obliquity, cos_obliquity)

        cusps[i + 1] = lon
        cusps[i + 7] = (lon + 180.0) % 360.0

    return cusps


def _houses_equal_mc(asc: float, mc: float) -> List[float]:
    """
    Equal houses from MC (code 'D').

    All houses are exactly 30° each, with the MC as the 10th house cusp.
    Unlike Equal from Ascendant, the Ascendant is NOT a house cusp in this system.

    Algorithm:
        H10 = MC
        H11 = MC + 30°
        H12 = MC + 60°
        H1 = MC + 90° (= MC - 270°)
        ... and so on

    Properties:
        - Works at all latitudes
        - MC is exactly on the 10th house cusp
        - Ascendant does NOT coincide with the 1st house cusp (stored in ascmc)

    Args:
        asc: Ascendant longitude in degrees (unused, kept for API compatibility)
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    # H10 = MC, then each house is 30° apart
    # H1 = MC - 270° = MC + 90°
    for i in range(1, 13):
        # House 10 is at MC, so offset = (i - 10) * 30
        cusps[i] = (mc + (i - 10) * 30.0) % 360.0
    return cusps


def _sunshine_arc_to_ecliptic(
    arc_offset_ra: float,
    is_diurnal: bool,
    armc: float,
    lat: float,
    sun_dec: float,
    sin_lat: float,
    cos_lat: float,
    cos_sun_dec: float,
    tan_sun_dec: float,
    sin_obliquity: float,
    cos_obliquity: float,
) -> float:
    """
    Map one trisected semi-arc offset to an ecliptic longitude.

    The Sunshine system divides the Sun's diurnal and nocturnal semi-arcs into
    thirds.  For each trisection point, this function uses spherical triangle
    geometry to find the corresponding ecliptic cusp (Makransky 1988,
    "Primary Directions").

    Steps:
        1. Convert the RA-arc offset to a great-circle arc length (xhs).
        2. Solve the spherical triangle {pole, arc midpoint, horizon point}
           via the spherical law of cosines to find side c.
        3. Apply the spherical law of sines to get the zenith distance (zd).
        4. Derive the equatorial RA intersection and the pole height of the
           house circle.
        5. Project to the ecliptic via the rising-point formula.

    Args:
        arc_offset_ra: RA offset along the semi-arc (degrees).
        is_diurnal: True for diurnal semi-arc cusps (houses 8-12 side),
                    False for nocturnal semi-arc cusps (houses 2-6 side).
        armc: Right Ascension of the Midheaven (degrees).
        lat: Geographic latitude (degrees).
        sun_dec: Sun's declination (degrees).
        sin_lat: sin(φ) — pre-computed for speed.
        cos_lat: cos(φ) — pre-computed for speed.
        cos_sun_dec: cos(δ☉) — pre-computed for speed.
        tan_sun_dec: tan(δ☉) — pre-computed for speed.
        sin_obliquity: sin(ε) — pre-computed for speed.
        cos_obliquity: cos(ε) — pre-computed for speed.

    Returns:
        Ecliptic longitude of the cusp (degrees, 0-360).
    """
    # Step 1: great-circle arc subtended by the RA offset on the Sun's parallel
    # xhs = 2 · arcsin(cos(δ☉) · sin(arc_offset / 2))
    great_circle_arc = 2.0 * math.degrees(
        math.asin(cos_sun_dec * math.sin(math.radians(arc_offset_ra / 2.0)))
    )

    # Step 2a: interior angle α at the pole of the spherical triangle
    # cos(α) = tan(δ☉) · tan(great_circle_arc / 2)
    cos_alpha = max(
        -1.0, min(1.0, tan_sun_dec * math.tan(math.radians(great_circle_arc / 2.0)))
    )
    alpha = math.degrees(math.acos(cos_alpha))

    # Step 2b: triangle side and opening angle differ for diurnal vs nocturnal arc
    if is_diurnal:
        # Diurnal arc: triangle opening angle is supplement of α
        triangle_angle = 180.0 - alpha
        triangle_side = 90.0 - lat + sun_dec
    else:
        # Nocturnal arc: triangle opening angle equals α
        triangle_angle = alpha
        triangle_side = 90.0 - lat - sun_dec

    # Step 2c: third side c via spherical law of cosines
    cos_c = math.cos(math.radians(great_circle_arc)) * math.cos(
        math.radians(triangle_side)
    ) + math.sin(math.radians(great_circle_arc)) * math.sin(
        math.radians(triangle_side)
    ) * math.cos(math.radians(triangle_angle))
    cos_c = max(-1.0, min(1.0, cos_c))
    side_c = math.degrees(math.acos(cos_c))

    # Step 3: zenith distance via spherical law of sines
    if side_c < 1e-6:
        zenith_dist = 0.0
    else:
        sin_zd = max(
            -1.0,
            min(
                1.0,
                math.sin(math.radians(great_circle_arc))
                * math.sin(math.radians(triangle_angle))
                / math.sin(math.radians(side_c)),
            ),
        )
        zenith_dist = math.degrees(math.asin(sin_zd))

    # Step 4: equatorial RA of the house-circle / equator intersection
    equator_ra_offset = math.degrees(
        math.atan(cos_lat * math.tan(math.radians(zenith_dist)))
    )

    # Pole height of the house circle
    pole_height = math.degrees(math.asin(math.sin(math.radians(zenith_dist)) * sin_lat))

    # Step 5: sign/direction convention and final RA
    if is_diurnal:
        ra_intersection = (armc + equator_ra_offset) % 360.0
    else:
        pole_height = -pole_height
        ra_intersection = (equator_ra_offset + armc + 180.0) % 360.0

    # Project equatorial intersection onto the ecliptic
    return _ra_to_ecliptic_longitude(
        ra_intersection, pole_height, sin_obliquity, cos_obliquity
    )


def _houses_sunshine(
    armc: float, lat: float, eps: float, asc: float, mc: float, sun_dec: float
) -> List[float]:
    """
    Sunshine (Makransky) house system (code 'I').

    Invented by Bob Makransky and published in 1988 in "Primary Directions".
    The diurnal and nocturnal arcs of the Sun are divided into thirds, and
    great circles through these trisection points and the horizon's north/south
    points define the house boundaries.  Each great circle is projected onto
    the ecliptic via spherical triangle arc-trisection geometry.

    Note: Cusps 11, 12, 2, 3 are NOT in exact opposition to cusps 5, 6, 8, 9.

    Algorithm:
        1. Compute the Sun's ascensional difference (AD) from its declination
           and the geographic latitude.
        2. Derive the diurnal semi-arc (DSA = 90° + AD) and the nocturnal
           semi-arc (NSA = 90° − AD); trisect each to get RA offsets.
        3. For nocturnal cusps (2, 3, 5, 6) call _sunshine_arc_to_ecliptic()
           with is_diurnal=False.
        4. For diurnal cusps (8, 9, 11, 12) call _sunshine_arc_to_ecliptic()
           with is_diurnal=True.
        5. If the MC is below the horizon, reflect all intermediate cusps by
           180°.

    Args:
        armc: Right Ascension of the Midheaven (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        sun_dec: Sun's declination in degrees

    Returns:
        List of 13 house cusp longitudes
    """
    NEAR_ZERO = 1e-10

    # Detect MC-below-horizon condition.
    # mcdec = declination of the MC, derived from ARMC and obliquity.
    # When |lat - mcdec| > 90°, the MC is below the horizon.
    mcdec = math.degrees(
        math.atan(math.sin(math.radians(armc)) * math.tan(math.radians(eps)))
    )
    mc_under_horizon = abs(lat - mcdec) > 90

    cusps = [0.0] * 13

    # When MC is below the horizon, flip MC and IC by 180° to keep MC
    # above the horizon (analogous to Regiomontanus polar handling).
    # The ASC is already hemisphere-corrected by the dispatch.
    if mc_under_horizon:
        mc = (mc + 180.0) % 360.0

    cusps[1] = asc
    cusps[10] = mc
    cusps[4] = (mc + 180.0) % 360.0
    cusps[7] = (asc + 180.0) % 360.0

    # Ascensional difference: sin(AD) = tan(δ☉) · tan(φ)
    ad_arg = math.tan(math.radians(sun_dec)) * math.tan(math.radians(lat))
    if ad_arg >= 1.0:
        ascensional_diff = 90.0 - NEAR_ZERO
    elif ad_arg <= -1.0:
        ascensional_diff = -90.0 + NEAR_ZERO
    else:
        ascensional_diff = math.degrees(math.asin(ad_arg))

    nocturnal_semi_arc = 90.0 - ascensional_diff  # NSA
    diurnal_semi_arc = 90.0 + ascensional_diff  # DSA

    # Pre-compute repeated trig quantities
    sin_lat = math.sin(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))
    cos_sun_dec = math.cos(math.radians(sun_dec))
    tan_sun_dec = math.tan(math.radians(sun_dec))
    sin_obliquity = math.sin(math.radians(eps))
    cos_obliquity = math.cos(math.radians(eps))

    # RA offsets along each semi-arc (positive = toward MC, negative = toward IC)
    # NSA offsets (nocturnal cusps 2, 3, 5, 6)
    nsa_third = nocturnal_semi_arc / 3.0
    nsa_offsets = {
        2: -2.0 * nsa_third,
        3: -1.0 * nsa_third,
        5: 1.0 * nsa_third,
        6: 2.0 * nsa_third,
    }
    # DSA offsets (diurnal cusps 8, 9, 11, 12)
    dsa_third = diurnal_semi_arc / 3.0
    dsa_offsets = {
        8: -2.0 * dsa_third,
        9: -1.0 * dsa_third,
        11: 1.0 * dsa_third,
        12: 2.0 * dsa_third,
    }

    # Shared kwargs for the geometric helper
    _common = dict(
        armc=armc,
        lat=lat,
        sun_dec=sun_dec,
        sin_lat=sin_lat,
        cos_lat=cos_lat,
        cos_sun_dec=cos_sun_dec,
        tan_sun_dec=tan_sun_dec,
        sin_obliquity=sin_obliquity,
        cos_obliquity=cos_obliquity,
    )

    # --- Nocturnal semi-arc cusps (below the horizon arc) ---
    for house_index, ra_offset in nsa_offsets.items():
        cusps[house_index] = _sunshine_arc_to_ecliptic(
            arc_offset_ra=ra_offset,
            is_diurnal=False,
            **_common,
        )

    # --- Diurnal semi-arc cusps (above the horizon arc) ---
    for house_index, ra_offset in dsa_offsets.items():
        cusps[house_index] = _sunshine_arc_to_ecliptic(
            arc_offset_ra=ra_offset,
            is_diurnal=True,
            **_common,
        )

    # When MC is below the horizon, shift all intermediate cusps by 180°.
    # This keeps the MC above the horizon, consistent with the polar-circle
    # handling used for Regiomontanus and Campanus systems.
    if mc_under_horizon:
        for i in (2, 3, 5, 6, 8, 9, 11, 12):
            cusps[i] = (cusps[i] + 180.0) % 360.0

    return cusps


def _houses_sunshine_makransky(
    armc: float, lat: float, eps: float, asc: float, mc: float, sun_dec: float
) -> List[float]:
    """
    Sunshine house system — Makransky variant (code 'i').

    Bob Makransky's original algorithm from *Primary Directions* (1988).
    Differs from Treindl's ('I') solution in the geometric projection:
    uses prime-vertical / meridian decomposition with explicit handling
    of meridian distance > 90° (which occurs at high latitudes in winter).

    The implementation follows the book's published computational
    procedure step by step, using the standard primary-directions
    quantities: MD (meridian distance), ZD (zenith distance), AD
    (ascensional difference), the pole height of a house circle, and the
    book's intermediate arcs and quadrant case tables.

    Raises PolarCircleError if the Sun is circumpolar (reference parity).

    Args:
        armc: Right Ascension of the Midheaven (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees
        sun_dec: Sun's declination in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused)
    """
    _NEAR_ZERO = 1e-10

    # MC-below-horizon check (same as Treindl variant)
    acmc = asc - mc
    while acmc > 180:
        acmc -= 360
    while acmc < -180:
        acmc += 360
    mc_under_horizon = acmc < 0
    if mc_under_horizon:
        asc = (asc + 180.0) % 360.0

    cusps = [0.0] * 13
    cusps[1] = asc
    cusps[10] = mc
    cusps[4] = (mc + 180.0) % 360.0
    cusps[7] = (asc + 180.0) % 360.0

    # Ascensional difference.  When |tan(dec)*tan(lat)| >= 1 the Sun is
    # circumpolar (never rises or never sets), the diurnal semi-arcs do
    # not exist, and the reference raises rather than substituting a
    # different house system (verified: pyswisseph errors for 'i' at
    # exactly |lat|+|dec| >= 90, e.g. 66.56/23.44, 80/11.47).
    ad_arg = math.tan(math.radians(sun_dec)) * math.tan(math.radians(lat))
    if abs(ad_arg) >= 1.0:
        raise PolarCircleError(
            message=(
                f"Sunshine houses undefined: the Sun (declination "
                f"{sun_dec:.4f}) is circumpolar at latitude {lat:.4f}"
            ),
            latitude=lat,
            house_system="i",
        )
    ascensional_diff = math.degrees(math.asin(ad_arg))

    night_semi_arc = 90.0 - ascensional_diff
    day_semi_arc = 90.0 + ascensional_diff

    # Equal trisection of the Sun's semi-arcs: each intermediate cusp sits
    # at a third of the night arc (houses 2/3/5/6) or of the day arc
    # (houses 8/9/11/12) away from the relevant meridian.
    arc_offset = [0.0] * 13
    arc_offset[2] = -2.0 * night_semi_arc / 3.0
    arc_offset[3] = -1.0 * night_semi_arc / 3.0
    arc_offset[5] = 1.0 * night_semi_arc / 3.0
    arc_offset[6] = 2.0 * night_semi_arc / 3.0
    arc_offset[8] = -2.0 * day_semi_arc / 3.0
    arc_offset[9] = -1.0 * day_semi_arc / 3.0
    arc_offset[11] = 1.0 * day_semi_arc / 3.0
    arc_offset[12] = 2.0 * day_semi_arc / 3.0

    sin_lat = math.sin(math.radians(lat))
    cos_lat = math.cos(math.radians(lat))
    tan_lat = math.tan(math.radians(lat))
    tan_dec = math.tan(math.radians(sun_dec))
    sin_ecl = math.sin(math.radians(eps))

    south = lat < 0

    for house in range(1, 13):
        if (house - 1) % 3 == 0:
            continue  # skip cardinal cusps 1, 4, 7, 10

        merid_dist = abs(arc_offset[house])

        # RA of the house point on the equator (night houses count from the
        # lower meridian, day houses from the upper)
        if house <= 6:
            ra_house = (armc + 180.0 + arc_offset[house]) % 360.0
        else:
            ra_house = (armc + arc_offset[house]) % 360.0

        if south:
            ra_house = (180.0 + ra_house) % 360.0

        # Zenith distance of the house point, built from the book's three
        # successive right spherical triangles (Primary Directions, 1988).
        if abs(merid_dist - 90.0) < _NEAR_ZERO:
            # Degenerate case: MD = 90°, the house point stands on the
            # east/west point and the triangles collapse.
            zenith_dist = 90.0 - math.degrees(math.atan(sin_lat * tan_dec))
        else:
            # Triangle 1 — from the meridian along the equator: the arc
            # between the zenith's meridian and the house point's vertical.
            if merid_dist < 90.0:
                zd_equator = math.degrees(
                    math.atan(cos_lat * math.tan(math.radians(merid_dist)))
                )
            else:
                # MD > 90°: the same triangle read from the opposite
                # (east-point) side.
                zd_equator = math.degrees(
                    math.atan(math.tan(math.radians(merid_dist - 90.0)) / cos_lat)
                )

            # Triangle 2 — the meridian arc between the equator and the
            # foot of the house point's vertical.
            merid_arc = math.degrees(math.atan(tan_lat * math.cos(math.radians(merid_dist))))

            # Meridian arc from the house point's declination parallel:
            # night houses add the Sun's declination, day houses subtract.
            if house <= 6:
                merid_arc_dec = merid_arc + sun_dec
            else:
                merid_arc_dec = merid_arc - sun_dec

            # Triangle 3 — correction carrying the declination offset back
            # onto the vertical of the house point.
            zd_correction = math.degrees(
                math.atan(
                    sin_lat
                    * math.sin(math.radians(merid_dist))
                    * math.tan(math.radians(merid_arc_dec))
                )
            )

            zenith_dist = zd_equator + zd_correction

        # Pole height of the house circle through the house point
        pole_height = math.degrees(
            math.asin(math.sin(math.radians(zenith_dist)) * sin_lat)
        )

        # Ascensional difference of the house circle's pole: the RA shift
        # between the house point and the circle's equator crossing.
        ad_arg = tan_dec * math.tan(math.radians(pole_height))
        ad_arg = max(-1.0, min(1.0, ad_arg))
        ad_pole = math.degrees(math.asin(ad_arg))

        # RA of the house circle's equator crossing (eastern houses shift
        # against the pole's AD, western houses with it)
        if house <= 3 or house >= 11:
            ra_crossing = (ra_house - ad_pole) % 360.0
        else:
            ra_crossing = (ra_house + ad_pole) % 360.0

        # Project (ra_crossing, pole_height) onto the ecliptic, resolving
        # quadrants with the book's case tables.
        if abs(ra_crossing - 90.0) < _NEAR_ZERO:
            ecl_arc = math.degrees(
                math.atan(sin_ecl * math.tan(math.radians(pole_height)))
            )
            if house <= 3 or house >= 11:
                cusp_lon = 90.0 + ecl_arc
            else:
                cusp_lon = 90.0 - ecl_arc
        elif abs(ra_crossing - 270.0) < _NEAR_ZERO:
            ecl_arc = math.degrees(
                math.atan(sin_ecl * math.tan(math.radians(pole_height)))
            )
            if house <= 3 or house >= 11:
                cusp_lon = 270.0 - ecl_arc
            else:
                cusp_lon = 270.0 + ecl_arc
        else:
            # Inclination of the house circle's meridian plane at the equator
            pole_merid_angle = math.degrees(
                math.atan(
                    abs(
                        math.tan(math.radians(pole_height))
                        / math.cos(math.radians(ra_crossing))
                    )
                )
            )
            # Tilt it into the ecliptic frame (sign per quadrant and house side)
            if house <= 3 or house >= 11:
                if 90.0 < ra_crossing < 270.0:
                    ecl_angle = pole_merid_angle - eps
                else:
                    ecl_angle = pole_merid_angle + eps
            else:
                if 90.0 < ra_crossing < 270.0:
                    ecl_angle = pole_merid_angle + eps
                else:
                    ecl_angle = pole_merid_angle - eps

            if abs(ecl_angle - 90.0) < _NEAR_ZERO:
                cusp_lon = 90.0 if ra_crossing < 180.0 else 270.0
            else:
                ecl_arc = math.degrees(
                    math.atan(
                        abs(
                            math.cos(math.radians(pole_merid_angle))
                            * math.tan(math.radians(ra_crossing))
                            / math.cos(math.radians(ecl_angle))
                        )
                    )
                )
                if ra_crossing < 90.0:
                    cusp_lon = ecl_arc
                elif 90.0 < ra_crossing < 180.0:
                    cusp_lon = 180.0 - ecl_arc
                elif 180.0 < ra_crossing < 270.0:
                    cusp_lon = 180.0 + ecl_arc
                else:
                    cusp_lon = 360.0 - ecl_arc

                # Reflected case table for ecl_angle > 90° (Makransky p. 146)
                if ecl_angle > 90.0:
                    if ra_crossing < 90.0:
                        cusp_lon = 180.0 - ecl_arc
                    elif 90.0 < ra_crossing < 180.0:
                        cusp_lon = ecl_arc
                    elif 180.0 < ra_crossing < 270.0:
                        cusp_lon = 360.0 - ecl_arc
                    else:
                        cusp_lon = 180.0 + ecl_arc

            if south:
                cusp_lon = (cusp_lon + 180.0) % 360.0

        cusps[house] = cusp_lon % 360.0

    return cusps


def _houses_horizontal(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Horizontal (Azimuthal) house system (code 'H').

    The horizon is divided into 12 equal arcs of 30 degrees each, starting
    from the East Point. For each division point, a vertical circle (great
    circle through zenith and nadir) passes through that point and intersects
    the ecliptic. Those intersection points define the house cusps.

    Algorithm derived from the geometric definition described in:
    Ralph William Holden, 'The Elements of House Division' (1977), ch. 10.

    Geometric derivation (standard spherical trigonometry):
        Let phi' = co-latitude (90 deg - |phi|, signed).
        For each horizon division at angle theta from the meridian
        (theta = k * 30 deg, k = 1..5), the vertical circle's parameters are:

            pole_height = arcsin(sin(theta) * sin(phi'))
            ra_offset   = atan2(cos(theta), sin(theta) * cos(phi'))

        These follow from solving the spherical triangle (zenith, celestial
        pole, horizon division point) via the spherical law of sines/cosines.

        The ecliptic intersection is then computed via the rising-point formula
        (Smart, 'Textbook on Spherical Astronomy', ch. 3).

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are cusps)
    """
    cusps = [0.0] * 13
    _NEAR_ZERO = 1e-10

    # Co-latitude: complement of geographic latitude
    if lat > 0:
        co_lat = 90.0 - lat
    else:
        co_lat = -90.0 - lat

    # Clamp to avoid singularity at the equator (|co_lat| = 90 deg)
    if abs(abs(co_lat) - 90.0) < _NEAR_ZERO:
        co_lat = (90.0 - _NEAR_ZERO) if co_lat > 0 else (-90.0 + _NEAR_ZERO)

    # ARMC rotated by 180 deg (base reference for house calculations)
    armc_base = (armc + 180.0) % 360.0

    sin_colat = math.sin(math.radians(co_lat))
    cos_colat = math.cos(math.radians(co_lat))

    # Compute 5 primary cusps by looping over horizon division angles.
    # The horizon is divided at theta = k * 30 deg (k = 1..5) measured from
    # the meridian toward the East Point along the horizon.
    #
    # Cusp mapping:  theta=30 -> cusp 3,  theta=60 -> cusp 2,
    #                theta=90 -> cusp 1 (East Point),
    #                theta=120 -> cusp 12, theta=150 -> cusp 11
    _CUSP_MAP = [(3, 30.0), (2, 60.0), (1, 90.0), (12, 120.0), (11, 150.0)]

    for cusp_idx, theta_deg in _CUSP_MAP:
        theta_rad = math.radians(theta_deg)
        sin_theta = math.sin(theta_rad)
        cos_theta = math.cos(theta_rad)

        # Pole height of the vertical circle at this division angle
        pole_height = math.degrees(math.asin(sin_theta * sin_colat))

        # RA offset from the base RA position (ARMC + 270 deg)
        ra_offset = math.degrees(math.atan2(cos_theta, sin_theta * cos_colat))

        cusps[cusp_idx] = _calc_ascendant(
            armc_base + 90.0 + ra_offset, eps, lat, pole_height
        )

    # Polar circle handling: check Asc-MC orientation
    if abs(co_lat) >= 90.0 - eps:
        acmc_diff = (asc - mc + 540.0) % 360.0 - 180.0
        if acmc_diff < 0:
            asc = (asc + 180.0) % 360.0
            mc = (mc + 180.0) % 360.0
            for i in range(1, 13):
                if 4 <= i < 10:
                    continue
                cusps[i] = (cusps[i] + 180.0) % 360.0

    # Flip primary cusps by 180 deg (convention: cusps represent the opposite
    # hemisphere from the direct calculation)
    for i in [1, 2, 3, 11, 12]:
        cusps[i] = (cusps[i] + 180.0) % 360.0

    # Check Asc/DC orientation
    acmc_diff = (asc - mc + 540.0) % 360.0 - 180.0
    if acmc_diff < 0:
        asc = (asc + 180.0) % 360.0

    # Set MC/IC and derive opposite cusps (4-9 from 10-3)
    cusps[10] = mc
    cusps[4] = (mc + 180.0) % 360.0
    cusps[7] = (cusps[1] + 180.0) % 360.0
    cusps[8] = (cusps[2] + 180.0) % 360.0
    cusps[9] = (cusps[3] + 180.0) % 360.0
    cusps[5] = (cusps[11] + 180.0) % 360.0
    cusps[6] = (cusps[12] + 180.0) % 360.0

    return cusps


def _houses_natural_gradient(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    Natural Gradient house system ('N').

    'N' maps to "Equal houses with 0° Aries as cusp 1".
    This is effectively a Whole Sign system starting from 0° Aries.

    Mathematical Formula:
        λᵢ = (i-1) × 30°
        where i = 1, 2, ..., 12

    Result:
        House 1 = 0° Aries (0°)
        House 2 = 0° Taurus (30°)
        House 3 = 0° Gemini (60°)
        ... etc.

    Properties:
        - Independent of time and location
        - Works at all latitudes
        - Computationally exact

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees (unused)
        lat: Geographic latitude in degrees (unused)
        eps: True obliquity of ecliptic in degrees (unused)
        asc: Ascendant longitude in degrees (unused)
        mc: Midheaven longitude in degrees (unused)

    Returns:
        List of 13 house cusp longitudes
    """
    cusps = [0.0] * 13
    for i in range(1, 13):
        cusps[i] = ((i - 1) * 30.0) % 360.0
    return cusps


def _apc_cusp(house: int, lat_rad: float, eps_rad: float, armc_rad: float) -> float:
    """
    Ecliptic longitude of one APC (Ascendant-Parallel Circle) house cusp.

    The APC system divides the declination parallel through the Ascendant
    into equal arcs: the diurnal part (Asc over the meridian to Dsc) and the
    nocturnal part each into six. Position circles through the North and
    South points of the horizon pass through each division point; their
    intersections with the ecliptic define the cusps.

    System definition published by L. Knegt (WvA/Ram school, Netherlands),
    'Handleiding voor het berekenen van huizentabellen'.

    Derivation from spherical trigonometry
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Let φ be the geographic latitude, ε the obliquity, θ the ARMC.

    *Ascensional difference of the Ascendant.*  The Ascendant's declination
    δ_A and its ascensional difference K (the hour-angle offset between its
    rising point and the east point) satisfy the classical relations
    sin K = tan φ · tan δ_A and the rising condition RA_A = θ + 90° − K.
    Substituting the Ascendant's position on the ecliptic and eliminating
    RA_A yields the closed form

        K   = arctan( tan φ · tan ε · cos θ / (1 + tan φ · tan ε · sin θ) )
        δ_A = arctan( sin K / tan φ )

    (for φ → 0 the parallel through the Ascendant approaches the equator's
    pole distance, handled as a limit below).

    *Division points.*  On the declination parallel the Ascendant sits at
    hour angle −(90° + K) (east horizon). Walking the parallel by equal
    sixths of the nocturnal arc (90° − K each, houses 1-7 eastward below
    the horizon) or of the diurnal arc (90° + K each, houses 8-12) puts
    division point k at right ascension

        a_k = θ + 90° + K + k · (90° ∓ K) / 3.

    *Projection.*  The position circle through the horizon North/South
    points and the division point (a_k, δ_A) crosses the ecliptic at

        λ = atan2( tan δ_A · tan φ · sin θ + sin a_k,
                   cos ε · (tan δ_A · tan φ · cos θ + cos a_k)
                     + sin ε · tan φ · sin(θ − a_k) ).

    Args:
        house: House number (1-12)
        lat_rad: Geographic latitude in radians
        eps_rad: Obliquity of the ecliptic in radians
        armc_rad: ARMC in radians

    Returns:
        Ecliptic longitude of the house cusp in degrees
    """
    _NEAR_ZERO = 1e-10

    # Ascensional difference and declination of the Ascendant
    if abs(math.degrees(lat_rad)) > 90 - _NEAR_ZERO:
        asc_diff = 0.0
        asc_declination = 0.0
    else:
        asc_diff = math.atan(
            math.tan(lat_rad)
            * math.tan(eps_rad)
            * math.cos(armc_rad)
            / (1 + math.tan(lat_rad) * math.tan(eps_rad) * math.sin(armc_rad))
        )

        if abs(math.degrees(lat_rad)) < _NEAR_ZERO:
            # Equator limit: the parallel through the Ascendant degenerates;
            # use a near-pole declination with the Ascendant's sign.
            asc_declination = math.radians(90.0 - _NEAR_ZERO)
            if lat_rad < 0:
                asc_declination = -asc_declination
        else:
            asc_declination = math.atan(math.sin(asc_diff) / math.tan(lat_rad))

    # Step index along the parallel: houses 1-7 walk the nocturnal arc
    # eastward from the Ascendant, houses 8-12 the diurnal arc.
    if house < 8:
        below_horizon = True
        step = house - 1
    else:
        below_horizon = False
        step = house - 13

    # Right ascension of the division point on the Ascendant's parallel
    half_pi = math.pi / 2
    if below_horizon:
        ra_division = (
            asc_diff + armc_rad + half_pi + step * (half_pi - asc_diff) / 3
        )
    else:
        ra_division = (
            asc_diff + armc_rad + half_pi + step * (half_pi + asc_diff) / 3
        )

    ra_division = ra_division % (2 * math.pi)

    # Ecliptic crossing of the position circle through the division point
    longitude = math.atan2(
        math.tan(asc_declination) * math.tan(lat_rad) * math.sin(armc_rad)
        + math.sin(ra_division),
        math.cos(eps_rad)
        * (
            math.tan(asc_declination) * math.tan(lat_rad) * math.cos(armc_rad)
            + math.cos(ra_division)
        )
        + math.sin(eps_rad) * math.tan(lat_rad) * math.sin(armc_rad - ra_division),
    )

    return math.degrees(longitude) % 360.0


def _houses_apc(
    armc: float, lat: float, eps: float, asc: float, mc: float
) -> List[float]:
    """
    APC (Ascendant-Parallel Circle) house system (code 'Y').

    The parallel of declination passing through the Ascendant is divided
    into 6 equal arcs above and 6 below the horizon. Position circles
    through the North and South points of the horizon pass through each
    division point; their intersections with the ecliptic define the cusps.

    Algorithm derived from the geometric definition by L. Knegt
    (WvA/Ram school, Netherlands).

    Reference: Knegt, 'Handleiding voor het berekenen van huizentabellen'
    (Manual for computing house tables).

    Args:
        armc: Sidereal time at Greenwich (RAMC) in degrees
        lat: Geographic latitude in degrees
        eps: True obliquity of ecliptic in degrees
        asc: Ascendant longitude in degrees
        mc: Midheaven longitude in degrees

    Returns:
        List of 13 house cusp longitudes (index 0 unused, 1-12 are cusps)
    """
    cusps = [0.0] * 13

    # Per-cusp computation works in radians
    lat_rad = math.radians(lat)
    eps_rad = math.radians(eps)
    armc_rad = math.radians(armc)

    # Calculate all 12 house cusps
    for i in range(1, 13):
        cusps[i] = _apc_cusp(i, lat_rad, eps_rad, armc_rad)

    # The per-cusp formula is ill-conditioned for the MC near |lat| = 90;
    # the directly computed MC is exact, so always substitute it.
    cusps[10] = mc
    cusps[4] = (mc + 180.0) % 360.0

    # Within polar circle, handle horizon crossing
    if abs(lat) >= 90.0 - eps:
        if difdeg2n(asc, mc) < 0:
            # Flip all cusps by 180° (clockwise direction)
            asc = (asc + 180.0) % 360.0
            mc = (mc + 180.0) % 360.0
            for i in range(1, 13):
                cusps[i] = (cusps[i] + 180.0) % 360.0

    return cusps


@overload
def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Tuple[float, float],
    lon_or_hsys: Union[bytes, str],
    lat_body: float = ...,
) -> float: ...


@overload
def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str],
    lon_or_hsys: float,
    lat_body: float = ...,
) -> float: ...


@overload
def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Tuple[float, float],
    lon_or_hsys: None = ...,
    lat_body: float = ...,
) -> float: ...


@overload
def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str],
    lon_or_hsys: None = ...,
    lat_body: float = ...,
) -> float: ...


def _house_pos_pythonic(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str, Tuple[float, float]],
    lon_or_hsys: Optional[Union[float, bytes, str]] = None,
    lat_body: float = 0.0,
) -> float:
    """
    Determine in which house a celestial body is located.

    Returns a decimal value where the integer part is the house number (1-12)
    and the decimal part indicates the position within the house
    (0.0 = start of cusp, 0.999... = end of house, just before next cusp).

    This function supports two calling conventions:

    1. Reference API compatible (5 args):
       _house_pos_pythonic(armc, lat, obliquity, objcoord, hsys)
       where objcoord is a tuple (lon, lat_body) and hsys is bytes/str

    2. Extended signature (6 args):
       _house_pos_pythonic(armc, lat, obliquity, hsys, lon, lat_body)
       where hsys is int/bytes/str and lon/lat_body are separate floats

    Args:
        armc: Right Ascension of Medium Coeli (ARMC) in degrees (0-360)
        lat: Geographic latitude in degrees (positive North, negative South)
        obliquity: True obliquity of the ecliptic in degrees
        hsys_or_objcoord: Either house system (int/bytes/str) or objcoord tuple (lon, lat)
        lon_or_hsys: Either body longitude (float) or house system (bytes/str)
        lat_body: Ecliptic latitude of the body in degrees (only for 6-arg form)

    Returns:
        Decimal value where:
            - Integer part (1-12): House number
            - Decimal part (0.0-0.999...): Position within house

    Example:
        >>> # Sun at 15° Aries, Placidus houses, Rome (5-arg reference API form)
        >>> pos = _house_pos_pythonic(292.957, 41.9, 23.4393, (15.0, 0.0), b'P')
        >>> # Or 6-arg extended form:
        >>> pos = _house_pos_pythonic(292.957, 41.9, 23.4393, ord('P'), 15.0, 0.0)
        >>> house = int(pos)  # House number (e.g., 10)
        >>> position = pos - house  # Position within house (e.g., 0.5 = halfway)
    """
    _NEAR_ZERO = 1e-10
    # Tiny offset (~0.28 microdeg) to avoid bodies landing exactly on a cusp boundary
    CUSP_BOUNDARY_OFFSET = 1.0 / 3600000.0

    # Declare typed local variables for proper type narrowing
    lon: float
    hsys_char: str
    hsys_int: int

    # Detect which calling convention is used.  objcoord may be any
    # sequence (the reference accepts lists as well as tuples).
    if isinstance(hsys_or_objcoord, (tuple, list)):
        # 5-arg reference API form: (armc, lat, obliquity, objcoord, hsys)
        objcoord = hsys_or_objcoord
        lon = objcoord[0]
        lat_body = objcoord[1] if len(objcoord) > 1 else 0.0
        # hsys comes from lon_or_hsys (bytes/str), default to b"P"
        if isinstance(lon_or_hsys, bytes):
            hsys_char = chr(lon_or_hsys[0])
            hsys_int = lon_or_hsys[0]
        elif isinstance(lon_or_hsys, str):
            hsys_char = lon_or_hsys[0]
            hsys_int = ord(lon_or_hsys[0])
        else:
            # Default to Placidus
            hsys_char = "P"
            hsys_int = ord("P")
    else:
        # 6-arg form: (armc, lat, obliquity, hsys, lon, lat_body)
        # hsys comes from hsys_or_objcoord (int/bytes/str)
        if isinstance(hsys_or_objcoord, bytes):
            hsys_char = chr(hsys_or_objcoord[0])
            hsys_int = hsys_or_objcoord[0]
        elif isinstance(hsys_or_objcoord, str):
            hsys_char = hsys_or_objcoord[0]
            hsys_int = ord(hsys_or_objcoord[0])
        else:
            # int case
            hsys_char = chr(hsys_or_objcoord)
            hsys_int = hsys_or_objcoord
        # lon comes from lon_or_hsys: a float, or an objcoord sequence
        # in the hsys-first calling order
        if isinstance(lon_or_hsys, (tuple, list)):
            lon = float(lon_or_hsys[0])
            if len(lon_or_hsys) > 1:
                lat_body = float(lon_or_hsys[1])
        else:
            lon = float(lon_or_hsys) if lon_or_hsys is not None else 0.0

    # Normalize inputs
    lon = lon % 360.0
    eps = obliquity
    geolat = lat

    # Convert ecliptic coordinates to equatorial via _rotate_spherical_x_axis rotation
    # (rotation around x-axis by obliquity angle)
    # This is crucial for proper house position when body has non-zero latitude
    lon_rad = math.radians(lon)
    lat_body_rad = math.radians(lat_body)
    eps_rad = math.radians(eps)

    # Ecliptic to equatorial transformation
    # RA: tan(RA) = (sin(lon)*cos(eps) - tan(lat)*sin(eps)) / cos(lon)
    # Dec: sin(Dec) = sin(lat)*cos(eps) + cos(lat)*sin(lon)*sin(eps)
    sin_lon = math.sin(lon_rad)
    cos_lon = math.cos(lon_rad)
    sin_lat = math.sin(lat_body_rad)
    cos_lat = math.cos(lat_body_rad)
    sin_eps = math.sin(eps_rad)
    cos_eps = math.cos(eps_rad)

    # Right Ascension
    y_ra = sin_lon * cos_eps - math.tan(lat_body_rad) * sin_eps
    x_ra = cos_lon
    ra = math.degrees(math.atan2(y_ra, x_ra)) % 360.0

    # Declination
    sin_dec = sin_lat * cos_eps + cos_lat * sin_lon * sin_eps
    if sin_dec > 1.0:
        sin_dec = 1.0
    if sin_dec < -1.0:
        sin_dec = -1.0
    dec = math.degrees(math.asin(sin_dec))

    # Signed meridian distances of the body, in (-180, 180]:
    #   md_upper — RA distance from the upper meridian (culmination, RA = ARMC)
    #   md_lower — RA distance from the lower meridian (anti-culmination)
    md_upper = (ra - armc) % 360.0
    if md_upper >= 180:
        md_upper -= 360.0
    md_lower = (md_upper + 180.0) % 360.0
    if md_lower >= 180:
        md_lower -= 360.0

    # Handle different house systems
    if hsys_char == "P" or hsys_char == "G":
        # Placidus / Gauquelin house position.
        # The Placidus position of a body is the fraction of its own
        # semi-arc already traversed, mapped onto a 360° house circle
        # (Holden, 'The Elements of House Division', 1977):
        #
        #   above horizon: position = (md_upper / day_semi_arc + 3) * 90°
        #   below horizon: position = (md_lower / night_semi_arc + 1) * 90°

        if 90.0 - abs(dec) <= abs(geolat):
            # Circumpolar body: it never rises or sets, so no semi-arcs
            # exist.  A never-rising body is placed in the nocturnal half
            # by its half-weighted distance from lower culmination; a
            # never-setting body in the diurnal half from upper culmination.
            if dec * geolat < 0:
                pos_deg = (90.0 + md_lower / 2.0) % 360.0
            else:
                pos_deg = (270.0 + md_upper / 2.0) % 360.0
        else:
            # Ascensional difference of the body: sin(AD) = tan δ · tan φ
            sin_ad = math.tan(math.radians(dec)) * math.tan(math.radians(geolat))
            if abs(sin_ad) > 1.0:
                sin_ad = 1.0 if sin_ad > 0 else -1.0
            asc_diff = math.degrees(math.asin(sin_ad))

            # Above-horizon test.  The altitude h of the body satisfies
            #   sin h = sin φ sin δ + cos φ cos δ cos H;
            # dividing by cos φ cos δ (> 0) gives the equivalent sign test
            #   tan φ tan δ + cos H  =  sin(AD) + cos H  >=  0.
            horizon_test = sin_ad + math.cos(math.radians(md_upper))
            above_horizon = horizon_test >= 0

            day_semi_arc = 90.0 + asc_diff
            night_semi_arc = 90.0 - asc_diff

            if above_horizon:
                pos_deg = (md_upper / day_semi_arc + 3.0) * 90.0
            else:
                pos_deg = (md_lower / night_semi_arc + 1.0) * 90.0

            # Add small offset for cusp precision
            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0

        if hsys_char == "G":
            # Gauquelin sectors run clockwise
            pos_deg = 360.0 - pos_deg
            hpos = pos_deg / 10.0 + 1.0
        else:
            hpos = pos_deg / 30.0 + 1.0

        return hpos

    elif hsys_char == "R":
        # Regiomontanus house position.
        #
        # Regiomontanus divides the celestial equator into equal arcs by
        # position circles through the North and South points of the
        # horizon; a body's position is the equator point cut by its own
        # position circle.  Intersecting the plane of the circle through
        # the horizon N/S points and the body (meridian distance M,
        # declination δ) with the equator gives the closed form
        #
        #   tan x = -(tan φ tan δ + cos M) / sin M
        #
        # measured so that x = 270° on the upper meridian (M = 0) and
        # x = 90° on the lower (|M| = 180°) — those two degenerate cases
        # are handled directly.
        # Degenerate meridian cases: the equator intersection is 270
        # (house 10 side) or 90 (house 4 side) according to whether the
        # meridian point is above or below the horizon — at polar
        # latitudes the upper-meridian point can be below it (altitude
        # 90 - |lat - dec| < 0 ⟺ tan(lat) tan(dec) < -1) and the limit
        # of the closed form lands on the opposite branch (verified:
        # pyswisseph puts the MC degree in house 4 at armc 90,
        # lat -80).
        _tt = math.tan(math.radians(geolat)) * math.tan(math.radians(dec))
        if abs(md_upper) < _NEAR_ZERO:
            pos_deg = 270.0 if _tt > -1.0 else 90.0
        elif 180.0 - abs(md_upper) < _NEAR_ZERO:
            pos_deg = 90.0 if _tt < 1.0 else 270.0
        else:
            if 90.0 - abs(geolat) < _NEAR_ZERO:
                geolat = 90.0 - _NEAR_ZERO if geolat > 0 else -90.0 + _NEAR_ZERO
            if 90.0 - abs(dec) < _NEAR_ZERO:
                dec = 90.0 - _NEAR_ZERO if dec > 0 else -90.0 + _NEAR_ZERO

            pole_term = math.tan(math.radians(geolat)) * math.tan(
                math.radians(dec)
            ) + math.cos(math.radians(md_upper))
            pos_deg = math.degrees(
                math.atan(-pole_term / math.sin(math.radians(md_upper)))
            )
            if md_upper < 0:
                pos_deg += 180.0
            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0

        hpos = pos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char == "C":
        # Campanus house position.
        # The Campanus system divides the prime vertical into 12 equal
        # 30° arcs.  A body's position follows from rotating its
        # equatorial coordinates into the frame whose fundamental plane
        # is the prime vertical (pole = horizon North point):
        #
        # 1. Express the body as (meridian distance - 90°, declination),
        #    i.e. longitude measured from the east point of the equator.
        # 2. Rotate about the shared x-axis (the east point) by -latitude,
        #    which carries the equator onto the prime vertical.
        # 3. The longitude in the rotated frame is the Campanus position.
        #
        # Reference: Meeus, 'Astronomical Algorithms' 2nd ed., Ch. 13;
        # Montenbruck & Pfleger, par. 1.3 (x-axis rotation matrix).

        ha_east = md_upper - 90.0

        ha_east_rad = math.radians(ha_east)
        dec_rad = math.radians(dec)
        rot_angle = math.radians(-geolat)

        cos_rot = math.cos(rot_angle)
        sin_rot = math.sin(rot_angle)
        cos_ha = math.cos(ha_east_rad)
        sin_ha = math.sin(ha_east_rad)
        cos_dec = math.cos(dec_rad)
        sin_dec_b = math.sin(dec_rad)

        # Convert to Cartesian
        vec_x = cos_dec * cos_ha
        vec_y = cos_dec * sin_ha
        vec_z = sin_dec_b

        # Rotation around x-axis:
        #   x' = x;  y' = y cos θ + z sin θ;  z' = -y sin θ + z cos θ
        rot_x = vec_x
        rot_y = vec_y * cos_rot + vec_z * sin_rot

        # Convert back to spherical (longitude only needed)
        pos_deg = math.degrees(math.atan2(rot_y, rot_x))
        pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0
        hpos = pos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char in ["E", "A", "W", "V", "D", "N"]:
        # Equal-based systems - simple longitude-based calculation
        asc = _calc_ascendant((armc + 90.0) % 360.0, eps, geolat, geolat)
        mc = _armc_to_mc(armc, eps)

        if hsys_char == "D":
            pos_deg = (lon - mc - 90.0) % 360.0
        elif hsys_char == "V":
            pos_deg = (lon - asc + 15.0) % 360.0
        elif hsys_char == "W":
            pos_deg = (lon - asc + (asc % 30.0)) % 360.0
        elif hsys_char == "N":
            pos_deg = lon
        else:
            pos_deg = (lon - asc) % 360.0

        pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0
        hpos = pos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char == "X":
        # Meridian (axial rotation)
        hpos = ((md_upper - 90.0) % 360.0) / 30.0 + 1.0
        return hpos

    elif hsys_char == "M":
        # Morinus
        # Project the ecliptic longitude onto the equatorial frame via:
        # tan(ra_equiv) = tan(λ) / cos(ε)
        # atan2(sin(λ), cos(λ)·cos(ε)) gives the RA-equivalent in [-180, 180],
        # then % 360 maps to [0, 360).
        a = lon
        if abs(a - 90.0) > _NEAR_ZERO and abs(a - 270.0) > _NEAR_ZERO:
            hpos_deg = (
                math.degrees(
                    math.atan2(
                        math.sin(math.radians(a)), math.cos(math.radians(a)) * cos_eps
                    )
                )
                % 360.0
            )
        else:
            hpos_deg = 90.0 if abs(a - 90.0) <= _NEAR_ZERO else 270.0
        hpos_deg = (hpos_deg - armc - 90.0) % 360.0
        hpos = hpos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char == "K":
        # Koch (Birthplace / GOH) house position.
        # The Koch system defines every house circle by the *birthplace
        # horizon*: each cusp rises after the MC by an equal share of the
        # MC degree's own semi-arc.  A body's position is therefore the
        # share of the MC's semi-arc that has elapsed between the body's
        # rising (its meridian distance corrected by its own ascensional
        # difference) and the MC's rising (corrected by the MC's AD):
        #
        #   ad_body      = arcsin(tan φ · tan δ_body)
        #   ad_mc        = arcsin(tan ε · tan φ · sin ARMC)
        #                  (δ_MC from sin δ = sin ε · sin λ_MC, folded into
        #                   one expression via the MC's RA = ARMC)
        #   mc_semi_arc  = 90° + ad_mc
        #
        #   east  (md_upper >= 0): fraction = (md_upper - ad_body + ad_mc)
        #                                       / mc_semi_arc
        #                          position = (fraction - 1) · 90°
        #   west  (md_upper <  0): fraction = (md_upper + 180 + ad_body
        #                                       + ad_mc) / mc_semi_arc
        #                          position = (fraction + 1) · 90°
        #
        # Bodies whose fraction falls outside [0, 2] never cross the
        # birthplace horizon inside this quadrant pair — the Koch position
        # is undefined and 0.0 is returned (matching the reference API).

        tan_eps = math.tan(math.radians(eps))

        invalid = False

        # Body's ascensional difference, saturated at ±90° for declinations
        # that never rise / never set at this latitude.
        if 90.0 - geolat < dec or -90.0 - geolat > dec:
            ad_body = 90.0
        elif geolat - 90.0 > dec or geolat + 90.0 < dec:
            ad_body = -90.0
        else:
            ad_arg = math.tan(math.radians(geolat)) * math.tan(math.radians(dec))
            ad_arg = max(-1.0, min(1.0, ad_arg))
            ad_body = math.degrees(math.asin(ad_arg))
            if abs(ad_body) < 1e-12:
                # An ascensional difference of ~1e-15 deg is float noise
                # from a declination that is physically zero (body on
                # the equinoctial points); its sign must not decide the
                # validity boundary below.
                ad_body = 0.0

        # MC's ascensional difference
        ad_mc_arg = (
            tan_eps * math.tan(math.radians(geolat)) * math.sin(math.radians(armc))
        )
        if abs(ad_mc_arg) > 1.0:
            ad_mc_arg = 1.0 if ad_mc_arg > 0 else -1.0
        ad_mc = math.degrees(math.asin(ad_mc_arg))
        if abs(ad_mc) < 1e-12:
            # Same noise floor as ad_body: sin(180 deg) = 1.2e-16 makes
            # the MC's AD ~1e-15 deg at cardinal ARMC; physically zero.
            ad_mc = 0.0

        # MC's semi-arc
        mc_semi_arc = 90.0 + ad_mc
        if mc_semi_arc == 0.0:
            invalid = True

        pos_deg = 0.0
        if not invalid and abs(mc_semi_arc) > 0:
            # Small tolerance for floating-point boundary checks.
            # When a body is exactly on the MC or IC, the fraction can be
            # -1e-17 instead of 0.0 due to IEEE 754 rounding.
            # Both boundaries are strict, matching the reference: a
            # fraction of exactly 0.0 (body on the meridian with zero
            # AD) is valid, but -1e-17 from a cancelling md/AD pair or
            # a fraction at/above 2 is classified undefined.  Bodies
            # exactly on the validity edge are decided by the same
            # comparisons the reference applies.
            if md_upper >= 0:  # east
                arc_fraction = (md_upper - ad_body + ad_mc) / mc_semi_arc
                if arc_fraction >= 2.0 or arc_fraction < 0.0:
                    invalid = True
                else:
                    pos_deg = (arc_fraction - 1.0) * 90.0
            else:  # west
                md_from_lower = md_upper + 180.0
                arc_fraction = (md_from_lower + ad_body + ad_mc) / mc_semi_arc
                if arc_fraction >= 2.0 or arc_fraction < 0.0:
                    invalid = True
                else:
                    pos_deg = (arc_fraction + 1.0) * 90.0

        if invalid:
            # Koch position undefined in the circumpolar area — return 0.0
            # matching the reference API behavior
            return 0.0
        else:
            pos_deg = pos_deg % 360.0
            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0
            hpos = pos_deg / 30.0 + 1.0
            return hpos

    elif hsys_char == "T":
        # Topocentric (Polich-Page) house position.
        # For house position purposes, Topocentric uses the same semi-arc
        # method as Placidus — the difference between the systems is only
        # in how cusps are computed (different pole heights per house),
        # not in how a body is placed within those cusps.
        if 90.0 - abs(dec) <= abs(geolat):
            if dec * geolat < 0:
                pos_deg = (90.0 + md_lower / 2.0) % 360.0
            else:
                pos_deg = (270.0 + md_upper / 2.0) % 360.0
        else:
            sin_ad = math.tan(math.radians(dec)) * math.tan(math.radians(geolat))
            if abs(sin_ad) > 1.0:
                sin_ad = 1.0 if sin_ad > 0 else -1.0
            asc_diff = math.degrees(math.asin(sin_ad))
            horizon_test = sin_ad + math.cos(math.radians(md_upper))
            above_horizon = horizon_test >= 0
            day_semi_arc = 90.0 + asc_diff
            night_semi_arc = 90.0 - asc_diff
            if above_horizon:
                pos_deg = (md_upper / day_semi_arc + 3.0) * 90.0
            else:
                pos_deg = (md_lower / night_semi_arc + 1.0) * 90.0
            pos_deg = (pos_deg + CUSP_BOUNDARY_OFFSET) % 360.0
        hpos = pos_deg / 30.0 + 1.0
        return hpos

    elif hsys_char == "B":
        # Alcabitius house position.
        # Alcabitius cusps are equally spaced in RA within each quadrant.
        # For _house_pos_pythonic, we interpolate the body's RA between the cusp RAs.
        cusps, _ = houses_armc(armc, lat, obliquity, hsys_int)

        # Convert cusp ecliptic longitudes to RA
        cusp_ras = []
        for c in cusps:
            c_r = math.radians(c)
            y_c = math.sin(c_r) * cos_eps
            x_c = math.cos(c_r)
            cusp_ras.append(math.degrees(math.atan2(y_c, x_c)) % 360.0)

        # Find which RA interval contains the body's RA
        for i in range(12):
            ra_start = cusp_ras[i]
            ra_end = cusp_ras[(i + 1) % 12]

            diff_to_body = (ra - ra_start + 360.0) % 360.0
            interval_size = (ra_end - ra_start + 360.0) % 360.0

            if interval_size < 0.0001:
                interval_size = 30.0

            if diff_to_body < interval_size:
                house_num = i + 1
                fraction = diff_to_body / interval_size
                fraction = max(0.0, min(fraction, 0.9999999999))
                return float(house_num) + fraction

        return 1.0

    # Default fallback: use cusp-based method (ecliptic longitude interpolation)
    # This applies to H, F, I, i, J, S, L, Q, Y, and others
    cusps, ascmc = houses_armc(armc, lat, obliquity, hsys_int)

    # Find the house containing this longitude
    for i in range(12):
        cusp_start = cusps[i]
        cusp_end = cusps[(i + 1) % 12]

        diff_to_body = (lon - cusp_start + 360.0) % 360.0
        house_size = (cusp_end - cusp_start + 360.0) % 360.0

        if house_size < 1e-9:
            # Zero-width house (coincident cusps happen at extreme
            # latitudes, e.g. Pullen SD at lat 89 collapses houses 2-3).
            # A body can only be "in" it when it sits exactly on the
            # cusp; treating the empty house as 30 deg wide used to
            # swallow a whole sign of sky into it.
            if diff_to_body < 1e-9:
                return float(i + 1)
            continue

        if diff_to_body < house_size:
            house_num = i + 1
            fraction = diff_to_body / house_size
            fraction = max(0.0, min(fraction, 0.9999999999))
            return float(house_num) + fraction

    return 1.0


def _armc_to_mc(armc: float, eps: float) -> float:
    """Convert ARMC to MC longitude."""
    mc_rad = math.atan2(math.tan(math.radians(armc)), math.cos(math.radians(eps)))
    mc = math.degrees(mc_rad)
    if mc < 0:
        mc += 360.0
    if 90.0 < armc <= 270.0:
        if mc < 90.0 or mc > 270.0:
            mc += 180.0
    elif armc > 270.0:
        if mc < 270.0:
            mc += 180.0
    elif armc <= 90.0:
        if mc > 90.0:
            mc += 180.0
    return mc % 360.0


@overload
def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Tuple[float, float],
    lon_or_hsys: Union[bytes, str],
    lat_body: float = ...,
) -> float: ...


@overload
def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str],
    lon_or_hsys: float,
    lat_body: float = ...,
) -> float: ...


@overload
def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Tuple[float, float],
    lon_or_hsys: None = ...,
    lat_body: float = ...,
) -> float: ...


@overload
def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str],
    lon_or_hsys: None = ...,
    lat_body: float = ...,
) -> float: ...


def house_pos(
    armc: float,
    lat: float,
    obliquity: float,
    hsys_or_objcoord: Union[int, bytes, str, Tuple[float, float]],
    lon_or_hsys: Optional[Union[float, bytes, str]] = None,
    lat_body: float = 0.0,
) -> float:
    """
    Determine in which house a celestial body is located.

    This function supports two calling conventions:

    1. Reference API compatible (5 args):
       house_pos(armc, lat, obliquity, objcoord, hsys)
       where objcoord is a tuple (lon, lat_body) and hsys is bytes/str

    2. Extended signature (6 args):
       house_pos(armc, lat, obliquity, hsys, lon, lat_body)
       where hsys is int/bytes/str and lon/lat_body are separate floats

    Returns a decimal value where the integer part is the house number (1-12)
    and the decimal part indicates the position within the house.

    Args:
        armc: Right Ascension of Medium Coeli (ARMC) in degrees
        lat: Geographic latitude in degrees
        obliquity: True obliquity of the ecliptic in degrees
        hsys_or_objcoord: Either house system (int/bytes/str) or objcoord tuple
        lon_or_hsys: Either body longitude (float) or house system (bytes/str)
        lat_body: Body's ecliptic latitude in degrees (only for 6-arg form)

    Returns:
        Decimal value: integer part = house number, decimal part = position within house
    """
    # Delegate to _house_pos_pythonic which now supports both calling conventions
    # Use type: ignore since the overloads cover all valid combinations
    return _house_pos_pythonic(armc, lat, obliquity, hsys_or_objcoord, lon_or_hsys, lat_body)  # type: ignore[arg-type]


def _gauquelin_sector_from_rise_set(
    jd: float,
    planet: int,
    lat: float,
    lon: float,
    altitude: float,
    pressure: float,
    temperature: float,
    flags: int,
    method: int,
) -> float:
    """
    Calculate Gauquelin sector using actual rise/set times.

    This implements methods 2-5 which use real rise/set event times
    rather than the hour angle approximation.

    Args:
        jd: Julian Day in UT
        planet: Planet ID
        lat, lon: Observer location
        altitude, pressure, temperature: Atmospheric parameters
        flags: Calculation flags
        method: 2=disc center no refraction, 3=disc center with refraction,
                4=disc edge no refraction, 5=disc edge with refraction

    Returns:
        Sector position in range [1, 37)
    """
    from .eclipse import rise_trans
    from .constants import (
        CALC_RISE,
        CALC_SET,
        BIT_DISC_CENTER,
        BIT_NO_REFRACTION,
    )

    def fallback() -> float:
        """Fall back to method 0 (hour angle approximation)."""
        return _gauquelin_sector_pythonic(
            jd, planet, lat, lon, altitude, pressure, temperature, flags, 0
        )

    # Determine rise_trans flags based on method
    rsmi_flags = 0
    if method in (2, 3):  # Disc center
        rsmi_flags |= BIT_DISC_CENTER
    # Methods 4, 5 use disc edge (upper limb) which is the default

    if method in (2, 4):  # No refraction
        rsmi_flags |= BIT_NO_REFRACTION
    # Methods 3, 5 include refraction (default when flag not set)

    # Strategy: Find the next rise and next set after jd, then find the
    # previous rise and set by searching ~1.5 days before those.
    # Compare which previous event is more recent to determine if
    # planet is above or below horizon.

    try:
        # Find next rise after jd
        retflag, tret = rise_trans(
            jd,
            planet,
            CALC_RISE | rsmi_flags,
            [lon, lat, altitude],
            pressure,
            temperature,
            flags,
        )
        if retflag == -2:
            return fallback()  # Circumpolar
        jd_next_rise = tret[0]

        # Find next set after jd
        retflag, tret = rise_trans(
            jd,
            planet,
            CALC_SET | rsmi_flags,
            [lon, lat, altitude],
            pressure,
            temperature,
            flags,
        )
        if retflag == -2:
            return fallback()  # Circumpolar
        jd_next_set = tret[0]

        # Find previous rise before jd by searching before the next rise
        # (planets rise roughly once per day, so ~1.5 days before next rise
        # should find the previous rise)
        retflag, tret = rise_trans(
            jd_next_rise - 1.5,
            planet,
            CALC_RISE | rsmi_flags,
            [lon, lat, altitude],
            pressure,
            temperature,
            flags,
        )
        if retflag == -2:
            return fallback()
        jd_prev_rise = tret[0]
        # Verify this rise is actually before jd
        if jd_prev_rise >= jd:
            # Try searching earlier
            retflag, tret = rise_trans(
                jd_next_rise - 2.5,
                planet,
                CALC_RISE | rsmi_flags,
                [lon, lat, altitude],
                pressure,
                temperature,
                flags,
            )
            jd_prev_rise = tret[0]
            if retflag == -2 or jd_prev_rise >= jd:
                return fallback()

        # Find previous set before jd by searching before the next set
        retflag, tret = rise_trans(
            jd_next_set - 1.5,
            planet,
            CALC_SET | rsmi_flags,
            [lon, lat, altitude],
            pressure,
            temperature,
            flags,
        )
        if retflag == -2:
            return fallback()
        jd_prev_set = tret[0]
        # Verify this set is actually before jd
        if jd_prev_set >= jd:
            # Try searching earlier
            retflag, tret = rise_trans(
                jd_next_set - 2.5,
                planet,
                CALC_SET | rsmi_flags,
                [lon, lat, altitude],
                pressure,
                temperature,
                flags,
            )
            jd_prev_set = tret[0]
            if retflag == -2 or jd_prev_set >= jd:
                return fallback()

    except (IndexError, TypeError, ValueError, ArithmeticError):
        return fallback()

    # Determine if planet is above or below horizon:
    # If the most recent event before jd was a rise, planet is above horizon.
    # If the most recent event before jd was a set, planet is below horizon.
    if jd_prev_rise > jd_prev_set:
        # Planet rose more recently than it set -> above horizon (sectors 1-18)
        # Diurnal arc runs from jd_prev_rise to jd_next_set
        diurnal_arc = jd_next_set - jd_prev_rise
        if diurnal_arc <= 0:
            return fallback()

        elapsed = jd - jd_prev_rise
        fraction = elapsed / diurnal_arc
        # Clamp fraction to [0, 1] for safety
        fraction = max(0.0, min(1.0, fraction))

        # Sector 1 is at rise, sector 10 is at MC, sector 18+ is near set
        sector = 1.0 + fraction * 18.0
    else:
        # Planet set more recently than it rose -> below horizon (sectors 19-36)
        # Nocturnal arc runs from jd_prev_set to jd_next_rise
        nocturnal_arc = jd_next_rise - jd_prev_set
        if nocturnal_arc <= 0:
            return fallback()

        elapsed = jd - jd_prev_set
        fraction = elapsed / nocturnal_arc
        # Clamp fraction to [0, 1] for safety
        fraction = max(0.0, min(1.0, fraction))

        # Sector 19 is at set, sector 28 is at IC, sector 36 is near next rise
        sector = 19.0 + fraction * 18.0

    # Normalize to range [1, 37)
    if sector >= 37.0:
        sector -= 36.0
    elif sector < 1.0:
        sector += 36.0

    return sector


def _gauquelin_sector_pythonic(
    jd: float,
    planet: "int | str",
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    flags: int = 0,
    method: int = 0,
) -> float:
    """
    Calculate the Gauquelin sector (1-36) in which a planet is located.

    Gauquelin sectors are a 36-fold division of the diurnal and nocturnal arcs
    used in statistical astrology research by Michel Gauquelin. Sectors are
    numbered clockwise from the Ascendant:
    - Sector 1: Rising (Ascendant)
    - Sector 10: Upper culmination (MC)
    - Sector 19: Setting (Descendant)
    - Sector 28: Lower culmination (IC)

    Reference API compatible function (gauquelin_sector equivalent).

    Args:
        jd: Julian Day in Universal Time (UT)
        planet: Planet/body ID (int, e.g. SUN, MOON) or star name (str)
        lat: Geographic latitude in degrees (positive North)
        lon: Geographic longitude in degrees (positive East)
        altitude: Geographic altitude in meters above sea level (default 0.0)
        pressure: Atmospheric pressure in mbar (default 1013.25)
        temperature: Atmospheric temperature in degrees Celsius (default 15.0)
        flags: Calculation flags (FLG_SWIEPH, FLG_TOPOCTR, etc.)
            Note: The gauquelin_sector() wrapper defaults to
            FLG_SWIEPH | FLG_TOPOCTR (32770) per pyswisseph.
        method: Calculation method:
            - 0: with latitude
            - 1: without latitude
            - 2: from rising/setting times of disc center
            - 3: from rising/setting times of disc center, with refraction
            - 4: from rising/setting times of disc edge
            - 5: from rising/setting times of disc edge, with refraction

    Returns:
        Sector position as float in range [1, 37).
        Integer part is sector number (1-36), decimal part is position within sector.

    Note:
        Methods 2-5 use actual rise/set times calculated via rise_trans().
        These methods may be slower than methods 0-1 but provide more
        accurate sector positions based on real rise/set events.

        For circumpolar objects (that never rise or set at the given
        latitude), methods 2-5 fall back to method 0.

    Example:
        >>> sector = _gauquelin_sector_pythonic(2451545.0, MARS, 48.85, 2.35)
        >>> print(f"Mars is in sector {int(sector)}")
    """
    # Methods 2-5: Use actual rise/set times
    if method in (2, 3, 4, 5):
        if isinstance(planet, str):
            raise NotImplementedError(
                "Star names are not supported with methods 2-5 for _gauquelin_sector_pythonic"
            )
        return _gauquelin_sector_from_rise_set(
            jd, planet, lat, lon, altitude, pressure, temperature, flags, method
        )

    # Methods 0-1: Use _house_pos_pythonic with Gauquelin house system ('G')
    # This computes Gauquelin sectors from house position
    from .planets import calc_ut
    from .fixed_stars import fixstar_ut
    from .cache import get_true_obliquity

    ts = get_timescale()
    t = ts.ut1_jd(jd)

    # Calculate obliquity of ecliptic
    eps = get_true_obliquity(t.tt)

    # Calculate ARMC (sidereal time at location)
    gast = float(t.gast)  # in hours
    armc_deg = (gast * 15.0 + lon) % 360.0

    # Get body position - planet (int) or star (str).
    # The reference API computes the body topocentrically from the geopos
    # argument (default flags include FLG_TOPOCTR): set the observer for
    # the sub-call and restore the caller's topo afterwards.
    from . import state as _state
    from .state import get_topo, set_topo

    saved_topo = get_topo()
    try:
        if flags & FLG_TOPOCTR:
            set_topo(lon, lat, altitude)
        if isinstance(planet, str):
            pos, _star_name, retflag = fixstar_ut(planet, jd, flags | FLG_SPEED)
            planet_lon = float(pos[0])
            planet_lat = float(pos[1])
        else:
            pos, retflag = calc_ut(jd, planet, flags | FLG_SPEED)
            planet_lon = pos[0]
            planet_lat = pos[1]
    finally:
        if flags & FLG_TOPOCTR:
            # Restore the exact previous observer object (or None) so the
            # caller's topo state and observer-cache identity are untouched.
            with _state._STATE_LOCK:
                _state._TOPO = saved_topo

    # For method 1 (without latitude), ignore ecliptic latitude
    if method == 1:
        planet_lat = 0.0

    # Use _house_pos_pythonic with 'G' (Gauquelin sectors) to get the sector position
    # _house_pos_pythonic returns values in [1, 37) for Gauquelin
    sector = _house_pos_pythonic(armc_deg, lat, eps, (planet_lon, planet_lat), "G")

    return sector


def gauquelin_sector(
    tjdut: float,
    body: "int | str",
    method: int,
    geopos: tuple[float, float, float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    flags: int = FLG_SWIEPH | FLG_TOPOCTR,
) -> float:
    """
    Calculate the Gauquelin sector (1-36) in which a planet is located.

    Gauquelin sectors are a 36-fold division of the diurnal and nocturnal arcs
    used in statistical astrology research by Michel Gauquelin. Sectors are
    numbered clockwise from the Ascendant:
    - Sector 1: Rising (Ascendant)
    - Sector 10: Upper culmination (MC)
    - Sector 19: Setting (Descendant)
    - Sector 28: Lower culmination (IC)

    Reference API compatible function (gauquelin_sector equivalent).

    Args:
        jd: Julian Day in Universal Time (UT)
        planet: Planet/body ID (SUN, MOON, etc.)
        lat: Geographic latitude in degrees (positive North)
        lon: Geographic longitude in degrees (positive East)
        altitude: Geographic altitude in meters above sea level (default 0.0)
        pressure: Atmospheric pressure in mbar (default 1013.25)
        temperature: Atmospheric temperature in degrees Celsius (default 15.0)
        flags: Calculation flags (FLG_SWIEPH, FLG_TOPOCTR, etc.)
            Note: The gauquelin_sector() wrapper defaults to
            FLG_SWIEPH | FLG_TOPOCTR (32770) per pyswisseph.
        method: Calculation method:
            - 0: with latitude
            - 1: without latitude
            - 2: from rising/setting times of disc center of planet
            - 3: from rising/setting times of disc center, incl. refraction
            - 4: from rising/setting times of disk edge of planet
            - 5: from rising/setting times of disk edge, incl. refraction
        geopos: Tuple of (longitude, latitude, altitude):
            - longitude: geographic longitude in degrees (eastern positive)
            - latitude: geographic latitude in degrees (northern positive)
            - altitude: geographic altitude in meters above sea level
        atpress: Atmospheric pressure in mbar (if 0, default 1013.25 is used)
        attemp: Atmospheric temperature in degrees Celsius (if 0, default 15 is used)
        flags: Bit flags for ephemeris (FLG_SWIEPH, etc.)

    Returns:
        Sector position as float in range [1, 37).
        Integer part is sector number (1-36), decimal part is position within sector.

    Note:
        This function matches the gauquelin_sector() API signature.
        Sectors are numbered clockwise from the Ascendant:
        - Sector 1: Rising (Ascendant)
        - Sector 10: Upper culmination (MC)
        - Sector 19: Setting (Descendant)
        - Sector 28: Lower culmination (IC)

    Example:
        >>> geopos = (2.35, 48.85, 0.0)  # Paris: (lon, lat, alt)
        >>> sector = gauquelin_sector(2451545.0, MARS, 0, geopos)
    """
    lon, lat, altitude = geopos
    # Use defaults if 0 is passed (standard API behavior)
    pressure = atpress if atpress != 0.0 else 1013.25
    temperature = attemp if attemp != 0.0 else 15.0

    return _gauquelin_sector_pythonic(
        tjdut, body, lat, lon, altitude, pressure, temperature, flags, method
    )
