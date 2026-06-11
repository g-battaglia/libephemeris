"""
Solar and Lunar eclipse calculations for libephemeris.

Finds eclipse events and calculates their circumstances.

Functions:
- _sol_eclipse_when_glob_pythonic: Find next global solar eclipse
- _sol_eclipse_when_loc_pythonic: Find eclipse at specific location
- _sol_eclipse_where_pythonic: Calculate path of eclipse
- _sol_eclipse_how_pythonic: Eclipse circumstances at location
- _lun_eclipse_when_pythonic: Find next lunar eclipse

Algorithm:
    Solar eclipses occur at New Moon when Moon is near the ecliptic plane.
    1. Find next New Moon (Sun-Moon conjunction in longitude)
    2. Check lunar latitude - if |lat| < ~1.5° eclipse is possible
    3. Calculate eclipse magnitude and type based on distances
    4. Use Besselian elements for precise timing of phases

    Lunar eclipses occur at Full Moon when Moon is near a lunar node.
    1. Find next Full Moon (Sun-Moon opposition in longitude)
    2. Check Moon's distance from node - if close, eclipse is possible
    3. Calculate eclipse type (penumbral, partial, total) based on geometry
    4. Calculate phase times based on shadow cone sizes

References:
    - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
    - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
    - Reference documentation
"""

from __future__ import annotations

import math
from typing import Sequence, Tuple, Union

import numpy as np
from skyfield.errors import EphemerisRangeError
from .constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    FLG_SPEED,
    FLG_SWIEPH,
    FLG_EQUATORIAL,
    ECL_TOTAL,
    ECL_ANNULAR,
    ECL_PARTIAL,
    ECL_ANNULAR_TOTAL,
    ECL_CENTRAL,
    ECL_NONCENTRAL,
    ECL_ALLTYPES_SOLAR,
    ECL_ALLTYPES_LUNAR,
    ECL_PENUMBRAL,
    ECL_GRAZING,
    ECL_VISIBLE,
    ECL_MAX_VISIBLE,
    ECL_1ST_VISIBLE,
    ECL_2ND_VISIBLE,
    ECL_3RD_VISIBLE,
    ECL_4TH_VISIBLE,
)
from .planets import calc_ut
from .state import get_timescale


def _get_leb_reader_safe():
    """Return LEB reader or None. Wraps get_leb_reader for eclipse dispatch."""
    from .state import get_leb_reader

    return get_leb_reader()


def _is_leb_out_of_range(exc: BaseException) -> bool:
    """Return True if the exception is an LEB out-of-range error.

    Recognises the ValueError raised by LEBReader / LEB2Reader / fast_calc
    when a JD is probed outside the LEB file's coverage. Both bodies and
    nutation coverage produce messages containing "outside range" or
    "outside nutation range".
    """
    return isinstance(exc, (KeyError, ValueError)) and (
        "outside" in str(exc).lower()
    )


def _call_with_leb_skyfield_fallback(impl, *args, **kwargs):
    """Run ``impl(reader=...)``, retrying with ``reader=None`` on LEB out-of-range.

    Eclipse / occultation / rise-transit search loops define LEB-backed inner
    closures inside ``if reader is not None`` blocks. If a search probe lands
    outside the LEB file's coverage, _topo_ecliptic raises ValueError. Without
    this helper there was no fallback, so partial/custom LEB files shorter
    than the search window would surface the underlying ValueError to callers.

    Behaviour:
        - If no LEB reader is available, run ``impl(reader=None)`` directly.
        - Otherwise, run ``impl(reader=reader)`` first.
        - If it raises an LEB out-of-range error, retry with ``reader=None``
          to force the Skyfield path. Any other exception is re-raised.

    Standard LEB files (base 1849-2150, medium 1550-2650, extended ±13000) all
    accommodate the search windows in this module, so the fallback only fires
    for custom files. The retry restarts computation from scratch (no shared
    intermediate state), which is acceptable for the rare-fallback case.
    """
    reader = _get_leb_reader_safe()
    if reader is not None:
        try:
            return impl(*args, reader=reader, **kwargs)
        except (KeyError, ValueError) as exc:
            if not _is_leb_out_of_range(exc):
                raise
    return impl(*args, reader=None, **kwargs)


def _reraise_if_leb_range_error(exc: BaseException) -> None:
    """Re-raise ``exc`` when it signals LEB out-of-range coverage.

    Defensive ``except`` blocks inside eclipse implementations must not
    swallow LEB coverage errors: doing so converts an out-of-range probe
    into fabricated zeros/sentinels and prevents
    ``_call_with_leb_skyfield_fallback`` from retrying on the Skyfield
    path. Call this first in any broad except that would otherwise
    degrade gracefully.
    """
    if _is_leb_out_of_range(exc):
        raise exc


# Constants for eclipse calculations
SYNODIC_MONTH = 29.530588853  # Mean synodic month in days
LUNAR_NODE_PERIOD = 6798.38  # Lunar node regression period in days
ECLIPSE_LIMIT_SOLAR = 18.5  # Maximum elongation from node for solar eclipse (degrees)
EARTH_RADIUS_KM = 6378.137  # Earth equatorial radius in km (WGS84)

# Edge case thresholds for eclipse calculations
# These constants define limits for shallow/near-miss eclipses
SHALLOW_ECLIPSE_MAG_THRESHOLD = 0.01  # Minimum magnitude for reliable contact times
NEAR_MISS_GAMMA_MARGIN = 0.02  # Margin from gamma limit for edge case handling
MINIMUM_SEPARATION_FOR_LENS = 1e-10  # Minimum separation to avoid division by zero

# 1 AU in km (IAU 2012 definition)
_AU_KM = 149597870.7

# Physical equatorial radii of planets in km
# Sources: IAU 2015 nominal values, NASA Planetary Fact Sheet
# Used for dynamic angular radius computation (matching SE's pla_diam[] approach)
_PLANET_RADIUS_KM = {
    MERCURY: 2439.7,
    VENUS: 6051.8,
    MARS: 3396.2,
    JUPITER: 71492.0,
    SATURN: 60268.0,  # Disc only, excludes rings
    URANUS: 25559.0,
    NEPTUNE: 24764.0,
    PLUTO: 1188.3,
}


def _calc_planet_angular_radius(planet_id: int, dist_au: float) -> float:
    """Compute angular radius of a planet in degrees from its geocentric distance.

    Uses the exact geometric formula: angular_radius = arcsin(physical_radius / distance),
    matching the Swiss Ephemeris approach of computing angular sizes dynamically from
    physical radii and actual distance rather than using static lookup tables.

    Args:
        planet_id: Planet constant (MERCURY, VENUS, etc.)
        dist_au: Geocentric distance in AU.

    Returns:
        Angular radius in degrees. Returns 0.0001° (~0.4") for unknown bodies
        (point source approximation).
    """
    radius_km = _PLANET_RADIUS_KM.get(planet_id)
    if radius_km is None:
        return 0.0001  # Point source for unknown bodies
    dist_km = dist_au * _AU_KM
    if dist_km <= 0:
        return 0.0001
    return math.degrees(math.asin(min(1.0, radius_km / dist_km)))


def _is_shallow_eclipse(magnitude: float) -> bool:
    """
    Check if an eclipse is very shallow (magnitude close to zero).

    Very shallow eclipses have unreliable contact time calculations because
    the penumbra barely grazes Earth or the Moon barely enters the shadow.

    Args:
        magnitude: Eclipse magnitude (0 to ~1.5 for solar, 0 to ~2 for lunar)

    Returns:
        True if the eclipse is considered shallow (magnitude < threshold)
    """
    return magnitude < SHALLOW_ECLIPSE_MAG_THRESHOLD


def _is_near_miss_eclipse(gamma: float, gamma_limit: float = 1.55) -> bool:
    """
    Check if an eclipse is a near-miss (gamma very close to the eclipse limit).

    Near-miss eclipses have gamma values very close to the limit where eclipses
    cease to occur. These require special handling to avoid numerical instability.

    Args:
        gamma: Eclipse gamma parameter (distance of shadow axis from Earth center)
        gamma_limit: Maximum gamma for any eclipse visibility (default 1.55)

    Returns:
        True if gamma is within NEAR_MISS_GAMMA_MARGIN of the limit
    """
    return abs(gamma) > (gamma_limit - NEAR_MISS_GAMMA_MARGIN)


def _safe_acos(x: float) -> float:
    """
    Safe arccosine that clamps input to valid range [-1, 1].

    Numerical errors can cause values slightly outside the valid range,
    which would raise a domain error. This function clamps the input.

    Args:
        x: Input value (should be in [-1, 1])

    Returns:
        math.acos of the clamped value
    """
    return math.acos(max(-1.0, min(1.0, x)))


def _safe_sqrt(x: float) -> float:
    """
    Safe square root that returns 0 for negative inputs.

    Numerical errors can cause small negative values where zero is expected.
    This function returns 0 for negative inputs instead of raising an error.

    Args:
        x: Input value (should be >= 0)

    Returns:
        math.sqrt of x if x >= 0, else 0.0
    """
    return math.sqrt(max(0.0, x))


def _calculate_obscuration_safe(r_sun: float, r_moon: float, d: float) -> float:
    """
    Calculate eclipse obscuration with edge case handling.

    Obscuration is the fraction of the Sun's area covered by the Moon.
    This function handles edge cases like:
    - Zero separation (total/annular eclipse)
    - Separation equal to sum of radii (no eclipse)
    - Very small separations that could cause numerical issues

    Args:
        r_sun: Sun's angular radius in degrees
        r_moon: Moon's angular radius in degrees
        d: Center-to-center separation in degrees

    Returns:
        Obscuration as a fraction. For total eclipses where Moon is larger
        than Sun, this can exceed 1.0 (matching reference API behavior where
        obscuration = (r_moon/r_sun)^2 represents the disc area ratio).
    """
    # Handle edge case: no overlap
    if d >= r_sun + r_moon:
        return 0.0

    # Handle edge case: one disk entirely within the other
    if d <= abs(r_sun - r_moon):
        # Both total and annular: obscuration = disc area ratio
        # For total eclipses (r_moon >= r_sun), this is >= 1.0
        # For annular eclipses (r_moon < r_sun), this is < 1.0
        return (r_moon / r_sun) ** 2

    # Handle edge case: zero or near-zero separation
    if d < MINIMUM_SEPARATION_FOR_LENS:
        return (r_moon / r_sun) ** 2

    # Calculate lens-shaped intersection area
    # d1 is distance from Sun center to intersection chord
    # d2 is distance from Moon center to intersection chord
    d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
    d2 = d - d1

    # Validate that the geometric configuration is valid
    if abs(d1) > r_sun or abs(d2) > r_moon:
        # This shouldn't happen if we passed earlier checks, but handle gracefully
        return 0.0

    # Calculate areas using lens formula with safe operations
    area1 = r_sun * r_sun * _safe_acos(d1 / r_sun) - d1 * _safe_sqrt(
        r_sun * r_sun - d1 * d1
    )
    area2 = r_moon * r_moon * _safe_acos(d2 / r_moon) - d2 * _safe_sqrt(
        r_moon * r_moon - d2 * d2
    )

    intersection_area = area1 + area2
    sun_area = math.pi * r_sun * r_sun

    obscuration = intersection_area / sun_area
    return max(0.0, obscuration)


def _calculate_magnitude_safe(
    gamma: float, moon_sun_ratio: float, gamma_limit_partial: float = 1.55
) -> float:
    """
    Calculate eclipse magnitude with edge case handling.

    Handles shallow partial eclipses where gamma is close to the eclipse limit.

    Args:
        gamma: Eclipse gamma parameter (shadow axis distance from Earth center)
        moon_sun_ratio: Ratio of Moon's apparent diameter to Sun's
        gamma_limit_partial: Maximum gamma for partial eclipse visibility

    Returns:
        Eclipse magnitude (0 to ~1.5), clamped to valid range
    """
    # For very shallow eclipses (gamma near limit), magnitude approaches 0
    if _is_near_miss_eclipse(gamma, gamma_limit_partial):
        # Use linear interpolation near the edge for smooth transition
        remaining = gamma_limit_partial - abs(gamma)
        if remaining <= 0:
            return 0.0
        # Magnitude decreases linearly as gamma approaches limit
        magnitude = (remaining / NEAR_MISS_GAMMA_MARGIN) * SHALLOW_ECLIPSE_MAG_THRESHOLD
        return max(0.0, min(1.5, magnitude * moon_sun_ratio))

    # Standard magnitude calculation
    if abs(gamma) >= gamma_limit_partial:
        return 0.0

    magnitude = 1.0 - abs(gamma) / gamma_limit_partial
    magnitude = magnitude * moon_sun_ratio

    return max(0.0, min(1.5, magnitude))


def _validate_contact_time(
    jd_contact: float, jd_max: float, max_offset_days: float = 0.25
) -> float:
    """
    Validate a calculated contact time and return 0 if invalid.

    Contact times should be within a reasonable range of the eclipse maximum.
    This function checks if the contact time is valid.

    Args:
        jd_contact: Calculated contact time (Julian Day)
        jd_max: Eclipse maximum time (Julian Day)
        max_offset_days: Maximum allowed offset from maximum (default 6 hours)

    Returns:
        The contact time if valid, 0.0 if invalid
    """
    if jd_contact <= 0:
        return 0.0

    offset = abs(jd_contact - jd_max)
    if offset > max_offset_days:
        return 0.0

    return jd_contact


def _find_next_new_moon(jd_start: float) -> float:
    """
    Find the next New Moon (Sun-Moon conjunction) after jd_start.

    Uses iterative refinement to find exact moment of conjunction.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of next New Moon
    """
    # Get current positions
    sun_pos, _ = calc_ut(jd_start, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_start, MOON, FLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation (Moon - Sun), normalized to -180 to 180
    elongation = (moon_lon - sun_lon) % 360.0
    if elongation > 180:
        elongation -= 360

    # Estimate time to next conjunction
    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day (average Moon speed - Sun speed)

    # Time until next conjunction (elongation = 0)
    if elongation > 0:
        # Moon is ahead, wait for next cycle
        dt = (360.0 - elongation) / relative_speed
    else:
        # Moon is behind Sun
        dt = (-elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = calc_ut(jd_guess, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd_guess, MOON, FLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation
        diff = (moon_lon - sun_lon) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _find_previous_new_moon(jd_start: float) -> float:
    """
    Find the previous New Moon (Sun-Moon conjunction) before jd_start.

    Uses iterative refinement to find exact moment of conjunction.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of previous New Moon
    """
    # Get current positions
    sun_pos, _ = calc_ut(jd_start, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_start, MOON, FLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation (Moon - Sun), normalized to -180 to 180
    elongation = (moon_lon - sun_lon) % 360.0
    if elongation > 180:
        elongation -= 360

    # Estimate time to previous conjunction
    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day (average Moon speed - Sun speed)

    # Time since last conjunction (go backwards)
    if elongation > 0:
        # Moon is ahead of Sun, last conjunction was elongation/speed days ago
        dt = -elongation / relative_speed
    else:
        # Moon is behind Sun, last conjunction was (360 + elongation)/speed days ago
        dt = -(360.0 + elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = calc_ut(jd_guess, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd_guess, MOON, FLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation
        diff = (moon_lon - sun_lon) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _get_moon_node_distance(jd: float, moon_lon: float) -> float:
    """
    Calculate angular distance of Moon from nearest lunar node.

    Args:
        jd: Julian Day (UT)
        moon_lon: Moon's ecliptic longitude in degrees

    Returns:
        Absolute angular distance from nearest node in degrees
    """
    from .lunar import calc_mean_lunar_node

    # Get mean node longitude
    ts = get_timescale()
    t = ts.ut1_jd(jd)
    node_lon = calc_mean_lunar_node(t.tt)

    # Calculate distance to both nodes
    dist_to_north = abs((moon_lon - node_lon + 180) % 360 - 180)
    dist_to_south = abs((moon_lon - (node_lon + 180) + 180) % 360 - 180)

    return min(dist_to_north, dist_to_south)


def _calc_gamma(jd: float) -> float:
    """
    Calculate the gamma parameter (sqrt(x² + y²)) at given JD using Besselian elements.

    The gamma parameter is the minimum distance of the Moon's shadow axis from
    Earth's center, measured in Earth equatorial radii. During an eclipse:
    - gamma < 1: central eclipse (total or annular) is possible
    - gamma > 1.5: no eclipse visible from Earth

    Args:
        jd: Julian Day (UT)

    Returns:
        Gamma value in Earth equatorial radii
    """
    from .state import get_leb_reader

    AU_TO_KM = 149597870.7

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            from .fast_calc import _apparent_icrs_cartesian
            from .time_utils import deltat

            jd_tt = jd + deltat(jd)
            sun_pos = _apparent_icrs_cartesian(reader, jd_tt, SUN)
            moon_pos = _apparent_icrs_cartesian(reader, jd_tt, MOON)
            _leb_ok = True
        except (KeyError, ValueError):
            pass
    if not _leb_ok:
        from .state import get_planets, get_timescale

        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)
        earth_at = eph["earth"].at(t)
        sun_pos = earth_at.observe(eph["sun"]).apparent().position.au
        moon_pos = earth_at.observe(eph["moon"]).apparent().position.au

    # Shadow axis direction (Sun to Moon)
    shadow_dir = [
        moon_pos[0] - sun_pos[0],
        moon_pos[1] - sun_pos[1],
        moon_pos[2] - sun_pos[2],
    ]
    shadow_len = math.sqrt(shadow_dir[0] ** 2 + shadow_dir[1] ** 2 + shadow_dir[2] ** 2)
    shadow_unit = [shadow_dir[i] / shadow_len for i in range(3)]

    # x-axis perpendicular to shadow, in equatorial plane
    x_axis_raw = [-shadow_unit[1], shadow_unit[0], 0.0]
    x_axis_len = math.sqrt(x_axis_raw[0] ** 2 + x_axis_raw[1] ** 2)

    if x_axis_len < 1e-10:
        x_axis = [1.0, 0.0, 0.0]
    else:
        x_axis = [x_axis_raw[0] / x_axis_len, x_axis_raw[1] / x_axis_len, 0.0]

    # y-axis completes right-handed system
    y_axis = [
        shadow_unit[1] * x_axis[2] - shadow_unit[2] * x_axis[1],
        shadow_unit[2] * x_axis[0] - shadow_unit[0] * x_axis[2],
        shadow_unit[0] * x_axis[1] - shadow_unit[1] * x_axis[0],
    ]

    # Project Moon position onto fundamental plane
    moon_along_axis = sum(moon_pos[i] * shadow_unit[i] for i in range(3))
    moon_perp = [moon_pos[i] - moon_along_axis * shadow_unit[i] for i in range(3)]

    x_au = sum(moon_perp[i] * x_axis[i] for i in range(3))
    y_au = sum(moon_perp[i] * y_axis[i] for i in range(3))

    # Convert to Earth radii
    x_earth = (x_au * AU_TO_KM) / EARTH_RADIUS_KM
    y_earth = (y_au * AU_TO_KM) / EARTH_RADIUS_KM

    return math.sqrt(x_earth**2 + y_earth**2)


def _refine_solar_eclipse_maximum(
    jd_approx: float, search_range: float = 0.125
) -> float:
    """
    Refine the solar eclipse maximum time using Besselian elements.

    Uses golden section search to find the time when gamma (shadow axis distance
    from Earth center) is minimized. This gives sub-second precision for eclipse
    maximum timing.

    Args:
        jd_approx: Approximate Julian Day of eclipse maximum (e.g., from New Moon)
        search_range: Search range in days (default ±3 hours)

    Returns:
        Julian Day of refined eclipse maximum (precision < 1 second)
    """
    phi = (1 + math.sqrt(5)) / 2  # Golden ratio

    jd_low = jd_approx - search_range
    jd_high = jd_approx + search_range

    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    gamma_a = _calc_gamma(jd_a)
    gamma_b = _calc_gamma(jd_b)

    # Golden section search for minimum gamma
    for _ in range(60):  # Converges to ~1e-9 days (~0.1 ms)
        if gamma_a < gamma_b:
            jd_high = jd_b
            jd_b = jd_a
            gamma_b = gamma_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            gamma_a = _calc_gamma(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            gamma_a = gamma_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            gamma_b = _calc_gamma(jd_b)

        if jd_high - jd_low < 1e-8:  # ~0.86 ms precision
            break

    return (jd_low + jd_high) / 2


def _calc_penumbra_limit(jd: float) -> float:
    """
    Calculate l1 (penumbral shadow radius on fundamental plane) at given JD.

    Returns the radius in Earth radii where the penumbral shadow intersects
    the fundamental plane.
    """
    AU_TO_KM = 149597870.7
    SUN_RADIUS_KM = 696340.0
    MOON_RADIUS_KM = 1737.4

    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)
    sun_dist_au = sun_pos[2]
    moon_dist_au = moon_pos[2]

    sun_dist_km = sun_dist_au * AU_TO_KM
    moon_dist_km = moon_dist_au * AU_TO_KM

    # The penumbral shadow cone extends from the Moon toward Earth.
    # The penumbral half-angle f1 is defined such that the penumbral cone's
    # outer edge connects the far limb of the Sun to the near limb of the Moon.
    #
    # Using similar triangles:
    # tan(f1) = (SUN_RADIUS + MOON_RADIUS) / (sun_dist - moon_dist)
    # At the fundamental plane (distance moon_dist from Moon toward Earth),
    # the penumbra radius is: l1_km = moon_dist * tan(f1) + MOON_RADIUS

    # But for Besselian elements, l1 is measured from the shadow axis,
    # which is the penumbral radius at the fundamental plane.
    #
    # Standard formula for l1 (penumbral shadow radius at fundamental plane):
    # l1 = sin(f1) + z * tan(f1) where z is along shadow axis
    # At the fundamental plane, this simplifies to the geometry below.

    # Angular semi-diameter of Sun from Earth
    sun_angular_rad = SUN_RADIUS_KM / sun_dist_km  # radians

    # Angular semi-diameter of Moon from Earth
    moon_angular_rad = MOON_RADIUS_KM / moon_dist_km  # radians

    # Penumbral cone angular radius from shadow axis (sum of angular radii)
    f1_rad = sun_angular_rad + moon_angular_rad

    # At the fundamental plane, the penumbral radius in Earth radii:
    # The shadow axis passes through Earth's center at distance = 0
    # The penumbral limit at distance z from Moon is: r = z * tan(f1)
    # But we measure from shadow axis, so: l1 = moon_dist * f1 (for small angles)

    # More accurate: l1 = k1 * sec(f1) + z * tan(f1) where k1 is constant
    # For practical purposes with small angles:
    l1_km = moon_dist_km * f1_rad
    l1_earth_radii = l1_km / EARTH_RADIUS_KM

    return l1_earth_radii


def _calc_umbra_limit(jd: float) -> float:
    """
    Calculate l2 (umbral/antumbral shadow radius on fundamental plane) at given JD.

    Returns the radius in Earth radii. Negative for umbra (total eclipse),
    positive for antumbra (annular eclipse).
    """
    AU_TO_KM = 149597870.7
    SUN_RADIUS_KM = 696340.0
    MOON_RADIUS_KM = 1737.4

    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)
    sun_dist_au = sun_pos[2]
    moon_dist_au = moon_pos[2]

    sun_dist_km = sun_dist_au * AU_TO_KM
    moon_dist_km = moon_dist_au * AU_TO_KM

    # Umbral cone half-angle
    f2 = math.atan((SUN_RADIUS_KM - MOON_RADIUS_KM) / sun_dist_km)

    # Distance from Moon to umbra vertex
    umbra_vertex_km = MOON_RADIUS_KM / math.tan(f2) if f2 > 1e-10 else float("inf")

    # Umbral radius at fundamental plane
    if moon_dist_km < umbra_vertex_km:
        # Umbra reaches Earth - negative by convention
        l2_km = MOON_RADIUS_KM - moon_dist_km * math.tan(f2)
        l2_earth_radii = -l2_km / EARTH_RADIUS_KM
    else:
        # Antumbra - Moon shadow vertex is before Earth
        l2_km = moon_dist_km * math.tan(f2) - MOON_RADIUS_KM
        l2_earth_radii = l2_km / EARTH_RADIUS_KM

    return l2_earth_radii


def _find_contact_time_besselian(
    jd_max: float,
    target_gamma: float,
    search_before: bool,
    search_range: float = 0.1,
) -> float:
    """
    Find the time when gamma equals target_gamma using bisection.

    Used to find eclipse contact times where the shadow boundary crosses Earth.

    Args:
        jd_max: Julian Day of eclipse maximum
        target_gamma: Target gamma value to search for
        search_before: If True, search before maximum; if False, search after
        search_range: Search range in days from maximum

    Returns:
        Julian Day when gamma equals target_gamma, or 0 if not found
    """
    if search_before:
        jd_start = jd_max - search_range
        jd_end = jd_max
    else:
        jd_start = jd_max
        jd_end = jd_max + search_range

    # Check if solution exists in range
    gamma_start = _calc_gamma(jd_start)
    gamma_end = _calc_gamma(jd_end)
    gamma_max = _calc_gamma(jd_max)

    # For contacts before max: gamma decreases from start to max, then increases
    # For contacts after max: gamma increases from max to end
    # We're looking for where gamma crosses target_gamma

    if search_before:
        if gamma_start < target_gamma or gamma_max > target_gamma:
            # Check edge case where max gamma might be above target
            if gamma_max < target_gamma:
                return 0.0
    else:
        if gamma_max > target_gamma or gamma_end < target_gamma:
            if gamma_max < target_gamma:
                return 0.0

    # Binary search
    for _ in range(60):  # ~0.1 second precision
        jd_mid = (jd_start + jd_end) / 2
        gamma_mid = _calc_gamma(jd_mid)

        if abs(gamma_mid - target_gamma) < 1e-8:
            return jd_mid

        if search_before:
            # Before max: gamma decreases then potentially rises
            # We want the first crossing from above
            if gamma_mid > target_gamma:
                jd_start = jd_mid
            else:
                jd_end = jd_mid
        else:
            # After max: gamma increases
            if gamma_mid < target_gamma:
                jd_start = jd_mid
            else:
                jd_end = jd_mid

        if jd_end - jd_start < 1e-8:
            break

    return (jd_start + jd_end) / 2


def _calculate_eclipse_phases_besselian(
    jd_max: float, eclipse_type: int
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate times of eclipse phases using Besselian elements for high precision.

    This function achieves timing precision of better than 10 seconds by
    calculating exact contact times based on when the penumbral and umbral
    shadow boundaries cross Earth's limb.

    Phase indices (matching reference API tret array):
        [0]: Time of maximum eclipse
        [1]: Time when eclipse takes place at local apparent noon
        [2]: Eclipse begins (first penumbral contact with Earth, P1)
        [3]: Eclipse ends (last penumbral contact with Earth, P4)
        [4]: Totality/annularity begins anywhere on Earth (U1)
        [5]: Totality/annularity ends (U4)
        [6]: Center line begins
        [7]: Center line ends
        [8]: Annular-total eclipse becomes total
        [9]: Annular-total eclipse becomes annular again

    Args:
        jd_max: Julian Day of eclipse maximum (refined using Besselian elements)
        eclipse_type: Eclipse type flags

    Returns:
        Tuple of 10 floats with phase times (JD UT), matching reference API format
    """
    is_central = bool(eclipse_type & ECL_CENTRAL)
    is_total = bool(eclipse_type & ECL_TOTAL)
    is_annular = bool(eclipse_type & ECL_ANNULAR)

    # Get l1 (penumbral limit) and l2 (umbral limit) at maximum
    l1 = _calc_penumbra_limit(jd_max)
    l2 = _calc_umbra_limit(jd_max)
    gamma_max = _calc_gamma(jd_max)

    # For global eclipse, first/fourth contacts occur when gamma = 1 + l1
    # (penumbral shadow touches Earth's limb from outside)
    penumbral_limit = 1.0 + l1  # Earth radius + penumbra radius

    # Calculate first contact (penumbra first touches Earth)
    t_first_contact = _find_contact_time_besselian(
        jd_max, penumbral_limit, search_before=True, search_range=0.15
    )

    # Calculate fourth contact (penumbra last leaves Earth)
    t_fourth_contact = _find_contact_time_besselian(
        jd_max, penumbral_limit, search_before=False, search_range=0.15
    )

    # For total/annular eclipses, calculate U1/U4 (totality begin/end
    # anywhere on Earth)
    t_second_contact = 0.0
    t_third_contact = 0.0

    if (is_total or is_annular) and abs(l2) > 0:
        # The umbra/antumbra first and last touches Earth at EXTERIOR
        # tangency: shadow-axis distance = 1 + |l2|. (The previous
        # interior criterion 1 - |l2| — umbra entirely inside the disc —
        # was ~2-4 minutes late/early systematically and returned zeros
        # for grazing-central eclipses.)
        umbral_limit = 1.0 + abs(l2)

        if gamma_max < umbral_limit:
            t_second_contact = _find_contact_time_besselian(
                jd_max, umbral_limit, search_before=True, search_range=0.10
            )
            t_third_contact = _find_contact_time_besselian(
                jd_max, umbral_limit, search_before=False, search_range=0.10
            )

    # tret[1]: time when the eclipse takes place at local apparent noon —
    # the shadow axis crosses the x = 0 plane (Besselian x sign change)
    # within the eclipse window.
    t_noon = 0.0
    t_lo = t_first_contact if t_first_contact else jd_max - 0.15
    t_hi = t_fourth_contact if t_fourth_contact else jd_max + 0.15
    try:
        x_lo = calc_besselian_x(t_lo)
        x_hi = calc_besselian_x(t_hi)
        if x_lo * x_hi < 0:
            a, b, fa = t_lo, t_hi, x_lo
            for _ in range(60):
                mid = 0.5 * (a + b)
                fm = calc_besselian_x(mid)
                if fa * fm <= 0:
                    b = mid
                else:
                    a, fa = mid, fm
            t_noon = 0.5 * (a + b)
    except (KeyError, ValueError, ArithmeticError) as _exc:
        _reraise_if_leb_range_error(_exc)
        t_noon = 0.0

    # tret[6]/[7]: center-line begin/end — the shadow AXIS touches the
    # fundamental-plane disc while gamma < 1.
    t_cl_begin = 0.0
    t_cl_end = 0.0
    if (is_total or is_annular) and gamma_max < 1.0:
        t_cl_begin = _find_contact_time_besselian(
            jd_max, 1.0, search_before=True, search_range=0.10
        )
        t_cl_end = _find_contact_time_besselian(
            jd_max, 1.0, search_before=False, search_range=0.10
        )

    return (
        jd_max,  # [0] Time of maximum eclipse
        t_noon,  # [1] Eclipse at local apparent noon
        # No fabricated fallbacks: an unresolved contact stays 0.0 (the
        # reference returns 0 rather than an invented +-1 h estimate)
        t_first_contact if t_first_contact else 0.0,  # [2] Eclipse begin
        t_fourth_contact if t_fourth_contact else 0.0,  # [3] Eclipse end
        t_second_contact,  # [4] Totality begin (or 0)
        t_third_contact,  # [5] Totality end (or 0)
        t_cl_begin,  # [6] Center line begin
        t_cl_end,  # [7] Center line end
        0.0,  # [8] Annular-total becomes total
        0.0,  # [9] Annular-total becomes annular again
    )


def _calculate_eclipse_type_and_magnitude(
    jd: float,
) -> Tuple[int, float, float, float]:
    """
    Determine eclipse type and magnitude at maximum eclipse.

    Uses geometric calculations based on Sun-Moon-Earth distances
    and apparent angular sizes. Includes edge case handling for:
    - Very shallow partial eclipses (magnitude near 0)
    - Near-miss eclipses (gamma close to eclipse limit)
    - Grazing eclipses at extreme gamma values

    Args:
        jd: Julian Day of eclipse maximum (UT)

    Returns:
        Tuple of (eclipse_type_flags, magnitude, gamma, moon_sun_ratio)
        - eclipse_type_flags: Bitmask of ECL_* constants
        - magnitude: Eclipse magnitude (fraction of Sun's diameter covered)
        - gamma: Gamma parameter (distance of Moon shadow axis from Earth center)
        - moon_sun_ratio: Ratio of apparent Moon diameter to Sun diameter
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Distances in AU
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Ecliptic coordinates for separation calculation
    sun_lon, sun_lat = sun_pos[0], sun_pos[1]
    moon_lon, moon_lat = moon_pos[0], moon_pos[1]

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    sun_angular_radius = (959.63 / 3600.0) / sun_dist
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Ratio of apparent sizes - handle edge case of very small Sun radius
    if sun_angular_radius < 1e-10:
        moon_sun_ratio = 1.0  # Fallback for numerical stability
    else:
        moon_sun_ratio = moon_angular_radius / sun_angular_radius

    # Calculate gamma using Besselian elements (shadow axis distance from Earth center)
    gamma = _calc_gamma(jd)

    # Penumbral and umbral limits at maximum (Earth radii)
    l1 = _calc_penumbra_limit(jd)
    l2 = _calc_umbra_limit(jd)

    gamma_limit_partial = 1.0 + l1
    gamma_limit_central = max(0.0, 1.0 - abs(l2))

    # Angular separation between Sun and Moon centers (degrees)
    def _angular_separation(
        lon1: float, lat1: float, lon2: float, lat2: float
    ) -> float:
        lon1_rad = math.radians(lon1)
        lat1_rad = math.radians(lat1)
        lon2_rad = math.radians(lon2)
        lat2_rad = math.radians(lat2)
        dlon = lon2_rad - lon1_rad
        cos_sep = math.sin(lat1_rad) * math.sin(lat2_rad) + math.cos(
            lat1_rad
        ) * math.cos(lat2_rad) * math.cos(dlon)
        return math.degrees(_safe_acos(cos_sep))

    separation = _angular_separation(sun_lon, sun_lat, moon_lon, moon_lat)

    # Eclipse magnitude: fraction of Sun's diameter covered (can exceed 1.0 for total)
    sum_radii = sun_angular_radius + moon_angular_radius
    if separation >= sum_radii:
        magnitude = 0.0
    else:
        magnitude = (sum_radii - separation) / (2.0 * sun_angular_radius)
        magnitude = max(0.0, min(1.5, magnitude))

    # Determine eclipse type based on gamma and shadow limits
    eclipse_type = 0

    # Edge case: gamma beyond eclipse limit
    if abs(gamma) > gamma_limit_partial:
        # No eclipse
        return 0, 0.0, gamma, moon_sun_ratio

    # Edge case: near-miss eclipse (gamma very close to limit).
    # The reference returns PARTIAL|NONCENTRAL (retflag 18) for shallow
    # partials; the non-reference ECL_GRAZING bit is no longer leaked
    # into any public retflag.
    if _is_near_miss_eclipse(gamma, gamma_limit_partial):
        eclipse_type = ECL_PARTIAL | ECL_NONCENTRAL
        magnitude = _calculate_magnitude_safe(
            gamma, moon_sun_ratio, gamma_limit_partial
        )
        return eclipse_type, magnitude, gamma, moon_sun_ratio

    eclipse_type = ECL_PARTIAL | ECL_NONCENTRAL

    # Hybrid detection: the umbra/antumbra cone tip is very close to
    # Earth's surface, so the eclipse transitions between total and
    # annular along its path (|l2| small enough for a sign change).
    is_hybrid = abs(l2) < 0.002  # ~13 km, empirical threshold

    if abs(gamma) <= 1.0:
        # Central: the shadow axis intersects Earth.
        # l2 < 0: umbra (total); l2 > 0: antumbra (annular).
        if is_hybrid:
            eclipse_type = ECL_ANNULAR_TOTAL | ECL_CENTRAL
        elif l2 < 0:
            eclipse_type = ECL_TOTAL | ECL_CENTRAL
        else:
            eclipse_type = ECL_ANNULAR | ECL_CENTRAL
    elif abs(gamma) <= 1.0 + abs(l2) and abs(l2) > 0:
        # Noncentral total/annular (~1% of eclipses): the umbra/antumbra
        # touches Earth but the shadow axis misses it. The reference
        # classifies these as TOTAL|NONCENTRAL / ANNULAR|NONCENTRAL.
        if is_hybrid:
            eclipse_type = ECL_ANNULAR_TOTAL | ECL_NONCENTRAL
        elif l2 < 0:
            eclipse_type = ECL_TOTAL | ECL_NONCENTRAL
        else:
            eclipse_type = ECL_ANNULAR | ECL_NONCENTRAL

    return eclipse_type, magnitude, gamma, moon_sun_ratio


def _calculate_eclipse_phases(
    jd_max: float, eclipse_type: int
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate times of eclipse phases (contacts) for a global eclipse.

    Uses Besselian elements for high-precision timing (< 10 seconds).

    Phase indices (matching reference API tret array):
        [0]: Time of maximum eclipse
        [1]: Time when eclipse takes place at local apparent noon
        [2]: Eclipse begins (first penumbral contact with Earth, P1)
        [3]: Eclipse ends (last penumbral contact with Earth, P4)
        [4]: Totality/annularity begins anywhere on Earth (U1)
        [5]: Totality/annularity ends (U4)
        [6]: Center line begins
        [7]: Center line ends
        [8]: Annular-total eclipse becomes total
        [9]: Annular-total eclipse becomes annular again

    Args:
        jd_max: Julian Day of maximum eclipse
        eclipse_type: Eclipse type flags

    Returns:
        Tuple of 10 floats with phase times (JD UT), matching reference API format
    """
    # Use the high-precision Besselian element-based calculation
    return _calculate_eclipse_phases_besselian(jd_max, eclipse_type)


def sol_eclipse_max_time(
    jd_approx: float,
    lat: float | None = None,
    lon: float | None = None,
    altitude: float = 0.0,
    search_range: float = 0.125,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, float]:
    """
    Calculate the precise time of maximum eclipse when Sun-Moon separation is minimum.

    This function finds the exact moment when the angular separation between the
    Sun and Moon centers is at its minimum, which corresponds to the time of
    maximum eclipse. It uses golden section search for sub-second precision.

    For global eclipse maximum (no observer location given):
        Uses the gamma parameter (minimum distance of Moon's shadow axis from
        Earth's center) to find global maximum. This is the time when the
        eclipse is at its maximum considering all possible Earth observers.

    For local eclipse maximum (observer location given):
        Uses the topocentric angular separation between Sun and Moon as seen
        from the specified observer location. This accounts for parallax and
        gives the precise local maximum time.

    Args:
        jd_approx: Approximate Julian Day (UT) of eclipse, typically near New Moon
                   or obtained from _sol_eclipse_when_glob_pythonic()[0]
        lat: Observer latitude in degrees (positive = North, negative = South).
             If None, calculates global maximum.
        lon: Observer longitude in degrees (positive = East, negative = West).
             If None, calculates global maximum.
        altitude: Observer altitude in meters above sea level (default 0).
                  Only used if lat and lon are provided.
        search_range: Search range in days around jd_approx (default ±3 hours).
                      For best results, jd_approx should be within this range
                      of the actual eclipse maximum.
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple of (jd_maximum, min_separation) where:
            - jd_maximum: Julian Day (UT) of maximum eclipse with sub-second precision
            - min_separation: Minimum angular separation in degrees between Sun and
                              Moon centers at maximum. For global maximum, this is
                              the gamma value in Earth radii (not degrees).

    Raises:
        ValueError: If only one of lat/lon is provided (both or neither required)

    Precision:
        Sub-second precision (< 1 second) achieved through golden section search
        with convergence to ~0.1 milliseconds.

    Algorithm:
        Global maximum:
            Uses Besselian elements approach - finds when gamma (perpendicular
            distance of Moon's shadow axis from Earth's center) is minimum.
            This corresponds to the time when the eclipse is at maximum for
            observers along the central line.

        Local maximum:
            Calculates topocentric apparent positions of Sun and Moon at each
            iteration, finding when their angular separation is minimum as
            seen from the observer's location on Earth's surface.

    Example:
        >>> # Find precise global maximum time for April 8, 2024 eclipse
        >>> from libephemeris import julday, sol_eclipse_max_time, _sol_eclipse_when_glob_pythonic
        >>> jd_start = julday(2024, 3, 1, 0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start)
        >>> jd_max, gamma = sol_eclipse_max_time(times[0])
        >>> print(f"Maximum at JD {jd_max:.8f}, gamma = {gamma:.6f}")

        >>> # Find local maximum time from Dallas, Texas
        >>> dallas_lat, dallas_lon = 32.7767, -96.7970
        >>> jd_local_max, separation = sol_eclipse_max_time(
        ...     times[0], lat=dallas_lat, lon=dallas_lon
        ... )
        >>> print(f"Local max at JD {jd_local_max:.8f}, sep = {separation:.6f}°")

    References:
        - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Reference documentation on Besselian elements
    """
    # Validate that both lat and lon are provided or neither
    if (lat is None) != (lon is None):
        raise ValueError(
            "Both lat and lon must be provided, or neither (for global maximum)"
        )

    # Determine if we're calculating global or local maximum
    is_local = lat is not None and lon is not None

    if is_local:
        # Local maximum: find minimum Sun-Moon separation from observer location
        assert lat is not None and lon is not None
        return _calc_local_eclipse_max_time(jd_approx, lat, lon, altitude, search_range)
    else:
        # Global maximum: find minimum gamma using Besselian elements
        return _calc_global_eclipse_max_time(jd_approx, search_range)


def _calc_global_eclipse_max_time(
    jd_approx: float, search_range: float
) -> Tuple[float, float]:
    """
    Calculate global eclipse maximum time using gamma minimization.

    Args:
        jd_approx: Approximate Julian Day of eclipse
        search_range: Search range in days

    Returns:
        Tuple of (jd_maximum, gamma_at_maximum)
    """
    phi = (1 + math.sqrt(5)) / 2  # Golden ratio

    jd_low = jd_approx - search_range
    jd_high = jd_approx + search_range

    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    gamma_a = _calc_gamma(jd_a)
    gamma_b = _calc_gamma(jd_b)

    # Golden section search for minimum gamma
    for _ in range(60):  # Converges to ~1e-9 days (~0.1 ms)
        if gamma_a < gamma_b:
            jd_high = jd_b
            jd_b = jd_a
            gamma_b = gamma_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            gamma_a = _calc_gamma(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            gamma_a = gamma_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            gamma_b = _calc_gamma(jd_b)

        if jd_high - jd_low < 1e-8:  # ~0.86 ms precision
            break

    jd_max = (jd_low + jd_high) / 2
    gamma_at_max = _calc_gamma(jd_max)

    return jd_max, gamma_at_max


def _calc_local_eclipse_max_time(
    jd_approx: float,
    lat: float,
    lon: float,
    altitude: float,
    search_range: float,
) -> Tuple[float, float]:
    """
    Calculate local eclipse maximum time using Sun-Moon separation minimization.

    Args:
        jd_approx: Approximate Julian Day of eclipse
        lat: Observer latitude in degrees
        lon: Observer longitude in degrees
        altitude: Observer altitude in meters
        search_range: Search range in days

    Returns:
        Tuple of (jd_maximum, min_separation_degrees)
    """
    from .state import get_leb_reader

    reader = get_leb_reader()

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        # geopos for topo/azalt: (lon, lat, alt)
        geopos = (lon, lat, altitude)

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon at given JD."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, geopos)
            return angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon at given JD."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            return sun_app.separation_from(moon_app).degrees

    # Golden section search for minimum separation
    phi = (1 + math.sqrt(5)) / 2

    jd_low = jd_approx - search_range
    jd_high = jd_approx + search_range

    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    sep_a = _get_separation(jd_a)
    sep_b = _get_separation(jd_b)

    for _ in range(60):  # Converges to ~1e-9 days (~0.1 ms)
        if sep_a < sep_b:
            jd_high = jd_b
            jd_b = jd_a
            sep_b = sep_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            sep_a = _get_separation(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            sep_a = sep_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            sep_b = _get_separation(jd_b)

        if jd_high - jd_low < 1e-8:  # ~0.86 ms precision
            break

    jd_max = (jd_low + jd_high) / 2
    min_separation = _get_separation(jd_max)

    return jd_max, min_separation


# Alias for reference API compatibility
sol_eclipse_max_time = sol_eclipse_max_time


def _sol_eclipse_when_glob_pythonic(
    jd_start: float,
    flags: int = FLG_SWIEPH,
    eclipse_type: int = 0,
    search_direction: str = "bidirectional",
) -> Tuple[int, Tuple[float, ...]]:
    """
    Find the next global solar eclipse after a given date.

    Searches forward in time from jd_start to find the next solar eclipse.
    Can filter by eclipse type (total, annular, partial, hybrid).

    By default uses bidirectional search to check ±15 days around jd_start,
    ensuring that eclipses close to the start date are not missed due to
    the search algorithm starting slightly after a New Moon.

    Args:
        jd_start: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        eclipse_type: Filter for specific eclipse type(s), bitmask of:
            - ECL_TOTAL (4): Total eclipse
            - ECL_ANNULAR (8): Annular eclipse
            - ECL_PARTIAL (16): Partial eclipse
            - ECL_ANNULAR_TOTAL (32): Hybrid eclipse
            - 0: Any eclipse type (default)
        search_direction: Direction to search for eclipses:
            - "forward": Only search forward from jd_start (per the reference API)
            - "backward": Only search backward from jd_start
            - "bidirectional" (default): Check ±15 days to ensure no eclipse
              is missed, returning the closest eclipse to jd_start that is
              >= jd_start (unless only found in backward direction)

    Returns:
        Tuple containing (matching reference API format):
            - retflag: Eclipse type flags bitmask (ECL_* constants)
            - tret: Tuple of 10 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse
                [1]: Time of first contact (partial begins)
                [2]: Time of second contact (central phase begins, or 0)
                [3]: Time of third contact (central phase ends, or 0)
                [4]: Time of fourth contact (partial ends)
                [5]: Time of sunrise on central line (or 0)
                [6]: Time of sunset on central line (or 0)
                [7]: Time when annular-total eclipse starts (or 0)
                [8]: Time when annular-total eclipse ends (or 0)
                [9]: Reserved (0)

    Raises:
        RuntimeError: If no eclipse found within search limit

    Algorithm:
        1. For bidirectional search: check previous New Moon within 15 days
        2. Find next New Moon after jd_start
        3. Check if Moon is close enough to node for eclipse
        4. If not eclipse, advance to next New Moon
        5. Calculate eclipse type and magnitude
        6. If eclipse_type filter set, check if matches
        7. Return closest matching eclipse >= jd_start

    Precision:
        Eclipse times accurate to ~1 minute for most eclipses.
        For higher precision, use Besselian elements method.

    Example:
        >>> # Find next total solar eclipse after Jan 1, 2024
        >>> from libephemeris import julday, ECL_TOTAL
        >>> jd = julday(2024, 1, 1, 0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd, eclipse_type=ECL_TOTAL)
        >>> print(f"Total eclipse at JD {times[0]:.5f}")

    References:
        - Reference API: sol_eclipse_when_glob()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 20  # Maximum search range
    MAX_NEW_MOONS = int(MAX_SEARCH_YEARS * 12.4)  # ~12.4 lunations per year
    BIDIRECTIONAL_WINDOW = 15.0  # Days to check backward in bidirectional mode

    # Validate search_direction parameter
    valid_directions = ("forward", "backward", "bidirectional")
    if search_direction not in valid_directions:
        raise ValueError(
            f"search_direction must be one of {valid_directions}, got '{search_direction}'"
        )

    # If eclipse_type is 0, accept any type
    if eclipse_type == 0:
        eclipse_type = ECL_ALLTYPES_SOLAR

    def _check_new_moon_for_eclipse(
        jd_new_moon: float,
    ) -> Union[Tuple[int, Tuple[float, ...]], None]:
        """Check if a New Moon corresponds to a matching eclipse."""
        # Get Moon position at New Moon
        moon_pos, _ = calc_ut(jd_new_moon, MOON, flags | FLG_SPEED)
        moon_lon = moon_pos[0]

        # Check if close enough to ecliptic for eclipse
        node_dist = _get_moon_node_distance(jd_new_moon, moon_lon)

        if node_dist < ECLIPSE_LIMIT_SOLAR:
            # Possible eclipse - check magnitude
            ecl_type, magnitude, gamma, ratio = _calculate_eclipse_type_and_magnitude(
                jd_new_moon
            )

            if ecl_type != 0:
                # Refine eclipse maximum using Besselian elements
                jd_max_refined = _refine_solar_eclipse_maximum(jd_new_moon)

                # Re-classify at the refined maximum time for accurate gamma.
                # The preliminary classification at New Moon can differ from the
                # refined one because gamma changes between the approximate New
                # Moon instant and the true eclipse maximum (e.g. an eclipse
                # classified as PARTIAL at New Moon may actually be central at
                # the refined maximum when gamma drops below 1.0).
                ecl_type_refined, mag_refined, gamma_refined, ratio_refined = (
                    _calculate_eclipse_type_and_magnitude(jd_max_refined)
                )
                if ecl_type_refined != 0:
                    ecl_type = ecl_type_refined

                # Check if matches filter
                type_matches = (
                    (eclipse_type & ECL_TOTAL and ecl_type & ECL_TOTAL)
                    or (eclipse_type & ECL_ANNULAR and ecl_type & ECL_ANNULAR)
                    or (eclipse_type & ECL_PARTIAL and ecl_type & ECL_PARTIAL)
                    or (
                        eclipse_type & ECL_ANNULAR_TOTAL
                        and ecl_type & ECL_ANNULAR_TOTAL
                    )
                )

                if type_matches:
                    # Calculate phase times using high-precision Besselian method
                    times = _calculate_eclipse_phases(jd_max_refined, ecl_type)
                    return ecl_type, times

        return None

    # For bidirectional search, first check if there's a recent eclipse we might miss
    backward_result: Union[Tuple[int, Tuple[float, ...]], None] = None
    if search_direction == "bidirectional":
        # Check the previous New Moon (within ~15 days before jd_start)
        jd_prev_new_moon = _find_previous_new_moon(jd_start)

        # Only consider if within the bidirectional window
        if jd_start - jd_prev_new_moon <= BIDIRECTIONAL_WINDOW:
            result = _check_new_moon_for_eclipse(jd_prev_new_moon)
            if result is not None:
                # Check if eclipse maximum is at or after jd_start
                # (eclipse might have started before but maximum after)
                jd_max = result[1][0]
                if jd_max >= jd_start:
                    # This eclipse's maximum is after jd_start, use it
                    return result
                else:
                    # Store for comparison - eclipse maximum is before jd_start
                    # but we'll still prefer forward search unless nothing found
                    backward_result = result

    # Handle backward-only search
    if search_direction == "backward":
        # The conjunction in longitude can lie AFTER jd_start while the
        # eclipse maximum lies a few minutes BEFORE it: probe the next
        # conjunction first so that eclipse is not skipped.
        jd_next_nm = _find_next_new_moon(jd_start)
        if jd_next_nm - jd_start <= BIDIRECTIONAL_WINDOW:
            result = _check_new_moon_for_eclipse(jd_next_nm)
            if result is not None and result[1][0] < jd_start:
                return result

        jd = jd_start
        for _ in range(MAX_NEW_MOONS):
            jd_prev_new_moon = _find_previous_new_moon(jd)
            result = _check_new_moon_for_eclipse(jd_prev_new_moon)
            if result is not None:
                # Direction invariant: a backward search must return an
                # eclipse whose maximum precedes jd_start (refinement can
                # drag the maximum across the start time).
                if result[1][0] < jd_start:
                    return result
            # Go further back
            jd = jd_prev_new_moon - 1

        raise RuntimeError(
            f"No matching solar eclipse found within {MAX_SEARCH_YEARS} years before JD {jd_start}"
        )

    # Forward search (default behavior for "forward" and "bidirectional")
    jd = jd_start

    for _ in range(MAX_NEW_MOONS):
        # Find next New Moon
        jd_new_moon = _find_next_new_moon(jd)

        result = _check_new_moon_for_eclipse(jd_new_moon)
        if result is not None and result[1][0] <= jd_start:
            # Direction invariant: a forward search must return an eclipse
            # whose maximum follows jd_start. The conjunction-anchored
            # candidate can refine to a maximum a few minutes BEFORE the
            # start time; skip to the next lunation instead of returning
            # a past eclipse (or looping forever in tret[0]+eps scans).
            result = None
        if result is not None:
            # For bidirectional mode, we have a forward result
            # Check if we had a backward result that might be closer
            if search_direction == "bidirectional" and backward_result is not None:
                # Compare distances from jd_start
                backward_jd_max = backward_result[1][0]
                forward_jd_max = result[1][0]

                # If backward eclipse max is closer to jd_start (but before it),
                # and forward is far away, still prefer forward since it's >= jd_start
                # The bidirectional check already returned if backward was >= jd_start
                pass  # Forward result is >= jd_start, use it

            return result

        # Advance to next lunation
        jd = jd_new_moon + 25  # Skip ahead ~25 days to ensure we find next New Moon

    raise RuntimeError(
        f"No matching solar eclipse found within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


def sol_eclipse_when_glob(
    tjdut: float,
    flags: int = FLG_SWIEPH,
    ecltype: int = 0,
    backwards: bool = False,
) -> Tuple[int, Tuple[float, ...]]:
    """Find the next (or previous) global solar eclipse (pyswisseph-compatible).

    Wrapper around _sol_eclipse_when_glob_pythonic() matching pyswisseph signature.

    Args:
        tjdut: Julian Day (UT) to start search from.
        flags: Calculation flags (default FLG_SWIEPH).
        ecltype: Eclipse type filter bitmask (0 = any).
        backwards: If True, search backward in time.

    Returns:
        Tuple of (retflag, tret) matching pyswisseph.
    """
    direction = "backward" if backwards else "forward"
    return _sol_eclipse_when_glob_pythonic(
        tjdut, flags=flags, eclipse_type=ecltype, search_direction=direction
    )


def _calculate_local_eclipse_phases(
    jd_max_global: float,
    lat: float,
    lon: float,
    altitude: float,
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """Calculate eclipse phase times as seen from a specific location.

    Wrapper around :func:`_calculate_local_eclipse_phases_impl` that adds
    LEB→Skyfield fallback for partial/custom LEB files. See the impl
    docstring for the full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _calculate_local_eclipse_phases_impl,
        jd_max_global, lat, lon, altitude,
    )


def _calculate_local_eclipse_phases_impl(
    jd_max_global: float,
    lat: float,
    lon: float,
    altitude: float,
    *,
    reader=None,
) -> Tuple[float, float, float, float, float, float, float, float, float, float]:
    """
    Calculate eclipse phase times as seen from a specific location.

    Uses geometric calculations based on the Moon's shadow passing over
    the observer's location.

    Args:
        jd_max_global: Julian Day of global maximum eclipse
        lat: Observer latitude in degrees (positive = North)
        lon: Observer longitude in degrees (positive = East)
        altitude: Observer altitude in meters above sea level

    Returns:
        Tuple of 10 floats with local eclipse circumstances:
            [0]: Time of maximum eclipse (local)
            [1]: Time of first contact (partial begins)
            [2]: Time of second contact (total/annular begins, or 0)
            [3]: Time of third contact (total/annular ends, or 0)
            [4]: Time of fourth contact (partial ends)
            [5]: Eclipse magnitude at maximum
            [6]: Fraction of Sun's diameter obscured
            [7]: Eclipse obscuration (area)
            [8]: Azimuth of Sun at maximum eclipse
            [9]: Altitude of Sun at maximum eclipse
    """
    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        try:
            from .fast_calc import _topo_ecliptic
            from .time_utils import deltat
            _topo_ecliptic(reader, jd_max_global + deltat(jd_max_global),
                           jd_max_global, SUN, (lon, lat, altitude))
        except (KeyError, ValueError):
            reader = None

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        geopos = (lon, lat, altitude)

        def _get_sun_moon_separation(jd: float) -> float:
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, geopos)
            return angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])

        def _get_sun_altaz(jd: float) -> Tuple[float, float]:
            """Get Sun altitude and azimuth at given JD from observer location."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, geopos)
            sun_az, sun_alt_true, sun_alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, sun_pos[:3])
            # azalt returns SE convention azimuth (S=0)
            return sun_alt_true, sun_az

        def _get_distances(jd: float) -> Tuple[float, float]:
            """Get Sun and Moon distances at given JD."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, geopos)
            return sun_pos[2], moon_pos[2]
    else:
        from skyfield.api import wgs84

        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        def _get_sun_moon_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon at given JD."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            # Get apparent positions from observer
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            # Calculate angular separation
            sep = sun_app.separation_from(moon_app)
            return sep.degrees

        def _get_sun_altaz(jd: float) -> Tuple[float, float]:
            """Get Sun altitude and azimuth at given JD from observer location."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            alt, az, _ = sun_app.altaz()
            # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_distances(jd: float) -> Tuple[float, float]:
            """Get Sun and Moon distances at given JD."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
            return sun_app.distance().au, moon_app.distance().au

    # Search for local maximum around global maximum
    # The local maximum can differ from global by up to ~1 hour
    search_range = 3.0 / 24.0  # 3 hours in days

    # First check if Sun is above horizon at global maximum
    sun_alt, sun_az = _get_sun_altaz(jd_max_global)
    if sun_alt < -1.0:  # Sun below horizon (with some margin for refraction)
        # Sun not visible - return zeros to indicate no local eclipse
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Find local maximum (minimum separation)
    # Use golden section search for efficiency
    jd_low = jd_max_global - search_range
    jd_high = jd_max_global + search_range

    phi = (1 + math.sqrt(5)) / 2  # Golden ratio
    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    sep_a = _get_sun_moon_separation(jd_a)
    sep_b = _get_sun_moon_separation(jd_b)

    for _ in range(30):  # Converge to ~0.1 second precision
        if sep_a < sep_b:
            jd_high = jd_b
            jd_b = jd_a
            sep_b = sep_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            sep_a = _get_sun_moon_separation(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            sep_a = sep_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            sep_b = _get_sun_moon_separation(jd_b)

        if jd_high - jd_low < 1e-7:  # ~0.01 second
            break

    jd_local_max = (jd_low + jd_high) / 2
    min_separation = _get_sun_moon_separation(jd_local_max)

    # Check if Sun is visible at local maximum
    sun_alt_max, sun_az_max = _get_sun_altaz(jd_local_max)
    if sun_alt_max < -1.0:
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Get angular sizes of Sun and Moon at local maximum
    sun_dist_au, moon_dist_au = _get_distances(jd_local_max)

    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au  # degrees
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)  # degrees

    sun_diameter = 2 * sun_angular_radius
    moon_diameter = 2 * moon_angular_radius

    # Check if there's an eclipse at this location
    # Eclipse occurs if separation < sum of radii
    sum_radii = sun_angular_radius + moon_angular_radius
    if min_separation > sum_radii:
        # No eclipse visible from this location
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Calculate eclipse magnitude
    # Magnitude = fraction of Sun's diameter covered by Moon
    overlap = sum_radii - min_separation
    # Edge case: handle very small sun_diameter (shouldn't happen but be safe)
    if sun_diameter < MINIMUM_SEPARATION_FOR_LENS:
        magnitude = 0.0
    else:
        magnitude = overlap / sun_diameter
        magnitude = max(0.0, min(magnitude, 1.0 + moon_diameter / sun_diameter))

    # Calculate obscuration using safe function that handles edge cases
    r_sun = sun_angular_radius
    r_moon = moon_angular_radius
    d = min_separation  # center-to-center separation

    obscuration = _calculate_obscuration_safe(r_sun, r_moon, d)

    # Find contact times using bisection
    # First/fourth contact: separation = sum of radii
    # Second/third contact: separation = difference of radii (for central eclipse)

    def _find_contact_time(
        jd_start: float, jd_end: float, target_sep: float, is_increasing: bool
    ) -> float:
        """Find time when separation equals target, searching in given direction.

        Includes edge case handling for shallow eclipses where contact may
        not occur or search may fail to converge.
        """
        jd_low = jd_start
        jd_high = jd_end

        # Edge case: check if target separation is achievable in search range
        sep_start = _get_sun_moon_separation(jd_start)
        sep_end = _get_sun_moon_separation(jd_end)

        # For decreasing search: sep should go from above to below target
        # For increasing search: sep should go from below to above target
        if is_increasing:
            if sep_start > target_sep or sep_end < target_sep:
                # Check if minimum separation in range is above target
                sep_mid = _get_sun_moon_separation((jd_start + jd_end) / 2)
                if min(sep_start, sep_end, sep_mid) > target_sep:
                    return 0.0  # Contact not achievable
        else:
            if sep_start < target_sep or sep_end > target_sep:
                # Check if minimum separation in range is above target
                sep_mid = _get_sun_moon_separation((jd_start + jd_end) / 2)
                if min(sep_start, sep_end, sep_mid) > target_sep:
                    return 0.0  # Contact not achievable

        for _ in range(50):
            jd_mid = (jd_low + jd_high) / 2
            sep = _get_sun_moon_separation(jd_mid)

            if abs(sep - target_sep) < 1e-6:
                return jd_mid

            if is_increasing:
                # Looking for separation increasing through target
                if sep < target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid
            else:
                # Looking for separation decreasing through target
                if sep > target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid

            if jd_high - jd_low < 1e-8:
                break

        result = (jd_low + jd_high) / 2
        # Validate the result - for shallow eclipses, contact may not be precise
        final_sep = _get_sun_moon_separation(result)
        if abs(final_sep - target_sep) > 0.001:  # More than ~4 arcsec error
            return 0.0  # Contact time not reliable

        return result

    # Search range for contacts
    contact_search = 2.0 / 24.0  # 2 hours

    # First contact (separation decreasing to sum of radii)
    jd_first = _find_contact_time(
        jd_local_max - contact_search, jd_local_max, sum_radii, is_increasing=False
    )

    # Fourth contact (separation increasing from sum of radii)
    jd_fourth = _find_contact_time(
        jd_local_max, jd_local_max + contact_search, sum_radii, is_increasing=True
    )

    # Check if first and fourth contacts are valid (Sun above horizon)
    first_alt, _ = _get_sun_altaz(jd_first)
    fourth_alt, _ = _get_sun_altaz(jd_fourth)

    if first_alt < -1.0:
        jd_first = 0.0
    if fourth_alt < -1.0:
        jd_fourth = 0.0

    # Second and third contacts (for total/annular eclipses)
    jd_second = 0.0
    jd_third = 0.0

    diff_radii = abs(moon_angular_radius - sun_angular_radius)
    if min_separation < diff_radii:
        # Central eclipse at this location
        jd_second = _find_contact_time(
            jd_local_max - contact_search / 4,
            jd_local_max,
            diff_radii,
            is_increasing=False,
        )
        jd_third = _find_contact_time(
            jd_local_max,
            jd_local_max + contact_search / 4,
            diff_radii,
            is_increasing=True,
        )

        # Check visibility
        second_alt, _ = _get_sun_altaz(jd_second)
        third_alt, _ = _get_sun_altaz(jd_third)

        if second_alt < -1.0:
            jd_second = 0.0
        if third_alt < -1.0:
            jd_third = 0.0

    return (
        jd_local_max,  # [0] Time of maximum
        jd_first,  # [1] First contact
        jd_second,  # [2] Second contact
        jd_third,  # [3] Third contact
        jd_fourth,  # [4] Fourth contact
        magnitude,  # [5] Magnitude
        overlap / sun_diameter,  # [6] Fraction covered
        obscuration,  # [7] Obscuration
        sun_az_max,  # [8] Azimuth
        sun_alt_max,  # [9] Altitude
    )


def _determine_local_eclipse_type(
    magnitude: float, moon_sun_ratio: float, is_central: bool
) -> int:
    """
    Determine eclipse type flags for local observation.

    Args:
        magnitude: Eclipse magnitude at location
        moon_sun_ratio: Ratio of Moon's apparent diameter to Sun's
        is_central: Whether this is a central eclipse at location

    Returns:
        Eclipse type flags bitmask
    """
    if magnitude <= 0:
        return 0

    flags = ECL_VISIBLE

    if is_central:
        flags |= ECL_CENTRAL
        if moon_sun_ratio >= 1.0:
            flags |= ECL_TOTAL
        elif moon_sun_ratio > 0.99:
            flags |= ECL_ANNULAR_TOTAL
        else:
            flags |= ECL_ANNULAR
    else:
        flags |= ECL_PARTIAL

    return flags


def _sol_eclipse_when_loc_pythonic(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next solar eclipse visible from a specific location.

    Searches forward in time from jd_start to find the next solar eclipse
    visible from the observer's location. Returns the times of eclipse phases
    as seen from that position.

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 10 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse (local)
                [1]: Time of first contact (partial begins)
                [2]: Time of second contact (total/annular begins, or 0)
                [3]: Time of third contact (total/annular ends, or 0)
                [4]: Time of fourth contact (partial ends)
                [5]: Time of sunrise on central line (or 0, unused)
                [6]: Time of sunset on central line (or 0, unused)
                [7-9]: Reserved (0)
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Eclipse magnitude
                [1]: Fraction of Sun's diameter covered
                [2]: Fraction of Sun's area obscured (obscuration)
                [3]: Azimuth of Sun at maximum eclipse (degrees)
                [4]: Altitude of Sun at maximum eclipse (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of Sun (degrees)
                [7]: Saros series number (0, not implemented)
                [8]: Inex number (0, not implemented)
                [9-10]: Reserved (0)
            - retflag: Eclipse type flags bitmask (ECL_* constants)

    Raises:
        RuntimeError: If no eclipse visible from location within search limit

    Algorithm:
        1. Use _sol_eclipse_when_glob_pythonic to find next global eclipse
        2. Calculate local circumstances at observer's location
        3. If eclipse not visible from location, continue to next global eclipse
        4. Return local phase times and attributes

    Precision:
        Eclipse times accurate to ~1 minute. Contact times depend on
        accurate Moon and Sun ephemeris positions.

    Example:
        >>> # Find next eclipse visible from Rome
        >>> from libephemeris import julday, _sol_eclipse_when_loc_pythonic
        >>> jd = julday(2024, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> ecl_type, times, attr = _sol_eclipse_when_loc_pythonic(jd, rome_lat, rome_lon)
        >>> print(f"Eclipse maximum at JD {times[0]:.5f}")
        >>> print(f"Magnitude: {attr[0]:.3f}")

    References:
        - Reference API: sol_eclipse_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 50  # Maximum search range
    MAX_ECLIPSES = int(MAX_SEARCH_YEARS * 2.4)  # ~2.4 eclipses per year on average

    jd = jd_start

    for _ in range(MAX_ECLIPSES):
        # Find next global eclipse
        try:
            global_type, global_times = _sol_eclipse_when_glob_pythonic(jd, flags)
        except RuntimeError:
            raise RuntimeError(
                f"No solar eclipse visible from lat={lat}, lon={lon} "
                f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
            )

        jd_max_global = global_times[0]

        # Calculate local circumstances
        local_data = _calculate_local_eclipse_phases(jd_max_global, lat, lon, altitude)

        jd_local_max = local_data[0]
        magnitude = local_data[5]

        # Check if eclipse is visible from this location
        if jd_local_max > 0 and magnitude > 0:
            # Eclipse visible! Prepare return values

            # Get Moon/Sun topocentric distances at local maximum
            from .constants import FLG_TOPOCTR as _TOPOCTR_LOC
            from .state import get_topo as _gt_loc
            import libephemeris as _le_loc
            _saved_topo_loc = _gt_loc()
            _le_loc.set_topo(lon, lat, altitude)
            try:
                sun_pos, _ = calc_ut(jd_local_max, SUN, FLG_SPEED | _TOPOCTR_LOC)
                moon_pos, _ = calc_ut(jd_local_max, MOON, FLG_SPEED | _TOPOCTR_LOC)
            finally:
                from libephemeris import state as _st_loc
                if _saved_topo_loc is not None:
                    _le_loc.set_topo(_saved_topo_loc.longitude.degrees, _saved_topo_loc.latitude.degrees, _saved_topo_loc.elevation.m)
                else:
                    _st_loc._TOPO = None
            sun_dist_au = sun_pos[2]
            moon_dist_au = moon_pos[2]

            sun_diameter = 2 * (959.63 / 3600.0) / sun_dist_au
            moon_diameter = 2 * (932.56 / 3600.0) * (0.002569 / moon_dist_au)

            # Check for central eclipse at this location
            is_central = local_data[2] > 0 and local_data[3] > 0

            moon_sun_ratio = moon_diameter / sun_diameter

            # Determine eclipse type
            ecl_type = _determine_local_eclipse_type(
                magnitude, moon_sun_ratio, is_central
            )

            # Add visibility flags
            if local_data[1] > 0:
                ecl_type |= ECL_1ST_VISIBLE
            if local_data[2] > 0:
                ecl_type |= ECL_2ND_VISIBLE
            if local_data[3] > 0:
                ecl_type |= ECL_3RD_VISIBLE
            if local_data[4] > 0:
                ecl_type |= ECL_4TH_VISIBLE
            ecl_type |= ECL_MAX_VISIBLE

            # Prepare times tuple (10 elements)
            times = (
                local_data[0],  # [0] Maximum
                local_data[1],  # [1] First contact
                local_data[2],  # [2] Second contact
                local_data[3],  # [3] Third contact
                local_data[4],  # [4] Fourth contact
                0.0,  # [5] Sunrise on central line (not implemented)
                0.0,  # [6] Sunset on central line (not implemented)
                0.0,  # [7] Reserved
                0.0,  # [8] Reserved
                0.0,  # [9] Reserved
            )

            # Prepare attributes tuple (11 elements)
            attr = (
                local_data[5],  # [0] Magnitude
                local_data[6],  # [1] Fraction covered
                local_data[7],  # [2] Obscuration
                local_data[8],  # [3] Azimuth
                local_data[9],  # [4] Altitude
                moon_diameter,  # [5] Moon diameter
                sun_diameter,  # [6] Sun diameter
                0.0,  # [7] Saros (not implemented)
                0.0,  # [8] Inex (not implemented)
                0.0,  # [9] Reserved
                0.0,  # [10] Reserved
            )

            return ecl_type, times, attr

        # Eclipse not visible from this location, try next
        jd = jd_max_global + 25  # Skip ahead to find next eclipse

    raise RuntimeError(
        f"No solar eclipse visible from lat={lat}, lon={lon} "
        f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Legacy alias for original implementation
_sol_eclipse_when_loc_legacy = _sol_eclipse_when_loc_pythonic


def sol_eclipse_when_loc(
    tjdut: float,
    geopos: "Sequence[float]",
    flags: int = FLG_SWIEPH,
    backwards: bool = False,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find the next solar eclipse visible from a geographic location.

    Wrapper around :func:`_sol_eclipse_when_loc_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_when_loc_impl, tjdut, geopos, flags, backwards,
    )


def _sol_eclipse_when_loc_impl(
    tjdut: float,
    geopos: "Sequence[float]",
    flags: int = FLG_SWIEPH,
    backwards: bool = False,
    *,
    reader,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next solar eclipse visible from a specific geographic location.

    This function matches the reference API signature exactly. It searches forward
    (or backward if specified) in time from tjdut to find the next solar
    eclipse visible from the observer's location specified by geopos.

    Args:
        tjdut: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude_degrees, latitude_degrees, altitude_meters]
                NOTE: longitude comes first (this matches reference API convention)
        backwards: If True, search backward in time instead of forward

    Returns:
        Tuple containing:
            - tret: Tuple of 7 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse (local)
                [1]: Time of first contact (partial begins)
                [2]: Time of second contact (totality/annularity begins, or 0)
                [3]: Time of third contact (totality/annularity ends, or 0)
                [4]: Time of fourth contact (partial ends)
                [5]: Time of sunrise if eclipse at sunrise (or 0)
                [6]: Time of sunset if eclipse at sunset (or 0)
            - attr: Tuple of 8 floats with eclipse attributes:
                [0]: Fraction of solar diameter covered by Moon
                [1]: Ratio of lunar diameter to solar diameter
                [2]: Fraction of solar disc area covered (obscuration)
                [3]: Core shadow width in km (0 for partial eclipses)
                [4]: Azimuth of Sun at maximum eclipse (degrees)
                [5]: True altitude of Sun at maximum eclipse (degrees)
                [6]: Apparent altitude of Sun with refraction (degrees)
                [7]: Angular distance of Moon center from Sun center (degrees)
            - retflag: Eclipse type flags bitmask (ECL_* constants)
                ECL_VISIBLE if any part visible, combined with eclipse type

    Raises:
        RuntimeError: If no eclipse visible from location within search limit
        ValueError: If geopos has wrong length

    Algorithm:
        1. Use sol_eclipse_when_glob to find candidate global eclipses
        2. For each candidate, calculate topocentric positions using set_topo()
           and FLG_TOPOCTR flag
        3. Check if the eclipse is visible from the observer location
        4. Calculate contact times by solving for when lunar limb touches solar limb
        5. Calculate obscuration fraction as area covered / total solar disc area
        6. Verify Sun is above horizon at observer location during eclipse

    Precision:
        Eclipse times accurate to ~1-2 minutes. Contact times depend on
        accurate Moon and Sun ephemeris positions with topocentric corrections.

    Example:
        >>> # Find next eclipse visible from Dallas, TX (Apr 8, 2024 eclipse)
        >>> from libephemeris import sol_eclipse_when_loc, julday, FLG_SWIEPH
        >>> jd = julday(2024, 1, 1, 0)
        >>> dallas_geopos = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> ecl_type, tret, attr = sol_eclipse_when_loc(jd, FLG_SWIEPH, dallas_geopos)
        >>> print(f"Eclipse maximum at JD {tret[0]:.5f}")
        >>> print(f"Obscuration: {attr[2]:.3f}")

    References:
        - Reference API: sol_eclipse_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    from typing import Sequence

    # Validate geopos
    if len(geopos) < 3:
        raise ValueError("geopos must be a sequence of [longitude, latitude, altitude]")

    # Extract geographic position (longitude first, then latitude - reference API convention)
    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    MAX_SEARCH_YEARS = 50
    MAX_ECLIPSES = int(MAX_SEARCH_YEARS * 2.4)

    # NOTE: the observer position is threaded through explicitly (LEB path
    # geopos tuples / Skyfield wgs84.latlon); the reference API's
    # swe_sol_eclipse_when_loc does not touch the global set_topo state,
    # and neither do we (a previous set_topo() call here leaked the
    # eclipse observer into subsequent FLG_TOPOCTR calculations).

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        # geopos for topo/azalt: (lon, lat, alt)
        _geopos = (lon, lat, altitude)

        def _get_sun_moon_data(
            jd: float,
        ) -> Tuple[float, float, float, float, float, float, float]:
            """Get Sun-Moon separation and Sun position data at given JD."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _geopos)
            sep = angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])
            sun_az, sun_alt_true, sun_alt_app = azalt(jd, ECL2HOR, _geopos, 0, 0, sun_pos[:3])
            return (
                sep,
                sun_alt_true,
                sun_az,
                sun_pos[2],
                moon_pos[2],
                sun_alt_true,
                sun_az,
            )

        def _get_sun_altaz(jd: float) -> Tuple[float, float, float]:
            """Get Sun altitude, azimuth, and apparent altitude at given JD."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _geopos)
            sun_az, sun_alt_true, sun_alt_app = azalt(jd, ECL2HOR, _geopos, 0, 0, sun_pos[:3])
            # Calculate apparent altitude with refraction
            true_alt = sun_alt_true
            if true_alt > -1.0:
                if true_alt > 0:
                    refraction = 1.0 / math.tan(
                        math.radians(true_alt + 7.31 / (true_alt + 4.4))
                    )
                    apparent_alt = true_alt + refraction / 60.0
                else:
                    apparent_alt = true_alt + 0.58
            else:
                apparent_alt = true_alt
            return true_alt, sun_az, apparent_alt

        def _get_angular_sizes(jd: float) -> Tuple[float, float, float, float]:
            """Get angular radii and diameters of Sun and Moon."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _geopos)
            sun_angular_radius = (959.63 / 3600.0) / sun_pos[2]
            moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_pos[2])
            return (
                sun_angular_radius,
                moon_angular_radius,
                2 * sun_angular_radius,
                2 * moon_angular_radius,
            )

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _geopos)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _geopos)
            return angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])
    else:
        from skyfield.api import wgs84

        from .state import get_planets, get_timescale

        # Get ephemeris
        eph = get_planets()
        ts = get_timescale()

        # Create observer
        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]
        observer = wgs84.latlon(lat, lon, altitude)

        def _get_sun_moon_data(
            jd: float,
        ) -> Tuple[float, float, float, float, float, float, float]:
            """Get Sun-Moon separation and Sun position data at given JD."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            # Get apparent positions from observer (topocentric)
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            # Angular separation
            sep = sun_app.separation_from(moon_app).degrees

            # Sun altitude and azimuth
            sun_alt, sun_az, _ = sun_app.altaz()
            # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
            sun_az_se = (sun_az.degrees + 180.0) % 360.0

            # Distances
            sun_dist_au = sun_app.distance().au
            moon_dist_au = moon_app.distance().au

            return (
                sep,
                sun_alt.degrees,
                sun_az_se,
                sun_dist_au,
                moon_dist_au,
                sun_alt.degrees,
                sun_az_se,
            )

        def _get_sun_altaz(jd: float) -> Tuple[float, float, float]:
            """Get Sun altitude, azimuth, and apparent altitude at given JD."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            alt, az, _ = sun_app.altaz()

            # Calculate apparent altitude with refraction
            # Using standard atmospheric refraction formula
            true_alt = alt.degrees
            if true_alt > -1.0:
                # Bennett's formula for refraction
                if true_alt > 0:
                    refraction = 1.0 / math.tan(
                        math.radians(true_alt + 7.31 / (true_alt + 4.4))
                    )
                    apparent_alt = true_alt + refraction / 60.0
                else:
                    # Near horizon, use approximation
                    apparent_alt = true_alt + 0.58
            else:
                apparent_alt = true_alt

            # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
            return true_alt, (az.degrees + 180.0) % 360.0, apparent_alt

        def _get_angular_sizes(jd: float) -> Tuple[float, float, float, float]:
            """Get angular radii and diameters of Sun and Moon."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            sun_dist_au = sun_app.distance().au
            moon_dist_au = moon_app.distance().au

            # Angular radii in degrees
            sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
            moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

            return (
                sun_angular_radius,
                moon_angular_radius,
                2 * sun_angular_radius,
                2 * moon_angular_radius,
            )

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            return sun_app.separation_from(moon_app).degrees

    def _find_local_maximum(
        jd_center: float, search_range: float = 3.0 / 24.0
    ) -> float:
        """Find local eclipse maximum (minimum separation) using golden section search."""
        phi = (1 + math.sqrt(5)) / 2

        jd_low = jd_center - search_range
        jd_high = jd_center + search_range

        jd_a = jd_high - (jd_high - jd_low) / phi
        jd_b = jd_low + (jd_high - jd_low) / phi

        sep_a = _get_separation(jd_a)
        sep_b = _get_separation(jd_b)

        for _ in range(40):
            if sep_a < sep_b:
                jd_high = jd_b
                jd_b = jd_a
                sep_b = sep_a
                jd_a = jd_high - (jd_high - jd_low) / phi
                sep_a = _get_separation(jd_a)
            else:
                jd_low = jd_a
                jd_a = jd_b
                sep_a = sep_b
                jd_b = jd_low + (jd_high - jd_low) / phi
                sep_b = _get_separation(jd_b)

            if jd_high - jd_low < 1e-7:
                break

        return (jd_low + jd_high) / 2

    def _find_contact_time(
        jd_start_search: float,
        jd_end_search: float,
        target_sep: float,
        is_increasing: bool,
    ) -> float:
        """Find time when separation equals target using bisection."""
        jd_low = jd_start_search
        jd_high = jd_end_search

        for _ in range(50):
            jd_mid = (jd_low + jd_high) / 2
            sep = _get_separation(jd_mid)

            if abs(sep - target_sep) < 1e-6:
                return jd_mid

            if is_increasing:
                if sep < target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid
            else:
                if sep > target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid

            if jd_high - jd_low < 1e-8:
                break

        return (jd_low + jd_high) / 2

    def _calculate_obscuration(
        sun_radius: float, moon_radius: float, separation: float
    ) -> float:
        """Calculate fraction of solar disc area covered by Moon."""
        r_sun = sun_radius
        r_moon = moon_radius
        d = separation

        if d >= r_sun + r_moon:
            return 0.0
        elif d <= abs(r_sun - r_moon):
            # Both total and annular: obscuration = disc area ratio
            # For total eclipses (r_moon >= r_sun), this is >= 1.0
            # matching reference API behavior
            return (r_moon / r_sun) ** 2
        else:
            # Partial overlap - lens formula
            d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
            d2 = d - d1

            if abs(d1) <= r_sun and abs(d2) <= r_moon:
                cos_arg1 = max(-1, min(1, d1 / r_sun))
                cos_arg2 = max(-1, min(1, d2 / r_moon))

                area1 = r_sun * r_sun * math.acos(cos_arg1) - d1 * math.sqrt(
                    max(0, r_sun * r_sun - d1 * d1)
                )
                area2 = r_moon * r_moon * math.acos(cos_arg2) - d2 * math.sqrt(
                    max(0, r_moon * r_moon - d2 * d2)
                )
                return (area1 + area2) / (math.pi * r_sun * r_sun)
            else:
                return 0.0

    def _find_previous_new_moon(jd: float) -> float:
        """Find the previous New Moon before jd."""
        # Start searching backward
        jd_search = jd
        for _ in range(30):  # At most ~2 years back
            jd_nm = _find_next_new_moon(jd_search - 35)
            if jd_nm < jd:
                # Found a new moon before jd, but we want the most recent one
                # Keep searching forward from here until we pass jd
                while jd_nm < jd:
                    jd_candidate = jd_nm
                    jd_nm = _find_next_new_moon(jd_nm + 25)
                return jd_candidate
            jd_search -= 30
        return jd - 29.5  # Fallback

    # Main search loop
    jd = tjdut

    # For backward search, we need to find eclipses before tjdut
    if backwards:
        # Find the most recent global eclipse before tjdut
        # We need to search backward in time
        search_direction = -1
    else:
        search_direction = 1

    for _ in range(MAX_ECLIPSES):
        # Find next/previous global eclipse
        try:
            if backwards:
                # For backward search, find a new moon before current position
                # and check if there's an eclipse there
                # We use a simpler approach: search forward from an earlier date
                # and find eclipses until we get one before tjdut
                earlier_jd = jd - 400  # About 1 year before
                eclipses_found = []
                temp_jd = earlier_jd
                while temp_jd < jd:
                    try:
                        global_type, global_times = _sol_eclipse_when_glob_pythonic(
                            temp_jd, flags
                        )
                        if global_times[0] < jd:
                            eclipses_found.append((global_type, global_times))
                        temp_jd = global_times[0] + 25
                    except RuntimeError:
                        break
                if not eclipses_found:
                    jd -= 400
                    continue
                # Take the most recent eclipse before jd
                global_type, global_times = eclipses_found[-1]
            else:
                global_type, global_times = _sol_eclipse_when_glob_pythonic(jd, flags)
        except RuntimeError:
            raise RuntimeError(
                f"No solar eclipse visible from lon={lon}, lat={lat} "
                f"within {MAX_SEARCH_YEARS} years of JD {tjdut}"
            )

        jd_max_global = global_times[0]

        # Find local maximum at this location
        jd_local_max = _find_local_maximum(jd_max_global)

        # Get angular sizes at local maximum
        sun_radius, moon_radius, sun_diam, moon_diam = _get_angular_sizes(jd_local_max)
        min_separation = _get_separation(jd_local_max)

        # Check if eclipse is visible (Moon overlaps Sun)
        sum_radii = sun_radius + moon_radius
        if min_separation > sum_radii:
            # No eclipse visible from this location
            if backwards:
                jd = jd_max_global - 1
            else:
                jd = jd_max_global + 25
            continue

        # Check if Sun is above horizon at local max
        true_alt, sun_az, apparent_alt = _get_sun_altaz(jd_local_max)

        if apparent_alt < 0.0:
            # Sun is below horizon (accounting for refraction) at local
            # maximum, but the eclipse may still be visible during sunrise
            # or sunset. Check if the eclipse is ongoing (separation <
            # sum_radii) when the Sun rises above the horizon. Search for
            # a time when both conditions are met: Sun visible (apparent
            # altitude >= 0) AND eclipse still in progress.
            horizon_eclipse_found = False

            # Check times around the local max for Sun above horizon
            # while eclipse is still in progress
            for dt_minutes in range(-180, 181, 5):
                jd_test = jd_local_max + dt_minutes / 1440.0
                test_alt, test_az, test_app_alt = _get_sun_altaz(jd_test)
                if test_app_alt >= 0.0:
                    test_sep = _get_separation(jd_test)
                    if test_sep < sum_radii:
                        # Eclipse visible above horizon at this time!
                        # Use this as the effective local maximum for
                        # the visible portion of the eclipse
                        horizon_eclipse_found = True
                        # Find the best visible time (minimum separation
                        # while Sun is above horizon)
                        jd_best = jd_test
                        sep_best = test_sep
                        for dt2 in range(dt_minutes, min(dt_minutes + 120, 181), 2):
                            jd_t2 = jd_local_max + dt2 / 1440.0
                            alt_t2, _, app_alt_t2 = _get_sun_altaz(jd_t2)
                            if app_alt_t2 < 0.0:
                                break
                            sep_t2 = _get_separation(jd_t2)
                            if sep_t2 >= sum_radii:
                                break
                            if sep_t2 < sep_best:
                                sep_best = sep_t2
                                jd_best = jd_t2
                        # Update local max to the best visible time
                        jd_local_max = jd_best
                        min_separation = sep_best
                        true_alt, sun_az, apparent_alt = _get_sun_altaz(jd_local_max)
                        # Recalculate angular sizes at new time
                        sun_radius, moon_radius, sun_diam, moon_diam = (
                            _get_angular_sizes(jd_local_max)
                        )
                        sum_radii = sun_radius + moon_radius
                        break

            if not horizon_eclipse_found:
                # Eclipse truly not visible from this location
                if backwards:
                    jd = jd_max_global - 1
                else:
                    jd = jd_max_global + 25
                continue

        # Eclipse visible! Calculate all parameters

        # Calculate magnitude (fraction of Sun diameter covered)
        overlap = sum_radii - min_separation
        magnitude = overlap / sun_diam
        magnitude = max(0.0, min(magnitude, 1.0 + moon_diam / sun_diam))

        # Calculate obscuration (area fraction)
        obscuration = _calculate_obscuration(sun_radius, moon_radius, min_separation)

        # Ratio of diameters
        ratio = moon_diam / sun_diam

        # Calculate contact times
        contact_search_range = 2.0 / 24.0  # 2 hours

        # First contact (separation decreasing to sum of radii)
        jd_first = _find_contact_time(
            jd_local_max - contact_search_range,
            jd_local_max,
            sum_radii,
            is_increasing=False,
        )

        # Fourth contact (separation increasing from sum of radii)
        jd_fourth = _find_contact_time(
            jd_local_max,
            jd_local_max + contact_search_range,
            sum_radii,
            is_increasing=True,
        )

        # Second and third contacts (for central eclipses)
        diff_radii = abs(moon_radius - sun_radius)
        if min_separation < diff_radii:
            # Central eclipse at this location
            jd_second = _find_contact_time(
                jd_local_max - contact_search_range / 4,
                jd_local_max,
                diff_radii,
                is_increasing=False,
            )
            jd_third = _find_contact_time(
                jd_local_max,
                jd_local_max + contact_search_range / 4,
                diff_radii,
                is_increasing=True,
            )
        else:
            jd_second = 0.0
            jd_third = 0.0

        # Check visibility of each contact — save original geometric times
        # for sunrise/sunset search before zeroing invisible contacts.
        jd_first_geom = jd_first
        jd_fourth_geom = jd_fourth

        first_alt, _, _ = _get_sun_altaz(jd_first)
        fourth_alt, _, _ = _get_sun_altaz(jd_fourth)

        if first_alt < -1.0:
            jd_first = 0.0
        if fourth_alt < -1.0:
            jd_fourth = 0.0

        if jd_second > 0:
            second_alt, _, _ = _get_sun_altaz(jd_second)
            if second_alt < -1.0:
                jd_second = 0.0

        if jd_third > 0:
            third_alt, _, _ = _get_sun_altaz(jd_third)
            if third_alt < -1.0:
                jd_third = 0.0

        # Check for sunrise/sunset during eclipse using original geometric
        # contact times (before visibility filtering zeroed them out).
        jd_sunrise = 0.0
        jd_sunset = 0.0

        # Check if Sun rises during eclipse
        if jd_first_geom > 0 and jd_fourth_geom > 0:
            alt_first, _, _ = _get_sun_altaz(jd_first_geom)
            alt_fourth, _, _ = _get_sun_altaz(jd_fourth_geom)

            if alt_first < 0 < alt_fourth:
                # Sun rises during eclipse - find approximate time
                t_low, t_high = jd_first_geom, jd_fourth_geom
                for _ in range(30):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _, _ = _get_sun_altaz(t_mid)
                    if alt_mid < 0:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                    if t_high - t_low < 1e-7:
                        break
                jd_sunrise = (t_low + t_high) / 2

            elif alt_first > 0 > alt_fourth:
                # Sun sets during eclipse
                t_low, t_high = jd_first_geom, jd_fourth_geom
                for _ in range(30):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _, _ = _get_sun_altaz(t_mid)
                    if alt_mid > 0:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                    if t_high - t_low < 1e-7:
                        break
                jd_sunset = (t_low + t_high) / 2

        # Calculate core shadow width.
        # This is always computed regardless of whether the eclipse is central
        # at the observer's location — it represents the shadow cone geometry
        # projected through the observer's Sun altitude at maximum.
        moon_pos, _ = calc_ut(jd_local_max, MOON, flags | FLG_SPEED)
        sun_pos, _ = calc_ut(jd_local_max, SUN, flags | FLG_SPEED)

        moon_dist_au = moon_pos[2]
        sun_dist_au = sun_pos[2]

        sun_radius_km = 696340.0
        moon_radius_km = 1737.4
        sun_dist_km = sun_dist_au * 149597870.7
        moon_dist_km = moon_dist_au * 149597870.7

        # Umbra/antumbra cone geometry
        alpha = math.atan((sun_radius_km - moon_radius_km) / sun_dist_km)

        if ratio >= 1.0:
            # Total eclipse - umbra
            umbra_radius_km = moon_radius_km - moon_dist_km * math.tan(alpha)
            umbra_radius_km = max(0, umbra_radius_km)
            is_total_shadow = True
        else:
            # Annular - antumbra
            umbra_radius_km = moon_dist_km * math.tan(alpha) - moon_radius_km
            umbra_radius_km = max(0, abs(umbra_radius_km))
            is_total_shadow = False

        # Path width affected by Sun altitude
        if true_alt > 0:
            cos_alt = math.cos(math.radians(90 - true_alt))
            cos_alt = max(0.1, cos_alt)
            shadow_width_km = 2 * umbra_radius_km / cos_alt
        else:
            shadow_width_km = 0.0

        shadow_width_km = min(1000.0, shadow_width_km)

        # Sign convention: negative for total eclipses (umbra),
        # positive for annular eclipses (antumbra)
        if is_total_shadow:
            shadow_width_km = -shadow_width_km

        # Determine eclipse type flags
        ecl_type = ECL_VISIBLE

        is_central = jd_second > 0 and jd_third > 0
        if is_central:
            ecl_type |= ECL_CENTRAL
            if ratio >= 1.0:
                ecl_type |= ECL_TOTAL
            elif ratio > 0.99:
                ecl_type |= ECL_ANNULAR_TOTAL
            else:
                ecl_type |= ECL_ANNULAR
        else:
            ecl_type |= ECL_PARTIAL

        # Add visibility flags
        if jd_first > 0:
            ecl_type |= ECL_1ST_VISIBLE
        if jd_second > 0:
            ecl_type |= ECL_2ND_VISIBLE
        if jd_third > 0:
            ecl_type |= ECL_3RD_VISIBLE
        if jd_fourth > 0:
            ecl_type |= ECL_4TH_VISIBLE
        ecl_type |= ECL_MAX_VISIBLE

        # Prepare tret tuple (10 elements matching pyswisseph layout)
        tret = (
            jd_local_max,  # [0] Maximum eclipse
            jd_first,  # [1] First contact
            jd_second,  # [2] Second contact
            jd_third,  # [3] Third contact
            jd_fourth,  # [4] Fourth contact
            jd_sunrise,  # [5] Sunrise during eclipse
            jd_sunset,  # [6] Sunset during eclipse
            0.0,  # [7] Reserved
            0.0,  # [8] Reserved
            0.0,  # [9] Reserved
        )

        # Prepare attr tuple (20 elements matching pyswisseph layout)
        attr = (
            magnitude,  # [0] Fraction of solar diameter covered
            ratio,  # [1] Ratio of lunar to solar diameter
            obscuration,  # [2] Fraction of solar disc area covered
            shadow_width_km,  # [3] Core shadow width in km
            sun_az,  # [4] Azimuth of sun at maximum
            true_alt,  # [5] True altitude of sun at maximum
            apparent_alt,  # [6] Apparent altitude with refraction
            min_separation,  # [7] Angular distance Moon-Sun centers
            ratio
            if (ecl_type & (ECL_TOTAL | ECL_ANNULAR))
            else magnitude,  # [8] Magnitude acc. to NASA
            *_get_saros_info(jd_max_global, "solar"),  # [9] Saros, [10] member
            0.0,  # [11] Reserved
            0.0,  # [12] Reserved
            0.0,  # [13] Reserved
            0.0,  # [14] Reserved
            0.0,  # [15] Reserved
            0.0,  # [16] Reserved
            0.0,  # [17] Reserved
            0.0,  # [18] Reserved
            0.0,  # [19] Reserved
        )

        return ecl_type, tret, attr

    raise RuntimeError(
        f"No solar eclipse visible from lon={lon}, lat={lat} "
        f"within {MAX_SEARCH_YEARS} years of JD {tjdut}"
    )


def _sol_eclipse_where_pythonic(
    jd: float,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Calculate the geographic location where a solar eclipse is central at a given time.

    This is the legacy function signature. For reference API compatibility,
    use sol_eclipse_where().

    Args:
        jd: Julian Day (UT) of the moment to calculate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        See sol_eclipse_where() for full return specification.
    """
    return sol_eclipse_where(jd, flags)


def sol_eclipse_where(
    tjdut: float,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the geographic location where a solar eclipse is at maximum at a given time.

    This function determines where on Earth the Moon's shadow axis intersects
    the Earth's surface at the specified Julian Day, using the WGS84 ellipsoid
    for accurate geodetic calculations.

    Args:
        tjdut: Julian Day (UT) during a solar eclipse
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - retflag: Eclipse type flags bitmask (ECL_* constants)
                Returns 0 if no central eclipse at this time
            - geopos: Tuple of 10 floats with geographic position per reference API:
                [0]: Geographic longitude of central line (degrees, East positive)
                [1]: Geographic latitude of central line (degrees, North positive)
                [2]: Longitude of northern limit of umbra
                [3]: Latitude of northern limit of umbra
                [4]: Longitude of southern limit of umbra
                [5]: Latitude of southern limit of umbra
                [6]: Longitude of northern limit of penumbra
                [7]: Latitude of northern limit of penumbra
                [8]: Longitude of southern limit of penumbra
                [9]: Latitude of southern limit of penumbra
            - attr: Tuple of 20 floats with eclipse attributes per reference API:
                [0]: Fraction of solar diameter covered by Moon (magnitude)
                [1]: Ratio of lunar diameter to solar diameter
                [2]: Fraction of Sun's area obscured (obscuration)
                [3]: Width of totality/annularity path (km)
                [4]: Sun's azimuth at central line (degrees)
                [5]: True altitude of Sun at central line (degrees)
                [6]: Apparent altitude with refraction (degrees)
                [7]: Angular distance Moon center from Sun center (degrees)
                [8-19]: Reserved for future use

    Note:
        If the shadow axis misses the Earth (partial eclipse or no eclipse),
        the function finds the point of maximum penumbral coverage instead.
        Uses WGS84 ellipsoid: a=6378.137 km, b=6356.752 km, f=1/298.257223563

    Algorithm:
        1. Start with sub-lunar point as initial guess
        2. Iteratively refine to find point of minimum Sun-Moon separation
        3. This accounts for parallax and gives accurate central line position
        4. Calculate eclipse attributes at that location

    Example:
        >>> # Find central line location during April 8, 2024 eclipse
        >>> from libephemeris import julday, sol_eclipse_where, FLG_SWIEPH
        >>> jd = 2460409.26  # ~18:18 UTC during maximum
        >>> ecl_type, geopos, attr = sol_eclipse_where(jd, FLG_SWIEPH)
        >>> print(f"Central at lon={geopos[0]:.2f}, lat={geopos[1]:.2f}")

    References:
        - Reference API: sol_eclipse_where()
        - Meeus "Astronomical Algorithms" Ch. 54
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
    """
    from .state import get_leb_reader

    AU_TO_KM = 149597870.7  # AU to km conversion

    reader = get_leb_reader()

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        # Get geocentric Moon position for initial guess (via calc_ut)
        moon_eq, _ = calc_ut(tjdut, MOON, FLG_EQUATORIAL | FLG_SPEED)
        moon_ra_deg = moon_eq[0]  # degrees
        moon_dec_deg = moon_eq[1]

        # Calculate sub-lunar point as initial guess
        jd_tt = tjdut + deltat(tjdut)
        # Approximate GMST from JD (same approach as Skyfield)
        T = (jd_tt - 2451545.0) / 36525.0
        gmst_hours = (280.46061837 + 360.98564736629 * (tjdut - 2451545.0) + 0.000387933 * T * T) / 15.0
        init_lon = moon_ra_deg - gmst_hours * 15.0
        init_lon = ((init_lon + 180) % 360) - 180
        init_lat = moon_dec_deg

        def get_separation(lat_arg: float, lon_arg: float) -> float:
            """Get angular separation between Sun and Moon from observer location."""
            try:
                gp = (lon_arg, lat_arg, 0.0)
                jd_tt_loc = tjdut + deltat(tjdut)
                sun_pos = _topo_ecliptic(reader, jd_tt_loc, tjdut, SUN, gp)
                moon_pos = _topo_ecliptic(reader, jd_tt_loc, tjdut, MOON, gp)
                return angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])
            except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
                _reraise_if_leb_range_error(_exc)
                return 999.0

    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(tjdut)

        # Get Earth, Sun, Moon objects
        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        earth_at = earth.at(t)

        # Get geocentric Moon position for initial guess
        moon_geo = earth_at.observe(moon).apparent()
        moon_ra, moon_dec, _ = moon_geo.radec()

        # Calculate sub-lunar point as initial guess
        gmst = t.gmst  # Greenwich Mean Sidereal Time in hours
        init_lon = moon_ra.hours * 15.0 - gmst * 15.0
        init_lon = ((init_lon + 180) % 360) - 180
        init_lat = moon_dec.degrees

        def get_separation(lat_arg: float, lon_arg: float) -> float:
            """Get angular separation between Sun and Moon from observer location."""
            try:
                observer = wgs84.latlon(lat_arg, lon_arg, 0.0)
                observer_at = earth + observer
                sun_app = observer_at.at(t).observe(sun).apparent()
                moon_app = observer_at.at(t).observe(moon).apparent()
                return sun_app.separation_from(moon_app).degrees
            except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
                _reraise_if_leb_range_error(_exc)
                return 999.0  # Return large value on error

    # Function to calculate Sun-Moon separation at a given location
    # (defined above in each branch)

    # Gradient descent to find minimum separation (eclipse center)
    # This finds the point on Earth where Sun and Moon appear closest
    central_lat = init_lat
    central_lon = init_lon

    # Use a multi-scale search: start coarse, then refine
    for scale in [5.0, 1.0, 0.1, 0.01, 0.001]:
        best_sep = get_separation(central_lat, central_lon)

        # Search in a small grid around current best
        improved = True
        iterations = 0
        while improved and iterations < 50:
            improved = False
            iterations += 1

            for dlat, dlon in [
                (scale, 0),
                (-scale, 0),
                (0, scale),
                (0, -scale),
                (scale, scale),
                (scale, -scale),
                (-scale, scale),
                (-scale, -scale),
            ]:
                test_lat = max(-89.0, min(89.0, central_lat + dlat))
                test_lon = central_lon + dlon
                # Normalize longitude
                test_lon = ((test_lon + 180) % 360) - 180

                sep = get_separation(test_lat, test_lon)
                if sep < best_sep:
                    best_sep = sep
                    central_lat = test_lat
                    central_lon = test_lon
                    improved = True

    # Normalize longitude to -180 to +180
    central_lon = ((central_lon + 180) % 360) - 180
    is_central = True

    # Now calculate eclipse attributes at this location
    try:
        if reader is not None:
            _gp = (central_lon, central_lat, 0.0)
            _jd_tt = tjdut + deltat(tjdut)
            sun_topo = _topo_ecliptic(reader, _jd_tt, tjdut, SUN, _gp)
            moon_topo = _topo_ecliptic(reader, _jd_tt, tjdut, MOON, _gp)

            # Get Sun altitude and azimuth at central line
            sun_az_val, sun_alt_true, sun_alt_app_val = azalt(tjdut, ECL2HOR, _gp, 0, 0, sun_topo[:3])
            sun_altitude = sun_alt_true
            sun_azimuth = sun_az_val

            # Calculate apparent altitude with refraction
            if sun_altitude > -1.0:
                if sun_altitude > 0:
                    refraction = 1.0 / math.tan(
                        math.radians(sun_altitude + 7.31 / (sun_altitude + 4.4))
                    )
                    apparent_alt = sun_altitude + refraction / 60.0
                else:
                    apparent_alt = sun_altitude + 0.58
            else:
                apparent_alt = sun_altitude

            separation = angular_separation(sun_topo[0], sun_topo[1], moon_topo[0], moon_topo[1])
            local_sun_dist = sun_topo[2]
            local_moon_dist = moon_topo[2]
        else:
            observer = wgs84.latlon(central_lat, central_lon, 0.0)
            observer_at = earth + observer

            # Get apparent positions from observer (topocentric)
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            # Get Sun altitude and azimuth at central line
            sun_alt, sun_az, _ = sun_app.altaz()
            sun_altitude = sun_alt.degrees
            # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
            sun_azimuth = (sun_az.degrees + 180.0) % 360.0

            # Calculate apparent altitude with refraction
            if sun_altitude > -1.0:
                if sun_altitude > 0:
                    # Bennett's formula for refraction
                    refraction = 1.0 / math.tan(
                        math.radians(sun_altitude + 7.31 / (sun_altitude + 4.4))
                    )
                    apparent_alt = sun_altitude + refraction / 60.0
                else:
                    apparent_alt = sun_altitude + 0.58
            else:
                apparent_alt = sun_altitude

            # Calculate angular separation
            separation = sun_app.separation_from(moon_app).degrees

            # Get local distances
            local_sun_dist = sun_app.distance().au
            local_moon_dist = moon_app.distance().au

        # Local angular radii (in degrees)
        local_sun_radius = (959.63 / 3600.0) / local_sun_dist
        local_moon_radius = (932.56 / 3600.0) * (0.002569 / local_moon_dist)
        local_ratio = local_moon_radius / local_sun_radius

        # Calculate magnitude at central line
        sum_radii = local_sun_radius + local_moon_radius
        if separation >= sum_radii:
            magnitude = 0.0
        elif separation <= abs(local_moon_radius - local_sun_radius):
            # Total or annular
            magnitude = local_moon_radius / local_sun_radius
        else:
            overlap = sum_radii - separation
            magnitude = overlap / (2 * local_sun_radius)
        magnitude = max(0.0, min(1.5, magnitude))

        # Calculate obscuration (area fraction)
        if separation >= sum_radii:
            obscuration = 0.0
        elif separation <= abs(local_moon_radius - local_sun_radius):
            # Both total and annular: obscuration = disc area ratio
            # For total eclipses (moon >= sun), this is >= 1.0
            obscuration = (local_moon_radius / local_sun_radius) ** 2
        else:
            r_sun = local_sun_radius
            r_moon = local_moon_radius
            d = separation
            d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
            d2 = d - d1
            cos_arg1 = max(-1, min(1, d1 / r_sun))
            cos_arg2 = max(-1, min(1, d2 / r_moon))
            area1 = r_sun * r_sun * math.acos(cos_arg1) - d1 * math.sqrt(
                max(0, r_sun * r_sun - d1 * d1)
            )
            area2 = r_moon * r_moon * math.acos(cos_arg2) - d2 * math.sqrt(
                max(0, r_moon * r_moon - d2 * d2)
            )
            obscuration = (area1 + area2) / (math.pi * r_sun * r_sun)
        obscuration = max(0.0, obscuration)

        # Calculate path width (km)
        sun_radius_km = 696340.0
        moon_radius_km = 1737.4
        sun_dist_km = local_sun_dist * AU_TO_KM
        moon_dist_km = local_moon_dist * AU_TO_KM

        # Umbra cone semi-angle
        alpha = math.atan((sun_radius_km - moon_radius_km) / sun_dist_km)

        if local_ratio >= 1.0:
            # Total eclipse - umbra reaches Earth
            umbra_radius_km = moon_radius_km - moon_dist_km * math.tan(alpha)
            umbra_radius_km = max(0, umbra_radius_km)
            is_total_shadow = True
        else:
            # Annular eclipse - antumbra
            umbra_radius_km = moon_dist_km * math.tan(alpha) - moon_radius_km
            umbra_radius_km = max(0, abs(umbra_radius_km))
            is_total_shadow = False

        # Path width affected by Sun altitude
        if sun_altitude > 0:
            sin_alt = math.sin(math.radians(sun_altitude))
            sin_alt = max(0.1, sin_alt)
            path_width_km = 2 * umbra_radius_km / sin_alt
        else:
            path_width_km = 0.0

        path_width_km = min(1000.0, path_width_km)

        # Sign convention: negative for total eclipses (umbra),
        # positive for annular eclipses (antumbra)
        if is_total_shadow:
            path_width_km = -path_width_km

    except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
        _reraise_if_leb_range_error(_exc)
        # If calculation fails, return zeros (10-element geopos, 20-element attr)
        return (
            0,
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
        )

    # Determine eclipse type flags
    if magnitude <= 0:
        return (
            0,
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
            (
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
            ),
        )

    eclipse_type = 0
    if is_central and separation <= abs(local_moon_radius - local_sun_radius):
        eclipse_type |= ECL_CENTRAL
        if local_ratio >= 1.0:
            eclipse_type |= ECL_TOTAL
        elif local_ratio > 0.99:
            eclipse_type |= ECL_ANNULAR_TOTAL
        else:
            eclipse_type |= ECL_ANNULAR
    else:
        eclipse_type |= ECL_PARTIAL

    # Calculate umbra and penumbra limits using Besselian elements
    # This uses the same algorithm as calc_eclipse_northern_limit and
    # calc_eclipse_southern_limit but for a single point in time

    # WGS84 ellipsoid parameters
    EARTH_FLATTENING_CALC = 1.0 / 298.257223563

    # Get Besselian elements at this time
    x_bes = calc_besselian_x(tjdut, flags)
    y_bes = calc_besselian_y(tjdut, flags)
    d_bes = calc_besselian_d(tjdut, flags)
    mu_bes = calc_besselian_mu(tjdut, flags)
    l2_bes = calc_besselian_l2(tjdut, flags)  # umbra/antumbra radius
    l1_bes = calc_besselian_l1(tjdut, flags)  # penumbra radius

    # Calculate gamma - distance from shadow axis to Earth center
    gamma = math.sqrt(x_bes * x_bes + y_bes * y_bes)

    # Initialize limit coordinates to 0.0 (will be filled if applicable)
    umbra_n_lon, umbra_n_lat = 0.0, 0.0
    umbra_s_lon, umbra_s_lat = 0.0, 0.0
    penumbra_n_lon, penumbra_n_lat = 0.0, 0.0
    penumbra_s_lon, penumbra_s_lat = 0.0, 0.0

    def calc_limit_point(
        x_val: float,
        y_val: float,
        shadow_radius: float,
        d_val: float,
        mu_val: float,
        north: bool,
    ) -> Tuple[float, float]:
        """
        Calculate the geographic coordinates of a limit point.

        Args:
            x_val: Besselian x coordinate
            y_val: Besselian y coordinate
            shadow_radius: Shadow radius (l1 for penumbra, |l2| for umbra)
            d_val: Declination of shadow axis in degrees
            mu_val: Hour angle of shadow axis in degrees
            north: True for northern limit, False for southern limit

        Returns:
            Tuple of (longitude, latitude) in degrees, or (0.0, 0.0) if invalid
        """
        # Offset y by the shadow radius (north or south)
        if north:
            y_offset = y_val + shadow_radius
        else:
            y_offset = y_val - shadow_radius

        # Calculate gamma for the limit point
        gamma_limit = math.sqrt(x_val * x_val + y_offset * y_offset)

        # Check if this point is on Earth's surface
        max_gamma_surface = 1.0 + EARTH_FLATTENING_CALC
        if gamma_limit >= max_gamma_surface:
            return 0.0, 0.0

        # Calculate z-factor (height above fundamental plane)
        if gamma_limit > 0.9999:
            z_factor = 0.0
        else:
            z_factor = math.sqrt(max(0, 1.0 - gamma_limit * gamma_limit))

        d_rad = math.radians(d_val)
        sin_d = math.sin(d_rad)
        cos_d = math.cos(d_rad)

        # Calculate latitude using the y_offset component and shadow axis declination
        sin_lat = y_offset * cos_d + z_factor * sin_d

        # Clamp to valid range
        sin_lat = max(-1.0, min(1.0, sin_lat))
        lat = math.degrees(math.asin(sin_lat))

        # For longitude, use the hour angle and x displacement
        if abs(cos_d) > 0.001:
            cos_lat = math.cos(math.radians(lat))
            if cos_lat > 0.001:
                lon_offset = math.degrees(math.atan2(x_val, z_factor * cos_d))
            else:
                lon_offset = 0.0
        else:
            lon_offset = 0.0

        # The longitude of the limit point
        lon = -mu_val + lon_offset

        # Apply correction for Earth's oblateness
        lat_geodetic = lat * (
            1.0 + EARTH_FLATTENING_CALC * math.sin(math.radians(lat)) ** 2
        )

        # Normalize longitude to -180 to +180
        lon = ((lon + 180.0) % 360.0) - 180.0

        return lon, lat_geodetic

    # Calculate umbra limits if the umbral shadow touches Earth
    umbra_radius = abs(l2_bes)
    max_gamma_umbra = 1.0 + EARTH_FLATTENING_CALC + umbra_radius
    if gamma < max_gamma_umbra and umbra_radius > 0.001:
        umbra_n_lon, umbra_n_lat = calc_limit_point(
            x_bes, y_bes, umbra_radius, d_bes, mu_bes, north=True
        )
        umbra_s_lon, umbra_s_lat = calc_limit_point(
            x_bes, y_bes, umbra_radius, d_bes, mu_bes, north=False
        )

    # Calculate penumbra limits if the penumbral shadow touches Earth
    penumbra_radius = abs(l1_bes)
    max_gamma_penumbra = 1.0 + EARTH_FLATTENING_CALC + penumbra_radius
    if gamma < max_gamma_penumbra and penumbra_radius > 0.001:
        penumbra_n_lon, penumbra_n_lat = calc_limit_point(
            x_bes, y_bes, penumbra_radius, d_bes, mu_bes, north=True
        )
        penumbra_s_lon, penumbra_s_lat = calc_limit_point(
            x_bes, y_bes, penumbra_radius, d_bes, mu_bes, north=False
        )

    # Prepare return tuples (reference API format)
    # geopos: 10 floats per reference API documentation
    # [0] longitude central line
    # [1] latitude central line
    # [2] lon northern limit umbra
    # [3] lat northern limit umbra
    # [4] lon southern limit umbra
    # [5] lat southern limit umbra
    # [6] lon northern limit penumbra
    # [7] lat northern limit penumbra
    # [8] lon southern limit penumbra
    # [9] lat southern limit penumbra
    geopos = (
        central_lon,  # [0] longitude central line
        central_lat,  # [1] latitude central line
        umbra_n_lon,  # [2] lon northern limit umbra
        umbra_n_lat,  # [3] lat northern limit umbra
        umbra_s_lon,  # [4] lon southern limit umbra
        umbra_s_lat,  # [5] lat southern limit umbra
        penumbra_n_lon,  # [6] lon northern limit penumbra
        penumbra_n_lat,  # [7] lat northern limit penumbra
        penumbra_s_lon,  # [8] lon southern limit penumbra
        penumbra_s_lat,  # [9] lat southern limit penumbra
    )

    # attr: 20 floats per reference API documentation
    # [0] Fraction of solar diameter covered by Moon (magnitude)
    # [1] Ratio of lunar diameter to solar diameter
    # [2] Fraction of solar disc area obscured (obscuration)
    # [3] Width of totality/annularity path (km)
    # [4] Sun's azimuth at central line (degrees)
    # [5] True altitude of Sun at central line (degrees)
    # [6] Apparent altitude with refraction (degrees)
    # [7] Angular distance Moon center from Sun center (degrees)
    # [8-19] Reserved for future use (zeros)
    attr = (
        magnitude,  # [0] Fraction of solar diameter covered
        local_ratio,  # [1] Moon/Sun diameter ratio
        obscuration,  # [2] Fraction of solar disc area covered
        path_width_km,  # [3] Core shadow width in km
        sun_azimuth,  # [4] Azimuth of sun at maximum
        sun_altitude,  # [5] True altitude of sun
        apparent_alt,  # [6] Apparent altitude with refraction
        separation,  # [7] Angular distance Moon-Sun centers
        local_ratio
        if (eclipse_type & (ECL_TOTAL | ECL_ANNULAR))
        else magnitude,  # [8] Magnitude acc. to NASA
        *_get_saros_info(tjdut, "solar"),  # [9] Saros, [10] member
        0.0,  # [11] Reserved
        0.0,  # [12] Reserved
        0.0,  # [13] Reserved
        0.0,  # [14] Reserved
        0.0,  # [15] Reserved
        0.0,  # [16] Reserved
        0.0,  # [17] Reserved
        0.0,  # [18] Reserved
        0.0,  # [19] Reserved
    )

    return eclipse_type, geopos, attr


def _sol_eclipse_how_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate the circumstances of a solar eclipse at a specific location and time.

    This is the legacy function signature. For reference API compatibility,
    use sol_eclipse_how().

    Args:
        jd: Julian Day (UT) of the moment to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        See sol_eclipse_how() for full specification.
    """
    geopos = (lon, lat, altitude)
    return sol_eclipse_how(jd, geopos, flags)


def sol_eclipse_how(
    tjdut: float,
    geopos: Sequence[float],
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate the circumstances of a solar eclipse at a location.

    Wrapper around :func:`_sol_eclipse_how_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_how_impl, tjdut, geopos, flags,
    )


def _sol_eclipse_how_impl(
    tjdut: float,
    geopos: Sequence[float],
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate the circumstances of a solar eclipse at a specific location and time.

    This function does NOT search for eclipses - it assumes the caller knows
    an eclipse is happening at tjdut and just calculates the circumstances
    at that exact moment from the specified location.

    Args:
        tjdut: Julian Day (UT) of a specific time during an eclipse
        flags: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level

    Returns:
        Tuple containing:
            - attr: Tuple of 20 floats with eclipse attributes per reference API spec:
                [0]: Fraction of solar diameter covered by Moon (magnitude)
                [1]: Ratio of lunar diameter to solar diameter
                [2]: Fraction of solar disc area obscured (obscuration)
                [3]: Diameter of core shadow in km (0 for partial eclipses)
                [4]: Azimuth of Sun (degrees)
                [5]: True altitude of Sun (degrees)
                [6]: Apparent altitude of Sun with refraction (degrees)
                [7]: Elongation of Moon from Sun (degrees)
                [8]: Magnitude according to NASA (currently 0.0, reserved)
                [9]: Saros series number (currently 0.0, reserved)
                [10]: Saros series member number (currently 0.0, reserved)
                [11-19]: Reserved for future use
            - retflag: Eclipse type flags bitmask (ECL_* constants)
                Returns 0 if no eclipse is visible from this location at this time

    Note:
        This function is intended for use when you already know an eclipse is
        occurring (e.g., from sol_eclipse_when_glob or sol_eclipse_when_loc).
        For a random time when no eclipse is occurring, magnitude will be 0
        and retflag will be 0.

    Algorithm:
        1. Calculate Sun and Moon apparent positions from observer location
        2. Compute angular separation between Sun and Moon
        3. Calculate eclipse magnitude from overlap of disks
        4. Determine eclipse type based on whether Moon fully covers Sun
        5. Calculate obscuration (area ratio)
        6. Include atmospheric refraction in altitude calculations

    Precision:
        Eclipse magnitude accurate to ~0.001 for central eclipses.
        Topocentric parallax included in calculations.

    Example:
        >>> # Calculate eclipse circumstances at Dallas during April 8, 2024 eclipse
        >>> from libephemeris import sol_eclipse_how, FLG_SWIEPH
        >>> jd = 2460409.26  # During eclipse maximum
        >>> dallas_geopos = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> ecl_type, attr = sol_eclipse_how(jd, dallas_geopos, FLG_SWIEPH)
        >>> print(f"Obscuration: {attr[2]:.3f}")

    References:
        - Reference API: sol_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # Validate and extract geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        _gp = (lon, lat, altitude)

        try:
            jd_tt = tjdut + deltat(tjdut)
            sun_topo = _topo_ecliptic(reader, jd_tt, tjdut, SUN, _gp)
            moon_topo = _topo_ecliptic(reader, jd_tt, tjdut, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0, (
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            )

        sun_az_val, sun_alt_true, sun_alt_app = azalt(tjdut, ECL2HOR, _gp, 0, 0, sun_topo[:3])
        sun_altitude = sun_alt_true
        sun_azimuth = sun_az_val

        # Calculate apparent altitude with refraction
        if sun_altitude > -1.0:
            if sun_altitude > 0:
                refraction = 1.0 / math.tan(
                    math.radians(sun_altitude + 7.31 / (sun_altitude + 4.4))
                )
                apparent_alt = sun_altitude + refraction / 60.0
            else:
                apparent_alt = sun_altitude + 0.58
        else:
            apparent_alt = sun_altitude

        # If Sun is below horizon, no visible eclipse
        if sun_altitude < -1.0:
            return 0, (
                0.0, 0.0, 0.0, 0.0, sun_azimuth, sun_altitude, apparent_alt,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            )

        separation = angular_separation(sun_topo[0], sun_topo[1], moon_topo[0], moon_topo[1])
        sun_dist_au = sun_topo[2]
        moon_dist_au = moon_topo[2]
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        # Get Skyfield time
        t = ts.ut1_jd(tjdut)

        # Create observer position
        observer_at = earth + observer

        # Get apparent positions from observer
        try:
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # If calculation fails, return zeros (20 elements per reference API spec)
            return 0, (
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            )

        # Get Sun altitude and azimuth
        sun_alt, sun_az, _ = sun_app.altaz()
        sun_altitude = sun_alt.degrees
        # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
        sun_azimuth = (sun_az.degrees + 180.0) % 360.0

        # Calculate apparent altitude with refraction
        if sun_altitude > -1.0:
            if sun_altitude > 0:
                # Bennett's formula for refraction
                refraction = 1.0 / math.tan(
                    math.radians(sun_altitude + 7.31 / (sun_altitude + 4.4))
                )
                apparent_alt = sun_altitude + refraction / 60.0
            else:
                # Near horizon, use approximation
                apparent_alt = sun_altitude + 0.58
        else:
            apparent_alt = sun_altitude

        # If Sun is below horizon, no visible eclipse (20 elements per reference API spec)
        if sun_altitude < -1.0:  # Allow for refraction near horizon
            return 0, (
                0.0, 0.0, 0.0, 0.0, sun_azimuth, sun_altitude, apparent_alt,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            )

        # Calculate angular separation between Sun and Moon
        separation = sun_app.separation_from(moon_app).degrees

        # Get distances for angular size calculations
        sun_dist_au = sun_app.distance().au
        moon_dist_au = moon_app.distance().au

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

    sun_diameter = 2 * sun_angular_radius
    moon_diameter = 2 * moon_angular_radius
    ratio = moon_diameter / sun_diameter

    # Check if there's any eclipse (disks overlapping)
    sum_radii = sun_angular_radius + moon_angular_radius
    if separation >= sum_radii:
        # No eclipse - Sun and Moon too far apart (20 elements per reference API spec)
        return 0, (
            0.0,  # [0] magnitude
            ratio,  # [1] ratio
            0.0,  # [2] obscuration
            0.0,  # [3] shadow width
            sun_azimuth,  # [4] azimuth
            sun_altitude,  # [5] true altitude
            apparent_alt,  # [6] apparent altitude
            separation,  # [7] separation
            0.0,  # [8] magnitude acc. to NASA
            0.0,  # [9] saros series number
            0.0,  # [10] saros series member number
            0.0,  # [11] reserved
            0.0,  # [12] reserved
            0.0,  # [13] reserved
            0.0,  # [14] reserved
            0.0,  # [15] reserved
            0.0,  # [16] reserved
            0.0,  # [17] reserved
            0.0,  # [18] reserved
            0.0,  # [19] reserved
        )

    # Calculate eclipse magnitude
    # Magnitude = fraction of Sun's diameter covered by Moon
    overlap = sum_radii - separation
    magnitude = overlap / sun_diameter
    magnitude = max(0.0, min(magnitude, 1.0 + moon_diameter / sun_diameter))

    # Calculate obscuration (fraction of Sun's area covered)
    r_sun = sun_angular_radius
    r_moon = moon_angular_radius
    d = separation  # center-to-center separation

    if d >= r_sun + r_moon:
        obscuration = 0.0
    elif d <= abs(r_sun - r_moon):
        # One disk entirely within the other
        # Both total and annular: obscuration = disc area ratio
        # For total eclipses (r_moon >= r_sun), this is >= 1.0
        obscuration = (r_moon / r_sun) ** 2
    else:
        # Partial overlap - use lens formula
        # Area of intersection of two circles
        d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
        d2 = d - d1

        if abs(d1) <= r_sun and abs(d2) <= r_moon:
            # Ensure values are valid for acos
            cos_arg1 = max(-1, min(1, d1 / r_sun))
            cos_arg2 = max(-1, min(1, d2 / r_moon))

            area1 = r_sun * r_sun * math.acos(cos_arg1) - d1 * math.sqrt(
                max(0, r_sun * r_sun - d1 * d1)
            )
            area2 = r_moon * r_moon * math.acos(cos_arg2) - d2 * math.sqrt(
                max(0, r_moon * r_moon - d2 * d2)
            )
            intersection_area = area1 + area2
            sun_area = math.pi * r_sun * r_sun
            obscuration = intersection_area / sun_area
        else:
            obscuration = 0.0

    obscuration = max(0.0, obscuration)

    # Calculate core shadow width.
    # Always computed regardless of whether the eclipse is central at the
    # observer — represents the shadow cone geometry projected through the
    # observer's Sun altitude.
    shadow_width_km = 0.0

    AU_TO_KM = 149597870.7
    sun_radius_km = 696340.0
    moon_radius_km = 1737.4
    sun_dist_km = sun_dist_au * AU_TO_KM
    moon_dist_km = moon_dist_au * AU_TO_KM

    # Umbra cone semi-angle
    alpha = math.atan((sun_radius_km - moon_radius_km) / sun_dist_km)

    if ratio >= 1.0:
        # Total eclipse - umbra
        umbra_radius_km = moon_radius_km - moon_dist_km * math.tan(alpha)
        umbra_radius_km = max(0, umbra_radius_km)
        is_total_shadow = True
    else:
        # Annular - antumbra
        umbra_radius_km = moon_dist_km * math.tan(alpha) - moon_radius_km
        umbra_radius_km = max(0, abs(umbra_radius_km))
        is_total_shadow = False

    # Path width affected by Sun altitude
    if sun_altitude > 0:
        sin_alt = math.sin(math.radians(sun_altitude))
        sin_alt = max(0.1, sin_alt)
        shadow_width_km = 2 * umbra_radius_km / sin_alt
    else:
        shadow_width_km = 0.0

    shadow_width_km = min(1000.0, shadow_width_km)

    # Sign convention: negative for total eclipses (umbra),
    # positive for annular eclipses (antumbra)
    if is_total_shadow:
        shadow_width_km = -shadow_width_km

    # Determine eclipse type flags
    diff_radii = abs(moon_angular_radius - sun_angular_radius)
    eclipse_type = ECL_VISIBLE
    moon_sun_ratio = moon_angular_radius / sun_angular_radius

    if separation < diff_radii:
        # Central eclipse at this location
        eclipse_type |= ECL_CENTRAL
        if moon_sun_ratio >= 1.0:
            eclipse_type |= ECL_TOTAL
        elif moon_sun_ratio > 0.99:
            eclipse_type |= ECL_ANNULAR_TOTAL
        else:
            eclipse_type |= ECL_ANNULAR
    else:
        # Partial eclipse
        eclipse_type |= ECL_PARTIAL

    # Prepare attributes tuple (20 elements matching reference API format)
    attr = (
        magnitude,  # [0] Fraction of solar diameter covered
        ratio,  # [1] Ratio of lunar to solar diameter
        obscuration,  # [2] Fraction of solar disc area covered
        shadow_width_km,  # [3] Core shadow width in km
        sun_azimuth,  # [4] Azimuth of sun
        sun_altitude,  # [5] True altitude of sun
        apparent_alt,  # [6] Apparent altitude with refraction
        separation,  # [7] Angular distance Moon-Sun centers
        ratio
        if (eclipse_type & (ECL_TOTAL | ECL_ANNULAR))
        else magnitude,  # [8] Magnitude acc. to NASA
        *_get_saros_info(tjdut, "solar"),  # [9] Saros, [10] member
        0.0,  # [11] reserved
        0.0,  # [12] reserved
        0.0,  # [13] reserved
        0.0,  # [14] reserved
        0.0,  # [15] reserved
        0.0,  # [16] reserved
        0.0,  # [17] reserved
        0.0,  # [18] reserved
        0.0,  # [19] reserved
    )

    return eclipse_type, attr


def sol_eclipse_how_details(
    tjd_ut: float,
    geopos: Sequence[float],
    ifl: int = 0,
) -> dict:
    """Calculate comprehensive solar eclipse circumstances at a location.

    Wrapper around :func:`_sol_eclipse_how_details_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_how_details_impl, tjd_ut, geopos, ifl,
    )


def _sol_eclipse_how_details_impl(
    tjd_ut: float,
    geopos: Sequence[float],
    ifl: int = 0,
    *,
    reader=None,
) -> dict:
    """
    Calculate comprehensive solar eclipse circumstances at a specific location.

    This enhanced version of sol_eclipse_how returns all eclipse details
    including contact times, maximum obscuration, position angles, and Sun
    altitude/azimuth during all eclipse phases.

    Args:
        tjd_ut: Julian Day (UT) of a specific time during an eclipse.
                This should be a time when an eclipse is known to occur
                (e.g., from sol_eclipse_when_glob or sol_eclipse_when_loc).
        ifl: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level

    Returns:
        Dictionary with comprehensive eclipse information:
            'eclipse_type': int - Eclipse type flags bitmask (ECL_* constants)
            'is_visible': bool - Whether eclipse is visible from this location
            'is_total': bool - Whether eclipse is total at this location
            'is_annular': bool - Whether eclipse is annular at this location
            'is_partial': bool - Whether eclipse is partial at this location

            Contact times (Julian Day UT, or 0.0 if not applicable):
            'jd_c1': float - First contact (partial phase begins)
            'jd_c2': float - Second contact (totality/annularity begins)
            'jd_max': float - Maximum eclipse
            'jd_c3': float - Third contact (totality/annularity ends)
            'jd_c4': float - Fourth contact (partial phase ends)

            Eclipse magnitude and obscuration:
            'magnitude': float - Fraction of solar diameter covered by Moon
            'max_magnitude': float - Maximum magnitude during eclipse
            'obscuration': float - Fraction of solar disc area covered (0-1)
            'max_obscuration': float - Maximum obscuration during eclipse (0-1)
            'obscuration_percent': float - Obscuration as percentage (0-100)
            'max_obscuration_percent': float - Maximum obscuration as percentage
            'ratio': float - Ratio of lunar diameter to solar diameter
            'shadow_width_km': float - Core shadow width in km (0 for partial)

            Position angles (degrees from North, clockwise):
            'position_angle_c1': float - Position angle at first contact
            'position_angle_c2': float - Position angle at second contact (or 0)
            'position_angle_c3': float - Position angle at third contact (or 0)
            'position_angle_c4': float - Position angle at fourth contact

            Sun position at each phase:
            'sun_alt_c1': float - Sun altitude at first contact
            'sun_az_c1': float - Sun azimuth at first contact
            'sun_alt_c2': float - Sun altitude at second contact
            'sun_az_c2': float - Sun azimuth at second contact
            'sun_alt_max': float - Sun altitude at maximum eclipse
            'sun_az_max': float - Sun azimuth at maximum eclipse
            'sun_alt_c3': float - Sun altitude at third contact
            'sun_az_c3': float - Sun azimuth at third contact
            'sun_alt_c4': float - Sun altitude at fourth contact
            'sun_az_c4': float - Sun azimuth at fourth contact

            Duration information:
            'duration_partial_minutes': float - Total duration of partial phases
            'duration_total_minutes': float - Duration of totality (0 if not total)

            Angular sizes:
            'sun_angular_radius': float - Sun angular radius in degrees
            'moon_angular_radius': float - Moon angular radius in degrees
            'separation': float - Angular separation Sun-Moon at given time

    Algorithm:
        1. Find eclipse maximum for this location using golden section search
        2. Calculate contact times using bisection method
        3. Compute position angles from celestial mechanics
        4. Calculate Sun altitude/azimuth at each phase

    Precision:
        Contact times accurate to ~1 second.
        Position angles accurate to ~0.1 degrees.

    Example:
        >>> from libephemeris import sol_eclipse_how_details, FLG_SWIEPH
        >>> jd = 2460409.28  # During April 8, 2024 eclipse
        >>> dallas = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> details = sol_eclipse_how_details(jd, dallas, FLG_SWIEPH)
        >>> print(f"Max obscuration: {details['max_obscuration_percent']:.1f}%")
        >>> print(f"First contact: JD {details['jd_c1']:.5f}")
        >>> print(f"Duration of totality: {details['duration_total_minutes']:.1f} min")

    References:
        - Reference documentation
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # Validate and extract geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        _gp = (lon, lat, altitude)

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
            return angular_separation(sun_pos[0], sun_pos[1], moon_pos[0], moon_pos[1])

        def _get_sun_altaz(jd: float) -> tuple:
            """Get Sun altitude and azimuth at given JD."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            sun_az, sun_alt_true, sun_alt_app = azalt(jd, ECL2HOR, _gp, 0, 0, sun_pos[:3])
            return sun_alt_true, sun_az

        def _get_angular_sizes(jd: float) -> tuple:
            """Get angular radii of Sun and Moon."""
            jd_tt = jd + deltat(jd)
            sun_pos = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
            sun_angular_radius = (959.63 / 3600.0) / sun_pos[2]
            moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_pos[2])
            return sun_angular_radius, moon_angular_radius

        def _calc_position_angle(jd: float) -> float:
            """Calculate position angle of Moon relative to Sun (topocentric ICRS)."""
            from .constants import FLG_EQUATORIAL as _EQ_PA, FLG_J2000 as _J2K_PA
            from .constants import FLG_TOPOCTR as _TP_PA
            from .state import get_topo as _gt_pa
            import libephemeris as _le_pa
            _saved_topo_pa = _gt_pa()
            _le_pa.set_topo(_gp[0], _gp[1], _gp[2])
            try:
                _pa_flags = _EQ_PA | _J2K_PA | _TP_PA | FLG_SPEED
                sun_eq, _ = calc_ut(jd, SUN, _pa_flags)
                moon_eq, _ = calc_ut(jd, MOON, _pa_flags)
            finally:
                from libephemeris import state as _st_pa
                if _saved_topo_pa is not None:
                    _le_pa.set_topo(_saved_topo_pa.longitude.degrees, _saved_topo_pa.latitude.degrees, _saved_topo_pa.elevation.m)
                else:
                    _st_pa._TOPO = None
            sun_ra_rad = math.radians(sun_eq[0])
            sun_dec_rad = math.radians(sun_eq[1])
            moon_ra_rad = math.radians(moon_eq[0])
            moon_dec_rad = math.radians(moon_eq[1])
            delta_ra = moon_ra_rad - sun_ra_rad
            y = math.sin(delta_ra)
            x = math.cos(sun_dec_rad) * math.tan(moon_dec_rad) - math.sin(
                sun_dec_rad
            ) * math.cos(delta_ra)
            pa = math.degrees(math.atan2(y, x))
            if pa < 0:
                pa += 360.0
            return pa
    else:
        from skyfield.api import wgs84

        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        def _get_sun_moon_positions(jd: float):
            """Get Sun and Moon positions from observer."""
            t = ts.ut1_jd(jd)
            observer_at = earth + observer

            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()

            return sun_app, moon_app

        def _get_separation(jd: float) -> float:
            """Get angular separation between Sun and Moon."""
            sun_app, moon_app = _get_sun_moon_positions(jd)
            return sun_app.separation_from(moon_app).degrees

        def _get_sun_altaz(jd: float) -> tuple:
            """Get Sun altitude and azimuth at given JD."""
            sun_app, _ = _get_sun_moon_positions(jd)
            alt, az, _ = sun_app.altaz()
            # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_angular_sizes(jd: float) -> tuple:
            """Get angular radii of Sun and Moon."""
            sun_app, moon_app = _get_sun_moon_positions(jd)

            sun_dist_au = sun_app.distance().au
            moon_dist_au = moon_app.distance().au

            # Angular radii in degrees
            sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
            moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

            return sun_angular_radius, moon_angular_radius

        def _calc_position_angle(jd: float) -> float:
            """
            Calculate position angle of Moon relative to Sun.

            Position angle is measured from North through East (counterclockwise).
            Returns angle in degrees (0-360).
            """
            sun_app, moon_app = _get_sun_moon_positions(jd)

            # Get RA and Dec for both bodies
            sun_ra, sun_dec, _ = sun_app.radec()
            moon_ra, moon_dec, _ = moon_app.radec()

            sun_ra_rad = sun_ra.radians
            sun_dec_rad = sun_dec.radians
            moon_ra_rad = moon_ra.radians
            moon_dec_rad = moon_dec.radians

            # Position angle formula
            delta_ra = moon_ra_rad - sun_ra_rad

            y = math.sin(delta_ra)
            x = math.cos(sun_dec_rad) * math.tan(moon_dec_rad) - math.sin(
                sun_dec_rad
            ) * math.cos(delta_ra)

            pa = math.degrees(math.atan2(y, x))
            if pa < 0:
                pa += 360.0

            return pa

    def _find_local_maximum(
        jd_center: float, search_range: float = 3.0 / 24.0
    ) -> float:
        """Find local eclipse maximum (minimum separation) using golden section."""
        phi = (1 + math.sqrt(5)) / 2

        jd_low = jd_center - search_range
        jd_high = jd_center + search_range

        jd_a = jd_high - (jd_high - jd_low) / phi
        jd_b = jd_low + (jd_high - jd_low) / phi

        sep_a = _get_separation(jd_a)
        sep_b = _get_separation(jd_b)

        for _ in range(50):
            if sep_a < sep_b:
                jd_high = jd_b
                jd_b = jd_a
                sep_b = sep_a
                jd_a = jd_high - (jd_high - jd_low) / phi
                sep_a = _get_separation(jd_a)
            else:
                jd_low = jd_a
                jd_a = jd_b
                sep_a = sep_b
                jd_b = jd_low + (jd_high - jd_low) / phi
                sep_b = _get_separation(jd_b)

            if jd_high - jd_low < 1e-8:
                break

        return (jd_low + jd_high) / 2

    def _find_contact_time(
        jd_start: float,
        jd_end: float,
        target_sep: float,
        is_increasing: bool,
    ) -> float:
        """Find time when separation equals target using bisection."""
        jd_low = jd_start
        jd_high = jd_end

        for _ in range(60):
            jd_mid = (jd_low + jd_high) / 2
            sep = _get_separation(jd_mid)

            if abs(sep - target_sep) < 1e-7:
                return jd_mid

            if is_increasing:
                if sep < target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid
            else:
                if sep > target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid

            if jd_high - jd_low < 1e-9:
                break

        return (jd_low + jd_high) / 2

    # First, get basic eclipse info at the given time
    eclipse_type, attr = sol_eclipse_how(tjd_ut, geopos, ifl)

    # Initialize result dictionary
    result = {
        "eclipse_type": eclipse_type,
        "is_visible": bool(eclipse_type & ECL_VISIBLE),
        "is_total": bool(eclipse_type & ECL_TOTAL),
        "is_annular": bool(eclipse_type & ECL_ANNULAR),
        "is_partial": bool(eclipse_type & ECL_PARTIAL),
        # Contact times
        "jd_c1": 0.0,
        "jd_c2": 0.0,
        "jd_max": 0.0,
        "jd_c3": 0.0,
        "jd_c4": 0.0,
        # Magnitude and obscuration
        "magnitude": attr[0],
        "max_magnitude": 0.0,
        "obscuration": attr[2],
        "max_obscuration": 0.0,
        "obscuration_percent": attr[2] * 100.0,
        "max_obscuration_percent": 0.0,
        "ratio": attr[1],
        "shadow_width_km": attr[3],
        # Position angles
        "position_angle_c1": 0.0,
        "position_angle_c2": 0.0,
        "position_angle_c3": 0.0,
        "position_angle_c4": 0.0,
        # Sun positions
        "sun_alt_c1": 0.0,
        "sun_az_c1": 0.0,
        "sun_alt_c2": 0.0,
        "sun_az_c2": 0.0,
        "sun_alt_max": attr[5],
        "sun_az_max": attr[4],
        "sun_alt_c3": 0.0,
        "sun_az_c3": 0.0,
        "sun_alt_c4": 0.0,
        "sun_az_c4": 0.0,
        # Durations
        "duration_partial_minutes": 0.0,
        "duration_total_minutes": 0.0,
        # Angular sizes
        "sun_angular_radius": 0.0,
        "moon_angular_radius": 0.0,
        "separation": attr[7],
    }

    # If no eclipse is visible, return early
    if eclipse_type == 0:
        return result

    # Get angular sizes
    sun_r, moon_r = _get_angular_sizes(tjd_ut)
    result["sun_angular_radius"] = sun_r
    result["moon_angular_radius"] = moon_r

    sum_radii = sun_r + moon_r
    diff_radii = abs(sun_r - moon_r)

    # Find local maximum for this observer
    jd_local_max = _find_local_maximum(tjd_ut)
    result["jd_max"] = jd_local_max

    # Get maximum eclipse circumstances
    max_sep = _get_separation(jd_local_max)
    sun_alt_max, sun_az_max = _get_sun_altaz(jd_local_max)
    result["sun_alt_max"] = sun_alt_max
    result["sun_az_max"] = sun_az_max

    # Calculate maximum magnitude and obscuration
    sun_r_max, moon_r_max = _get_angular_sizes(jd_local_max)
    sun_diameter = 2 * sun_r_max
    moon_diameter = 2 * moon_r_max

    sum_radii_max = sun_r_max + moon_r_max
    overlap = sum_radii_max - max_sep
    max_magnitude = max(
        0.0, min(overlap / sun_diameter, 1.0 + moon_diameter / sun_diameter)
    )
    result["max_magnitude"] = max_magnitude
    result["ratio"] = moon_diameter / sun_diameter

    # Calculate maximum obscuration
    d = max_sep
    if d >= sun_r_max + moon_r_max:
        max_obscuration = 0.0
    elif d <= abs(sun_r_max - moon_r_max):
        # Both total and annular: obscuration = disc area ratio
        # For total eclipses (moon >= sun), this is >= 1.0
        max_obscuration = (moon_r_max / sun_r_max) ** 2
    else:
        d1 = (d * d + sun_r_max * sun_r_max - moon_r_max * moon_r_max) / (2 * d)
        d2 = d - d1

        if abs(d1) <= sun_r_max and abs(d2) <= moon_r_max:
            cos_arg1 = max(-1, min(1, d1 / sun_r_max))
            cos_arg2 = max(-1, min(1, d2 / moon_r_max))

            area1 = sun_r_max * sun_r_max * math.acos(cos_arg1) - d1 * math.sqrt(
                max(0, sun_r_max * sun_r_max - d1 * d1)
            )
            area2 = moon_r_max * moon_r_max * math.acos(cos_arg2) - d2 * math.sqrt(
                max(0, moon_r_max * moon_r_max - d2 * d2)
            )
            intersection_area = area1 + area2
            sun_area = math.pi * sun_r_max * sun_r_max
            max_obscuration = intersection_area / sun_area
        else:
            max_obscuration = 0.0

    max_obscuration = max(0.0, max_obscuration)
    result["max_obscuration"] = max_obscuration
    result["max_obscuration_percent"] = max_obscuration * 100.0

    # Update eclipse type based on maximum
    is_central = max_sep < diff_radii
    if is_central:
        result["is_partial"] = False
        if moon_r_max >= sun_r_max:
            result["is_total"] = True
            result["is_annular"] = False
        else:
            result["is_total"] = False
            result["is_annular"] = True
    else:
        result["is_partial"] = True
        result["is_total"] = False
        result["is_annular"] = False

    # Calculate contact times
    contact_search_range = 2.5 / 24.0  # 2.5 hours

    # First contact (separation decreasing to sum of radii)
    try:
        jd_c1 = _find_contact_time(
            jd_local_max - contact_search_range,
            jd_local_max,
            sum_radii_max,
            is_increasing=False,
        )
        result["jd_c1"] = jd_c1

        # Sun position at C1
        sun_alt_c1, sun_az_c1 = _get_sun_altaz(jd_c1)
        result["sun_alt_c1"] = sun_alt_c1
        result["sun_az_c1"] = sun_az_c1

        # Position angle at C1
        result["position_angle_c1"] = _calc_position_angle(jd_c1)
    except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
        _reraise_if_leb_range_error(_exc)
        pass

    # Fourth contact (separation increasing from sum of radii)
    try:
        jd_c4 = _find_contact_time(
            jd_local_max,
            jd_local_max + contact_search_range,
            sum_radii_max,
            is_increasing=True,
        )
        result["jd_c4"] = jd_c4

        # Sun position at C4
        sun_alt_c4, sun_az_c4 = _get_sun_altaz(jd_c4)
        result["sun_alt_c4"] = sun_alt_c4
        result["sun_az_c4"] = sun_az_c4

        # Position angle at C4
        result["position_angle_c4"] = _calc_position_angle(jd_c4)
    except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
        _reraise_if_leb_range_error(_exc)
        pass

    # Second and third contacts (for central eclipses only)
    diff_radii_max = abs(sun_r_max - moon_r_max)
    if is_central and max_sep < diff_radii_max:
        # Second contact
        try:
            jd_c2 = _find_contact_time(
                jd_local_max - contact_search_range / 4,
                jd_local_max,
                diff_radii_max,
                is_increasing=False,
            )
            result["jd_c2"] = jd_c2

            # Sun position at C2
            sun_alt_c2, sun_az_c2 = _get_sun_altaz(jd_c2)
            result["sun_alt_c2"] = sun_alt_c2
            result["sun_az_c2"] = sun_az_c2

            # Position angle at C2
            result["position_angle_c2"] = _calc_position_angle(jd_c2)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            pass

        # Third contact
        try:
            jd_c3 = _find_contact_time(
                jd_local_max,
                jd_local_max + contact_search_range / 4,
                diff_radii_max,
                is_increasing=True,
            )
            result["jd_c3"] = jd_c3

            # Sun position at C3
            sun_alt_c3, sun_az_c3 = _get_sun_altaz(jd_c3)
            result["sun_alt_c3"] = sun_alt_c3
            result["sun_az_c3"] = sun_az_c3

            # Position angle at C3
            result["position_angle_c3"] = _calc_position_angle(jd_c3)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            pass

    # Calculate durations
    if result["jd_c1"] > 0 and result["jd_c4"] > 0:
        result["duration_partial_minutes"] = (
            (result["jd_c4"] - result["jd_c1"]) * 24 * 60
        )

    if result["jd_c2"] > 0 and result["jd_c3"] > 0:
        result["duration_total_minutes"] = (result["jd_c3"] - result["jd_c2"]) * 24 * 60

    return result


def _sol_eclipse_how_details_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> dict:
    """
    Calculate comprehensive solar eclipse circumstances at a specific location.

    This is the legacy function signature. For reference API-style arguments,
    use sol_eclipse_how_details().

    Args:
        jd: Julian Day (UT) of the moment to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Dictionary with comprehensive eclipse information.
        See sol_eclipse_how_details() for full specification.
    """
    geopos = (lon, lat, altitude)
    return sol_eclipse_how_details(jd, geopos, flags)


# =============================================================================
# LUNAR ECLIPSE FUNCTIONS
# =============================================================================

# Constants for lunar eclipse calculations
ECLIPSE_LIMIT_LUNAR = 18.5  # Maximum elongation from node for lunar eclipse (degrees)
# Note: penumbral eclipses can occur up to ~18° from a node (observed max 18.02°).
# The previous value of 12.0° was too restrictive and missed shallow penumbral eclipses.


def _find_next_full_moon(jd_start: float) -> float:
    """
    Find the next Full Moon (Sun-Moon opposition) after jd_start.

    Uses iterative refinement to find exact moment of opposition.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of next Full Moon
    """
    # Get current positions
    sun_pos, _ = calc_ut(jd_start, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_start, MOON, FLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation from opposition (Moon - Sun - 180°)
    elongation = (moon_lon - sun_lon - 180.0) % 360.0
    if elongation > 180:
        elongation -= 360

    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day

    # Time until next opposition (elongation from 180° = 0)
    if elongation > 0:
        # Past opposition, wait for next cycle
        dt = (360.0 - elongation) / relative_speed
    else:
        dt = (-elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = calc_ut(jd_guess, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd_guess, MOON, FLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation from opposition
        diff = (moon_lon - sun_lon - 180.0) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _find_previous_full_moon(jd_start: float) -> float:
    """
    Find the previous Full Moon (Sun-Moon opposition) before jd_start.

    Uses iterative refinement to find exact moment of opposition.

    Args:
        jd_start: Julian Day (UT) to start search from

    Returns:
        Julian Day of previous Full Moon
    """
    # Get current positions
    sun_pos, _ = calc_ut(jd_start, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_start, MOON, FLG_SPEED)

    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]

    # Calculate elongation from opposition (Moon - Sun - 180°)
    elongation = (moon_lon - sun_lon - 180.0) % 360.0
    if elongation > 180:
        elongation -= 360

    # Moon gains ~12.2° per day on Sun
    relative_speed = 12.190749  # degrees/day

    # Time since last opposition (go backwards)
    if elongation > 0:
        # Past opposition, last opposition was elongation/speed days ago
        dt = -elongation / relative_speed
    else:
        # Before opposition, last opposition was (360 + elongation)/speed days ago
        dt = -(360.0 + elongation) / relative_speed

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(20):
        sun_pos, _ = calc_ut(jd_guess, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd_guess, MOON, FLG_SPEED)

        sun_lon = sun_pos[0]
        moon_lon = moon_pos[0]
        sun_speed = sun_pos[3]
        moon_speed = moon_pos[3]

        # Elongation from opposition
        diff = (moon_lon - sun_lon - 180.0) % 360.0
        if diff > 180:
            diff -= 360

        # Convergence check (< 0.1 arcsec)
        if abs(diff) < 1e-5:
            return jd_guess

        # Newton-Raphson step
        rel_speed = moon_speed - sun_speed
        if abs(rel_speed) < 0.1:
            rel_speed = 12.19

        jd_guess -= diff / rel_speed

    return jd_guess


def _calculate_lunar_eclipse_type_and_magnitude(
    jd: float,
) -> Tuple[int, float, float, float, float, float]:
    """
    Determine lunar eclipse type and magnitude at maximum eclipse.

    Uses geometric calculations based on Earth's shadow cones at the Moon's distance.
    Includes edge case handling for:
    - Very shallow penumbral eclipses (penumbral magnitude near 0)
    - Very shallow umbral eclipses (umbral magnitude near 0)
    - Division by zero protection

    Args:
        jd: Julian Day of eclipse maximum (UT)

    Returns:
        Tuple of (eclipse_type_flags, magnitude_umbral, magnitude_penumbral,
                  gamma, penumbra_radius, umbra_radius)
        - eclipse_type_flags: Bitmask of ECL_* constants
        - magnitude_umbral: Eclipse magnitude (fraction of Moon in umbra)
        - magnitude_penumbral: Penumbral eclipse magnitude
        - gamma: Gamma parameter (Moon's distance from shadow axis in Earth radii)
        - penumbra_radius: Penumbral shadow radius at Moon distance (degrees)
        - umbra_radius: Umbral shadow radius at Moon distance (degrees)
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Distances in AU
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Moon's ecliptic latitude (how far from the ecliptic)
    moon_lat = moon_pos[1]

    # Constants for shadow calculations
    # Sun angular radius at 1 AU: 959.63 arcsec
    # Earth radius: 6378.137 km = 4.2635e-5 AU
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_AU = 4.2635e-5
    EARTH_RADIUS_KM = 6378.137

    # Sun's angular semi-diameter at actual distance (in degrees)
    sun_semidiameter = (SUN_RADIUS_ARCSEC / 3600.0) / sun_dist

    # Moon's angular semi-diameter (932.56 arcsec at mean distance 0.002569 AU)
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Earth's angular semi-diameter as seen from Moon (in degrees)
    earth_semidiameter = math.degrees(math.atan(EARTH_RADIUS_AU / moon_dist))

    # Calculate shadow cone sizes at Moon's distance
    # Using the geometric relationship of similar triangles

    # Sun's angular radius as seen from Earth
    sun_angular_radius_from_earth = sun_semidiameter

    # Penumbra radius at Moon's distance (Earth + Sun shadow)
    # Penumbra outer edge: Earth appears larger, adding Sun's angular size
    penumbra_radius = earth_semidiameter + sun_angular_radius_from_earth

    # Umbra radius at Moon's distance (Earth minus Sun shadow)
    # The umbra is smaller because of the Sun's finite size
    # Using more accurate shadow cone calculation
    # Shadow cone angle from geometry
    sun_dist_km = sun_dist * 149597870.7
    moon_dist_km = moon_dist * 149597870.7
    sun_radius_km = (SUN_RADIUS_ARCSEC / 206265.0) * sun_dist_km

    # Umbra cone semi-angle - protect against division by zero
    denom = sun_dist_km
    if abs(denom) < 1e-10:
        umbra_cone_angle = 0.0
    else:
        umbra_cone_angle = math.atan((sun_radius_km - EARTH_RADIUS_KM) / denom)
    # Umbra radius at Moon distance
    umbra_radius_km = EARTH_RADIUS_KM - moon_dist_km * math.tan(umbra_cone_angle)

    if umbra_radius_km > 0:
        umbra_radius = math.degrees(math.atan(umbra_radius_km / moon_dist_km))
    else:
        # Umbra doesn't reach Moon (antumbra) - extremely rare for lunar eclipses
        umbra_radius = 0.0

    # Apply atmospheric enlargement of Earth's shadow (Danjon factor).
    # Earth's atmosphere refracts sunlight, effectively enlarging the shadow
    # by approximately 1/85 (~1.18%). This correction is standard in eclipse
    # magnitude calculations and is necessary for correct classification of
    # borderline total/partial eclipses.
    # Reference: Danjon, A. (1951), "Les éclipses de Lune"
    # Also: Meeus, "Astronomical Algorithms", Ch. 54
    _SHADOW_ENLARGEMENT = 1.0 + 1.0 / 85.0
    penumbra_radius *= _SHADOW_ENLARGEMENT
    umbra_radius *= _SHADOW_ENLARGEMENT

    # Moon's distance from the shadow axis (in degrees).
    # The shadow centre is at the anti-solar point (Sun + 180° in longitude,
    # negated latitude).  At eclipse maximum the longitude term is nearly
    # zero and |moon_lat| alone is sufficient, but at non-maximum times
    # (e.g. moonrise/moonset) the longitude offset can be significant
    # and must be included.
    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]
    shadow_lon = (sun_lon + 180.0) % 360.0

    delta_lon = moon_lon - shadow_lon
    if delta_lon > 180.0:
        delta_lon -= 360.0
    elif delta_lon < -180.0:
        delta_lon += 360.0

    # Small-angle approximation for angular separation on the ecliptic
    # (cos correction for longitude at the Moon's latitude)
    cos_lat = math.cos(math.radians(moon_lat))
    moon_distance_from_axis = math.sqrt((delta_lon * cos_lat) ** 2 + moon_lat**2)

    # Gamma: Moon's distance from shadow axis in Earth radii
    # Edge case: protect against zero earth_semidiameter
    if abs(earth_semidiameter) < 1e-10:
        gamma = 0.0
    else:
        gamma = moon_distance_from_axis / earth_semidiameter

    # Calculate eclipse magnitudes
    # Edge case: protect against zero moon_semidiameter
    if abs(moon_semidiameter) < 1e-10:
        return 0, 0.0, 0.0, gamma, penumbra_radius, umbra_radius

    # Penumbral magnitude: how far Moon penetrates into penumbra
    penumbral_mag = (penumbra_radius + moon_semidiameter - moon_distance_from_axis) / (
        2 * moon_semidiameter
    )

    # Umbral magnitude: how far Moon penetrates into umbra
    umbral_mag = (umbra_radius + moon_semidiameter - moon_distance_from_axis) / (
        2 * moon_semidiameter
    )

    # Determine eclipse type
    eclipse_type = 0

    # Edge case: no eclipse (penumbral magnitude <= 0)
    if penumbral_mag <= 0:
        # No eclipse
        return 0, 0.0, 0.0, gamma, penumbra_radius, umbra_radius

    # Edge case: shallow penumbral eclipse
    if penumbral_mag > 0 and penumbral_mag < SHALLOW_ECLIPSE_MAG_THRESHOLD:
        # Very shallow penumbral eclipse - mark as grazing
        # Shallow grazing penumbral: the reference reports plain PENUMBRAL
        eclipse_type = ECL_PENUMBRAL
        penumbral_mag = max(0.0, penumbral_mag)
        return eclipse_type, 0.0, penumbral_mag, gamma, penumbra_radius, umbra_radius

    if umbral_mag <= 0:
        # Penumbral only
        eclipse_type = ECL_PENUMBRAL
        penumbral_mag = max(0.0, penumbral_mag)
        return eclipse_type, 0.0, penumbral_mag, gamma, penumbra_radius, umbra_radius

    # Edge case: shallow umbral (partial) eclipse
    if umbral_mag > 0 and umbral_mag < SHALLOW_ECLIPSE_MAG_THRESHOLD:
        # Very shallow partial umbral eclipse - mark as grazing
        # Shallow grazing partial: the reference reports plain PARTIAL
        eclipse_type = ECL_PARTIAL
        umbral_mag = max(0.0, min(1.0, umbral_mag))
        penumbral_mag = max(0.0, penumbral_mag)
        return (
            eclipse_type,
            umbral_mag,
            penumbral_mag,
            gamma,
            penumbra_radius,
            umbra_radius,
        )

    if umbral_mag >= 1.0:
        # Total umbral eclipse
        eclipse_type = ECL_TOTAL
        umbral_mag = max(0.0, umbral_mag)
        penumbral_mag = max(0.0, penumbral_mag)
    else:
        # Partial umbral eclipse
        eclipse_type = ECL_PARTIAL
        umbral_mag = max(0.0, min(1.0, umbral_mag))
        penumbral_mag = max(0.0, penumbral_mag)

    return eclipse_type, umbral_mag, penumbral_mag, gamma, penumbra_radius, umbra_radius


def _refine_lunar_eclipse_maximum(
    jd_full_moon: float, search_range: float = 0.5
) -> float:
    """
    Refine the lunar eclipse maximum time from Full Moon approximation.

    The eclipse maximum occurs when the Moon is closest to Earth's shadow
    center (the anti-Sun point), not exactly at Full Moon (opposition in
    longitude). This function uses golden section search to find the true
    maximum by minimizing the angular separation.

    Args:
        jd_full_moon: Julian Day of Full Moon (initial approximation)
        search_range: Search range in days (±0.5 days covers the offset)

    Returns:
        Julian Day of true eclipse maximum (when Moon is closest to shadow axis)
    """

    def _get_shadow_separation(jd: float) -> float:
        """Calculate angular separation between Moon and shadow center."""
        sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
        moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

        # Shadow center is the anti-Sun point (180° from Sun)
        shadow_lon = (sun_pos[0] + 180.0) % 360.0
        shadow_lat = -sun_pos[1]  # Opposite latitude

        # Calculate angular separation using spherical law of cosines
        d_lon = math.radians(moon_pos[0] - shadow_lon)
        lat1 = math.radians(shadow_lat)
        lat2 = math.radians(moon_pos[1])

        cos_sep = math.sin(lat1) * math.sin(lat2) + math.cos(lat1) * math.cos(
            lat2
        ) * math.cos(d_lon)
        cos_sep = max(-1.0, min(1.0, cos_sep))
        return math.degrees(math.acos(cos_sep))

    # Golden section search to minimize separation
    phi = (1 + math.sqrt(5)) / 2

    jd_low = jd_full_moon - search_range
    jd_high = jd_full_moon + search_range

    jd_a = jd_high - (jd_high - jd_low) / phi
    jd_b = jd_low + (jd_high - jd_low) / phi

    sep_a = _get_shadow_separation(jd_a)
    sep_b = _get_shadow_separation(jd_b)

    for _ in range(50):
        if sep_a < sep_b:
            jd_high = jd_b
            jd_b = jd_a
            sep_b = sep_a
            jd_a = jd_high - (jd_high - jd_low) / phi
            sep_a = _get_shadow_separation(jd_a)
        else:
            jd_low = jd_a
            jd_a = jd_b
            sep_a = sep_b
            jd_b = jd_low + (jd_high - jd_low) / phi
            sep_b = _get_shadow_separation(jd_b)

        if jd_high - jd_low < 1e-8:  # ~0.86 ms precision
            break

    return (jd_low + jd_high) / 2


def _calculate_lunar_eclipse_phases(
    jd_max: float,
    eclipse_type: int,
    moon_semidiameter: float,
    umbra_radius: float,
    penumbra_radius: float,
    moon_lat_at_max: float,
) -> Tuple[float, float, float, float, float, float, float, float]:
    """
    Calculate times of lunar eclipse phases (contacts).

    Phase indices (pyswisseph-compatible 10-element format):
        [0]: Time of maximum eclipse
        [1]: Reserved (0)
        [2]: Time of partial eclipse beginning (Moon enters umbra)
        [3]: Time of partial eclipse ending (Moon leaves umbra)
        [4]: Time of total eclipse beginning (or 0 if not total)
        [5]: Time of total eclipse ending (or 0 if not total)
        [6]: Time of penumbral eclipse beginning
        [7]: Time of penumbral eclipse ending
        [8]: Reserved (0)
        [9]: Reserved (0)

    Args:
        jd_max: Julian Day of maximum eclipse
        eclipse_type: Eclipse type flags
        moon_semidiameter: Moon's angular semi-diameter in degrees
        umbra_radius: Umbral shadow radius in degrees
        penumbra_radius: Penumbral shadow radius in degrees
        moon_lat_at_max: Moon's ecliptic latitude at maximum in degrees

    Returns:
        Tuple of 8 floats with phase times (JD UT)
    """
    # Actually calculate from Moon's speed minus Sun's speed
    sun_pos, _ = calc_ut(jd_max, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd_max, MOON, FLG_SPEED)

    # Speed of Moon relative to shadow (in longitude)
    relative_speed = abs(moon_pos[3] - sun_pos[3])  # degrees/day
    if relative_speed < 0.1:
        relative_speed = 0.55  # fallback to typical value

    # Distance from shadow axis at maximum
    y = abs(moon_lat_at_max)

    # Calculate half-duration of each phase using geometry
    # The Moon moves through the shadow at an angle
    # Half-duration = sqrt(R² - y²) / speed where R is radius and y is impact parameter

    def calc_half_duration(radius: float) -> float:
        """Calculate half-duration for given shadow radius.

        Includes edge case handling for shallow eclipses where
        the Moon barely touches the shadow boundary.
        """
        r_total = radius + moon_semidiameter
        if y >= r_total:
            return 0.0
        # Calculate chord half-length using safe sqrt
        chord_sq = r_total * r_total - y * y
        # Edge case: very shallow eclipse where chord is nearly zero
        if chord_sq < 1e-20:
            return 0.0
        half_chord = _safe_sqrt(chord_sq)
        return half_chord / relative_speed

    def calc_total_half_duration(radius: float) -> float:
        """Calculate half-duration of total phase (Moon fully inside).

        Includes edge case handling for nearly-total eclipses where
        the Moon barely fits inside the umbra.
        """
        r_inner = radius - moon_semidiameter
        if r_inner <= 0 or y >= r_inner:
            return 0.0
        chord_sq = r_inner * r_inner - y * y
        # Edge case: nearly-total eclipse where chord is nearly zero
        if chord_sq < 1e-20:
            return 0.0
        half_chord = _safe_sqrt(chord_sq)
        return half_chord / relative_speed

    # Penumbral phase times
    penumbral_half_dur = calc_half_duration(penumbra_radius)
    t_pen_begin = jd_max - penumbral_half_dur
    t_pen_end = jd_max + penumbral_half_dur

    # Partial (umbral) phase times
    partial_half_dur = calc_half_duration(umbra_radius)
    if partial_half_dur > 0:
        t_partial_begin = jd_max - partial_half_dur
        t_partial_end = jd_max + partial_half_dur
    else:
        t_partial_begin = 0.0
        t_partial_end = 0.0

    # Total phase times
    if eclipse_type & ECL_TOTAL:
        total_half_dur = calc_total_half_duration(umbra_radius)
        if total_half_dur > 0:
            t_total_begin = jd_max - total_half_dur
            t_total_end = jd_max + total_half_dur
        else:
            t_total_begin = 0.0
            t_total_end = 0.0
    else:
        t_total_begin = 0.0
        t_total_end = 0.0

    return (
        jd_max,  # [0] Maximum eclipse
        0.0,  # [1] Reserved
        t_partial_begin,  # [2] Partial begins (Moon enters umbra)
        t_partial_end,  # [3] Partial ends (Moon leaves umbra)
        t_total_begin,  # [4] Total begins (or 0 if not total)
        t_total_end,  # [5] Total ends (or 0 if not total)
        t_pen_begin,  # [6] Penumbral begins
        t_pen_end,  # [7] Penumbral ends
        0.0,  # [8] Reserved
        0.0,  # [9] Reserved
    )


def _lun_eclipse_when_pythonic(
    jd_start: float,
    flags: int = FLG_SWIEPH,
    eclipse_type: int = 0,
    backwards: bool = False,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Find the next (or previous) lunar eclipse from a given date.

    Searches forward (or backward when ``backwards`` is True) in time from
    jd_start. Can filter by eclipse type (total, partial, penumbral).

    Args:
        jd_start: Julian Day (UT) to start search from
        flags: Calculation flags (FLG_SWIEPH, etc.)
        eclipse_type: Filter for specific eclipse type(s), bitmask of:
            - ECL_TOTAL (4): Total lunar eclipse
            - ECL_PARTIAL (16): Partial lunar eclipse
            - ECL_PENUMBRAL (64): Penumbral lunar eclipse
            - 0: Any eclipse type (default)

    Returns:
        Tuple containing (matching reference API format):
            - retflag: Eclipse type flags bitmask (ECL_* constants)
            - tret: Tuple of 10 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse
                [1]: Reserved (0)
                [2]: Time of partial eclipse beginning (Moon enters umbra)
                [3]: Time of partial eclipse ending (Moon leaves umbra)
                [4]: Time of total eclipse beginning (or 0 if not total)
                [5]: Time of total eclipse ending (or 0 if not total)
                [6]: Time of penumbral eclipse beginning
                [7]: Time of penumbral eclipse ending
                [8]: Reserved (0)
                [9]: Reserved (0)

    Raises:
        RuntimeError: If no eclipse found within search limit

    Algorithm:
        1. Find next Full Moon after jd_start
        2. Check if Moon is close enough to node for eclipse
        3. If not eclipse, advance to next Full Moon
        4. Calculate eclipse type and magnitude
        5. If eclipse_type filter set, check if matches
        6. Calculate phase times

    Precision:
        Eclipse times accurate to ~1 minute for most eclipses.

    Example:
        >>> # Find next total lunar eclipse after Jan 1, 2024
        >>> from libephemeris import julday, ECL_TOTAL
        >>> jd = julday(2024, 1, 1, 0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd, eclipse_type=ECL_TOTAL)
        >>> print(f"Total lunar eclipse at JD {times[0]:.5f}")

    References:
        - Reference API: lun_eclipse_when()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 20  # Maximum search range
    MAX_FULL_MOONS = int(MAX_SEARCH_YEARS * 12.4)  # ~12.4 lunations per year

    # Mask out non-lunar eclipse type bits (CENTRAL=1, NONCENTRAL=2,
    # ANNULAR=8, ANNULAR_TOTAL=32) — pyswisseph ignores these for lunar
    # eclipses and treats them as "any type".
    eclipse_type = eclipse_type & ECL_ALLTYPES_LUNAR

    # If eclipse_type is 0, accept any type
    if eclipse_type == 0:
        eclipse_type = ECL_ALLTYPES_LUNAR

    jd = jd_start

    for _ in range(MAX_FULL_MOONS):
        # Find next (or previous) Full Moon
        if backwards:
            jd_full_moon = _find_previous_full_moon(jd)
        else:
            jd_full_moon = _find_next_full_moon(jd)

        # Get Moon position at Full Moon
        moon_pos, _ = calc_ut(jd_full_moon, MOON, flags | FLG_SPEED)
        moon_lon = moon_pos[0]
        moon_lat = moon_pos[1]

        # Check if close enough to ecliptic for eclipse
        # Lunar eclipse possible if Moon is near a node
        node_dist = _get_moon_node_distance(jd_full_moon, moon_lon)

        if node_dist < ECLIPSE_LIMIT_LUNAR:
            # Refine the eclipse maximum time from Full Moon approximation
            # Eclipse maximum occurs when Moon is closest to shadow axis,
            # not exactly at Full Moon (opposition in longitude)
            jd_max = _refine_lunar_eclipse_maximum(jd_full_moon)

            # Direction invariant: a forward search must return an eclipse
            # after jd_start, a backward search one before it (refinement
            # can drag the maximum across the start time).
            if (not backwards and jd_max <= jd_start) or (
                backwards and jd_max >= jd_start
            ):
                jd = jd_full_moon + (-25 if backwards else 25)
                continue

            # Possible eclipse - check magnitude at refined maximum time
            (
                ecl_type,
                umbral_mag,
                penumbral_mag,
                gamma,
                penumbra_radius,
                umbra_radius,
            ) = _calculate_lunar_eclipse_type_and_magnitude(jd_max)

            if ecl_type != 0:
                # Eclipse found (including grazing — pyswisseph reports
                # borderline eclipses as partial, so we must not skip them)
                type_matches = (
                    (eclipse_type & ECL_TOTAL and ecl_type & ECL_TOTAL)
                    or (eclipse_type & ECL_PARTIAL and ecl_type & ECL_PARTIAL)
                    or (eclipse_type & ECL_PENUMBRAL and ecl_type & ECL_PENUMBRAL)
                )

                if type_matches:
                    # Get moon position at refined maximum for phase calculations
                    moon_pos_max, _ = calc_ut(jd_max, MOON, flags | FLG_SPEED)
                    moon_dist = moon_pos_max[2]
                    moon_lat_max = moon_pos_max[1]
                    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

                    # Calculate phase times
                    times = _calculate_lunar_eclipse_phases(
                        jd_max,
                        ecl_type,
                        moon_semidiameter,
                        umbra_radius,
                        penumbra_radius,
                        moon_lat_max,
                    )
                    # Strip internal GRAZING flag — pyswisseph
                    # doesn't include it in the returned type
                    return ecl_type & ~ECL_GRAZING, times

        # Advance (or retreat) to the next lunation: ~25 days keeps us
        # safely within one synodic month of the next Full Moon.
        jd = jd_full_moon + (-25 if backwards else 25)

    raise RuntimeError(
        f"No matching lunar eclipse found within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


def lun_eclipse_when(
    tjdut: float,
    flags: int = FLG_SWIEPH,
    ecltype: int = 0,
    backwards: bool = False,
) -> Tuple[int, Tuple[float, ...]]:
    """Find the next (or previous) lunar eclipse globally (pyswisseph-compatible).

    Wrapper around _lun_eclipse_when_pythonic() matching pyswisseph signature.

    Args:
        tjdut: Julian Day (UT) to start search from.
        flags: Calculation flags (default FLG_SWIEPH).
        ecltype: Eclipse type filter bitmask (0 = any).
        backwards: If True, search backward in time.

    Returns:
        Tuple of (retflag, tret) matching pyswisseph.
    """
    return _lun_eclipse_when_pythonic(
        tjdut, flags=flags, eclipse_type=ecltype, backwards=backwards
    )


def _lun_eclipse_when_loc_pythonic(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    backwards: bool = False,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next lunar eclipse visible from a specific location.

    Searches forward in time from jd_start to find the next lunar eclipse
    where the Moon is above the horizon during at least part of the eclipse.
    The Moon must be above the horizon for the eclipse to be visible.

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 10 floats with eclipse phase times (JD UT):
                [0]: Time of maximum eclipse (local visibility)
                [1]: Time of partial eclipse beginning (Moon enters umbra)
                [2]: Time of total eclipse beginning (or 0 if not total)
                [3]: Time of total eclipse ending (or 0 if not total)
                [4]: Time of partial eclipse ending (Moon leaves umbra)
                [5]: Time of penumbral eclipse beginning
                [6]: Time of penumbral eclipse ending
                [7]: Time of moonrise (if Moon rises during eclipse, else 0)
                [8]: Time of moonset (if Moon sets during eclipse, else 0)
                [9]: Reserved (0)
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Umbral magnitude
                [1]: Penumbral magnitude
                [2]: Reserved (0)
                [3]: Azimuth of Moon at maximum eclipse (degrees)
                [4]: Altitude of Moon at maximum eclipse (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of umbral shadow (degrees)
                [7]: Apparent diameter of penumbral shadow (degrees)
                [8]: Saros series number (0, not implemented)
                [9-10]: Reserved (0)
            - retflag: Eclipse type flags bitmask (ECL_* constants)
                Includes visibility flags (ECL_VISIBLE, ECL_*_VISIBLE)

    Raises:
        RuntimeError: If no eclipse visible from location within search limit

    Algorithm:
        1. Use _lun_eclipse_when_pythonic to find next global lunar eclipse
        2. Calculate Moon's altitude at observer's location during eclipse
        3. If Moon is below horizon during entire eclipse, continue to next
        4. Return local phase times, visibility flags, and attributes

    Precision:
        Eclipse times accurate to ~1 minute. Visibility depends on
        accurate horizon calculations.

    Example:
        >>> # Find next lunar eclipse visible from Rome
        >>> from libephemeris import julday, _lun_eclipse_when_loc_pythonic
        >>> jd = julday(2024, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> ecl_type, times, attr = _lun_eclipse_when_loc_pythonic(jd, rome_lat, rome_lon)
        >>> print(f"Eclipse maximum at JD {times[0]:.5f}")
        >>> print(f"Moon altitude: {attr[4]:.1f}°")

    References:
        - Reference API: lun_eclipse_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    MAX_SEARCH_YEARS = 50  # Maximum search range
    MAX_ECLIPSES = int(
        MAX_SEARCH_YEARS * 2.4
    )  # ~2.4 lunar eclipses per year on average

    # reader is provided by the caller (None forces Skyfield path)

    # Geometric center altitude threshold for moonrise/moonset.
    # Atmospheric refraction near the horizon lifts the Moon's apparent
    # position by roughly 34'.  The reference API considers the Moon
    # "visible" until refraction can no longer bring its centre above
    # the geometric horizon, which occurs at a geometric centre altitude
    # of approximately -0.36°.  This empirical value matches the
    # reference API's moonrise/moonset times to within a few seconds
    # across a wide range of eclipses and observer locations.
    _MOON_HORIZON_ALT = -0.36

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        _gp = (lon, lat, altitude)

        def _get_moon_altaz(jd: float) -> Tuple[float, float]:
            """Get Moon geometric altitude and azimuth at given JD."""
            jd_tt = jd + deltat(jd)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
            moon_az, moon_alt_true, moon_alt_app = azalt(jd, ECL2HOR, _gp, 0, 0, moon_pos[:3])
            return moon_alt_true, moon_az

        def _get_moon_apparent_alt(jd: float) -> float:
            """Get Moon apparent altitude with atmospheric refraction."""
            jd_tt = jd + deltat(jd)
            moon_pos = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
            moon_az, moon_alt_true, moon_alt_app = azalt(jd, ECL2HOR, _gp, 1013.25, 10.0, moon_pos[:3])
            return moon_alt_app
    else:
        from .state import get_planets, get_timescale
        from skyfield.api import wgs84

        # Get ephemeris
        eph = get_planets()
        ts = get_timescale()

        # Create observer
        earth = eph["earth"]
        moon_obj = eph["moon"]
        observer = wgs84.latlon(lat, lon, altitude)

        def _get_moon_altaz(jd: float) -> Tuple[float, float]:
            """Get Moon geometric altitude and azimuth at given JD.

            Returns:
                Tuple of (geometric_altitude_degrees, azimuth_degrees).
                Azimuth uses the reference API convention (S=0, W=90).
            """
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            moon_app = observer_at.at(t).observe(moon_obj).apparent()
            alt, az, _ = moon_app.altaz()
            # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_moon_apparent_alt(jd: float) -> float:
            """Get Moon apparent altitude with atmospheric refraction.

            Uses standard atmosphere (T=10°C, P=1013.25 mbar) for
            atmospheric refraction correction, matching the conditions
            used by the reference API.

            Returns:
                Apparent altitude in degrees (refraction-corrected).
            """
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            moon_app = observer_at.at(t).observe(moon_obj).apparent()
            alt, _, _ = moon_app.altaz(temperature_C=10.0, pressure_mbar=1013.25)
            return alt.degrees

    def _is_moon_visible(jd: float, min_alt: float = -1.0) -> bool:
        """Check if Moon is above horizon (accounting for refraction).

        Uses the geometric centre altitude threshold that corresponds
        to the refraction-corrected horizon.
        """
        alt, _ = _get_moon_altaz(jd)
        return alt > max(min_alt, _MOON_HORIZON_ALT)

    jd = jd_start

    for _ in range(MAX_ECLIPSES):
        # Find next global lunar eclipse
        try:
            global_type, global_times = _lun_eclipse_when_pythonic(
                jd, flags, backwards=backwards
            )
        except RuntimeError:
            raise RuntimeError(
                f"No lunar eclipse visible from lat={lat}, lon={lon} "
                f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
            )

        jd_max = global_times[0]
        jd_partial_begin = global_times[2]
        jd_partial_end = global_times[3]
        jd_total_begin = global_times[4]
        jd_total_end = global_times[5]
        jd_pen_begin = global_times[6]
        jd_pen_end = global_times[7]

        # Check if Moon is visible during any part of the eclipse
        # Sample at key phases: penumbral begin/end, partial begin/end, maximum
        check_times = [jd_max]
        if jd_pen_begin > 0:
            check_times.append(jd_pen_begin)
        if jd_pen_end > 0:
            check_times.append(jd_pen_end)
        if jd_partial_begin > 0:
            check_times.append(jd_partial_begin)
        if jd_partial_end > 0:
            check_times.append(jd_partial_end)
        if jd_total_begin > 0:
            check_times.append(jd_total_begin)
        if jd_total_end > 0:
            check_times.append(jd_total_end)

        # Also sample at several intermediate points
        eclipse_start = jd_pen_begin if jd_pen_begin > 0 else jd_max - 2 / 24
        eclipse_end = jd_pen_end if jd_pen_end > 0 else jd_max + 2 / 24
        for i in range(5):
            t = eclipse_start + (eclipse_end - eclipse_start) * i / 4
            if t not in check_times:
                check_times.append(t)

        moon_visible = any(_is_moon_visible(t) for t in check_times)

        if not moon_visible:
            # Eclipse not visible from this location, try next
            jd = jd_max + (-25 if backwards else 25)
            continue

        # Eclipse visible! Calculate local circumstances

        # Find moonrise/moonset during eclipse first, as we need these
        # both for the local maximum determination and for the output.
        eclipse_start = jd_pen_begin if jd_pen_begin > 0 else jd_max - 3 / 24
        eclipse_end = jd_pen_end if jd_pen_end > 0 else jd_max + 3 / 24

        moonrise_time = 0.0
        moonset_time = 0.0

        if eclipse_start > 0 and eclipse_end > 0:
            alt_start, _ = _get_moon_altaz(eclipse_start)
            alt_end, _ = _get_moon_altaz(eclipse_end)

            # Use refraction-corrected horizon threshold for rise/set.
            horizon = _MOON_HORIZON_ALT

            if alt_start < horizon < alt_end:
                # Moon rises during eclipse — binary search for crossing
                t_low, t_high = eclipse_start, eclipse_end
                for _ in range(50):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _ = _get_moon_altaz(t_mid)
                    if alt_mid < horizon:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                    if t_high - t_low < 1e-9:
                        break
                moonrise_time = (t_low + t_high) / 2

            elif alt_start > horizon > alt_end:
                # Moon sets during eclipse — binary search for crossing
                t_low, t_high = eclipse_start, eclipse_end
                for _ in range(50):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _ = _get_moon_altaz(t_mid)
                    if alt_mid > horizon:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                    if t_high - t_low < 1e-9:
                        break
                moonset_time = (t_low + t_high) / 2

        # Determine the local maximum time. If Moon is below the
        # refraction-corrected horizon at the global maximum, use the
        # moonrise or moonset time as the local max (matching reference
        # API behavior — reports the horizon crossing as the "maximum"
        # for partially-visible eclipses).
        moon_alt_at_max, _ = _get_moon_altaz(jd_max)

        jd_local_max = jd_max
        if moon_alt_at_max < _MOON_HORIZON_ALT:
            # Moon is below horizon at global max. Use moonrise/moonset
            # as local max, whichever is closest to the global max.
            if moonset_time > 0 and moonrise_time > 0:
                # Both rise and set during eclipse — use whichever is
                # closer to the global max
                if abs(moonset_time - jd_max) < abs(moonrise_time - jd_max):
                    jd_local_max = moonset_time
                else:
                    jd_local_max = moonrise_time
            elif moonset_time > 0:
                jd_local_max = moonset_time
            elif moonrise_time > 0:
                jd_local_max = moonrise_time
            # else: fall back to global max (shouldn't happen since
            # moon_visible check passed above)

        # Get Moon position at local maximum
        moon_pos, _ = calc_ut(jd_local_max, MOON, FLG_SPEED)
        moon_dist = moon_pos[2]

        # Get eclipse type and magnitude at local maximum
        (
            ecl_type_flags,
            umbral_mag,
            penumbral_mag,
            gamma,
            penumbra_radius,
            umbra_radius,
        ) = _calculate_lunar_eclipse_type_and_magnitude(jd_local_max)

        # Get Moon's alt/az at local maximum
        moon_alt, moon_az = _get_moon_altaz(jd_local_max)
        moon_app_alt = _get_moon_apparent_alt(jd_local_max)

        # Calculate apparent diameters
        # Moon semi-diameter: 932.56 arcsec at mean distance 0.002569 AU
        moon_diameter = 2 * (932.56 / 3600.0) * (0.002569 / moon_dist)
        umbra_diameter = 2 * umbra_radius
        penumbra_diameter = 2 * penumbra_radius

        # Determine which phases are visible
        ecl_type = ecl_type_flags

        # Check visibility at each phase
        if jd_pen_begin > 0 and _is_moon_visible(jd_pen_begin):
            ecl_type |= ECL_1ST_VISIBLE
        if jd_partial_begin > 0 and _is_moon_visible(jd_partial_begin):
            ecl_type |= ECL_2ND_VISIBLE
        if jd_partial_end > 0 and _is_moon_visible(jd_partial_end):
            ecl_type |= ECL_3RD_VISIBLE
        if jd_pen_end > 0 and _is_moon_visible(jd_pen_end):
            ecl_type |= ECL_4TH_VISIBLE
        if _is_moon_visible(jd_max):
            ecl_type |= ECL_MAX_VISIBLE
        ecl_type |= ECL_VISIBLE

        # Filter contact times by Moon visibility — zero out phases
        # where the Moon is below the horizon (not visible to observer)
        local_partial_begin = jd_partial_begin
        local_partial_end = jd_partial_end
        local_total_begin = jd_total_begin
        local_total_end = jd_total_end
        local_pen_begin = jd_pen_begin
        local_pen_end = jd_pen_end

        if local_pen_begin > 0 and not _is_moon_visible(local_pen_begin):
            local_pen_begin = 0.0
        if local_pen_end > 0 and not _is_moon_visible(local_pen_end):
            local_pen_end = 0.0
        if local_partial_begin > 0 and not _is_moon_visible(local_partial_begin):
            local_partial_begin = 0.0
        if local_partial_end > 0 and not _is_moon_visible(local_partial_end):
            local_partial_end = 0.0
        if local_total_begin > 0 and not _is_moon_visible(local_total_begin):
            local_total_begin = 0.0
        if local_total_end > 0 and not _is_moon_visible(local_total_end):
            local_total_end = 0.0

        # Prepare times tuple (10 elements matching pyswisseph layout)
        times = (
            jd_local_max,  # [0] Maximum (local visibility)
            0.0,  # [1] Reserved
            local_partial_begin,  # [2] Partial begins (enters umbra)
            local_partial_end,  # [3] Partial ends (leaves umbra)
            local_total_begin,  # [4] Total begins
            local_total_end,  # [5] Total ends
            local_pen_begin,  # [6] Penumbral begins
            local_pen_end,  # [7] Penumbral ends
            moonrise_time,  # [8] Moonrise during eclipse
            moonset_time,  # [9] Moonset during eclipse
        )

        # Prepare attributes tuple (20 elements matching pyswisseph layout)
        attr = (
            umbral_mag,  # [0] Umbral magnitude
            penumbral_mag,  # [1] Penumbral magnitude
            0.0,  # [2] Reserved
            0.0,  # [3] Reserved
            moon_az,  # [4] Azimuth of Moon at maximum
            moon_alt,  # [5] True altitude of Moon at maximum
            moon_app_alt,  # [6] Apparent altitude (with refraction)
            0.0,  # [7] Distance from opposition (degrees)
            umbral_mag,  # [8] Umbral magnitude (equals [0])
            *_get_saros_info(global_times[0], "lunar"),  # [9] Saros, [10] member
            0.0,  # [11] Reserved
            0.0,  # [12] Reserved
            0.0,  # [13] Reserved
            0.0,  # [14] Reserved
            0.0,  # [15] Reserved
            0.0,  # [16] Reserved
            0.0,  # [17] Reserved
            0.0,  # [18] Reserved
            0.0,  # [19] Reserved
        )

        return ecl_type, times, attr

    raise RuntimeError(
        f"No lunar eclipse visible from lat={lat}, lon={lon} "
        f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


def lun_eclipse_when_loc(
    tjdut: float,
    geopos: "Sequence[float]",
    flags: int = FLG_SWIEPH,
    backwards: bool = False,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find the next lunar eclipse visible from a geographic position (pyswisseph-compatible).

    Wrapper around _lun_eclipse_when_loc_pythonic() matching pyswisseph signature.

    Args:
        tjdut: Julian Day (UT) to start search from.
        geopos: Sequence of [longitude, latitude, altitude].
        flags: Calculation flags (default FLG_SWIEPH).
        backwards: If True, search backward in time.

    Returns:
        Tuple of (retflag, tret, attr) matching pyswisseph.
    """
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    return _call_with_leb_skyfield_fallback(
        _lun_eclipse_when_loc_pythonic, tjdut, lat, lon, altitude, flags, backwards,
    )


def _lun_eclipse_how_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate the circumstances of a lunar eclipse at a specific location and time.

    This function determines how a lunar eclipse appears from a given geographic
    position at a specific Julian Day. Unlike _lun_eclipse_when_loc_pythonic which finds the
    next eclipse, this function calculates the eclipse magnitude, Moon position,
    and other circumstances for a known eclipse time.

    Args:
        jd: Julian Day (UT) of the moment to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - retflag: Eclipse type flags bitmask (ECL_* constants)
                Returns 0 if no eclipse is occurring at this time
                Includes ECL_VISIBLE if Moon is above horizon
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Umbral eclipse magnitude (fraction of Moon's diameter in umbra)
                [1]: Penumbral eclipse magnitude
                [2]: Reserved (0)
                [3]: Azimuth of Moon at the given time (degrees)
                [4]: Altitude of Moon at the given time (degrees)
                [5]: Apparent diameter of Moon (degrees)
                [6]: Apparent diameter of umbral shadow at Moon's distance (degrees)
                [7]: Apparent diameter of penumbral shadow at Moon's distance (degrees)
                [8]: Saros series number (0, not implemented)
                [9]: Reserved (0)
                [10]: Reserved (0)

    Note:
        This function is intended for use when you already know an eclipse is
        occurring (e.g., from _lun_eclipse_when_pythonic or _lun_eclipse_when_loc_pythonic).
        For a random time when no eclipse is occurring, magnitude will be 0
        and retflag will be 0.

    Algorithm:
        1. Calculate Moon's apparent position from observer location
        2. Calculate Earth's shadow cone geometry at Moon's distance
        3. Compute umbral and penumbral magnitudes
        4. Determine eclipse type based on Moon's penetration into shadows
        5. Return attributes including Moon's altitude and azimuth

    Precision:
        Eclipse magnitude accurate to ~0.001 for central eclipses.
        Topocentric parallax included in calculations.

    Example:
        >>> # Calculate eclipse circumstances at Rome during a lunar eclipse
        >>> from libephemeris import julday, _lun_eclipse_how_pythonic
        >>> jd = julday(2022, 5, 16, 4.2)  # During May 2022 eclipse
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> ecl_type, attr = _lun_eclipse_how_pythonic(jd, rome_lat, rome_lon)
        >>> print(f"Umbral magnitude: {attr[0]:.3f}")
        >>> print(f"Moon altitude: {attr[4]:.1f}°")

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        _gp = (lon, lat, altitude)

        try:
            jd_tt = jd + deltat(jd)
            moon_topo = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0, tuple([0.0] * 20)

        moon_az_val, moon_alt_true, moon_alt_app = azalt(jd, ECL2HOR, _gp, 0, 0, moon_topo[:3])
        moon_altitude = moon_alt_true
        moon_azimuth = moon_az_val
        moon_dist_au = moon_topo[2]
    else:
        from skyfield.api import wgs84

        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Moon object
        moon = eph["moon"]
        earth = eph["earth"]

        # Get Skyfield time
        t = ts.ut1_jd(jd)

        # Create observer position
        observer_at = earth + observer

        # Get Moon apparent position from observer
        try:
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # If calculation fails, return zeros (20 elements)
            return 0, tuple([0.0] * 20)

        # Get Moon altitude and azimuth
        moon_alt, moon_az, _ = moon_app.altaz()
        moon_altitude = moon_alt.degrees
        # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
        moon_azimuth = (moon_az.degrees + 180.0) % 360.0

        # Get Moon's distance for angular size calculation
        moon_dist_au = moon_app.distance().au

    # Calculate eclipse geometry using the same calculations as _lun_eclipse_when_pythonic
    (
        ecl_type_flags,
        umbral_mag,
        penumbral_mag,
        gamma,
        penumbra_radius,
        umbra_radius,
    ) = _calculate_lunar_eclipse_type_and_magnitude(jd)

    # Calculate apparent diameters
    # Moon semi-diameter: 932.56 arcsec at mean distance 0.002569 AU
    moon_diameter = 2 * (932.56 / 3600.0) * (0.002569 / moon_dist_au)
    umbra_diameter = 2 * umbra_radius
    penumbra_diameter = 2 * penumbra_radius

    # Determine eclipse type flags
    eclipse_type = 0

    if penumbral_mag <= 0 and umbral_mag <= 0:
        # No eclipse - Moon too far from Earth's shadow
        return 0, (
            0.0,  # [0] umbral magnitude
            0.0,  # [1] penumbral magnitude
            0.0,  # [2] reserved
            0.0,  # [3] reserved
            moon_azimuth,  # [4] Moon azimuth
            moon_altitude,  # [5] true altitude
            moon_altitude,  # [6] apparent altitude
            0.0,  # [7] distance from opposition
            0.0,  # [8] eclipse magnitude
            0.0,  # [9] Saros series number
            0.0,  # [10] Saros series member
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,  # [11-19] Reserved
        )

    # There is an eclipse - set type flags
    eclipse_type = ecl_type_flags

    # Check if Moon is above horizon
    if moon_altitude > -1.0:  # Allow for refraction near horizon
        eclipse_type |= ECL_VISIBLE
        eclipse_type |= ECL_MAX_VISIBLE

    # Prepare attributes tuple (20 elements matching pyswisseph layout)
    attr = (
        max(0.0, umbral_mag),  # [0] Umbral magnitude
        max(0.0, penumbral_mag),  # [1] Penumbral magnitude
        0.0,  # [2] Reserved
        0.0,  # [3] Reserved
        moon_azimuth,  # [4] Azimuth of Moon
        moon_altitude,  # [5] True altitude of Moon
        moon_altitude,  # [6] Apparent altitude (approx)
        0.0,  # [7] Distance from opposition (degrees)
        max(0.0, umbral_mag),  # [8] Eclipse magnitude (equals [0])
        *_get_saros_info(jd, "lunar"),  # [9] Saros, [10] member
        0.0,  # [11] Reserved
        0.0,  # [12] Reserved
        0.0,  # [13] Reserved
        0.0,  # [14] Reserved
        0.0,  # [15] Reserved
        0.0,  # [16] Reserved
        0.0,  # [17] Reserved
        0.0,  # [18] Reserved
        0.0,  # [19] Reserved
    )

    return eclipse_type, attr


def lun_eclipse_how(
    tjdut: float,
    geopos: Sequence[float],
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate detailed circumstances of a lunar eclipse from a location.

    Wrapper around :func:`_lun_eclipse_how_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _lun_eclipse_how_impl, tjdut, geopos, flags,
    )


def _lun_eclipse_how_impl(
    tjdut: float,
    geopos: Sequence[float],
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate detailed circumstances of a lunar eclipse from a specific location.

    This function matches the reference API signature exactly. It calculates
    the eclipse circumstances at a specific Julian Day and geographic location,
    including magnitudes, Moon position, and visibility flags.

    Unlike solar eclipses, lunar eclipses look the same from everywhere on
    the night side of Earth, but we need to check if the Moon is above
    the horizon from the observer's location.

    Args:
        tjdut: Julian Day (UT) during a lunar eclipse
        flags: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive) - NOTE: longitude first!
            - latitude in degrees (North positive)
            - altitude in meters above sea level

    Returns:
        Tuple containing:
            - retflag: Eclipse type flags bitmask combined with visibility flags:
                Returns eclipse type (ECL_TOTAL, ECL_PARTIAL, ECL_PENUMBRAL)
                combined with ECL_VISIBLE (128) if Moon is above horizon.
                Returns 0 if no eclipse is occurring at this time.
            - attr: Tuple of 11 floats with eclipse attributes:
                [0]: Umbral eclipse magnitude (fraction of Moon diameter
                     covered by umbra; >1 means Moon fully in umbra)
                [1]: Penumbral eclipse magnitude
                [2]: Reserved (0)
                [3]: Reserved (0)
                [4]: Azimuth of Moon at tjdut (degrees)
                [5]: True altitude of Moon at tjdut (degrees)
                [6]: Apparent altitude of Moon with atmospheric refraction (degrees)
                [7]: Distance of Moon center from shadow axis (in Moon radii)
                [8]: Eclipse type at this moment (ECL_TOTAL, ECL_PARTIAL,
                     ECL_PENUMBRAL, or 0)
                [9]: Apparent diameter of Moon (degrees)
                [10]: Apparent diameter of umbral shadow (degrees)

    Note:
        This function is intended for use when you already know an eclipse is
        occurring (e.g., from lun_eclipse_when). For a random time when
        no eclipse is occurring, magnitude will be 0 and retflag will be 0.

    Algorithm:
        1. Get Moon's apparent position from observer location (topocentric)
        2. Calculate Earth's shadow cone geometry at Moon's distance
        3. Compute umbral and penumbral magnitudes
        4. Calculate Moon's distance from shadow axis in Moon radii
        5. Determine eclipse type based on penetration into shadow
        6. Apply atmospheric refraction for apparent altitude
        7. Set visibility flags if Moon is above horizon

    Precision:
        Eclipse magnitude accurate to ~0.01.
        Moon altitude accurate to within ~1 degree.

    Example:
        >>> # Calculate eclipse circumstances at Los Angeles during Nov 8, 2022 eclipse
        >>> from libephemeris import lun_eclipse_how, FLG_SWIEPH
        >>> jd = 2459892.4  # Maximum of Nov 8, 2022 total lunar eclipse
        >>> la_geopos = [-118.24, 34.05, 0]  # lon, lat, alt
        >>> ecl_type, attr = lun_eclipse_how(jd, la_geopos, FLG_SWIEPH)
        >>> print(f"Umbral magnitude: {attr[0]:.3f}")
        >>> print(f"Moon altitude: {attr[5]:.1f}°")

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # Validate and extract geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        _gp = (lon, lat, altitude)

        try:
            jd_tt = tjdut + deltat(tjdut)
            moon_topo = _topo_ecliptic(reader, jd_tt, tjdut, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0, tuple([0.0] * 20)

        moon_az_val, moon_alt_true, moon_alt_app = azalt(tjdut, ECL2HOR, _gp, 0, 0, moon_topo[:3])
        moon_altitude = moon_alt_true
        moon_azimuth = moon_az_val
        moon_dist_au = moon_topo[2]

        # Calculate apparent altitude with atmospheric refraction
        if moon_altitude > -1.0:
            if moon_altitude > 0:
                refraction = 1.0 / math.tan(
                    math.radians(moon_altitude + 7.31 / (moon_altitude + 4.4))
                )
                apparent_altitude = moon_altitude + refraction / 60.0
            else:
                apparent_altitude = moon_altitude + 0.58
        else:
            apparent_altitude = moon_altitude
    else:
        from skyfield.api import wgs84

        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Moon object
        moon = eph["moon"]
        earth = eph["earth"]

        # Get Skyfield time
        t = ts.ut1_jd(tjdut)

        # Create observer position
        observer_at = earth + observer

        # Get Moon apparent position from observer (topocentric)
        try:
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # If calculation fails, return zeros (20 elements)
            return 0, tuple([0.0] * 20)

        # Get Moon altitude and azimuth
        moon_alt, moon_az, _ = moon_app.altaz()
        moon_altitude = moon_alt.degrees
        # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
        moon_azimuth = (moon_az.degrees + 180.0) % 360.0

        # Get Moon's distance for angular size calculation
        moon_dist_au = moon_app.distance().au

        # Calculate apparent altitude with atmospheric refraction
        # Bennett's formula for atmospheric refraction
        if moon_altitude > -1.0:
            if moon_altitude > 0:
                # Standard atmospheric refraction formula
                refraction = 1.0 / math.tan(
                    math.radians(moon_altitude + 7.31 / (moon_altitude + 4.4))
                )
                apparent_altitude = moon_altitude + refraction / 60.0
            else:
                # Near horizon, use approximation
                # Refraction at horizon is about 34 arcminutes
                apparent_altitude = moon_altitude + 0.58
        else:
            apparent_altitude = moon_altitude

    # Calculate eclipse geometry using the same calculations as _lun_eclipse_when_pythonic
    (
        ecl_type_flags,
        umbral_mag,
        penumbral_mag,
        gamma,
        penumbra_radius,
        umbra_radius,
    ) = _calculate_lunar_eclipse_type_and_magnitude(tjdut)

    # Calculate apparent diameters
    # Moon semi-diameter: 932.56 arcsec at mean distance 0.002569 AU
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist_au)
    moon_diameter = 2 * moon_semidiameter
    umbra_diameter = 2 * umbra_radius
    penumbra_diameter = 2 * penumbra_radius

    # Calculate distance from shadow axis in Moon radii
    # Get Moon's ecliptic latitude at this time
    moon_pos, _ = calc_ut(tjdut, MOON, flags | FLG_SPEED)
    moon_lat = moon_pos[1]

    # The center distance is how far the Moon center is from the shadow axis
    # This is approximately the Moon's ecliptic latitude
    # expressed in Moon radii: center_distance = moon_lat / moon_semidiameter
    center_distance_radii = abs(moon_lat) / moon_semidiameter

    # Determine eclipse type flags and current phase type
    eclipse_type = 0
    current_phase_type = 0

    if penumbral_mag <= 0 and umbral_mag <= 0:
        # No eclipse - Moon too far from Earth's shadow
        return 0, (
            0.0,  # [0] umbral magnitude
            0.0,  # [1] penumbral magnitude
            0.0,  # [2] reserved
            0.0,  # [3] reserved
            moon_azimuth,  # [4] Moon azimuth
            moon_altitude,  # [5] true altitude
            apparent_altitude,  # [6] apparent altitude
            center_distance_radii,  # [7] distance from opposition
            0.0,  # [8] eclipse magnitude
            0.0,  # [9] Saros series number
            0.0,  # [10] Saros series member
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,  # [11-19] Reserved
        )

    # There is an eclipse - set type flags
    eclipse_type = ecl_type_flags

    # Determine current phase type at this moment
    if umbral_mag > 1.0:
        current_phase_type = ECL_TOTAL
    elif umbral_mag > 0:
        current_phase_type = ECL_PARTIAL
    elif penumbral_mag > 0:
        current_phase_type = ECL_PENUMBRAL

    # Check if Moon is above horizon and set visibility flags
    # Moon must be above horizon for eclipse to be visible
    if moon_altitude > -1.0:  # Allow for refraction near horizon
        eclipse_type |= ECL_VISIBLE
        # Also indicate which phase is visible
        if current_phase_type:
            eclipse_type |= ECL_MAX_VISIBLE

    # Prepare attributes tuple (20 elements matching pyswisseph layout)
    attr = (
        max(0.0, umbral_mag),  # [0] Umbral magnitude
        max(0.0, penumbral_mag),  # [1] Penumbral magnitude
        0.0,  # [2] Reserved
        0.0,  # [3] Reserved
        moon_azimuth,  # [4] Azimuth of Moon
        moon_altitude,  # [5] True altitude of Moon
        apparent_altitude,  # [6] Apparent altitude with refraction
        center_distance_radii,  # [7] Distance from opposition (degrees)
        max(0.0, umbral_mag),  # [8] Eclipse magnitude (equals [0])
        *_get_saros_info(tjdut, "lunar"),  # [9] Saros, [10] member
        0.0,  # [11] Reserved
        0.0,  # [12] Reserved
        0.0,  # [13] Reserved
        0.0,  # [14] Reserved
        0.0,  # [15] Reserved
        0.0,  # [16] Reserved
        0.0,  # [17] Reserved
        0.0,  # [18] Reserved
        0.0,  # [19] Reserved
    )

    return eclipse_type, attr


def _kernel_jd_bounds(eph) -> Tuple[float, float]:
    """Return (earliest_start_jd, latest_end_jd) across all loaded SPK segments.

    Used to clamp occultation searches to the loaded DE kernel's coverage so the
    batch scan never requests JDs outside it. Split kernels (e.g. DE441, whose
    every body is broken into two segments at JD 2440432.5 / 1969-07-30) expose
    several end_jd values; the kernel's true upper bound is the *maximum* end_jd.
    Taking the minimum would clamp the search to the split point (~1969) and drop
    every later event.

    Args:
        eph: Loaded ephemeris object containing SPK segments.

    Returns:
        Tuple of (minimum segment start_jd, maximum segment end_jd).
    """
    starts = [seg.spk_segment.start_jd for seg in eph.segments]
    ends = [seg.spk_segment.end_jd for seg in eph.segments]
    return min(starts), max(ends)


def lun_occult_when_glob(
    tjdut: float,
    body: "int | str",
    flags: int = FLG_SWIEPH,
    ecltype: int = 0,
    backwards: bool = False,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Find the next lunar occultation of a planet or fixed star globally (UT).

    A lunar occultation occurs when the Moon passes in front of (occults)
    a planet or star as seen from Earth. This function searches forward
    (or backward) in time to find the next such event globally (somewhere on Earth).

    This function matches the reference lun_occult_when_glob() API exactly.

    Args:
        tjdut: Julian Day (UT) to start search from
        body: Planet identifier (int) or star name (str).
            Use MERCURY, VENUS, etc. for planets.
            Use a star name string (e.g., "Regulus") for star occultations.
        flags: Calculation flags (FLG_SWIEPH, FLG_J2000, etc.)
        ecltype: Bit flags for eclipse type filter (0 = any type):
            ECL_CENTRAL, ECL_NONCENTRAL, ECL_TOTAL,
            ECL_ANNULAR, ECL_PARTIAL, ECL_ANNULAR_TOTAL
        backwards: If True, search backward in time; if False, search forward.

    Returns:
        Tuple containing:
            - retflags: Returned bit flags (int):
                - 0 if no occultation found
                - ECL_TOTAL or ECL_ANNULAR or ECL_PARTIAL
                - ECL_CENTRAL, ECL_NONCENTRAL
            - tret: Tuple of 10 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation
                [1]: Time when occultation takes place at local apparent noon
                [2]: Time of occultation begin
                [3]: Time of occultation end
                [4]: Time of totality begin
                [5]: Time of totality end
                [6]: Time of center line begin
                [7]: Time of center line end
                [8]: Time when annular-total occultation becomes total
                [9]: Time when annular-total occultation becomes annular again

    Raises:
        ValueError: If neither planet ID nor star name is specified

    Algorithm:
        1. Calculate Moon's position and the target body's position
        2. Find the next conjunction in right ascension or longitude
        3. At conjunction, check angular separation
        4. If separation < Moon's angular radius + target's angular radius,
           an occultation occurs
        5. Refine timing using Newton-Raphson iteration
        6. Calculate contact times

    Precision:
        Occultation times accurate to ~1 minute for most events.

    Example:
        >>> # Find next occultation of Regulus by the Moon
        >>> from libephemeris import julday
        >>> jd = julday(2024, 1, 1, 0)
        >>> retflags, tret = lun_occult_when_glob(jd, "Regulus", FLG_SWIEPH, 0)
        >>> print(f"Occultation at JD {tret[0]:.5f}")

        >>> # Find next occultation of Venus by the Moon
        >>> retflags, tret = lun_occult_when_glob(jd, VENUS, FLG_SWIEPH, 0)
        >>> print(f"Venus occultation at JD {tret[0]:.5f}")

    References:
        - Reference API: lun_occult_when_glob()
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from .state import get_leb_reader, get_planets, get_timescale
    from .fixed_stars import fixstar_ut
    from .constants import (
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
    )
    from .planets import _PLANET_MAP, get_planet_target

    # Interpret body parameter: int = planet ID, str = star name
    if isinstance(body, str):
        planet = 0
        star_name = body
    else:
        planet = body
        star_name = ""

    if planet == 0 and not star_name:
        raise ValueError(
            "Either planet ID or star name must be specified for occultation search"
        )

    jd_start = tjdut

    MAX_SEARCH_YEARS = 20
    MAX_ITERATIONS = int(MAX_SEARCH_YEARS * 365)  # Check roughly daily

    MOON_MEAN_ANGULAR_RADIUS_DEG = 932.56 / 3600.0

    reader = get_leb_reader()
    if reader is not None:
        from .constants import FLG_EQUATORIAL

        def _get_moon_position(jd: float) -> Tuple[float, float, float, float]:
            moon_eq, _ = calc_ut(jd, MOON, FLG_EQUATORIAL | FLG_SPEED)
            dist_au = moon_eq[2]
            moon_angular_radius = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist_au)
            return moon_eq[0], moon_eq[1], dist_au, moon_angular_radius
    else:
        eph = get_planets()
        ts = get_timescale()
        earth = eph["earth"]
        moon_body = eph["moon"]

        def _get_moon_position(jd: float) -> Tuple[float, float, float, float]:
            t = ts.ut1_jd(jd)
            moon_app = earth.at(t).observe(moon_body).apparent()
            ra, dec, dist = moon_app.radec(epoch="date")
            moon_angular_radius = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist.au)
            return ra.hours * 15.0, dec.degrees, dist.au, moon_angular_radius

    if reader is not None:
        def _get_target_position(jd: float) -> Tuple[float, float, float]:
            if planet == 0:
                from .fixed_stars import fixstar_ut
                pos, _, _ = fixstar_ut(star_name, jd, FLG_EQUATORIAL | FLG_SPEED)
                return pos[0], pos[1], 0.0001
            else:
                if planet not in _PLANET_MAP:
                    raise ValueError(f"illegal planet number {planet}.")
                target_eq, _ = calc_ut(jd, planet, FLG_EQUATORIAL | FLG_SPEED)
                angular_radius = _calc_planet_angular_radius(planet, target_eq[2])
                return target_eq[0], target_eq[1], angular_radius
    else:
        def _get_target_position(jd: float) -> Tuple[float, float, float]:
            if planet == 0:
                from .fixed_stars import FIXED_STARS
                from .fixed_stars import _resolve_star_id
                from skyfield.api import Star as SkyfieldStar

                star_id, err, _ = _resolve_star_id(star_name)
                if err is not None:
                    raise ValueError(err)
                star = FIXED_STARS[star_id]
                t = ts.ut1_jd(jd)
                star_obj = SkyfieldStar(
                    ra_hours=star.ra_j2000 / 15.0,
                    dec_degrees=star.dec_j2000,
                    ra_mas_per_year=star.pm_ra * 1000.0,
                    dec_mas_per_year=star.pm_dec * 1000.0,
                )
                star_app = earth.at(t).observe(star_obj).apparent()
                ra, dec, _ = star_app.radec(epoch="date")
                return ra.hours * 15.0, dec.degrees, 0.0001
            else:
                if planet not in _PLANET_MAP:
                    raise ValueError(f"illegal planet number {planet}.")
                target_name = _PLANET_MAP[planet]
                target = get_planet_target(eph, target_name)
                t = ts.ut1_jd(jd)
                target_app = earth.at(t).observe(target).apparent()
                ra, dec, dist = target_app.radec(epoch="date")
                angular_radius = _calc_planet_angular_radius(planet, dist.au)
                return ra.hours * 15.0, dec.degrees, angular_radius

    def _angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
        """Calculate angular separation between two points (in degrees)."""
        # Convert to radians
        ra1_r = math.radians(ra1)
        dec1_r = math.radians(dec1)
        ra2_r = math.radians(ra2)
        dec2_r = math.radians(dec2)

        # Haversine formula for spherical distance
        cos_sep = math.sin(dec1_r) * math.sin(dec2_r) + math.cos(dec1_r) * math.cos(
            dec2_r
        ) * math.cos(ra1_r - ra2_r)
        cos_sep = max(-1.0, min(1.0, cos_sep))
        return math.degrees(math.acos(cos_sep))

    def _check_occultation(jd: float) -> Tuple[bool, float, float, float]:
        """
        Check if occultation occurs at given time.

        For global occultation detection, we use an expanded threshold that
        accounts for lunar parallax (~1 degree). An occultation might be
        visible from somewhere on Earth even if not from the geocenter.

        Returns: (is_occultation, separation, moon_radius, target_radius)
        """
        moon_ra, moon_dec, moon_dist, moon_radius = _get_moon_position(jd)
        target_ra, target_dec, target_radius = _get_target_position(jd)

        separation = _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

        # Occultation threshold: Moon radius + target radius + parallax margin
        # Lunar parallax can be up to ~1 degree, so we add margin to detect
        # global occultations that might only be visible from some Earth locations
        LUNAR_PARALLAX_MARGIN = 1.0  # degrees
        threshold = moon_radius + target_radius + LUNAR_PARALLAX_MARGIN

        return separation <= threshold, separation, moon_radius, target_radius

    def _find_minimum_separation(jd_start_search: float, jd_end_search: float) -> float:
        """Find time of minimum separation using golden section search.

        Uses increased iterations for sub-second precision needed for
        accurate occultation timing.
        """
        phi = (1 + math.sqrt(5)) / 2  # Golden ratio

        a = jd_start_search
        b = jd_end_search

        c = b - (b - a) / phi
        d = a + (b - a) / phi

        def get_sep(jd: float) -> float:
            moon_ra, moon_dec, _, _ = _get_moon_position(jd)
            target_ra, target_dec, _ = _get_target_position(jd)
            return _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

        fc = get_sep(c)
        fd = get_sep(d)

        for _ in range(60):  # Increased iterations for sub-second precision
            if fc < fd:
                b = d
                d = c
                fd = fc
                c = b - (b - a) / phi
                fc = get_sep(c)
            else:
                a = c
                c = d
                fc = fd
                d = a + (b - a) / phi
                fd = get_sep(d)

            if b - a < 1e-7:  # ~0.01 second precision
                break

        return (a + b) / 2

    def _get_separation_at_jd(jd: float) -> float:
        """Get angular separation between Moon and target at given JD."""
        moon_ra, moon_dec, _, _ = _get_moon_position(jd)
        target_ra, target_dec, _ = _get_target_position(jd)
        return _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

    def _find_contact_time_bisection(
        jd_start: float,
        jd_end: float,
        target_sep: float,
        is_increasing: bool,
    ) -> float:
        """Find time when separation equals target using bisection search.

        This is more accurate than chord-based estimates for contact times.

        Args:
            jd_start: Start of search range
            jd_end: End of search range
            target_sep: Target angular separation (threshold)
            is_increasing: True if looking for separation increasing through target

        Returns:
            Julian day when separation equals target, or 0.0 if not found
        """
        jd_low = jd_start
        jd_high = jd_end

        # Check if target separation is achievable in search range
        sep_start = _get_separation_at_jd(jd_start)
        sep_end = _get_separation_at_jd(jd_end)

        if is_increasing:
            # For increasing: sep should go from below to above target
            if sep_start > target_sep and sep_end > target_sep:
                return 0.0  # Never crosses target
        else:
            # For decreasing: sep should go from above to below target
            if sep_start < target_sep and sep_end < target_sep:
                return 0.0  # Never crosses target

        for _ in range(50):  # Converge to ~0.01 second precision
            jd_mid = (jd_low + jd_high) / 2
            sep = _get_separation_at_jd(jd_mid)

            if abs(sep - target_sep) < 1e-7:  # ~0.4 arcsec precision
                return jd_mid

            if is_increasing:
                # Looking for separation increasing through target
                if sep < target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid
            else:
                # Looking for separation decreasing through target
                if sep > target_sep:
                    jd_low = jd_mid
                else:
                    jd_high = jd_mid

            if jd_high - jd_low < 1e-8:  # ~0.001 second precision
                break

        result = (jd_low + jd_high) / 2
        # Validate the result
        final_sep = _get_separation_at_jd(result)
        if abs(final_sep - target_sep) > 0.0001:  # ~0.4 arcsec tolerance
            return 0.0  # Contact time not reliable

        return result

    def _calculate_contact_times(
        jd_max: float, min_sep: float, moon_radius: float, target_radius: float
    ) -> Tuple[float, float, float, float]:
        """Calculate contact times for the occultation using bisection search.

        Uses precise bisection search to find contact times, similar to solar
        eclipse calculations. This is more accurate than chord-based estimates.
        """
        # First/fourth contact: separation = moon_radius + target_radius
        outer_threshold = moon_radius + target_radius

        # Second/third contact: separation = |moon_radius - target_radius|
        # (for total occultation when target is inside Moon's disk)
        inner_threshold = abs(moon_radius - target_radius)

        # Search range for contacts (2 hours before/after maximum)
        contact_search = 2.0 / 24.0

        jd_first = 0.0
        jd_fourth = 0.0
        jd_second = 0.0
        jd_third = 0.0

        if min_sep < outer_threshold:
            # Find first contact (separation decreasing to threshold before max)
            jd_first = _find_contact_time_bisection(
                jd_max - contact_search, jd_max, outer_threshold, is_increasing=False
            )

            # Find fourth contact (separation increasing to threshold after max)
            jd_fourth = _find_contact_time_bisection(
                jd_max, jd_max + contact_search, outer_threshold, is_increasing=True
            )

        # For total occultation (target fully inside Moon)
        if min_sep < inner_threshold:
            # Find second contact (separation decreasing to inner threshold)
            jd_second = _find_contact_time_bisection(
                jd_max - contact_search, jd_max, inner_threshold, is_increasing=False
            )

            # Find third contact (separation increasing to inner threshold)
            jd_third = _find_contact_time_bisection(
                jd_max, jd_max + contact_search, inner_threshold, is_increasing=True
            )

        return jd_first, jd_second, jd_third, jd_fourth

    def _calculate_contact_times_global(
        jd_max: float,
        min_sep: float,
        moon_radius: float,
        target_radius: float,
        parallax: float,
    ) -> Tuple[float, float, float, float]:
        """Calculate contact times for global occultations accounting for parallax.

        For global occultations, the geocentric separation may be larger than
        moon_radius + target_radius, but an occultation still occurs from some
        Earth location due to lunar parallax. The contact times represent when
        the occultation begins/ends at the optimal viewing location.

        Contact occurs when geo_sep = threshold + parallax.

        Args:
            jd_max: Time of minimum geocentric separation
            min_sep: Minimum geocentric separation (degrees)
            moon_radius: Moon's angular radius (degrees)
            target_radius: Target's angular radius (degrees)
            parallax: Lunar horizontal parallax (degrees)

        Returns:
            Tuple of (jd_first, jd_second, jd_third, jd_fourth)
        """
        outer_threshold = moon_radius + target_radius
        inner_threshold = abs(moon_radius - target_radius)
        contact_search = 6.0 / 24.0  # 6 hours search range for global events

        jd_first = 0.0
        jd_fourth = 0.0
        jd_second = 0.0
        jd_third = 0.0

        # Contact times: geo_sep = threshold + parallax
        target_sep_outer = outer_threshold + parallax
        target_sep_inner = inner_threshold + parallax

        def find_contact(
            start_jd: float, end_jd: float, target_sep: float, is_increasing: bool
        ) -> float:
            """Find contact time using bisection on geocentric separation."""
            jd_low, jd_high = start_jd, end_jd

            for _ in range(60):
                jd_mid = (jd_low + jd_high) / 2
                moon_ra, moon_dec, _, _ = _get_moon_position(jd_mid)
                target_ra, target_dec, _ = _get_target_position(jd_mid)
                geo_sep = _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

                if abs(geo_sep - target_sep) < 1e-7:
                    return jd_mid

                if is_increasing:
                    if geo_sep < target_sep:
                        jd_low = jd_mid
                    else:
                        jd_high = jd_mid
                else:
                    if geo_sep > target_sep:
                        jd_low = jd_mid
                    else:
                        jd_high = jd_mid

                if jd_high - jd_low < 1e-8:
                    break

            return (jd_low + jd_high) / 2

        # Find first contact (before max, separation decreasing)
        jd_first = find_contact(
            jd_max - contact_search, jd_max, target_sep_outer, is_increasing=False
        )

        # Find fourth contact (after max, separation increasing)
        jd_fourth = find_contact(
            jd_max, jd_max + contact_search, target_sep_outer, is_increasing=True
        )

        # For total occultation, also find inner contacts
        if min_sep < inner_threshold + parallax:
            jd_second = find_contact(
                jd_max - contact_search, jd_max, target_sep_inner, is_increasing=False
            )
            jd_third = find_contact(
                jd_max, jd_max + contact_search, target_sep_inner, is_increasing=True
            )

        return jd_first, jd_second, jd_third, jd_fourth

    # --- Vectorized batch search (~20-50x speedup) ---
    # The batch search uses Skyfield vectorized calls (earth.at(t_array)).
    # When reader is None (Skyfield-only), eph/ts/earth/moon_body are
    # already defined from the dispatch above.  When reader is not None,
    # we need them for the batch scan.
    if reader is not None:
        # Batch scan requires Skyfield; load only if not already loaded
        from .state import get_planets as _gp_batch
        eph = _gp_batch()
        ts = get_timescale()
        earth = eph["earth"]
        moon_body = eph["moon"]

    # Pre-resolve target once for batch computation
    if planet != 0 and planet not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {planet}.")

    if planet == 0:
        from .fixed_stars import FIXED_STARS as _FS_B
        from .fixed_stars import _resolve_star_id as _rsi_b
        from skyfield.api import Star as _SkyfieldStarBatch

        _b_star_id, _b_err, _ = _rsi_b(star_name)
        if _b_err is not None:
            raise ValueError(_b_err)
        _b_star = _FS_B[_b_star_id]
        # Build a Skyfield Star object for vectorized epoch-of-date positions
        _b_star_sf = _SkyfieldStarBatch(
            ra_hours=_b_star.ra_j2000 / 15.0,
            dec_degrees=_b_star.dec_j2000,
            ra_mas_per_year=_b_star.pm_ra * 1000.0,
            dec_mas_per_year=_b_star.pm_dec * 1000.0,
        )
        _b_target = None
    else:
        _b_star = None
        # Use the plain Skyfield ephemeris body for vectorized batch calls.
        # get_planet_target() returns custom wrappers (_CobCorrectedTarget,
        # _SpkCenterTarget) whose _observe_from_bcrs doesn't support array
        # times.  The COB correction is < 0.01" — irrelevant for the coarse
        # candidate scan.  Fine refinement later uses the scalar functions.
        from .planets import _PLANET_FALLBACK

        _planet_name_sf = _PLANET_MAP[planet]
        try:
            _b_target = eph[_planet_name_sf]
        except KeyError:
            _b_target = eph[_PLANET_FALLBACK.get(_planet_name_sf, _planet_name_sf)]

    def _batch_separations(
        jd_arr: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Compute Moon-target angular separations for an array of JDs.

        Returns ``(separations, moon_radii, target_radii, moon_dist_au)``
        as numpy arrays of the same length as *jd_arr*.
        """
        jds = np.asarray(jd_arr, dtype=np.float64)
        t = ts.ut1_jd(jds)

        # Moon geocentric apparent position
        moon_app = earth.at(t).observe(moon_body).apparent()
        m_ra_obj, m_dec_obj, m_dist_obj = moon_app.radec(epoch="date")
        m_ra = np.asarray(m_ra_obj.hours) * 15.0
        m_dec = np.asarray(m_dec_obj.degrees)
        m_dist = np.asarray(m_dist_obj.au)
        m_radii = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / m_dist)

        # Target position
        if _b_star is not None:
            # Use Skyfield Star for proper precession to epoch of date,
            # matching the Moon's epoch="date" frame.
            star_app = earth.at(t).observe(_b_star_sf).apparent()
            t_ra_o, t_dec_o, _ = star_app.radec(epoch="date")
            t_ra = np.asarray(t_ra_o.hours) * 15.0
            t_dec = np.asarray(t_dec_o.degrees)
            t_radii = np.full(len(jds), 0.0001)
        else:
            t_app = earth.at(t).observe(_b_target).apparent()
            t_ra_o, t_dec_o, t_dist_o = t_app.radec(epoch="date")
            t_ra = np.asarray(t_ra_o.hours) * 15.0
            t_dec = np.asarray(t_dec_o.degrees)
            t_dist = np.asarray(t_dist_o.au)
            t_radii = np.array(
                [_calc_planet_angular_radius(planet, float(d)) for d in t_dist]
            )

        # Vectorized haversine angular separation
        ra1, dec1 = np.radians(m_ra), np.radians(m_dec)
        ra2, dec2 = np.radians(t_ra), np.radians(t_dec)
        cos_sep = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(
            ra1 - ra2
        )
        np.clip(cos_sep, -1.0, 1.0, out=cos_sep)
        seps = np.degrees(np.arccos(cos_sep))

        return seps, m_radii, t_radii, m_dist

    # Coarse-to-fine batch search
    _SCAN_STEP = 0.5  # days  (Moon moves ~13°/day; 0.5d catches all passes)
    _SCAN_CHUNK = 1000  # JDs per vectorised call
    _CANDIDATE_DEG = 5.0  # wide net for coarse candidates

    direction = -1.0 if backwards else 1.0
    search_limit = MAX_SEARCH_YEARS * 365.25
    jd_cursor = jd_start

    # Clamp search to the loaded DE kernel range to avoid EphemerisRangeError
    _eph_jd_min, _eph_jd_max = _kernel_jd_bounds(eph)

    while abs(jd_cursor - jd_start) < search_limit:
        chunk_end = jd_cursor + direction * _SCAN_CHUNK * _SCAN_STEP
        # Clamp to search boundary
        if direction > 0:
            chunk_end = min(chunk_end, jd_start + search_limit, _eph_jd_max)
        else:
            chunk_end = max(chunk_end, jd_start - search_limit, _eph_jd_min)

        chunk_jds = np.arange(jd_cursor, chunk_end, direction * _SCAN_STEP)
        if len(chunk_jds) == 0:
            break

        seps, m_radii, t_radii, m_dists = _batch_separations(chunk_jds)

        # Candidate selection: separation below per-point occultation
        # threshold OR below a fixed wide threshold
        occ_thresh = m_radii + t_radii + 1.0  # + parallax margin
        candidate_mask = seps < np.maximum(occ_thresh, _CANDIDATE_DEG)

        if not np.any(candidate_mask):
            jd_cursor = float(chunk_jds[-1]) + direction * _SCAN_STEP
            continue

        # Group adjacent candidate indices into distinct events
        cand_idx = np.where(candidate_mask)[0]
        events: list[list[int]] = []
        cur_group: list[int] = [int(cand_idx[0])]
        for ci in range(1, len(cand_idx)):
            if cand_idx[ci] - cand_idx[ci - 1] <= 5:
                cur_group.append(int(cand_idx[ci]))
            else:
                events.append(cur_group)
                cur_group = [int(cand_idx[ci])]
        events.append(cur_group)

        for event_group in events:
            # Best coarse point in this group
            g_seps = seps[event_group]
            best_i = event_group[int(np.argmin(g_seps))]
            jd_near = float(chunk_jds[best_i])

            # Fine batch refinement around candidate (±0.5 day, 0.02 step)
            fine_jds = np.arange(jd_near - 0.5, jd_near + 0.5, 0.02)
            fine_seps, _, _, _ = _batch_separations(fine_jds)
            jd_refined = float(fine_jds[int(np.argmin(fine_seps))])

            # Precise minimum via golden-section search (scalar calls)
            jd_max = _find_minimum_separation(jd_refined - 0.25, jd_refined + 0.25)

            # Reject events whose maximum falls before the search start
            # (forward) or after (backward).  The fine refinement window
            # (±0.5 day) can drag jd_max past jd_start, producing results
            # that predate the caller's requested start time.
            if direction > 0 and jd_max < jd_start:
                continue
            if direction < 0 and jd_max > jd_start:
                continue

            # --- Verification (identical to original logic) ---
            moon_ra, moon_dec, moon_dist, moon_r = _get_moon_position(jd_max)
            target_ra, target_dec, target_r = _get_target_position(jd_max)
            min_sep = _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

            moon_dist_km = moon_dist * 149597870.7
            LUNAR_PARALLAX = math.degrees(math.asin(6378.0 / moon_dist_km))

            if min_sep >= (moon_r + target_r + LUNAR_PARALLAX):
                continue  # not a true occultation

            effective_sep = max(0, min_sep - LUNAR_PARALLAX)

            jd_first, jd_second, jd_third, jd_fourth = _calculate_contact_times_global(
                jd_max, min_sep, moon_r, target_r, LUNAR_PARALLAX
            )

            grazing_threshold = 0.9 * moon_r
            is_grazing = effective_sep > grazing_threshold

            if effective_sep < abs(moon_r - target_r):
                ecl_type = ECL_ANNULAR if target_r > moon_r else ECL_TOTAL
            else:
                ecl_type = ECL_PARTIAL


            tret = (
                jd_max,
                0.0,
                jd_first,
                jd_fourth,
                jd_second,
                jd_third,
                0.0,
                0.0,
                0.0,
                0.0,
            )
            return ecl_type, tret

        jd_cursor = float(chunk_jds[-1]) + direction * _SCAN_STEP

    target_desc = star_name if planet == 0 else f"planet {planet}"
    raise RuntimeError(
        f"No lunar occultation of {target_desc} found within "
        f"{MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for reference API compatibility
lun_occult_when_glob = lun_occult_when_glob


def _lun_occult_when_loc_pythonic(
    tjdut: float,
    ipl: int,
    starname: str,
    lat: float,
    lon: float,
    alt: float = 0.0,
    flags: int = FLG_SWIEPH,
    backwards: bool = False,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next lunar occultation visible from a specific location.

    A lunar occultation occurs when the Moon passes in front of (occults)
    a planet or star as seen from Earth. This function searches forward
    (or backward) in time to find the next occultation visible from a
    specific geographic location, where both the Moon and the target are
    above the horizon.

    This function matches the reference lun_occult_when_loc() API exactly.

    Args:
        tjdut: Julian Day (UT) to start search from
        ipl: Planet identifier (int). Use MERCURY, VENUS, etc.
            Use 0 if searching for a star occultation.
        starname: Star name (str) for star occultations, or empty string
            if searching for a planet occultation (e.g., "Regulus").
        lat: Geographic latitude in degrees (North positive)
        lon: Geographic longitude in degrees (East positive)
        alt: Altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)
        backwards: If True, search backward in time; if False, search forward.

    Returns:
        Tuple containing:
            - retflag: Occultation type flags bitmask (ECL_* constants)
                0: No occultation found
                ECL_TOTAL: Total occultation (body fully behind Moon)
                ECL_PARTIAL: Partial occultation (body partially behind Moon)
                ECL_VISIBLE: Occultation visible from location
                ECL_MAX_VISIBLE: Maximum visible from location
                ECL_1ST_VISIBLE: First contact visible
                ECL_2ND_VISIBLE: Second contact visible
                ECL_3RD_VISIBLE: Third contact visible
                ECL_4TH_VISIBLE: Fourth contact visible
            - tret: Tuple of 10 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation (minimum separation)
                [1]: Time of first contact (occultation begins)
                [2]: Time of second contact (full occultation begins, or 0)
                [3]: Time of third contact (full occultation ends, or 0)
                [4]: Time of fourth contact (occultation ends)
                [5-9]: Reserved (0)
            - attr: Tuple of 20 floats with occultation attributes:
                [0]: Fraction of target diameter covered by Moon (magnitude)
                [1]: Ratio of lunar diameter to target diameter
                [2]: Fraction of target disc covered by Moon (obscuration)
                [3]: Diameter of core shadow in km (0 for stars)
                [4]: Azimuth of target at maximum occultation (degrees)
                [5]: True altitude of target above horizon at maximum (degrees)
                [6]: Apparent altitude of target above horizon at maximum (degrees)
                [7]: Angular separation (elongation) at maximum (degrees)
                [8-19]: Reserved (0)

    Raises:
        ValueError: If neither planet nor star_name is specified

    Algorithm:
        1. Use lun_occult_when_glob to find next global occultation
        2. Calculate Moon and target altitude at observer's location
        3. If both bodies are below horizon during entire occultation,
           continue to next global occultation
        4. Return local phase times, visibility flags, and attributes

    Precision:
        Occultation times accurate to ~1 minute for most events.
        Visibility depends on accurate horizon calculations.

    Example:
        >>> # Find next occultation of Regulus visible from Rome
        >>> from libephemeris import julday, _lun_occult_when_loc_pythonic
        >>> jd = julday(2017, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> ocl_type, times, attr = _lun_occult_when_loc_pythonic(jd, 0, "Regulus",
        ...                                              rome_lat, rome_lon)
        >>> print(f"Occultation maximum at JD {times[0]:.5f}")
        >>> print(f"Moon altitude: {attr[5]:.1f}°")

    References:
        - Reference API: lun_occult_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from skyfield.api import wgs84

    from .constants import (
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
    )
    from .fixed_stars import FIXED_STARS, _resolve_star_id
    from .planets import _PLANET_MAP, get_planet_target
    from .state import get_planets, get_timescale

    planet = ipl
    star_name = starname
    altitude = alt

    if planet == 0 and not star_name:
        raise ValueError(
            "Either planet ID or star name must be specified for occultation search"
        )

    jd_start = tjdut

    MAX_SEARCH_YEARS = 50
    MAX_OCCULTATIONS = int(MAX_SEARCH_YEARS * 20)

    MOON_MEAN_ANGULAR_RADIUS_DEG = 932.56 / 3600.0

    # reader is provided by the caller (None forces Skyfield path)
    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        geopos = (lon, lat, altitude)

        def _get_moon_altaz(jd: float) -> Tuple[float, float]:
            jd_tt = jd + deltat(jd)
            pos = _topo_ecliptic(reader, jd_tt, jd, MOON, geopos, FLG_SPEED)
            az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, pos[:3])
            return alt_true, az

        def _get_target_altaz(jd: float) -> Tuple[float, float]:
            jd_tt = jd + deltat(jd)
            if planet == 0:
                from .fixed_stars import fixstar_ut
                star_pos, _, _ = fixstar_ut(star_name, jd, FLG_SPEED)
                az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, star_pos[:3])
                return alt_true, az
            else:
                pos = _topo_ecliptic(reader, jd_tt, jd, planet, geopos, FLG_SPEED)
                az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, pos[:3])
                return alt_true, az

        def _is_visible(jd: float, min_alt: float = -1.0) -> bool:
            moon_alt, _ = _get_moon_altaz(jd)
            target_alt, _ = _get_target_altaz(jd)
            return moon_alt > min_alt and target_alt > min_alt

        def _get_moon_angular_info(jd: float) -> Tuple[float, float, float]:
            moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)
            dist = moon_pos[2]
            moon_diameter = 2 * MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist)
            return moon_diameter, dist, MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist)

        def _get_target_angular_info(jd: float) -> float:
            if planet == 0:
                return 0.0001
            target_pos, _ = calc_ut(jd, planet, FLG_SPEED)
            dist = target_pos[2]
            return 2 * _calc_planet_angular_radius(planet, dist)

        def _angular_separation_at_jd(jd_check: float) -> float:
            """Topocentric angular separation between Moon and target (LEB path).

            Uses observer's actual position to account for lunar parallax,
            which can be up to ~1 degree and significantly affects timing.
            """
            jd_tt = jd_check + deltat(jd_check)
            moon_ecl = _topo_ecliptic(
                reader, jd_tt, jd_check, MOON, geopos, FLG_SPEED
            )
            if planet == 0:
                from .fixed_stars import fixstar_ut
                star_pos, _, _ = fixstar_ut(star_name, jd_check, FLG_SPEED)
                target_lon, target_lat = star_pos[0], star_pos[1]
            else:
                target_ecl = _topo_ecliptic(
                    reader, jd_tt, jd_check, planet, geopos, FLG_SPEED
                )
                target_lon, target_lat = target_ecl[0], target_ecl[1]
            return angular_separation(
                moon_ecl[0], moon_ecl[1], target_lon, target_lat
            )
    else:
        eph = get_planets()
        ts = get_timescale()
        earth = eph["earth"]
        moon_body = eph["moon"]
        observer = wgs84.latlon(lat, lon, altitude)

        def _get_moon_altaz(jd: float) -> Tuple[float, float]:
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            moon_app = observer_at.at(t).observe(moon_body).apparent()
            alt, az, _ = moon_app.altaz()
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_target_altaz(jd: float) -> Tuple[float, float]:
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            if planet == 0:
                star_id, err, _ = _resolve_star_id(star_name)
                if err is not None:
                    raise ValueError(err)
                star = FIXED_STARS[star_id]
                t_years = (jd - 2451545.0) / 365.25
                ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
                from skyfield.api import Star
                star_obj = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
                star_app = observer_at.at(t).observe(star_obj).apparent()
                alt, az, _ = star_app.altaz()
                return alt.degrees, (az.degrees + 180.0) % 360.0
            else:
                if planet not in _PLANET_MAP:
                    raise ValueError(f"illegal planet number {planet}.")
                target_name = _PLANET_MAP[planet]
                target = get_planet_target(eph, target_name)
                target_app = observer_at.at(t).observe(target).apparent()
                alt, az, _ = target_app.altaz()
                return alt.degrees, (az.degrees + 180.0) % 360.0

        def _is_visible(jd: float, min_alt: float = -1.0) -> bool:
            moon_alt, _ = _get_moon_altaz(jd)
            target_alt, _ = _get_target_altaz(jd)
            return moon_alt > min_alt and target_alt > min_alt

        def _get_moon_angular_info(jd: float) -> Tuple[float, float, float]:
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            moon_app = observer_at.at(t).observe(moon_body).apparent()
            dist = moon_app.distance().au
            moon_diameter = 2 * MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist)
            return moon_diameter, dist, MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist)

        def _get_target_angular_info(jd: float) -> float:
            if planet == 0:
                return 0.0001
            if planet not in _PLANET_MAP:
                return 0.0001
            target_name = _PLANET_MAP[planet]
            target = get_planet_target(eph, target_name)
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            target_app = observer_at.at(t).observe(target).apparent()
            dist = target_app.distance().au
            return 2 * _calc_planet_angular_radius(planet, dist)

        def _angular_separation_at_jd(jd_check: float) -> float:
            """Topocentric angular separation between Moon and target (Skyfield path).

            Uses observer's actual position to account for lunar parallax,
            which can be up to ~1 degree and significantly affects timing.
            """
            t = ts.ut1_jd(jd_check)
            observer_at = earth + observer
            moon_app = observer_at.at(t).observe(moon_body).apparent()

            if planet == 0:
                star_id, err, _ = _resolve_star_id(star_name)
                if err is not None:
                    raise ValueError(err)
                star = FIXED_STARS[star_id]
                t_years = (jd_check - 2451545.0) / 365.25
                ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0

                from skyfield.api import Star

                star_obj = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
                target_app = observer_at.at(t).observe(star_obj).apparent()
            else:
                target_name = _PLANET_MAP[planet]
                target = get_planet_target(eph, target_name)
                target_app = observer_at.at(t).observe(target).apparent()

            sep = moon_app.separation_from(target_app)
            return sep.degrees

    jd = jd_start

    for _ in range(MAX_OCCULTATIONS):
        # Find next global occultation
        try:
            _body = star_name if star_name else planet
            global_type, global_times = lun_occult_when_glob(
                jd, _body, flags, 0, backwards
            )
        except RuntimeError:
            raise RuntimeError(
                f"No lunar occultation of {'star ' + star_name if planet == 0 else 'planet ' + str(planet)} "
                f"visible from lat={lat}, lon={lon} "
                f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
            )
        except EphemerisRangeError:
            raise RuntimeError(
                f"No lunar occultation of {'star ' + star_name if planet == 0 else 'planet ' + str(planet)} "
                f"visible from lat={lat}, lon={lon} "
                f"within ephemeris range (search from JD {jd_start})"
            )

        # Extract times using reference API indices:
        # [0]: time of maximum, [2]: occultation begin, [3]: occultation end
        # [4]: totality begin, [5]: totality end
        jd_max = global_times[0]
        jd_first = global_times[2]  # occultation begin
        jd_second = global_times[4]  # totality begin
        jd_third = global_times[5]  # totality end
        jd_fourth = global_times[3]  # occultation end

        # Check if occultation is visible from this location
        # Visibility criteria for local occultation observation:
        # - The maximum must have both Moon and target above horizon
        # - At least one contact must be visible for the event to be "locally visible"
        max_visible = _is_visible(jd_max)

        # Check visibility at contact points
        first_visible = jd_first > 0 and _is_visible(jd_first)
        fourth_visible = jd_fourth > 0 and _is_visible(jd_fourth)

        # Sample at key phases for visibility flags
        check_times = [jd_max]
        if jd_first > 0:
            check_times.append(jd_first)
        if jd_fourth > 0:
            check_times.append(jd_fourth)
        if jd_second > 0:
            check_times.append(jd_second)
        if jd_third > 0:
            check_times.append(jd_third)

        # Also sample at intermediate points
        if jd_first > 0 and jd_fourth > 0:
            duration = jd_fourth - jd_first
            for i in range(5):
                t = jd_first + duration * i / 4
                if t not in check_times:
                    check_times.append(t)

        # Check if any phase is visible (for visibility flags)
        any_visible = any(_is_visible(t) for t in check_times)

        # Visibility criterion: event is observable if the fourth contact is visible
        # This ensures the complete occultation (ingress and egress) is observable
        # If fourth contact is not visible (e.g., Moon/planet sets during occultation),
        # skip this event and search for the next one
        if not fourth_visible:
            # Fourth contact not visible - incomplete event, try next
            jd = jd_max + (-1 if backwards else 1)  # Step to the next event
            continue

        # Occultation visible! Calculate local circumstances with topocentric parallax

        # Find local maximum using golden section search with topocentric positions
        def _find_local_maximum(
            jd_center: float, search_range: float = 3.0 / 24.0
        ) -> float:
            """Find local occultation maximum using golden section search."""
            phi = (1 + math.sqrt(5)) / 2

            a = jd_center - search_range
            b = jd_center + search_range

            c = b - (b - a) / phi
            d = a + (b - a) / phi

            fc = _angular_separation_at_jd(c)
            fd = _angular_separation_at_jd(d)

            for _ in range(60):  # High precision iterations
                if fc < fd:
                    b = d
                    d = c
                    fd = fc
                    c = b - (b - a) / phi
                    fc = _angular_separation_at_jd(c)
                else:
                    a = c
                    c = d
                    fc = fd
                    d = a + (b - a) / phi
                    fd = _angular_separation_at_jd(d)

                if b - a < 1e-7:  # ~0.01 second precision
                    break

            return (a + b) / 2

        # Find contact time using bisection with topocentric positions
        def _find_local_contact_time(
            jd_start_search: float,
            jd_end_search: float,
            target_sep: float,
            is_increasing: bool,
        ) -> float:
            """Find time when topocentric separation equals target using bisection."""
            jd_low = jd_start_search
            jd_high = jd_end_search

            # Check if target separation is achievable in search range
            sep_start = _angular_separation_at_jd(jd_start_search)
            sep_end = _angular_separation_at_jd(jd_end_search)

            if is_increasing:
                if sep_start > target_sep and sep_end > target_sep:
                    return 0.0
            else:
                if sep_start < target_sep and sep_end < target_sep:
                    return 0.0

            for _ in range(50):
                jd_mid = (jd_low + jd_high) / 2
                sep = _angular_separation_at_jd(jd_mid)

                if abs(sep - target_sep) < 1e-7:
                    return jd_mid

                if is_increasing:
                    if sep < target_sep:
                        jd_low = jd_mid
                    else:
                        jd_high = jd_mid
                else:
                    if sep > target_sep:
                        jd_low = jd_mid
                    else:
                        jd_high = jd_mid

                if jd_high - jd_low < 1e-8:
                    break

            result = (jd_low + jd_high) / 2
            final_sep = _angular_separation_at_jd(result)
            if abs(final_sep - target_sep) > 0.0001:
                return 0.0

            return result

        # Recalculate local maximum using topocentric positions
        # This accounts for lunar parallax (~1 degree) which affects timing
        jd_max_local = _find_local_maximum(jd_max)

        # Early check: if local maximum is not visible, skip to next event
        if not _is_visible(jd_max_local):
            jd = jd_max + (-1 if backwards else 1)
            continue

        # Get Moon's alt/az at local maximum
        moon_alt, moon_az = _get_moon_altaz(jd_max_local)
        target_alt, target_az = _get_target_altaz(jd_max_local)

        # Get apparent diameters
        moon_diameter, moon_dist, moon_radius = _get_moon_angular_info(jd_max_local)
        target_diameter = _get_target_angular_info(jd_max_local)
        target_radius = target_diameter / 2

        # Get minimum separation at local maximum
        min_separation = _angular_separation_at_jd(jd_max_local)

        # Calculate contact times using topocentric positions
        outer_threshold = moon_radius + target_radius
        inner_threshold = abs(moon_radius - target_radius)
        contact_search = 2.0 / 24.0  # 2 hours search range

        # Early check: if topocentric separation is too large, no occultation at this location
        if min_separation >= outer_threshold:
            jd = jd_max + (-1 if backwards else 1)
            continue

        jd_first_local = 0.0
        jd_fourth_local = 0.0
        jd_second_local = 0.0
        jd_third_local = 0.0

        if min_separation < outer_threshold:
            # Find first contact (separation decreasing to threshold before max)
            jd_first_local = _find_local_contact_time(
                jd_max_local - contact_search,
                jd_max_local,
                outer_threshold,
                is_increasing=False,
            )

            # Find fourth contact (separation increasing to threshold after max)
            jd_fourth_local = _find_local_contact_time(
                jd_max_local,
                jd_max_local + contact_search,
                outer_threshold,
                is_increasing=True,
            )

        # For total occultation
        if min_separation < inner_threshold:
            jd_second_local = _find_local_contact_time(
                jd_max_local - contact_search,
                jd_max_local,
                inner_threshold,
                is_increasing=False,
            )
            jd_third_local = _find_local_contact_time(
                jd_max_local,
                jd_max_local + contact_search,
                inner_threshold,
                is_increasing=True,
            )

        # Calculate fraction covered at maximum
        if min_separation < abs(moon_radius - target_radius):
            fraction_covered = 1.0
        elif min_separation < moon_radius + target_radius:
            overlap = (moon_radius + target_radius) - min_separation
            fraction_covered = min(1.0, overlap / target_diameter)
        else:
            fraction_covered = 0.0

        # Determine occultation type
        ecl_type = global_type

        # If there's no actual occultation at this location, skip to next global event
        if jd_first_local <= 0 or jd_fourth_local <= 0:
            jd = jd_max + (-1 if backwards else 1)
            continue

        # Check visibility at each phase using local times
        if jd_first_local > 0 and _is_visible(jd_first_local):
            ecl_type |= ECL_1ST_VISIBLE
        if jd_second_local > 0 and _is_visible(jd_second_local):
            ecl_type |= ECL_2ND_VISIBLE
        if jd_third_local > 0 and _is_visible(jd_third_local):
            ecl_type |= ECL_3RD_VISIBLE
        if jd_fourth_local > 0 and _is_visible(jd_fourth_local):
            ecl_type |= ECL_4TH_VISIBLE
        if _is_visible(jd_max_local):
            ecl_type |= ECL_MAX_VISIBLE
        ecl_type |= ECL_VISIBLE

        # Require the fourth contact to be visible for a complete observation
        # If fourth contact is not visible at this location, skip to next global event
        fourth_visible = jd_fourth_local > 0 and _is_visible(jd_fourth_local)
        if not fourth_visible:
            jd = jd_max + (-1 if backwards else 1)
            continue

        # Find moonrise/moonset during occultation
        moonrise_time = 0.0
        moonset_time = 0.0

        if jd_first_local > 0 and jd_fourth_local > 0:
            alt_begin, _ = _get_moon_altaz(jd_first_local)
            alt_end, _ = _get_moon_altaz(jd_fourth_local)
            if alt_begin < 0 < alt_end:
                # Moon rises during occultation - find approximate time
                t_low, t_high = jd_first_local, jd_fourth_local
                for _ in range(30):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _ = _get_moon_altaz(t_mid)
                    if alt_mid < 0:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                moonrise_time = (t_low + t_high) / 2
            elif alt_begin > 0 > alt_end:
                # Moon sets during occultation
                t_low, t_high = jd_first_local, jd_fourth_local
                for _ in range(30):
                    t_mid = (t_low + t_high) / 2
                    alt_mid, _ = _get_moon_altaz(t_mid)
                    if alt_mid > 0:
                        t_low = t_mid
                    else:
                        t_high = t_mid
                moonset_time = (t_low + t_high) / 2

        # Prepare times tuple (10 elements) using LOCAL topocentric times
        times = (
            jd_max_local,  # [0] Maximum (topocentric)
            jd_first_local,  # [1] First contact (topocentric)
            jd_second_local,  # [2] Second contact (topocentric)
            jd_third_local,  # [3] Third contact (topocentric)
            jd_fourth_local,  # [4] Fourth contact (topocentric)
            0.0,  # [5] Reserved
            0.0,  # [6] Reserved
            moonrise_time,  # [7] Moonrise during occultation
            moonset_time,  # [8] Moonset during occultation
            0.0,  # [9] Reserved
        )

        # Prepare attributes tuple (20 elements, reference API compatible)
        # Calculate ratio of diameters
        diameter_ratio = moon_diameter / target_diameter if target_diameter > 0 else 0.0

        attr = (
            fraction_covered,  # [0] Fraction of target diameter covered (magnitude)
            diameter_ratio,  # [1] Ratio of lunar diameter to target diameter
            fraction_covered,  # [2] Fraction of disc covered (obscuration, same as magnitude for point sources)
            0.0,  # [3] Diameter of core shadow in km (0 for stars)
            target_az,  # [4] Azimuth of target at maximum
            target_alt,  # [5] True altitude of target above horizon
            target_alt,  # [6] Apparent altitude (same as true for simplicity)
            min_separation,  # [7] Angular separation (elongation) at maximum
            0.0,  # [8] Reserved
            0.0,  # [9] Reserved
            0.0,  # [10] Reserved
            0.0,  # [11] Reserved
            0.0,  # [12] Reserved
            0.0,  # [13] Reserved
            0.0,  # [14] Reserved
            0.0,  # [15] Reserved
            0.0,  # [16] Reserved
            0.0,  # [17] Reserved
            0.0,  # [18] Reserved
            0.0,  # [19] Reserved
        )

        return ecl_type, times, attr

    target_desc = star_name if planet == 0 else f"planet {planet}"
    raise RuntimeError(
        f"No lunar occultation of {target_desc} visible from lat={lat}, lon={lon} "
        f"within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


def lun_occult_when_loc(
    tjdut: float,
    body: "int | str",
    geopos: "Sequence[float]",
    flags: int = FLG_SWIEPH,
    backwards: bool = False,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next lunar occultation visible from a specific location.

    This function matches the reference lun_occult_when_loc() API exactly.

    A lunar occultation occurs when the Moon passes in front of (occults)
    a planet or star as seen from Earth. This function searches forward
    (or backward) in time to find the next occultation visible from a specific
    geographic location, where both the Moon and the target are above the horizon.

    Args:
        tjdut: Julian Day (UT) to start search from
        body: Planet identifier (int) or star name (str)
            For planets: MERCURY, VENUS, MARS, JUPITER, SATURN, etc.
            For stars: e.g., "Regulus", "Spica", "Aldebaran"
        geopos: Sequence of [longitude_degrees, latitude_degrees, altitude_meters]
                NOTE: longitude comes first (this matches reference API convention)
        flags: Calculation flags (FLG_SWIEPH, etc.)
        backwards: If True, search backward in time instead of forward

    Returns:
        Tuple containing:
            - retflags: Occultation type flags bitmask (int):
                ECL_TOTAL: Total occultation (body fully behind Moon)
                ECL_PARTIAL: Partial occultation (body partially behind Moon)
                ECL_VISIBLE: Occultation visible from location
                ECL_MAX_VISIBLE: Maximum visible from location
                ECL_1ST_VISIBLE: First contact visible
                ECL_4TH_VISIBLE: Fourth contact visible
            - tret: Tuple of 10 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation (minimum separation)
                [1]: Time of first contact (occultation begins)
                [2]: Time of second contact (full occultation begins, or 0)
                [3]: Time of third contact (full occultation ends, or 0)
                [4]: Time of fourth contact (occultation ends)
                [5]: Reserved (0)
                [6]: Reserved (0)
                [7]: Time of moonrise (if Moon rises during occultation, else 0)
                [8]: Time of moonset (if Moon sets during occultation, else 0)
                [9]: Reserved (0)
            - attr: Tuple of 20 floats with occultation attributes:
                [0]: Fraction of target diameter covered by Moon (magnitude)
                [1]: Ratio of lunar diameter to target diameter
                [2]: Fraction of target disc covered by Moon (obscuration)
                [3]: Diameter of core shadow in km (0 for stars)
                [4]: Azimuth of target at maximum occultation (degrees)
                [5]: True altitude of target above horizon at maximum (degrees)
                [6]: Apparent altitude of target above horizon at maximum (degrees)
                [7]: Angular separation (elongation) at maximum (degrees)
                [8-19]: Reserved (0)

    Raises:
        RuntimeError: If no occultation visible from location within search limit
        ValueError: If body is invalid or geopos has wrong length

    Example:
        >>> # Find next occultation of Regulus visible from Rome
        >>> from libephemeris import julday, lun_occult_when_loc, FLG_SWIEPH
        >>> jd = julday(2017, 1, 1, 0)
        >>> rome_geopos = [12.4964, 41.9028, 0]  # lon, lat, alt
        >>> retflags, tret, attr = lun_occult_when_loc(jd, "Regulus", rome_geopos)
        >>> print(f"Occultation maximum at JD {tret[0]:.5f}")
        >>> print(f"Moon altitude: {attr[5]:.1f}°")

    References:
        - Reference API: lun_occult_when_loc()
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from typing import Sequence

    # Validate geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    # Extract geographic position (reference API uses lon, lat, alt order)
    lon = geopos[0]
    lat = geopos[1]
    altitude = geopos[2]

    # Determine if body is planet ID or star name
    if isinstance(body, str):
        planet = 0
        star_name = body
    else:
        planet = body
        star_name = ""

    # Call the internal implementation with LEB→Skyfield fallback
    ecl_type, times, attr = _call_with_leb_skyfield_fallback(
        _lun_occult_when_loc_pythonic,
        tjdut, planet, star_name, lat, lon, altitude, flags, backwards,
    )

    # Return in reference API order: (retflags, tret, attr)
    return ecl_type, times, attr


def _lun_occult_where_pythonic(
    tjdut: float,
    body: Union[int, str],
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Calculate where on Earth a lunar occultation is visible at a given time.

    This function determines where on Earth the lunar occultation of a planet
    or star is visible at the specified Julian Day. It returns the geographic
    coordinates of the central line (where the occultation is most central)
    and attributes about the occultation geometry.

    This function matches the reference lun_occult_where() API exactly.

    Args:
        tjdut: Julian Day (UT) of the moment to calculate
        body: Planet ID (int) or star name (str) to check for occultation.
            For planets use MERCURY, VENUS, etc.
            For stars use the star name as a string (e.g. "Regulus").
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - retflag: Occultation type flags bitmask (ECL_* constants)
                Returns 0 if no occultation at this time
            - geopos: Tuple of 10 floats with geographic positions:
                [0]: Geographic longitude of central occultation (degrees, East+)
                [1]: Geographic latitude of central occultation (degrees, North+)
                [2]: Geographic longitude of northern limit (degrees)
                [3]: Geographic latitude of northern limit (degrees)
                [4]: Geographic longitude of southern limit (degrees)
                [5]: Geographic latitude of southern limit (degrees)
                [6]: Geographic longitude of northern penumbra limit (degrees)
                [7]: Geographic latitude of northern penumbra limit (degrees)
                [8]: Geographic longitude of southern penumbra limit (degrees)
                [9]: Geographic latitude of southern penumbra limit (degrees)
            - attr: Tuple of 20 floats with occultation attributes:
                [0]: Fraction of target covered by Moon (magnitude)
                [1]: Ratio of lunar diameter to target diameter
                [2]: Fraction of disc covered (obscuration)
                [3]: Diameter of core shadow in km
                [4]: Moon's azimuth at central line (degrees)
                [5]: True altitude of Moon at central line (degrees)
                [6]: Apparent altitude of Moon (with refraction)
                [7]: Angular distance Moon center from target center
                [8-19]: Reserved for future use

    Raises:
        ValueError: If neither planet ID nor star name is specified

    Note:
        If there is no occultation at the given time (Moon too far from the
        target), geopos will contain zeros and retflag will be 0.

    Algorithm:
        1. Calculate Moon and target positions at the given time
        2. Check if angular separation is small enough for occultation
        3. Use gradient descent to find point of minimum Moon-target separation
        4. Calculate path limits based on occultation geometry
        5. Calculate attributes at the central line location

    Example:
        >>> # Find where Regulus occultation is visible (reference API style)
        >>> from libephemeris import julday, _lun_occult_where_pythonic
        >>> jd = julday(2017, 6, 28, 10.0)  # During a known occultation
        >>> retflag, geopos, attr = _lun_occult_where_pythonic(jd, "Regulus")
        >>> print(f"Central line at lon={geopos[0]:.2f}, lat={geopos[1]:.2f}")
        >>>
        >>> # Or using planet ID for planet occultation
        >>> retflag, geopos, attr = _lun_occult_where_pythonic(jd, VENUS)

    References:
        - Reference API: lun_occult_where()
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from skyfield.api import wgs84

    from .constants import (
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
    )
    from .fixed_stars import FIXED_STARS, _resolve_star_id
    from .planets import _PLANET_MAP, get_planet_target
    from .state import get_planets, get_timescale

    # Handle the body parameter - can be int (planet ID) or str (star name)
    if isinstance(body, str):
        planet = 0
        star_name = body
    else:
        planet = body
        star_name = ""

    if planet == 0 and not star_name:
        raise ValueError("Planet ID or star name must be specified for occultation")

    jd = tjdut

    MOON_MEAN_ANGULAR_RADIUS_DEG = 932.56 / 3600.0
    zero_geopos = (0.0,) * 10
    zero_attr = (0.0,) * 20

    # reader is provided by the caller (None forces Skyfield path)
    _reader = reader
    if _reader is not None:
        from .constants import FLG_EQUATORIAL

        def _get_moon_geocentric(jd_calc: float) -> Tuple[float, float, float, float]:
            moon_eq, _ = calc_ut(jd_calc, MOON, FLG_EQUATORIAL | FLG_SPEED)
            dist_au = moon_eq[2]
            moon_ar = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist_au)
            return moon_eq[0], moon_eq[1], dist_au, moon_ar

        def _get_target_position(jd_calc: float) -> Tuple[float, float, float]:
            if planet == 0:
                from .fixed_stars import fixstar_ut
                pos, _, _ = fixstar_ut(star_name, jd_calc, FLG_EQUATORIAL | FLG_SPEED)
                return pos[0], pos[1], 0.0001
            else:
                if planet not in _PLANET_MAP:
                    raise ValueError(f"illegal planet number {planet}.")
                tgt_eq, _ = calc_ut(jd_calc, planet, FLG_EQUATORIAL | FLG_SPEED)
                ar = _calc_planet_angular_radius(planet, tgt_eq[2])
                return tgt_eq[0], tgt_eq[1], ar
    else:
        eph = get_planets()
        ts = get_timescale()
        earth = eph["earth"]
        moon_body = eph["moon"]

        def _get_moon_geocentric(jd_calc: float) -> Tuple[float, float, float, float]:
            t = ts.ut1_jd(jd_calc)
            moon_app = earth.at(t).observe(moon_body).apparent()
            ra, dec, dist = moon_app.radec(epoch="date")
            moon_ar = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / dist.au)
            return ra.hours * 15.0, dec.degrees, dist.au, moon_ar

        def _get_target_position(jd_calc: float) -> Tuple[float, float, float]:
            if planet == 0:
                star_id, err, _ = _resolve_star_id(star_name)
                if err is not None:
                    raise ValueError(err)
                star = FIXED_STARS[star_id]
                t_years = (jd_calc - 2451545.0) / 365.25
                ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
                return ra_deg, dec_deg, 0.0001
            else:
                if planet not in _PLANET_MAP:
                    raise ValueError(f"illegal planet number {planet}.")
                target_name = _PLANET_MAP[planet]
                target = get_planet_target(eph, target_name)
                t = ts.ut1_jd(jd_calc)
                target_app = earth.at(t).observe(target).apparent()
                ra, dec, dist = target_app.radec(epoch="date")
                angular_radius = _calc_planet_angular_radius(planet, dist.au)
                return ra.hours * 15.0, dec.degrees, angular_radius

    def _angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
        """Calculate angular separation between two points (in degrees)."""
        # Convert to radians
        ra1_r = math.radians(ra1)
        dec1_r = math.radians(dec1)
        ra2_r = math.radians(ra2)
        dec2_r = math.radians(dec2)

        # Haversine formula for spherical distance
        cos_sep = math.sin(dec1_r) * math.sin(dec2_r) + math.cos(dec1_r) * math.cos(
            dec2_r
        ) * math.cos(ra1_r - ra2_r)
        cos_sep = max(-1.0, min(1.0, cos_sep))
        return math.degrees(math.acos(cos_sep))

    # Get Moon and target positions
    moon_ra, moon_dec, moon_dist, moon_radius = _get_moon_geocentric(jd)
    target_ra, target_dec, target_radius = _get_target_position(jd)

    # Calculate angular separation
    separation = _angular_separation(moon_ra, moon_dec, target_ra, target_dec)

    # Check if occultation is occurring (somewhere on Earth)
    # Use parallax margin to detect occultations visible from some Earth locations
    # even if not from geocenter. Lunar parallax can be up to ~1 degree.
    LUNAR_PARALLAX_MARGIN = 1.0  # degrees
    occultation_threshold = moon_radius + target_radius + LUNAR_PARALLAX_MARGIN
    if separation > occultation_threshold:
        # No occultation at this time
        return 0, zero_geopos, zero_attr

    # Occultation is occurring - calculate where on Earth it's visible

    # Calculate GMST for initial guess
    _ts_ow = get_timescale()
    _t_ow = _ts_ow.ut1_jd(jd)
    gmst = _t_ow.gmst
    gmst_deg = gmst * 15.0

    init_lon = moon_ra - gmst_deg
    init_lon = ((init_lon + 180) % 360) - 180
    init_lat = moon_dec

    if _reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import angular_separation as _ang_sep_ow

        def get_separation(lat: float, lon: float) -> float:
            try:
                _gp = (lon, lat, 0.0)
                jd_tt = jd + deltat(jd)
                moon_pos = _topo_ecliptic(_reader, jd_tt, jd, MOON, _gp, FLG_SPEED)
                if planet == 0:
                    from .fixed_stars import fixstar_ut
                    tgt_pos, _, _ = fixstar_ut(star_name, jd, FLG_SPEED)
                else:
                    tgt_pos = _topo_ecliptic(_reader, jd_tt, jd, planet, _gp, FLG_SPEED)
                return _ang_sep_ow(moon_pos[0], moon_pos[1], tgt_pos[0], tgt_pos[1])
            except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
                _reraise_if_leb_range_error(_exc)
                return 999.0
    else:
        ts = _ts_ow
        t = _t_ow
        earth = eph["earth"]
        moon_body = eph["moon"]

        def _get_target_for_observer(observer_at, time_obj):
            if planet == 0:
                star_id, err, _ = _resolve_star_id(star_name)
                if err is not None:
                    raise ValueError(err)
                star = FIXED_STARS[star_id]
                t_years = (jd - 2451545.0) / 365.25
                ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
                from skyfield.api import Star
                star_obj = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
                return observer_at.at(time_obj).observe(star_obj).apparent()
            else:
                target_name = _PLANET_MAP[planet]
                target = get_planet_target(eph, target_name)
                return observer_at.at(time_obj).observe(target).apparent()

        def get_separation(lat: float, lon: float) -> float:
            try:
                observer = wgs84.latlon(lat, lon, 0.0)
                observer_at = earth + observer
                moon_app = observer_at.at(t).observe(moon_body).apparent()
                target_app = _get_target_for_observer(observer_at, t)
                return moon_app.separation_from(target_app).degrees
            except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
                _reraise_if_leb_range_error(_exc)
                return 999.0

    # Gradient descent to find minimum separation (occultation center)
    # This finds the point on Earth where Moon and target appear closest
    central_lat = init_lat
    central_lon = init_lon

    # Use a multi-scale search: start coarse, then refine
    for scale in [5.0, 1.0, 0.1, 0.01, 0.001]:
        best_sep = get_separation(central_lat, central_lon)

        # Search in a small grid around current best
        improved = True
        iterations = 0
        while improved and iterations < 50:
            improved = False
            iterations += 1

            for dlat, dlon in [
                (scale, 0),
                (-scale, 0),
                (0, scale),
                (0, -scale),
                (scale, scale),
                (scale, -scale),
                (-scale, scale),
                (-scale, -scale),
            ]:
                test_lat = max(-89.0, min(89.0, central_lat + dlat))
                test_lon = central_lon + dlon
                # Normalize longitude
                test_lon = ((test_lon + 180) % 360) - 180

                sep = get_separation(test_lat, test_lon)
                if sep < best_sep:
                    best_sep = sep
                    central_lat = test_lat
                    central_lon = test_lon
                    improved = True

    # Normalize longitude to -180 to +180
    central_lon = ((central_lon + 180) % 360) - 180

    # Clamp latitude to valid range
    central_lat = max(-90.0, min(90.0, central_lat))

    # Calculate occultation path width and limits
    # The occultation is visible where the Moon-target angular distance
    # as seen from the observer is less than the sum of their angular radii

    # Moon's parallax affects where the occultation is visible
    # The parallax is approximately: Moon distance in Earth radii ~ 60
    # Angular parallax ~ asin(1/60) ~ 0.95 degrees
    # Note: parallax calculation available for future refinement

    # Path width estimation based on Moon's angular diameter
    # The path width is approximately 2 * Moon_radius / cos(elevation) in degrees
    # Converting to km: 1 degree ~ 111 km at equator
    path_width_deg = 2 * moon_radius
    path_width_km = path_width_deg * 111.0  # Approximate at equator

    # Calculate northern and southern limits
    # These are approximately +/- Moon's angular radius in latitude from the central line
    # adjusted for the Moon's parallax effect
    north_lat = central_lat + moon_radius
    south_lat = central_lat - moon_radius

    # Clamp to valid latitude range
    north_lat = max(-90.0, min(90.0, north_lat))
    south_lat = max(-90.0, min(90.0, south_lat))

    # Northern and southern limits have approximately the same longitude
    # as the central line (simplified calculation)
    north_lon = central_lon
    south_lon = central_lon

    # Calculate sunrise/sunset limits
    # These are where the Moon is at the horizon during the occultation
    # Simplified: approximately ±90° in longitude from central line
    sunrise_lon = ((central_lon - 90) + 180) % 360 - 180
    sunset_lon = ((central_lon + 90) + 180) % 360 - 180
    sunrise_lat = central_lat
    sunset_lat = central_lat

    # Calculate attributes at central line
    try:
        if _reader is not None:
            from .time_utils import deltat
            from .utils import azalt, ECL2HOR

            _gp_c = (central_lon, central_lat, 0.0)
            _jd_tt_c = jd + deltat(jd)
            _moon_c = _topo_ecliptic(_reader, _jd_tt_c, jd, MOON, _gp_c, FLG_SPEED)
            _moon_az_c, moon_altitude, apparent_alt = azalt(jd, ECL2HOR, _gp_c, 1013.25, 15.0, _moon_c[:3])
            moon_azimuth = _moon_az_c

            if planet == 0:
                from .fixed_stars import fixstar_ut
                _tgt_c, _, _ = fixstar_ut(star_name, jd, FLG_SPEED)
            else:
                _tgt_c = _topo_ecliptic(_reader, _jd_tt_c, jd, planet, _gp_c, FLG_SPEED)
            local_separation = _ang_sep_ow(_moon_c[0], _moon_c[1], _tgt_c[0], _tgt_c[1])
            local_moon_dist = _moon_c[2]
        else:
            observer = wgs84.latlon(central_lat, central_lon, 0.0)
            observer_at = earth + observer

            moon_app = observer_at.at(t).observe(moon_body).apparent()
            target_app = _get_target_for_observer(observer_at, t)

            moon_alt, moon_az, _ = moon_app.altaz()
            moon_altitude = moon_alt.degrees
            moon_azimuth = (moon_az.degrees + 180.0) % 360.0

            if moon_altitude > -1.0:
                if moon_altitude > 0:
                    refraction = 1.0 / math.tan(
                        math.radians(moon_altitude + 7.31 / (moon_altitude + 4.4))
                    )
                    apparent_alt = moon_altitude + refraction / 60.0
                else:
                    apparent_alt = moon_altitude + 0.58
            else:
                apparent_alt = moon_altitude

            local_separation = moon_app.separation_from(target_app).degrees
            local_moon_dist = moon_app.distance().au
        local_moon_radius = MOON_MEAN_ANGULAR_RADIUS_DEG * (0.002569 / local_moon_dist)
        local_moon_diameter = 2 * local_moon_radius
        local_target_diameter = 2 * target_radius

        # Fraction covered (magnitude)
        if local_separation < abs(local_moon_radius - target_radius):
            fraction_covered = 1.0
        elif local_separation < local_moon_radius + target_radius:
            overlap = (local_moon_radius + target_radius) - local_separation
            fraction_covered = (
                min(1.0, overlap / (2 * target_radius)) if target_radius > 0 else 1.0
            )
        else:
            fraction_covered = 0.0

        # Ratio of lunar diameter to target diameter
        if local_target_diameter > 0:
            diameter_ratio = local_moon_diameter / local_target_diameter
        else:
            diameter_ratio = 999.0  # Moon much larger than star

    except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
        _reraise_if_leb_range_error(_exc)
        # If calculation fails, use defaults
        moon_azimuth = 0.0
        moon_altitude = 0.0
        apparent_alt = 0.0
        local_separation = separation
        local_moon_diameter = 2 * moon_radius
        local_target_diameter = 2 * target_radius
        fraction_covered = 0.0
        diameter_ratio = 0.0
        local_moon_radius = moon_radius

    # Determine occultation type based on local separation at central line
    if local_separation < abs(local_moon_radius - target_radius):
        if target_radius > local_moon_radius:
            # Target larger than Moon (very rare for occultations)
            eclipse_type = ECL_ANNULAR
        else:
            # Total occultation
            eclipse_type = ECL_TOTAL
    else:
        # Partial occultation
        eclipse_type = ECL_PARTIAL

    # Prepare return tuples
    geopos = (
        central_lon,  # [0] Central longitude
        central_lat,  # [1] Central latitude
        north_lon,  # [2] Northern limit longitude
        north_lat,  # [3] Northern limit latitude
        south_lon,  # [4] Southern limit longitude
        south_lat,  # [5] Southern limit latitude
        sunrise_lon,  # [6] Penumbra north limit longitude
        sunrise_lat,  # [7] Penumbra north limit latitude
        sunset_lon,  # [8] Penumbra south limit longitude
        sunset_lat,  # [9] Penumbra south limit latitude
    )

    attr = (
        fraction_covered,  # [0] Fraction covered (magnitude)
        diameter_ratio,  # [1] Ratio of lunar diameter to target diameter
        fraction_covered,  # [2] Fraction of disc covered (obscuration)
        path_width_km,  # [3] Diameter of core shadow in km
        moon_azimuth,  # [4] Moon azimuth at central line
        moon_altitude,  # [5] True altitude of Moon
        apparent_alt,  # [6] Apparent altitude (with refraction)
        local_separation,  # [7] Angular distance Moon center from target
        0.0,  # [8] Reserved
        0.0,  # [9] Reserved
        0.0,  # [10] Reserved
        0.0,  # [11] Reserved
        0.0,  # [12] Reserved
        0.0,  # [13] Reserved
        0.0,  # [14] Reserved
        0.0,  # [15] Reserved
        0.0,  # [16] Reserved
        0.0,  # [17] Reserved
        0.0,  # [18] Reserved
        0.0,  # [19] Reserved
    )

    return eclipse_type, geopos, attr


# Alias for internal API compatibility
_lun_occult_where_internal = _lun_occult_where_pythonic


def lun_occult_where(
    tjdut: float,
    body: "int | str" = 0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Calculate where on Earth a lunar occultation is visible (pyswisseph-compatible).

    Wrapper around _lun_occult_where_pythonic() matching pyswisseph's signature.

    Args:
        tjdut: Julian Day (UT) of the moment to calculate.
        body: Planet identifier (int) or star name (str).
        flags: Calculation flags (default FLG_SWIEPH).

    Returns:
        Tuple of (retflag, geopos, attr) matching pyswisseph.
    """
    return _call_with_leb_skyfield_fallback(
        _lun_occult_where_internal, tjdut, body, flags,
    )


# =============================================================================
# RISE, SET, AND TRANSIT CALCULATIONS
# =============================================================================


def rise_trans(
    tjdut: float,
    body: "Union[int, str]",
    rsmi: int,
    geopos: Sequence[float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate rise, set, or transit time for a celestial body.

    Wrapper around :func:`_rise_trans_impl` that adds LEB→Skyfield fallback
    for partial/custom LEB files. See the impl docstring for the full API.
    """
    return _call_with_leb_skyfield_fallback(
        _rise_trans_impl, tjdut, body, rsmi, geopos, atpress, attemp, flags,
    )


def _rise_trans_impl(
    tjdut: float,
    body: "Union[int, str]",
    rsmi: int,
    geopos: Sequence[float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate rise, set, or transit time for a celestial body.

    This function finds the next occurrence of a rising, setting, or transit
    (meridian passage) event for a celestial body as seen from a specific
    geographic location.

    Args:
        tjdut: Julian Day (UT) to start search from
        body: Planet/body ID (int, e.g. SUN, MOON) or fixed star
            name (str, e.g. "Sirius", "Regulus")
        rsmi: Event type and calculation flags (bitmask):
            - CALC_RISE (1): Rise time (body crossing horizon going up)
            - CALC_SET (2): Set time (body crossing horizon going down)
            - CALC_MTRANSIT (4): Upper meridian transit (culmination)
            - CALC_ITRANSIT (8): Lower meridian transit (anti-culmination)
            Additional flags (OR with event type):
            - BIT_DISC_CENTER (256): Use disc center instead of upper limb
            - BIT_DISC_BOTTOM (8192): Use lower limb of disc
            - BIT_NO_REFRACTION (512): Ignore atmospheric refraction
            - BIT_CIVIL_TWILIGHT (1024): Sun at -6 degrees
            - BIT_NAUTIC_TWILIGHT (2048): Sun at -12 degrees
            - BIT_ASTRO_TWILIGHT (4096): Sun at -18 degrees
        geopos: Geographic position as sequence [lon, lat, alt]:
            - 0: geographic longitude in degrees (eastern positive)
            - 1: geographic latitude in degrees (northern positive)
            - 2: geographic altitude in meters above sea level
        atpress: Atmospheric pressure in mbar/hPa for refraction (default 0.0)
        attemp: Atmospheric temperature in degrees Celsius (default 0.0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - res: integer flag (0 = event found, -2 = circumpolar)
            - tret: tuple of 10 floats, of which tret[0] = JD of event

    Raises:
        ValueError: If invalid planet ID or star name

    Note:
        For circumpolar objects (always above or below horizon at the given
        latitude), the function returns (-2, tret) with tret[0] = 0.0.
        For transits, circumpolar objects still have valid transit times.

    Algorithm:
        1. For transits: Find when body crosses the local meridian
           (Local Sidereal Time = body's Right Ascension)
        2. For rise/set: Find when body's altitude crosses the horizon
           accounting for refraction and disc size
        3. Uses Newton-Raphson iteration for precise timing

    Precision:
        Rise/set times accurate to ~1 minute for Sun/Moon (due to refraction
        uncertainty), better than 1 second for transit times.

    Example:
        >>> from libephemeris import julday, rise_trans, SUN, CALC_RISE
        >>> jd = julday(2024, 6, 21, 0)
        >>> # Find sunrise at Rome (geopos = [lon, lat, alt])
        >>> res, tret = rise_trans(jd, SUN, CALC_RISE, [12.5, 41.9, 0])
        >>> print(f"Sunrise at JD {tret[0]:.5f}")
        >>> # Find Sirius rise time
        >>> res, tret = rise_trans(jd, "Sirius", CALC_RISE, [12.5, 41.9, 0])

    References:
        - Reference API: rise_trans()
        - Meeus "Astronomical Algorithms" Ch. 15 (Rise, Set, Transit)
    """
    from skyfield.api import wgs84

    from .constants import (
        CALC_RISE,
        CALC_SET,
        CALC_MTRANSIT,
        CALC_ITRANSIT,
        BIT_DISC_CENTER,
        BIT_DISC_BOTTOM,
        BIT_NO_REFRACTION,
        BIT_CIVIL_TWILIGHT,
        BIT_NAUTIC_TWILIGHT,
        BIT_ASTRO_TWILIGHT,
    )
    from .planets import _PLANET_MAP, get_planet_target
    from .state import get_planets, get_timescale

    # Unpack geopos: [longitude, latitude, altitude]
    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2]) if len(geopos) > 2 else 0.0

    # Map parameter names for internal use
    jd_start = tjdut
    pressure = atpress
    temperature = attemp
    is_fixed_star = isinstance(body, str)

    # Extract event type from rsmi (lower bits)
    event_type = rsmi & 0x0F  # First 4 bits for event type

    # Validate event type
    if event_type not in (
        CALC_RISE,
        CALC_SET,
        CALC_MTRANSIT,
        CALC_ITRANSIT,
    ):
        raise ValueError(
            "Invalid event type in rsmi: %d. Use CALC_RISE, CALC_SET, CALC_MTRANSIT, or CALC_ITRANSIT"
            % rsmi
        )

    # reader is provided by the caller (None forces Skyfield path)
    _reader_rt = reader
    _use_leb_rt = _reader_rt is not None

    if _use_leb_rt:
        try:
            from . import fast_calc as _fc_rt
            from .time_utils import deltat as _sd_rt
            if not is_fixed_star:
                _fc_rt.fast_calc_ut(_reader_rt, jd_start, body, FLG_SPEED)
        except (KeyError, ValueError):
            _use_leb_rt = False

    if _use_leb_rt:
        from .constants import FLG_EQUATORIAL

        geopos_leb = (lon, lat, altitude)
        if is_fixed_star:
            planet = -1
        else:
            planet = body
            if planet not in _PLANET_MAP:
                raise ValueError("illegal planet number %d." % planet)
        ts = get_timescale()
        earth = None
        target = None
        observer = None
    else:
        eph = get_planets()
        ts = get_timescale()
        earth = eph["earth"]
        if is_fixed_star:
            from skyfield.api import Star as SkyfieldStar
            from .fixed_stars import FIXED_STARS, _resolve_star_id

            star_id, err, _ = _resolve_star_id(body)
            if err is not None:
                raise ValueError(err)
            star_entry = FIXED_STARS[star_id]
            t_years = (jd_start - 2451545.0) / 365.25
            ra_deg = star_entry.ra_j2000 + (star_entry.pm_ra * t_years) / 3600.0
            dec_deg = star_entry.dec_j2000 + (star_entry.pm_dec * t_years) / 3600.0
            target = SkyfieldStar(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
            planet = -1
        else:
            planet = body
            if planet not in _PLANET_MAP:
                raise ValueError("illegal planet number %d." % planet)
            target_name = _PLANET_MAP[planet]
            target = get_planet_target(eph, target_name)
        observer = wgs84.latlon(lat, lon, altitude)

    # Handle atpress=0: estimate pressure from altitude (pyswisseph convention).
    # When atpress is 0, Swiss Ephemeris estimates pressure using the
    # barometric formula. We replicate this behavior.
    if pressure == 0.0:
        # Barometric formula: P = 1013.25 * (1 - 2.25577e-5 * h)^5.25588
        pressure = 1013.25 * (1.0 - 2.25577e-5 * altitude) ** 5.25588
        if temperature == 0.0:
            # Standard temperature lapse rate: T = 15 - 0.0065 * h
            temperature = 15.0 - 0.0065 * altitude

    # Determine the altitude threshold for rise/set
    if rsmi & BIT_CIVIL_TWILIGHT:
        horizon_alt = -6.0
        use_refraction = False  # Twilight uses geometric position
    elif rsmi & BIT_NAUTIC_TWILIGHT:
        horizon_alt = -12.0
        use_refraction = False
    elif rsmi & BIT_ASTRO_TWILIGHT:
        horizon_alt = -18.0
        use_refraction = False
    else:
        horizon_alt = 0.0
        use_refraction = not (rsmi & BIT_NO_REFRACTION)

    # Calculate horizon refraction per IAU standard
    # Standard horizon refraction is 34 arcminutes (0.5667°), the conventional
    # value used in sunrise/sunset calculations (IAU/Meeus convention).
    # The Bennett (1982) formula gives ~29' but the traditional 34' is retained
    # for API compatibility with established ephemeris software.
    if use_refraction:
        # Standard horizon refraction: 34 arcminutes = 0.5667 degrees
        # Apply pressure and temperature correction (standard: 1010 mbar, 10°C = 283K)
        base_refraction = 34.0 / 60.0  # degrees
        pressure_factor = pressure / 1010.0
        temperature_factor = 283.0 / (273.0 + temperature)
        correction_factor = pressure_factor * temperature_factor
        refraction = base_refraction * correction_factor
    else:
        refraction = 0.0

    # Calculate semi-diameter dynamically based on body distance
    # Calculate actual angular size at observation time
    # We need to compute this at an initial estimate time
    if _use_leb_rt:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        _AU_KM = 149597870.7

        def _get_semi_diameter(jd: float, body: int) -> float:
            if body == SUN:
                sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
                return (959.63 / sun_pos[2]) / 3600.0
            elif body == MOON:
                moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)
                dist_km = moon_pos[2] * _AU_KM
                return math.degrees(math.atan(1737.4 / dist_km))
            return 0.0

        def _get_body_altaz(jd: float) -> Tuple[float, float]:
            jd_tt = jd + deltat(jd)
            if is_fixed_star:
                from .fixed_stars import fixstar_ut
                star_pos, _, _ = fixstar_ut(body, jd, FLG_SPEED)
                az, alt_true, _ = azalt(jd, ECL2HOR, geopos_leb, 0, 0, star_pos[:3])
            else:
                pos = _topo_ecliptic(_reader_rt, jd_tt, jd, planet, geopos_leb, FLG_SPEED)
                az, alt_true, _ = azalt(jd, ECL2HOR, geopos_leb, 0, 0, pos[:3])
            return alt_true, az

        def _get_body_ra_dec(jd: float) -> Tuple[float, float]:
            if is_fixed_star:
                from .fixed_stars import fixstar_ut
                pos, _, _ = fixstar_ut(body, jd, FLG_EQUATORIAL | FLG_SPEED)
                return pos[0] / 15.0, pos[1]
            else:
                eq, _ = calc_ut(jd, planet, FLG_EQUATORIAL | FLG_SPEED)
                return eq[0] / 15.0, eq[1]
    else:
        def _get_semi_diameter(jd: float, body: int) -> float:
            if body == SUN:
                t = ts.ut1_jd(jd)
                sun_pos = earth.at(t).observe(eph["sun"])
                distance_au = sun_pos.distance().au
                return (959.63 / distance_au) / 3600.0
            elif body == MOON:
                t = ts.ut1_jd(jd)
                moon_pos = earth.at(t).observe(eph["moon"])
                distance_km = moon_pos.distance().km
                return math.degrees(math.atan(1737.4 / distance_km))
            return 0.0

        def _get_body_altaz(jd: float) -> Tuple[float, float]:
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            body_app = observer_at.at(t).observe(target).apparent()
            alt, az, _ = body_app.altaz()
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_body_ra_dec(jd: float) -> Tuple[float, float]:
            t = ts.ut1_jd(jd)
            body_app = earth.at(t).observe(target).apparent()
            ra, dec, _ = body_app.radec(epoch="date")
            return ra.hours, dec.degrees

    initial_sd = _get_semi_diameter(jd_start, planet)

    is_twilight = bool(
        rsmi & (BIT_CIVIL_TWILIGHT | BIT_NAUTIC_TWILIGHT | BIT_ASTRO_TWILIGHT)
    )
    if is_twilight or (rsmi & BIT_DISC_CENTER):
        disc_correction = 0.0
    elif rsmi & BIT_DISC_BOTTOM:
        disc_correction = initial_sd
    else:
        disc_correction = -initial_sd

    target_altitude = horizon_alt - refraction + disc_correction

    if event_type in (CALC_MTRANSIT, CALC_ITRANSIT):
        if _use_leb_rt:
            return _calculate_transit_leb(
                jd_start, lon, event_type, _get_body_ra_dec, CALC_ITRANSIT,
            )
        return _calculate_transit(
            jd_start, lat, lon, event_type, ts, earth, target, observer,
            CALC_ITRANSIT,
        )

    return _calculate_rise_set(
        jd_start, lat, lon, event_type, target_altitude, ts, earth, target,
        observer, _get_body_altaz, _get_body_ra_dec, CALC_RISE, CALC_SET,
        rsmi,
    )


def _make_tret(jd_event: float = 0.0) -> Tuple[float, ...]:
    """Build a 10-element tret tuple for rise_trans return value.

    pyswisseph returns a 10-element tuple where only index 0 is meaningful
    (the Julian Day of the event). Remaining elements are always 0.0.
    """
    return (jd_event, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def _calculate_transit_leb(
    jd_start: float,
    lon: float,
    event_type: int,
    get_ra_dec,
    CALC_ITRANSIT: int,
) -> Tuple[int, Tuple[float, ...]]:
    """LEB version of _calculate_transit using a callable for RA/Dec."""

    def _get_body_hour_angle(jd: float) -> float:
        ra_hours, _ = get_ra_dec(jd)
        ra_deg = ra_hours * 15.0
        t = get_timescale().ut1_jd(jd)
        # Apparent sidereal time: the RA is apparent (equinox of date with
        # nutation), so GAST is the consistent pairing — GMST left a ~1.1 s
        # equation-of-equinoxes error in transit times.
        gast = t.gast
        lst = gast + lon / 15.0
        lst_deg = lst * 15.0
        ha = (lst_deg - ra_deg) % 360.0
        if ha > 180:
            ha -= 360
        return ha

    ha = _get_body_hour_angle(jd_start)
    if event_type == CALC_ITRANSIT:
        ha = (ha - 180) % 360
        if ha > 180:
            ha -= 360

    sidereal_day = 0.99726957
    if ha > 0:
        dt = (360.0 - ha) / 360.0 * sidereal_day
    else:
        dt = (-ha) / 360.0 * sidereal_day

    jd_guess = jd_start + dt
    for _ in range(30):
        ha = _get_body_hour_angle(jd_guess)
        if event_type == CALC_ITRANSIT:
            ha = (ha - 180) % 360
            if ha > 180:
                ha -= 360
        if abs(ha) < 1.0 / 3600.0:
            return 0, _make_tret(jd_guess)
        jd_guess -= ha / (360.0 / sidereal_day)

    return 0, _make_tret(jd_guess)


def _calculate_transit(
    jd_start: float,
    lat: float,
    lon: float,
    event_type: int,
    ts,
    earth,
    target,
    observer,
    CALC_ITRANSIT: int,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate meridian transit time."""
    # Transit occurs when body's RA = Local Sidereal Time
    # LST = GMST + longitude (in hours)
    # So transit occurs when GMST = RA - longitude/15

    def _get_body_hour_angle(jd: float) -> float:
        """Calculate body's hour angle (in degrees). 0 = on meridian."""
        t = ts.ut1_jd(jd)

        # Get body's RA (epoch of date)
        body_app = earth.at(t).observe(target).apparent()
        ra, _, _ = body_app.radec(epoch="date")
        ra_deg = ra.hours * 15.0  # Convert to degrees

        # Get Local APPARENT Sidereal Time (apparent RA pairs with GAST;
        # GMST left a ~1.1 s equation-of-equinoxes error in transits)
        gast = t.gast  # in hours
        lst = gast + lon / 15.0  # in hours
        lst_deg = lst * 15.0  # Convert to degrees

        # Hour angle = LST - RA (normalized to -180 to +180)
        ha = (lst_deg - ra_deg) % 360.0
        if ha > 180:
            ha -= 360
        return ha

    # Get initial hour angle
    ha = _get_body_hour_angle(jd_start)

    # Adjust for lower transit (HA = 180 instead of 0)
    if event_type == CALC_ITRANSIT:
        ha = (ha - 180) % 360
        if ha > 180:
            ha -= 360

    # Estimate time to next transit
    # Earth rotates 360° in ~23h56m (sidereal day)
    # So hour angle changes ~15°/hour = 360°/day (approximately)
    sidereal_day = 0.99726957  # days

    if ha > 0:
        # Body is west of meridian, wait for it to come around
        dt = (360.0 - ha) / 360.0 * sidereal_day
    else:
        # Body is east of meridian, it will cross soon
        dt = (-ha) / 360.0 * sidereal_day

    jd_guess = jd_start + dt

    # Newton-Raphson refinement
    for _ in range(30):
        ha = _get_body_hour_angle(jd_guess)

        # Adjust for lower transit
        if event_type == CALC_ITRANSIT:
            ha = (ha - 180) % 360
            if ha > 180:
                ha -= 360

        # Check convergence (< 1 arcsecond in HA = < 0.07 seconds in time)
        if abs(ha) < 1.0 / 3600.0:
            return 0, _make_tret(jd_guess)

        # Rate of change of HA: approximately 360°/sidereal day
        ha_rate = 360.0 / sidereal_day  # degrees/day

        # Newton-Raphson step
        jd_guess -= ha / ha_rate

    # If we get here, convergence failed but return best estimate
    return 0, _make_tret(jd_guess)


def _calculate_rise_set(
    jd_start: float,
    lat: float,
    lon: float,
    event_type: int,
    target_altitude: float,
    ts,
    earth,
    target,
    observer,
    _get_body_altaz,
    _get_body_ra_dec,
    CALC_RISE: int,
    CALC_SET: int,
    rsmi: int,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate rise or set time using bisection and Newton-Raphson."""

    # Get current altitude and declination
    alt, _ = _get_body_altaz(jd_start)
    _, dec = _get_body_ra_dec(jd_start)

    # Check for circumpolar objects
    # An object is circumpolar if its altitude never crosses the target altitude
    # Calculate the altitude at upper and lower transit

    # Maximum altitude = 90 - |lat - dec|
    # Minimum altitude = -(90 - |lat + dec|) = |lat + dec| - 90
    # For circumpolar: min_alt > target_altitude (never sets)
    # For never rises: max_alt < target_altitude

    # For fast-moving bodies (esp. Moon ~13°/day), declination changes
    # significantly during the search window. We estimate the declination
    # speed and add a margin to avoid false circumpolar detection.
    max_search = 2.0  # Search up to 2 days ahead
    _, dec2 = _get_body_ra_dec(jd_start + 1.0 / 24.0)
    dec_speed = abs(dec2 - dec) * 24.0  # degrees per day
    margin = dec_speed * max_search  # total possible dec change in search window

    # Calculate the altitude at upper and lower transit
    if lat >= 0:  # Northern hemisphere
        max_alt = 90.0 - abs(lat - dec)  # At upper transit
        min_alt = dec - (90.0 - lat)  # At lower transit
    else:  # Southern hemisphere
        max_alt = 90.0 - abs(lat - dec)
        min_alt = -dec - (90.0 + lat)

    # Check circumpolar conditions with margin for declination change
    is_circumpolar_above = (min_alt - margin) > target_altitude  # Never sets
    is_circumpolar_below = (max_alt + margin) < target_altitude  # Never rises

    if event_type == CALC_RISE and is_circumpolar_below:
        return -2, _make_tret()  # Never rises
    if event_type == CALC_SET and is_circumpolar_above:
        return -2, _make_tret()  # Never sets
    if event_type == CALC_RISE and is_circumpolar_above:
        return -2, _make_tret()  # Always above horizon, no rise
    if event_type == CALC_SET and is_circumpolar_below:
        return -2, _make_tret()  # Always below horizon, no set

    # Determine search step based on how close to grazing conditions.
    # For the Moon at high latitudes, rise and set can be very close together
    # (< 30 min), so 1-hour steps would miss the crossing entirely.
    # Use finer steps when altitude extremes are close to the target.
    grazing_margin = min(abs(max_alt - target_altitude), abs(min_alt - target_altitude))
    if grazing_margin < 2.0:
        # Near-grazing: use 10-minute steps to catch brief crossings
        search_step = 10.0 / 1440.0  # 10 minutes
    elif grazing_margin < 5.0:
        # Moderately close: use 20-minute steps
        search_step = 20.0 / 1440.0  # 20 minutes
    else:
        search_step = 1.0 / 24.0  # 1 hour steps

    # Find bracket where altitude crosses target
    jd = jd_start
    alt_prev, _ = _get_body_altaz(jd)

    max_search = 2.0  # Search up to 2 days ahead
    jd_cross_start = None
    jd_cross_end = None

    for i in range(int(max_search / search_step) + 1):
        jd = jd_start + i * search_step
        alt, _ = _get_body_altaz(jd)

        # Check for crossing
        if event_type == CALC_RISE:
            # Rising: altitude going from below to above target
            if alt_prev < target_altitude and alt >= target_altitude:
                jd_cross_start = jd - search_step
                jd_cross_end = jd
                break
        else:  # CALC_SET
            # Setting: altitude going from above to below target
            if alt_prev >= target_altitude and alt < target_altitude:
                jd_cross_start = jd - search_step
                jd_cross_end = jd
                break

        alt_prev = alt

    if jd_cross_start is None:
        # No crossing found in search range - might be circumpolar or rare case
        return -2, _make_tret()

    # jd_cross_end is always set when jd_cross_start is set
    assert jd_cross_end is not None

    # Bisection to narrow down the crossing
    # Use tighter tolerance for better precision
    # 0.0001° ≈ 0.36 arcsec ≈ 0.02 seconds of time
    for _ in range(50):  # More iterations for tighter convergence
        jd_mid = (jd_cross_start + jd_cross_end) / 2
        alt_mid, _ = _get_body_altaz(jd_mid)

        if (
            abs(alt_mid - target_altitude) < 0.0001
        ):  # ~0.36 arcsec (improved from 0.001°)
            return 0, _make_tret(jd_mid)

        if event_type == CALC_RISE:
            if alt_mid < target_altitude:
                jd_cross_start = jd_mid
            else:
                jd_cross_end = jd_mid
        else:  # CALC_SET
            if alt_mid >= target_altitude:
                jd_cross_start = jd_mid
            else:
                jd_cross_end = jd_mid

    return 0, _make_tret((jd_cross_start + jd_cross_end) / 2)


# Alias for reference API compatibility
rise_trans = rise_trans


def rise_trans_true_hor(
    tjdut: float,
    body: "Union[int, str]",
    rsmi: int,
    geopos: Sequence[float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    horhgt: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...]]:
    """Calculate rise, set, or transit with custom horizon altitude.

    Wrapper around :func:`_rise_trans_true_hor_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _rise_trans_true_hor_impl,
        tjdut, body, rsmi, geopos, atpress, attemp, horhgt, flags,
    )


def _rise_trans_true_hor_impl(
    tjdut: float,
    body: "Union[int, str]",
    rsmi: int,
    geopos: Sequence[float],
    atpress: float = 0.0,
    attemp: float = 0.0,
    horhgt: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Calculate rise, set, or transit time with a custom horizon altitude.

    This function is similar to rise_trans() but allows specifying a custom
    horizon altitude. This is useful for locations with mountains or buildings
    that occlude the real horizon.

    Args:
        tjdut: Julian Day (UT) to start search from
        body: Planet/body ID (int, e.g. SUN, MOON) or fixed star
            name (str, e.g. "Sirius", "Regulus")
        rsmi: Event type and calculation flags (bitmask):
            - CALC_RISE (1): Rise time (body crossing horizon going up)
            - CALC_SET (2): Set time (body crossing horizon going down)
            - CALC_MTRANSIT (4): Upper meridian transit (culmination)
            - CALC_ITRANSIT (8): Lower meridian transit (anti-culmination)
            Additional flags (OR with event type):
            - BIT_DISC_CENTER (256): Use disc center instead of upper limb
            - BIT_DISC_BOTTOM (8192): Use lower limb of disc
            - BIT_NO_REFRACTION (512): Ignore atmospheric refraction
        geopos: Geographic position as sequence [lon, lat, alt]:
            - 0: geographic longitude in degrees (eastern positive)
            - 1: geographic latitude in degrees (northern positive)
            - 2: geographic altitude in meters above sea level
        atpress: Atmospheric pressure in mbar/hPa for refraction (default 0.0)
        attemp: Atmospheric temperature in degrees Celsius (default 0.0)
        horhgt: Height of local horizon in degrees (default 0.0).
            Positive values mean the horizon is elevated (e.g., mountains),
            so rise times will be later and set times earlier.
            Negative values mean the horizon is depressed (e.g., observer
            on a mountain).
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - res: integer flag (0 = event found, -2 = circumpolar)
            - tret: tuple of 10 floats, of which tret[0] = JD of event

    Raises:
        ValueError: If invalid planet ID or star name

    Note:
        For circumpolar objects (always above or below horizon at the given
        latitude), the function returns (-2, tret) with tret[0] = 0.0.
        For transits, circumpolar objects still have valid transit times.

        Twilight flags (BIT_CIVIL_TWILIGHT, BIT_NAUTIC_TWILIGHT,
        BIT_ASTRO_TWILIGHT) are NOT supported in this function since
        the horhgt parameter already specifies the target altitude.

    Algorithm:
        1. For transits: Find when body crosses the local meridian
           (Local Sidereal Time = body's Right Ascension)
        2. For rise/set: Find when body's altitude crosses the custom horizon
           accounting for refraction and disc size
        3. Uses Newton-Raphson iteration for precise timing

    Precision:
        Rise/set times accurate to ~1 minute for Sun/Moon (due to refraction
        uncertainty), better than 1 second for transit times.

    Example:
        >>> from libephemeris import julday, rise_trans_true_hor, SUN, CALC_RISE
        >>> jd = julday(2024, 6, 21, 0)
        >>> # Find sunrise with mountains at 5 deg above horizon at Rome
        >>> res, tret = rise_trans_true_hor(jd, SUN, CALC_RISE,
        ...                                [12.5, 41.9, 0], horhgt=5.0)
        >>> print(f"Sunrise at JD {tret[0]:.5f}")

    References:
        - Reference API: rise_trans_true_hor()
        - Meeus "Astronomical Algorithms" Ch. 15 (Rise, Set, Transit)
    """
    from skyfield.api import wgs84

    from .constants import (
        CALC_RISE,
        CALC_SET,
        CALC_MTRANSIT,
        CALC_ITRANSIT,
        BIT_DISC_CENTER,
        BIT_DISC_BOTTOM,
        BIT_NO_REFRACTION,
    )
    from .planets import _PLANET_MAP, get_planet_target
    from .state import get_planets, get_timescale

    # Unpack geopos: [longitude, latitude, altitude]
    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2]) if len(geopos) > 2 else 0.0

    # Map parameter names for internal use
    jd_start = tjdut
    pressure = atpress
    temperature = attemp
    horizon_altitude = horhgt
    is_fixed_star = isinstance(body, str)

    # Extract event type from rsmi (lower bits)
    event_type = rsmi & 0x0F  # First 4 bits for event type

    # Validate event type
    if event_type not in (
        CALC_RISE,
        CALC_SET,
        CALC_MTRANSIT,
        CALC_ITRANSIT,
    ):
        raise ValueError(
            "Invalid event type in rsmi: %d. Use CALC_RISE, CALC_SET, CALC_MTRANSIT, or CALC_ITRANSIT"
            % rsmi
        )

    # reader is provided by the caller (None forces Skyfield path)
    _reader_rt = reader
    _use_leb_rt = _reader_rt is not None

    if _use_leb_rt:
        try:
            from . import fast_calc as _fc_rt
            from .time_utils import deltat as _sd_rt
            if not is_fixed_star:
                _fc_rt.fast_calc_ut(_reader_rt, jd_start, body, FLG_SPEED)
        except (KeyError, ValueError):
            _use_leb_rt = False

    if _use_leb_rt:
        from .constants import FLG_EQUATORIAL

        geopos_leb = (lon, lat, altitude)
        if is_fixed_star:
            planet = -1
        else:
            planet = body
            if planet not in _PLANET_MAP:
                raise ValueError("illegal planet number %d." % planet)
        ts = get_timescale()
        earth = None
        target = None
        observer = None
    else:
        eph = get_planets()
        ts = get_timescale()
        earth = eph["earth"]
        if is_fixed_star:
            from skyfield.api import Star as SkyfieldStar
            from .fixed_stars import FIXED_STARS, _resolve_star_id

            star_id, err, _ = _resolve_star_id(body)
            if err is not None:
                raise ValueError(err)
            star_entry = FIXED_STARS[star_id]
            t_years = (jd_start - 2451545.0) / 365.25
            ra_deg = star_entry.ra_j2000 + (star_entry.pm_ra * t_years) / 3600.0
            dec_deg = star_entry.dec_j2000 + (star_entry.pm_dec * t_years) / 3600.0
            target = SkyfieldStar(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
            planet = -1
        else:
            planet = body
            if planet not in _PLANET_MAP:
                raise ValueError("illegal planet number %d." % planet)
            target_name = _PLANET_MAP[planet]
            target = get_planet_target(eph, target_name)
        observer = wgs84.latlon(lat, lon, altitude)

    # Handle atpress=0: estimate pressure from altitude (pyswisseph convention).
    # When atpress is 0, Swiss Ephemeris estimates pressure using the
    # barometric formula. We replicate this behavior.
    if pressure == 0.0:
        # Barometric formula: P = 1013.25 * (1 - 2.25577e-5 * h)^5.25588
        pressure = 1013.25 * (1.0 - 2.25577e-5 * altitude) ** 5.25588
        if temperature == 0.0:
            # Standard temperature lapse rate: T = 15 - 0.0065 * h
            temperature = 15.0 - 0.0065 * altitude

    # Clamp negative horizon altitudes to 0.0 for API compatibility.
    # Confirmed by independent testing: negative horizon_altitude values
    # produce identical results to horizon_altitude=0.0 in the reference
    # implementation.
    effective_horizon = max(0.0, horizon_altitude)

    # Determine refraction using Bennett (1982) formula
    # R = 1/tan(h + 7.31/(h+4.4)) arcminutes, where h = apparent altitude in degrees
    # For true_hor, refraction depends on the custom horizon altitude, not a flat value.
    # Apply pressure/temperature correction: factor = (P/1010) * (283/(273+T))
    if rsmi & BIT_NO_REFRACTION:
        refraction = 0.0
    else:
        # Bennett (1982) refraction at the horizon altitude
        # Use the effective horizon altitude for the refraction estimate
        h = effective_horizon
        r_arcmin = 1.0 / math.tan(math.radians(h + 7.31 / (h + 4.4)))
        # Apply pressure/temperature correction (standard: 1010 mbar, 10°C = 283K)
        pressure_factor = pressure / 1010.0
        temperature_factor = 283.0 / (273.0 + temperature)
        refraction = max(0.0, r_arcmin / 60.0 * pressure_factor * temperature_factor)

    # Use the effective (clamped) horizon altitude
    horizon_alt = effective_horizon

    # Account for disc semi-diameter
    # Sun: ~16 arcmin, Moon: ~16 arcmin (varies)
    if rsmi & BIT_DISC_CENTER:
        disc_correction = 0.0
    elif rsmi & BIT_DISC_BOTTOM:
        # Lower limb: add semi-diameter (rises later, sets earlier)
        if planet == SUN:
            disc_correction = 16.0 / 60.0  # degrees
        elif planet == MOON:
            disc_correction = 16.0 / 60.0
        else:
            disc_correction = 0.0
    else:
        # Default: upper limb (subtract semi-diameter)
        if planet == SUN:
            disc_correction = -16.0 / 60.0  # degrees
        elif planet == MOON:
            disc_correction = -16.0 / 60.0
        else:
            disc_correction = 0.0

    # Effective horizon altitude (negative means below geometric horizon)
    target_altitude = horizon_alt - refraction + disc_correction

    if _use_leb_rt:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        def _get_body_altaz(jd: float) -> Tuple[float, float]:
            jd_tt = jd + deltat(jd)
            if is_fixed_star:
                from .fixed_stars import fixstar_ut
                star_pos, _, _ = fixstar_ut(body, jd, FLG_SPEED)
                az, alt_true, _ = azalt(jd, ECL2HOR, geopos_leb, 0, 0, star_pos[:3])
            else:
                pos = _topo_ecliptic(_reader_rt, jd_tt, jd, planet, geopos_leb, FLG_SPEED)
                az, alt_true, _ = azalt(jd, ECL2HOR, geopos_leb, 0, 0, pos[:3])
            return alt_true, az

        def _get_body_ra_dec(jd: float) -> Tuple[float, float]:
            if is_fixed_star:
                from .fixed_stars import fixstar_ut
                pos, _, _ = fixstar_ut(body, jd, FLG_EQUATORIAL | FLG_SPEED)
                return pos[0] / 15.0, pos[1]
            else:
                eq, _ = calc_ut(jd, planet, FLG_EQUATORIAL | FLG_SPEED)
                return eq[0] / 15.0, eq[1]
    else:
        def _get_body_altaz(jd: float) -> Tuple[float, float]:
            t = ts.ut1_jd(jd)
            observer_at = earth + observer
            body_app = observer_at.at(t).observe(target).apparent()
            alt, az, _ = body_app.altaz()
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_body_ra_dec(jd: float) -> Tuple[float, float]:
            t = ts.ut1_jd(jd)
            body_app = earth.at(t).observe(target).apparent()
            ra, dec, _ = body_app.radec(epoch="date")
            return ra.hours, dec.degrees

    if event_type in (CALC_MTRANSIT, CALC_ITRANSIT):
        if _use_leb_rt:
            return _calculate_transit_leb(
                jd_start, lon, event_type, _get_body_ra_dec, CALC_ITRANSIT,
            )
        return _calculate_transit(
            jd_start, lat, lon, event_type, ts, earth, target, observer,
            CALC_ITRANSIT,
        )

    return _calculate_rise_set(
        jd_start, lat, lon, event_type, target_altitude, ts, earth, target,
        observer, _get_body_altaz, _get_body_ra_dec, CALC_RISE, CALC_SET,
        rsmi,
    )


# Alias for reference API compatibility
rise_trans_true_hor = rise_trans_true_hor


# =============================================================================
# HELIACAL RISING AND SETTING CALCULATIONS
# =============================================================================


def heliacal_ut(
    jd_start: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SUN,
    event_type: int = 1,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, int]:
    """
    Calculate heliacal rising or setting time for a celestial body.

    Heliacal events are the first/last visibility of a celestial body
    at dawn or dusk. These were fundamental for ancient calendars:
    - Heliacal rising: First morning visibility after a period of invisibility
    - Heliacal setting: Last evening visibility before becoming invisible

    Args:
        jd_start: Julian Day (UT) to start search from
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        humidity: Relative humidity 0.0-1.0 for atmospheric extinction (default 0.5)
        body: Planet/body ID (MERCURY, VENUS, MARS, JUPITER, SATURN,
              or fixed stars). Note: SUN and MOON are not valid for heliacal events.
        event_type: Type of heliacal event:
            - HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - EVENING_FIRST (3): First evening visibility (after superior conjunction)
            - MORNING_LAST (4): Last morning visibility (before superior conjunction)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - jd_event: Julian Day (UT) of the heliacal event, or 0.0 if not found
            - retflag: Return flag (event_type on success, negative on error)

    Raises:
        ValueError: If invalid body ID or event_type

    Algorithm:
        The algorithm searches for the moment when:
        1. The body is at a specific altitude above the horizon (arcus visionis)
        2. The Sun is at twilight position (typically -6° to -12° below horizon)
        3. The body's apparent magnitude is brighter than the sky's limiting magnitude

        For heliacal rising (morning first):
        - Search forward for when the body first becomes visible at dawn
        - Body must be above horizon while Sun is still below
        - Sky must be dark enough for the body to be seen

        For heliacal setting (evening last):
        - Search forward for when the body is last visible at dusk
        - Body must be above horizon while Sun is setting
        - Sky brightness must not overwhelm the body's light

    Historical Note:
        Heliacal risings were crucial for ancient calendars. The heliacal
        rising of Sirius marked the Egyptian new year and predicted the
        Nile flood. Babylonians used heliacal events to track planetary
        positions without modern instruments.

    Example:
        >>> from libephemeris import julday, heliacal_ut, VENUS, HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Find next heliacal rising of Venus from Rome
        >>> jd_event, flag = heliacal_ut(jd, 41.9, 12.5, body=VENUS,
        ...                              event_type=HELIACAL_RISING)
        >>> print(f"Heliacal rising at JD {jd_event:.5f}")

    References:
        - Reference API: heliacal_ut()
        - Schoch "Planets in Mesopotamian Astral Science"
        - Ptolemy's criteria for heliacal visibility
    """
    from .constants import (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    )
    from .planets import _PLANET_MAP, get_planet_target, pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST."
        )

    # Sun and Moon are not valid for heliacal events
    if body == SUN:
        raise ValueError("SUN is not valid for heliacal calculations")
    if body == MOON:
        raise ValueError("MOON is not valid for heliacal calculations")

    # Validate body
    if body not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {body}.")

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Get celestial bodies
    target_name = _PLANET_MAP[body]
    target = get_planet_target(eph, target_name)
    sun = eph["sun"]
    earth = eph["earth"]

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Note: Arcus visionis (minimum altitude for visibility) and
    # sun altitude thresholds are used in the visibility check functions.
    # Typical values: arcus visionis ~7° for bright objects, sun altitude ~-8°.

    def _get_altitudes(jd: float) -> Tuple[float, float, float]:
        """Get Sun altitude, body altitude, and body azimuth at given JD."""
        t = ts.ut1_jd(jd)
        observer_at = earth + observer

        # Sun position
        sun_app = observer_at.at(t).observe(sun).apparent()
        sun_alt, _, _ = sun_app.altaz()

        # Body position
        body_app = observer_at.at(t).observe(target).apparent()
        body_alt, body_az, _ = body_app.altaz()

        # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
        return sun_alt.degrees, body_alt.degrees, (body_az.degrees + 180.0) % 360.0

    def _get_elongation(jd: float) -> float:
        """Get the elongation of body from Sun in degrees."""
        try:
            pheno = pheno_ut(jd, body, flags)
            return pheno[2]  # Elongation
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # Fallback: calculate elongation manually
            t = ts.ut1_jd(jd)
            sun_app = earth.at(t).observe(sun).apparent()
            body_app = earth.at(t).observe(target).apparent()
            return body_app.separation_from(sun_app).degrees

    def _get_body_magnitude(jd: float) -> float:
        """Get the visual magnitude of the body."""
        try:
            pheno = pheno_ut(jd, body, flags)
            return pheno[4]  # Visual magnitude
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0  # Default to bright magnitude

    def _calculate_limiting_magnitude(sun_alt: float, body_alt: float) -> float:
        """
        Calculate the limiting magnitude based on sky brightness.

        The limiting magnitude depends on:
        - Sun's altitude below horizon (darker = higher limit)
        - Body's altitude (extinction increases near horizon)
        - Atmospheric conditions (humidity, pressure)

        This is a simplified model based on Schaefer (1990).
        """
        # Sky brightness model (simplified)
        # When Sun is below -18°, sky is fully dark (limiting mag ~6.5)
        # When Sun is at -6° (civil twilight), limiting mag is ~3
        # When Sun is at 0° (horizon), limiting mag is ~-2

        if sun_alt >= 0:
            return -2.0  # Daylight - essentially nothing visible
        elif sun_alt >= -6:
            # Civil twilight
            return -2.0 + (sun_alt + 6) * (-3.0 - (-2.0)) / (-6)
        elif sun_alt >= -12:
            # Nautical twilight
            return 3.0 + (sun_alt + 12) * (5.0 - 3.0) / 6
        elif sun_alt >= -18:
            # Astronomical twilight
            return 5.0 + (sun_alt + 18) * (6.5 - 5.0) / 6
        else:
            # Full darkness
            return 6.5

        # Atmospheric extinction correction (simplified)
        # Add extinction based on altitude (airmass)
        # At low altitudes, extinction can be 0.5-1.0 magnitudes

    def _is_body_visible(jd: float) -> Tuple[bool, float, float, float]:
        """
        Check if body is visible at given time.

        Returns: (is_visible, sun_alt, body_alt, elongation)
        """
        sun_alt, body_alt, body_az = _get_altitudes(jd)
        elongation = _get_elongation(jd)

        # Body must be above minimum altitude
        if body_alt < 0:
            return False, sun_alt, body_alt, elongation

        # Sun must be below horizon
        if sun_alt > 0:
            return False, sun_alt, body_alt, elongation

        # Check limiting magnitude vs body magnitude
        limiting_mag = _calculate_limiting_magnitude(sun_alt, body_alt)
        body_mag = _get_body_magnitude(jd)

        # Body is visible if its magnitude is brighter (lower) than limiting magnitude
        # Also require sufficient elongation from Sun
        min_elongation = 10.0  # Minimum elongation for visibility (degrees)

        is_visible = (body_mag <= limiting_mag) and (elongation >= min_elongation)

        return is_visible, sun_alt, body_alt, elongation

    def _find_twilight_time(jd: float, sun_target_alt: float, rising: bool) -> float:
        """
        Find when Sun crosses target altitude (morning or evening).

        Args:
            jd: Starting JD
            sun_target_alt: Target Sun altitude (negative for below horizon)
            rising: True for morning (Sun rising), False for evening (Sun setting)
        """
        # Search within one day
        for _ in range(50):
            sun_alt, _, _ = _get_altitudes(jd)

            # Check if we're at the right phase of day
            if rising:
                # Morning: looking for Sun rising through target altitude
                if abs(sun_alt - sun_target_alt) < 0.01:
                    return jd
            else:
                # Evening: looking for Sun setting through target altitude
                if abs(sun_alt - sun_target_alt) < 0.01:
                    return jd

            # Adjust time based on Sun's position
            # Sun moves ~15°/hour = 0.25°/minute = 360°/day
            sun_rate = 360.0 / 1.0  # degrees per day
            diff = sun_target_alt - sun_alt

            # Crude estimate of time adjustment
            dt = diff / sun_rate

            # Limit step size
            dt = max(-0.1, min(0.1, dt))

            if abs(dt) < 1e-6:
                return jd

            jd += dt

        return jd

    def _search_heliacal_rising(jd_start: float) -> float:
        """
        Search for heliacal rising (morning first visibility).

        The body becomes visible in the morning before sunrise after
        a period of being hidden in the Sun's glare.
        """
        max_days = 400  # Search up to more than a year

        # Use coarse search first (check once per day at approximate twilight)
        for day in range(max_days):
            jd_day = jd_start + day

            # Check at approximately 5 AM local time (rough morning twilight estimate)
            # This is a coarse check - we'll refine if we find visibility
            for hour in [4, 5, 6]:  # Check a few morning hours
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                # Look for conditions: Sun below horizon but not too deep, body above
                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        # Found a potential heliacal rising
                        # Refine the time using binary search
                        return _refine_heliacal_time(jd_check, is_morning=True)

        return 0.0  # Not found

    def _search_heliacal_setting(jd_start: float) -> float:
        """
        Search for heliacal setting (evening last visibility).

        The body is last visible in the evening after sunset before
        becoming hidden in the Sun's glare.
        """
        max_days = 400

        # First, find when the body is visible in the evening
        first_visible_jd = 0.0

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few evening hours (coarse search)
            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                # Look for conditions: Sun below horizon, body above
                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        first_visible_jd = jd_check
                        break

            if first_visible_jd > 0:
                break

        if first_visible_jd == 0:
            return 0.0

        # Now find the LAST day of visibility
        last_visible_jd = first_visible_jd
        start_day = int(first_visible_jd - jd_start)

        for day in range(start_day, start_day + 120):  # Check up to 120 days ahead
            jd_day = jd_start + day
            found_visible = False

            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0
                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3 and body_alt > 1:
                    visible, _, _, _ = _is_body_visible(jd_check)
                    if visible:
                        last_visible_jd = jd_check
                        found_visible = True
                        break

            if not found_visible and day > start_day:
                # No visibility found today - return the last visible day
                return _refine_heliacal_time(last_visible_jd, is_morning=False)

        # If still visible after 120 days, return the last checked
        return _refine_heliacal_time(last_visible_jd, is_morning=False)

    def _search_evening_first(jd_start: float) -> float:
        """
        Search for evening first visibility (after superior conjunction).

        The body appears in the evening sky for the first time after
        passing behind the Sun (superior conjunction for inferior planets,
        or conjunction for superior planets).
        """
        max_days = 400
        was_invisible = False

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few evening hours (coarse search)
            for hour in [18, 19, 20]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)
                    if not visible and body_alt < 5:
                        was_invisible = True
                    elif visible and was_invisible:
                        # First evening visibility after being invisible
                        return _refine_heliacal_time(jd_check, is_morning=False)

        return 0.0

    def _search_morning_last(jd_start: float) -> float:
        """
        Search for morning last visibility (before superior conjunction).

        The body is last visible in the morning sky before passing
        behind the Sun.
        """
        max_days = 400
        last_visible_jd = 0.0
        found_visible = False

        for day in range(max_days):
            jd_day = jd_start + day

            # Check at a few morning hours (coarse search)
            for hour in [4, 5, 6]:
                jd_check = jd_day + hour / 24.0

                sun_alt, body_alt, _ = _get_altitudes(jd_check)

                if -15 < sun_alt < -3:
                    visible, s_alt, b_alt, elong = _is_body_visible(jd_check)

                    if visible:
                        last_visible_jd = jd_check
                        found_visible = True
                    elif found_visible and not visible:
                        # Found the transition to invisibility
                        return _refine_heliacal_time(last_visible_jd, is_morning=True)

        return last_visible_jd if last_visible_jd > 0 else 0.0

    def _refine_heliacal_time(jd_approx: float, is_morning: bool) -> float:
        """
        Refine the heliacal event time using binary search.

        Find the exact moment when the body becomes just visible/invisible.
        """
        # Use binary search to find the transition point
        jd_low = jd_approx - 0.1  # ~2.4 hours before
        jd_high = jd_approx + 0.1  # ~2.4 hours after

        for _ in range(30):  # ~30 iterations gives very high precision
            jd_mid = (jd_low + jd_high) / 2

            visible, sun_alt, body_alt, elong = _is_body_visible(jd_mid)

            # For heliacal rising: looking for first visibility
            # For heliacal setting: looking for last visibility
            if is_morning:
                # Morning: visibility increases as time progresses (Sun rises)
                # But we want first visibility, so search backwards
                if visible:
                    jd_high = jd_mid  # Look earlier
                else:
                    jd_low = jd_mid  # Look later
            else:
                # Evening: visibility decreases as time progresses (sky darkens)
                # For last visibility, we want the last moment still visible
                if visible:
                    jd_low = jd_mid  # Look later for last visibility
                else:
                    jd_high = jd_mid  # Look earlier

            if jd_high - jd_low < 1e-6:  # ~0.1 second precision
                break

        return (jd_low + jd_high) / 2

    # Main search logic based on event type
    if event_type == HELIACAL_RISING:
        jd_event = _search_heliacal_rising(jd_start)
    elif event_type == HELIACAL_SETTING:
        jd_event = _search_heliacal_setting(jd_start)
    elif event_type == EVENING_FIRST:
        jd_event = _search_evening_first(jd_start)
    elif event_type == MORNING_LAST:
        jd_event = _search_morning_last(jd_start)
    else:
        jd_event = 0.0

    if jd_event > 0:
        return jd_event, event_type
    else:
        return 0.0, -1  # Not found


# Note: heliacal_ut is defined in heliacal.py with the pyswisseph-compatible
# API signature (returns 3-tuple). The internal heliacal_ut here returns
# (jd_event, retflag) and is used by heliacal_ut internally.


def heliacal_pheno_ut(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    pressure: float = 1013.25,
    temperature: float = 15.0,
    humidity: float = 0.5,
    body: int = SUN,
    event_type: int = 1,
    flags: int = FLG_SWIEPH,
) -> Tuple[Tuple[float, ...], int]:
    """
    Provides data relevant for the calculation of heliacal risings and settings.

    This function calculates detailed phenomena associated with heliacal events,
    including altitudes, azimuths, arcus visionis, magnitude, visibility times,
    and other parameters used in heliacal visibility calculations.

    Args:
        jd: Julian Day (UT) for the calculation
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        pressure: Atmospheric pressure in mbar/hPa for refraction (default 1013.25)
        temperature: Temperature in Celsius for refraction (default 15)
        humidity: Relative humidity 0.0-1.0 for atmospheric extinction (default 0.5)
        body: Planet/body ID (MERCURY, VENUS, MARS, JUPITER, SATURN, etc.)
        event_type: Type of heliacal event:
            - HELIACAL_RISING (1): Morning first visibility (heliacal rising)
            - HELIACAL_SETTING (2): Evening last visibility (heliacal setting)
            - EVENING_FIRST (3): First evening visibility (after superior conjunction)
            - MORNING_LAST (4): Last morning visibility (before superior conjunction)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - dret: Tuple of 50 floats with heliacal phenomena data:
                - 0: AltO [deg] topocentric altitude of object (unrefracted)
                - 1: AppAltO [deg] apparent altitude of object (refracted)
                - 2: GeoAltO [deg] geocentric altitude of object
                - 3: AziO [deg] azimuth of object
                - 4: AltS [deg] topocentric altitude of Sun
                - 5: AziS [deg] azimuth of Sun
                - 6: TAVact [deg] actual topocentric arcus visionis
                - 7: ARCVact [deg] actual (geocentric) arcus visionis
                - 8: DAZact [deg] actual difference between object's and sun's azimuth
                - 9: ARCLact [deg] actual longitude difference between object and sun
                - 10: kact [-] extinction coefficient
                - 11: minTAV [deg] smallest topocentric arcus visionis
                - 12: TfirstVR [JDN] first time object is visible, according to VR
                - 13: TbVR [JDN] optimum time the object is visible, according to VR
                - 14: TlastVR [JDN] last time object is visible, according to VR
                - 15: TbYallop [JDN] best time the object is visible, according to Yallop
                - 16: WMoon [deg] crescent width of Moon
                - 17: qYal [-] q-test value of Yallop
                - 18: qCrit [-] q-test criterion of Yallop
                - 19: ParO [deg] parallax of object
                - 20: Magn [-] magnitude of object
                - 21: RiseO [JDN] rise/set time of object
                - 22: RiseS [JDN] rise/set time of Sun
                - 23: Lag [JDN] rise/set time of object minus rise/set time of Sun
                - 24: TvisVR [JDN] visibility duration
                - 25: LMoon [deg] crescent length of Moon
                - 26-49: Reserved for future use
            - retflag: Return flag (flags on success, negative on error)

    Raises:
        ValueError: If invalid body ID or event_type

    Example:
        >>> from libephemeris import julday, heliacal_pheno_ut, VENUS, HELIACAL_RISING
        >>> jd = julday(2024, 1, 1, 0)
        >>> # Get heliacal phenomena for Venus at Rome
        >>> dret, flag = heliacal_pheno_ut(jd, 41.9, 12.5, body=VENUS,
        ...                                 event_type=HELIACAL_RISING)
        >>> print(f"Object altitude: {dret[0]:.2f}°, Sun altitude: {dret[4]:.2f}°")

    References:
        - Reference API: heliacal_pheno_ut()
        - Schoch "Planets in Mesopotamian Astral Science"
    """
    import math

    from .constants import (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    )
    from .planets import _PLANET_MAP, get_planet_target, pheno_ut
    from .state import get_planets, get_timescale
    from skyfield.api import wgs84

    # Validate event type
    if event_type not in (
        HELIACAL_RISING,
        HELIACAL_SETTING,
        EVENING_FIRST,
        MORNING_LAST,
    ):
        raise ValueError(
            f"Invalid event_type: {event_type}. Use HELIACAL_RISING, "
            "HELIACAL_SETTING, EVENING_FIRST, or MORNING_LAST."
        )

    # Validate body
    if body not in _PLANET_MAP:
        raise ValueError(f"illegal planet number {body}.")

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Get celestial bodies
    target_name = _PLANET_MAP[body]
    target = get_planet_target(eph, target_name)
    sun = eph["sun"]
    earth = eph["earth"]

    # Create observer location
    observer = wgs84.latlon(lat, lon, altitude)

    # Initialize result array with 50 zeros
    dret = [0.0] * 50

    # Calculate positions at the given time
    t = ts.ut1_jd(jd)
    observer_at = earth + observer

    # Calculate Sun position
    sun_app = observer_at.at(t).observe(sun).apparent()
    sun_alt_topo, sun_az, _ = sun_app.altaz()
    sun_alt_deg = sun_alt_topo.degrees
    # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
    sun_az_deg = (sun_az.degrees + 180.0) % 360.0

    # Calculate body position
    body_app = observer_at.at(t).observe(target).apparent()
    body_alt_topo, body_az, body_dist = body_app.altaz()
    body_alt_deg = body_alt_topo.degrees
    # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
    body_az_deg = (body_az.degrees + 180.0) % 360.0

    # Get geocentric altitude (without refraction or topocentric correction)
    body_geo = earth.at(t).observe(target).apparent()
    body_geo_ra, body_geo_dec, body_geo_dist = body_geo.radec()

    # Calculate geocentric altitude using hour angle
    # First get the local sidereal time
    gast = t.gast
    lst = gast + lon / 15.0  # Local sidereal time in hours
    ra_hours = body_geo_ra.hours
    ha = (lst - ra_hours) * 15.0  # Hour angle in degrees

    # Geocentric altitude calculation
    dec_rad = math.radians(body_geo_dec.degrees)
    lat_rad = math.radians(lat)
    ha_rad = math.radians(ha)

    sin_alt = math.sin(lat_rad) * math.sin(dec_rad) + math.cos(lat_rad) * math.cos(
        dec_rad
    ) * math.cos(ha_rad)
    sin_alt = max(-1.0, min(1.0, sin_alt))
    geo_alt_deg = math.degrees(math.asin(sin_alt))

    # Calculate atmospheric refraction
    # Use simplified formula: R = 1.02 / tan(h + 10.3/(h + 5.11)) in arcminutes
    if body_alt_deg > -1:
        refraction = 1.02 / math.tan(
            math.radians(body_alt_deg + 10.3 / (body_alt_deg + 5.11))
        )
        refraction /= 60.0  # Convert arcminutes to degrees
    else:
        refraction = 0.5  # Near horizon approximation

    # Apparent altitude (with refraction)
    app_alt_deg = body_alt_deg + refraction

    # Calculate arcus visionis (altitude difference between body and Sun)
    # Topocentric arcus visionis
    tav_act = body_alt_deg - sun_alt_deg

    # Geocentric arcus visionis
    arcv_act = geo_alt_deg - sun_alt_deg

    # Azimuth difference
    daz_act = body_az_deg - sun_az_deg
    # Normalize to -180 to +180
    while daz_act > 180:
        daz_act -= 360
    while daz_act < -180:
        daz_act += 360

    # Get elongation (longitude difference) from Sun using pheno
    try:
        pheno = pheno_ut(jd, body, flags)
        elongation = pheno[2]  # Elongation
        magnitude = pheno[4]  # Visual magnitude
        phase_angle = pheno[1]  # Phase angle
    except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
        _reraise_if_leb_range_error(_exc)
        # Calculate elongation manually
        sun_geo = earth.at(t).observe(sun).apparent()
        elongation = body_geo.separation_from(sun_geo).degrees
        magnitude = 0.0
        phase_angle = 0.0

    # Calculate extinction coefficient (simplified Schaefer model)
    # k = k_rayleigh + k_aerosol + k_ozone
    # Simplified: k increases with humidity and decreases with altitude
    # Typical values: 0.2 (excellent) to 0.5 (poor)
    pressure_factor = pressure / 1013.25
    humidity_factor = 1.0 + 0.5 * humidity
    altitude_factor = math.exp(-altitude / 8500.0)  # Scale height ~8.5 km
    k_act = 0.25 * pressure_factor * humidity_factor * altitude_factor

    # Minimum topocentric arcus visionis (simplified)
    # Depends on magnitude and atmospheric conditions
    # Brighter objects need less arcus visionis
    if magnitude < 0:
        min_tav = 5.0 + 0.5 * magnitude  # Bright objects
    elif magnitude < 2:
        min_tav = 6.0 + magnitude * 0.5
    else:
        min_tav = 7.0 + magnitude * 0.7

    # Parallax of object (in degrees)
    # Parallax = arcsin(Earth_radius / distance)
    # Earth radius ~ 6371 km, distance in AU (1 AU ~ 149,597,870.7 km)
    earth_radius_au = 6371.0 / 149597870.7
    if body_geo_dist.au > 0:
        parallax = math.degrees(math.asin(earth_radius_au / body_geo_dist.au))
    else:
        parallax = 0.0

    # Calculate rise/set times for object and Sun
    # For morning events, we look for rise times; for evening, set times
    is_morning = event_type in (HELIACAL_RISING, MORNING_LAST)

    # Use a simple estimate for rise/set time based on altitude
    # Rise/set occurs when altitude crosses 0 (corrected for refraction)
    # Time rate: approximately 1 degree of altitude per 4 minutes (at mid-latitudes)
    rise_set_correction = -0.833  # Standard refraction + semidiameter for Sun

    # Estimate object rise/set time
    if body_alt_deg != 0:
        # Rough estimate: time to rise/set based on altitude rate
        # Altitude rate ~ 15° cos(lat) per hour at the horizon
        alt_rate = 15.0 * math.cos(lat_rad)  # degrees per hour
        if alt_rate > 0:
            if is_morning:
                # Object rising: how long until it crosses horizon
                time_to_horizon = (
                    body_alt_deg - rise_set_correction
                ) / alt_rate  # hours
                rise_o = jd - time_to_horizon / 24.0
            else:
                # Object setting
                time_to_horizon = (
                    body_alt_deg - rise_set_correction
                ) / alt_rate  # hours
                rise_o = jd + time_to_horizon / 24.0
        else:
            rise_o = jd
    else:
        rise_o = jd

    # Estimate Sun rise/set time
    if sun_alt_deg != 0:
        alt_rate = 15.0 * math.cos(lat_rad)
        if alt_rate > 0:
            if is_morning:
                time_to_horizon = (sun_alt_deg - rise_set_correction) / alt_rate
                rise_s = jd - time_to_horizon / 24.0
            else:
                time_to_horizon = (sun_alt_deg - rise_set_correction) / alt_rate
                rise_s = jd + time_to_horizon / 24.0
        else:
            rise_s = jd
    else:
        rise_s = jd

    # Lag time (object rise - Sun rise)
    lag = rise_o - rise_s

    # Visibility duration estimate (simplified)
    # Based on how long the object is above horizon while Sun is below
    if sun_alt_deg < -6 and body_alt_deg > 0:
        # Object visible during civil twilight or darker
        # Estimate visibility window
        tvis_vr = abs(sun_alt_deg + 6) / 15.0 / 24.0  # days
    else:
        tvis_vr = 0.0

    # For Moon-specific calculations
    w_moon = 0.0  # Crescent width
    l_moon = 0.0  # Crescent length
    illumination = 0.0

    if body == MOON:
        # Calculate Moon phase and crescent geometry
        try:
            moon_pheno = pheno_ut(jd, MOON, flags)
            phase = moon_pheno[0]  # Phase 0-1
            illumination = phase * 100.0  # Percentage

            # Crescent width approximation (Danjon's formula)
            # W = 15 * (1 - cos(phase_angle/2)) arcminutes (approximate)
            if len(moon_pheno) > 1:
                pa_rad = math.radians(moon_pheno[1])  # Phase angle in radians
                w_moon = 15.0 * (1 - math.cos(pa_rad / 2)) / 60.0  # In degrees

            # Crescent length (semicircle approximation)
            # L ≈ pi * D / 2 where D is diameter
            moon_diameter = 0.5  # Approximately 0.5 degrees
            l_moon = math.pi * moon_diameter / 2
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            pass

    # Yallop q-test (for lunar crescent visibility)
    # q = (ARCV - (11.8371 - 6.3226*W + 0.7319*W^2 - 0.1018*W^3)) / 10
    # Simplified version
    if body == MOON:
        w = w_moon * 60.0  # Convert to arcminutes for formula
        q_criterion = 11.8371 - 6.3226 * w + 0.7319 * w**2 - 0.1018 * w**3
        q_yallop = (arcv_act - q_criterion) / 10.0
        q_crit = q_criterion
    else:
        q_yallop = 0.0
        q_crit = 0.0

    # First, best, and last visibility times (VR = visibility rule)
    # Simplified: based on Sun altitude thresholds
    # First visibility: Sun at about -10°
    # Best visibility: Sun at about -8°
    # Last visibility: Sun at about -6°
    sun_rate = 15.0 * math.cos(lat_rad)  # degrees per hour
    if sun_rate > 0:
        if is_morning:
            t_first_vr = jd + (sun_alt_deg + 10) / sun_rate / 24.0
            t_best_vr = jd + (sun_alt_deg + 8) / sun_rate / 24.0
            t_last_vr = jd + (sun_alt_deg + 6) / sun_rate / 24.0
        else:
            t_first_vr = jd + (-6 - sun_alt_deg) / sun_rate / 24.0
            t_best_vr = jd + (-8 - sun_alt_deg) / sun_rate / 24.0
            t_last_vr = jd + (-10 - sun_alt_deg) / sun_rate / 24.0
    else:
        t_first_vr = jd
        t_best_vr = jd
        t_last_vr = jd

    # Best time according to Yallop (for Moon)
    t_b_yallop = t_best_vr  # Use same as best VR for simplicity

    # Fill in the result array
    dret[0] = body_alt_deg  # AltO - topocentric altitude (unrefracted)
    dret[1] = app_alt_deg  # AppAltO - apparent altitude (refracted)
    dret[2] = geo_alt_deg  # GeoAltO - geocentric altitude
    dret[3] = body_az_deg  # AziO - azimuth of object
    dret[4] = sun_alt_deg  # AltS - topocentric altitude of Sun
    dret[5] = sun_az_deg  # AziS - azimuth of Sun
    dret[6] = tav_act  # TAVact - topocentric arcus visionis
    dret[7] = arcv_act  # ARCVact - geocentric arcus visionis
    dret[8] = daz_act  # DAZact - azimuth difference
    dret[9] = elongation  # ARCLact - elongation from Sun
    dret[10] = k_act  # kact - extinction coefficient
    dret[11] = min_tav  # minTAV - minimum topocentric arcus visionis
    dret[12] = t_first_vr  # TfirstVR - first visibility time
    dret[13] = t_best_vr  # TbVR - best visibility time
    dret[14] = t_last_vr  # TlastVR - last visibility time
    dret[15] = t_b_yallop  # TbYallop - best time according to Yallop
    dret[16] = w_moon  # WMoon - crescent width
    dret[17] = q_yallop  # qYal - Yallop q-test value
    dret[18] = q_crit  # qCrit - Yallop criterion
    dret[19] = parallax  # ParO - parallax of object
    dret[20] = magnitude  # Magn - magnitude
    dret[21] = rise_o  # RiseO - rise/set time of object
    dret[22] = rise_s  # RiseS - rise/set time of Sun
    dret[23] = lag  # Lag - time difference
    dret[24] = tvis_vr  # TvisVR - visibility duration
    dret[25] = l_moon  # LMoon - crescent length
    dret[26] = phase_angle  # CVAact (using phase angle)
    dret[27] = illumination  # Illum - illumination percentage
    # dret[28] onwards are reserved, already 0.0

    return tuple(dret), flags


# =============================================================================
# VISUAL LIMITING MAGNITUDE CALCULATIONS
# =============================================================================


def vis_limit_mag(
    jd: float,
    geopos: tuple,
    atmo: tuple,
    observer: tuple,
    objname: str,
    flags: int = FLG_SWIEPH,
) -> Tuple[float, Tuple[float, ...]]:
    """
    Calculate the limiting visual magnitude for observing a celestial body.

    This function determines whether a celestial body (planet, star, etc.)
    is visible given the current sky brightness conditions, atmospheric
    parameters, and observer characteristics. It returns both the visibility
    status and detailed information about the observation conditions.

    Args:
        jd: Julian Day (UT) for the observation time
        geopos: Geographic position as a sequence:
            - [0]: Geographic longitude in degrees (east positive)
            - [1]: Geographic latitude in degrees (north positive)
            - [2]: Altitude above sea level in meters
        atmo: Atmospheric conditions as a sequence:
            - [0]: Atmospheric pressure in mbar/hPa
            - [1]: Atmospheric temperature in degrees Celsius
            - [2]: Relative humidity in percent (0-100)
            - [3]: If >= 1: Meteorological Range in km
                   If 0-1: Total atmospheric coefficient (ktot)
                   If 0: Compute ktot from other parameters
        observer: Observer characteristics as a sequence:
            - [0]: Age of observer in years (default 36)
            - [1]: Snellen ratio of observer's eyes (default 1.0 = normal)
            For optical instruments (when HELFLAG_OPTICAL_PARAMS is set):
            - [2]: 0 = monocular, 1 = binocular
            - [3]: Telescope magnification (0 = naked eye)
            - [4]: Optical aperture (telescope diameter) in mm
            - [5]: Optical transmission coefficient
        objname: Name of the object to observe. Can be:
            - Planet name (e.g., "Venus", "Mars", "Jupiter")
            - Fixed star name (e.g., "Sirius", "Aldebaran")
            - Planet number as string (e.g., "2" for Venus)
        flags: Calculation flags combining ephemeris and heliacal flags:
            - FLG_SWIEPH, FLG_JPLEPH, etc. for ephemeris
            - HELFLAG_OPTICAL_PARAMS: Use optical instrument parameters
            - HELFLAG_NO_DETAILS: Skip detailed calculations
            - HELFLAG_VISLIM_DARK: Assume Sun at nadir (dark sky)
            - HELFLAG_VISLIM_NOMOON: Exclude Moon's brightness contribution

    Returns:
        Tuple containing:
            - result: Visibility status code:
                - (-2): Object is below horizon
                - (0): OK, photopic vision (bright conditions)
                - (1): OK, scotopic vision (dark conditions)
                - (2): OK, near limit between photopic/scotopic
            - dret: Tuple of 8 floats with observation details:
                - [0]: Limiting visual magnitude (object visible if mag < this)
                - [1]: Altitude of object in degrees
                - [2]: Azimuth of object in degrees
                - [3]: Altitude of Sun in degrees
                - [4]: Azimuth of Sun in degrees
                - [5]: Altitude of Moon in degrees
                - [6]: Azimuth of Moon in degrees
                - [7]: Magnitude of object

    Raises:
        ValueError: If objname is empty or invalid

    Algorithm:
        The limiting magnitude calculation is based on Schaefer's model (1990)
        which considers:
        1. Sky brightness from Sun, Moon, zodiacal light, airglow
        2. Atmospheric extinction based on airmass and conditions
        3. Observer's eye adaptation (scotopic vs photopic)
        4. Optional optical instrument characteristics

        The sky background brightness varies with:
        - Sun altitude (twilight contribution)
        - Moon altitude and phase (moonlight)
        - Atmospheric scattering

    Example:
        >>> from libephemeris import julday, vis_limit_mag
        >>> jd = julday(2024, 8, 15, 22.0)
        >>> # Rome location
        >>> geopos = (12.5, 41.9, 0)
        >>> # Standard atmosphere
        >>> atmo = (1013.25, 15.0, 50.0, 0.0)
        >>> # Normal observer
        >>> observer = (36, 1.0)
        >>> result, dret = vis_limit_mag(jd, geopos, atmo, observer, "Venus")
        >>> if dret[0] > dret[7]:
        ...     print("Venus is visible")
        ... else:
        ...     print("Venus is not visible")

    References:
        - Reference API: vis_limit_mag()
        - Schaefer, B.E. (1990) "Telescopic Limiting Magnitudes"
        - Schaefer, B.E. (1993) "Astronomy and the Limits of Vision"
    """
    import math
    from .planets import _PLANET_MAP, get_planet_target, pheno_ut
    from .fixed_stars import fixstar2_ut, fixstar2_mag
    from .state import get_planets, get_timescale
    from .constants import (
        SUN,
        MOON,
        HELFLAG_VISLIM_DARK,
        HELFLAG_VISLIM_NOMOON,
        HELFLAG_BELOW_HORIZON,
        HELFLAG_PHOTOPIC,
        HELFLAG_SCOTOPIC,
        HELFLAG_MIXED,
    )
    from skyfield.api import wgs84

    if not objname:
        raise ValueError("objname cannot be empty")

    # Parse geographic position
    lon = geopos[0] if len(geopos) > 0 else 0.0
    lat = geopos[1] if len(geopos) > 1 else 0.0
    alt_m = geopos[2] if len(geopos) > 2 else 0.0

    # Parse atmospheric conditions
    pressure = atmo[0] if len(atmo) > 0 else 1013.25
    temperature = atmo[1] if len(atmo) > 1 else 15.0
    humidity_pct = atmo[2] if len(atmo) > 2 else 50.0
    met_range = atmo[3] if len(atmo) > 3 else 0.0

    # Parse observer data
    observer_age = observer[0] if len(observer) > 0 else 36.0
    snellen_ratio = observer[1] if len(observer) > 1 else 1.0

    # Get ephemeris and timescale
    eph = get_planets()
    ts = get_timescale()

    # Create observer location
    obs_location = wgs84.latlon(lat, lon, alt_m)
    earth = eph["earth"]
    sun = eph["sun"]
    moon = eph["moon"]

    t = ts.ut1_jd(jd)
    observer_at = earth + obs_location

    # Calculate Sun position
    sun_app = observer_at.at(t).observe(sun).apparent()
    sun_alt_deg, sun_az_deg, _ = sun_app.altaz()
    sun_alt = sun_alt_deg.degrees
    sun_az = sun_az_deg.degrees

    # Calculate Moon position
    moon_app = observer_at.at(t).observe(moon).apparent()
    moon_alt_deg, moon_az_deg, _ = moon_app.altaz()
    moon_alt = moon_alt_deg.degrees
    moon_az = moon_az_deg.degrees

    # Determine if objname is a planet ID or name
    body_id = None
    is_fixed_star = False

    # Try parsing as integer (planet ID)
    try:
        body_id = int(objname)
    except ValueError:
        # Try to find planet by name
        name_upper = objname.upper().strip()
        planet_names = {
            "SUN": SUN,
            "MOON": MOON,
            "MERCURY": 2,
            "VENUS": 3,
            "MARS": 4,
            "JUPITER": 5,
            "SATURN": 6,
            "URANUS": 7,
            "NEPTUNE": 8,
            "PLUTO": 9,
        }
        if name_upper in planet_names:
            body_id = planet_names[name_upper]
        else:
            # Assume it's a fixed star
            is_fixed_star = True

    # Calculate object position and magnitude
    obj_alt = 0.0
    obj_az = 0.0
    obj_mag = 0.0

    if is_fixed_star:
        # Fixed star calculation
        try:
            star_result, star_name_out, retflag = fixstar2_ut(
                objname, jd, flags & 0xFF
            )

            # star_result is (lon, lat, dist, lon_speed, lat_speed, dist_speed)
            # We need to convert ecliptic to horizontal
            star_lon = star_result[0]
            star_lat = star_result[1]

            # Get star magnitude
            try:
                star_mag_val, _star_name_mag = fixstar2_mag(objname)
                obj_mag = star_mag_val
            except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
                _reraise_if_leb_range_error(_exc)
                obj_mag = 2.0  # Default magnitude if not found

            # Convert ecliptic to equatorial then to horizontal
            # Simplified: use azalt function if available
            from .utils import azalt, ECL2HOR

            hor_result = azalt(
                jd,
                ECL2HOR,
                (lon, lat, alt_m),
                pressure,
                temperature,
                (star_lon, star_lat, 1.0),
            )
            obj_az = hor_result[0]
            obj_alt = hor_result[1]

        except ValueError:
            raise
        except (ValueError, ArithmeticError, IndexError) as e:
            # Star not found or other error
            raise ValueError(f"could not find star name {objname.lower()}: {e}")
    else:
        # Planet calculation
        if body_id is None:
            raise ValueError(f"Unknown object: {objname}")

        # Get planet name from _PLANET_MAP
        if body_id in _PLANET_MAP:
            target_name = _PLANET_MAP[body_id]
            target = get_planet_target(eph, target_name)

            # Calculate position
            body_app = observer_at.at(t).observe(target).apparent()
            body_alt_deg, body_az_deg, _ = body_app.altaz()
            obj_alt = body_alt_deg.degrees
            obj_az = body_az_deg.degrees

            # Get magnitude from pheno
            try:
                pheno_result = pheno_ut(jd, body_id, flags)
                obj_mag = pheno_result[4]  # Visual magnitude
            except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
                _reraise_if_leb_range_error(_exc)
                obj_mag = 0.0  # Default bright
        else:
            raise ValueError(f"illegal planet number {body_id}.")

    # Check if object is below horizon
    if obj_alt < 0:
        # Apply HELFLAG options before returning
        use_dark_sky_early = bool(flags & HELFLAG_VISLIM_DARK)
        exclude_moon_early = bool(flags & HELFLAG_VISLIM_NOMOON)
        ret_sun_alt = -90.0 if use_dark_sky_early else sun_alt
        ret_moon_alt = -90.0 if exclude_moon_early else moon_alt
        dret = (
            0.0,
            obj_alt,
            obj_az,
            ret_sun_alt,
            sun_az,
            ret_moon_alt,
            moon_az,
            obj_mag,
        )
        return float(HELFLAG_BELOW_HORIZON), dret

    # Apply HELFLAG options
    use_dark_sky = bool(flags & HELFLAG_VISLIM_DARK)
    exclude_moon = bool(flags & HELFLAG_VISLIM_NOMOON)

    if use_dark_sky:
        sun_alt = -90.0  # Assume Sun at nadir

    if exclude_moon:
        moon_alt = -90.0  # Assume Moon at nadir

    # Calculate atmospheric extinction coefficient
    # Simplified Schaefer model
    def calc_extinction_coeff(pressure_mbar, temp_c, humidity_pct, met_range_km):
        """Calculate total atmospheric extinction coefficient."""
        # Rayleigh scattering (molecular)
        k_rayleigh = 0.1451 * (pressure_mbar / 1013.25)

        # Aerosol scattering
        if met_range_km >= 1:
            # Meteorological range given
            k_aerosol = 3.912 / met_range_km - 0.106
        elif 0 < met_range_km < 1:
            # ktot given directly
            return met_range_km
        else:
            # Estimate from humidity
            # Higher humidity = more scattering
            k_aerosol = 0.1 + 0.2 * (humidity_pct / 100.0)

        # Ozone absorption (small contribution)
        k_ozone = 0.016

        return k_rayleigh + k_aerosol + k_ozone

    k_ext = calc_extinction_coeff(pressure, temperature, humidity_pct, met_range)

    # Calculate airmass
    def calc_airmass(altitude_deg):
        """Calculate airmass using Kasten & Young formula."""
        if altitude_deg <= 0:
            return 40.0  # Maximum airmass at/below horizon
        z = 90.0 - altitude_deg  # Zenith angle
        z_rad = math.radians(z)
        # Kasten & Young (1989)
        airmass = 1.0 / (math.cos(z_rad) + 0.50572 * (96.07995 - z) ** (-1.6364))
        return min(airmass, 40.0)

    airmass = calc_airmass(obj_alt)

    # Calculate sky brightness contribution from Sun
    def calc_sun_brightness_contribution(sun_altitude):
        """Calculate sky brightness factor from Sun position."""
        if sun_altitude >= 0:
            return 1e10  # Daylight - essentially infinite brightness
        elif sun_altitude >= -6:
            # Civil twilight
            return 10 ** (4.0 + sun_altitude * 0.4)
        elif sun_altitude >= -12:
            # Nautical twilight
            return 10 ** (1.6 + (sun_altitude + 6) * 0.2)
        elif sun_altitude >= -18:
            # Astronomical twilight
            return 10 ** (0.4 + (sun_altitude + 12) * 0.1)
        else:
            # Full night
            return 1.0

    # Calculate sky brightness contribution from Moon
    def calc_moon_brightness_contribution(moon_altitude, jd_ut):
        """Calculate sky brightness factor from Moon position and phase."""
        if moon_altitude <= 0:
            return 0.0

        # Get Moon phase
        try:
            moon_pheno = pheno_ut(jd_ut, MOON, flags & 0xFF)
            phase_angle = moon_pheno[0]  # Phase angle in degrees
            illumination = (1 + math.cos(math.radians(phase_angle))) / 2.0
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            illumination = 0.5  # Assume half moon

        # Moon brightness increases with altitude and illumination
        moon_factor = illumination * math.sin(math.radians(moon_altitude))
        return moon_factor * 100.0

    sun_contrib = calc_sun_brightness_contribution(sun_alt)
    moon_contrib = (
        calc_moon_brightness_contribution(moon_alt, jd) if not exclude_moon else 0.0
    )

    # Total sky brightness (arbitrary units)
    total_sky_brightness = sun_contrib + moon_contrib + 1.0  # 1.0 = natural sky glow

    # Calculate limiting magnitude
    # Base limiting magnitude for perfect conditions: ~6.5 for naked eye
    base_limit_mag = 6.5

    # Adjust for observer characteristics
    # Age reduces sensitivity (approximately 0.05 mag per decade over 30)
    age_factor = max(0, (observer_age - 30) / 10.0) * 0.05

    # Snellen ratio adjustment (better eyes = higher limit)
    snellen_factor = -2.5 * math.log10(snellen_ratio) if snellen_ratio > 0 else 0

    # Airmass extinction
    airmass_extinction = k_ext * airmass

    # Sky brightness reduction of limiting magnitude
    if total_sky_brightness > 1:
        sky_reduction = 2.5 * math.log10(total_sky_brightness)
    else:
        sky_reduction = 0.0

    # Calculate final limiting magnitude
    limiting_mag = base_limit_mag - age_factor - snellen_factor - sky_reduction

    # Apply extinction to object magnitude (object appears fainter)
    apparent_obj_mag = obj_mag + airmass_extinction

    # Determine vision type based on sky brightness
    if sun_alt >= -6:
        # Photopic vision (daylight/bright twilight)
        vision_type = HELFLAG_PHOTOPIC
    elif sun_alt >= -12:
        # Mixed/mesopic vision
        vision_type = HELFLAG_MIXED
    else:
        # Scotopic vision (night)
        vision_type = HELFLAG_SCOTOPIC

    # Build result tuple
    dret = (
        limiting_mag,  # 0: Limiting visual magnitude
        obj_alt,  # 1: Altitude of object
        obj_az,  # 2: Azimuth of object
        sun_alt,  # 3: Altitude of Sun
        sun_az,  # 4: Azimuth of Sun
        moon_alt,  # 5: Altitude of Moon
        moon_az,  # 6: Azimuth of Moon
        apparent_obj_mag,  # 7: Apparent magnitude of object (with extinction)
    )

    return float(vision_type), dret


# =============================================================================
# BESSELIAN ELEMENTS
# =============================================================================

# Earth's equatorial radius in km (WGS84)
EARTH_RADIUS_KM = 6378.137


from dataclasses import dataclass


@dataclass
class BesselianElements:
    """
    Besselian elements for a solar eclipse at a given reference time.

    Besselian elements are the fundamental quantities used to calculate
    the circumstances of a solar eclipse for any location on Earth.
    They describe the geometry of the Moon's shadow relative to Earth
    in a standardized coordinate system called the fundamental plane.

    The fundamental plane passes through Earth's center and is perpendicular
    to the Moon-Sun line (shadow axis). The coordinate system has:
    - x-axis: pointing east (increasing right ascension)
    - y-axis: pointing north
    - z-axis: along shadow axis toward the Moon

    Attributes:
        t0: Reference time (Julian Day, UT) for these elements.
            Elements are exact at this time; derivatives allow interpolation.
        x: Shadow axis x-coordinate on fundamental plane (Earth radii).
            Positive = shadow axis is east of Earth's center.
        y: Shadow axis y-coordinate on fundamental plane (Earth radii).
            Positive = shadow axis is north of Earth's equator.
        d: Declination of the shadow axis (degrees).
            The angle between the shadow axis and the equatorial plane.
        l1: Radius of the penumbral shadow cone on fundamental plane (Earth radii).
            Observers within this radius see at least a partial eclipse.
        l2: Radius of the umbral/antumbral shadow cone (Earth radii).
            Positive for annular eclipses (antumbra), negative for total (umbra).
        mu: Greenwich hour angle of the shadow axis (degrees).
            The angle measured westward from Greenwich meridian to the
            sub-solar point on the fundamental plane.
        dx_dt: Rate of change of x (Earth radii per hour).
        dy_dt: Rate of change of y (Earth radii per hour).
        dd_dt: Rate of change of d (degrees per hour).
        dl1_dt: Rate of change of l1 (Earth radii per hour).
        dl2_dt: Rate of change of l2 (Earth radii per hour).
        dmu_dt: Rate of change of mu (degrees per hour).

    Example:
        >>> from libephemeris import julday, BesselianElements
        >>> # Create elements manually for illustration
        >>> elements = BesselianElements(
        ...     t0=julday(2024, 4, 8, 18.0),
        ...     x=0.3145, y=0.2731, d=7.5821,
        ...     l1=0.5436, l2=-0.0047, mu=89.1234,
        ...     dx_dt=0.5123, dy_dt=0.1456, dd_dt=0.0012,
        ...     dl1_dt=-0.0001, dl2_dt=-0.0001, dmu_dt=15.0041
        ... )
        >>> print(f"Shadow at x={elements.x:.4f}, y={elements.y:.4f}")

    Note:
        All angular rates (dd_dt, dmu_dt) are in degrees per hour.
        All linear rates (dx_dt, dy_dt, dl1_dt, dl2_dt) are in Earth radii per hour.
        Elements can be interpolated to nearby times t using:
            x(t) ≈ x + dx_dt * (t - t0) * 24  (where t is in days)

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """

    t0: float  # Reference time (Julian Day, UT)

    # Besselian elements
    x: float  # Shadow x-coordinate (Earth radii)
    y: float  # Shadow y-coordinate (Earth radii)
    d: float  # Declination of shadow axis (degrees)
    l1: float  # Penumbral cone radius (Earth radii)
    l2: float  # Umbral cone radius (Earth radii)
    mu: float  # Greenwich hour angle (degrees)

    # Time derivatives (rates of change per hour)
    dx_dt: float  # Rate of change of x (Earth radii/hour)
    dy_dt: float  # Rate of change of y (Earth radii/hour)
    dd_dt: float  # Rate of change of d (degrees/hour)
    dl1_dt: float  # Rate of change of l1 (Earth radii/hour)
    dl2_dt: float  # Rate of change of l2 (Earth radii/hour)
    dmu_dt: float  # Rate of change of mu (degrees/hour)


def calc_besselian_x(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian x coordinate for a solar eclipse.

    The Besselian x coordinate is the x-component (positive eastward) of the
    Moon's shadow axis intersection with the fundamental plane, measured in
    Earth equatorial radii. The fundamental plane is the plane through Earth's
    center perpendicular to the Moon-Sun line (the shadow axis direction).

    The x-axis points east in the direction of increasing right ascension,
    perpendicular to the shadow axis. This coordinate, along with y, describes
    where the Moon's shadow cone axis pierces the fundamental plane.

    Args:
        jd: Julian Day (UT) at which to calculate the Besselian x coordinate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The Besselian x coordinate in Earth equatorial radii.
        Positive values indicate the shadow axis is east of Earth's center,
        negative values indicate west.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Convert to 3D Cartesian coordinates (AU)
        3. Compute the shadow axis direction unit vector (Sun to Moon direction)
        4. Compute the fundamental plane axes:
           - z-axis: along shadow axis (from Sun toward Moon)
           - x-axis: perpendicular to z, in equatorial plane (toward increasing RA)
           - y-axis: completes right-handed system (toward north)
        5. Project Moon's position onto the fundamental plane
        6. Extract x-component and convert to Earth radii

    Mathematical basis:
        The shadow axis is defined as the line from the Sun's center through
        the Moon's center. The fundamental plane is perpendicular to this axis
        and passes through Earth's center. The Besselian x coordinate is the
        east-west displacement (in Earth radii) of where this axis pierces
        the fundamental plane.

    Precision:
        Accurate to ~0.0001 Earth radii (~0.6 km) for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_x
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> x = calc_besselian_x(jd)
        >>> print(f"Besselian x = {x:.4f} Earth radii")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    from .state import get_leb_reader

    # AU to km conversion
    AU_TO_KM = 149597870.7

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            from .fast_calc import _apparent_icrs_cartesian
            from .time_utils import deltat

            jd_tt = jd + deltat(jd)
            sun_pos = _apparent_icrs_cartesian(reader, jd_tt, SUN)
            moon_pos = _apparent_icrs_cartesian(reader, jd_tt, MOON)
            _leb_ok = True
        except (KeyError, ValueError):
            pass
    if not _leb_ok:
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)

        # Get Earth, Sun, Moon
        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        # Get geocentric positions in AU (ICRS/J2000 equatorial frame)
        earth_at = earth.at(t)

        # Geocentric apparent positions
        sun_apparent = earth_at.observe(sun).apparent()
        moon_apparent = earth_at.observe(moon).apparent()

        # Get Cartesian positions in AU
        sun_pos = sun_apparent.position.au  # [x, y, z] in AU
        moon_pos = moon_apparent.position.au  # [x, y, z] in AU

    # Shadow axis direction: from Sun toward Moon
    # This is the direction of the shadow cone axis
    shadow_dir = [
        moon_pos[0] - sun_pos[0],
        moon_pos[1] - sun_pos[1],
        moon_pos[2] - sun_pos[2],
    ]

    # Normalize to get unit vector
    shadow_len = math.sqrt(shadow_dir[0] ** 2 + shadow_dir[1] ** 2 + shadow_dir[2] ** 2)
    shadow_unit = [
        shadow_dir[0] / shadow_len,
        shadow_dir[1] / shadow_len,
        shadow_dir[2] / shadow_len,
    ]

    # The fundamental plane passes through Earth's center and is perpendicular
    # to the shadow axis. We need to define an x-axis in this plane pointing
    # eastward (in the direction of increasing RA).

    # The conventional Besselian x-axis is perpendicular to the shadow axis
    # and lies in the equatorial plane (or as close to it as possible).
    # We use cross product: x_axis = k × shadow_axis (where k = [0, 0, 1])
    # This gives a vector pointing eastward, perpendicular to shadow axis.

    # Cross product of [0, 0, 1] × shadow_unit
    x_axis_raw = [-shadow_unit[1], shadow_unit[0], 0.0]

    # Normalize x-axis
    x_axis_len = math.sqrt(x_axis_raw[0] ** 2 + x_axis_raw[1] ** 2)

    # Handle edge case where shadow axis is nearly parallel to z-axis
    if x_axis_len < 1e-10:
        x_axis = [1.0, 0.0, 0.0]
    else:
        x_axis = [x_axis_raw[0] / x_axis_len, x_axis_raw[1] / x_axis_len, 0.0]

    # Project Moon's position onto the fundamental plane
    # Dot product of Moon position with shadow axis unit vector
    moon_along_axis = (
        moon_pos[0] * shadow_unit[0]
        + moon_pos[1] * shadow_unit[1]
        + moon_pos[2] * shadow_unit[2]
    )

    # The perpendicular component (projection onto fundamental plane)
    # is the Moon position minus the along-axis component
    moon_perp = [
        moon_pos[0] - moon_along_axis * shadow_unit[0],
        moon_pos[1] - moon_along_axis * shadow_unit[1],
        moon_pos[2] - moon_along_axis * shadow_unit[2],
    ]

    # The x coordinate is the dot product with the x-axis
    x_au = (
        moon_perp[0] * x_axis[0] + moon_perp[1] * x_axis[1] + moon_perp[2] * x_axis[2]
    )

    # Convert from AU to Earth radii
    x_km = x_au * AU_TO_KM
    x_earth_radii = x_km / EARTH_RADIUS_KM

    return x_earth_radii


def calc_besselian_y(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian y coordinate for a solar eclipse.

    The Besselian y coordinate is the y-component (positive northward) of the
    Moon's shadow axis intersection with the fundamental plane, measured in
    Earth equatorial radii. The fundamental plane is the plane through Earth's
    center perpendicular to the Moon-Sun line (the shadow axis direction).

    The y-axis points north, perpendicular to both the shadow axis and the
    x-axis. This coordinate, along with x, describes where the Moon's shadow
    cone axis pierces the fundamental plane.

    Args:
        jd: Julian Day (UT) at which to calculate the Besselian y coordinate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The Besselian y coordinate in Earth equatorial radii.
        Positive values indicate the shadow axis is north of Earth's center,
        negative values indicate south.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Convert to 3D Cartesian coordinates (AU)
        3. Compute the shadow axis direction unit vector (Sun to Moon direction)
        4. Compute the fundamental plane axes:
           - z-axis: along shadow axis (from Sun toward Moon)
           - x-axis: perpendicular to z, in equatorial plane (toward increasing RA)
           - y-axis: completes right-handed system (toward north)
        5. Project Moon's position onto the fundamental plane
        6. Extract y-component and convert to Earth radii

    Mathematical basis:
        The shadow axis is defined as the line from the Sun's center through
        the Moon's center. The fundamental plane is perpendicular to this axis
        and passes through Earth's center. The Besselian y coordinate is the
        north-south displacement (in Earth radii) of where this axis pierces
        the fundamental plane.

    Precision:
        Accurate to ~0.0001 Earth radii (~0.6 km) for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_y
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> y = calc_besselian_y(jd)
        >>> print(f"Besselian y = {y:.4f} Earth radii")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    from .state import get_leb_reader

    # AU to km conversion
    AU_TO_KM = 149597870.7

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            from .fast_calc import _apparent_icrs_cartesian
            from .time_utils import deltat

            jd_tt = jd + deltat(jd)
            sun_pos = _apparent_icrs_cartesian(reader, jd_tt, SUN)
            moon_pos = _apparent_icrs_cartesian(reader, jd_tt, MOON)
            _leb_ok = True
        except (KeyError, ValueError):
            pass
    if not _leb_ok:
        from .state import get_planets, get_timescale

        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)

        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        earth_at = earth.at(t)

        sun_apparent = earth_at.observe(sun).apparent()
        moon_apparent = earth_at.observe(moon).apparent()

        sun_pos = sun_apparent.position.au
        moon_pos = moon_apparent.position.au

    # Shadow axis direction: from Sun toward Moon
    shadow_dir = [
        moon_pos[0] - sun_pos[0],
        moon_pos[1] - sun_pos[1],
        moon_pos[2] - sun_pos[2],
    ]

    # Normalize to get unit vector
    shadow_len = math.sqrt(shadow_dir[0] ** 2 + shadow_dir[1] ** 2 + shadow_dir[2] ** 2)
    shadow_unit = [
        shadow_dir[0] / shadow_len,
        shadow_dir[1] / shadow_len,
        shadow_dir[2] / shadow_len,
    ]

    x_axis_raw = [-shadow_unit[1], shadow_unit[0], 0.0]

    # Normalize x-axis
    x_axis_len = math.sqrt(x_axis_raw[0] ** 2 + x_axis_raw[1] ** 2)

    # Handle edge case where shadow axis is nearly parallel to z-axis
    # (extremely rare for solar eclipses)
    if x_axis_len < 1e-10:
        # Shadow axis nearly vertical - use y-axis as x instead
        x_axis = [1.0, 0.0, 0.0]
    else:
        x_axis = [x_axis_raw[0] / x_axis_len, x_axis_raw[1] / x_axis_len, 0.0]

    # Compute y-axis to complete right-handed system: y = shadow_unit × x_axis
    # This gives a vector pointing toward the north celestial pole in the
    # fundamental plane (perpendicular to both shadow axis and x-axis)
    y_axis = [
        shadow_unit[1] * x_axis[2] - shadow_unit[2] * x_axis[1],
        shadow_unit[2] * x_axis[0] - shadow_unit[0] * x_axis[2],
        shadow_unit[0] * x_axis[1] - shadow_unit[1] * x_axis[0],
    ]

    # Project Moon's position onto the fundamental plane
    # The Moon's position in the fundamental plane coordinate system:
    # - Find how far Moon is along the shadow axis from Earth's center
    # - The projection onto the fundamental plane gives us the displacement

    # Dot product of Moon position with shadow axis unit vector
    # This gives the component along the shadow axis
    moon_along_axis = (
        moon_pos[0] * shadow_unit[0]
        + moon_pos[1] * shadow_unit[1]
        + moon_pos[2] * shadow_unit[2]
    )

    # The perpendicular component (projection onto fundamental plane)
    # is the Moon position minus the along-axis component
    moon_perp = [
        moon_pos[0] - moon_along_axis * shadow_unit[0],
        moon_pos[1] - moon_along_axis * shadow_unit[1],
        moon_pos[2] - moon_along_axis * shadow_unit[2],
    ]

    # The y coordinate is the dot product with the y-axis
    y_au = (
        moon_perp[0] * y_axis[0] + moon_perp[1] * y_axis[1] + moon_perp[2] * y_axis[2]
    )

    # Convert from AU to Earth radii
    y_km = y_au * AU_TO_KM
    y_earth_radii = y_km / EARTH_RADIUS_KM

    return y_earth_radii


def calc_besselian_d(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian element d (declination) for a solar eclipse.

    The Besselian element d is the declination of the Moon's shadow axis,
    measured in degrees. In the standard Besselian element convention, the
    shadow axis direction is defined as pointing from the shadow cone vertex
    toward Earth. Since the Sun, Moon, and shadow axis are nearly collinear
    during an eclipse, d closely tracks the Sun's declination.

    The declination d, together with the hour angle mu, defines the orientation
    of the shadow axis in equatorial coordinates.

    Args:
        jd: Julian Day (UT) at which to calculate the declination
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The declination d in degrees.
        Positive values indicate the shadow axis points north of the equator,
        negative values indicate south.
        Range: -90 to +90 degrees.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Convert to 3D Cartesian coordinates (AU)
        3. Compute the shadow axis direction (from Moon toward Sun, which is
           the conventional Besselian axis direction toward the observer)
        4. Normalize and extract the z-component (toward north celestial pole)
        5. Compute declination as arcsin(z) and convert to degrees

    Mathematical basis:
        The shadow axis in Besselian elements points from the shadow cone
        vertex toward Earth. For a unit vector in equatorial coordinates,
        the declination is arcsin(z) where z is the component toward the
        north celestial pole. During a solar eclipse, this declination
        closely matches the Sun's geocentric declination.

    Precision:
        Accurate to better than 0.001 degrees for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_d
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> d = calc_besselian_d(jd)
        >>> print(f"Besselian d = {d:.4f} degrees")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    from .state import get_leb_reader

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            from .fast_calc import _apparent_icrs_cartesian
            from .time_utils import deltat

            jd_tt = jd + deltat(jd)
            sun_pos = _apparent_icrs_cartesian(reader, jd_tt, SUN)
            moon_pos = _apparent_icrs_cartesian(reader, jd_tt, MOON)
            _leb_ok = True
        except (KeyError, ValueError):
            pass
    if not _leb_ok:
        from .state import get_planets, get_timescale

        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)

        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        earth_at = earth.at(t)

        sun_apparent = earth_at.observe(sun).apparent()
        moon_apparent = earth_at.observe(moon).apparent()

        sun_pos = sun_apparent.position.au
        moon_pos = moon_apparent.position.au

    # Shadow axis direction: from Sun toward Moon
    shadow_dir = [
        moon_pos[0] - sun_pos[0],
        moon_pos[1] - sun_pos[1],
        moon_pos[2] - sun_pos[2],
    ]

    # Normalize to get unit vector
    shadow_len = math.sqrt(shadow_dir[0] ** 2 + shadow_dir[1] ** 2 + shadow_dir[2] ** 2)

    # In Besselian element convention, d is the declination of the shadow axis
    # pointing FROM the Moon TOWARD Earth (the direction the shadow travels).
    # The shadow_dir computed above (Moon - Sun) points in the opposite direction,
    # so we negate the z-component when computing the declination.
    shadow_unit_z = -shadow_dir[2] / shadow_len

    # The declination is the angle from the equatorial plane
    # For a unit vector in equatorial coordinates, dec = arcsin(z)
    # Clamp to [-1, 1] to avoid numerical issues with arcsin
    shadow_unit_z = max(-1.0, min(1.0, shadow_unit_z))
    d_rad = math.asin(shadow_unit_z)

    # Convert to degrees
    d_deg = math.degrees(d_rad)

    return d_deg


def calc_besselian_l1(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian element l1 (penumbral shadow radius) for a solar eclipse.

    The Besselian element l1 is the radius of the penumbral shadow cone where it
    intersects the fundamental plane, measured in Earth equatorial radii. The
    fundamental plane is the plane through Earth's center perpendicular to the
    Moon-Sun line (the shadow axis direction).

    The penumbral cone is formed by the external tangent lines from the Sun's limb
    to the Moon's limb. Observers within this radius on the fundamental plane
    would see some portion of the Sun obscured by the Moon.

    Args:
        jd: Julian Day (UT) at which to calculate l1
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The penumbral shadow radius l1 in Earth equatorial radii.
        Always positive. Typical values range from 0.5 to 0.6 Earth radii
        during solar eclipses.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Calculate the Sun-Moon distance
        3. Calculate the penumbral cone semi-angle f1 using:
           sin(f1) = (R_sun + R_moon) / D_sun-moon
        4. Calculate the penumbral cone vertex distance from the Sun:
           k1 = R_sun / sin(f1)
        5. Calculate the distance from the vertex to Earth's center (fundamental plane):
           d = D_sun - k1
        6. The penumbral radius at the fundamental plane is:
           l1 = d * tan(f1)
        7. Convert to Earth radii

    Mathematical basis:
        The penumbral cone's semi-angle f1 is determined by the external
        tangent geometry. The cone vertex is located behind the Sun at a
        distance k1 = R_sun / sin(f1) from the Sun center. The penumbral
        shadow radius at any plane perpendicular to the shadow axis is
        simply the distance from that plane to the vertex multiplied by
        tan(f1).

        For the fundamental plane (at Earth's center), the distance from
        the vertex is approximately D_sun - k1, giving:
            l1 = (D_sun - k1) * tan(f1)

    Physical constants used:
        - Sun radius: 696,000 km (IAU nominal)
        - Moon radius: 1,737.4 km (mean radius)
        - Earth radius: 6,378.137 km (WGS84 equatorial)
        - AU: 149,597,870.7 km

    Precision:
        Accurate to ~0.0001 Earth radii (~0.6 km) for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_l1
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> l1 = calc_besselian_l1(jd)
        >>> print(f"Besselian l1 = {l1:.4f} Earth radii")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    from .state import get_leb_reader

    # Physical constants
    AU_TO_KM = 149597870.7  # km per AU
    SUN_RADIUS_KM = 696000.0  # IAU nominal solar radius
    MOON_RADIUS_KM = 1737.4  # Mean lunar radius

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            from .fast_calc import _apparent_icrs_cartesian
            from .time_utils import deltat

            jd_tt = jd + deltat(jd)
            sun_pos = _apparent_icrs_cartesian(reader, jd_tt, SUN)
            moon_pos = _apparent_icrs_cartesian(reader, jd_tt, MOON)
            _leb_ok = True
        except (KeyError, ValueError):
            pass
    if not _leb_ok:
        from .state import get_planets, get_timescale

        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)

        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        earth_at = earth.at(t)

        sun_apparent = earth_at.observe(sun).apparent()
        moon_apparent = earth_at.observe(moon).apparent()

        sun_pos = sun_apparent.position.au
        moon_pos = moon_apparent.position.au

    # Calculate distances
    # Sun distance from Earth (geocentric)
    sun_dist_au = math.sqrt(sun_pos[0] ** 2 + sun_pos[1] ** 2 + sun_pos[2] ** 2)

    # Moon distance from Earth (geocentric)
    moon_dist_au = math.sqrt(moon_pos[0] ** 2 + moon_pos[1] ** 2 + moon_pos[2] ** 2)

    # Shadow axis direction: from Sun toward Moon
    shadow_dir = [
        moon_pos[0] - sun_pos[0],
        moon_pos[1] - sun_pos[1],
        moon_pos[2] - sun_pos[2],
    ]

    # Distance from Sun to Moon (in AU)
    sun_moon_dist_au = math.sqrt(
        shadow_dir[0] ** 2 + shadow_dir[1] ** 2 + shadow_dir[2] ** 2
    )

    # Convert to km
    sun_dist_km = sun_dist_au * AU_TO_KM
    sun_moon_dist_km = sun_moon_dist_au * AU_TO_KM

    # Calculate the penumbral cone semi-angle f1
    # The penumbral cone is formed by external tangent lines from Sun's limb to Moon's limb
    # sin(f1) = (R_sun + R_moon) / D_sun-moon
    sum_radii = SUN_RADIUS_KM + MOON_RADIUS_KM
    sin_f1 = sum_radii / sun_moon_dist_km

    # Clamp to valid range for asin (numerical safety)
    sin_f1 = max(-1.0, min(1.0, sin_f1))
    f1_rad = math.asin(sin_f1)
    tan_f1 = math.tan(f1_rad)

    # The penumbral cone vertex is located at distance k1 from the Sun center
    # along the shadow axis (toward the Moon):
    # k1 = R_sun / sin(f1)
    k1_km = SUN_RADIUS_KM / sin_f1

    # The distance from the penumbral vertex to Earth's center (fundamental plane)
    # is the Sun's distance minus the vertex offset
    vertex_to_earth_km = sun_dist_km - k1_km

    # The penumbral shadow radius at the fundamental plane is:
    # l1 = vertex_to_earth * tan(f1)
    l1_km = vertex_to_earth_km * tan_f1

    # Convert to Earth radii
    l1_earth_radii = l1_km / EARTH_RADIUS_KM

    return l1_earth_radii


def calc_besselian_l2(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian element l2 (umbral/antumbral shadow radius) for a solar eclipse.

    The Besselian element l2 is the radius of the umbral (or antumbral) shadow cone
    where it intersects the fundamental plane, measured in Earth equatorial radii.
    The fundamental plane is the plane through Earth's center perpendicular to the
    Moon-Sun line (the shadow axis direction).

    The umbral cone is formed by the internal tangent lines from the Sun's limb
    to the Moon's limb. The sign of l2 indicates the eclipse type:
    - l2 > 0: Umbral cone apex is beyond Earth (total eclipse possible)
    - l2 < 0: Umbral cone apex is between Moon and Earth (annular eclipse)

    Args:
        jd: Julian Day (UT) at which to calculate l2
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The umbral/antumbral shadow radius l2 in Earth equatorial radii.
        Positive for total eclipses (umbral shadow), negative for annular
        eclipses (antumbral shadow). Typical absolute values range from
        0.003 to 0.05 Earth radii during solar eclipses.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Calculate the Sun-Moon distance
        3. Calculate the umbral cone semi-angle f2 using:
           sin(f2) = (R_sun - R_moon) / D_sun-moon
        4. Calculate the umbral cone vertex distance from the Moon:
           c2 = R_moon / sin(f2)  (distance from Moon to umbral vertex)
        5. Calculate whether vertex is beyond or before Earth:
           If Moon distance > c2: umbral (total), l2 > 0
           If Moon distance < c2: antumbral (annular), l2 < 0
        6. The umbral radius at the fundamental plane is:
           l2 = ±|d_vertex_to_earth| * tan(f2)
        7. Convert to Earth radii

    Mathematical basis:
        The umbral cone's semi-angle f2 is determined by the internal
        tangent geometry. The cone vertex is located at distance c2 from
        the Moon center along the shadow axis (toward Earth):
        c2 = R_moon / sin(f2)

        For the fundamental plane (at Earth's center), the distance from
        the vertex determines the shadow radius. If the vertex is beyond
        Earth (umbral), l2 is positive. If the vertex is between Moon and
        Earth (antumbral), l2 is negative, indicating the shadow has
        diverged from the cone and forms an antumbral zone.

    Physical constants used:
        - Sun radius: 696,000 km (IAU nominal)
        - Moon radius: 1,737.4 km (mean radius)
        - Earth radius: 6,378.137 km (WGS84 equatorial)
        - AU: 149,597,870.7 km

    Precision:
        Accurate to ~0.0001 Earth radii (~0.6 km) for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions.

    Example:
        >>> from libephemeris import julday, calc_besselian_l2
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> l2 = calc_besselian_l2(jd)
        >>> print(f"Besselian l2 = {l2:.4f} Earth radii")
        >>> if l2 > 0:
        ...     print("Total eclipse (umbral shadow)")
        ... else:
        ...     print("Annular eclipse (antumbral shadow)")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    from .state import get_leb_reader

    # Physical constants
    AU_TO_KM = 149597870.7  # km per AU
    SUN_RADIUS_KM = 696000.0  # IAU nominal solar radius
    MOON_RADIUS_KM = 1737.4  # Mean lunar radius

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            from .fast_calc import _apparent_icrs_cartesian
            from .time_utils import deltat

            jd_tt = jd + deltat(jd)
            sun_pos = _apparent_icrs_cartesian(reader, jd_tt, SUN)
            moon_pos = _apparent_icrs_cartesian(reader, jd_tt, MOON)
            _leb_ok = True
        except (KeyError, ValueError):
            pass
    if not _leb_ok:
        from .state import get_planets, get_timescale

        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)

        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        earth_at = earth.at(t)

        sun_apparent = earth_at.observe(sun).apparent()
        moon_apparent = earth_at.observe(moon).apparent()

        sun_pos = sun_apparent.position.au
        moon_pos = moon_apparent.position.au

    # Calculate distances
    # Moon distance from Earth (geocentric)
    moon_dist_au = math.sqrt(moon_pos[0] ** 2 + moon_pos[1] ** 2 + moon_pos[2] ** 2)

    # Shadow axis direction: from Sun toward Moon
    shadow_dir = [
        moon_pos[0] - sun_pos[0],
        moon_pos[1] - sun_pos[1],
        moon_pos[2] - sun_pos[2],
    ]

    # Distance from Sun to Moon (in AU)
    sun_moon_dist_au = math.sqrt(
        shadow_dir[0] ** 2 + shadow_dir[1] ** 2 + shadow_dir[2] ** 2
    )

    # Convert to km
    moon_dist_km = moon_dist_au * AU_TO_KM
    sun_moon_dist_km = sun_moon_dist_au * AU_TO_KM

    # Calculate the umbral cone semi-angle f2
    # The umbral cone is formed by internal tangent lines from Sun's limb to Moon's limb
    # sin(f2) = (R_sun - R_moon) / D_sun-moon
    diff_radii = SUN_RADIUS_KM - MOON_RADIUS_KM
    sin_f2 = diff_radii / sun_moon_dist_km

    # Clamp to valid range for asin (numerical safety)
    sin_f2 = max(-1.0, min(1.0, sin_f2))
    f2_rad = math.asin(sin_f2)
    tan_f2 = math.tan(f2_rad)

    # The umbral cone vertex is located at distance c2 from the Moon center
    # along the shadow axis (toward Earth):
    # c2 = R_moon / sin(f2)
    c2_km = MOON_RADIUS_KM / sin_f2

    # The distance from Earth's center to the umbral vertex
    # c2 is the distance from Moon to the umbral cone vertex (toward Earth)
    # For total eclipse: vertex is beyond Earth, so c2 > moon_dist → l2 > 0
    # For annular eclipse: vertex is between Moon and Earth, so c2 < moon_dist → l2 < 0
    vertex_to_earth_km = c2_km - moon_dist_km

    # The umbral/antumbral shadow radius at the fundamental plane is:
    # l2 = vertex_to_earth * tan(f2)
    # Positive l2 = umbral shadow (total eclipse)
    # Negative l2 = antumbral shadow (annular eclipse)
    l2_km = vertex_to_earth_km * tan_f2

    # Convert to Earth radii
    l2_earth_radii = l2_km / EARTH_RADIUS_KM

    return l2_earth_radii


def calc_besselian_mu(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the Besselian element mu (Greenwich hour angle) for a solar eclipse.

    The Besselian element mu is the Greenwich hour angle of the shadow axis,
    measured in degrees. It represents the angle between the Greenwich meridian
    and the hour circle containing the shadow axis, measured westward from
    Greenwich along the celestial equator.

    Together with the declination d, mu defines the orientation of the shadow
    axis in the equatorial coordinate system. As Earth rotates, mu increases
    at approximately 15 degrees per hour (360 degrees per sidereal day).

    Args:
        jd: Julian Day (UT) at which to calculate mu
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The Greenwich hour angle mu in degrees.
        Range: 0 to 360 degrees, measured westward from Greenwich.
        Values increase with time as Earth rotates eastward.

    Algorithm:
        1. Get geocentric equatorial positions of Moon and Sun
        2. Convert to 3D Cartesian coordinates (AU)
        3. Compute the shadow axis direction (from Moon toward Sun, the
           conventional Besselian axis direction toward the observer)
        4. Calculate the right ascension (RA) of the shadow axis from
           its x,y components: RA = atan2(y, x)
        5. Calculate Greenwich Apparent Sidereal Time (GAST)
        6. Compute mu = GAST - RA, normalized to 0-360 degrees

    Mathematical basis:
        The hour angle H of a celestial object is defined as:
            H = LST - RA (Local Sidereal Time minus Right Ascension)

        For the Greenwich meridian, this becomes:
            mu = GAST - RA_shadow

        where RA_shadow is the right ascension of the shadow axis direction.
        The shadow axis in Besselian elements points from the Moon toward
        the Sun (toward the observer on Earth).

    Physical interpretation:
        mu indicates where the shadow axis crosses the celestial equator
        relative to the Greenwich meridian. When mu = 0, the shadow axis
        is on the Greenwich meridian. As Earth rotates, mu increases by
        about 15 degrees per hour.

    Precision:
        Accurate to better than 0.01 degrees for typical eclipses.
        Uses full precision Moon and Sun ephemeris positions and
        accurate sidereal time calculation.

    Example:
        >>> from libephemeris import julday, calc_besselian_mu
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> mu = calc_besselian_mu(jd)
        >>> print(f"Besselian mu = {mu:.4f} degrees")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    from .state import get_leb_reader

    reader = get_leb_reader()
    _leb_ok = False
    if reader is not None:
        try:
            from .fast_calc import _apparent_icrs_cartesian
            from .time_utils import deltat

            jd_tt = jd + deltat(jd)
            sun_pos = _apparent_icrs_cartesian(reader, jd_tt, SUN)
            moon_pos = _apparent_icrs_cartesian(reader, jd_tt, MOON)
            _leb_ok = True
        except (KeyError, ValueError):
            pass
    if not _leb_ok:
        from .state import get_planets, get_timescale

        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)

        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        earth_at = earth.at(t)

        sun_apparent = earth_at.observe(sun).apparent()
        moon_apparent = earth_at.observe(moon).apparent()

        sun_pos = sun_apparent.position.au
        moon_pos = moon_apparent.position.au

    # Shadow axis direction: from Moon toward Sun
    # This is the conventional Besselian axis direction (toward Earth/observer)
    # Same convention as used in calc_besselian_d
    shadow_dir = [
        sun_pos[0] - moon_pos[0],
        sun_pos[1] - moon_pos[1],
        sun_pos[2] - moon_pos[2],
    ]

    # Calculate the right ascension of the shadow axis
    # RA = atan2(y, x) in the equatorial coordinate system
    # The x-axis points toward the vernal equinox, y-axis at RA=90 degrees
    ra_rad = math.atan2(shadow_dir[1], shadow_dir[0])
    ra_deg = math.degrees(ra_rad)

    # Normalize RA to 0-360 degrees
    if ra_deg < 0:
        ra_deg += 360.0

    # Calculate Greenwich Apparent Sidereal Time (GAST)
    from .state import get_timescale as _get_ts_mu

    ts_mu = _get_ts_mu()
    t_mu = ts_mu.ut1_jd(jd)
    gast = t_mu.gast  # Greenwich Apparent Sidereal Time in hours

    # Convert GAST from hours to degrees (1 hour = 15 degrees)
    gast_deg = gast * 15.0

    # Calculate Greenwich hour angle: mu = GAST - RA
    mu_deg = gast_deg - ra_deg

    # Normalize to 0-360 degrees
    mu_deg = mu_deg % 360.0
    if mu_deg < 0:
        mu_deg += 360.0

    return mu_deg


# =============================================================================
# BESSELIAN ELEMENT TIME DERIVATIVES
# =============================================================================
# These functions calculate the hourly rates of change for each Besselian
# element. These derivatives are essential for interpolating element values
# during an eclipse and for calculating eclipse contact times.


def calc_besselian_dx_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element x (dx/dt).

    This function computes the hourly rate of change of the Besselian x
    coordinate, which represents how fast the shadow axis is moving eastward
    (or westward if negative) across the fundamental plane.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dx/dt in Earth radii per hour.
        Positive values indicate eastward motion, negative indicates westward.

    Algorithm:
        Uses symmetric numerical differentiation:
        dx/dt ≈ (x(t+h) - x(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-6 Earth radii per hour for typical
        eclipse conditions. The small time step minimizes truncation error
        while the symmetric difference eliminates first-order errors.

    Example:
        >>> from libephemeris import julday, calc_besselian_dx_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dx_dt = calc_besselian_dx_dt(jd)
        >>> print(f"dx/dt = {dx_dt:.6f} Earth radii/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    # This provides good balance between truncation and round-off errors
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate x at t-h and t+h
    x_minus = calc_besselian_x(jd - h, flags)
    x_plus = calc_besselian_x(jd + h, flags)

    # Symmetric difference quotient
    # dx/dt in Earth radii per day
    dx_dt_per_day = (x_plus - x_minus) / (2 * h)

    # Convert to Earth radii per hour (divide by 24)
    dx_dt_per_hour = dx_dt_per_day / 24.0

    return dx_dt_per_hour


def calc_besselian_dy_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element y (dy/dt).

    This function computes the hourly rate of change of the Besselian y
    coordinate, which represents how fast the shadow axis is moving northward
    (or southward if negative) across the fundamental plane.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dy/dt in Earth radii per hour.
        Positive values indicate northward motion, negative indicates southward.

    Algorithm:
        Uses symmetric numerical differentiation:
        dy/dt ≈ (y(t+h) - y(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-6 Earth radii per hour for typical
        eclipse conditions. The small time step minimizes truncation error
        while the symmetric difference eliminates first-order errors.

    Example:
        >>> from libephemeris import julday, calc_besselian_dy_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dy_dt = calc_besselian_dy_dt(jd)
        >>> print(f"dy/dt = {dy_dt:.6f} Earth radii/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate y at t-h and t+h
    y_minus = calc_besselian_y(jd - h, flags)
    y_plus = calc_besselian_y(jd + h, flags)

    # Symmetric difference quotient
    # dy/dt in Earth radii per day
    dy_dt_per_day = (y_plus - y_minus) / (2 * h)

    # Convert to Earth radii per hour
    dy_dt_per_hour = dy_dt_per_day / 24.0

    return dy_dt_per_hour


def calc_besselian_dd_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element d (dd/dt).

    This function computes the hourly rate of change of the shadow axis
    declination. During an eclipse, the declination changes very slowly
    (typically less than 1 degree per hour) as it tracks the Sun's motion.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dd/dt in degrees per hour.
        Positive values indicate the shadow axis is moving northward,
        negative indicates southward.

    Algorithm:
        Uses symmetric numerical differentiation:
        dd/dt ≈ (d(t+h) - d(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-5 degrees per hour for typical
        eclipse conditions.

    Example:
        >>> from libephemeris import julday, calc_besselian_dd_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dd_dt = calc_besselian_dd_dt(jd)
        >>> print(f"dd/dt = {dd_dt:.6f} degrees/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate d at t-h and t+h
    d_minus = calc_besselian_d(jd - h, flags)
    d_plus = calc_besselian_d(jd + h, flags)

    # Symmetric difference quotient
    # dd/dt in degrees per day
    dd_dt_per_day = (d_plus - d_minus) / (2 * h)

    # Convert to degrees per hour
    dd_dt_per_hour = dd_dt_per_day / 24.0

    return dd_dt_per_hour


def calc_besselian_dl1_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element l1 (dl1/dt).

    This function computes the hourly rate of change of the penumbral
    shadow radius. The penumbral radius changes very slowly during an
    eclipse as the Earth-Moon-Sun geometry evolves.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dl1/dt in Earth radii per hour.
        The value is typically very small (order of 1e-4 to 1e-3).

    Algorithm:
        Uses symmetric numerical differentiation:
        dl1/dt ≈ (l1(t+h) - l1(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-7 Earth radii per hour for typical
        eclipse conditions.

    Example:
        >>> from libephemeris import julday, calc_besselian_dl1_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dl1_dt = calc_besselian_dl1_dt(jd)
        >>> print(f"dl1/dt = {dl1_dt:.8f} Earth radii/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate l1 at t-h and t+h
    l1_minus = calc_besselian_l1(jd - h, flags)
    l1_plus = calc_besselian_l1(jd + h, flags)

    # Symmetric difference quotient
    # dl1/dt in Earth radii per day
    dl1_dt_per_day = (l1_plus - l1_minus) / (2 * h)

    # Convert to Earth radii per hour
    dl1_dt_per_hour = dl1_dt_per_day / 24.0

    return dl1_dt_per_hour


def calc_besselian_dl2_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element l2 (dl2/dt).

    This function computes the hourly rate of change of the umbral/antumbral
    shadow radius. The umbral radius changes very slowly during an eclipse
    as the Earth-Moon-Sun geometry evolves.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dl2/dt in Earth radii per hour.
        The value is typically very small (order of 1e-5 to 1e-4).

    Algorithm:
        Uses symmetric numerical differentiation:
        dl2/dt ≈ (l2(t+h) - l2(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.

    Precision:
        Accurate to approximately 1e-8 Earth radii per hour for typical
        eclipse conditions.

    Example:
        >>> from libephemeris import julday, calc_besselian_dl2_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dl2_dt = calc_besselian_dl2_dt(jd)
        >>> print(f"dl2/dt = {dl2_dt:.8f} Earth radii/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate l2 at t-h and t+h
    l2_minus = calc_besselian_l2(jd - h, flags)
    l2_plus = calc_besselian_l2(jd + h, flags)

    # Symmetric difference quotient
    # dl2/dt in Earth radii per day
    dl2_dt_per_day = (l2_plus - l2_minus) / (2 * h)

    # Convert to Earth radii per hour
    dl2_dt_per_hour = dl2_dt_per_day / 24.0

    return dl2_dt_per_hour


def calc_besselian_dmu_dt(jd: float, flags: int = FLG_SWIEPH) -> float:
    """
    Calculate the time derivative of Besselian element mu (dmu/dt).

    This function computes the hourly rate of change of the Greenwich hour
    angle. Since Earth rotates at approximately 15 degrees per hour, dmu/dt
    is typically close to this value, with small variations due to the
    changing right ascension of the shadow axis.

    The derivative is calculated using symmetric numerical differentiation
    (central difference method) for improved accuracy, with special handling
    for the 0/360 degree wraparound.

    Args:
        jd: Julian Day (UT) at which to calculate the derivative
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        The rate of change dmu/dt in degrees per hour.
        Expected to be approximately 15 degrees/hour (Earth's rotation rate)
        with small variations (typically ±0.04 degrees/hour).

    Algorithm:
        Uses symmetric numerical differentiation:
        dmu/dt ≈ (mu(t+h) - mu(t-h)) / (2h)
        where h = 1/1440 days (1 minute) for optimal accuracy.
        Special handling ensures correct results near the 0/360 boundary.

    Physical interpretation:
        dmu/dt ≈ 15.04 deg/hour (Earth's sidereal rotation rate)
              - d(RA_shadow)/dt (shadow axis RA rate of change)

        The deviation from exactly 15°/hour is due to the apparent
        motion of the Sun (and hence the shadow axis) in right ascension.

    Precision:
        Accurate to approximately 1e-4 degrees per hour for typical
        eclipse conditions.

    Example:
        >>> from libephemeris import julday, calc_besselian_dmu_dt
        >>> # April 8, 2024 total solar eclipse near maximum
        >>> jd = julday(2024, 4, 8, 18.3)  # ~18:18 UTC
        >>> dmu_dt = calc_besselian_dmu_dt(jd)
        >>> print(f"dmu/dt = {dmu_dt:.4f} degrees/hour")

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time step for numerical differentiation: 1 minute in days
    h = 1.0 / 1440.0  # 1 minute = 1/1440 days

    # Calculate mu at t-h and t+h
    mu_minus = calc_besselian_mu(jd - h, flags)
    mu_plus = calc_besselian_mu(jd + h, flags)

    # Handle wraparound at 0/360 degrees
    # mu increases with time (typically ~15 deg/hour)
    # If mu_minus is near 360 and mu_plus is near 0, add 360 to mu_plus
    # If mu_plus is near 360 and mu_minus is near 0, add 360 to mu_minus
    delta_mu = mu_plus - mu_minus

    if delta_mu < -180.0:
        # mu wrapped from ~360 to ~0 going forward
        delta_mu += 360.0
    elif delta_mu > 180.0:
        # mu wrapped from ~0 to ~360 going backward (very rare)
        delta_mu -= 360.0

    # dmu/dt in degrees per day
    dmu_dt_per_day = delta_mu / (2 * h)

    # Convert to degrees per hour
    dmu_dt_per_hour = dmu_dt_per_day / 24.0

    return dmu_dt_per_hour


def interpolate_besselian_elements(
    elements: BesselianElements, target_jd: float
) -> BesselianElements:
    """
    Interpolate Besselian elements to any time during an eclipse.

    Uses first-order Taylor series expansion to interpolate Besselian elements
    from their reference time (t0) to a target time. This is accurate for
    times within approximately ±1-2 hours of the reference time.

    The interpolation formula for each element is:
        element(t) = element(t0) + d_element_dt * (t - t0) * 24

    where t and t0 are in Julian Days, and derivatives are per hour.

    Args:
        elements: BesselianElements object containing values and derivatives
                  at the reference time t0.
        target_jd: Target Julian Day (UT) at which to interpolate elements.

    Returns:
        A new BesselianElements object with interpolated values at target_jd.
        The new object has:
        - t0 set to target_jd (the new reference time)
        - All elements (x, y, d, l1, l2, mu) interpolated to target_jd
        - All derivatives (dx_dt, dy_dt, etc.) copied from the original
          (derivatives change slowly and are valid for nearby times)

    Example:
        >>> from libephemeris import (
        ...     julday, BesselianElements, interpolate_besselian_elements
        ... )
        >>> # Elements at reference time for April 8, 2024 eclipse
        >>> elements = BesselianElements(
        ...     t0=julday(2024, 4, 8, 18.0),
        ...     x=0.3145, y=0.2731, d=7.5821,
        ...     l1=0.5436, l2=-0.0047, mu=89.1234,
        ...     dx_dt=0.5123, dy_dt=0.1456, dd_dt=0.0012,
        ...     dl1_dt=-0.0001, dl2_dt=-0.0001, dmu_dt=15.0041
        ... )
        >>> # Interpolate to 30 minutes later
        >>> later_jd = julday(2024, 4, 8, 18.5)
        >>> interpolated = interpolate_besselian_elements(elements, later_jd)
        >>> print(f"x at t0+0.5h: {interpolated.x:.4f}")

    Note:
        - The mu (Greenwich hour angle) is normalized to [0, 360) degrees
          after interpolation to handle wraparound.
        - For high accuracy over longer time spans, compute fresh elements
          using the individual calc_besselian_* functions.
        - Accuracy degrades quadratically with time from reference:
          - Within ±15 min: error < 0.001 Earth radii
          - Within ±1 hour: error < 0.01 Earth radii
          - Beyond ±2 hours: consider recomputing elements

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Time difference in hours (derivatives are per hour)
    dt_hours = (target_jd - elements.t0) * 24.0

    # Interpolate each element using first-order Taylor series
    x_new = elements.x + elements.dx_dt * dt_hours
    y_new = elements.y + elements.dy_dt * dt_hours
    d_new = elements.d + elements.dd_dt * dt_hours
    l1_new = elements.l1 + elements.dl1_dt * dt_hours
    l2_new = elements.l2 + elements.dl2_dt * dt_hours
    mu_new = elements.mu + elements.dmu_dt * dt_hours

    # Normalize mu to [0, 360) degrees
    mu_new = mu_new % 360.0

    # Return new BesselianElements with interpolated values
    # Derivatives are copied as they change slowly over short time spans
    return BesselianElements(
        t0=target_jd,
        x=x_new,
        y=y_new,
        d=d_new,
        l1=l1_new,
        l2=l2_new,
        mu=mu_new,
        dx_dt=elements.dx_dt,
        dy_dt=elements.dy_dt,
        dd_dt=elements.dd_dt,
        dl1_dt=elements.dl1_dt,
        dl2_dt=elements.dl2_dt,
        dmu_dt=elements.dmu_dt,
    )


def calc_eclipse_first_contact_c1(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of first external contact (C1) for a solar eclipse.

    First contact (C1) is the moment when the Moon's disk first externally
    touches the Sun's disk, marking the beginning of a solar eclipse. At this
    instant, the penumbral shadow cone first touches Earth's surface.

    This function uses Besselian elements to precisely calculate C1. The
    condition for C1 is when gamma (the distance of the shadow axis from
    Earth's center) equals 1 + l1 (Earth radius plus penumbral cone radius),
    occurring before eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function searches backward from
                this time to find C1.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of first contact C1. Returns 0.0 if C1 cannot be
        determined (which would indicate the input time is not near a valid
        solar eclipse).

    Algorithm:
        1. Calculate l1 (penumbral radius) at eclipse maximum
        2. Compute target gamma = 1 + l1 (condition for penumbra touching Earth)
        3. Use binary search to find when gamma equals this target before maximum
        4. The search proceeds from (jd_max - search_range) to jd_max

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the gamma value converges to within 1e-8 Earth radii,
        which corresponds to approximately 0.06 km or 0.04 seconds of time.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, calc_eclipse_first_contact_c1
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate first contact
        >>> jd_c1 = calc_eclipse_first_contact_c1(jd_max)
        >>> print(f"First contact C1: JD {jd_c1:.6f}")

    Note:
        - C1 is also known as "first contact" or "partial phase begins"
        - For local eclipse circumstances (at a specific observer location),
          use _sol_eclipse_when_loc_pythonic() which returns contact times in its result
        - The returned time is for the global eclipse (when penumbra first
          touches Earth anywhere), not for a specific observer location

    See Also:
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_besselian_x, calc_besselian_y: Individual Besselian element functions

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Get l1 (penumbral radius) at maximum eclipse
    l1 = _calc_penumbra_limit(jd_max)

    # For global eclipse, first contact occurs when gamma = 1 + l1
    # (penumbral shadow touches Earth's limb from outside)
    penumbral_limit = 1.0 + l1  # Earth radius + penumbra radius

    # Calculate first contact (penumbra first touches Earth)
    # Search backward from maximum with a range of ~3.6 hours
    t_first_contact = _find_contact_time_besselian(
        jd_max, penumbral_limit, search_before=True, search_range=0.15
    )

    return t_first_contact


def calc_eclipse_second_contact_c2(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of second contact (C2) for a solar eclipse.

    Second contact (C2) is the moment when totality or annularity begins -
    when the Moon is fully inside (total eclipse) or outside (annular eclipse)
    the Sun's disk. At this instant, the umbral (total) or antumbral (annular)
    shadow first touches Earth's surface.

    This function uses Besselian elements to precisely calculate C2. The
    condition for C2 is when gamma (the distance of the shadow axis from
    Earth's center) equals 1 - |l2| (Earth radius minus umbral/antumbral cone
    radius), occurring before eclipse maximum.

    Note: C2 only exists for central eclipses (total or annular). For partial
    eclipses, this function returns 0.0.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function searches backward from
                this time to find C2.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of second contact C2. Returns 0.0 if:
        - The eclipse is not a central eclipse (total or annular)
        - C2 cannot be determined (gamma at maximum exceeds umbral limit)
        - The input time is not near a valid solar eclipse

    Algorithm:
        1. Calculate l2 (umbral/antumbral radius) at eclipse maximum
        2. Compute target gamma = 1 - |l2| (condition for umbra touching Earth)
        3. Check if gamma at maximum is less than umbral limit (central phase possible)
        4. Use binary search to find when gamma equals this target before maximum
        5. The search proceeds from (jd_max - search_range) to jd_max

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the gamma value converges to within 1e-8 Earth radii,
        which corresponds to approximately 0.06 km or 0.04 seconds of time.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, calc_eclipse_second_contact_c2
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate second contact
        >>> jd_c2 = calc_eclipse_second_contact_c2(jd_max)
        >>> print(f"Second contact C2: JD {jd_c2:.6f}")

    Note:
        - C2 is also known as "totality begins" or "annularity begins"
        - For partial eclipses, C2 does not exist (returns 0.0)
        - The duration of totality/annularity is (C3 - C2)
        - For local eclipse circumstances (at a specific observer location),
          use _sol_eclipse_when_loc_pythonic() which returns contact times in its result
        - The returned time is for the global eclipse (when umbra/antumbra first
          touches Earth anywhere), not for a specific observer location

    See Also:
        - calc_eclipse_first_contact_c1: Calculate first contact (eclipse begins)
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_besselian_l2: Calculate umbral/antumbral shadow radius

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Get l2 (umbral/antumbral radius) at maximum eclipse
    l2 = _calc_umbra_limit(jd_max)

    # Get gamma at maximum to check if central phase is possible
    gamma_max = _calc_gamma(jd_max)

    # The umbra/antumbra first touches Earth at EXTERIOR tangency:
    # shadow-axis distance = 1 + |l2| (consistent with tret[4] of
    # sol_eclipse_when_glob and with NASA's U1 convention).
    umbral_limit = 1.0 + abs(l2)

    # Check if a total/annular phase is possible anywhere on Earth
    if gamma_max >= umbral_limit:
        # No umbral contact - eclipse is partial only
        return 0.0

    # Calculate second contact (umbra/antumbra first touches Earth)
    # Search backward from maximum with a range of ~2.4 hours
    # (umbral contact typically occurs 0.5-1.5 hours before maximum)
    t_second_contact = _find_contact_time_besselian(
        jd_max, umbral_limit, search_before=True, search_range=0.10
    )

    return t_second_contact


def calc_eclipse_third_contact_c3(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of third contact (C3) for a solar eclipse.

    Third contact (C3) is the moment when totality or annularity ends -
    when the Moon's umbral (total) or antumbral (annular) shadow last touches
    Earth's surface. At this instant, the central phase of the eclipse ends
    and the partial phase resumes.

    This function uses Besselian elements to precisely calculate C3. The
    condition for C3 is when gamma (the distance of the shadow axis from
    Earth's center) equals 1 - |l2| (Earth radius minus umbral/antumbral cone
    radius), occurring after eclipse maximum.

    Note: C3 only exists for central eclipses (total or annular). For partial
    eclipses, this function returns 0.0.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function searches forward from
                this time to find C3.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of third contact C3. Returns 0.0 if:
        - The eclipse is not a central eclipse (total or annular)
        - C3 cannot be determined (gamma at maximum exceeds umbral limit)
        - The input time is not near a valid solar eclipse

    Algorithm:
        1. Calculate l2 (umbral/antumbral radius) at eclipse maximum
        2. Compute target gamma = 1 - |l2| (condition for umbra leaving Earth)
        3. Check if gamma at maximum is less than umbral limit (central phase possible)
        4. Use binary search to find when gamma equals this target after maximum
        5. The search proceeds from jd_max to (jd_max + search_range)

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the gamma value converges to within 1e-8 Earth radii,
        which corresponds to approximately 0.06 km or 0.04 seconds of time.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, calc_eclipse_third_contact_c3
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate third contact
        >>> jd_c3 = calc_eclipse_third_contact_c3(jd_max)
        >>> print(f"Third contact C3: JD {jd_c3:.6f}")

    Note:
        - C3 is also known as "totality ends" or "annularity ends"
        - For partial eclipses, C3 does not exist (returns 0.0)
        - The duration of totality/annularity is (C3 - C2)
        - For local eclipse circumstances (at a specific observer location),
          use _sol_eclipse_when_loc_pythonic() which returns contact times in its result
        - The returned time is for the global eclipse (when umbra/antumbra last
          leaves Earth anywhere), not for a specific observer location

    See Also:
        - calc_eclipse_first_contact_c1: Calculate first contact (eclipse begins)
        - calc_eclipse_second_contact_c2: Calculate second contact (totality begins)
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_besselian_l2: Calculate umbral/antumbral shadow radius

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Get l2 (umbral/antumbral radius) at maximum eclipse
    l2 = _calc_umbra_limit(jd_max)

    # Get gamma at maximum to check if central phase is possible
    gamma_max = _calc_gamma(jd_max)

    # For central eclipse, third contact occurs when gamma = 1 - |l2|
    # (umbral/antumbral shadow leaves Earth's limb from inside)
    # Exterior tangency, consistent with tret[5]/U4 (see C2 helper)
    umbral_limit = 1.0 + abs(l2)

    # Check if central phase is possible (gamma at max must be less than umbral limit)
    if gamma_max >= umbral_limit:
        # No central phase - eclipse is partial only
        return 0.0

    # Calculate third contact (umbra/antumbra last leaves Earth)
    # Search forward from maximum with a range of ~2.4 hours
    # (umbral contact typically occurs 0.5-1.5 hours after maximum)
    t_third_contact = _find_contact_time_besselian(
        jd_max, umbral_limit, search_before=False, search_range=0.10
    )

    return t_third_contact


def calc_eclipse_fourth_contact_c4(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of fourth external contact (C4) for a solar eclipse.

    Fourth contact (C4) is the moment when the Moon's disk completely separates
    from the Sun's disk externally, marking the end of a solar eclipse. At this
    instant, the penumbral shadow cone last touches Earth's surface.

    This function uses Besselian elements to precisely calculate C4. The
    condition for C4 is when gamma (the distance of the shadow axis from
    Earth's center) equals 1 + l1 (Earth radius plus penumbral cone radius),
    occurring after eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function searches forward from
                this time to find C4.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of fourth contact C4. Returns 0.0 if C4 cannot be
        determined (which would indicate the input time is not near a valid
        solar eclipse).

    Algorithm:
        1. Calculate l1 (penumbral radius) at eclipse maximum
        2. Compute target gamma = 1 + l1 (condition for penumbra leaving Earth)
        3. Use binary search to find when gamma equals this target after maximum
        4. The search proceeds from jd_max to (jd_max + search_range)

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the gamma value converges to within 1e-8 Earth radii,
        which corresponds to approximately 0.06 km or 0.04 seconds of time.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, calc_eclipse_fourth_contact_c4
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate fourth contact
        >>> jd_c4 = calc_eclipse_fourth_contact_c4(jd_max)
        >>> print(f"Fourth contact C4: JD {jd_c4:.6f}")

    Note:
        - C4 is also known as "fourth contact" or "partial phase ends"
        - C4 marks the end of the eclipse when the penumbra completely leaves Earth
        - For local eclipse circumstances (at a specific observer location),
          use _sol_eclipse_when_loc_pythonic() which returns contact times in its result
        - The returned time is for the global eclipse (when penumbra last
          leaves Earth anywhere), not for a specific observer location
        - The total duration of the eclipse is (C4 - C1)

    See Also:
        - calc_eclipse_first_contact_c1: Calculate first contact (eclipse begins)
        - calc_eclipse_second_contact_c2: Calculate second contact (totality begins)
        - calc_eclipse_third_contact_c3: Calculate third contact (totality ends)
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_besselian_l1: Calculate penumbral shadow radius

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Get l1 (penumbral radius) at maximum eclipse
    l1 = _calc_penumbra_limit(jd_max)

    # For global eclipse, fourth contact occurs when gamma = 1 + l1
    # (penumbral shadow leaves Earth's limb from inside)
    penumbral_limit = 1.0 + l1  # Earth radius + penumbra radius

    # Calculate fourth contact (penumbra completely leaves Earth)
    # Search forward from maximum with a range of ~3.6 hours
    # (penumbral contact typically occurs 1.5-3 hours after maximum)
    t_fourth_contact = _find_contact_time_besselian(
        jd_max, penumbral_limit, search_before=False, search_range=0.15
    )

    return t_fourth_contact


def _calc_lunar_eclipse_penumbral_separation(jd: float) -> float:
    """
    Calculate the separation between Moon's nearest limb and penumbral shadow edge.

    Returns the angular separation in degrees between the Moon's nearest limb
    and the penumbral shadow boundary. The penumbral contact occurs when this
    value equals zero (Moon's limb touches penumbra edge).

    Args:
        jd: Julian Day (UT)

    Returns:
        Separation in degrees. Positive means Moon is outside penumbra,
        negative means Moon's limb is inside penumbra.
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Extract coordinates
    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]
    moon_lat = moon_pos[1]
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Shadow center is opposite the Sun (anti-Sun point on ecliptic)
    shadow_lon = (sun_lon + 180.0) % 360.0

    # Calculate Moon's angular distance from shadow axis
    # The shadow axis passes through Earth opposite the Sun
    dlon = moon_lon - shadow_lon
    if dlon > 180:
        dlon -= 360
    if dlon < -180:
        dlon += 360

    # Angular distance from Moon center to shadow axis (degrees)
    # Using spherical geometry for accuracy
    moon_lat_rad = math.radians(moon_lat)
    distance_from_axis = math.sqrt(dlon**2 * math.cos(moon_lat_rad) ** 2 + moon_lat**2)

    # Constants for shadow geometry
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_AU = 4.2635e-5

    # Sun's angular semi-diameter at actual distance (in degrees)
    sun_semidiameter = (SUN_RADIUS_ARCSEC / 3600.0) / sun_dist

    # Moon's angular semi-diameter
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Earth's angular semi-diameter as seen from Moon (in degrees)
    earth_semidiameter = math.degrees(math.atan(EARTH_RADIUS_AU / moon_dist))

    # Penumbra radius at Moon's distance
    penumbra_radius = earth_semidiameter + sun_semidiameter

    # Apply Danjon atmospheric shadow enlargement (same factor used in
    # _calculate_lunar_eclipse_type_and_magnitude for consistency)
    _SHADOW_ENLARGEMENT = 1.0 + 1.0 / 85.0
    penumbra_radius *= _SHADOW_ENLARGEMENT

    # Separation: positive = Moon outside penumbra, negative = Moon inside
    # Contact occurs when Moon's nearest limb touches penumbra edge
    # At P1/P4: distance_from_axis - moon_semidiameter = penumbra_radius
    # So: separation = distance_from_axis - moon_semidiameter - penumbra_radius = 0
    separation = distance_from_axis - moon_semidiameter - penumbra_radius

    return separation


def _find_lunar_penumbral_contact_time(
    jd_max: float,
    search_before: bool,
    search_range: float = 0.20,
) -> float:
    """
    Find the time of lunar eclipse penumbral contact using binary search.

    Uses iterative refinement to find when Moon's limb touches the penumbral
    shadow boundary.

    Args:
        jd_max: Julian Day of eclipse maximum
        search_before: If True, search for P1 (before maximum), else P4 (after)
        search_range: Search range in days from maximum (default 0.20 = ~4.8 hours)

    Returns:
        Julian Day of the contact, or 0.0 if not found
    """
    # Set up search bounds
    if search_before:
        jd_low = jd_max - search_range
        jd_high = jd_max
    else:
        jd_low = jd_max
        jd_high = jd_max + search_range

    # Verify that a contact exists in the search range
    sep_low = _calc_lunar_eclipse_penumbral_separation(jd_low)
    sep_high = _calc_lunar_eclipse_penumbral_separation(jd_high)

    # For P1: sep_low should be positive (outside) and sep_high negative (inside)
    # For P4: sep_low should be negative (inside) and sep_high positive (outside)
    if search_before:
        # P1: looking for transition from positive to negative
        if sep_low < 0:
            # Moon already in penumbra at start of search - no P1 in range
            return 0.0
        if sep_high > 0:
            # Moon still outside penumbra at maximum - no eclipse
            return 0.0
    else:
        # P4: looking for transition from negative to positive
        if sep_low > 0:
            # Moon already outside penumbra at maximum
            return 0.0
        if sep_high < 0:
            # Moon still in penumbra at end of search
            return 0.0

    # Binary search for the contact time (when separation = 0)
    for _ in range(60):  # ~60 iterations gives sub-second precision
        jd_mid = (jd_low + jd_high) / 2
        sep_mid = _calc_lunar_eclipse_penumbral_separation(jd_mid)

        if search_before:
            # P1: separation goes from positive to negative
            if sep_mid > 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid
        else:
            # P4: separation goes from negative to positive
            if sep_mid < 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid

        # Check convergence (precision ~0.05 seconds)
        if jd_high - jd_low < 5e-7:
            break

    return (jd_low + jd_high) / 2


def calc_lunar_eclipse_penumbral_first_contact_p1(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of penumbral first contact (P1) for a lunar eclipse.

    P1 is the moment when the Moon's leading limb first enters Earth's
    penumbral shadow, marking the beginning of the lunar eclipse. At this
    instant, a very subtle shading begins on the Moon's eastern limb, though
    it is typically not visible to the naked eye until the Moon penetrates
    deeper into the penumbra.

    This function uses iterative refinement with geometric calculations to
    precisely determine P1. The condition for P1 is when the angular
    separation between the Moon's limb and the penumbral shadow edge equals
    zero, occurring before eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches backward from this time to find P1.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of penumbral first contact P1. Returns 0.0 if P1
        cannot be determined (which would indicate the input time is not
        near a valid lunar eclipse).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine penumbral shadow radius at Moon's distance
        3. Use binary search to find when Moon's limb touches penumbra edge
        4. The search proceeds from (jd_max - search_range) to jd_max

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_when_pythonic
        >>> from libephemeris import calc_lunar_eclipse_penumbral_first_contact_p1
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate penumbral first contact
        >>> jd_p1 = calc_lunar_eclipse_penumbral_first_contact_p1(jd_max)
        >>> print(f"Penumbral first contact P1: JD {jd_p1:.6f}")

    Note:
        - P1 is also known as "penumbral eclipse begins" or "first penumbral contact"
        - The penumbral phase (P1 to P4) is the total duration of the eclipse
        - For penumbral-only eclipses, P1 and P4 are the only contact points
        - The penumbral shading is typically visible only after the Moon has
          penetrated significantly into the penumbra (~70% or more)

    See Also:
        - calc_lunar_eclipse_penumbral_fourth_contact_p4: Calculate P4 (eclipse ends)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times
        - _lun_eclipse_when_loc_pythonic: Get local lunar eclipse circumstances

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate P1 (penumbral first contact)
    # Search backward from maximum with a range of ~4.8 hours
    # (penumbral contact typically occurs 2-3 hours before maximum)
    t_p1 = _find_lunar_penumbral_contact_time(
        jd_max, search_before=True, search_range=0.20
    )

    return t_p1


def calc_lunar_eclipse_penumbral_fourth_contact_p4(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of penumbral fourth contact (P4) for a lunar eclipse.

    P4 is the moment when the Moon's trailing limb completely exits Earth's
    penumbral shadow, marking the end of the lunar eclipse. At this instant,
    the last trace of penumbral shading leaves the Moon's western limb,
    although the shading would have been imperceptible for some time before.

    This function uses iterative refinement with geometric calculations to
    precisely determine P4. The condition for P4 is when the angular
    separation between the Moon's limb and the penumbral shadow edge equals
    zero, occurring after eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches forward from this time to find P4.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of penumbral fourth contact P4. Returns 0.0 if P4
        cannot be determined (which would indicate the input time is not
        near a valid lunar eclipse).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine penumbral shadow radius at Moon's distance
        3. Use binary search to find when Moon's limb exits penumbra edge
        4. The search proceeds from jd_max to (jd_max + search_range)

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_when_pythonic
        >>> from libephemeris import calc_lunar_eclipse_penumbral_fourth_contact_p4
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate penumbral fourth contact
        >>> jd_p4 = calc_lunar_eclipse_penumbral_fourth_contact_p4(jd_max)
        >>> print(f"Penumbral fourth contact P4: JD {jd_p4:.6f}")

    Note:
        - P4 is also known as "penumbral eclipse ends" or "last penumbral contact"
        - The total duration of the eclipse (P1 to P4) can be several hours
        - P4 marks the end of any visible or instrumental eclipse effects
        - The duration (P4 - P1) represents the complete penumbral phase

    See Also:
        - calc_lunar_eclipse_penumbral_first_contact_p1: Calculate P1 (eclipse begins)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times
        - _lun_eclipse_when_loc_pythonic: Get local lunar eclipse circumstances

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate P4 (penumbral fourth contact)
    # Search forward from maximum with a range of ~4.8 hours
    # (penumbral contact typically occurs 2-3 hours after maximum)
    t_p4 = _find_lunar_penumbral_contact_time(
        jd_max, search_before=False, search_range=0.20
    )

    return t_p4


def _calc_lunar_eclipse_umbral_outer_separation(jd: float) -> float:
    """
    Calculate the separation between Moon's nearest limb and umbral shadow edge.

    Returns the angular separation in degrees between the Moon's nearest limb
    and the umbral shadow boundary. The umbral contacts U1/U4 occur when this
    value equals zero (Moon's limb touches umbra edge).

    Args:
        jd: Julian Day (UT)

    Returns:
        Separation in degrees. Positive means Moon is outside umbra,
        negative means Moon's limb is inside umbra.
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Extract coordinates
    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]
    moon_lat = moon_pos[1]
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Shadow center is opposite the Sun (anti-Sun point on ecliptic)
    shadow_lon = (sun_lon + 180.0) % 360.0

    # Calculate Moon's angular distance from shadow axis
    dlon = moon_lon - shadow_lon
    if dlon > 180:
        dlon -= 360
    if dlon < -180:
        dlon += 360

    # Angular distance from Moon center to shadow axis (degrees)
    moon_lat_rad = math.radians(moon_lat)
    distance_from_axis = math.sqrt(dlon**2 * math.cos(moon_lat_rad) ** 2 + moon_lat**2)

    # Constants for shadow geometry
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_KM = 6378.137

    # Moon's angular semi-diameter
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Calculate umbra radius at Moon's distance
    sun_dist_km = sun_dist * 149597870.7
    moon_dist_km = moon_dist * 149597870.7
    sun_radius_km = (SUN_RADIUS_ARCSEC / 206265.0) * sun_dist_km

    # Umbra cone semi-angle
    umbra_cone_angle = math.atan((sun_radius_km - EARTH_RADIUS_KM) / sun_dist_km)
    # Umbra radius at Moon distance
    umbra_radius_km = EARTH_RADIUS_KM - moon_dist_km * math.tan(umbra_cone_angle)

    if umbra_radius_km > 0:
        umbra_radius = math.degrees(math.atan(umbra_radius_km / moon_dist_km))
    else:
        # Umbra doesn't reach Moon - no umbral eclipse possible
        umbra_radius = 0.0

    # Separation: positive = Moon outside umbra, negative = Moon inside
    # At U1/U4: distance_from_axis - moon_semidiameter = umbra_radius
    # So: separation = distance_from_axis - moon_semidiameter - umbra_radius = 0
    separation = distance_from_axis - moon_semidiameter - umbra_radius

    return separation


def _calc_lunar_eclipse_umbral_inner_separation(jd: float) -> float:
    """
    Calculate the separation for totality contacts (U2/U3).

    Returns the angular separation in degrees. At U2/U3, the Moon's trailing/leading
    limb touches the umbra edge (entire Moon inside umbra).

    Args:
        jd: Julian Day (UT)

    Returns:
        Separation in degrees. Positive means Moon is not fully inside umbra,
        negative means Moon is completely inside umbra.
    """
    # Get positions
    sun_pos, _ = calc_ut(jd, SUN, FLG_SPEED)
    moon_pos, _ = calc_ut(jd, MOON, FLG_SPEED)

    # Extract coordinates
    sun_lon = sun_pos[0]
    moon_lon = moon_pos[0]
    moon_lat = moon_pos[1]
    sun_dist = sun_pos[2]
    moon_dist = moon_pos[2]

    # Shadow center is opposite the Sun
    shadow_lon = (sun_lon + 180.0) % 360.0

    # Calculate Moon's angular distance from shadow axis
    dlon = moon_lon - shadow_lon
    if dlon > 180:
        dlon -= 360
    if dlon < -180:
        dlon += 360

    # Angular distance from Moon center to shadow axis
    moon_lat_rad = math.radians(moon_lat)
    distance_from_axis = math.sqrt(dlon**2 * math.cos(moon_lat_rad) ** 2 + moon_lat**2)

    # Constants for shadow geometry
    SUN_RADIUS_ARCSEC = 959.63
    EARTH_RADIUS_KM = 6378.137

    # Moon's angular semi-diameter
    moon_semidiameter = (932.56 / 3600.0) * (0.002569 / moon_dist)

    # Calculate umbra radius at Moon's distance
    sun_dist_km = sun_dist * 149597870.7
    moon_dist_km = moon_dist * 149597870.7
    sun_radius_km = (SUN_RADIUS_ARCSEC / 206265.0) * sun_dist_km

    # Umbra cone semi-angle
    umbra_cone_angle = math.atan((sun_radius_km - EARTH_RADIUS_KM) / sun_dist_km)
    # Umbra radius at Moon distance
    umbra_radius_km = EARTH_RADIUS_KM - moon_dist_km * math.tan(umbra_cone_angle)

    if umbra_radius_km > 0:
        umbra_radius = math.degrees(math.atan(umbra_radius_km / moon_dist_km))
    else:
        # No umbral eclipse possible
        umbra_radius = 0.0

    # For U2/U3 (totality contacts): Moon is fully inside umbra when
    # distance_from_axis + moon_semidiameter = umbra_radius
    # (farthest edge of Moon touches umbra edge)
    separation = distance_from_axis + moon_semidiameter - umbra_radius

    return separation


def _find_lunar_umbral_outer_contact_time(
    jd_max: float,
    search_before: bool,
    search_range: float = 0.15,
) -> float:
    """
    Find the time of lunar eclipse umbral outer contact (U1 or U4) using binary search.

    Args:
        jd_max: Julian Day of eclipse maximum
        search_before: If True, search for U1 (before maximum), else U4 (after)
        search_range: Search range in days from maximum

    Returns:
        Julian Day of the contact, or 0.0 if not found
    """
    # Set up search bounds
    if search_before:
        jd_low = jd_max - search_range
        jd_high = jd_max
    else:
        jd_low = jd_max
        jd_high = jd_max + search_range

    # Verify that a contact exists in the search range
    sep_low = _calc_lunar_eclipse_umbral_outer_separation(jd_low)
    sep_high = _calc_lunar_eclipse_umbral_outer_separation(jd_high)

    if search_before:
        # U1: looking for transition from positive to negative
        if sep_low < 0:
            # Moon already in umbra at start of search - expand search
            jd_low = jd_max - search_range * 1.5
            sep_low = _calc_lunar_eclipse_umbral_outer_separation(jd_low)
            if sep_low < 0:
                return 0.0
        if sep_high > 0:
            # Moon still outside umbra at maximum - no umbral eclipse
            return 0.0
    else:
        # U4: looking for transition from negative to positive
        if sep_low > 0:
            # Moon already outside umbra at maximum
            return 0.0
        if sep_high < 0:
            # Moon still in umbra at end of search - expand search
            jd_high = jd_max + search_range * 1.5
            sep_high = _calc_lunar_eclipse_umbral_outer_separation(jd_high)
            if sep_high < 0:
                return 0.0

    # Binary search for the contact time (when separation = 0)
    for _ in range(60):  # ~60 iterations gives sub-second precision
        jd_mid = (jd_low + jd_high) / 2
        sep_mid = _calc_lunar_eclipse_umbral_outer_separation(jd_mid)

        if search_before:
            # U1: separation goes from positive to negative
            if sep_mid > 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid
        else:
            # U4: separation goes from negative to positive
            if sep_mid < 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid

        # Check convergence (precision ~0.05 seconds)
        if jd_high - jd_low < 5e-7:
            break

    return (jd_low + jd_high) / 2


def _find_lunar_umbral_inner_contact_time(
    jd_max: float,
    search_before: bool,
    search_range: float = 0.10,
) -> float:
    """
    Find the time of lunar eclipse totality contact (U2 or U3) using binary search.

    Args:
        jd_max: Julian Day of eclipse maximum
        search_before: If True, search for U2 (totality begins), else U3 (totality ends)
        search_range: Search range in days from maximum

    Returns:
        Julian Day of the contact, or 0.0 if not found (no total eclipse)
    """
    # First check if this is a total eclipse by checking separation at maximum
    sep_max = _calc_lunar_eclipse_umbral_inner_separation(jd_max)
    if sep_max > 0:
        # Not a total eclipse - Moon never fully enters umbra
        return 0.0

    # Set up search bounds
    if search_before:
        jd_low = jd_max - search_range
        jd_high = jd_max
    else:
        jd_low = jd_max
        jd_high = jd_max + search_range

    # Verify that a contact exists in the search range
    sep_low = _calc_lunar_eclipse_umbral_inner_separation(jd_low)
    sep_high = _calc_lunar_eclipse_umbral_inner_separation(jd_high)

    if search_before:
        # U2: looking for transition from positive to negative
        if sep_low < 0:
            # Moon already fully in umbra at start - expand search
            jd_low = jd_max - search_range * 1.5
            sep_low = _calc_lunar_eclipse_umbral_inner_separation(jd_low)
            if sep_low < 0:
                jd_low = jd_max - search_range * 2.0
                sep_low = _calc_lunar_eclipse_umbral_inner_separation(jd_low)
                if sep_low < 0:
                    return 0.0
        if sep_high > 0:
            # Moon not fully in umbra at maximum - not a total eclipse
            return 0.0
    else:
        # U3: looking for transition from negative to positive
        if sep_low > 0:
            # Moon not fully in umbra at maximum - not a total eclipse
            return 0.0
        if sep_high < 0:
            # Moon still fully in umbra at end - expand search
            jd_high = jd_max + search_range * 1.5
            sep_high = _calc_lunar_eclipse_umbral_inner_separation(jd_high)
            if sep_high < 0:
                jd_high = jd_max + search_range * 2.0
                sep_high = _calc_lunar_eclipse_umbral_inner_separation(jd_high)
                if sep_high < 0:
                    return 0.0

    # Binary search for the contact time (when separation = 0)
    for _ in range(60):
        jd_mid = (jd_low + jd_high) / 2
        sep_mid = _calc_lunar_eclipse_umbral_inner_separation(jd_mid)

        if search_before:
            # U2: separation goes from positive to negative
            if sep_mid > 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid
        else:
            # U3: separation goes from negative to positive
            if sep_mid < 0:
                jd_low = jd_mid
            else:
                jd_high = jd_mid

        # Check convergence
        if jd_high - jd_low < 5e-7:
            break

    return (jd_low + jd_high) / 2


def calc_lunar_eclipse_umbral_first_contact_u1(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of umbral first contact (U1) for a lunar eclipse.

    U1 is the moment when the Moon's leading limb first enters Earth's
    umbral shadow, marking the beginning of the partial phase of the
    lunar eclipse. At this instant, the eastern edge of the Moon starts
    to darken noticeably as it enters the umbra.

    This function uses iterative refinement with geometric calculations to
    precisely determine U1. The condition for U1 is when the angular
    separation between the Moon's limb and the umbral shadow edge equals
    zero, occurring before eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches backward from this time to find U1.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of umbral first contact U1. Returns 0.0 if U1
        cannot be determined (which indicates either a penumbral-only eclipse
        or the input time is not near a valid lunar eclipse).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine umbral shadow radius at Moon's distance
        3. Use binary search to find when Moon's limb touches umbra edge
        4. The search proceeds from (jd_max - search_range) to jd_max

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_when_pythonic
        >>> from libephemeris import calc_lunar_eclipse_umbral_first_contact_u1
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate umbral first contact
        >>> jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max)
        >>> print(f"Umbral first contact U1: JD {jd_u1:.6f}")

    Note:
        - U1 is also known as "partial eclipse begins" or "first umbral contact"
        - U1 only occurs for partial and total lunar eclipses, not penumbral-only
        - The umbral phase (U1 to U4) is the most dramatic part of the eclipse
        - U1 marks when the Moon begins to visibly darken

    See Also:
        - calc_lunar_eclipse_umbral_fourth_contact_u4: Calculate U4 (partial ends)
        - calc_lunar_eclipse_umbral_second_contact_u2: Calculate U2 (totality begins)
        - calc_lunar_eclipse_umbral_third_contact_u3: Calculate U3 (totality ends)
        - calc_lunar_eclipse_penumbral_first_contact_p1: Calculate P1 (penumbral begins)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate U1 (umbral first contact)
    # Search backward from maximum with a range of ~3.6 hours
    # (umbral contact typically occurs 1-2 hours before maximum)
    t_u1 = _find_lunar_umbral_outer_contact_time(
        jd_max, search_before=True, search_range=0.15
    )

    return t_u1


def calc_lunar_eclipse_umbral_second_contact_u2(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of umbral second contact (U2) for a lunar eclipse.

    U2 is the moment when the Moon's trailing limb enters Earth's umbral
    shadow, marking the beginning of totality. At this instant, the entire
    Moon is within the umbra and the total phase of the eclipse begins.
    The Moon appears deeply red or copper-colored during totality.

    This function uses iterative refinement with geometric calculations to
    precisely determine U2. The condition for U2 is when the angular
    separation between the Moon's trailing limb and the umbral shadow edge
    equals zero, occurring before eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches backward from this time to find U2.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of umbral second contact U2. Returns 0.0 if U2
        cannot be determined (which indicates a partial-only or penumbral-only
        eclipse where totality does not occur).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine umbral shadow radius at Moon's distance
        3. Use binary search to find when Moon fully enters umbra
        4. The search proceeds from (jd_max - search_range) to jd_max

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_when_pythonic
        >>> from libephemeris import calc_lunar_eclipse_umbral_second_contact_u2
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate umbral second contact (totality begins)
        >>> jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max)
        >>> print(f"Totality begins U2: JD {jd_u2:.6f}")

    Note:
        - U2 only occurs for total lunar eclipses
        - Returns 0.0 for partial or penumbral eclipses
        - The duration (U3 - U2) is the length of totality
        - During totality, the Moon appears red due to refracted sunlight

    See Also:
        - calc_lunar_eclipse_umbral_third_contact_u3: Calculate U3 (totality ends)
        - calc_lunar_eclipse_umbral_first_contact_u1: Calculate U1 (partial begins)
        - calc_lunar_eclipse_umbral_fourth_contact_u4: Calculate U4 (partial ends)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate U2 (totality begins)
    # Search backward from maximum with a range of ~2.4 hours
    t_u2 = _find_lunar_umbral_inner_contact_time(
        jd_max, search_before=True, search_range=0.10
    )

    return t_u2


def calc_lunar_eclipse_umbral_third_contact_u3(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of umbral third contact (U3) for a lunar eclipse.

    U3 is the moment when the Moon's leading limb exits Earth's umbral
    shadow, marking the end of totality. At this instant, the western
    edge of the Moon begins to brighten as it emerges from the umbra,
    and the total phase of the eclipse ends.

    This function uses iterative refinement with geometric calculations to
    precisely determine U3. The condition for U3 is when the angular
    separation between the Moon's leading limb and the umbral shadow edge
    equals zero, occurring after eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches forward from this time to find U3.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of umbral third contact U3. Returns 0.0 if U3
        cannot be determined (which indicates a partial-only or penumbral-only
        eclipse where totality does not occur).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine umbral shadow radius at Moon's distance
        3. Use binary search to find when Moon begins exiting umbra
        4. The search proceeds from jd_max to (jd_max + search_range)

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_when_pythonic
        >>> from libephemeris import calc_lunar_eclipse_umbral_third_contact_u3
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate umbral third contact (totality ends)
        >>> jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max)
        >>> print(f"Totality ends U3: JD {jd_u3:.6f}")

    Note:
        - U3 only occurs for total lunar eclipses
        - Returns 0.0 for partial or penumbral eclipses
        - The duration (U3 - U2) is the length of totality
        - U3 marks when the Moon begins to brighten

    See Also:
        - calc_lunar_eclipse_umbral_second_contact_u2: Calculate U2 (totality begins)
        - calc_lunar_eclipse_umbral_first_contact_u1: Calculate U1 (partial begins)
        - calc_lunar_eclipse_umbral_fourth_contact_u4: Calculate U4 (partial ends)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate U3 (totality ends)
    # Search forward from maximum with a range of ~2.4 hours
    t_u3 = _find_lunar_umbral_inner_contact_time(
        jd_max, search_before=False, search_range=0.10
    )

    return t_u3


def calc_lunar_eclipse_umbral_fourth_contact_u4(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the time of umbral fourth contact (U4) for a lunar eclipse.

    U4 is the moment when the Moon's trailing limb completely exits Earth's
    umbral shadow, marking the end of the partial phase of the lunar eclipse.
    At this instant, the last portion of the Moon emerges from the umbra and
    returns to its normal brightness (though still in the penumbra).

    This function uses iterative refinement with geometric calculations to
    precisely determine U4. The condition for U4 is when the angular
    separation between the Moon's trailing limb and the umbral shadow edge
    equals zero, occurring after eclipse maximum.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function searches forward from this time to find U4.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Julian Day (UT) of umbral fourth contact U4. Returns 0.0 if U4
        cannot be determined (which indicates a penumbral-only eclipse
        or the input time is not near a valid lunar eclipse).

    Algorithm:
        1. Calculate Moon's position and Earth's shadow geometry
        2. Determine umbral shadow radius at Moon's distance
        3. Use binary search to find when Moon's limb exits umbra edge
        4. The search proceeds from jd_max to (jd_max + search_range)

    Precision:
        The calculation achieves timing precision better than 1 second by
        iterating until the time value converges to within ~0.05 seconds.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_when_pythonic
        >>> from libephemeris import calc_lunar_eclipse_umbral_fourth_contact_u4
        >>> from libephemeris import ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate umbral fourth contact
        >>> jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max)
        >>> print(f"Umbral fourth contact U4: JD {jd_u4:.6f}")

    Note:
        - U4 is also known as "partial eclipse ends" or "last umbral contact"
        - U4 only occurs for partial and total lunar eclipses, not penumbral-only
        - The duration (U4 - U1) is the length of the partial/umbral phase
        - After U4, the Moon is only in the penumbra until P4

    See Also:
        - calc_lunar_eclipse_umbral_first_contact_u1: Calculate U1 (partial begins)
        - calc_lunar_eclipse_umbral_second_contact_u2: Calculate U2 (totality begins)
        - calc_lunar_eclipse_umbral_third_contact_u3: Calculate U3 (totality ends)
        - calc_lunar_eclipse_penumbral_fourth_contact_p4: Calculate P4 (eclipse ends)
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate U4 (umbral fourth contact)
    # Search forward from maximum with a range of ~3.6 hours
    t_u4 = _find_lunar_umbral_outer_contact_time(
        jd_max, search_before=False, search_range=0.15
    )

    return t_u4


def calc_solar_eclipse_duration(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the duration of totality or annularity for a solar eclipse.

    For central solar eclipses (total or annular), this function calculates
    the duration of the central phase - the time between second contact (C2)
    and third contact (C3). This represents the length of time the umbral
    (total) or antumbral (annular) shadow touches Earth's surface.

    For non-central eclipses (partial only), the function returns 0.0 since
    there is no totality or annularity phase.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _sol_eclipse_when_glob_pythonic()
                or _sol_eclipse_when_loc_pythonic(). The function calculates C2 and C3
                relative to this time.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Duration of totality/annularity in minutes. Returns 0.0 if:
        - The eclipse is not a central eclipse (total or annular)
        - C2 or C3 cannot be determined
        - The input time is not near a valid solar eclipse

    Algorithm:
        1. Calculate second contact C2 (when central phase begins)
        2. Calculate third contact C3 (when central phase ends)
        3. Return (C3 - C2) converted from days to minutes

    Precision:
        The calculation achieves timing precision better than 1 second,
        resulting in duration precision of approximately ±0.03 minutes (±2 seconds).

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic
        >>> from libephemeris import calc_solar_eclipse_duration, ECL_TOTAL
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate duration of totality
        >>> duration = calc_solar_eclipse_duration(jd_max)
        >>> print(f"Duration of totality: {duration:.2f} minutes")

    Note:
        - This is the global duration (time umbra/antumbra is on Earth anywhere)
        - Local duration at a specific observer location is typically shorter
        - For local eclipse duration, use _sol_eclipse_when_loc_pythonic() which provides
          contact times for the specific location
        - Total solar eclipses can have central durations up to ~7.5 minutes
        - Annular eclipses can have central durations up to ~12 minutes
        - The global duration is always longer than any local duration

    See Also:
        - calc_eclipse_second_contact_c2: Calculate C2 (central phase begins)
        - calc_eclipse_third_contact_c3: Calculate C3 (central phase ends)
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and all phase times
        - _sol_eclipse_when_loc_pythonic: Get local eclipse circumstances including contacts
        - calc_lunar_eclipse_total_duration: Duration of lunar eclipse totality
        - calc_lunar_eclipse_umbral_duration: Duration of lunar eclipse umbral phase

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate C2 (second contact - central phase begins)
    jd_c2 = calc_eclipse_second_contact_c2(jd_max, flags=flags)

    # Check if central phase exists
    if jd_c2 == 0.0:
        return 0.0

    # Calculate C3 (third contact - central phase ends)
    jd_c3 = calc_eclipse_third_contact_c3(jd_max, flags=flags)

    # Check if C3 was found
    if jd_c3 == 0.0:
        return 0.0

    # Calculate duration in minutes (C3 - C2) * 24 hours * 60 minutes
    duration_days = jd_c3 - jd_c2
    duration_minutes = duration_days * 24.0 * 60.0

    return duration_minutes


def calc_lunar_eclipse_total_duration(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the duration of totality for a total lunar eclipse.

    For total lunar eclipses, this function calculates the duration of
    totality - the time between umbral second contact (U2) and umbral
    third contact (U3). This represents the length of time the entire
    Moon is within Earth's umbral shadow.

    For partial or penumbral lunar eclipses, the function returns 0.0
    since totality does not occur.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function calculates U2 and U3 relative to this time.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Duration of totality in minutes. Returns 0.0 if:
        - The eclipse is not a total lunar eclipse
        - U2 or U3 cannot be determined
        - The input time is not near a valid lunar eclipse

    Algorithm:
        1. Calculate umbral second contact U2 (when totality begins)
        2. Calculate umbral third contact U3 (when totality ends)
        3. Return (U3 - U2) converted from days to minutes

    Precision:
        The calculation achieves timing precision better than 1 second,
        resulting in duration precision of approximately ±0.03 minutes (±2 seconds).

    Example:
        >>> from libephemeris import julday, _lun_eclipse_when_pythonic
        >>> from libephemeris import calc_lunar_eclipse_total_duration, ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate duration of totality
        >>> duration = calc_lunar_eclipse_total_duration(jd_max)
        >>> print(f"Duration of totality: {duration:.2f} minutes")

    Note:
        - U2 marks when the Moon fully enters the umbra (totality begins)
        - U3 marks when the Moon begins exiting the umbra (totality ends)
        - During totality, the Moon appears red/copper due to refracted sunlight
        - Total lunar eclipse durations can range from about 14 to 107 minutes
        - Unlike solar eclipses, lunar eclipse totality is visible from the
          entire night side of Earth

    See Also:
        - calc_lunar_eclipse_umbral_second_contact_u2: Calculate U2 (totality begins)
        - calc_lunar_eclipse_umbral_third_contact_u3: Calculate U3 (totality ends)
        - calc_lunar_eclipse_umbral_duration: Duration of umbral/partial phase
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate U2 (umbral second contact - totality begins)
    jd_u2 = calc_lunar_eclipse_umbral_second_contact_u2(jd_max, flags=flags)

    # Check if totality exists
    if jd_u2 == 0.0:
        return 0.0

    # Calculate U3 (umbral third contact - totality ends)
    jd_u3 = calc_lunar_eclipse_umbral_third_contact_u3(jd_max, flags=flags)

    # Check if U3 was found
    if jd_u3 == 0.0:
        return 0.0

    # Calculate duration in minutes (U3 - U2) * 24 hours * 60 minutes
    duration_days = jd_u3 - jd_u2
    duration_minutes = duration_days * 24.0 * 60.0

    return duration_minutes


def calc_lunar_eclipse_umbral_duration(
    jd_max: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the duration of the umbral (partial) phase for a lunar eclipse.

    For lunar eclipses that enter Earth's umbral shadow (total or partial
    eclipses), this function calculates the duration of the umbral phase -
    the time between umbral first contact (U1) and umbral fourth contact (U4).
    This represents the length of time any part of the Moon is within the
    Earth's umbral shadow.

    For penumbral-only lunar eclipses, the function returns 0.0 since
    the Moon does not enter the umbral shadow.

    Args:
        jd_max: Julian Day (UT) of eclipse maximum. This should be the time
                of greatest eclipse, which can be obtained from _lun_eclipse_when_pythonic().
                The function calculates U1 and U4 relative to this time.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        Duration of umbral phase in minutes. Returns 0.0 if:
        - The eclipse is a penumbral-only eclipse
        - U1 or U4 cannot be determined
        - The input time is not near a valid lunar eclipse

    Algorithm:
        1. Calculate umbral first contact U1 (when partial phase begins)
        2. Calculate umbral fourth contact U4 (when partial phase ends)
        3. Return (U4 - U1) converted from days to minutes

    Precision:
        The calculation achieves timing precision better than 1 second,
        resulting in duration precision of approximately ±0.03 minutes (±2 seconds).

    Example:
        >>> from libephemeris import julday, _lun_eclipse_when_pythonic
        >>> from libephemeris import calc_lunar_eclipse_umbral_duration, ECL_TOTAL
        >>> # Find the November 8, 2022 total lunar eclipse
        >>> jd_start = julday(2022, 10, 1, 0.0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate duration of umbral phase
        >>> duration = calc_lunar_eclipse_umbral_duration(jd_max)
        >>> print(f"Duration of umbral phase: {duration:.2f} minutes")

    Note:
        - U1 marks when the Moon first enters the umbra (partial phase begins)
        - U4 marks when the Moon fully exits the umbra (partial phase ends)
        - This duration includes totality (if it occurs) plus both partial phases
        - Umbral phase durations can range from about 24 to 236 minutes
        - Unlike solar eclipses, lunar eclipse umbral phase is visible from
          the entire night side of Earth

    See Also:
        - calc_lunar_eclipse_umbral_first_contact_u1: Calculate U1 (partial begins)
        - calc_lunar_eclipse_umbral_fourth_contact_u4: Calculate U4 (partial ends)
        - calc_lunar_eclipse_total_duration: Duration of totality phase only
        - _lun_eclipse_when_pythonic: Find next lunar eclipse and all phase times

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Lunar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
    """
    # Calculate U1 (umbral first contact - partial phase begins)
    jd_u1 = calc_lunar_eclipse_umbral_first_contact_u1(jd_max, flags=flags)

    # Check if umbral phase exists
    if jd_u1 == 0.0:
        return 0.0

    # Calculate U4 (umbral fourth contact - partial phase ends)
    jd_u4 = calc_lunar_eclipse_umbral_fourth_contact_u4(jd_max, flags=flags)

    # Check if U4 was found
    if jd_u4 == 0.0:
        return 0.0

    # Calculate duration in minutes (U4 - U1) * 24 hours * 60 minutes
    duration_days = jd_u4 - jd_u1
    duration_minutes = duration_days * 24.0 * 60.0

    return duration_minutes


# =============================================================================
# ECLIPSE PATH WIDTH CALCULATION
# =============================================================================


def calc_eclipse_path_width(
    jd: float,
    lat: float | None = None,
    lon: float | None = None,
    flags: int = FLG_SWIEPH,
) -> float:
    """Calculate the width of the central eclipse path at a location.

    Wrapper around :func:`_calc_eclipse_path_width_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _calc_eclipse_path_width_impl, jd, lat, lon, flags,
    )


def _calc_eclipse_path_width_impl(
    jd: float,
    lat: float | None = None,
    lon: float | None = None,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> float:
    """
    Calculate the width of the path of totality or annularity for a central solar eclipse.

    For central solar eclipses (total, annular, or hybrid), this function calculates
    the width of the umbral (total) or antumbral (annular) shadow path on Earth's
    surface at a specific time or geographic location.

    The path width is the perpendicular distance across the central line within
    which observers will see totality or annularity. Outside this band, the
    eclipse appears partial.

    Args:
        jd: Julian Day (UT) at which to calculate the path width. This should be
            a time during a central solar eclipse when the shadow is on Earth.
        lat: Observer latitude in degrees (positive = North, negative = South).
             If provided along with lon, calculates path width at this location.
             If None, calculates path width at the central line for time jd.
        lon: Observer longitude in degrees (positive = East, negative = West).
             If provided along with lat, calculates path width at this location.
             If None, calculates path width at the central line for time jd.
        flags: Calculation flags (FLG_SWIEPH by default). Controls which
               ephemeris to use for the underlying calculations.

    Returns:
        The path width in kilometers. Returns 0.0 if:
        - The eclipse is not a central eclipse (partial only)
        - The shadow does not touch Earth at the given time
        - No valid eclipse is occurring at the given time
        - The Sun is below the horizon at the specified location

    Algorithm:
        1. Calculate the Besselian elements l2 (umbral/antumbral cone radius)
        2. Get the Sun altitude at the central line (or specified location)
        3. Calculate the geometric shadow width using the formula:
           width = 2 * |l2| * Earth_radius / sin(sun_altitude)
        4. Apply corrections for Earth's curvature and observer elevation

    Physical interpretation:
        - The umbral/antumbral shadow is a cone with its vertex near the Moon
        - Where this cone intersects Earth's surface, it creates an elliptical
          shadow footprint
        - The path width is the minor axis of this ellipse perpendicular to
          the direction of shadow motion
        - Path width is minimum when the Sun is at zenith and increases as
          the Sun is lower in the sky

    Precision:
        Accurate to approximately ±1 km for typical eclipses when the Sun is
        more than 10° above the horizon. Less accurate near sunrise/sunset
        due to grazing geometry.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, ECL_TOTAL
        >>> from libephemeris import calc_eclipse_path_width
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start, eclipse_type=ECL_TOTAL)
        >>> jd_max = times[0]
        >>> # Calculate path width at eclipse maximum
        >>> width = calc_eclipse_path_width(jd_max)
        >>> print(f"Path width at maximum: {width:.1f} km")

        >>> # Calculate path width at a specific location
        >>> dallas_lat, dallas_lon = 32.7767, -96.7970
        >>> width_dallas = calc_eclipse_path_width(jd_max, lat=dallas_lat, lon=dallas_lon)
        >>> print(f"Path width at Dallas: {width_dallas:.1f} km")

    Note:
        - Path width for total eclipses typically ranges from 0 to ~270 km
        - Path width for annular eclipses typically ranges from 0 to ~380 km
        - Very narrow eclipses (< 10 km) indicate the shadow axis barely
          grazes Earth's surface
        - At the limits of the central path (sunrise/sunset), the path width
          can exceed 1000 km due to grazing geometry
        - The calculation assumes a spherical Earth with WGS84 equatorial radius

    See Also:
        - sol_eclipse_where: Find central line coordinates and path limits
        - calc_besselian_l2: Calculate the umbral/antumbral cone radius
        - _sol_eclipse_when_glob_pythonic: Find next solar eclipse and phase times
        - _sol_eclipse_how_pythonic: Get eclipse circumstances at a specific location

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac (2013), Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    # reader is provided by the caller (None forces Skyfield path)

    # Constants
    AU_TO_KM = 149597870.7
    SUN_RADIUS_KM = 696340.0
    MOON_RADIUS_KM = 1737.4

    # Determine the location for calculation
    if lat is not None and lon is not None:
        # Use specified location
        calc_lat = lat
        calc_lon = lon
    else:
        # Find central line location at this time using sol_eclipse_where
        ecl_type, geopos_result, attr = sol_eclipse_where(jd, flags)
        calc_lon = geopos_result[0]
        calc_lat = geopos_result[1]

        # Check if this is a central eclipse
        if not (ecl_type & ECL_CENTRAL):
            return 0.0

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        _gp = (calc_lon, calc_lat, 0.0)

        try:
            jd_tt = jd + deltat(jd)
            sun_topo = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_topo = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        sun_az_val, sun_alt_true, sun_alt_app = azalt(jd, ECL2HOR, _gp, 0, 0, sun_topo[:3])
        sun_altitude = sun_alt_true

        if sun_altitude <= 0:
            return 0.0

        sun_dist_au = sun_topo[2]
        moon_dist_au = moon_topo[2]
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()
        t = ts.ut1_jd(jd)

        # Get celestial bodies
        earth = eph["earth"]
        sun = eph["sun"]
        moon = eph["moon"]

        # Create observer at the calculation location
        try:
            observer = wgs84.latlon(calc_lat, calc_lon, 0.0)
            observer_at = earth + observer
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        # Get Sun and Moon apparent positions from the observer
        try:
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        # Get Sun altitude
        sun_alt, sun_az, _ = sun_app.altaz()
        sun_altitude = sun_alt.degrees

        # If Sun is below horizon, no path width
        if sun_altitude <= 0:
            return 0.0

        # Get distances
        sun_dist_au = sun_app.distance().au
        moon_dist_au = moon_app.distance().au

    sun_dist_km = sun_dist_au * AU_TO_KM
    moon_dist_km = moon_dist_au * AU_TO_KM

    # Calculate angular radii for eclipse type determination
    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au  # degrees
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)  # degrees

    moon_sun_ratio = moon_angular_radius / sun_angular_radius

    # Check if Sun and Moon are close enough for central eclipse
    if reader is not None:
        from .utils import angular_separation as _ang_sep_pw

        separation = _ang_sep_pw(sun_topo[0], sun_topo[1], moon_topo[0], moon_topo[1])
    else:
        separation = sun_app.separation_from(moon_app).degrees
    diff_radii = abs(moon_angular_radius - sun_angular_radius)

    if separation > diff_radii:
        # Not a central eclipse at this location/time
        return 0.0

    # Calculate the umbral/antumbral cone geometry
    # The umbra cone semi-angle (alpha) is the angle at the Moon
    # between the cone edge and the shadow axis
    alpha = math.atan((SUN_RADIUS_KM - MOON_RADIUS_KM) / sun_dist_km)

    # Calculate the shadow radius at Earth's surface
    if moon_sun_ratio >= 1.0:
        # Total eclipse - umbra reaches Earth
        # The umbra radius decreases along the shadow axis
        umbra_radius_km = MOON_RADIUS_KM - moon_dist_km * math.tan(alpha)
        if umbra_radius_km <= 0:
            return 0.0  # Umbra vertex is before Earth
    else:
        # Annular eclipse - antumbra
        # The antumbra radius increases along the shadow axis
        umbra_radius_km = moon_dist_km * math.tan(alpha) - MOON_RADIUS_KM
        if umbra_radius_km < 0:
            umbra_radius_km = abs(umbra_radius_km)

    # Calculate the path width
    # The shadow hits Earth at an angle determined by the Sun's altitude
    # Path width = 2 * umbra_radius / sin(sun_altitude)
    sin_alt = math.sin(math.radians(sun_altitude))

    # Protect against very small values near horizon (grazing eclipses)
    if sin_alt < 0.05:
        sin_alt = 0.05  # Limit to ~3 degrees to avoid unrealistic values

    path_width_km = 2.0 * umbra_radius_km / sin_alt

    # Apply a reasonable upper limit for grazing eclipses
    path_width_km = min(path_width_km, 1500.0)

    return max(0.0, path_width_km)


# Alias for reference API compatibility
calc_eclipse_path_width = calc_eclipse_path_width


# =============================================================================
# SAROS AND INEX SERIES CALCULATION
# =============================================================================

# Saros cycle period in days (223 synodic months)
SAROS_CYCLE_DAYS = 6585.3211  # More precise: 6585.3211 days

# Inex cycle period in days (358 synodic months, ~28.945 years)
# The Inex cycle relates eclipses that occur at nearly the same longitude
# but shifted by one lunar node (ascending to descending or vice versa).
INEX_CYCLE_DAYS = 10571.9509  # 358 * 29.530588853 = 10571.9509 days

# Reference solar eclipses with known Saros series and member numbers
# Comprehensive coverage of all ~25 active solar Saros series (2015-2028)
# Data from NASA's Five Millennium Canon of Solar Eclipses (Espenak & Meeus)
# Format: (JD of eclipse maximum, Saros series number, member number)
_SOLAR_SAROS_REFERENCES = (
    (2457101.907, 120, 61),
    (2457278.788, 125, 54),
    (2457456.581, 130, 52),
    (2457632.880, 135, 39),
    (2457811.120, 140, 29),
    (2457987.268, 145, 22),
    (2458165.369, 150, 17),
    (2458312.626, 117, 69),
    (2458489.571, 122, 58),
    (2458667.308, 127, 58),
    (2458843.721, 132, 46),
    (2459021.778, 137, 36),
    (2459198.176, 142, 23),
    (2459375.946, 147, 23),
    (2459552.815, 152, 13),
    (2459700.362, 119, 66),
    (2459877.958, 124, 55),
    (2460054.678, 129, 52),
    (2460232.250, 134, 44),
    (2460409.262, 139, 30),
    (2460586.281, 144, 17),
    (2460763.950, 149, 21),
    (2460940.321, 154, 7),
    (2461089.008, 121, 61),
    (2461265.240, 126, 48),
)

# Reference lunar eclipses with known Saros series and member numbers
# Comprehensive coverage of all ~30 active lunar Saros series (2015-2030)
# Data from NASA's Five Millennium Canon of Lunar Eclipses (Espenak & Meeus)
# Format: (JD of eclipse maximum, Saros series number, member number)
_LUNAR_SAROS_REFERENCES = (
    (2457117.000, 132, 30),
    (2457293.616, 137, 26),
    (2457470.991, 142, 18),
    (2457648.288, 147, 8),
    (2457795.530, 114, 59),
    (2457973.264, 119, 61),
    (2458150.062, 124, 49),
    (2458327.348, 129, 38),
    (2458504.717, 134, 27),
    (2458681.396, 139, 21),
    (2458859.299, 144, 16),
    (2459006.309, 111, 67),
    (2459183.905, 116, 58),
    (2459360.971, 121, 55),
    (2459537.877, 126, 45),
    (2459715.675, 131, 34),
    (2459891.958, 136, 20),
    (2460070.224, 141, 24),
    (2460246.343, 146, 11),
    (2460394.801, 113, 64),
    (2460571.614, 118, 52),
    (2460748.791, 123, 53),
    (2460926.258, 128, 41),
    (2461102.982, 133, 27),
    (2461280.676, 138, 29),
    (2461457.467, 143, 18),
    (2461605.169, 110, 72),
    (2461782.676, 115, 58),
    (2461959.264, 120, 58),
    (2462137.203, 125, 49),
)


def _get_saros_info(
    jd_eclipse: float,
    eclipse_type: str = "solar",
) -> Tuple[float, float]:
    """Return (saros_series, saros_member) for an eclipse as floats.

    Internal helper used to populate attr[9] and attr[10] in eclipse
    attribute tuples. Returns (0.0, 0.0) if no match is found.

    Args:
        jd_eclipse: Julian Day (UT) of the eclipse maximum.
        eclipse_type: "solar" or "lunar".

    Returns:
        Tuple of (saros_series_number, saros_member_number) as floats.
    """
    if eclipse_type.lower() == "solar":
        references = _SOLAR_SAROS_REFERENCES
    elif eclipse_type.lower() == "lunar":
        references = _LUNAR_SAROS_REFERENCES
    else:
        return (0.0, 0.0)

    best_series = 0
    best_member = 0
    best_residual = float("inf")

    for ref_jd, ref_saros, ref_member in references:
        delta = jd_eclipse - ref_jd
        n_cycles = round(delta / SAROS_CYCLE_DAYS)
        residual = abs(delta - n_cycles * SAROS_CYCLE_DAYS)

        if residual < best_residual:
            best_residual = residual
            best_series = ref_saros
            best_member = ref_member + n_cycles

    # Sanity check: residual should be < ~1 day for a valid eclipse
    if best_residual > 2.0:
        return (0.0, 0.0)

    return (float(best_series), float(best_member))


def get_saros_number(
    jd_eclipse: float,
    eclipse_type: str = "solar",
) -> int:
    """
    Determine the Saros series number for an eclipse.

    The Saros cycle is a period of approximately 6585.3211 days (18 years,
    11 days, 8 hours) after which the Sun, Moon, and lunar nodes return
    to nearly the same relative geometry, causing similar eclipses to recur.

    Each Saros series is numbered sequentially. Solar eclipses have series
    numbered approximately 1-180 (though not all are active at any given time),
    and lunar eclipses similarly have their own numbering.

    Args:
        jd_eclipse: Julian Day (UT) of the eclipse maximum. This should be
                    obtained from _sol_eclipse_when_glob_pythonic() or _lun_eclipse_when_pythonic().
        eclipse_type: Type of eclipse - either "solar" or "lunar".
                      Defaults to "solar".

    Returns:
        The Saros series number (integer). Solar Saros numbers typically range
        from about 1 to 180, and lunar Saros numbers similarly.

    Algorithm:
        1. For the given eclipse JD, calculate the difference from known
           reference eclipses with documented Saros series numbers.
        2. Determine how many Saros cycles separate the input eclipse from
           the reference.
        3. Each Saros cycle advances the series number by 1 (earlier eclipses
           have lower numbers in the same series, later have higher).

    Precision:
        The function is accurate for eclipses within several centuries of
        the reference eclipses. Very distant eclipses (more than ~500 years
        from reference) may have small errors due to the slight variation
        in the Saros period.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, get_saros_number
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 3, 1, 0.0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start)
        >>> jd_max = times[0]
        >>> saros = get_saros_number(jd_max, "solar")
        >>> print(f"Saros series: {saros}")  # Should print 139

    Note:
        - Each Saros series begins with partial eclipses near one pole,
          progresses through central (total/annular) eclipses, and ends
          with partial eclipses near the other pole.
        - A typical Saros series contains 70-80 eclipses over ~1200-1400 years.
        - The same Saros number for solar and lunar eclipses are different
          series (they use separate numbering systems).

    References:
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - van Gent, R.H. "A Catalogue of Eclipse Cycles"
    """
    if eclipse_type.lower() not in ("solar", "lunar"):
        raise ValueError(
            f"eclipse_type must be 'solar' or 'lunar', got '{eclipse_type}'"
        )
    saros_series, _ = _get_saros_info(jd_eclipse, eclipse_type)
    return int(saros_series)


# =============================================================================
# INEX SERIES CALCULATION
# =============================================================================

# Reference solar eclipses with known Inex series numbers
# The Inex cycle connects eclipses at the same solar longitude
# but at opposite lunar nodes (ascending <-> descending).
# Series numbers from van Gent's eclipse cycle catalog.
# Format: (JD of eclipse maximum, Inex series number)
_SOLAR_INEX_REFERENCES = (
    # Inex 50 - includes 2024 total solar eclipse
    (2460409.296, 50),  # 2024-Apr-08 total solar eclipse (Inex 50)
    # Inex 49 - includes 2017 total solar eclipse
    (2457987.767, 49),  # 2017-Aug-21 total solar eclipse (Inex 49)
    # Inex 51 - includes 2020 total solar eclipse
    (2459203.604, 51),  # 2020-Dec-14 total solar eclipse (Inex 51)
    # Inex 48 - includes 2019 total solar eclipse
    (2458675.604, 48),  # 2019-Jul-02 total solar eclipse (Inex 48)
    # Inex 52 - includes 2021 annular solar eclipse
    (2459369.971, 52),  # 2021-Jun-10 annular solar eclipse (Inex 52)
)

# Reference lunar eclipses with known Inex series numbers
# Format: (JD of eclipse maximum, Inex series number)
_LUNAR_INEX_REFERENCES = (
    # Inex 50 - includes 2022 total lunar eclipse
    (2459891.459, 50),  # 2022-Nov-08 total lunar eclipse (Inex 50)
    # Inex 49 - includes 2022 total lunar eclipse
    (2459695.604, 49),  # 2022-May-16 total lunar eclipse (Inex 49)
    # Inex 48 - includes 2021 total lunar eclipse
    (2459356.917, 48),  # 2021-May-26 total lunar eclipse (Inex 48)
    # Inex 51 - includes 2018 total lunar eclipse
    (2458310.835, 51),  # 2018-Jul-27 total lunar eclipse (Inex 51)
    # Inex 47 - includes 2019 total lunar eclipse
    (2458497.459, 47),  # 2019-Jan-21 total lunar eclipse (Inex 47)
)


def get_inex_number(
    jd_eclipse: float,
    eclipse_type: str = "solar",
) -> int:
    """
    Determine the Inex series number for an eclipse.

    The Inex cycle is a period of approximately 10571.95 days (358 synodic
    months, ~28.945 years) after which eclipses occur at the same solar
    longitude but at the opposite lunar node.

    While the Saros cycle connects eclipses with similar geometry (same node,
    similar latitude), the Inex cycle connects eclipses at the same longitude
    but at different nodes. Together, Saros and Inex form a two-dimensional
    numbering system for eclipse series.

    Args:
        jd_eclipse: Julian Day (UT) of the eclipse maximum. This should be
                    obtained from _sol_eclipse_when_glob_pythonic() or _lun_eclipse_when_pythonic().
        eclipse_type: Type of eclipse - either "solar" or "lunar".
                      Defaults to "solar".

    Returns:
        The Inex series number (integer).

    Algorithm:
        1. For the given eclipse JD, calculate the difference from known
           reference eclipses with documented Inex series numbers.
        2. Determine how many Inex cycles separate the input eclipse from
           the reference.
        3. Find the reference that gives the smallest residual (closest
           match to an integer number of cycles).

    Precision:
        The function is accurate for eclipses within several centuries of
        the reference eclipses. Very distant eclipses may have small errors
        due to slight variations in the Inex period over long time spans.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, get_inex_number
        >>> # Find the April 8, 2024 total solar eclipse
        >>> jd_start = julday(2024, 3, 1, 0.0)
        >>> ecl_type, times = _sol_eclipse_when_glob_pythonic(jd_start)
        >>> jd_max = times[0]
        >>> inex = get_inex_number(jd_max, "solar")
        >>> print(f"Inex series: {inex}")  # Should print 50

    Note:
        - The Inex series number combined with the Saros series number
          uniquely identifies each eclipse in the long-term eclipse cycle.
        - Each Inex series contains eclipses that occur at the same solar
          longitude but alternate between ascending and descending nodes.
        - Unlike Saros, which relates eclipses of similar appearance, Inex
          relates eclipses at the same celestial longitude.

    References:
        - van Gent, R.H. "A Catalogue of Eclipse Cycles"
        - Meeus, J. "Mathematical Astronomy Morsels"
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
    """
    if eclipse_type.lower() == "solar":
        references = _SOLAR_INEX_REFERENCES
    elif eclipse_type.lower() == "lunar":
        references = _LUNAR_INEX_REFERENCES
    else:
        raise ValueError(
            f"eclipse_type must be 'solar' or 'lunar', got '{eclipse_type}'"
        )

    # Use the first reference eclipse as our anchor point
    ref_jd, ref_inex = references[0]

    # Calculate how many Inex cycles separate the input from the reference
    delta_days = jd_eclipse - ref_jd

    # Number of complete Inex cycles (can be positive or negative)
    cycles = round(delta_days / INEX_CYCLE_DAYS)

    # Find which reference gives the closest match
    best_match_series = ref_inex
    best_match_residual = abs(delta_days - cycles * INEX_CYCLE_DAYS)

    for ref_jd_test, ref_inex_test in references:
        delta = jd_eclipse - ref_jd_test
        cycles_test = round(delta / INEX_CYCLE_DAYS)
        residual = abs(delta - cycles_test * INEX_CYCLE_DAYS)

        if residual < best_match_residual:
            best_match_residual = residual
            best_match_series = ref_inex_test

    return best_match_series


def calc_eclipse_central_line(
    jd_start: float,
    jd_end: float,
    step_minutes: float = 1.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], Tuple[float, ...]]:
    """
    Calculate the geographic coordinates of points along the central line of a solar eclipse.

    The central line is the path traced on Earth's surface by the intersection of
    the Moon's shadow axis with the Earth. Along this line, observers see the
    eclipse at its maximum (either total or annular, depending on eclipse type).

    This function computes a series of (latitude, longitude) points along the
    central line at specified time intervals using the same algorithm as
    _sol_eclipse_where_pythonic() for consistent, accurate results.

    Args:
        jd_start: Julian Day (UT) to start calculating the central line.
                  Should be during a solar eclipse (ideally from _sol_eclipse_when_glob_pythonic).
        jd_end: Julian Day (UT) to end the calculation.
                Should be after jd_start and during the same eclipse.
        step_minutes: Time step in minutes between calculated points (default 1.0).
                      Smaller values give a more detailed path but take longer.
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing three tuples:
            - times: Tuple of Julian Day (UT) values for each point
            - latitudes: Tuple of geographic latitudes in degrees (North positive)
            - longitudes: Tuple of geographic longitudes in degrees (East positive)

        Points where the shadow axis doesn't intersect Earth's surface
        are omitted from the results.

    Algorithm:
        For each time step, uses sol_eclipse_where() to find the central
        line position. This ensures consistency with _sol_eclipse_where_pythonic() by
        using the same iterative algorithm to find the point of minimum
        Sun-Moon angular separation on Earth's surface.

    Precision:
        Geographic coordinates accurate to approximately 0.001 degrees (~100 m)
        for points along the central line, consistent with _sol_eclipse_where_pythonic().

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, calc_eclipse_central_line
        >>> # Find April 8, 2024 total solar eclipse
        >>> jd = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times_ecl = _sol_eclipse_when_glob_pythonic(jd)
        >>> jd_c1, jd_c4 = times_ecl[1], times_ecl[4]  # First and fourth contacts
        >>> # Calculate central line path
        >>> times, lats, lons = calc_eclipse_central_line(jd_c1, jd_c4, step_minutes=5.0)
        >>> for i in range(len(times)):
        ...     print(f"JD {times[i]:.5f}: lat={lats[i]:.2f}°, lon={lons[i]:.2f}°")

    Note:
        - The function only returns points where the central eclipse is visible on Earth.
        - For partial-only eclipses, an empty tuple is returned.
        - The central line is only defined for central eclipses (total, annular, or hybrid).
        - Points near the beginning and end may be at extreme latitudes as the shadow
          enters and exits Earth.

    See Also:
        - _sol_eclipse_when_glob_pythonic: Find the next solar eclipse
        - _sol_eclipse_where_pythonic: Find central line point at a specific time
        - sol_eclipse_where: Underlying function for central line calculation

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac, Ch. 11
    """
    times_list: list[float] = []
    latitudes_list: list[float] = []
    longitudes_list: list[float] = []

    # Convert step to days
    step_days = step_minutes / (24.0 * 60.0)

    # Iterate through time range
    jd = jd_start
    while jd <= jd_end:
        # Use sol_eclipse_where for consistent results
        # This function uses iterative refinement to find the exact central line position
        retflag, geopos, attr = sol_eclipse_where(jd, flags)

        # Check if we have a valid central eclipse point
        # retflag > 0 indicates an eclipse is happening
        # geopos[0] and geopos[1] contain the central line coordinates
        if retflag > 0:
            lon = geopos[0]
            lat = geopos[1]

            # Check if we have valid coordinates (not zeros or placeholder values)
            # Also verify we're at a central eclipse point using Besselian gamma
            x = calc_besselian_x(jd, flags)
            y = calc_besselian_y(jd, flags)
            gamma = math.sqrt(x * x + y * y)

            # Include point if gamma indicates central eclipse is on Earth's surface
            # gamma < ~1.5 allows for Earth's oblateness and shadow geometry
            if gamma < 1.5 and (abs(lon) > 0.0001 or abs(lat) > 0.0001):
                times_list.append(jd)
                latitudes_list.append(lat)
                longitudes_list.append(lon)

        jd += step_days

    return tuple(times_list), tuple(latitudes_list), tuple(longitudes_list)


# Alias for reference API naming convention
calc_eclipse_central_line = calc_eclipse_central_line


def calc_eclipse_northern_limit(
    jd_start: float,
    jd_end: float,
    step_minutes: float = 1.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], Tuple[float, ...]]:
    """
    Calculate the geographic coordinates of points along the northern limit of the umbral shadow.

    The northern limit is the northern boundary of the path of totality or annularity.
    Along this line, observers see the edge of the umbral (or antumbral) shadow passing
    overhead - they are at the extreme northern edge of where totality/annularity is visible.

    This function computes a series of (latitude, longitude) points along the northern
    limit at specified time intervals, using Besselian elements for accurate calculations.

    Args:
        jd_start: Julian Day (UT) to start calculating the northern limit.
                  Should be during a central solar eclipse (from _sol_eclipse_when_glob_pythonic).
        jd_end: Julian Day (UT) to end the calculation.
                Should be after jd_start and during the same eclipse.
        step_minutes: Time step in minutes between calculated points (default 1.0).
                      Smaller values give a more detailed path but take longer.
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing three tuples:
            - times: Tuple of Julian Day (UT) values for each point
            - latitudes: Tuple of geographic latitudes in degrees (North positive)
            - longitudes: Tuple of geographic longitudes in degrees (East positive)

        Points where the shadow axis doesn't intersect Earth's surface (gamma > 1)
        or where the umbral shadow doesn't touch Earth are omitted from the results.

    Algorithm:
        For each time step:
        1. Calculate Besselian elements (x, y, d, mu, l2) using the library functions
        2. Calculate gamma = sqrt(x² + y²), the distance from shadow axis to Earth center
        3. Calculate l2, the umbral/antumbral cone radius at the fundamental plane
        4. If gamma + |l2| < ~1.5 (shadow touches Earth):
           - Offset the central line position northward by the umbral radius
           - Account for the shadow cone geometry and Earth's curvature
           - Convert to geographic coordinates using spherical trigonometry
        5. Collect all valid points into the result tuples

    Precision:
        Geographic coordinates accurate to approximately 0.1 degrees (~10 km)
        for points along the northern limit. Accuracy is best at mid-eclipse
        and decreases near the ends of the path where grazing geometry occurs.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, calc_eclipse_northern_limit
        >>> # Find April 8, 2024 total solar eclipse
        >>> jd = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times_ecl = _sol_eclipse_when_glob_pythonic(jd)
        >>> jd_c1, jd_c4 = times_ecl[1], times_ecl[4]  # First and fourth contacts
        >>> # Calculate northern limit path
        >>> times, lats, lons = calc_eclipse_northern_limit(jd_c1, jd_c4, step_minutes=5.0)
        >>> for i in range(len(times)):
        ...     print(f"JD {times[i]:.5f}: lat={lats[i]:.2f}°, lon={lons[i]:.2f}°")

    Note:
        - The function only returns points where the umbral/antumbral shadow touches Earth.
        - For partial-only eclipses, an empty tuple is returned.
        - The northern limit is only defined for central eclipses (total, annular, or hybrid).
        - The northern limit latitude is always greater than or equal to the central line
          latitude at the same time (in the northern hemisphere sense).
        - Near the ends of the eclipse path, the northern limit may curve sharply as the
          shadow enters or exits Earth's surface at a grazing angle.

    See Also:
        - calc_eclipse_central_line: Calculate the central line of the eclipse
        - calc_eclipse_path_width: Calculate the width of the shadow path
        - _sol_eclipse_where_pythonic: Find central line coordinates at a specific time
        - calc_besselian_l2: Calculate the umbral/antumbral cone radius

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac, Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    # WGS84 ellipsoid parameters
    EARTH_FLATTENING = 1.0 / 298.257223563

    times_list: list[float] = []
    latitudes_list: list[float] = []
    longitudes_list: list[float] = []

    # Convert step to days
    step_days = step_minutes / (24.0 * 60.0)

    # Iterate through time range
    jd = jd_start
    while jd <= jd_end:
        # Get Besselian elements at this time
        x = calc_besselian_x(jd, flags)
        y = calc_besselian_y(jd, flags)
        d = calc_besselian_d(jd, flags)
        mu = calc_besselian_mu(jd, flags)
        l2 = calc_besselian_l2(jd, flags)

        # Calculate gamma - distance from shadow axis to Earth center
        gamma = math.sqrt(x * x + y * y)

        # The umbral radius in Earth radii (absolute value)
        # l2 is negative for umbra (total), positive for antumbra (annular)
        umbra_radius = abs(l2)

        # For the northern limit to touch Earth, we need gamma + umbra_radius < ~1.5
        # (accounting for Earth's oblateness)
        max_gamma = 1.0 + EARTH_FLATTENING + umbra_radius
        if gamma < max_gamma and umbra_radius > 0.001:
            # Calculate the position of the northern limit
            # The northern limit is offset from the shadow axis in the +y direction
            # (north) by the umbral radius

            # Convert d (declination of shadow axis) and mu (hour angle) to radians
            d_rad = math.radians(d)
            mu_rad = math.radians(mu)

            # The perpendicular distance from Earth center to shadow axis is gamma
            # The y-coordinate points north in the fundamental plane

            # For the northern limit, we offset y by the umbral radius
            # The direction perpendicular to the shadow motion and toward north
            # depends on the shadow axis orientation

            # In the fundamental plane coordinate system:
            # x = east-west displacement (positive east)
            # y = north-south displacement (positive north)
            # The umbra has radius l2, so the northern edge is at y + |l2|

            y_north = y + umbra_radius

            # Calculate gamma for the northern limit point
            gamma_north = math.sqrt(x * x + y_north * y_north)

            # Check if this point is on Earth's surface
            max_gamma_surface = 1.0 + EARTH_FLATTENING
            if gamma_north < max_gamma_surface:
                # Calculate z-factor (height above fundamental plane)
                if gamma_north > 0.9999:
                    z_factor = 0.0
                else:
                    z_factor = math.sqrt(max(0, 1.0 - gamma_north * gamma_north))

                sin_d = math.sin(d_rad)
                cos_d = math.cos(d_rad)

                # Calculate latitude using the y_north component and shadow axis declination
                sin_lat = y_north * cos_d + z_factor * sin_d

                # Clamp to valid range
                sin_lat = max(-1.0, min(1.0, sin_lat))
                lat = math.degrees(math.asin(sin_lat))

                # For longitude, use the hour angle and x displacement
                if abs(cos_d) > 0.001:
                    cos_lat = math.cos(math.radians(lat))
                    if cos_lat > 0.001:
                        lon_offset = math.degrees(math.atan2(x, z_factor * cos_d))
                    else:
                        lon_offset = 0.0
                else:
                    lon_offset = 0.0

                # The longitude of the northern limit point
                lon = -mu + lon_offset

                # Apply correction for Earth's oblateness
                lat_geodetic = lat * (
                    1.0 + EARTH_FLATTENING * math.sin(math.radians(lat)) ** 2
                )

                # Normalize longitude to -180 to +180
                lon = ((lon + 180.0) % 360.0) - 180.0

                # Store this point
                times_list.append(jd)
                latitudes_list.append(lat_geodetic)
                longitudes_list.append(lon)

        jd += step_days

    return tuple(times_list), tuple(latitudes_list), tuple(longitudes_list)


# Alias for reference API naming convention
calc_eclipse_northern_limit = calc_eclipse_northern_limit


def calc_eclipse_southern_limit(
    jd_start: float,
    jd_end: float,
    step_minutes: float = 1.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[Tuple[float, ...], Tuple[float, ...], Tuple[float, ...]]:
    """
    Calculate the geographic coordinates of points along the southern limit of the umbral shadow.

    The southern limit is the southern boundary of the path of totality or annularity.
    Along this line, observers see the edge of the umbral (or antumbral) shadow passing
    overhead - they are at the extreme southern edge of where totality/annularity is visible.

    This function computes a series of (latitude, longitude) points along the southern
    limit at specified time intervals, using Besselian elements for accurate calculations.

    Args:
        jd_start: Julian Day (UT) to start calculating the southern limit.
                  Should be during a central solar eclipse (from _sol_eclipse_when_glob_pythonic).
        jd_end: Julian Day (UT) to end the calculation.
                Should be after jd_start and during the same eclipse.
        step_minutes: Time step in minutes between calculated points (default 1.0).
                      Smaller values give a more detailed path but take longer.
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing three tuples:
            - times: Tuple of Julian Day (UT) values for each point
            - latitudes: Tuple of geographic latitudes in degrees (North positive)
            - longitudes: Tuple of geographic longitudes in degrees (East positive)

        Points where the shadow axis doesn't intersect Earth's surface (gamma > 1)
        or where the umbral shadow doesn't touch Earth are omitted from the results.

    Algorithm:
        For each time step:
        1. Calculate Besselian elements (x, y, d, mu, l2) using the library functions
        2. Calculate gamma = sqrt(x² + y²), the distance from shadow axis to Earth center
        3. Calculate l2, the umbral/antumbral cone radius at the fundamental plane
        4. If gamma + |l2| < ~1.5 (shadow touches Earth):
           - Offset the central line position southward by the umbral radius
           - Account for the shadow cone geometry and Earth's curvature
           - Convert to geographic coordinates using spherical trigonometry
        5. Collect all valid points into the result tuples

    Precision:
        Geographic coordinates accurate to approximately 0.1 degrees (~10 km)
        for points along the southern limit. Accuracy is best at mid-eclipse
        and decreases near the ends of the path where grazing geometry occurs.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_when_glob_pythonic, calc_eclipse_southern_limit
        >>> # Find April 8, 2024 total solar eclipse
        >>> jd = julday(2024, 1, 1, 0.0)
        >>> ecl_type, times_ecl = _sol_eclipse_when_glob_pythonic(jd)
        >>> jd_c1, jd_c4 = times_ecl[1], times_ecl[4]  # First and fourth contacts
        >>> # Calculate southern limit path
        >>> times, lats, lons = calc_eclipse_southern_limit(jd_c1, jd_c4, step_minutes=5.0)
        >>> for i in range(len(times)):
        ...     print(f"JD {times[i]:.5f}: lat={lats[i]:.2f}°, lon={lons[i]:.2f}°")

    Note:
        - The function only returns points where the umbral/antumbral shadow touches Earth.
        - For partial-only eclipses, an empty tuple is returned.
        - The southern limit is only defined for central eclipses (total, annular, or hybrid).
        - The southern limit latitude is always less than or equal to the central line
          latitude at the same time (in the northern hemisphere sense).
        - Near the ends of the eclipse path, the southern limit may curve sharply as the
          shadow enters or exits Earth's surface at a grazing angle.

    See Also:
        - calc_eclipse_central_line: Calculate the central line of the eclipse
        - calc_eclipse_northern_limit: Calculate the northern limit of the eclipse
        - calc_eclipse_path_width: Calculate the width of the shadow path
        - _sol_eclipse_where_pythonic: Find central line coordinates at a specific time
        - calc_besselian_l2: Calculate the umbral/antumbral cone radius

    References:
        - Meeus, J. "Astronomical Algorithms", Ch. 54 (Solar Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
        - Explanatory Supplement to the Astronomical Almanac, Ch. 11
        - Chauvenet's "Manual of Spherical and Practical Astronomy", Vol. 1
    """
    # WGS84 ellipsoid parameters
    EARTH_FLATTENING = 1.0 / 298.257223563

    times_list: list[float] = []
    latitudes_list: list[float] = []
    longitudes_list: list[float] = []

    # Convert step to days
    step_days = step_minutes / (24.0 * 60.0)

    # Iterate through time range
    jd = jd_start
    while jd <= jd_end:
        # Get Besselian elements at this time
        x = calc_besselian_x(jd, flags)
        y = calc_besselian_y(jd, flags)
        d = calc_besselian_d(jd, flags)
        mu = calc_besselian_mu(jd, flags)
        l2 = calc_besselian_l2(jd, flags)

        # Calculate gamma - distance from shadow axis to Earth center
        gamma = math.sqrt(x * x + y * y)

        # The umbral radius in Earth radii (absolute value)
        # l2 is negative for umbra (total), positive for antumbra (annular)
        umbra_radius = abs(l2)

        # For the southern limit to touch Earth, we need gamma + umbra_radius < ~1.5
        # (accounting for Earth's oblateness)
        max_gamma = 1.0 + EARTH_FLATTENING + umbra_radius
        if gamma < max_gamma and umbra_radius > 0.001:
            # Calculate the position of the southern limit
            # The southern limit is offset from the shadow axis in the -y direction
            # (south) by the umbral radius

            # Convert d (declination of shadow axis) and mu (hour angle) to radians
            d_rad = math.radians(d)
            mu_rad = math.radians(mu)

            # The perpendicular distance from Earth center to shadow axis is gamma
            # The y-coordinate points north in the fundamental plane

            # For the southern limit, we offset y by the umbral radius in the
            # negative direction (southward)
            # The direction perpendicular to the shadow motion and toward south
            # depends on the shadow axis orientation

            # In the fundamental plane coordinate system:
            # x = east-west displacement (positive east)
            # y = north-south displacement (positive north)
            # The umbra has radius l2, so the southern edge is at y - |l2|

            y_south = y - umbra_radius

            # Calculate gamma for the southern limit point
            gamma_south = math.sqrt(x * x + y_south * y_south)

            # Check if this point is on Earth's surface
            max_gamma_surface = 1.0 + EARTH_FLATTENING
            if gamma_south < max_gamma_surface:
                # Calculate z-factor (height above fundamental plane)
                if gamma_south > 0.9999:
                    z_factor = 0.0
                else:
                    z_factor = math.sqrt(max(0, 1.0 - gamma_south * gamma_south))

                sin_d = math.sin(d_rad)
                cos_d = math.cos(d_rad)

                # Calculate latitude using the y_south component and shadow axis declination
                sin_lat = y_south * cos_d + z_factor * sin_d

                # Clamp to valid range
                sin_lat = max(-1.0, min(1.0, sin_lat))
                lat = math.degrees(math.asin(sin_lat))

                # For longitude, use the hour angle and x displacement
                if abs(cos_d) > 0.001:
                    cos_lat = math.cos(math.radians(lat))
                    if cos_lat > 0.001:
                        lon_offset = math.degrees(math.atan2(x, z_factor * cos_d))
                    else:
                        lon_offset = 0.0
                else:
                    lon_offset = 0.0

                # The longitude of the southern limit point
                lon = -mu + lon_offset

                # Apply correction for Earth's oblateness
                lat_geodetic = lat * (
                    1.0 + EARTH_FLATTENING * math.sin(math.radians(lat)) ** 2
                )

                # Normalize longitude to -180 to +180
                lon = ((lon + 180.0) % 360.0) - 180.0

                # Store this point
                times_list.append(jd)
                latitudes_list.append(lat_geodetic)
                longitudes_list.append(lon)

        jd += step_days

    return tuple(times_list), tuple(latitudes_list), tuple(longitudes_list)


# Alias for reference API naming convention
calc_eclipse_southern_limit = calc_eclipse_southern_limit


def _sol_eclipse_magnitude_at_loc_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> float:
    """
    Calculate the eclipse magnitude at a specific geographic location and time.

    Eclipse magnitude is defined as the fraction of the Sun's diameter that
    is covered by the Moon. This is a simplified convenience function that
    returns just the magnitude value without additional eclipse attributes.

    Args:
        jd: Julian Day (UT) of the time to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Eclipse magnitude as a float:
            - 0.0 if no eclipse is visible (Sun and Moon not overlapping)
            - 0.0 to 1.0 for partial eclipses (fraction of diameter covered)
            - 1.0 for the moment of totality in a total eclipse
            - > 1.0 for total eclipses (the excess indicates how much larger
              the Moon appears than the Sun)

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous magnitude at the given time. To find eclipse events,
        use _sol_eclipse_when_glob_pythonic() or _sol_eclipse_when_loc_pythonic() first.

        If the Sun is below the horizon at the observer's location, the
        magnitude will be 0.0 since the eclipse is not visible.

    Algorithm:
        1. Calculate topocentric Sun and Moon apparent positions
        2. Compute angular separation between centers
        3. Compute angular radii of Sun and Moon
        4. Calculate overlap and express as fraction of Sun's diameter

    Precision:
        Magnitude accurate to ~0.001 for central eclipses.
        Topocentric parallax is included in calculations.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_magnitude_at_loc_pythonic
        >>> # Calculate magnitude during April 8, 2024 eclipse from Dallas
        >>> jd = 2460409.28  # During eclipse
        >>> dallas_lat, dallas_lon = 32.7767, -96.797
        >>> magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, dallas_lat, dallas_lon)
        >>> print(f"Magnitude: {magnitude:.4f}")

        >>> # Check magnitude at a location outside the eclipse path
        >>> london_lat, london_lon = 51.5074, -0.1278
        >>> magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, london_lat, london_lon)
        >>> print(f"London magnitude: {magnitude:.4f}")  # Will be 0.0

    References:
        - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
        - Reference documentation
    """
    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        _gp = (lon, lat, altitude)

        try:
            jd_tt = jd + deltat(jd)
            sun_topo = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_topo = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        sun_az_val, sun_alt_true, sun_alt_app = azalt(jd, ECL2HOR, _gp, 0, 0, sun_topo[:3])
        sun_altitude = sun_alt_true

        if sun_altitude < -1.0:
            return 0.0

        separation = angular_separation(sun_topo[0], sun_topo[1], moon_topo[0], moon_topo[1])
        sun_dist_au = sun_topo[2]
        moon_dist_au = moon_topo[2]
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        # Get Skyfield time
        t = ts.ut1_jd(jd)

        # Create observer position
        observer_at = earth + observer

        # Get apparent positions from observer (topocentric)
        try:
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # If calculation fails, return 0 magnitude
            return 0.0

        # Get Sun altitude to check visibility
        sun_alt, _, _ = sun_app.altaz()
        sun_altitude = sun_alt.degrees

        # If Sun is below horizon, no visible eclipse
        if sun_altitude < -1.0:  # Allow for refraction near horizon
            return 0.0

        # Calculate angular separation between Sun and Moon
        separation = sun_app.separation_from(moon_app).degrees

        # Get distances for angular size calculations
        sun_dist_au = sun_app.distance().au
        moon_dist_au = moon_app.distance().au

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

    sun_diameter = 2 * sun_angular_radius
    moon_diameter = 2 * moon_angular_radius

    # Check if there's any eclipse (disks overlapping)
    sum_radii = sun_angular_radius + moon_angular_radius
    if separation >= sum_radii:
        # No eclipse - Sun and Moon too far apart
        return 0.0

    # Calculate eclipse magnitude
    # Magnitude = fraction of Sun's diameter covered by Moon
    overlap = sum_radii - separation
    magnitude = overlap / sun_diameter

    # Clamp to valid range (can exceed 1.0 for total eclipse)
    magnitude = max(0.0, min(magnitude, 1.0 + moon_diameter / sun_diameter))

    return magnitude


def sol_eclipse_magnitude_at_loc(
    tjd_ut: float,
    geopos: Sequence[float],
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the eclipse magnitude at a specific geographic location and time.

    This function matches the reference API naming convention. It is a convenience
    wrapper that returns just the eclipse magnitude (fraction of solar diameter
    covered by the Moon) without the full attribute array.

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level

    Returns:
        Eclipse magnitude as a float:
            - 0.0 if no eclipse is visible
            - 0.0 to 1.0 for partial eclipses
            - >= 1.0 for total eclipses

    Raises:
        ValueError: If geopos has wrong length

    Example:
        >>> from libephemeris import sol_eclipse_magnitude_at_loc, FLG_SWIEPH
        >>> # Calculate magnitude during April 8, 2024 eclipse from Dallas
        >>> jd = 2460409.28
        >>> dallas_geopos = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> magnitude = sol_eclipse_magnitude_at_loc(jd, FLG_SWIEPH, dallas_geopos)
        >>> print(f"Magnitude: {magnitude:.4f}")

    References:
        - Reference API: sol_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # Validate geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    # Extract geographic position (longitude first, then latitude - reference API convention)
    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_magnitude_at_loc_pythonic,
        tjd_ut, lat, lon, altitude, ifl,
    )


def _sol_eclipse_obscuration_at_loc_pythonic(
    jd: float,
    lat: float,
    lon: float,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> float:
    """
    Calculate the eclipse obscuration at a specific geographic location and time.

    Eclipse obscuration is defined as the fraction of the Sun's disc area that
    is covered by the Moon. This differs from eclipse magnitude, which is the
    fraction of the Sun's diameter covered.

    The relationship between magnitude and obscuration is non-linear:
    - Obscuration = 0 when magnitude = 0 (no eclipse)
    - Obscuration < magnitude for partial eclipses
    - Obscuration = 1 when magnitude >= 1 (total eclipse)
    - For annular eclipses, obscuration = (moon_radius/sun_radius)^2

    Args:
        jd: Julian Day (UT) of the time to calculate
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Eclipse obscuration as a float:
            - 0.0 if no eclipse is visible (Sun and Moon not overlapping)
            - 0.0 to 1.0 for partial eclipses (fraction of Sun's area covered)
            - 1.0 for total eclipse (Moon completely covers Sun)
            - (moon_radius/sun_radius)^2 for annular eclipse (Moon inside Sun)

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous obscuration at the given time. To find eclipse events,
        use _sol_eclipse_when_glob_pythonic() or _sol_eclipse_when_loc_pythonic() first.

        If the Sun is below the horizon at the observer's location, the
        obscuration will be 0.0 since the eclipse is not visible.

    Algorithm:
        1. Calculate topocentric Sun and Moon apparent positions
        2. Compute angular separation between centers
        3. Compute angular radii of Sun and Moon
        4. Calculate the intersection area of two overlapping circles
        5. Express as fraction of Sun's total disc area

        The intersection area formula uses the lens/vesica piscis formula:
        For two circles with radii r1, r2 and center separation d, the
        intersection area involves the sum of two circular segments.

    Precision:
        Obscuration accurate to ~0.001 for central eclipses.
        Topocentric parallax is included in calculations.

    Example:
        >>> from libephemeris import julday, _sol_eclipse_obscuration_at_loc_pythonic
        >>> # Calculate obscuration during April 8, 2024 eclipse from Dallas
        >>> jd = 2460409.28  # During eclipse
        >>> dallas_lat, dallas_lon = 32.7767, -96.797
        >>> obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, dallas_lat, dallas_lon)
        >>> print(f"Obscuration: {obscuration:.4f}")

        >>> # Compare with magnitude
        >>> from libephemeris import _sol_eclipse_magnitude_at_loc_pythonic
        >>> magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, dallas_lat, dallas_lon)
        >>> print(f"Magnitude: {magnitude:.4f}, Obscuration: {obscuration:.4f}")

    References:
        - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
        - Reference documentation
        - Intersection of two circles formula (computational geometry)
    """
    # reader is provided by the caller (None forces Skyfield path)

    if reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR, angular_separation

        _gp = (lon, lat, altitude)

        try:
            jd_tt = jd + deltat(jd)
            sun_topo = _topo_ecliptic(reader, jd_tt, jd, SUN, _gp)
            moon_topo = _topo_ecliptic(reader, jd_tt, jd, MOON, _gp)
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            return 0.0

        sun_az_val, sun_alt_true, sun_alt_app = azalt(jd, ECL2HOR, _gp, 0, 0, sun_topo[:3])
        sun_altitude = sun_alt_true

        if sun_altitude < -1.0:
            return 0.0

        separation = angular_separation(sun_topo[0], sun_topo[1], moon_topo[0], moon_topo[1])
        sun_dist_au = sun_topo[2]
        moon_dist_au = moon_topo[2]
    else:
        from skyfield.api import wgs84
        from .state import get_planets, get_timescale

        # Get ephemeris and timescale
        eph = get_planets()
        ts = get_timescale()

        # Create observer location
        observer = wgs84.latlon(lat, lon, altitude)

        # Get Sun and Moon objects
        sun = eph["sun"]
        moon = eph["moon"]
        earth = eph["earth"]

        # Get Skyfield time
        t = ts.ut1_jd(jd)

        # Create observer position
        observer_at = earth + observer

        # Get apparent positions from observer (topocentric)
        try:
            sun_app = observer_at.at(t).observe(sun).apparent()
            moon_app = observer_at.at(t).observe(moon).apparent()
        except (KeyError, ValueError, ArithmeticError, IndexError) as _exc:
            _reraise_if_leb_range_error(_exc)
            # If calculation fails, return 0 obscuration
            return 0.0

        # Get Sun altitude to check visibility
        sun_alt, _, _ = sun_app.altaz()
        sun_altitude = sun_alt.degrees

        # If Sun is below horizon, no visible eclipse
        if sun_altitude < -1.0:  # Allow for refraction near horizon
            return 0.0

        # Calculate angular separation between Sun and Moon
        separation = sun_app.separation_from(moon_app).degrees

        # Get distances for angular size calculations
        sun_dist_au = sun_app.distance().au
        moon_dist_au = moon_app.distance().au

    # Angular radii (in degrees)
    # Sun: mean radius 959.63" at 1 AU
    # Moon: mean radius 932.56" at mean distance (0.002569 AU)
    sun_angular_radius = (959.63 / 3600.0) / sun_dist_au
    moon_angular_radius = (932.56 / 3600.0) * (0.002569 / moon_dist_au)

    r_sun = sun_angular_radius
    r_moon = moon_angular_radius
    d = separation  # center-to-center separation

    # Calculate obscuration (fraction of Sun's area covered)
    if d >= r_sun + r_moon:
        # No overlap - no eclipse
        return 0.0
    elif d <= abs(r_sun - r_moon):
        # One disk entirely within the other
        # Both total and annular: obscuration = disc area ratio
        # For total eclipses (r_moon >= r_sun), this is >= 1.0
        return (r_moon / r_sun) ** 2
    else:
        # Partial overlap - use lens formula for intersection of two circles
        # The intersection area is the sum of two circular segments
        #
        # For two circles with radii r1, r2 and center separation d:
        # d1 = distance from center of circle 1 to the chord
        # d2 = distance from center of circle 2 to the chord
        # d1 = (d^2 + r1^2 - r2^2) / (2*d)
        # d2 = d - d1
        #
        # Each segment area = r^2 * arccos(d_i/r) - d_i * sqrt(r^2 - d_i^2)

        d1 = (d * d + r_sun * r_sun - r_moon * r_moon) / (2 * d)
        d2 = d - d1

        if abs(d1) <= r_sun and abs(d2) <= r_moon:
            # Ensure values are valid for acos
            cos_arg1 = max(-1, min(1, d1 / r_sun))
            cos_arg2 = max(-1, min(1, d2 / r_moon))

            # Segment areas
            area1 = r_sun * r_sun * math.acos(cos_arg1) - d1 * math.sqrt(
                max(0, r_sun * r_sun - d1 * d1)
            )
            area2 = r_moon * r_moon * math.acos(cos_arg2) - d2 * math.sqrt(
                max(0, r_moon * r_moon - d2 * d2)
            )

            intersection_area = area1 + area2
            sun_area = math.pi * r_sun * r_sun
            obscuration = intersection_area / sun_area
        else:
            obscuration = 0.0

    return max(0.0, obscuration)


def sol_eclipse_obscuration_at_loc(
    tjd_ut: float,
    geopos: Sequence[float],
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the eclipse obscuration at a specific geographic location and time.

    This function matches the reference API naming convention. It is a convenience
    wrapper that returns just the eclipse obscuration (fraction of solar disc
    area covered by the Moon) without the full attribute array.

    Obscuration differs from magnitude:
    - Magnitude: fraction of Sun's DIAMETER covered
    - Obscuration: fraction of Sun's AREA covered

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)
        geopos: Sequence of [longitude, latitude, altitude]:
            - longitude in degrees (East positive)
            - latitude in degrees (North positive)
            - altitude in meters above sea level

    Returns:
        Eclipse obscuration as a float:
            - 0.0 if no eclipse is visible
            - 0.0 to 1.0 for partial/annular eclipses
            - 1.0 for total eclipse

    Raises:
        ValueError: If geopos has wrong length

    Example:
        >>> from libephemeris import sol_eclipse_obscuration_at_loc, FLG_SWIEPH
        >>> # Calculate obscuration during April 8, 2024 eclipse from Dallas
        >>> jd = 2460409.28
        >>> dallas_geopos = [-96.797, 32.7767, 0]  # lon, lat, alt
        >>> obs = sol_eclipse_obscuration_at_loc(jd, FLG_SWIEPH, dallas_geopos)
        >>> print(f"Obscuration: {obs:.4f}")

    References:
        - Reference API: sol_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    # Validate geopos
    if len(geopos) < 3:
        raise ValueError("geopos must have at least 3 elements: [lon, lat, alt]")

    # Extract geographic position (lon first, lat second - reference API convention)
    lon = float(geopos[0])
    lat = float(geopos[1])
    altitude = float(geopos[2])

    return _call_with_leb_skyfield_fallback(
        _sol_eclipse_obscuration_at_loc_pythonic,
        tjd_ut, lat, lon, altitude, ifl,
    )


def _lun_eclipse_umbral_magnitude_pythonic(
    jd: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the umbral magnitude for a lunar eclipse at a specific time.

    Umbral magnitude is defined as the fraction of the Moon's diameter that
    is within Earth's umbral (dark) shadow. This is a simplified convenience
    function that returns just the umbral magnitude value.

    Unlike solar eclipse magnitude (which depends on observer location), lunar
    eclipse magnitude is the same for all observers who can see the Moon,
    since the Moon physically enters Earth's shadow.

    Args:
        jd: Julian Day (UT) of the time to calculate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Umbral magnitude as a float:
            - 0.0 if the Moon is not in the umbral shadow (penumbral-only
              eclipse or no eclipse)
            - 0.0 to 1.0 for partial umbral eclipses (fraction of Moon's
              diameter within umbra)
            - >= 1.0 for total lunar eclipses (Moon fully within umbra;
              values > 1.0 indicate how deep the Moon is in the shadow)

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous umbral magnitude at the given time. To find eclipse
        events, use _lun_eclipse_when_pythonic() first.

        For penumbral-only eclipses, this function returns 0.0 since the
        Moon has not entered the umbra. Use _lun_eclipse_how_pythonic() to get the
        penumbral magnitude in such cases.

    Algorithm:
        1. Calculate Moon's ecliptic position
        2. Calculate Earth's umbral shadow cone geometry at Moon's distance
        3. Determine how much of the Moon's diameter is within the umbra

    Precision:
        Magnitude accurate to ~0.01 for typical eclipses.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_umbral_magnitude_pythonic, _lun_eclipse_when_pythonic
        >>> # First find a lunar eclipse
        >>> jd_start = julday(2022, 5, 1, 0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start)
        >>> jd_max = times[0]  # Time of maximum eclipse
        >>> # Calculate umbral magnitude at maximum
        >>> umbral_mag = _lun_eclipse_umbral_magnitude_pythonic(jd_max)
        >>> print(f"Umbral magnitude: {umbral_mag:.4f}")

        >>> # Check magnitude at a random time (no eclipse)
        >>> jd_no_eclipse = julday(2022, 6, 1, 12.0)
        >>> mag = _lun_eclipse_umbral_magnitude_pythonic(jd_no_eclipse)
        >>> print(f"Magnitude: {mag:.4f}")  # Will be 0.0

    References:
        - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
        - Reference documentation
    """
    # Use the existing calculation function
    (
        ecl_type_flags,
        umbral_mag,
        penumbral_mag,
        gamma,
        penumbra_radius,
        umbra_radius,
    ) = _calculate_lunar_eclipse_type_and_magnitude(jd)

    # Return the umbral magnitude, clamped to non-negative
    return max(0.0, umbral_mag)


def lun_eclipse_umbral_magnitude(
    tjd_ut: float,
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the umbral magnitude for a lunar eclipse at a specific time.

    This function matches the reference API naming convention. It is a convenience
    function that returns just the umbral magnitude (fraction of Moon's diameter
    within Earth's umbral shadow).

    Unlike solar eclipses, lunar eclipse magnitude does not depend on observer
    location - the Moon physically enters Earth's shadow, so the magnitude is
    the same for all observers who can see the Moon.

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Umbral magnitude as a float:
            - 0.0 if Moon is not in umbra (penumbral-only or no eclipse)
            - 0.0 to 1.0 for partial umbral eclipses
            - >= 1.0 for total lunar eclipses

    Example:
        >>> from libephemeris import lun_eclipse_umbral_magnitude, FLG_SWIEPH
        >>> # Calculate umbral magnitude during Nov 8, 2022 total lunar eclipse
        >>> jd = 2459892.0  # During eclipse
        >>> umbral_mag = lun_eclipse_umbral_magnitude(jd, FLG_SWIEPH)
        >>> print(f"Umbral magnitude: {umbral_mag:.4f}")

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    return _lun_eclipse_umbral_magnitude_pythonic(tjd_ut, ifl)


def _lun_eclipse_penumbral_magnitude_pythonic(
    jd: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the penumbral magnitude for a lunar eclipse at a specific time.

    Penumbral magnitude is defined as the fraction of the Moon's diameter that
    is within Earth's penumbral (partial) shadow. This is a simplified convenience
    function that returns just the penumbral magnitude value.

    Unlike solar eclipse magnitude (which depends on observer location), lunar
    eclipse magnitude is the same for all observers who can see the Moon,
    since the Moon physically enters Earth's shadow.

    Args:
        jd: Julian Day (UT) of the time to calculate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Penumbral magnitude as a float:
            - 0.0 if the Moon is not in the penumbral shadow (no eclipse)
            - 0.0 to 1.0 for penumbral-only eclipses or during penumbral
              phases of umbral eclipses
            - > 1.0 when the Moon is fully within the penumbra (and possibly
              also within the umbra); values > 1.0 indicate how deep the Moon
              is in the penumbral shadow

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous penumbral magnitude at the given time. To find eclipse
        events, use _lun_eclipse_when_pythonic() first.

        For all lunar eclipses (penumbral, partial, or total), the penumbral
        magnitude will be > 0 during the eclipse. For penumbral-only eclipses,
        use this function to determine the eclipse magnitude.

    Algorithm:
        1. Calculate Moon's ecliptic position
        2. Calculate Earth's penumbral shadow cone geometry at Moon's distance
        3. Determine how much of the Moon's diameter is within the penumbra

    Precision:
        Magnitude accurate to ~0.01 for typical eclipses.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_penumbral_magnitude_pythonic, _lun_eclipse_when_pythonic
        >>> from libephemeris import ECL_PENUMBRAL
        >>> # First find a penumbral lunar eclipse
        >>> jd_start = julday(2020, 1, 1, 0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start, eclipse_type=ECL_PENUMBRAL)
        >>> jd_max = times[0]  # Time of maximum eclipse
        >>> # Calculate penumbral magnitude at maximum
        >>> penumbral_mag = _lun_eclipse_penumbral_magnitude_pythonic(jd_max)
        >>> print(f"Penumbral magnitude: {penumbral_mag:.4f}")

        >>> # Check magnitude at a random time (no eclipse)
        >>> jd_no_eclipse = julday(2022, 6, 1, 12.0)
        >>> mag = _lun_eclipse_penumbral_magnitude_pythonic(jd_no_eclipse)
        >>> print(f"Magnitude: {mag:.4f}")  # Will be 0.0

    References:
        - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
        - Reference documentation
    """
    # Use the existing calculation function
    (
        ecl_type_flags,
        umbral_mag,
        penumbral_mag,
        gamma,
        penumbra_radius,
        umbra_radius,
    ) = _calculate_lunar_eclipse_type_and_magnitude(jd)

    # Return the penumbral magnitude, clamped to non-negative
    return max(0.0, penumbral_mag)


def lun_eclipse_penumbral_magnitude(
    tjd_ut: float,
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the penumbral magnitude for a lunar eclipse at a specific time.

    This function matches the reference API naming convention. It is a convenience
    function that returns just the penumbral magnitude (fraction of Moon's diameter
    within Earth's penumbral shadow).

    Unlike solar eclipses, lunar eclipse magnitude does not depend on observer
    location - the Moon physically enters Earth's shadow, so the magnitude is
    the same for all observers who can see the Moon.

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Penumbral magnitude as a float:
            - 0.0 if Moon is not in penumbra (no eclipse)
            - 0.0 to 1.0 for penumbral-only eclipses (partial immersion)
            - > 1.0 when Moon is fully within the penumbra

    Example:
        >>> from libephemeris import lun_eclipse_penumbral_magnitude, FLG_SWIEPH
        >>> # Calculate penumbral magnitude during Jan 10, 2020 penumbral lunar eclipse
        >>> jd = 2458859.0  # During eclipse
        >>> penumbral_mag = lun_eclipse_penumbral_magnitude(jd, FLG_SWIEPH)
        >>> print(f"Penumbral magnitude: {penumbral_mag:.4f}")

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    return _lun_eclipse_penumbral_magnitude_pythonic(tjd_ut, ifl)


def _lun_eclipse_gamma_pythonic(
    jd: float,
    flags: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the gamma parameter for a lunar eclipse at a specific time.

    Gamma is the distance of the Moon's center from Earth's shadow axis,
    measured in Earth radii. It is a fundamental parameter for characterizing
    lunar eclipses, indicating how centrally the Moon passes through the shadow.

    Unlike solar eclipse gamma (which uses the fundamental plane), lunar eclipse
    gamma is measured perpendicular to the shadow axis in the plane containing
    the Moon. A gamma of 0 means the Moon passes exactly through the center of
    Earth's shadow.

    The sign of gamma indicates which side of the shadow axis the Moon passes:
        - Positive gamma: Moon passes north of shadow axis
        - Negative gamma: Moon passes south of shadow axis

    Args:
        jd: Julian Day (UT) of the time to calculate
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Gamma value as a float:
            - 0.0: Moon center is on the shadow axis (most central eclipse)
            - |gamma| < ~0.25: Deep total eclipse (Moon well within umbra)
            - |gamma| < ~0.75: Total eclipse possible
            - |gamma| < ~1.0: Partial umbral eclipse possible
            - |gamma| < ~1.5: Penumbral eclipse possible
            - |gamma| > ~1.5: No eclipse (Moon misses the shadow)

    Note:
        This function does NOT search for eclipses - it calculates the
        instantaneous gamma at the given time. To find eclipse events,
        use _lun_eclipse_when_pythonic() first.

        The gamma parameter is useful for:
        - Classifying eclipse centrality
        - Predicting eclipse magnitude
        - Analyzing Saros series patterns

    Algorithm:
        1. Calculate Moon's ecliptic latitude
        2. Calculate Earth's angular semi-diameter as seen from Moon
        3. Compute gamma as the ratio of lunar latitude to Earth's angular radius

    Precision:
        Gamma accurate to ~0.001 for typical eclipses.

    Example:
        >>> from libephemeris import julday, _lun_eclipse_gamma_pythonic, _lun_eclipse_when_pythonic
        >>> # First find a lunar eclipse
        >>> jd_start = julday(2022, 5, 1, 0)
        >>> ecl_type, times = _lun_eclipse_when_pythonic(jd_start)
        >>> jd_max = times[0]  # Time of maximum eclipse
        >>> # Calculate gamma at maximum
        >>> gamma = _lun_eclipse_gamma_pythonic(jd_max)
        >>> print(f"Gamma: {gamma:.4f}")

        >>> # Check gamma at a random time (far from eclipse)
        >>> jd_no_eclipse = julday(2022, 6, 1, 12.0)
        >>> gamma = _lun_eclipse_gamma_pythonic(jd_no_eclipse)
        >>> print(f"Gamma: {gamma:.4f}")  # Will be large (no eclipse)

    References:
        - Meeus "Astronomical Algorithms" Ch. 54 (Eclipses)
        - Espenak & Meeus "Five Millennium Canon of Lunar Eclipses"
        - Reference documentation
    """
    # Use the existing calculation function
    (
        ecl_type_flags,
        umbral_mag,
        penumbral_mag,
        gamma,
        penumbra_radius,
        umbra_radius,
    ) = _calculate_lunar_eclipse_type_and_magnitude(jd)

    return gamma


def lun_eclipse_gamma(
    tjd_ut: float,
    ifl: int = FLG_SWIEPH,
) -> float:
    """
    Calculate the gamma parameter for a lunar eclipse at a specific time.

    This function matches the reference API naming convention. Gamma represents
    the distance of the Moon's center from Earth's shadow axis, measured in
    Earth radii.

    Unlike solar eclipse gamma, lunar eclipse gamma does not depend on observer
    location - it is a geometric property of the Moon's passage through Earth's
    shadow.

    Args:
        tjd_ut: Julian Day (UT) of the time to calculate
        ifl: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Gamma value as a float:
            - Positive: Moon passes north of shadow axis
            - Negative: Moon passes south of shadow axis
            - |gamma| ~ 0: Most central eclipse
            - |gamma| > ~1.5: No eclipse

    Example:
        >>> from libephemeris import lun_eclipse_gamma, FLG_SWIEPH
        >>> # Calculate gamma during Nov 8, 2022 total lunar eclipse
        >>> jd = 2459892.0  # During eclipse
        >>> gamma = lun_eclipse_gamma(jd, FLG_SWIEPH)
        >>> print(f"Gamma: {gamma:.4f}")

    References:
        - Reference API: lun_eclipse_how()
        - Meeus "Astronomical Algorithms" Ch. 54
    """
    return _lun_eclipse_gamma_pythonic(tjd_ut, ifl)


def planet_occult_when_glob(
    tjdut: float,
    occulting_planet: int,
    occulted_planet: int = 0,
    starname: str = "",
    flags: int = FLG_SWIEPH,
    direction: int = 0,
) -> Tuple[int, Tuple[float, ...]]:
    """
    Find the next planetary occultation globally (UT).

    A planetary occultation occurs when one planet passes in front of (occults)
    another planet or star as seen from Earth. This is different from lunar
    occultations - here the occulting body is a planet (e.g., Venus, Jupiter).

    Planetary occultations are rare events. Mutual occultations between planets
    typically occur only a few times per century.

    Args:
        tjdut: Julian Day (UT) to start search from
        occulting_planet: Planet ID of the occulting (foreground) planet.
            Use VENUS, MARS, JUPITER, etc.
        occulted_planet: Planet ID of the occulted (background) planet.
            Use 0 if searching for a star occultation.
        starname: Star name (str). Use empty string "" if searching for a planet.
        flags: Calculation flags (FLG_SWIEPH, etc.)
        direction: Search direction. 0 or positive = forward in time,
                   negative = backward in time.

    Returns:
        Tuple containing:
            - retflags: Returned bit flags (int):
                - 0 if no occultation found
                - ECL_TOTAL or ECL_PARTIAL
            - tret: Tuple of 10 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation
                [1]: Reserved (0)
                [2]: Time of occultation begin
                [3]: Time of occultation end
                [4]: Time of totality begin (if total, else 0)
                [5]: Time of totality end (if total, else 0)
                [6-9]: Reserved (0)

    Raises:
        RuntimeError: If no occultation found within search limit
        ValueError: If neither occulted_planet nor starname is specified

    Algorithm:
        1. Calculate both planets' positions over time
        2. Find conjunctions in right ascension
        3. At conjunction, check angular separation
        4. If separation < occulting planet's angular radius + occulted body's radius,
           an occultation occurs
        5. Refine timing using golden section search
        6. Calculate contact times

    Historical Events:
        - 1818 Dec 3: Venus occulted Jupiter
        - 2065 Nov 22: Venus will occult Jupiter
        - 2123 Sep 14: Venus will occult Jupiter

    Example:
        >>> # Find next planetary occultation of Jupiter by Venus
        >>> from libephemeris import julday, planet_occult_when_glob, VENUS, JUPITER
        >>> jd = julday(2060, 1, 1, 0)
        >>> retflags, tret = planet_occult_when_glob(jd, VENUS, JUPITER)
        >>> print(f"Occultation at JD {tret[0]:.5f}")

        >>> # Find next occultation of Regulus by Venus
        >>> retflags, tret = planet_occult_when_glob(jd, VENUS, 0, "Regulus")

    References:
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
        - Herald, D. & Sinnott, R. "Planetary Occultations"
    """
    from .state import get_planets, get_timescale
    from .fixed_stars import FIXED_STARS, _resolve_star_id
    from .constants import (
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
        SUN,
        MOON,
    )
    from .planets import _PLANET_MAP, get_planet_target

    if occulted_planet == 0 and not starname:
        raise ValueError(
            "Either occulted_planet ID or starname must be specified for occultation search"
        )

    # Validate occulting planet - can't be Sun or Moon (use other functions for those)
    if occulting_planet == SUN:
        raise ValueError("Sun cannot be an occulting body - use _sol_eclipse_when_glob_pythonic")
    if occulting_planet == MOON:
        raise ValueError("Moon cannot be an occulting body - use lun_occult_when_glob")
    if occulting_planet not in _PLANET_MAP:
        raise ValueError(f"Invalid occulting planet ID: {occulting_planet}")

    # Validate occulted planet if specified
    if occulted_planet != 0:
        if occulted_planet not in _PLANET_MAP:
            raise ValueError(f"Invalid occulted planet ID: {occulted_planet}")
        if occulted_planet == occulting_planet:
            raise ValueError("Occulting and occulted planets cannot be the same")

    jd_start = tjdut

    # Planetary occultations are very rare - search up to 150 years
    MAX_SEARCH_YEARS = 150
    MAX_ITERATIONS = int(MAX_SEARCH_YEARS * 365)  # Check roughly daily

    from .constants import FLG_EQUATORIAL

    def _get_planet_position(
        jd: float, planet_id: int
    ) -> Tuple[float, float, float, float]:
        """Get planet's geocentric RA, Dec, distance, and angular radius."""

        if planet_id not in _PLANET_MAP:
            raise ValueError(f"Invalid planet ID: {planet_id}")

        # Use calc_ut (LEB-aware) for equatorial positions
        planet_eq, _ = calc_ut(jd, planet_id, FLG_EQUATORIAL | FLG_SPEED)
        ra_deg = planet_eq[0]
        dec_deg = planet_eq[1]
        dist_au = planet_eq[2]

        angular_radius = _calc_planet_angular_radius(planet_id, dist_au)

        return ra_deg, dec_deg, dist_au, angular_radius

    def _get_target_position(jd: float) -> Tuple[float, float, float]:
        """Get target (occulted) body's geocentric RA, Dec, and angular radius."""
        if occulted_planet == 0:
            # Fixed star
            star_id, err, _ = _resolve_star_id(starname)
            if err is not None:
                raise ValueError(err)

            star = FIXED_STARS[star_id]

            # Time from J2000.0
            t_years = (jd - 2451545.0) / 365.25

            # Apply proper motion
            ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
            dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0

            # Stars have negligible angular radius
            return ra_deg, dec_deg, 0.0001  # ~0.4 arcsec for point source
        else:
            # Planet
            ra, dec, dist, angular_radius = _get_planet_position(jd, occulted_planet)
            return ra, dec, angular_radius

    def _angular_separation(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
        """Calculate angular separation between two points (in degrees)."""
        ra1_r = math.radians(ra1)
        dec1_r = math.radians(dec1)
        ra2_r = math.radians(ra2)
        dec2_r = math.radians(dec2)

        cos_sep = math.sin(dec1_r) * math.sin(dec2_r) + math.cos(dec1_r) * math.cos(
            dec2_r
        ) * math.cos(ra1_r - ra2_r)
        cos_sep = max(-1.0, min(1.0, cos_sep))
        return math.degrees(math.acos(cos_sep))

    def _check_occultation(jd: float) -> Tuple[bool, float, float, float]:
        """
        Check if occultation occurs at given time.

        Returns: (is_occultation, separation, occulting_radius, target_radius)
        """
        occ_ra, occ_dec, occ_dist, occ_radius = _get_planet_position(
            jd, occulting_planet
        )
        target_ra, target_dec, target_radius = _get_target_position(jd)

        separation = _angular_separation(occ_ra, occ_dec, target_ra, target_dec)

        # Occultation threshold: occulting planet radius + target radius
        threshold = occ_radius + target_radius

        return separation <= threshold, separation, occ_radius, target_radius

    def _find_minimum_separation(jd_start_search: float, jd_end_search: float) -> float:
        """Find time of minimum separation using golden section search."""
        phi = (1 + math.sqrt(5)) / 2  # Golden ratio

        a = jd_start_search
        b = jd_end_search

        c = b - (b - a) / phi
        d = a + (b - a) / phi

        def get_sep(jd: float) -> float:
            occ_ra, occ_dec, _, _ = _get_planet_position(jd, occulting_planet)
            target_ra, target_dec, _ = _get_target_position(jd)
            return _angular_separation(occ_ra, occ_dec, target_ra, target_dec)

        fc = get_sep(c)
        fd = get_sep(d)

        for _ in range(50):  # Converge to ~0.1 second precision
            if fc < fd:
                b = d
                d = c
                fd = fc
                c = b - (b - a) / phi
                fc = get_sep(c)
            else:
                a = c
                c = d
                fc = fd
                d = a + (b - a) / phi
                fd = get_sep(d)

            if b - a < 1e-6:  # ~0.1 second precision
                break

        return (a + b) / 2

    def _calculate_contact_times(
        jd_max: float, min_sep: float, occ_radius: float, target_radius: float
    ) -> Tuple[float, float, float, float]:
        """Calculate contact times for the occultation."""
        outer_threshold = occ_radius + target_radius
        inner_threshold = abs(occ_radius - target_radius)

        # Estimate relative angular speed
        dt_test = 1.0 / 24.0  # 1 hour
        occ_ra1, occ_dec1, _, _ = _get_planet_position(
            jd_max - dt_test, occulting_planet
        )
        occ_ra2, occ_dec2, _, _ = _get_planet_position(
            jd_max + dt_test, occulting_planet
        )
        target_ra1, target_dec1, _ = _get_target_position(jd_max - dt_test)
        target_ra2, target_dec2, _ = _get_target_position(jd_max + dt_test)

        # Relative motion
        d_ra = (occ_ra2 - occ_ra1) - (target_ra2 - target_ra1)
        d_dec = (occ_dec2 - occ_dec1) - (target_dec2 - target_dec1)
        relative_speed = math.sqrt(d_ra**2 + d_dec**2) / (2 * dt_test)

        if relative_speed < 0.001:  # Very slow motion, use fallback
            relative_speed = 0.5  # degrees/day fallback

        # Contact times based on chord geometry
        if min_sep < outer_threshold:
            half_chord_outer = math.sqrt(max(0, outer_threshold**2 - min_sep**2))
            dt_outer = half_chord_outer / relative_speed

            jd_first = jd_max - dt_outer
            jd_fourth = jd_max + dt_outer
        else:
            jd_first = 0.0
            jd_fourth = 0.0

        # For total occultation
        if min_sep < inner_threshold:
            half_chord_inner = math.sqrt(max(0, inner_threshold**2 - min_sep**2))
            dt_inner = half_chord_inner / relative_speed

            jd_second = jd_max - dt_inner
            jd_third = jd_max + dt_inner
        else:
            jd_second = 0.0
            jd_third = 0.0

        return jd_first, jd_second, jd_third, jd_fourth

    # Main search loop
    jd = jd_start
    step = 1.0  # Check every day initially

    for iteration in range(MAX_ITERATIONS):
        try:
            is_occ, sep, occ_r, target_r = _check_occultation(jd)
        except (ValueError, ArithmeticError, IndexError) as e:
            # Handle ephemeris range errors or other calculation errors
            if "ephemeris" in str(e).lower() or "range" in str(e).lower():
                # Exceeded ephemeris range - stop searching
                break
            raise

        if is_occ:
            # Found an occultation! Refine timing
            jd_max = _find_minimum_separation(jd - 1.0, jd + 1.0)

            # Verify it's still an occultation
            is_occ_refined, min_sep, occ_r, target_r = _check_occultation(jd_max)

            if is_occ_refined:
                # Calculate contact times
                jd_first, jd_second, jd_third, jd_fourth = _calculate_contact_times(
                    jd_max, min_sep, occ_r, target_r
                )

                # Determine occultation type
                # Grazing threshold: when the target passes within the outer 10%
                # of the occulting planet's disc
                grazing_threshold = 0.9 * occ_r
                is_grazing = min_sep > grazing_threshold

                if min_sep < abs(occ_r - target_r):
                    ecl_type = ECL_TOTAL
                else:
                    ecl_type = ECL_PARTIAL


                tret = (
                    jd_max,  # [0] Time of maximum
                    0.0,  # [1] Reserved
                    jd_first,  # [2] Begin
                    jd_fourth,  # [3] End
                    jd_second,  # [4] Totality begin
                    jd_third,  # [5] Totality end
                    0.0,  # [6] Reserved
                    0.0,  # [7] Reserved
                    0.0,  # [8] Reserved
                    0.0,  # [9] Reserved
                )

                return ecl_type, tret

        # Adaptive step based on angular separation
        if sep < 0.5:  # Very close - within half a degree
            step = 0.01  # Check every ~15 minutes
        elif sep < 2.0:  # Close - within 2 degrees
            step = 0.1
        elif sep < 10.0:
            step = 0.5
        else:
            step = 1.0

        jd += step

        if jd > jd_start + MAX_SEARCH_YEARS * 365.25:
            break

    if occulted_planet == 0:
        target_desc = starname
    else:
        target_desc = f"planet {occulted_planet}"
    occ_desc = f"planet {occulting_planet}"

    raise RuntimeError(
        f"No planetary occultation of {target_desc} by {occ_desc} found within "
        f"{MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for reference API compatibility
planet_occult_when_glob = planet_occult_when_glob


def planet_occult_when_loc(
    jd_start: float,
    occulting_planet: int,
    occulted_planet: int = 0,
    star_name: str = "",
    lat: float = 0.0,
    lon: float = 0.0,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """Find the next planetary occultation visible from a specific location.

    Wrapper around :func:`_planet_occult_when_loc_impl` that adds LEB→Skyfield
    fallback for partial/custom LEB files. See the impl docstring for the
    full API contract.
    """
    return _call_with_leb_skyfield_fallback(
        _planet_occult_when_loc_impl,
        jd_start, occulting_planet, occulted_planet, star_name,
        lat, lon, altitude, flags,
    )


def _planet_occult_when_loc_impl(
    jd_start: float,
    occulting_planet: int,
    occulted_planet: int = 0,
    star_name: str = "",
    lat: float = 0.0,
    lon: float = 0.0,
    altitude: float = 0.0,
    flags: int = FLG_SWIEPH,
    *,
    reader=None,
) -> Tuple[int, Tuple[float, ...], Tuple[float, ...]]:
    """
    Find the next planetary occultation visible from a specific location.

    A planetary occultation occurs when one planet passes in front of (occults)
    another planet or star. This function searches forward in time to find the
    next occultation visible from a specific geographic location, where both
    the occulting and occulted bodies are above the horizon.

    Args:
        jd_start: Julian Day (UT) to start search from
        occulting_planet: Planet ID of the occulting (foreground) planet.
            Use VENUS, MARS, JUPITER, etc.
        occulted_planet: Planet ID of the occulted (background) planet.
            Set to 0 if searching for a fixed star occultation.
        star_name: Name of fixed star to check (e.g. "Regulus", "Spica").
            Only used if occulted_planet is 0.
        lat: Observer latitude in degrees (positive = North, negative = South)
        lon: Observer longitude in degrees (positive = East, negative = West)
        altitude: Observer altitude in meters above sea level (default 0)
        flags: Calculation flags (FLG_SWIEPH, etc.)

    Returns:
        Tuple containing:
            - times: Tuple of 10 floats with occultation phase times (JD UT):
                [0]: Time of maximum occultation (minimum separation)
                [1]: Time of first contact (occultation begins)
                [2]: Time of second contact (full occultation begins, or 0)
                [3]: Time of third contact (full occultation ends, or 0)
                [4]: Time of fourth contact (occultation ends)
                [5-9]: Reserved (0)
            - attr: Tuple of 20 floats with occultation attributes:
                [0]: Fraction of target covered (magnitude)
                [1]: Ratio of occulting to occulted angular diameter
                [2]: Fraction of disc covered (obscuration)
                [3]: Reserved (0)
                [4]: Azimuth of occulted body at maximum (degrees)
                [5]: True altitude of occulted body at maximum (degrees)
                [6]: Apparent altitude (with refraction)
                [7]: Angular separation at maximum (degrees)
                [8-19]: Reserved (0)
            - retflag: Occultation type flags (ECL_* constants)

    Raises:
        RuntimeError: If no occultation visible from location within search limit
        ValueError: If neither occulted_planet nor star_name is specified

    Example:
        >>> # Find next occultation of Regulus by Venus visible from Rome
        >>> from libephemeris import julday, planet_occult_when_loc, VENUS
        >>> jd = julday(2020, 1, 1, 0)
        >>> rome_lat, rome_lon = 41.9028, 12.4964
        >>> ecl_type, times, attr = planet_occult_when_loc(
        ...     jd, VENUS, 0, "Regulus", rome_lat, rome_lon
        ... )
        >>> print(f"Occultation max at JD {times[0]:.5f}")

    References:
        - Meeus "Astronomical Algorithms" Ch. 9 (Angular Separation)
    """
    from skyfield.api import wgs84

    from .constants import (
        MERCURY,
        VENUS,
        MARS,
        JUPITER,
        SATURN,
        URANUS,
        NEPTUNE,
        PLUTO,
        SUN,
        MOON,
    )
    from .fixed_stars import FIXED_STARS, _resolve_star_id
    from .planets import _PLANET_MAP, get_planet_target
    from .state import get_planets, get_timescale

    if occulted_planet == 0 and not star_name:
        raise ValueError(
            "Either occulted_planet ID or star_name must be specified for occultation search"
        )

    # Validate planets
    if occulting_planet == SUN:
        raise ValueError("Sun cannot be an occulting body - use _sol_eclipse_when_loc_pythonic")
    if occulting_planet == MOON:
        raise ValueError("Moon cannot be an occulting body - use _lun_occult_when_loc_pythonic")
    if occulting_planet not in _PLANET_MAP:
        raise ValueError(f"Invalid occulting planet ID: {occulting_planet}")

    if occulted_planet != 0:
        if occulted_planet not in _PLANET_MAP:
            raise ValueError(f"Invalid occulted planet ID: {occulted_planet}")
        if occulted_planet == occulting_planet:
            raise ValueError("Occulting and occulted planets cannot be the same")

    MAX_SEARCH_YEARS = 150
    MAX_GLOBAL_SEARCHES = 100

    # reader is provided by the caller (None forces Skyfield path)
    _reader = reader
    if _reader is not None:
        from .fast_calc import _topo_ecliptic
        from .time_utils import deltat
        from .utils import azalt, ECL2HOR

        geopos = (lon, lat, altitude)

        def _get_body_altitude(jd: float, planet_id: int) -> Tuple[float, float]:
            jd_tt = jd + deltat(jd)
            pos = _topo_ecliptic(_reader, jd_tt, jd, planet_id, geopos, FLG_SPEED)
            az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, pos[:3])
            return alt_true, az

        def _get_target_altitude(jd: float) -> Tuple[float, float]:
            jd_tt = jd + deltat(jd)
            if occulted_planet == 0:
                from .fixed_stars import fixstar_ut
                star_pos, _, _ = fixstar_ut(star_name, jd, FLG_SPEED)
                az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, star_pos[:3])
            else:
                pos = _topo_ecliptic(_reader, jd_tt, jd, occulted_planet, geopos, FLG_SPEED)
                az, alt_true, _alt_app = azalt(jd, ECL2HOR, geopos, 0, 0, pos[:3])
            return alt_true, az
    else:
        eph = get_planets()
        ts = get_timescale()
        earth = eph["earth"]
        observer = wgs84.latlon(lat, lon, altitude)
        observer_at = earth + observer

        def _get_body_altitude(jd: float, planet_id: int) -> Tuple[float, float]:
            from .planets import get_planet_target
            target_name = _PLANET_MAP[planet_id]
            target = get_planet_target(eph, target_name)
            t = ts.ut1_jd(jd)
            target_app = observer_at.at(t).observe(target).apparent()
            alt, az, _ = target_app.altaz()
            return alt.degrees, (az.degrees + 180.0) % 360.0

        def _get_target_altitude(jd: float) -> Tuple[float, float]:
            from .planets import get_planet_target
            t = ts.ut1_jd(jd)
            if occulted_planet == 0:
                star_id, err, _ = _resolve_star_id(star_name)
                if err is not None:
                    raise ValueError(err)
                star = FIXED_STARS[star_id]
                t_years = (jd - 2451545.0) / 365.25
                ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
                from skyfield.api import Star
                star_obj = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
                target_app = observer_at.at(t).observe(star_obj).apparent()
            else:
                target_name = _PLANET_MAP[occulted_planet]
                target = get_planet_target(eph, target_name)
                target_app = observer_at.at(t).observe(target).apparent()
            alt, az, _ = target_app.altaz()
            # Convert Skyfield navigational azimuth (N=0) to SE convention (S=0)
            return alt.degrees, (az.degrees + 180.0) % 360.0

    def _is_visible_at_location(jd: float) -> bool:
        """Check if both bodies are above horizon at given time."""
        occ_alt, _ = _get_body_altitude(jd, occulting_planet)
        target_alt, _ = _get_target_altitude(jd)

        # Both bodies must be above horizon (with some margin for twilight)
        MIN_ALT = -0.5  # Allow slightly below horizon for refraction
        return occ_alt > MIN_ALT and target_alt > MIN_ALT

    current_jd = jd_start

    for search_count in range(MAX_GLOBAL_SEARCHES):
        try:
            # Find next global occultation
            retflags, tret = planet_occult_when_glob(
                current_jd, occulting_planet, occulted_planet, star_name, flags, 0
            )

            jd_max = tret[0]
            jd_begin = tret[2]
            jd_end = tret[3]

            # Check if any part of the occultation is visible from this location
            # Sample multiple times during the occultation
            visible = False
            visible_times = []

            if jd_begin > 0 and jd_end > 0:
                sample_count = 10
                for i in range(sample_count + 1):
                    sample_jd = jd_begin + (jd_end - jd_begin) * i / sample_count
                    if _is_visible_at_location(sample_jd):
                        visible = True
                        visible_times.append(sample_jd)

            # Also check at maximum
            if _is_visible_at_location(jd_max):
                visible = True
                visible_times.append(jd_max)

            if visible:
                occ_alt, occ_az = _get_body_altitude(jd_max, occulting_planet)
                target_alt, target_az = _get_target_altitude(jd_max)

                if _reader is not None:
                    from .fast_calc import _topo_ecliptic as _tp_occ
                    from .time_utils import deltat as _sd_occ
                    from .utils import angular_separation as _ang_sep_occ

                    _jd_tt_occ = jd_max + _sd_occ(jd_max)
                    _gp_occ = (lon, lat, altitude)
                    occ_pos = _tp_occ(_reader, _jd_tt_occ, jd_max, occulting_planet, _gp_occ, FLG_SPEED)
                    if occulted_planet == 0:
                        from .fixed_stars import fixstar_ut
                        tgt_pos, _, _ = fixstar_ut(star_name, jd_max, FLG_SPEED)
                        target_radius = 0.0001
                    else:
                        tgt_pos = _tp_occ(_reader, _jd_tt_occ, jd_max, occulted_planet, _gp_occ, FLG_SPEED)
                        target_radius = _calc_planet_angular_radius(occulted_planet, tgt_pos[2])
                    separation = _ang_sep_occ(occ_pos[0], occ_pos[1], tgt_pos[0], tgt_pos[1])
                    occ_dist = occ_pos[2]
                    occ_radius = _calc_planet_angular_radius(occulting_planet, occ_dist)
                else:
                    from .planets import get_planet_target

                    t = ts.ut1_jd(jd_max)
                    target_name = _PLANET_MAP[occulting_planet]
                    occ_body = get_planet_target(eph, target_name)
                    occ_app = observer_at.at(t).observe(occ_body).apparent()
                    if occulted_planet == 0:
                        star_id, _, _ = _resolve_star_id(star_name)
                        star = FIXED_STARS[star_id]
                        t_years = (jd_max - 2451545.0) / 365.25
                        ra_deg = star.ra_j2000 + (star.pm_ra * t_years) / 3600.0
                        dec_deg = star.dec_j2000 + (star.pm_dec * t_years) / 3600.0
                        from skyfield.api import Star
                        target_body = Star(ra_hours=ra_deg / 15.0, dec_degrees=dec_deg)
                    else:
                        target_body = get_planet_target(eph, _PLANET_MAP[occulted_planet])
                    target_app = observer_at.at(t).observe(target_body).apparent()
                    separation = occ_app.separation_from(target_app).degrees
                    occ_dist = occ_app.distance().au
                    occ_radius = _calc_planet_angular_radius(occulting_planet, occ_dist)
                    if occulted_planet == 0:
                        target_radius = 0.0001
                    else:
                        target_dist = target_app.distance().au
                        target_radius = _calc_planet_angular_radius(
                            occulted_planet, target_dist
                        )

                if target_radius > 0:
                    magnitude = max(0, 1.0 - separation / (occ_radius + target_radius))
                    ratio = (
                        occ_radius / target_radius if target_radius > 0.0001 else 999.0
                    )
                    obscuration = magnitude  # Simplified
                else:
                    magnitude = 1.0 if separation < occ_radius else 0.0
                    ratio = 999.0
                    obscuration = magnitude

                # Refraction correction (simplified)
                apparent_alt = (
                    target_alt
                    + (
                        1.02
                        / math.tan(
                            math.radians(target_alt + 10.3 / (target_alt + 5.11))
                        )
                    )
                    / 60.0
                    if target_alt > -1
                    else target_alt
                )

                times = (
                    jd_max,  # [0] Maximum
                    jd_begin if jd_begin > 0 else 0.0,  # [1] First contact
                    tret[4] if tret[4] > 0 else 0.0,  # [2] Second contact
                    tret[5] if tret[5] > 0 else 0.0,  # [3] Third contact
                    jd_end if jd_end > 0 else 0.0,  # [4] Fourth contact
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,  # [5-9] Reserved
                )

                attr = (
                    magnitude,  # [0] Magnitude
                    ratio,  # [1] Diameter ratio
                    obscuration,  # [2] Obscuration
                    0.0,  # [3] Reserved
                    target_az,  # [4] Azimuth
                    target_alt,  # [5] True altitude
                    apparent_alt,  # [6] Apparent altitude
                    separation,  # [7] Separation
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,  # [8-19] Reserved
                )

                return retflags, times, attr

            # Not visible from this location - continue search after this event
            current_jd = jd_end + 1.0 if jd_end > 0 else jd_max + 1.0

        except RuntimeError:
            # No more global occultations found
            break

    if occulted_planet == 0:
        target_desc = star_name
    else:
        target_desc = f"planet {occulted_planet}"
    occ_desc = f"planet {occulting_planet}"

    raise RuntimeError(
        f"No planetary occultation of {target_desc} by {occ_desc} visible from "
        f"lat={lat}, lon={lon} found within {MAX_SEARCH_YEARS} years of JD {jd_start}"
    )


# Alias for reference API compatibility
planet_occult_when_loc = planet_occult_when_loc
