# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Fast calculation pipeline using precomputed .leb binary ephemeris data.

This module reimplements the Skyfield pipeline using LEBReader as the data
source, handling coordinate transforms, light-time correction, aberration,
and flag dispatch.

Three pipelines:
    A: ICRS barycentric bodies (Sun-Pluto, Earth, Chiron, Ceres-Vesta)
    B: Ecliptic direct bodies (lunar nodes, Lilith variants)
    C: Independently sourced heliocentric fictitious bodies

Provenance:
    Stored coefficients are generated from the registered NASA JPL/IAU source
    pipeline and evaluated with the documented LEB Chebyshev representation.
    Light-time, deflection, aberration, precession/nutation, frame, center, and
    flag ordering follow the public IAU/IERS and JPL conventions cited in
    ``docs/leb/algorithms.md``. Dispatch order, cache boundaries, and finite-
    difference steps are project engineering and are documented beside their
    use. This module contains no independently fitted astronomical table.
"""

from __future__ import annotations

import math
import threading
from functools import lru_cache
from typing import TYPE_CHECKING, Callable, Dict, Optional, Tuple

from .constants import (
    EARTH,
    FICT_OFFSET,
    INTP_APOG,
    INTP_PERG,
    JUPITER,
    MARS,
    MEAN_APOG,
    MEAN_NODE,
    MERCURY,
    MOON,
    NEPTUNE,
    OSCU_APOG,
    PLUTO,
    SATURN,
    SUN,
    TRUE_NODE,
    URANUS,
    VENUS,
    WALDEMATH,
    FLG_BARYCTR,
    FLG_EQUATORIAL,
    FLG_HELCTR,
    FLG_ICRS,
    FLG_J2000,
    FLG_MOSEPH,
    FLG_NOABERR,
    FLG_NOGDEFL,
    FLG_NONUT,
    FLG_RADIANS,
    FLG_SIDEREAL,
    FLG_SPEED,
    FLG_SPEED3,
    FLG_TOPOCTR,
    FLG_TRUEPOS,
    FLG_XYZ,
)

from .leb_format import (
    COORD_ECLIPTIC,
    COORD_HELIO_ECL,
    COORD_ICRS_BARY,
    COORD_ICRS_BARY_SYSTEM,
)
from .precession_vondrak import (
    method_b_accumulated_precession,
    vondrak_mean_obliquity_deg,
    vondrak_pn_matrix,
    vondrak_precession_matrix,
)
from .ayanamsha_definitions import (
    AYANAMSHA_DEFINING,
    DYNAMIC_AYANAMSHA_MODES,
)

# These planet IDs were historically importable from ``libephemeris.fast_calc``.
# Keep the module-level surface even though the current dispatcher no longer
# needs each name directly.
_LEGACY_REEXPORTED_PLANET_IDS = (
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
)

if TYPE_CHECKING:
    from typing import Protocol, Union

    from .leb2_reader import LEB2Reader
    from .leb_composite import CompositeLEBReader, TieredLEBReader
    from .leb_reader import LEBReader

    LEBReaderLike = Union[LEBReader, LEB2Reader, CompositeLEBReader, TieredLEBReader]

    class _DeflectorSource(Protocol):
        """Structural type for anything that can supply deflector states.

        Only ``eval_body`` is required. The parameter is positional-only so
        both the LEB readers (``jd``) and the Skyfield adapter (``jd_tt``)
        satisfy it regardless of their keyword name.
        """

        def eval_body(
            self, body_id: int, jd: float, /
        ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]: ...

# =============================================================================
# CONSTANTS
# =============================================================================

# Speed of light in AU/day (NOVAS/Skyfield ``C_AUDAY``; adopted for
# cross-backend bit-identity). c = 299792458 m/s exactly (BIPM/CODATA SI) with
# the IAU 1976 System au (1.49597870691e11 m) reproduces 173.1446326846693.
C_LIGHT_AU_DAY = 173.1446326846693  # Speed of light in AU/day
J2000 = 2451545.0  # J2000.0 epoch in JD

# SPICE/ECLIPJ2000 obliquity (IAU 1976, 84381.448"). Used ONLY as the
# internal bridge frame for reviewed analytical heliocentric orbits: the
# stored/native theory frame and the Earth-vector rotation
# must use the same fixed standard frame on both backends.
OBLIQUITY_J2000_DEG = 23.4392911
OBLIQUITY_J2000_RAD = math.radians(OBLIQUITY_J2000_DEG)

# IAU 2006 mean obliquity at J2000 (Capitaine et al. 2003; SOFA/ERFA obl06).
# All user-visible J2000 ecliptic/equatorial frame conversions use this
# independently standardized value.
OBLIQUITY_J2000_REF_DEG = 84381.406 / 3600.0
OBLIQUITY_J2000_REF_RAD = math.radians(OBLIQUITY_J2000_REF_DEG)


def _build_icrs_to_j2000_frames() -> Tuple[
    Tuple[Tuple[float, float, float], ...], Tuple[Tuple[float, float, float], ...]
]:
    """Build ICRS -> mean-equator/ecliptic J2000 frame matrices.

    ERFA ``bp06`` supplies the IAU 2006 frame-bias matrix. The ecliptic matrix
    then applies the IAU 2006 mean obliquity at J2000.
    """
    import erfa

    rb = erfa.bp06(2451545.0, 0.0)[0]
    bias = tuple((float(row[0]), float(row[1]), float(row[2])) for row in rb)
    ce = math.cos(OBLIQUITY_J2000_REF_RAD)
    se = math.sin(OBLIQUITY_J2000_REF_RAD)
    rx = ((1.0, 0.0, 0.0), (0.0, ce, se), (0.0, -se, ce))
    ecl = tuple(
        (
            sum(rx[i][k] * bias[k][0] for k in range(3)),
            sum(rx[i][k] * bias[k][1] for k in range(3)),
            sum(rx[i][k] * bias[k][2] for k in range(3)),
        )
        for i in range(3)
    )
    return bias, ecl


_FRAME_BIAS_J2000, _ICRS_TO_ECL_J2000_REF = _build_icrs_to_j2000_frames()

# IAU 2006 mean obliquity polynomial coefficients (arcseconds)
# eps = 84381.406 - 46.836769*T - 0.0001831*T^2 + 0.00200340*T^3 ...
_OBLIQUITY_COEFFS = (
    84381.406,
    -46.836769,
    -0.0001831,
    0.00200340,
    -0.000000576,
    -0.0000000434,
)


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================


def _mean_obliquity_iau2006(jd_tt: float) -> float:
    """Compute IAU 2006 mean obliquity in degrees.

    Args:
        jd_tt: Julian Day in TT.

    Returns:
        Mean obliquity in degrees.
    """
    T = (jd_tt - J2000) / 36525.0
    eps_arcsec = _OBLIQUITY_COEFFS[0]
    T_power = T
    for i in range(1, len(_OBLIQUITY_COEFFS)):
        eps_arcsec += _OBLIQUITY_COEFFS[i] * T_power
        T_power *= T
    return eps_arcsec / 3600.0


def _vec3_sub(a: Tuple[float, ...], b: Tuple[float, ...]) -> Tuple[float, float, float]:
    """Subtract two 3-vectors."""
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _vec3_dist(v: Tuple[float, float, float]) -> float:
    """Euclidean distance of a 3-vector."""
    return math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])


def _cartesian_to_spherical(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Convert Cartesian to spherical (lon, lat, dist) in degrees and AU.

    Returns:
        (longitude_deg, latitude_deg, distance_au)
    """
    dist = math.sqrt(x * x + y * y + z * z)
    if dist == 0.0:
        return (0.0, 0.0, 0.0)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(max(-1.0, min(1.0, z / dist))))
    return (lon, lat, dist)


def _cartesian_velocity_to_spherical(
    x: float,
    y: float,
    z: float,
    vx: float,
    vy: float,
    vz: float,
) -> Tuple[float, float, float]:
    """Convert Cartesian velocity to spherical velocity.

    Given position (x, y, z) and velocity (vx, vy, vz), compute the time
    derivatives of (longitude, latitude, distance) using the standard
    Cartesian-to-spherical Jacobian.

    Args:
        x, y, z: Position in AU (Cartesian).
        vx, vy, vz: Velocity in AU/day (Cartesian).

    Returns:
        (dlon_deg_day, dlat_deg_day, ddist_au_day)
    """
    r_xy_sq = x * x + y * y
    r_sq = r_xy_sq + z * z
    r = math.sqrt(r_sq)
    r_xy = math.sqrt(r_xy_sq)

    if r == 0.0 or r_xy == 0.0:
        # At the pole or origin — angular rates undefined
        ddist = math.sqrt(vx * vx + vy * vy + vz * vz) if r == 0.0 else 0.0
        return (0.0, 0.0, ddist)

    dlon_rad = (x * vy - y * vx) / r_xy_sq  # rad/day
    dlat_rad = (vz * r_xy_sq - z * (x * vx + y * vy)) / (r_sq * r_xy)  # rad/day
    ddist = (x * vx + y * vy + z * vz) / r  # AU/day

    dlon_deg = math.degrees(dlon_rad)  # deg/day
    dlat_deg = math.degrees(dlat_rad)  # deg/day

    return (dlon_deg, dlat_deg, ddist)


def _spherical_to_cartesian_with_velocity(
    lon: float,
    lat: float,
    dist: float,
    dlon: float,
    dlat: float,
    ddist: float,
) -> Tuple[float, float, float, float, float, float]:
    """Convert spherical position+velocity to Cartesian (the inverse Jacobian).

    Shared by the FLG_XYZ post-processing in both the LEB/Skyfield path
    (``fast_calc``) and the Horizons backend so the two stay identical.

    Args:
        lon, lat: Longitude/latitude in DEGREES.
        dist: Radial distance in AU.
        dlon, dlat: Longitude/latitude speed in DEGREES/day.
        ddist: Radial speed in AU/day.

    Returns:
        (x, y, z, vx, vy, vz) — position in AU, velocity in AU/day.  The
        velocity is exactly zero when all input rates are zero.
    """
    lon_r = math.radians(lon)
    lat_r = math.radians(lat)
    cos_lat = math.cos(lat_r)
    sin_lat = math.sin(lat_r)
    cos_lon = math.cos(lon_r)
    sin_lon = math.sin(lon_r)
    x = dist * cos_lat * cos_lon
    y = dist * cos_lat * sin_lon
    z = dist * sin_lat
    if dlon != 0.0 or dlat != 0.0 or ddist != 0.0:
        dlon_r = math.radians(dlon)
        dlat_r = math.radians(dlat)
        vx = (
            -dist * cos_lat * sin_lon * dlon_r
            - dist * sin_lat * cos_lon * dlat_r
            + cos_lat * cos_lon * ddist
        )
        vy = (
            dist * cos_lat * cos_lon * dlon_r
            - dist * sin_lat * sin_lon * dlat_r
            + cos_lat * sin_lon * ddist
        )
        vz = dist * cos_lat * dlat_r + sin_lat * ddist
    else:
        vx = vy = vz = 0.0
    # Cast to native floats: callers may pass numpy scalars (e.g. star vectors),
    # and several return this helper directly, so honour the native-float contract.
    return (float(x), float(y), float(z), float(vx), float(vy), float(vz))


_SPEED3_STEP_DAYS = 0.0001


def _apply_central_speed3_stencil(
    current: Tuple[float, float, float, float, float, float],
    previous: Tuple[float, float, float, float, float, float],
    following: Tuple[float, float, float, float, float, float],
    iflag: int,
) -> Tuple[float, float, float, float, float, float]:
    """Rebuild SPEED3 rates with a standard centered finite difference.

    The clean-room implementation uses ``(P(t+h) - P(t-h)) / (2*h)``.
    Longitude is unwrapped around the current position for spherical output;
    Cartesian output differentiates each axis independently.
    """
    if iflag & FLG_XYZ:
        velocity = tuple(
            (following[index] - previous[index]) / (2.0 * _SPEED3_STEP_DAYS)
            for index in range(3)
        )
    else:
        half_cycle = math.pi if iflag & FLG_RADIANS else 180.0
        full_cycle = 2.0 * half_cycle
        previous_delta = (
            previous[0] - current[0] + half_cycle
        ) % full_cycle - half_cycle
        following_delta = (
            following[0] - current[0] + half_cycle
        ) % full_cycle - half_cycle
        velocity = (
            (following_delta - previous_delta) / (2.0 * _SPEED3_STEP_DAYS),
            (following[1] - previous[1]) / (2.0 * _SPEED3_STEP_DAYS),
            (following[2] - previous[2]) / (2.0 * _SPEED3_STEP_DAYS),
        )
    return (
        float(current[0]),
        float(current[1]),
        float(current[2]),
        float(velocity[0]),
        float(velocity[1]),
        float(velocity[2]),
    )


def _rotate_equatorial_to_ecliptic(
    x: float, y: float, z: float, eps_rad: float
) -> Tuple[float, float, float]:
    """Rotate from equatorial to ecliptic frame.

    Args:
        x, y, z: Equatorial Cartesian coordinates.
        eps_rad: Obliquity in radians.

    Returns:
        (x_ecl, y_ecl, z_ecl) in ecliptic frame.
    """
    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    return (
        x,
        y * cos_eps + z * sin_eps,
        -y * sin_eps + z * cos_eps,
    )


def _rotate_icrs_to_ecliptic_j2000(
    x: float, y: float, z: float
) -> Tuple[float, float, float]:
    """Rotate ICRS to IAU 2006 mean-ecliptic J2000 coordinates."""
    m = _ICRS_TO_ECL_J2000_REF
    return (
        m[0][0] * x + m[0][1] * y + m[0][2] * z,
        m[1][0] * x + m[1][1] * y + m[1][2] * z,
        m[2][0] * x + m[2][1] * y + m[2][2] * z,
    )


def _rotate_icrs_to_j2000_mean_equator(
    x: float, y: float, z: float
) -> Tuple[float, float, float]:
    """Rotate ICRS to the mean equator/equinox of J2000 (frame bias only)."""
    m = _FRAME_BIAS_J2000
    return (
        m[0][0] * x + m[0][1] * y + m[0][2] * z,
        m[1][0] * x + m[1][1] * y + m[1][2] * z,
        m[2][0] * x + m[2][1] * y + m[2][2] * z,
    )


def _apply_aberration(
    geo: Tuple[float, float, float],
    earth_vel: Tuple[float, float, float],
    light_time: float = 0.0,
) -> Tuple[float, float, float]:
    """Apply relativistic aberration to a geometric position vector.

    Special-relativistic stellar aberration (Explanatory Supplement to the
    Astronomical Almanac, 3rd ed. 2013, section 7.2.3; Kaplan 2005, USNO
    Circular 179, section 6). Uses the full special-relativistic formula
    implemented by Skyfield's ``add_aberration()``, including the Lorentz
    factor gamma.

    When *light_time* is zero or negative, falls back to a first-order
    approximation (classical Bradley aberration) for backward
    compatibility with callers that do not track light-time.

    Args:
        geo: Light-time-corrected geocentric position vector (ICRS, AU).
        earth_vel: Earth barycentric velocity vector (ICRS, AU/day).
        light_time: One-way light travel time in days (observer→target).

    Returns:
        Aberrated position vector (same frame, AU).
    """
    dist = _vec3_dist(geo)
    if dist == 0.0:
        return geo

    # ── Full special-relativistic aberration (Skyfield formula) ──────────
    if light_time > 0.0:
        p1mag = light_time * C_LIGHT_AU_DAY  # distance in AU
        vemag = math.sqrt(earth_vel[0] ** 2 + earth_vel[1] ** 2 + earth_vel[2] ** 2)
        if vemag == 0.0 or p1mag == 0.0:
            return geo

        beta = vemag / C_LIGHT_AU_DAY
        dot = geo[0] * earth_vel[0] + geo[1] * earth_vel[1] + geo[2] * earth_vel[2]
        cosd = dot / (p1mag * vemag)
        gammai = math.sqrt(1.0 - beta * beta)  # inverse Lorentz factor
        p = beta * cosd
        q = (1.0 + p / (1.0 + gammai)) * light_time
        r = 1.0 + p
        if abs(r) < 1e-30:
            r = 1.0

        return (
            (gammai * geo[0] + q * earth_vel[0]) / r,
            (gammai * geo[1] + q * earth_vel[1]) / r,
            (gammai * geo[2] + q * earth_vel[2]) / r,
        )

    # ── Fallback: first-order Bradley aberration ────────────────────────
    ux = geo[0] / dist
    uy = geo[1] / dist
    uz = geo[2] / dist

    vx = earth_vel[0] / C_LIGHT_AU_DAY
    vy = earth_vel[1] / C_LIGHT_AU_DAY
    vz = earth_vel[2] / C_LIGHT_AU_DAY

    dot = ux * vx + uy * vy + uz * vz

    ax = ux + vx - ux * dot
    ay = uy + vy - uy * dot
    az = uz + vz - uz * dot

    a_dist = math.sqrt(ax * ax + ay * ay + az * az)
    if a_dist == 0.0:
        return geo

    return (ax / a_dist * dist, ay / a_dist * dist, az / a_dist * dist)


# Gravitational deflection constants (matching Skyfield's relativity module so
# the deflection is bit-identical to the Skyfield backend).
# Heliocentric gravitational constant GM_sun (m^3/s^2). This is the value
# carried by the NOVAS/Skyfield relativity model (equal to Skyfield's ``GS``);
# it is NOT the DE440 GM_sun of Park, Folkner, Williams & Boggs (2021), AJ 161,
# 105 (1.327124400413e20). It is retained here only to reproduce Skyfield's
# deflection bit-for-bit.
_GS = 1.32712440017987e20  # heliocentric gravitational constant (m^3/s^2)
_C_MS = 299792458.0  # speed of light (m/s), exact SI (BIPM/CODATA)
_AU_M = 149597870700  # 1 AU in metres (IAU 2012 Resolution B2)

# Deflector reciprocal masses (solar mass / deflector mass). Sun-to-planet mass
# ratios of the order used in the Astronomical Almanac (Selected Astronomical
# Constants), matching the deflector set of the Skyfield relativity model.
_DEFLECTORS: Tuple[Tuple[int, float], ...] = (
    (SUN, 1.0),  # Sun
    (5, 1047.3486),  # Jupiter barycenter
    (6, 3497.898),  # Saturn barycenter
)

# Body IDs in the LEB file for deflector barycenters
# (Jupiter=5, Saturn=6 map to body_ids 5, 6 in LEB — their barycentric ICRS)
_DEFLECTOR_LEB_IDS = {SUN: SUN, 5: 5, 6: 6}

# =============================================================================
# CENTER-OF-BODY (COB) CORRECTION FOR SYSTEM BARYCENTERS
# =============================================================================

# Barycenter names for outer planets (used by get_cob_offset / SPK centers)
_SYSTEM_BARY_NAMES: Dict[int, str] = {
    5: "jupiter barycenter",
    6: "saturn barycenter",
    7: "uranus barycenter",
    8: "neptune barycenter",
    9: "pluto barycenter",
}

# Bodies stored as a system barycentre (COORD_ICRS_BARY_SYSTEM): the giant
# planets that can carry a runtime JPL center-of-body offset
# (barycentre -> planet centre).
_SYSTEM_BARY_BODIES = frozenset(_SYSTEM_BARY_NAMES)


def _apply_cob_correction(
    pos: Tuple[float, float, float],
    ipl: int,
    jd_tt: float,
) -> Tuple[float, float, float]:
    """Apply center-of-body correction to a system barycenter position.

    Converts system barycenter ICRS position to planet center ICRS position
    by adding the offset from a JPL planet-center SPK segment when one is
    available. If no such segment covers the epoch, the system barycentre is
    returned unchanged.

    Args:
        pos: System barycenter ICRS position (x, y, z) in AU.
        ipl: Body ID (5=Jupiter, 6=Saturn, 7=Uranus, 8=Neptune, 9=Pluto).
        jd_tt: Julian Day TT at which to evaluate the COB offset.

    Returns:
        Planet center ICRS position (x, y, z) in AU.
    """
    from .cache import get_cached_time_tt
    from .planets import _PLANET_CENTER_NAIF_IDS
    from .state import get_planet_center_segment

    bary_name = _SYSTEM_BARY_NAMES.get(ipl)
    if bary_name is None:
        return pos  # This body does not need a system-to-center offset.

    t = get_cached_time_tt(jd_tt)

    # Map body_id to planet name for NAIF lookup
    planet_name = {5: "jupiter", 6: "saturn", 7: "uranus", 8: "neptune", 9: "pluto"}[
        ipl
    ]

    # Try SPK center offset first (high precision)
    if planet_name in _PLANET_CENTER_NAIF_IDS:
        naif_id = _PLANET_CENTER_NAIF_IDS[planet_name]
        seg = get_planet_center_segment(naif_id, jd_tt)
        if seg is not None:
            try:
                offset_pos = seg.at(t).position.au
                return (
                    pos[0] + float(offset_pos[0]),
                    pos[1] + float(offset_pos[1]),
                    pos[2] + float(offset_pos[2]),
                )
            except (ValueError, ArithmeticError):
                pass

    # A system barycentre is the explicit fallback when no independently
    # sourced JPL center segment is available for this epoch.
    return pos


def _apply_gravitational_deflection(
    geo: Tuple[float, float, float],
    earth_bary: Tuple[float, float, float],
    jd_tt: float,
    light_time: float,
    reader: "_DeflectorSource",
) -> Tuple[float, float, float]:
    """Apply PPN gravitational light deflection by Sun, Jupiter, Saturn.

    Post-Newtonian (PPN, gamma = 1) light deflection in the field of each
    deflector (Explanatory Supplement to the Astronomical Almanac, 3rd ed.
    2013, section 7.2.4; Kaplan 2005, USNO Circular 179, section 6). Matches
    Skyfield's ``apparent(deflectors=(10, 599, 699))`` formula. Uses LEB data
    for deflector positions (barycentric ICRS).

    For the Sun, deflection is the dominant correction (~max 1.75" at
    the limb, typically 0.01–4" for planets).  Jupiter and Saturn add
    ~0.01" near their limbs.

    Args:
        geo: Light-time-corrected geocentric ICRS position (AU).
        earth_bary: Earth ICRS barycentric position (AU).
        jd_tt: Observation Julian Day in TT.
        light_time: Light travel time to target (days).
        reader: Open LEBReader for deflector positions.

    Returns:
        Deflection-corrected geocentric ICRS position (AU).
    """
    result = list(geo)
    pmag = math.sqrt(result[0] ** 2 + result[1] ** 2 + result[2] ** 2)
    if pmag == 0.0:
        return geo

    for defl_body_id, rmass in _DEFLECTORS:
        # 1. Deflector barycentric position at observation time
        try:
            defl_pos, _ = reader.eval_body(defl_body_id, jd_tt)
        except (KeyError, ValueError):
            continue

        # 2. Deflector relative to observer
        gpv = (
            defl_pos[0] - earth_bary[0],
            defl_pos[1] - earth_bary[1],
            defl_pos[2] - earth_bary[2],
        )

        # 3. Unit vector to target
        phat = (result[0] / pmag, result[1] / pmag, result[2] / pmag)

        # 4. Light-time difference: when photon passed closest to deflector
        dlt = (phat[0] * gpv[0] + phat[1] * gpv[1] + phat[2] * gpv[2]) / C_LIGHT_AU_DAY

        # 5. Clamp and compute time at closest approach
        tclose_offset = max(0.0, min(dlt, light_time))
        tclose_jd = jd_tt - tclose_offset

        # 6. Deflector position at closest approach
        try:
            defl_close, _ = reader.eval_body(defl_body_id, tclose_jd)
        except (KeyError, ValueError):
            continue

        # 7. pe = observer - deflector (observer relative to deflector)
        pe = (
            earth_bary[0] - defl_close[0],
            earth_bary[1] - defl_close[1],
            earth_bary[2] - defl_close[2],
        )

        # 8. pq = target relative to deflector (from observer frame)
        pq = (result[0] + pe[0], result[1] + pe[1], result[2] + pe[2])

        qmag = math.sqrt(pq[0] ** 2 + pq[1] ** 2 + pq[2] ** 2)
        emag = math.sqrt(pe[0] ** 2 + pe[1] ** 2 + pe[2] ** 2)

        # Skip a deflector the observer is inside of or grazing (emag below
        # ~2 solar radii). This happens for planet-centric observers at the
        # deflecting body itself (e.g. calc_pctr from Jupiter/Saturn, where
        # pe is just the centre-vs-system-barycentre offset, ~1e-6 AU, well
        # inside the planet). The far-field PPN formula is singular as
        # emag -> 0 and produced a spurious ~10" "deflection" of Mars seen
        # from Jupiter — physically, the deflection of light arriving at the
        # deflector's own centre vanishes. The threshold cannot affect any
        # supported geocentric/heliocentric observer (Earth is >= 0.98 AU
        # from the Sun and >= 4 AU from Jupiter/Saturn).
        if qmag == 0.0 or emag < 0.01:
            continue

        qhat = (pq[0] / qmag, pq[1] / qmag, pq[2] / qmag)
        ehat = (pe[0] / emag, pe[1] / emag, pe[2] / emag)

        pdotq = phat[0] * qhat[0] + phat[1] * qhat[1] + phat[2] * qhat[2]
        qdote = qhat[0] * ehat[0] + qhat[1] * ehat[1] + qhat[2] * ehat[2]
        edotp = ehat[0] * phat[0] + ehat[1] * phat[1] + ehat[2] * phat[2]

        # Skip if object is on the line-of-sight to the deflector
        if abs(edotp) > 0.99999999999:
            continue

        fac1 = 2.0 * _GS / (_C_MS * _C_MS * emag * _AU_M * rmass)
        fac2 = 1.0 + qdote
        if abs(fac2) < 1e-30:
            continue

        coeff = fac1 / fac2 * pmag
        result[0] += coeff * (pdotq * ehat[0] - edotp * qhat[0])
        result[1] += coeff * (pdotq * ehat[1] - edotp * qhat[1])
        result[2] += coeff * (pdotq * ehat[2] - edotp * qhat[2])

    return (result[0], result[1], result[2])


@lru_cache(maxsize=64)
def _get_skyfield_frame_data(
    jd_tt: float,
) -> Tuple[
    Tuple[Tuple[float, float, float], ...],
    float,
    float,
    float,
]:
    """Get precession-nutation matrix and nutation angles for the Skyfield path.

    The nutation angles still come from Skyfield's own IAU 2006/2000A model (so
    the modern nutation, and the LEB-vs-Skyfield comparison tests, are
    unchanged), but the precession and of-date mean obliquity now come from the
    Vondrák 2011 long-term model via ``vondrak_pn_matrix`` -- the same single
    source the LEB path uses, so the two backends cannot diverge at remote
    epochs. (Previously this read Skyfield's ``t.M`` IAU 2006 BPN matrix.)

    Args:
        jd_tt: Julian Day in TT.

    Returns:
        (pn_mat, dpsi, deps, eps_true_rad) where:
        - pn_mat: 3x3 ICRS→true-equatorial-of-date rotation as nested tuples
        - dpsi: nutation in longitude (radians, IAU 2006/2000A)
        - deps: nutation in obliquity (radians, IAU 2006/2000A)
        - eps_true_rad: true obliquity (radians)
    """
    from .cache import get_cached_time_tt

    t = get_cached_time_tt(jd_tt)

    # Nutation angles from Skyfield's own computation (unchanged).
    dpsi, deps = t._nutation_angles_radians
    dpsi = float(dpsi)
    deps = float(deps)

    # Vondrák 2011 long-term precession (via erfa) + Skyfield's nutation,
    # built into the ICRS->true-equator-of-date matrix.
    pn_mat, eps_true_rad = vondrak_pn_matrix(jd_tt, dpsi, deps)

    return pn_mat, dpsi, deps, eps_true_rad


@lru_cache(maxsize=64)
def _get_precession_matrix(
    jd_tt: float,
) -> Tuple[Tuple[float, float, float], ...]:
    """Get the ICRS->mean-equator-of-date rotation matrix (no nutation).

    Used for the public sidereal+equatorial convention, which is expressed on
    the mean equator.
    Now sourced from the Vondrák 2011 long-term precession (via erfa, ICRS frame
    bias included) -- the same single source as the LEB path -- so both backends
    agree and remain valid at remote epochs. (Previously this used Skyfield's
    ``mean_equator_and_equinox_of_date`` IAU 2006 frame.)

    The ``@lru_cache`` is retained so ``cache.clear_caches()`` can call
    ``_get_precession_matrix.cache_clear()`` (the underlying Vondrák matrix is
    state-independent, but the cache-clearing contract must hold).

    Args:
        jd_tt: Julian Day in TT.

    Returns:
        3x3 rotation matrix as nested tuples.
    """
    return vondrak_precession_matrix(jd_tt)


# =============================================================================
# PURE-LEB FRAME DATA (no Skyfield dependency)
# =============================================================================
# IAU 2006 Fukushima-Williams precession polynomials (IERS 2010 Table 5.1)
# combined with LEB-stored IAU 2006/2000A nutation to build the full
# bias-precession-nutation matrix without any Skyfield or erfa calls.

_AS2R = math.pi / (180.0 * 3600.0)  # arcseconds to radians
_J2000 = 2451545.0
_CENTURY = 36525.0


def _iau2006_precession_angles(
    jd_tt: float,
) -> Tuple[float, float, float, float]:
    """IAU 2006 Fukushima-Williams precession angles.

    Equivalent to erfa.pfw06. Returns (gamb, phib, psib, epsa) in radians.
    """
    T = (jd_tt - _J2000) / _CENTURY
    gamb = (
        -0.052928
        + (
            10.556378
            + (0.4932044 + (-0.00031238 + (-0.000002788 + 0.0000000260 * T) * T) * T)
            * T
        )
        * T
    ) * _AS2R
    phib = (
        84381.412819
        + (
            -46.811016
            + (0.0511268 + (0.00053289 + (-0.000000440 - 0.0000000176 * T) * T) * T) * T
        )
        * T
    ) * _AS2R
    psib = (
        -0.041775
        + (
            5038.481484
            + (1.5584175 + (-0.00018522 + (-0.000026452 - 0.0000000148 * T) * T) * T)
            * T
        )
        * T
    ) * _AS2R
    epsa = (
        84381.406
        + (
            -46.836769
            + (-0.0001831 + (0.00200340 + (-0.000000576 - 0.0000000434 * T) * T) * T)
            * T
        )
        * T
    ) * _AS2R
    return gamb, phib, psib, epsa


def _fw2m(
    gamb: float, phib: float, psi: float, eps: float
) -> Tuple[Tuple[float, float, float], ...]:
    """Fukushima-Williams angles to BPN rotation matrix.

    Computes R = R1(-eps) x R3(-psi) x R1(phib) x R3(gamb).
    Equivalent to erfa.fw2m / SOFA iauFw2m.
    """
    sg = math.sin(gamb)
    cg = math.cos(gamb)
    sp = math.sin(phib)
    cp = math.cos(phib)
    ss = math.sin(psi)
    cs = math.cos(psi)
    se = math.sin(eps)
    ce = math.cos(eps)

    # BA = R1(phib) x R3(gamb)
    -sg * cp
    cg * cp
    ba20 = sg * sp
    ba21 = -cg * sp

    # CBA = R3(-psi) x BA
    cba00 = cs * cg + ss * sg * cp
    cba01 = cs * sg - ss * cg * cp
    cba02 = -ss * sp
    cba10 = ss * cg - cs * sg * cp
    cba11 = ss * sg + cs * cg * cp
    cba12 = cs * sp

    # DCBA = R1(-eps) x CBA
    return (
        (cba00, cba01, cba02),
        (ce * cba10 - se * ba20, ce * cba11 - se * ba21, ce * cba12 - se * cp),
        (se * cba10 + ce * ba20, se * cba11 + ce * ba21, se * cba12 + ce * cp),
    )


# Thread-local cache for _get_leb_frame_data (keyed by exact jd_tt float).
# Using a plain dict + maxsize check is faster than @lru_cache for this
# because the reader argument makes lru_cache less efficient.
_leb_frame_cache: Dict[Tuple[int, float], Tuple] = {}
_LEB_FRAME_CACHE_MAX = 64


def _get_leb_frame_data(
    reader: "LEBReaderLike",
    jd_tt: float,
) -> Tuple[
    Tuple[Tuple[float, float, float], ...],
    float,
    float,
    float,
]:
    """Get precession-nutation matrix and nutation angles from LEB data.

    Skyfield-free replacement for _get_skyfield_frame_data(): the nutation comes
    from the LEB-stored Chebyshev coefficients, and the precession from the
    Vondrák 2011 long-term model via ``vondrak_pn_matrix`` (see
    ``precession_vondrak.py`` for the erfa/Vondrák provenance). Long-term-valid,
    unlike the IAU 2006 polynomials previously used here.

    Args:
        reader: LEBReader or LEB2Reader with nutation data.
        jd_tt: Julian Day in TT.

    Returns:
        Same (pn_mat, dpsi, deps, eps_true_rad) as _get_skyfield_frame_data.
    """
    # Keyed by (reader identity, jd): per-context readers may carry
    # different nutation tables, and a reader swap must not serve stale
    # frame data.
    _cache_key = (id(reader), jd_tt)
    cached = _leb_frame_cache.get(_cache_key)
    if cached is not None:
        return cached

    # Nutation from LEB Chebyshev (~1.5 µs)
    dpsi, deps = reader.eval_nutation(jd_tt)

    # Vondrák 2011 long-term precession (via erfa) combined with the LEB
    # nutation above into the ICRS->true-equator-of-date matrix.
    pn_mat, eps_true_rad = vondrak_pn_matrix(jd_tt, dpsi, deps)

    result = (pn_mat, dpsi, deps, eps_true_rad)

    # Cache with bounded size
    if len(_leb_frame_cache) >= _LEB_FRAME_CACHE_MAX:
        _leb_frame_cache.clear()
    _leb_frame_cache[_cache_key] = result

    return result


def _get_leb_precession_matrix(
    jd_tt: float,
) -> Tuple[Tuple[float, float, float], ...]:
    """Get the ICRS->mean-equator-of-date precession matrix (no nutation).

    Skyfield-free; uses the Vondrák 2011 long-term precession (via erfa,
    frame bias included) -- the same source as _get_precession_matrix().
    """
    return vondrak_precession_matrix(jd_tt)


# Thread-local reference to the active reader for frame data dispatch.
# Set by _fast_calc_core() and _pipeline_icrs() at the start of each
# calculation. Thread-local so two threads using different readers (e.g.
# per-context readers) never evaluate frame data from each other\'s
# reader. A global generation counter invalidates every thread\'s slot
# when state.close() releases the reader (via _reset_active_reader()).
_active_local = threading.local()
_active_generation: int = 0


def _reset_active_reader() -> None:
    """Invalidate all threads\' active-reader references.

    Called by state.close() to avoid stale references to a closed
    LEB reader (whose mmap has been released).
    """
    global _active_generation
    _active_generation += 1


def _set_active_reader(reader: "LEBReaderLike") -> None:
    """Bind the active reader for this thread\'s frame-data dispatch."""
    _active_local.gen = _active_generation
    _active_local.reader = reader
    _active_local.has_nutation = (
        hasattr(reader, "has_nutation") and reader.has_nutation()
    )


def _active_has_nutation() -> bool:
    """Return whether this thread's current LEB reader contains nutation data."""
    return getattr(_active_local, "gen", -1) == _active_generation and getattr(
        _active_local, "has_nutation", False
    )


def _frame_data(
    jd_tt: float,
) -> Tuple[Tuple[Tuple[float, float, float], ...], float, float, float]:
    """Get frame data from LEB (fast) or Skyfield (fallback).

    Uses the thread-local active reader set by _fast_calc_core().
    """
    if _active_has_nutation():
        return _get_leb_frame_data(_active_local.reader, jd_tt)
    return _get_skyfield_frame_data(jd_tt)


def _prec_matrix(
    jd_tt: float,
) -> Tuple[Tuple[float, float, float], ...]:
    """Get precession matrix from pure Python (fast) or Skyfield (fallback)."""
    if _active_has_nutation():
        return _get_leb_precession_matrix(jd_tt)
    return _get_precession_matrix(jd_tt)


def _mat3_vec3(
    mat: Tuple[Tuple[float, float, float], ...],
    vec: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """Multiply a 3x3 matrix by a 3-vector."""
    return (
        mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2],
        mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2],
        mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2],
    )


def _cotrans(lon: float, lat: float, eps: float) -> Tuple[float, float]:
    """Coordinate transform between ecliptic and equatorial.

    Negative eps: ecliptic -> equatorial.
    Positive eps: equatorial -> ecliptic.

    Args:
        lon: Longitude/RA in degrees.
        lat: Latitude/Dec in degrees.
        eps: Obliquity in degrees (sign determines direction).

    Returns:
        (transformed_lon, transformed_lat) in degrees.
    """
    eps_rad = math.radians(-eps)  # Ecliptic-to-equatorial rotation sign.
    lon_rad = math.radians(lon)
    lat_rad = math.radians(lat)

    cos_eps = math.cos(eps_rad)
    sin_eps = math.sin(eps_rad)
    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)

    # Spherical rotation
    x = cos_lat * cos_lon
    y = cos_lat * sin_lon * cos_eps - sin_lat * sin_eps
    z = cos_lat * sin_lon * sin_eps + sin_lat * cos_eps

    new_lon = math.degrees(math.atan2(y, x)) % 360.0
    r = math.sqrt(x * x + y * y + z * z)
    new_lat = math.degrees(math.asin(max(-1.0, min(1.0, z / r)))) if r > 0 else 0.0

    return new_lon, new_lat


def _precess_ecliptic(
    lon: float, lat: float, from_jd: float, to_jd: float
) -> Tuple[float, float]:
    """Precess ecliptic coordinates between two epochs.

    Uses the astrometry module's precession.

    Args:
        lon, lat: Ecliptic coordinates in degrees.
        from_jd, to_jd: Julian Days (TT).

    Returns:
        (precessed_lon, precessed_lat) in degrees.
    """
    from .astrometry import _precess_ecliptic as _pe

    return _pe(lon, lat, from_jd, to_jd)


# =============================================================================
# AYANAMSHA
# =============================================================================

# IAU 2006 general precession polynomial (arcsec/century)
_PREC_COEFFS = (5028.796195, 1.1054348, 0.00007964, -0.000023857, -0.0000000383)


def _general_precession_deg(jd_tt: float) -> float:
    """IAU 2006 general precession in longitude since J2000, in degrees."""
    centuries = (jd_tt - J2000) / 36525.0
    power = centuries
    arcseconds = 0.0
    for coefficient in _PREC_COEFFS:
        arcseconds += coefficient * power
        power *= centuries
    return arcseconds / 3600.0


def _general_precession_rate_deg_day(jd_tt: float) -> float:
    """General-precession rate dP/dT at ``jd_tt``, in degrees/day.

    The longitude-speed corrections (ayanamsa drift for sidereal SPEED output
    and the of-date→J2000 equinox-motion frame term) both subtract this rate
    from ``dlon``. ``_PREC_COEFFS`` are arcsec/century, so dP/dT =
    c0 + 2*c1*T + ... (only the first two terms matter at the required
    precision); convert to deg/day with / 3600 / 36525.

    Shared by both subtraction sites in ``_fast_calc_core`` and by the matching
    sites in ``horizons_backend._to_ecliptic_output`` so the LEB/Skyfield and
    Horizons backends stay identical.
    """
    T = (jd_tt - J2000) / 36525.0
    prec_rate_arcsec_cy = _PREC_COEFFS[0] + 2.0 * _PREC_COEFFS[1] * T
    return prec_rate_arcsec_cy / (3600.0 * 36525.0)


# Modes whose live catalog/frame geometry is evaluated by
# planets._calc_ayanamsa(). Formula modes remain local to the LEB path.
_DELEGATED_AYANAMSHA_MODES = DYNAMIC_AYANAMSHA_MODES

# Modes whose ayanamsha drift follows a live catalog direction instead of the
# IAU general-precession polynomial. Horizons imports this private set for the
# same derivative choice.
_STAR_BASED_MODES = DYNAMIC_AYANAMSHA_MODES


def _calc_ayanamsa_from_leb(
    reader: "LEBReaderLike",
    jd_tt: float,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
) -> float:
    """Compute ayanamsa for the LEB fast path.

    Epoch/formula modes and SIDM_USER are evaluated locally from the shared
    defining table. Live catalog and galactic modes delegate to
    ``planets._calc_ayanamsa``; the planet itself remains on the LEB path.

    Args:
        reader: Open LEBReader instance.
        jd_tt: Julian Day in TT.
        sid_mode: Sidereal mode ID (if None, reads from global state).
        sid_t0: Reference epoch JD for custom ayanamsha.
        sid_ayan_t0: Ayanamsha value at reference epoch (degrees).

    Returns:
        Ayanamsa in degrees.

    """
    if sid_mode is None:
        # Read via the locked getter so a concurrent EphemerisContext swap
        # (which installs its sidereal config into these globals under
        # _CONTEXT_SWAP_LOCK for the whole calc) cannot leak into a
        # module-level sidereal request. Reading the globals directly here
        # bypassed the lock and could return another context's ayanamsha.
        from .state import get_sid_mode

        sid_mode, sid_t0, sid_ayan_t0 = get_sid_mode(full=True)

    mode = sid_mode if sid_mode is not None else 0

    if mode in _DELEGATED_AYANAMSHA_MODES:
        # Late imports avoid a cycle. _calc_ayanamsa accepts UT, so invert
        # TT approximately with the same Delta T model used in the forward
        # conversion.
        from .planets import _calc_ayanamsa
        from .time_utils import deltat

        jd_ut = jd_tt - deltat(jd_tt)
        return _calc_ayanamsa(jd_ut, mode)

    if mode == 255:
        # SIDM_USER: sidereal zero point fixed on the mean ecliptic of t0.
        # ayanamsha(t) = ayan_t0 + accumP_B(t; t0), Method-B long-term
        # precession (Vondrák ecliptic frame). t0 is used directly as TT (no
        # J2000 sentinel). Mirrors planets._calc_ayanamsa's SIDM_USER branch.
        ayan_t0 = sid_ayan_t0 if sid_ayan_t0 is not None else 0.0
        t0 = sid_t0 if sid_t0 is not None else J2000
        mean_aya = ayan_t0 + method_b_accumulated_precession(jd_tt, t0)
    elif mode in AYANAMSHA_DEFINING:
        # Registered defining pair propagated with Method-B on the defining
        # epoch's mean ecliptic. The per-mode evidence status (primary,
        # geometric, secondary, or project convention) is intentionally kept
        # in docs/reference/ayanamsha.md; a numeric pair alone is not asserted
        # to be an author's published statement.
        defining_value, defining_epoch = AYANAMSHA_DEFINING[mode]
        mean_aya = defining_value + method_b_accumulated_precession(
            jd_tt, defining_epoch
        )
    else:
        raise KeyError(mode)

    # Mean ayanamsa (without nutation) - get_ayanamsa_ut() returns mean value
    return mean_aya % 360.0


# =============================================================================
# TOPOCENTRIC HELPERS
# =============================================================================

_EARTH_OMEGA_RAD_S = 7.2921150e-5  # Earth angular velocity (rad/s)
_SEC_PER_DAY = 86400.0


def _topocentric_offset(
    geopos: Tuple[float, float, float],
    jd_tt: float,
    jd_ut1: float,
    reader: "LEBReaderLike",
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Compute observer ICRS position and velocity offset from geocenter.

    Uses ERFA for WGS84→ITRF and celestial-to-terrestrial matrix.
    Returns (pos_au, vel_au_day) both as 3-tuples in ICRS.

    Args:
        geopos: (longitude_deg, latitude_deg, altitude_m)
        jd_tt: Julian Day TT
        jd_ut1: Julian Day UT1
        reader: LEB reader (for frame data)
    """
    import erfa

    lon_rad = math.radians(geopos[0])
    lat_rad = math.radians(geopos[1])
    alt_m = geopos[2]

    # WGS84 geodetic → ITRF geocentric (metres)
    xyz_itrf_m = erfa.gd2gc(1, lon_rad, lat_rad, alt_m)  # 1 = WGS84

    # Celestial-to-terrestrial matrix (polar motion xp=yp=0)
    c2t = erfa.c2t06a(
        2451545.0, jd_tt - 2451545.0, 2451545.0, jd_ut1 - 2451545.0, 0.0, 0.0
    )

    # Terrestrial → celestial (transpose)
    t2c_0 = (float(c2t[0][0]), float(c2t[1][0]), float(c2t[2][0]))
    t2c_1 = (float(c2t[0][1]), float(c2t[1][1]), float(c2t[2][1]))
    t2c_2 = (float(c2t[0][2]), float(c2t[1][2]), float(c2t[2][2]))

    # Position: ITRF → ICRS (AU)
    x_m = float(xyz_itrf_m[0])
    y_m = float(xyz_itrf_m[1])
    z_m = float(xyz_itrf_m[2])

    pos_au = (
        (t2c_0[0] * x_m + t2c_0[1] * y_m + t2c_0[2] * z_m) / _AU_M,
        (t2c_1[0] * x_m + t2c_1[1] * y_m + t2c_1[2] * z_m) / _AU_M,
        (t2c_2[0] * x_m + t2c_2[1] * y_m + t2c_2[2] * z_m) / _AU_M,
    )

    # Velocity: Earth rotation cross product in ITRF, then rotate to ICRS
    omega = _EARTH_OMEGA_RAD_S
    vx_itrf = -omega * y_m
    vy_itrf = omega * x_m
    vz_itrf = 0.0

    vel_au_day = (
        (t2c_0[0] * vx_itrf + t2c_0[1] * vy_itrf + t2c_0[2] * vz_itrf)
        * _SEC_PER_DAY
        / _AU_M,
        (t2c_1[0] * vx_itrf + t2c_1[1] * vy_itrf + t2c_1[2] * vz_itrf)
        * _SEC_PER_DAY
        / _AU_M,
        (t2c_2[0] * vx_itrf + t2c_2[1] * vy_itrf + t2c_2[2] * vz_itrf)
        * _SEC_PER_DAY
        / _AU_M,
    )

    return pos_au, vel_au_day


def _topo_ecliptic(
    reader: "LEBReaderLike",
    jd_tt: float,
    jd_ut1: float,
    ipl: int,
    geopos: Tuple[float, float, float],
    iflag: int = 0,
) -> Tuple[float, float, float, float, float, float]:
    """Compute topocentric ecliptic position without mutating global state.

    Runs _pipeline_icrs with observer = Earth + topocentric offset.
    Returns (lon, lat, dist, dlon, dlat, ddist) in ecliptic of date.

    Args:
        reader: Open LEB reader.
        jd_tt: Julian Day TT.
        jd_ut1: Julian Day UT1.
        ipl: Body ID (SUN, MOON, etc.)
        geopos: (longitude_deg, latitude_deg, altitude_m)
        iflag: Additional flags (FLG_SPEED, etc.). FLG_TOPOCTR is implied.
    """
    obs_pos, obs_vel = _topocentric_offset(geopos, jd_tt, jd_ut1, reader)

    want_velocity = bool(iflag & FLG_SPEED)

    # Check if body is a system barycenter (needs COB correction)
    _is_sys_bary = False
    if reader.has_body(ipl):
        _body = reader._bodies[ipl]
        _is_sys_bary = _body.coord_type == COORD_ICRS_BARY_SYSTEM

    result = _pipeline_icrs(
        reader,
        jd_tt,
        ipl,
        iflag & ~FLG_TOPOCTR,
        want_velocity=want_velocity,
        is_system_bary=_is_sys_bary,
        topo_offset=(obs_pos, obs_vel),
    )

    if want_velocity:
        return result  # type: ignore[return-value]
    return result[0], result[1], result[2], 0.0, 0.0, 0.0  # type: ignore[return-value]


def _apparent_icrs_cartesian(
    reader: "LEBReaderLike",
    jd_tt: float,
    ipl: int,
) -> Tuple[float, float, float]:
    """Geocentric apparent ICRS Cartesian position (post-aberration/deflection).

    Runs _pipeline_icrs steps 1-6 and returns the geo[] vector in ICRS
    without any frame rotation. Equivalent to Skyfield's .position.au
    on an apparent() object.

    Used by eclipse.py Pattern E for Besselian shadow geometry.
    """
    result = _pipeline_icrs(
        reader,
        jd_tt,
        ipl,
        FLG_EQUATORIAL | FLG_J2000,
        want_velocity=False,
        want_xyz=True,
    )
    return result[0], result[1], result[2]  # type: ignore[return-value]


# =============================================================================
# PIPELINE A: ICRS BARYCENTRIC BODIES
# =============================================================================

# Central-difference half-step (days) for the Pipeline-A velocity. The velocity
# is built analytically from the exact Chebyshev derivative for the body's
# geometric motion (no float64-jd quantization), then the aberration,
# deflection and frame-rotation contributions are recovered by differencing
# those smooth corrections over ±_VEL_H on analytically-extrapolated inputs.
# Because the body's geometric vector is extrapolated analytically (not by
# re-evaluating the body ephemeris at jd±h), _VEL_H does NOT suffer the ULP
# quantization that a naive position central difference would. It can therefore
# be kept small to bound the O(h**2) truncation of the fast-moving Moon's terms
# — dominated by the topocentric diurnal parallax rate — while leaving the
# planets exact. The lower bound is set by the smooth position-level inputs that
# ARE re-evaluated (frame matrices, COB and topocentric offsets): their ULP
# noise is amplified by 1/(2h). The 0.0005-day half-step keeps both truncation
# and round-off terms small; it is a numerical differentiation parameter, not
# an astronomical fitted constant.
_VEL_H = 0.0005


def _eval_body_center_state(
    reader: "LEBReaderLike",
    ipl: int,
    jd_tt: float,
    is_system_bary: bool,
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """Evaluate the requested body center and its barycentric state.

    LEB records for outer planets can store a system barycentre. When a JPL
    planet-center segment covers the epoch, add that segment's position and a
    centered derivative of the same offset. Otherwise retain the system
    barycentre as an explicit data-availability fallback.
    """
    pos, vel = reader.eval_body(ipl, jd_tt)
    if not is_system_bary:
        return pos, vel

    center_pos = _apply_cob_correction(pos, ipl, jd_tt)
    if center_pos == pos:
        return pos, vel

    zero = (0.0, 0.0, 0.0)
    offset_before = _apply_cob_correction(zero, ipl, jd_tt - _VEL_H)
    offset_after = _apply_cob_correction(zero, ipl, jd_tt + _VEL_H)
    scale = 1.0 / (2.0 * _VEL_H)
    center_vel = (
        vel[0] + (offset_after[0] - offset_before[0]) * scale,
        vel[1] + (offset_after[1] - offset_before[1]) * scale,
        vel[2] + (offset_after[2] - offset_before[2]) * scale,
    )
    return center_pos, center_vel


def _nutation_rate_deg_day(jd_tt: float) -> float:
    """Central-difference rate of nutation in longitude dΔψ/dt (deg/day).

    Reads the nutation from ``_frame_data`` (the active LEB reader's stored
    IAU 2006/2000A series, or Skyfield's own model as fallback) so the rate is
    consistent with the Δψ this backend adds to positions. Shares the ``_VEL_H``
    step with the sidereal ayanamsha-rate correction in ``_fast_calc_core`` so
    the two dΔψ/dt terms cancel to the intrinsic rate exactly.
    """
    _, dpsi_m, _, _ = _frame_data(jd_tt - _VEL_H)
    _, dpsi_p, _, _ = _frame_data(jd_tt + _VEL_H)
    return math.degrees(dpsi_p - dpsi_m) / (2.0 * _VEL_H)


def _frame_transform(
    geo: Tuple[float, float, float],
    jd_tt: float,
    iflag: int,
    want_xyz: bool,
) -> Tuple[float, float, float]:
    """Rotate an apparent ICRS vector into the requested output frame.

    Shared by the Pipeline-A position and velocity paths so that the reported
    speed is the exact derivative of the reported position (the frame rotation
    is evaluated at the sample epoch, capturing the precession/nutation
    equinox motion). Returns spherical (lon_deg, lat_deg, dist_au) unless
    *want_xyz*, in which case the Cartesian components of the rotated vector
    are returned.
    """
    if (iflag & FLG_EQUATORIAL) and (iflag & FLG_J2000):
        # Mean equator/equinox of J2000: ICRS rotated by IAU 2006 frame bias.
        v = _rotate_icrs_to_j2000_mean_equator(geo[0], geo[1], geo[2])
    elif (iflag & FLG_EQUATORIAL) and (iflag & FLG_SIDEREAL):
        # SID+EQ: mean equator of date (P matrix, no nutation).
        v = _mat3_vec3(_prec_matrix(jd_tt), geo)
    elif iflag & FLG_EQUATORIAL:
        # NONUT: mean equator (P matrix); else true equator (PNM).
        if iflag & FLG_NONUT:
            _eq_mat = _prec_matrix(jd_tt)
        else:
            _eq_mat, _, _, _ = _frame_data(jd_tt)
        v = _mat3_vec3(_eq_mat, geo)
    elif iflag & FLG_J2000:
        v = _rotate_icrs_to_ecliptic_j2000(geo[0], geo[1], geo[2])
    else:
        # TRUE ECLIPTIC OF DATE (default). FLG_NONUT: mean ecliptic (no nutation).
        if iflag & FLG_NONUT:
            _rot_mat = _prec_matrix(jd_tt)
            eps_rad = math.radians(vondrak_mean_obliquity_deg(jd_tt))
        else:
            _rot_mat, _, _, eps_rad = _frame_data(jd_tt)
        geo_eq = _mat3_vec3(_rot_mat, geo)
        v = _rotate_equatorial_to_ecliptic(geo_eq[0], geo_eq[1], geo_eq[2], eps_rad)

    if want_xyz:
        return v[0], v[1], v[2]
    return _cartesian_to_spherical(v[0], v[1], v[2])


def _pipeline_icrs(
    reader: "LEBReaderLike",
    jd_tt: float,
    ipl: int,
    iflag: int,
    want_velocity: bool = False,
    is_system_bary: bool = False,
    topo_offset: Optional[
        Tuple[Tuple[float, float, float], Tuple[float, float, float]]
    ] = None,
    want_xyz: bool = False,
    topo_offset_fn: Optional[
        Callable[[float], Tuple[Tuple[float, float, float], Tuple[float, float, float]]]
    ] = None,
) -> Tuple[float, ...]:
    """Pipeline A: compute ecliptic coordinates for ICRS barycentric bodies.

    Handles both planet-center bodies (COORD_ICRS_BARY) and system-barycenter
    bodies (COORD_ICRS_BARY_SYSTEM). For system barycenters, a JPL
    center-of-body segment is applied at runtime when available.

    When want_velocity is False (default), returns (lon, lat, dist).
    When want_velocity is True, returns (lon, lat, dist, dlon, dlat, ddist)
    where the velocity components are analytically derived from the Chebyshev
    polynomial derivatives, transformed through the same rotation matrices
    as the position.

    Args:
        reader: Open LEBReader.
        jd_tt: Julian Day TT.
        ipl: Body ID.
        iflag: Flags.
        want_velocity: Whether to compute velocity.
        is_system_bary: If True, stored data is system barycenter; apply COB.

    Returns:
        (lon_deg, lat_deg, dist_au) or
        (lon_deg, lat_deg, dist_au, dlon_deg_day, dlat_deg_day, ddist_au_day)
    """
    # Bind the active reader for _frame_data / _prec_matrix dispatch.
    # Callers like _topo_ecliptic enter this pipeline without going through
    # _fast_calc_core, so the module-level reference must be refreshed here
    # to avoid using a stale (closed) reader after state.close().
    _set_active_reader(reader)

    # 1. Evaluate the public body center, retaining the system barycentre only
    # when no JPL center segment covers this epoch.
    target_pos, target_vel = _eval_body_center_state(reader, ipl, jd_tt, is_system_bary)

    # Pre-initialize velocity variables to satisfy type checker.
    # These are always set before use when want_velocity=True.
    geo_vel: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    dlon = dlat = ddist = 0.0

    # 2. Observer selection (defer Earth fetch for helio/bary)
    if iflag & FLG_HELCTR:
        sun_pos, sun_vel = reader.eval_body(SUN, jd_tt)
        observer = sun_pos
        observer_vel = sun_vel
        earth_vel = (0.0, 0.0, 0.0)  # Not needed (aberration skipped)
    elif iflag & FLG_BARYCTR:
        observer = (0.0, 0.0, 0.0)
        observer_vel = (0.0, 0.0, 0.0)
        earth_vel = (0.0, 0.0, 0.0)  # Not needed (aberration skipped)
    else:
        # Geocentric (default) — need Earth position and velocity
        earth_pos, earth_vel = reader.eval_body(EARTH, jd_tt)
        if topo_offset is not None:
            topo_pos, topo_vel = topo_offset
            observer = (
                earth_pos[0] + topo_pos[0],
                earth_pos[1] + topo_pos[1],
                earth_pos[2] + topo_pos[2],
            )
            observer_vel = (
                earth_vel[0] + topo_vel[0],
                earth_vel[1] + topo_vel[1],
                earth_vel[2] + topo_vel[2],
            )
            earth_vel = observer_vel
        else:
            observer = earth_pos
            observer_vel = earth_vel

    # 3. Geometric vector
    geo = _vec3_sub(target_pos, observer)

    # 4. Light-time correction (unless FLG_TRUEPOS). The same public body
    # center is evaluated at every retarded epoch.
    retarded_vel = target_vel
    lt = 0.0
    if not (iflag & FLG_TRUEPOS):
        # Solve |target(t-lt) - observer(t)| = c*lt by fixed-point iteration.
        # For HELCTR the observer is the Sun at the observation epoch; for
        # BARYCTR it is the SSB.
        retarded_pos = target_pos
        for _ in range(3):  # Fixed-point iterations
            dist = _vec3_dist(geo)
            if dist == 0.0:
                break
            lt = dist / C_LIGHT_AU_DAY
            retarded_pos, retarded_vel = _eval_body_center_state(
                reader, ipl, jd_tt - lt, is_system_bary
            )
            geo = _vec3_sub(retarded_pos, observer)

    if want_velocity:
        # Exact geometric velocity of the light-time-corrected vector, from the
        # Chebyshev derivative (no float64-jd quantization). The aberration,
        # deflection and frame-rotation contributions are added below.
        geo_vel = _vec3_sub(retarded_vel, observer_vel)
        geo_geom = geo  # pre-deflection, pre-aberration apparent vector at jd_tt

        # Light-time retardation rate: the apparent place is target(t - lt(t)),
        # whose exact derivative carries a factor (1 - d(lt)/dt) on the target's
        # retarded velocity. Differentiating the propagation equation gives
        #   lt' = p̂·(r' - o') / (c + p̂·r')
        # from differentiating |target(t-lt) - observer(t)| = c·lt, where p̂ is
        # the unit vector along the observer-to-target propagation distance.
        if not (iflag & FLG_TRUEPOS) and lt > 0.0:
            _dg = _vec3_dist(geo_geom)
            if _dg > 0.0:
                _px, _py, _pz = (
                    geo_geom[0] / _dg,
                    geo_geom[1] / _dg,
                    geo_geom[2] / _dg,
                )
                _pr = (
                    _px * retarded_vel[0]
                    + _py * retarded_vel[1]
                    + _pz * retarded_vel[2]
                )
                _pv = _px * geo_vel[0] + _py * geo_vel[1] + _pz * geo_vel[2]
                _ltdot = _pv / (C_LIGHT_AU_DAY + _pr)
                geo_vel = (
                    retarded_vel[0] * (1.0 - _ltdot) - observer_vel[0],
                    retarded_vel[1] * (1.0 - _ltdot) - observer_vel[1],
                    retarded_vel[2] * (1.0 - _ltdot) - observer_vel[2],
                )

    # 5. Gravitational deflection by Sun, Jupiter, Saturn (finite-source PPN
    #    point-mass reduction documented in the IERS/Skyfield source chain).
    #    The correction can be arcsecond-scale near solar conjunction; its
    #    exact magnitude depends on finite observer/deflector/target geometry.
    #    Skipped for helio/bary/truepos/nogdefl and for the Moon (negligible at
    #    ~0.0026 AU, deflection < 0.000001"). NOABERR deliberately does NOT
    #    skip deflection: the public flag contract assigns that bit only to
    #    aberration (FLG_ASTROMETRIC combines NOABERR and NOGDEFL).
    _do_defl = not (iflag & (FLG_NOGDEFL | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS))
    _do_defl = _do_defl and ipl != MOON and lt > 0.0
    _defl_delta: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    if _do_defl:
        _pre_defl = geo
        geo = _apply_gravitational_deflection(geo, observer, jd_tt, lt, reader)
        if want_velocity:
            _defl_delta = _vec3_sub(geo, _pre_defl)

    # 6. Aberration (full special-relativistic, matching Skyfield).
    _do_aber = not (iflag & (FLG_NOABERR | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS))
    if _do_aber:
        geo = _apply_aberration(geo, earth_vel, lt)

    # 6. Coordinate transform (position).
    #
    # Sidereal handling for Pipeline A (ICRS bodies):
    #   SID+EQ: use mean equator (P matrix, no nutation) — NO ayanamsha subtraction
    #   SID+EQ+J2K: same as non-sidereal EQ+J2K (ICRS ≡ J2000 equatorial)
    #   SID only: output ecliptic-of-date, ayanamsha subtracted in _fast_calc_core
    #   SID+J2K: output J2000 ecliptic, ayanamsha subtracted in _fast_calc_core
    _want_xyz = want_xyz or bool(iflag & FLG_XYZ)

    lon_deg, lat_deg, dist = _frame_transform(geo, jd_tt, iflag, _want_xyz)

    if not want_velocity:
        return lon_deg, lat_deg, dist

    # 7. Velocity: the reported speed must be the exact time-derivative of the
    #    reported position. The analytic geometric velocity (geo_vel) omits the
    #    time-derivatives of the aberration/deflection corrections (up to
    #    ~4"/day for the Moon) and of the of-date frame rotation
    #    (precession + nutation of the equinox, ~0.09"/day). Recover the full
    #    derivative by central-differencing the apparent place through the same
    #    frame transform at jd_tt ± _VEL_H, extrapolating the body's geometric
    #    vector analytically (so no ephemeris re-evaluation at jd±h and hence no
    #    ULP quantization) while re-forming the smooth aberration/deflection/
    #    frame terms at each sample epoch. This mirrors the Skyfield backend,
    #    which central-differences its full apparent-place pipeline.
    #
    # Topocentric observer at the sample epochs. Linear extrapolation of the
    # observer (via geo_vel, which carries observer_vel) misses the diurnal
    # Earth-rotation curvature of the topocentre; recompute the true offset at
    # jd_tt ± _VEL_H so the parallax and diurnal-aberration rates are exact.
    _tofs_m: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = (
        None
    )
    _tofs_p: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = (
        None
    )
    if topo_offset is not None and topo_offset_fn is not None:
        _tofs_m = topo_offset_fn(jd_tt - _VEL_H)
        _tofs_p = topo_offset_fn(jd_tt + _VEL_H)

    ev_m: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    ev_p: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    if _do_aber:
        # Observer (Earth [+ topocentre]) velocity at the sample epochs drives
        # the annual (and diurnal, for topocentric) aberration rate.
        _, _evm = reader.eval_body(EARTH, jd_tt - _VEL_H)
        _, _evp = reader.eval_body(EARTH, jd_tt + _VEL_H)
        if topo_offset is not None:
            _tvm = _tofs_m[1] if _tofs_m is not None else topo_offset[1]
            _tvp = _tofs_p[1] if _tofs_p is not None else topo_offset[1]
            ev_m = (_evm[0] + _tvm[0], _evm[1] + _tvm[1], _evm[2] + _tvm[2])
            ev_p = (_evp[0] + _tvp[0], _evp[1] + _tvp[1], _evp[2] + _tvp[2])
        else:
            ev_m, ev_p = _evm, _evp

    def _apparent_at(
        s: float,
        ev_s: Tuple[float, float, float],
        tofs_s: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]],
    ) -> Tuple[float, float, float]:
        g = (
            geo_geom[0] + geo_vel[0] * s,
            geo_geom[1] + geo_vel[1] * s,
            geo_geom[2] + geo_vel[2] * s,
        )
        if tofs_s is not None and topo_offset is not None:
            # Replace the linearly-extrapolated topocentre with the true one:
            # subtract (topo(t+s) - topo(t) - topo_vel·s), the diurnal curvature.
            _tp = tofs_s[0]
            _t0 = topo_offset[0]
            _tv = topo_offset[1]
            g = (
                g[0] - (_tp[0] - _t0[0] - _tv[0] * s),
                g[1] - (_tp[1] - _t0[1] - _tv[1] * s),
                g[2] - (_tp[2] - _t0[2] - _tv[2] * s),
            )
        # Light-time to the target at the sample epoch (the extrapolated
        # geometric vector has magnitude c·lt_s). The aberration formula
        # normalises by this distance, so a frozen lt would bias it.
        lt_s = _vec3_dist(g) / C_LIGHT_AU_DAY
        if _do_defl:
            if ipl == SUN:
                # Sun target: the Sun is its own dominant "deflector", so the
                # PPN geometry is degenerate (q ~ the Sun's barycentric motion
                # over lt, a ~1e-7 AU vector). Re-evaluating it on the
                # extrapolated vector produces a huge spurious rate (~8"/day);
                # physically the term is essentially constant over ±_VEL_H, so
                # freeze the position-path delta. (The degenerate case is
                # handled here on its own — the non-degenerate deflectors below
                # are NOT frozen.)
                g = (
                    g[0] + _defl_delta[0],
                    g[1] + _defl_delta[1],
                    g[2] + _defl_delta[2],
                )
            else:
                # Deflectors AND observer evaluated at the sample epoch. Near
                # solar conjunction the impact-parameter rate is dominated by
                # the Sun's apparent motion (~0.96 deg/day), not by the
                # target's: freezing the deflectors at jd_tt biased slow
                # targets by ~0.2-0.3"/day at elongation < 2 deg (and ~20"/day
                # through a limb transit). The observer is extrapolated
                # linearly — ample for the smooth pe = observer - deflector
                # geometry over ±_VEL_H.
                obs_s = (
                    observer[0] + observer_vel[0] * s,
                    observer[1] + observer_vel[1] * s,
                    observer[2] + observer_vel[2] * s,
                )
                g = _apply_gravitational_deflection(g, obs_s, jd_tt + s, lt_s, reader)
        if _do_aber:
            g = _apply_aberration(g, ev_s, lt_s)
        return g

    w_m = _frame_transform(
        _apparent_at(-_VEL_H, ev_m, _tofs_m), jd_tt - _VEL_H, iflag, _want_xyz
    )
    w_p = _frame_transform(
        _apparent_at(_VEL_H, ev_p, _tofs_p), jd_tt + _VEL_H, iflag, _want_xyz
    )
    inv2h = 1.0 / (2.0 * _VEL_H)
    if _want_xyz:
        dlon = (w_p[0] - w_m[0]) * inv2h
        dlat = (w_p[1] - w_m[1]) * inv2h
        ddist = (w_p[2] - w_m[2]) * inv2h
    else:
        _dl = (w_p[0] - w_m[0]) % 360.0
        if _dl > 180.0:
            _dl -= 360.0
        dlon = _dl * inv2h
        dlat = (w_p[1] - w_m[1]) * inv2h
        ddist = (w_p[2] - w_m[2]) * inv2h

    return lon_deg, lat_deg, dist, dlon, dlat, ddist


# =============================================================================
# PIPELINE B: ECLIPTIC DIRECT BODIES
# =============================================================================


# Geocentric gravitational parameter of the Earth-Moon system in
# AU^3/day^2: GM_sun (Gaussian constant squared) over the Sun/(E+M)
# mass ratio.
_GM_EARTH_MOON = 0.01720209895**2 / 328900.5596


def _osculating_node_radius(reader: "LEBReaderLike", jd_tt: float) -> float:
    """Osculating-ellipse radius at the Moon's ascending node.

    Reconstructs the geocentric lunar state vector from the LEB Moon and
    Earth channels (barycentric ICRS Chebyshev states), rotates it to the
    ecliptic of date, computes the instantaneous two-body elements, and
    evaluates ``r = p / (1 + e*cos(nu_node))`` from that independently sourced
    state.
    """
    moon_pos, moon_vel = reader.eval_body(MOON, jd_tt)
    earth_pos, earth_vel = reader.eval_body(EARTH, jd_tt)
    pn_mat, _, _, eps_rad = _frame_data(jd_tt)
    ce, se = math.cos(eps_rad), math.sin(eps_rad)

    def _to_ecl(vec):
        x, y, z = _mat3_vec3(pn_mat, vec)
        return (x, y * ce + z * se, -y * se + z * ce)

    r = _to_ecl(tuple(moon_pos[i] - earth_pos[i] for i in range(3)))
    v = _to_ecl(tuple(moon_vel[i] - earth_vel[i] for i in range(3)))
    hx = r[1] * v[2] - r[2] * v[1]
    hy = r[2] * v[0] - r[0] * v[2]
    hz = r[0] * v[1] - r[1] * v[0]
    h2 = hx * hx + hy * hy + hz * hz
    p = h2 / _GM_EARTH_MOON
    r_mag = math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)
    v2 = v[0] ** 2 + v[1] ** 2 + v[2] ** 2
    r_dot_v = r[0] * v[0] + r[1] * v[1] + r[2] * v[2]
    c1 = v2 / _GM_EARTH_MOON - 1.0 / r_mag
    c2 = r_dot_v / _GM_EARTH_MOON
    ex = c1 * r[0] - c2 * v[0]
    ey = c1 * r[1] - c2 * v[1]
    ez = c1 * r[2] - c2 * v[2]
    e_mag = math.sqrt(ex * ex + ey * ey + ez * ez)
    # Ascending node direction n = k x h; cos(nu_node) = n_hat . e_hat
    nx, ny = -hy, hx
    n_mag = math.hypot(nx, ny)
    if n_mag <= 0.0 or e_mag <= 0.0:
        raise ValueError("degenerate lunar state")
    cos_nu = (nx * ex + ny * ey) / (n_mag * e_mag)
    return p / (1.0 + e_mag * max(-1.0, min(1.0, cos_nu)))


def _pipeline_ecliptic(
    reader: "LEBReaderLike",
    jd_tt: float,
    ipl: int,
    iflag: int,
) -> Tuple[float, float, float, float, float, float]:
    """Pipeline B: evaluate ecliptic-direct bodies.

    Returns:
        (lon, lat, dist, dlon, dlat, ddist)
    """
    # Public frame convention: SIDEREAL+EQUATORIAL and J2000/NONUT use the
    # mean (nutation-free) ecliptic; other of-date requests use the true one.
    _sid_eq = bool(iflag & FLG_SIDEREAL) and bool(iflag & FLG_EQUATORIAL)
    _want_nonut_b = bool(iflag & (FLG_NONUT | FLG_J2000)) or _sid_eq

    _clean_model_body = ipl in (MEAN_NODE, MEAN_APOG, INTP_APOG, INTP_PERG)

    def _nutation_rate(jd: float) -> float:
        """Return the local longitude-nutation rate in degrees per day."""
        step = 0.001
        _, dpsi_minus, _, _ = _frame_data(jd - step)
        _, dpsi_plus, _, _ = _frame_data(jd + step)
        return math.degrees(dpsi_plus - dpsi_minus) / (2.0 * step)

    def _ofdate_state(
        jd: float,
    ) -> Tuple[float, float, float, float, float, float]:
        """Return a reported of-date state for an ecliptic-direct body.

        Mean and interpolated apsides always come from the packaged clean-room
        models, never from their legacy LEB channels.  True-node and
        osculating-apogee positions retain the JPL-derived LEB path.
        """
        use_speed3_state = bool(iflag & FLG_SPEED3)
        use_speed3_velocity = bool(iflag & FLG_SPEED) and (
            use_speed3_state or bool(iflag & FLG_TOPOCTR)
        )
        if ipl == MEAN_NODE:
            from .lunar import calc_mean_lunar_node_state

            lo, la, di, dlo, dla, ddi = calc_mean_lunar_node_state(
                jd, speed3=use_speed3_velocity
            )
        elif ipl == MEAN_APOG:
            from .lunar import calc_mean_lilith_state

            lo, la, di, dlo, dla, ddi = calc_mean_lilith_state(
                jd, speed3=use_speed3_velocity
            )
        elif ipl in (INTP_APOG, INTP_PERG):
            from .lunar import calc_interpolated_apse_state

            lo, la, di, dlo, dla, ddi = calc_interpolated_apse_state(
                jd, ipl, speed3=use_speed3_state
            )
            if use_speed3_velocity and not use_speed3_state:
                speed3 = calc_interpolated_apse_state(jd, ipl, speed3=True)
                dlo, dla, ddi = speed3[3:]
        else:
            (lo, la, di), _ = reader.eval_body(ipl, jd)
            dlo = dla = ddi = 0.0

        if ipl == TRUE_NODE:
            # The direct channel stores an angular-momentum proxy rather than
            # an osculating node radius. Reconstruct the geocentric lunar state
            # and evaluate the two-body ellipse at the ascending node.
            try:
                di = _osculating_node_radius(reader, jd)
            except (KeyError, ValueError, ArithmeticError):
                pass  # keep the stored proxy if the Moon channel is missing

        # Nutation in longitude (Δψ) handling. Mean bodies are stored on the
        # mean ecliptic (add Δψ for the true ecliptic of date); true/osculating
        # bodies include it (strip Δψ for a nutation-free frame). FLG_NONUT,
        # FLG_J2000 and SIDEREAL+EQUATORIAL select the nutation-free frame.
        if ipl in (MEAN_NODE, MEAN_APOG):
            if not _want_nonut_b:
                _, dpsi_rad, _, _ = _frame_data(jd)
                lo = (lo + math.degrees(dpsi_rad)) % 360.0
                dlo += _nutation_rate(jd)
        elif ipl in (TRUE_NODE, OSCU_APOG, INTP_APOG, INTP_PERG):
            if _want_nonut_b:
                _, dpsi_rad, _, _ = _frame_data(jd)
                lo = (lo - math.degrees(dpsi_rad)) % 360.0
                if _clean_model_body:
                    dlo -= _nutation_rate(jd)

        return lo, la, di, dlo, dla, ddi

    lon, lat, dist, dlon, dlat, ddist = _ofdate_state(jd_tt)

    # Velocity: the reported speed is the central difference of the reported
    # of-date position, over the SAME window the Skyfield backend uses, so the
    # two backends agree per body:
    #   * OscuApog's osculating curve oscillates sub-daily -> tight 0.05 d window
    #   * every other ecliptic-direct point -> 0.5 d
    # Because _ofdate_pos already carries the distance/latitude overrides and
    # Δψ, dlon/dlat/ddist are automatically the true derivatives of the
    # reported position — including the nutation rate (matching the ±Δψ on the
    # position) and the osculating-node/3-harmonic-latitude motion (fixing the
    # previously stale ddist/dlat that were left as the retired channel's
    # Chebyshev derivative). At an ephemeris edge the ±dt samples may fall out
    # of range; return zero speed there, matching the Skyfield boundary path.
    if (iflag & FLG_SPEED) and not _clean_model_body:
        _dt = (
            0.0001
            if iflag & (FLG_SPEED3 | FLG_TOPOCTR)
            else (0.05 if ipl in (OSCU_APOG, TRUE_NODE) else 0.5)
        )
        try:
            lo_m, la_m, di_m, _, _, _ = _ofdate_state(jd_tt - _dt)
            lo_p, la_p, di_p, _, _, _ = _ofdate_state(jd_tt + _dt)
            _dl = lo_p - lo_m
            if _dl > 180.0:
                _dl -= 360.0
            elif _dl < -180.0:
                _dl += 360.0
            dlon = _dl / (2.0 * _dt)
            dlat = (la_p - la_m) / (2.0 * _dt)
            ddist = (di_p - di_m) / (2.0 * _dt)
        except (KeyError, ValueError, ArithmeticError):
            dlon = dlat = ddist = 0.0
    elif not (iflag & FLG_SPEED):
        dlon = dlat = ddist = 0.0

    # Coordinate transforms for ecliptic-direct bodies.
    # Input coords are always ecliptic of date.
    #
    # FLG_J2000 is honored uniformly for every ecliptic-direct body,
    # including the non-mean lunar points under SIDEREAL.
    _effective_j2000 = bool(iflag & FLG_J2000)

    if (iflag & FLG_EQUATORIAL) and _effective_j2000:
        # J2000 equatorial: precess ecliptic-of-date -> J2000 ecliptic,
        # then rotate J2000 ecliptic -> J2000 equatorial.
        eps = OBLIQUITY_J2000_REF_DEG

        def _ecl_date_to_eq_j2000(
            lo: float, la: float, src_epoch: float
        ) -> tuple[float, float]:
            # The source of-date epoch varies between the two velocity samples
            # so the precession captures the ~50.29"/yr equinox drift; the
            # reported J2000 speed is then the true derivative of the reported
            # J2000 position, consistently across the two backends.
            lo_j, la_j = _precess_ecliptic(lo, la, src_epoch, J2000)
            return _cotrans(lo_j, la_j, -eps)

        # Velocity via finite difference on original ecliptic coords
        dt_step = 0.001  # days
        eq_now_lon, eq_now_lat = _ecl_date_to_eq_j2000(lon, lat, jd_tt)
        eq_fwd_lon, eq_fwd_lat = _ecl_date_to_eq_j2000(
            lon + dlon * dt_step, lat + dlat * dt_step, jd_tt + dt_step
        )
        d_eq_lon = eq_fwd_lon - eq_now_lon
        if d_eq_lon > 180.0:
            d_eq_lon -= 360.0
        elif d_eq_lon < -180.0:
            d_eq_lon += 360.0
        dlon = d_eq_lon / dt_step
        dlat = (eq_fwd_lat - eq_now_lat) / dt_step
        lon = eq_now_lon
        lat = eq_now_lat

    elif iflag & FLG_EQUATORIAL:
        # Equatorial of date: rotate ecliptic-of-date → equatorial-of-date.
        # A fixed-obliquity rotation (no equinox drift to capture here).
        # Sidereal/NONUT modes use mean obliquity (no nutation).
        if (iflag & FLG_SIDEREAL) or (iflag & FLG_NONUT):
            eps = vondrak_mean_obliquity_deg(jd_tt)
        else:
            _, _, deps, _ = _frame_data(jd_tt)
            eps_mean = vondrak_mean_obliquity_deg(jd_tt)
            eps = eps_mean + math.degrees(deps)

        # Velocity via finite difference on original ecliptic coords
        dt_step = 0.001  # days
        eq_now_lon, eq_now_lat = _cotrans(lon, lat, -eps)
        eq_fwd_lon, eq_fwd_lat = _cotrans(
            lon + dlon * dt_step, lat + dlat * dt_step, -eps
        )
        d_eq_lon = eq_fwd_lon - eq_now_lon
        if d_eq_lon > 180.0:
            d_eq_lon -= 360.0
        elif d_eq_lon < -180.0:
            d_eq_lon += 360.0
        dlon = d_eq_lon / dt_step
        dlat = (eq_fwd_lat - eq_now_lat) / dt_step
        lon = eq_now_lon
        lat = eq_now_lat

    elif _effective_j2000:
        # J2000 ecliptic: precess from ecliptic of date to J2000. The source
        # of-date epoch varies between the two velocity samples so the reported
        # J2000 speed is the true derivative of the reported J2000 position —
        # it differs from the of-date speed by the general-precession rate
        # (~0.14"/day), as required by the time-dependent frame rotation.
        dt_step = 0.001  # days
        j_now_lon, j_now_lat = _precess_ecliptic(lon, lat, jd_tt, J2000)
        j_fwd_lon, j_fwd_lat = _precess_ecliptic(
            lon + dlon * dt_step, lat + dlat * dt_step, jd_tt + dt_step, J2000
        )
        d_j_lon = j_fwd_lon - j_now_lon
        if d_j_lon > 180.0:
            d_j_lon -= 360.0
        elif d_j_lon < -180.0:
            d_j_lon += 360.0
        dlon = d_j_lon / dt_step
        dlat = (j_fwd_lat - j_now_lat) / dt_step
        lon = j_now_lon
        lat = j_now_lat

    return lon, lat, dist, dlon, dlat, ddist


# =============================================================================
# PIPELINE C: HELIOCENTRIC BODIES
# =============================================================================


def _pipeline_helio(
    reader: "LEBReaderLike",
    jd_tt: float,
    ipl: int,
    iflag: int,
) -> Tuple[float, float, float, float, float, float]:
    """Pipeline C: evaluate reviewed heliocentric ecliptic bodies.

    LEB stores these as heliocentric J2000 ecliptic (lon, lat, dist).
    Default output is geocentric ecliptic of date, consistent with the
    Skyfield path in planets.py.

    Returns:
        (lon, lat, dist, dlon, dlat, ddist)
    """
    is_helio = bool(iflag & FLG_HELCTR)

    if is_helio:
        # Heliocentric: LEB data is already heliocentric J2000 ecliptic
        (lon, lat, dist), (dlon, dlat, ddist) = reader.eval_body(ipl, jd_tt)
    else:
        # Geocentric: convert heliocentric J2000 ecliptic → geocentric J2000
        cos_e = math.cos(OBLIQUITY_J2000_RAD)
        sin_e = math.sin(OBLIQUITY_J2000_RAD)

        def _geo_j2000(jd: float) -> Tuple[float, float, float]:
            """Geocentric J2000 ecliptic position from LEB data."""
            (h_lon, h_lat, h_dist), _ = reader.eval_body(ipl, jd)
            # Heliocentric J2000 ecliptic → Cartesian
            h_lon_r = math.radians(h_lon)
            h_lat_r = math.radians(h_lat)
            cl = math.cos(h_lat_r)
            xh = h_dist * cl * math.cos(h_lon_r)
            yh = h_dist * cl * math.sin(h_lon_r)
            zh = h_dist * math.sin(h_lat_r)
            # Earth heliocentric: bary(Earth) - bary(Sun) in ICRS
            # Apply light-time correction by evaluating Earth at the retarded
            # epoch while the solar observer remains at observation time.
            earth_pos, _ = reader.eval_body(EARTH, jd)
            sun_pos, _ = reader.eval_body(SUN, jd)
            ex0 = earth_pos[0] - sun_pos[0]
            ey0 = earth_pos[1] - sun_pos[1]
            ez0 = earth_pos[2] - sun_pos[2]
            r_es = math.sqrt(ex0 * ex0 + ey0 * ey0 + ez0 * ez0)
            lt_es = r_es / C_LIGHT_AU_DAY
            earth_ret, _ = reader.eval_body(EARTH, jd - lt_es)
            ex = earth_ret[0] - sun_pos[0]
            ey = earth_ret[1] - sun_pos[1]
            ez = earth_ret[2] - sun_pos[2]
            # Rotate ICRS (≈J2000 equatorial) → J2000 ecliptic
            earth_ecl_x = ex
            earth_ecl_y = ey * cos_e + ez * sin_e
            earth_ecl_z = -ey * sin_e + ez * cos_e
            # Geocentric = body_helio - earth_helio
            xg = xh - earth_ecl_x
            yg = yh - earth_ecl_y
            zg = zh - earth_ecl_z
            rg = math.sqrt(xg * xg + yg * yg + zg * zg)
            lon_g = math.degrees(math.atan2(yg, xg)) % 360.0
            sin_b = max(-1.0, min(1.0, zg / rg)) if rg > 0 else 0.0
            lat_g = math.degrees(math.asin(sin_b))
            return lon_g, lat_g, rg

        lon, lat, dist = _geo_j2000(jd_tt)
        # Velocity via central difference (matching Skyfield path: dt=1.0 day)
        dt_v = 1.0
        prev = _geo_j2000(jd_tt - dt_v)
        nxt = _geo_j2000(jd_tt + dt_v)
        dlon = (nxt[0] - prev[0]) / (2.0 * dt_v)
        # Unwrap a 0/360 boundary crossing: the longitudes are in [0, 360), so a
        # crossing is a ~360 deg jump in the raw difference; after dividing by
        # 2*dt_v the threshold and correction carry the same 1/(2*dt_v) factor.
        # Valid only because these bodies are slow (<<90 deg/day at dt_v=1): a
        # real speed above 180/(2*dt_v) cannot occur, so it can only be a wrap.
        # Do NOT reuse this for fast bodies -- it would misread real motion.
        if dlon > 180.0 / (2.0 * dt_v):
            dlon -= 360.0 / (2.0 * dt_v)
        elif dlon < -180.0 / (2.0 * dt_v):
            dlon += 360.0 / (2.0 * dt_v)
        dlat = (nxt[1] - prev[1]) / (2.0 * dt_v)
        ddist = (nxt[2] - prev[2]) / (2.0 * dt_v)

    # Position is now J2000 ecliptic (helio or geo).
    is_j2000 = bool(iflag & FLG_J2000)
    is_equatorial = bool(iflag & FLG_EQUATORIAL)

    # For any of-date output (J2000 flag not set), precess the J2000 position
    # AND velocity to the ecliptic of date. Carrying the velocity through the
    # SAME (epoch-varying) precession — instead of leaving it frozen in the
    # J2000 frame — makes the of-date speed the true time derivative of the
    # of-date position: it picks up the ~50.29"/yr general-precession rate,
    # matching the Skyfield backend (_precess_ecliptic_state). Without this,
    # the derivative would omit the precession contribution.
    if not is_j2000:
        dt_step = 0.001  # days
        now_lon, now_lat = _precess_ecliptic(lon, lat, J2000, jd_tt)
        fwd_lon, fwd_lat = _precess_ecliptic(
            lon + dlon * dt_step, lat + dlat * dt_step, J2000, jd_tt + dt_step
        )
        d_lon = fwd_lon - now_lon
        if d_lon > 180.0:
            d_lon -= 360.0
        elif d_lon < -180.0:
            d_lon += 360.0
        dlon = d_lon / dt_step
        dlat = (fwd_lat - now_lat) / dt_step
        lon = now_lon
        lat = now_lat

    # Rotate the mean ecliptic of date onto the TRUE ecliptic of date by adding
    # nutation in longitude (Δψ), so these bodies are reported on the true
    # ecliptic like every other of-date body. The velocity picks up the
    # nutation-in-longitude
    # rate dΔψ/dt so the reported speed stays the exact derivative of the
    # reported longitude. Skipped for the nutation-free frames: FLG_J2000 (the
    # J2000 ecliptic carries no nutation), FLG_NONUT (mean ecliptic), and
    # SIDEREAL+EQUATORIAL (rotated to the equator with the mean obliquity below,
    # so the longitude must stay mean too). For plain SIDEREAL the true ayanamsha
    # (mean + Δψ) is subtracted downstream in _fast_calc_core, which then also
    # removes dΔψ/dt from the speed — the Δψ added here is what makes that
    # cancellation land on the intrinsic sidereal rate.
    _sid_eq = bool(iflag & FLG_SIDEREAL) and is_equatorial
    if not is_j2000 and not (iflag & FLG_NONUT) and not _sid_eq:
        _, dpsi_rad, _, _ = _frame_data(jd_tt)
        lon = (lon + math.degrees(dpsi_rad)) % 360.0
        if iflag & FLG_SPEED:
            dlon += _nutation_rate_deg_day(jd_tt)

    # For EQ+J2000, both backends rotate the J2000 ecliptic to the equator with
    # the J2000 obliquity — a consistent J2000 equatorial frame, same as every
    # other body class under the shared IAU 2006 frame definition.
    if is_equatorial and is_j2000:
        # EQ+J2000: J2000 ecliptic → J2000 equatorial (fixed J2000 obliquity).
        eps = OBLIQUITY_J2000_REF_DEG

        dt_step = 0.001
        eq_now_lon, eq_now_lat = _cotrans(lon, lat, -eps)
        eq_fwd_lon, eq_fwd_lat = _cotrans(
            lon + dlon * dt_step, lat + dlat * dt_step, -eps
        )
        d_eq_lon = eq_fwd_lon - eq_now_lon
        if d_eq_lon > 180.0:
            d_eq_lon -= 360.0
        elif d_eq_lon < -180.0:
            d_eq_lon += 360.0
        dlon = d_eq_lon / dt_step
        dlat = (eq_fwd_lat - eq_now_lat) / dt_step
        lon = eq_now_lon
        lat = eq_now_lat

    elif is_equatorial:
        # Equatorial of date: coords (position and velocity) are now on the
        # ecliptic of date (precessed above); rotate ecliptic-of-date →
        # equatorial-of-date with the of-date obliquity (a fixed-obliquity
        # rotation). Sidereal/NONUT modes use mean obliquity (no nutation).
        if (iflag & FLG_SIDEREAL) or (iflag & FLG_NONUT):
            eps = vondrak_mean_obliquity_deg(jd_tt)
        else:
            _, _, deps, _ = _frame_data(jd_tt)
            eps_mean = vondrak_mean_obliquity_deg(jd_tt)
            eps = eps_mean + math.degrees(deps)

        dt_step = 0.001
        eq_now_lon, eq_now_lat = _cotrans(lon, lat, -eps)
        eq_fwd_lon, eq_fwd_lat = _cotrans(
            lon + dlon * dt_step, lat + dlat * dt_step, -eps
        )
        d_eq_lon = eq_fwd_lon - eq_now_lon
        if d_eq_lon > 180.0:
            d_eq_lon -= 360.0
        elif d_eq_lon < -180.0:
            d_eq_lon += 360.0
        dlon = d_eq_lon / dt_step
        dlat = (eq_fwd_lat - eq_now_lat) / dt_step
        lon = eq_now_lon
        lat = eq_now_lat

    # else: J2000 ecliptic (already in frame) or of-date ecliptic (precessed
    # above) — position and velocity are final.

    return lon, lat, dist, dlon, dlat, ddist


# =============================================================================
# ENTRY POINTS
# =============================================================================


def _escalate_sealed_range_miss(
    reader: "LEBReaderLike", ipl: int, jd_tt: float, err: Exception
) -> Optional[Exception]:
    """Name the REQUESTED body when neither it nor its LEB support covers jd.

    In sealed ``leb`` mode a curated minor body requested outside its own LEB
    window falls through to the declared Keplerian model, which reads Sun/Earth
    states from the same LEB. When the support states are ALSO out of range the
    Keplerian reduction fails on the internal Sun lookup and would surface an
    ``EphemerisRangeError`` naming ``Body 0`` (the support Sun) instead of the
    body the caller actually asked for. Detect that both the requested body and
    the Sun are out of range and return a typed error naming the requested body
    with its own stored window.

    Returns ``None`` -- letting the original ``ValueError``/``KeyError``
    propagate so calc_ut()/calc() take the normal out-of-range → Keplerian
    fall-through -- when the support states still cover the date (the Keplerian
    model can run), when the request was not a body-range miss (e.g. an
    unsupported-flag ``KeyError`` on an in-range body), or when the reader does
    not expose per-body coverage.
    """
    from .state import get_calc_mode

    if get_calc_mode() != "leb":
        return None
    # Core bodies (Sun..mean Apogee, Earth) are handled downstream by
    # planets._raise_leb_range_miss, which emits the canonical sealed-mode
    # "does not silently substitute" contract. Only the curated minor bodies
    # that fall through to the Keplerian model (and there read the internal Sun)
    # can be mislabelled as "Body 0", so restrict escalation to non-core ids.
    from .planets import _LEB_CORE_BODY_IDS

    if ipl in _LEB_CORE_BODY_IDS:
        return None
    coverage_fn = getattr(reader, "body_coverage", None)
    if coverage_fn is None:
        return None
    body_cov = coverage_fn(ipl)
    if body_cov is None:
        # Body absent from the reader entirely: not a date-range miss for this
        # body. Leave the original error for the caller's missing-body policy.
        return None
    b_start, b_end = float(body_cov[0]), float(body_cov[1])
    if b_start <= jd_tt <= b_end:
        # The requested body IS covered here; the failure was something else
        # (unsupported flag, corrupt section). Do not relabel it.
        return None
    sun_cov = coverage_fn(SUN)
    if sun_cov is not None and float(sun_cov[0]) <= jd_tt <= float(sun_cov[1]):
        # Keplerian support (Sun/Earth) is available: allow the documented
        # out-of-range → declared-local-model fall-through to proceed.
        return None
    from .exceptions import EphemerisRangeError

    return EphemerisRangeError(
        message=(
            f"Body {ipl} at JD {jd_tt:.6f} is outside active LEB coverage range "
            f"[{b_start:.6f}, {b_end:.6f}], and no in-range LEB support state is "
            "available for the declared local model."
        ),
        requested_jd=jd_tt,
        start_jd=b_start,
        end_jd=b_end,
        body_id=ipl,
    )


def fast_calc_ut(
    reader: "LEBReaderLike",
    tjd_ut: float,
    ipl: int,
    iflag: int,
    *,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
    topo: Optional[Tuple[float, float, float]] = None,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Fast equivalent of calc_ut() using precomputed .leb data.

    Args:
        reader: An open LEBReader instance.
        tjd_ut: Julian Day in Universal Time (UT1).
        ipl: Planet/body ID (SUN, MOON, etc.).
        iflag: Calculation flags (FLG_SPEED, FLG_HELCTR, etc.).
        sid_mode: Sidereal mode override (for thread-safe context calls).
        sid_t0: Sidereal reference epoch override.
        sid_ayan_t0: Sidereal ayanamsha-at-epoch override.

    Returns:
        Same as calc_ut(): ((lon, lat, dist, dlon, dlat, ddist), iflag)

    Raises:
        KeyError: If the body/flag is not handled by this direct LEB reducer;
            the caller applies the active mode's source policy.
        ValueError: If JD is outside the .leb file's range.
    """
    # FLG_ICRS is not handled by this direct reducer. The caller's mode-aware
    # vector resolver performs it without crossing the LEB source boundary.
    if iflag & FLG_ICRS:
        raise KeyError("FLG_ICRS not supported in LEB mode")

    # Strip FLG_MOSEPH (always ignored)
    iflag = iflag & ~FLG_MOSEPH

    # FLG_TOPOCTR: an explicit topo override (context calls) wins over
    # the global state set via set_topo().
    topo_offset = None
    if iflag & FLG_TOPOCTR:
        if topo is not None:
            topo_geopos = (float(topo[0]), float(topo[1]), float(topo[2]))
        else:
            from .state import get_topo

            _topo_obj = get_topo()
            if _topo_obj is None:
                raise ValueError("set_topo() must be called before FLG_TOPOCTR")
            topo_geopos = (
                float(_topo_obj.longitude.degrees),
                float(_topo_obj.latitude.degrees),
                float(_topo_obj.elevation.m),
            )

    # Snapshot sidereal state once at entry (thread-safe)
    if sid_mode is None and (iflag & FLG_SIDEREAL):
        # Read via the locked getter so a concurrent EphemerisContext swap
        # (which installs its sidereal config into these globals under
        # _CONTEXT_SWAP_LOCK for the whole calc) cannot leak into a
        # module-level sidereal request. Reading the globals directly here
        # bypassed the lock and could return another context's ayanamsha.
        from .state import get_sid_mode

        sid_mode, sid_t0, sid_ayan_t0 = get_sid_mode(full=True)

    # Delta-T conversion: UT -> TT
    from .time_utils import deltat

    delta_t = deltat(tjd_ut)
    jd_tt = tjd_ut + delta_t

    topo_offset_fn = None
    if iflag & FLG_TOPOCTR:
        _tgp = topo_geopos  # type: ignore[possibly-undefined]
        topo_offset = _topocentric_offset(_tgp, jd_tt, tjd_ut, reader)
        if iflag & FLG_SPEED:
            # Provider used by the velocity central difference to recover the
            # observer's diurnal (Earth-rotation) motion at jd_tt ± _VEL_H.
            def topo_offset_fn(
                jd_sample: float, _dt: float = jd_tt - tjd_ut
            ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
                return _topocentric_offset(_tgp, jd_sample, jd_sample - _dt, reader)

    try:
        return _fast_calc_core(
            reader,
            jd_tt,
            tjd_ut,
            ipl,
            iflag,
            sid_mode=sid_mode,
            sid_t0=sid_t0,
            sid_ayan_t0=sid_ayan_t0,
            topo_offset=topo_offset,
            topo_offset_fn=topo_offset_fn,
        )
    except (KeyError, ValueError) as err:
        escalated = _escalate_sealed_range_miss(reader, ipl, jd_tt, err)
        if escalated is not None:
            raise escalated from err
        raise


def fast_calc_tt(
    reader: "LEBReaderLike",
    tjd_tt: float,
    ipl: int,
    iflag: int,
    *,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
    topo: Optional[Tuple[float, float, float]] = None,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Fast equivalent of calc() using precomputed .leb data.

    Args:
        reader: An open LEBReader instance.
        tjd_tt: Julian Day in Terrestrial Time (TT).
        ipl: Planet/body ID.
        iflag: Calculation flags.
        sid_mode: Sidereal mode override (for thread-safe context calls).
        sid_t0: Sidereal reference epoch override.
        sid_ayan_t0: Sidereal ayanamsha-at-epoch override.

    Returns:
        Same as calc(): ((lon, lat, dist, dlon, dlat, ddist), iflag)

    Raises:
        KeyError: If body is not in the .leb file.
        ValueError: If JD is outside the .leb file's range.
    """
    if iflag & FLG_ICRS:
        raise KeyError("FLG_ICRS not supported in LEB mode")

    iflag = iflag & ~FLG_MOSEPH

    # FLG_TOPOCTR: an explicit topo override (context calls) wins over
    # the global state set via set_topo().
    topo_offset = None
    if iflag & FLG_TOPOCTR:
        if topo is not None:
            topo_geopos = (float(topo[0]), float(topo[1]), float(topo[2]))
        else:
            from .state import get_topo

            _topo_obj = get_topo()
            if _topo_obj is None:
                raise ValueError("set_topo() must be called before FLG_TOPOCTR")
            topo_geopos = (
                float(_topo_obj.longitude.degrees),
                float(_topo_obj.latitude.degrees),
                float(_topo_obj.elevation.m),
            )

    if sid_mode is None and (iflag & FLG_SIDEREAL):
        # Read via the locked getter so a concurrent EphemerisContext swap
        # (which installs its sidereal config into these globals under
        # _CONTEXT_SWAP_LOCK for the whole calc) cannot leak into a
        # module-level sidereal request. Reading the globals directly here
        # bypassed the lock and could return another context's ayanamsha.
        from .state import get_sid_mode

        sid_mode, sid_t0, sid_ayan_t0 = get_sid_mode(full=True)

    # TT -> UT via the canonical deltat(), NOT the ΔT table baked into the LEB
    # file (which can be stale by several seconds at historical dates, e.g.
    # ~9 s near 1719 in the medium tier). Using the same deltat() as
    # fast_calc_ut() keeps the two entry points coherent (a TT jd and its
    # UT partner give the same topocentric / sidereal result). deltat() takes
    # UT; invert with one fixed-point step (ΔT varies slowly, so it converges
    # to well under a millisecond). Fall back to the file table only if the
    # canonical deltat() is unavailable.
    from .time_utils import deltat

    try:
        _dt_days = deltat(tjd_tt)
        tjd_ut = tjd_tt - deltat(tjd_tt - _dt_days)
    except (KeyError, ValueError):
        tjd_ut = tjd_tt - reader.delta_t(tjd_tt)

    topo_offset_fn = None
    if iflag & FLG_TOPOCTR:
        _tgp = topo_geopos  # type: ignore[possibly-undefined]
        topo_offset = _topocentric_offset(_tgp, tjd_tt, tjd_ut, reader)
        if iflag & FLG_SPEED:
            # Provider used by the velocity central difference to recover the
            # observer's diurnal (Earth-rotation) motion at jd_tt ± _VEL_H.
            def topo_offset_fn(
                jd_sample: float, _dt: float = tjd_tt - tjd_ut
            ) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
                return _topocentric_offset(_tgp, jd_sample, jd_sample - _dt, reader)

    try:
        return _fast_calc_core(
            reader,
            tjd_tt,
            tjd_ut,
            ipl,
            iflag,
            sid_mode=sid_mode,
            sid_t0=sid_t0,
            sid_ayan_t0=sid_ayan_t0,
            topo_offset=topo_offset,
            topo_offset_fn=topo_offset_fn,
        )
    except (KeyError, ValueError) as err:
        escalated = _escalate_sealed_range_miss(reader, ipl, tjd_tt, err)
        if escalated is not None:
            raise escalated from err
        raise


def _fast_calc_core(
    reader: "LEBReaderLike",
    jd_tt: float,
    jd_ut: float,
    ipl: int,
    iflag: int,
    *,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
    topo_offset: Optional[
        Tuple[Tuple[float, float, float], Tuple[float, float, float]]
    ] = None,
    topo_offset_fn: Optional[
        Callable[[float], Tuple[Tuple[float, float, float], Tuple[float, float, float]]]
    ] = None,
) -> Tuple[Tuple[float, float, float, float, float, float], int]:
    """Core fast calculation logic shared by fast_calc_ut and fast_calc_tt.

    Args:
        reader: Open LEBReader.
        jd_tt: Julian Day TT (for position computation).
        jd_ut: Julian Day UT (for sidereal ayanamsa).
        ipl: Body ID.
        iflag: Flags.
        sid_mode: Sidereal mode (already snapshot at entry point).
        sid_t0: Sidereal reference epoch.
        sid_ayan_t0: Sidereal ayanamsha at reference epoch.

    Returns:
        ((lon, lat, dist, dlon, dlat, ddist), iflag)
    """
    # Set active reader for frame data dispatch (avoids threading reader
    # through every coordinate-transform helper).
    _set_active_reader(reader)

    output_iflag = iflag
    if FICT_OFFSET <= ipl <= WALDEMATH:
        # Never source a fictitious-body position from a persisted LEB channel.
        # Legacy files may contain obsolete coefficients, while all public
        # IDs 40--58 now have runtime analytical dispatch. A KeyError is the
        # normal LEB miss signal and lets calc/calc_ut continue to that path.
        raise KeyError(
            f"Fictitious body {ipl} is calculated from its runtime analytical model"
        )
    # FLG_J2000 is honored uniformly for every lunar point, including the
    # non-mean ones under SIDEREAL: ayanamsha and ecliptic-plane precession
    # are independent, composable rotations, so the deferred sid+j2k pipeline
    # applies to all ecliptic-direct bodies alike.

    # Clean-room compatibility models do not depend on legacy LEB channels.
    # They remain available even when a custom LEB omits those body records.
    _external_ecliptic_model = ipl in (
        MEAN_NODE,
        MEAN_APOG,
        INTP_APOG,
        INTP_PERG,
    )
    if not _external_ecliptic_model and not reader.has_body(ipl):
        raise KeyError(f"Body {ipl} not in LEB file")

    # Moon-derived abstract points have no heliocentric or barycentric
    # geometric definition; use the public zero sentinel for those centers.
    if ipl in (MEAN_NODE, TRUE_NODE, MEAN_APOG, OSCU_APOG, INTP_APOG, INTP_PERG) and (
        iflag & (FLG_HELCTR | FLG_BARYCTR)
    ):
        return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0), output_iflag

    coord_type = (
        COORD_ECLIPTIC if _external_ecliptic_model else reader._bodies[ipl].coord_type
    )

    _pipeline_a = False
    _pipeline_c = False

    # Dispatch to appropriate pipeline based on coordinate type.
    # When XYZ+SIDEREAL, force Pipeline A to return spherical so the
    # shared sidereal correction (longitude subtraction) can run first;
    # XYZ conversion then happens in the post-processing stage.
    _xyz_sid = bool(iflag & FLG_XYZ) and bool(iflag & FLG_SIDEREAL)
    _pipe_iflag = (iflag & ~FLG_XYZ) if _xyz_sid else iflag

    if coord_type == COORD_ICRS_BARY:
        _pipeline_a = not _xyz_sid
        if iflag & FLG_SPEED:
            lon, lat, dist, dlon, dlat, ddist = _pipeline_icrs(
                reader,
                jd_tt,
                ipl,
                _pipe_iflag,
                want_velocity=True,
                topo_offset=topo_offset,
                topo_offset_fn=topo_offset_fn,
            )
        else:
            lon, lat, dist = _pipeline_icrs(
                reader,
                jd_tt,
                ipl,
                _pipe_iflag,
                topo_offset=topo_offset,
            )
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    elif coord_type == COORD_ICRS_BARY_SYSTEM:
        _pipeline_a = not _xyz_sid
        if iflag & FLG_SPEED:
            lon, lat, dist, dlon, dlat, ddist = _pipeline_icrs(
                reader,
                jd_tt,
                ipl,
                _pipe_iflag,
                want_velocity=True,
                is_system_bary=True,
                topo_offset=topo_offset,
                topo_offset_fn=topo_offset_fn,
            )
        else:
            lon, lat, dist = _pipeline_icrs(
                reader,
                jd_tt,
                ipl,
                _pipe_iflag,
                is_system_bary=True,
                topo_offset=topo_offset,
            )
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    elif coord_type == COORD_ECLIPTIC:
        if (iflag & FLG_TOPOCTR) and ipl not in (
            MEAN_NODE,
            TRUE_NODE,
            MEAN_APOG,
            OSCU_APOG,
            INTP_APOG,
            INTP_PERG,
        ):
            raise KeyError("FLG_TOPOCTR not supported for ecliptic-direct bodies")
        # Pipeline B: ecliptic direct
        lon, lat, dist, dlon, dlat, ddist = _pipeline_ecliptic(
            reader, jd_tt, ipl, iflag
        )
        if not (iflag & FLG_SPEED):
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    elif coord_type == COORD_HELIO_ECL:
        if iflag & FLG_TOPOCTR:
            raise KeyError("FLG_TOPOCTR not supported for heliocentric bodies")
        # Pipeline C: heliocentric
        _pipeline_c = True
        lon, lat, dist, dlon, dlat, ddist = _pipeline_helio(reader, jd_tt, ipl, iflag)
        if not (iflag & FLG_SPEED):
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    else:
        raise ValueError(f"Unknown coord_type {coord_type}")

    # ICRS (Pipeline-A) bodies requested as XYZ+SIDEREAL are forced through the
    # spherical path (FLG_XYZ stripped via _pipe_iflag) so the sidereal longitude
    # subtraction can run; the Cartesian conversion then happens at the end. At
    # that point dlon is still a spherical longitude rate exactly like the
    # Pipeline-A spherical path, so the J2000 / ayanamsa-drift speed gates below
    # must treat them identically — _pipeline_a was set False only to route the
    # XYZ conversion. Without this, the XYZ sidereal-J2000 velocity omits the two
    # precession-rate terms the equivalent spherical output applies, leaving the
    # XYZ and spherical results (which must be exact coordinate transforms of one
    # another) inconsistent by ~2x the general-precession rate in dlon.
    _xyz_sid_pipea = _xyz_sid and coord_type in (
        COORD_ICRS_BARY,
        COORD_ICRS_BARY_SYSTEM,
    )

    # Sidereal correction
    #
    # Public API convention: EQUATORIAL output does not subtract ayanamsha
    # from right ascension for any body type:
    #   - Pipeline A (ICRS): already outputs equatorial coordinates using the
    #     mean equator (P matrix) when sidereal is set.
    #   - Pipeline B/C (ecliptic/helio): the pipeline handles equatorial
    #     conversion internally; sidereal correction is not applied.
    # So skip sidereal correction entirely when EQUATORIAL is requested.
    #
    # For J2000 output, coordinates are already in J2000 ecliptic (no nutation),
    # so we use mean ayanamsha; for ecliptic of date we use true ayanamsha.
    _skip_sidereal = bool(iflag & FLG_EQUATORIAL)

    if (iflag & FLG_SIDEREAL) and not _skip_sidereal:
        try:
            mean_aya = _calc_ayanamsa_from_leb(
                reader,
                jd_tt,
                sid_mode=sid_mode,
                sid_t0=sid_t0,
                sid_ayan_t0=sid_ayan_t0,
            )
            # J2000 ecliptic has no nutation component → mean ayanamsha.
            # Ecliptic of date includes nutation → true ayanamsha (mean + Δψ).
            #
            # FLG_J2000 is honored for every lunar point, so any body with the
            # J2000 bit set uses the mean ayanamsha after precession.
            _eff_mean_aya = bool(iflag & (FLG_J2000 | FLG_NONUT))
            if _eff_mean_aya:
                aya = mean_aya
            else:
                try:
                    _, dpsi_sid, _, _ = _frame_data(jd_tt)
                    nutation_deg = math.degrees(dpsi_sid)
                    aya = mean_aya + nutation_deg
                except (KeyError, ValueError) as _nut_exc:
                    if _external_ecliptic_model:
                        raise KeyError(
                            "LEB frame data do not cover this clean-model epoch"
                        ) from _nut_exc
                    from .logging_config import get_logger

                    # Mean-only ayanamsha is a ~9-17" silent precision loss;
                    # never hide it, and never eat programming errors.
                    get_logger().warning(
                        "sidereal: nutation unavailable (%s); using mean ayanamsha",
                        _nut_exc,
                    )
                    aya = mean_aya

            lon = (lon - aya) % 360.0

            # Sidereal speed correction (ayanamsa drift): the ayanamsha drifts
            # ~50"/yr, so its time-derivative (~the general-precession rate) is
            # removed from dlon for sidereal SPEED output. _PREC_COEFFS are
            # arcsec/century: dP/dT = c0 + 2*c1*T + ...; convert to deg/day with
            # / 3600 / 36525.
            #
            # Fires for every sidereal SPEED request — including FLG_J2000. The
            # velocity coming out of the pipeline is the true derivative of the
            # tropical position in the output frame (of-date ecliptic, or J2000
            # ecliptic under FLG_J2000). The sidereal longitude subtracts the
            # ayanamsha, so its rate subtracts the ayanamsha drift (~ the general
            # precession rate) in BOTH frames: in the of-date frame the tropical
            # rate carries precession which is (partly) cancelled; in the J2000
            # frame the tropical rate has no equinox motion but the mean
            # ayanamsha still drifts, so the term is the full ~p. The
            # FLG_SPEED-absent case must be excluded: dlon was zeroed and stays
            # 0.0.
            if iflag & FLG_SPEED:
                if sid_mode in _STAR_BASED_MODES:
                    # Star-anchored modes: ayanamsha = anchor-star apparent
                    # longitude − C, so its drift is the anchor star's
                    # apparent-longitude drift (dominated by the annual-aberration
                    # derivative), NOT the IAU 2006 precession polynomial. The
                    # position correction already delegates the ayanamsha value
                    # to Skyfield for these modes; do the same for the rate by
                    # central-differencing the actual ayanamsha (mean or true per
                    # the active flags), as the Skyfield path's
                    # _apply_sidereal_correction does. A fixed precession rate
                    # cannot represent the apparent motion of the anchor.
                    def _aya_of(_jd: float) -> float:
                        _m = _calc_ayanamsa_from_leb(
                            reader,
                            _jd,
                            sid_mode=sid_mode,
                            sid_t0=sid_t0,
                            sid_ayan_t0=sid_ayan_t0,
                        )
                        if _eff_mean_aya:
                            return _m
                        try:
                            _, _dpsi, _, _ = _frame_data(_jd)
                            return _m + math.degrees(_dpsi)
                        except (KeyError, ValueError):
                            return _m

                    # The ayanamsha is reported mod 360, so the two samples can
                    # straddle the 0/360 wrap; unwrap the delta into
                    # (-180, 180] before dividing or the raw -360 difference
                    # reports a ~1.6e7 deg/day speed spike.
                    _dt_aya = 1.0 / 86400.0
                    _d_aya = (
                        _aya_of(jd_tt + _dt_aya) - _aya_of(jd_tt - _dt_aya) + 180.0
                    ) % 360.0 - 180.0
                    dlon -= _d_aya / (2.0 * _dt_aya)
                else:
                    dlon -= _general_precession_rate_deg_day(jd_tt)

                    # Of-date sidereal longitude subtracts the true ayanamsha
                    # (mean + Δψ), so its derivative must also subtract dΔψ/dt.
                    # Mean-ayanamsha frames (J2000/NONUT) carry no nutation and
                    # are excluded. This keeps speed equal to the derivative of
                    # the reported position for every ecliptic-direct body.
                    if not _eff_mean_aya:
                        try:
                            dlon -= _nutation_rate_deg_day(jd_tt)
                        except (KeyError, ValueError):
                            pass

        except KeyError:
            # Star-based sidereal mode, fall back
            raise

    # NOTE: Pipeline-A J2000 spherical output needs NO extra precession-rate
    # subtraction here. _pipeline_icrs now central-differences the apparent
    # place through the J2000 frame transform, so the reported dlon is already
    # the exact derivative of the reported J2000 longitude (the J2000 frame is
    # fixed, so there is no equinox motion to remove). The former subtraction
    # here double-counted it against XYZ output and, together with the sidereal
    # drift above, against SID+J2000.

    # FLG_XYZ post-processing for Pipeline B/C (spherical → Cartesian).
    # Pipeline A handles XYZ internally via _pipeline_icrs.
    if (iflag & FLG_XYZ) and not _pipeline_a:
        result = _spherical_to_cartesian_with_velocity(
            lon, lat, dist, dlon, dlat, ddist
        )
    elif (iflag & FLG_RADIANS) and not (iflag & FLG_XYZ):
        # FLG_RADIANS: convert angular output from degrees to radians.
        _d2r = math.pi / 180.0
        result = (
            lon * _d2r,
            lat * _d2r,
            dist,
            dlon * _d2r,
            dlat * _d2r,
            ddist,
        )
    else:
        result = (lon, lat, dist, dlon, dlat, ddist)

    # Mean lunar points apply SPEED3 after every requested frame and output
    # transformation. Reusing this core with velocity bits removed gives the
    # exact reported positions at t±h without recursively rebuilding a speed.
    if (iflag & FLG_SPEED3) and ipl in (MEAN_NODE, MEAN_APOG):
        sample_iflag = iflag & ~(FLG_SPEED | FLG_SPEED3)
        previous, _ = _fast_calc_core(
            reader,
            jd_tt - _SPEED3_STEP_DAYS,
            jd_ut - _SPEED3_STEP_DAYS,
            ipl,
            sample_iflag,
            sid_mode=sid_mode,
            sid_t0=sid_t0,
            sid_ayan_t0=sid_ayan_t0,
            topo_offset=topo_offset,
            topo_offset_fn=topo_offset_fn,
        )
        following, _ = _fast_calc_core(
            reader,
            jd_tt + _SPEED3_STEP_DAYS,
            jd_ut + _SPEED3_STEP_DAYS,
            ipl,
            sample_iflag,
            sid_mode=sid_mode,
            sid_t0=sid_t0,
            sid_ayan_t0=sid_ayan_t0,
            topo_offset=topo_offset,
            topo_offset_fn=topo_offset_fn,
        )
        result = _apply_central_speed3_stencil(result, previous, following, iflag)

    return result, output_iflag
