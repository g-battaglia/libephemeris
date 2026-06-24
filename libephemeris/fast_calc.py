# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Fast calculation pipeline using precomputed .leb binary ephemeris data.

This module reimplements the Skyfield pipeline using LEBReader as the data
source, handling coordinate transforms, light-time correction, aberration,
and flag dispatch.

Three pipelines:
    A: ICRS barycentric bodies (Sun-Pluto, Earth, Chiron, Ceres-Vesta)
    B: Ecliptic direct bodies (lunar nodes, Lilith variants)
    C: Heliocentric bodies (Uranians, Transpluto)
"""

from __future__ import annotations

import math
import threading
from functools import lru_cache
from typing import TYPE_CHECKING, Dict, Optional, Tuple

from .constants import (
    EARTH,
    INTP_APOG,
    INTP_PERG,
    MEAN_APOG,
    MEAN_NODE,
    MOON,
    OSCU_APOG,
    SUN,
    TRUE_NODE,
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
    FLG_TOPOCTR,
    FLG_TRUEPOS,
    FLG_XYZ,
    _MOON_MEAN_DIST_AU,
    _MOON_MEAN_APOG_DIST_AU,
)
from .leb_format import (
    COORD_ECLIPTIC,
    COORD_HELIO_ECL,
    COORD_ICRS_BARY,
    COORD_ICRS_BARY_SYSTEM,
)
from .precession_vondrak import (
    vondrak_mean_obliquity_deg,
    vondrak_pn_matrix,
    vondrak_precession_matrix,
)

if TYPE_CHECKING:
    from typing import Union

    from .leb2_reader import LEB2Reader
    from .leb_composite import CompositeLEBReader
    from .leb_reader import LEBReader

    LEBReaderLike = Union[LEBReader, LEB2Reader, CompositeLEBReader]

# =============================================================================
# CONSTANTS
# =============================================================================

C_LIGHT_AU_DAY = 173.1446326846693  # Speed of light in AU/day
J2000 = 2451545.0  # J2000.0 epoch in JD
OBLIQUITY_J2000_DEG = 23.4392911  # Mean obliquity at J2000 (degrees)
OBLIQUITY_J2000_RAD = math.radians(OBLIQUITY_J2000_DEG)

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
    return (x, y, z, vx, vy, vz)


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
    """Rotate ICRS (≈equatorial J2000) to ecliptic J2000."""
    return _rotate_equatorial_to_ecliptic(x, y, z, OBLIQUITY_J2000_RAD)


def _apply_aberration(
    geo: Tuple[float, float, float],
    earth_vel: Tuple[float, float, float],
    light_time: float = 0.0,
) -> Tuple[float, float, float]:
    """Apply relativistic aberration to a geometric position vector.

    Uses the full special-relativistic formula matching Skyfield's
    ``add_aberration()``.  This includes the Lorentz factor γ for
    exact agreement with the reference pipeline (< 0.001 mas).

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


# Gravitational deflection constants (matching Skyfield's relativity module)
_GS = 1.32712440017987e20  # heliocentric gravitational constant (m^3/s^2)
_C_MS = 299792458.0  # speed of light (m/s)
_AU_M = 149597870700  # 1 AU in metres

# Deflector reciprocal masses (solar mass / deflector mass)
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


def _apply_cob_correction(
    pos: Tuple[float, float, float],
    ipl: int,
    jd_tt: float,
) -> Tuple[float, float, float]:
    """Apply center-of-body correction to a system barycenter position.

    Converts system barycenter ICRS position to planet center ICRS position
    by adding the COB offset. Uses SPK planet_centers segments where available
    (high precision), falling back to analytical moon-theory COB corrections.

    This matches the behavior of _SpkCenterTarget in planets.py.

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
        return pos  # Not an outer planet, no COB needed

    t = get_cached_time_tt(jd_tt)

    # Map body_id to planet name for NAIF lookup
    planet_name = {5: "jupiter", 6: "saturn", 7: "uranus", 8: "neptune", 9: "pluto"}[
        ipl
    ]

    # Try SPK center offset first (high precision)
    if planet_name in _PLANET_CENTER_NAIF_IDS:
        naif_id = _PLANET_CENTER_NAIF_IDS[planet_name]
        seg = get_planet_center_segment(naif_id)
        if seg is not None:
            try:
                offset_pos = seg.at(t).position.au
                return (
                    pos[0] + float(offset_pos[0]),
                    pos[1] + float(offset_pos[1]),
                    pos[2] + float(offset_pos[2]),
                )
            except (ValueError, ArithmeticError):
                pass  # Fall through to analytical COB

    # Fallback: analytical COB from moon theories
    from .moon_theories import get_cob_offset

    offset = get_cob_offset(bary_name, t)
    return (
        pos[0] + offset[0],
        pos[1] + offset[1],
        pos[2] + offset[2],
    )


def _apply_gravitational_deflection(
    geo: Tuple[float, float, float],
    earth_bary: Tuple[float, float, float],
    jd_tt: float,
    light_time: float,
    reader: "LEBReaderLike",
) -> Tuple[float, float, float]:
    """Apply PPN gravitational light deflection by Sun, Jupiter, Saturn.

    Matches Skyfield's ``apparent(deflectors=(10, 599, 699))`` formula.
    Uses LEB data for deflector positions (barycentric ICRS).

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

        if qmag == 0.0 or emag == 0.0:
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

    Used for sidereal+equatorial output, where the reference ephemeris uses the mean equator.
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
    return (
        getattr(_active_local, "gen", -1) == _active_generation
        and getattr(_active_local, "has_nutation", False)
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
    eps_rad = math.radians(-eps)  # Negate to match the reference convention
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


# Ayanamsha J2000 offsets for formula-based sidereal modes (degrees).
# Star-based / galactic modes have (0.0) placeholders and require Skyfield.
# These values mirror the ayanamsha_data dict in planets._calc_ayanamsa().
_AYANAMSHA_J2000: Dict[int, float] = {
    0: 24.740300,  # SIDM_FAGAN_BRADLEY
    1: 23.857092,  # SIDM_LAHIRI
    2: 27.815753,  # SIDM_DELUCE
    3: 22.410791,  # SIDM_RAMAN
    4: 20.057541,  # SIDM_USHASHASHI
    5: 23.760240,  # SIDM_KRISHNAMURTI
    6: 28.359679,  # SIDM_DJWHAL_KHUL
    7: 22.478803,  # SIDM_YUKTESHWAR
    8: 22.762137,  # SIDM_JN_BHASIN
    9: 23.533640,  # SIDM_BABYL_KUGLER1
    10: 24.933640,  # SIDM_BABYL_KUGLER2
    11: 25.783640,  # SIDM_BABYL_KUGLER3
    12: 24.733640,  # SIDM_BABYL_HUBER
    13: 24.522528,  # SIDM_BABYL_ETPSC
    14: 24.758924,  # SIDM_ALDEBARAN_15TAU
    15: 20.247788,  # SIDM_HIPPARCHOS
    16: 19.992959,  # SIDM_SASSANIAN
    18: 0.0,  # SIDM_J2000
    19: 1.396581,  # SIDM_J1900
    20: 0.698370,  # SIDM_B1950
    21: 20.895059,  # SIDM_SURYASIDDHANTA
    22: 20.680425,  # SIDM_SURYASIDDHANTA_MSUN
    23: 20.895060,  # SIDM_ARYABHATA
    24: 20.657427,  # SIDM_ARYABHATA_MSUN
    25: 20.103388,  # SIDM_SS_REVATI
    26: 23.005763,  # SIDM_SS_CITRA
    37: 20.575847,  # SIDM_ARYABHATA_522
    38: 24.615753,  # SIDM_BABYL_BRITTON
    41: 25.000019,  # SIDM_GALEQU_FIORENZA
}

# Star-based modes that cannot be computed without Skyfield
_STAR_BASED_MODES = frozenset({17, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 39, 40, 42})


def _calc_ayanamsa_from_leb(
    reader: "LEBReaderLike",
    jd_tt: float,
    sid_mode: Optional[int] = None,
    sid_t0: Optional[float] = None,
    sid_ayan_t0: Optional[float] = None,
) -> float:
    """Compute ayanamsa using only .leb data (no Skyfield).

    For formula-based ayanamsha modes only. Star-based and galactic modes
    require Skyfield and will raise KeyError (triggering fallback).

    Args:
        reader: Open LEBReader instance.
        jd_tt: Julian Day in TT.
        sid_mode: Sidereal mode ID (if None, reads from global state).
        sid_t0: Reference epoch JD for custom ayanamsha.
        sid_ayan_t0: Ayanamsha value at reference epoch (degrees).

    Returns:
        Ayanamsa in degrees.

    Raises:
        KeyError: If the active sidereal mode requires Skyfield.
    """
    if sid_mode is None:
        from .state import _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0

        sid_mode = _SIDEREAL_MODE
        sid_t0 = _SIDEREAL_T0
        sid_ayan_t0 = _SIDEREAL_AYAN_T0

    T = (jd_tt - J2000) / 36525.0

    # IAU 2006 general precession
    precession_arcsec = 0.0
    T_power = T
    for coeff in _PREC_COEFFS:
        precession_arcsec += coeff * T_power
        T_power *= T

    # Get reference offset for active sidereal mode
    # Default Fagan/Bradley (0), matching the reference API
    mode = sid_mode if sid_mode is not None else 0

    if mode in _STAR_BASED_MODES:
        raise KeyError(f"Star-based sidereal mode {mode} requires Skyfield")

    if mode == 255:
        # SIDM_USER: custom user-defined ayanamsha
        ayan_t0 = sid_ayan_t0 if sid_ayan_t0 is not None else 0.0
        t0 = sid_t0 if sid_t0 is not None else J2000
        T0 = (t0 - J2000) / 36525.0

        # Delta precession from user epoch to current epoch
        def _prec(Tc: float) -> float:
            return sum(c * Tc ** (i + 1) for i, c in enumerate(_PREC_COEFFS))

        delta_prec = _prec(T) - _prec(T0)
        mean_aya = ayan_t0 + delta_prec / 3600.0
    elif mode in _AYANAMSHA_J2000:
        aya_j2000 = _AYANAMSHA_J2000[mode]
        mean_aya = aya_j2000 + precession_arcsec / 3600.0
    else:
        raise KeyError(f"Unknown sidereal mode {mode}")

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
    c2t = erfa.c2t06a(2451545.0, jd_tt - 2451545.0, 2451545.0, jd_ut1 - 2451545.0, 0.0, 0.0)

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
        (t2c_0[0] * vx_itrf + t2c_0[1] * vy_itrf + t2c_0[2] * vz_itrf) * _SEC_PER_DAY / _AU_M,
        (t2c_1[0] * vx_itrf + t2c_1[1] * vy_itrf + t2c_1[2] * vz_itrf) * _SEC_PER_DAY / _AU_M,
        (t2c_2[0] * vx_itrf + t2c_2[1] * vy_itrf + t2c_2[2] * vz_itrf) * _SEC_PER_DAY / _AU_M,
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


def _pipeline_icrs(
    reader: "LEBReaderLike",
    jd_tt: float,
    ipl: int,
    iflag: int,
    want_velocity: bool = False,
    is_system_bary: bool = False,
    topo_offset: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
    want_xyz: bool = False,
) -> Tuple[float, ...]:
    """Pipeline A: compute ecliptic coordinates for ICRS barycentric bodies.

    Handles both planet-center bodies (COORD_ICRS_BARY) and system-barycenter
    bodies (COORD_ICRS_BARY_SYSTEM). For system barycenters, COB correction
    is applied at runtime to match Skyfield's _SpkCenterTarget behavior.

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

    # 1. Get body position (and velocity if needed)
    target_pos, target_vel = reader.eval_body(ipl, jd_tt)

    # 1b. For system barycenters, apply COB only for TRUEPOS (no light-time).
    #     For normal path, COB is deferred until after light-time iteration
    #     to match Skyfield's _SpkCenterTarget._observe_from_bcrs() behavior:
    #     iterate light-time on barycenter, apply COB once at retarded time.
    if is_system_bary and (iflag & FLG_TRUEPOS):
        target_pos = _apply_cob_correction(target_pos, ipl, jd_tt)

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

    # 4. Light-time correction (unless FLG_TRUEPOS)
    #    For system barycenters: iterate on raw barycenter positions (smooth),
    #    then apply COB once after convergence at OBSERVER time. This matches
    #    Skyfield's _SpkCenterTarget._observe_from_bcrs() which:
    #      1. Calls barycenter._observe_from_bcrs(observer) to iterate light-time
    #         on the barycenter, returning (pos, vel, t, light_time) where
    #         t = observer.t (the OBSERVATION time, not retarded time).
    #      2. Evaluates center_segment.at(t) at OBSERVER time, not retarded time.
    #    We must match this: COB offset evaluated at jd_tt (observer time).
    retarded_vel = target_vel
    lt = 0.0
    if not (iflag & FLG_TRUEPOS):
        for _ in range(3):  # Fixed-point iterations
            dist = _vec3_dist(geo)
            if dist == 0.0:
                break
            lt = dist / C_LIGHT_AU_DAY
            retarded_pos, retarded_vel = reader.eval_body(ipl, jd_tt - lt)
            geo = _vec3_sub(retarded_pos, observer)

        # Apply COB correction at OBSERVER time (jd_tt), matching Skyfield's
        # _SpkCenterTarget._observe_from_bcrs() which evaluates the center
        # segment at observer.t, not at the retarded time.
        if is_system_bary and lt > 0.0:
            retarded_pos_cob = _apply_cob_correction(
                (geo[0] + observer[0], geo[1] + observer[1], geo[2] + observer[2]),
                ipl,
                jd_tt,
            )
            geo = _vec3_sub(retarded_pos_cob, observer)

    if want_velocity:
        geo_vel = _vec3_sub(retarded_vel, observer_vel)

    # 5. Gravitational deflection by Sun, Jupiter, Saturn (PPN formula).
    #    Dominant correction: up to ~4" for Saturn near the Sun's limb.
    #    Skipped for helio/bary/truepos/nogdefl and for the Moon (negligible at
    #    ~0.0026 AU, deflection < 0.000001"). NOABERR deliberately does NOT
    #    skip deflection: the reference API disables only aberration with it
    #    (FLG_ASTROMETRIC = NOABERR|NOGDEFL).
    if not (iflag & (FLG_NOGDEFL | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS)):
        if ipl != MOON and lt > 0.0:
            geo = _apply_gravitational_deflection(geo, observer, jd_tt, lt, reader)

    # 6. Aberration (full special-relativistic, matching Skyfield).
    #    For velocity: the aberration correction depends on Earth's velocity
    #    which changes slowly (~0.017 deg/day²). The velocity component of
    #    aberration is ~1e-8 deg/day — negligible. We skip it.
    if not (iflag & (FLG_NOABERR | FLG_HELCTR | FLG_BARYCTR | FLG_TRUEPOS)):
        geo = _apply_aberration(geo, earth_vel, lt)

    # 6. Coordinate transform — apply the same transform to velocity
    #
    # Sidereal handling for Pipeline A (ICRS bodies):
    #   SID+EQ: use mean equator (P matrix, no nutation) — NO ayanamsha subtraction
    #   SID+EQ+J2K: same as non-sidereal EQ+J2K (ICRS ≡ J2000 equatorial)
    #   SID only: output ecliptic-of-date, ayanamsha subtracted in _fast_calc_core
    #   SID+J2K: output J2000 ecliptic, ayanamsha subtracted in _fast_calc_core
    _is_sidereal = bool(iflag & FLG_SIDEREAL)
    _want_nonut = bool(iflag & FLG_NONUT)

    _want_xyz = want_xyz or bool(iflag & FLG_XYZ)

    if (iflag & FLG_EQUATORIAL) and (iflag & FLG_J2000):
        # ICRS J2000 equatorial -- geo is already in this frame
        if _want_xyz:
            lon_deg, lat_deg, dist = geo[0], geo[1], geo[2]
        else:
            lon_deg, lat_deg, dist = _cartesian_to_spherical(geo[0], geo[1], geo[2])
        if want_velocity:
            if _want_xyz:
                dlon, dlat, ddist = geo_vel[0], geo_vel[1], geo_vel[2]
            else:
                dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                    geo[0], geo[1], geo[2],
                    geo_vel[0], geo_vel[1], geo_vel[2],
                )

    elif (iflag & FLG_EQUATORIAL) and _is_sidereal:
        p_mat = _prec_matrix(jd_tt)
        geo_eq = _mat3_vec3(p_mat, geo)
        if _want_xyz:
            lon_deg, lat_deg, dist = geo_eq[0], geo_eq[1], geo_eq[2]
        else:
            lon_deg, lat_deg, dist = _cartesian_to_spherical(
                geo_eq[0], geo_eq[1], geo_eq[2]
            )
        if want_velocity:
            vel_eq = _mat3_vec3(p_mat, geo_vel)
            if _want_xyz:
                dlon, dlat, ddist = vel_eq[0], vel_eq[1], vel_eq[2]
            else:
                dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                    geo_eq[0], geo_eq[1], geo_eq[2],
                    vel_eq[0], vel_eq[1], vel_eq[2],
                )

    elif iflag & FLG_EQUATORIAL:
        # NONUT: use P matrix (mean equator), else PNM (true equator)
        if _want_nonut:
            _eq_mat = _prec_matrix(jd_tt)
        else:
            _eq_mat, _, _, _ = _frame_data(jd_tt)
        geo_eq = _mat3_vec3(_eq_mat, geo)
        if _want_xyz:
            lon_deg, lat_deg, dist = geo_eq[0], geo_eq[1], geo_eq[2]
        else:
            lon_deg, lat_deg, dist = _cartesian_to_spherical(
                geo_eq[0], geo_eq[1], geo_eq[2]
            )
        if want_velocity:
            vel_eq = _mat3_vec3(_eq_mat, geo_vel)
            if _want_xyz:
                dlon, dlat, ddist = vel_eq[0], vel_eq[1], vel_eq[2]
            else:
                dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                    geo_eq[0], geo_eq[1], geo_eq[2],
                    vel_eq[0], vel_eq[1], vel_eq[2],
                )

    elif iflag & FLG_J2000:
        ecl = _rotate_icrs_to_ecliptic_j2000(geo[0], geo[1], geo[2])
        if _want_xyz:
            lon_deg, lat_deg, dist = ecl[0], ecl[1], ecl[2]
        else:
            lon_deg, lat_deg, dist = _cartesian_to_spherical(ecl[0], ecl[1], ecl[2])
        if want_velocity:
            vel_ecl = _rotate_icrs_to_ecliptic_j2000(
                geo_vel[0], geo_vel[1], geo_vel[2]
            )
            if _want_xyz:
                dlon, dlat, ddist = vel_ecl[0], vel_ecl[1], vel_ecl[2]
            else:
                dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                    ecl[0], ecl[1], ecl[2],
                    vel_ecl[0], vel_ecl[1], vel_ecl[2],
                )

    else:
        # TRUE ECLIPTIC OF DATE (default) -- most common path
        # FLG_NONUT: mean ecliptic (P matrix only, no nutation)
        if _want_nonut:
            _rot_mat = _prec_matrix(jd_tt)
            eps_rad = math.radians(vondrak_mean_obliquity_deg(jd_tt))
        else:
            _rot_mat, _, _, eps_rad = _frame_data(jd_tt)

        geo_eq = _mat3_vec3(_rot_mat, geo)

        ecl = _rotate_equatorial_to_ecliptic(
            geo_eq[0], geo_eq[1], geo_eq[2], eps_rad
        )
        if _want_xyz:
            lon_deg, lat_deg, dist = ecl[0], ecl[1], ecl[2]
        else:
            lon_deg, lat_deg, dist = _cartesian_to_spherical(
                ecl[0], ecl[1], ecl[2]
            )

        if want_velocity:
            vel_eq = _mat3_vec3(_rot_mat, geo_vel)
            vel_ecl = _rotate_equatorial_to_ecliptic(
                vel_eq[0], vel_eq[1], vel_eq[2], eps_rad
            )
            if _want_xyz:
                dlon, dlat, ddist = vel_ecl[0], vel_ecl[1], vel_ecl[2]
            else:
                dlon, dlat, ddist = _cartesian_velocity_to_spherical(
                    ecl[0], ecl[1], ecl[2],
                    vel_ecl[0], vel_ecl[1], vel_ecl[2],
                )

    if want_velocity:
        return lon_deg, lat_deg, dist, dlon, dlat, ddist
    return lon_deg, lat_deg, dist


# =============================================================================
# PIPELINE B: ECLIPTIC DIRECT BODIES
# =============================================================================


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
    (lon, lat, dist), (dlon, dlat, ddist) = reader.eval_body(ipl, jd_tt)

    # Override distance/latitude for lunar analytical bodies.
    # LEB stores pre-computed Chebyshev coefficients that may use older
    # models.  These overrides compute values analytically at runtime
    # for better match with the reference ephemeris.
    if ipl == MEAN_NODE:
        # LEB stores 0 for dist; the reference ephemeris returns mean distance constant.
        dist = _MOON_MEAN_DIST_AU
    elif ipl == TRUE_NODE:
        # KNOWN LIMITATION (see divergences.md 7): the LEB-stored value is an
        # h_mag proxy (~0.0015 AU), NOT the osculating node radius the reference ephemeris and
        # the Skyfield backend return (~0.0024 AU). Recomputing it standalone
        # would require reconstructing the geocentric-ecliptic Moon state inside
        # this pipeline; the longitude/latitude (the used quantities) are correct,
        # so the distance proxy is left as-is and documented. Use the Skyfield
        # backend if the True Node distance / FLG_XYZ is needed.
        pass  # keep dist proxy from LEB (lon/lat correct; distance documented)
    elif ipl == MEAN_APOG:
        # LEB stores old 5.145°·sin(ω) latitude model (max error ~20").
        # Override with 3-harmonic model fitted to the reference ephemeris output:
        #   lat = 5.1490449082·sin(ω) + 0.0034412113·sin(3ω)
        # where ω = (apogee_lon - node_lon).
        # Max residual vs the reference ephemeris: ~1.3", RMS: ~0.7".
        from .lunar import calc_mean_lunar_node

        node_lon = calc_mean_lunar_node(jd_tt)
        omega_rad = math.radians((lon - node_lon) % 360.0)
        lat = 5.1490449082 * math.sin(omega_rad) + 0.0034412113 * math.sin(
            3.0 * omega_rad
        )
        dist = _MOON_MEAN_APOG_DIST_AU

    # Nutation in longitude (dpsi) handling.
    # When SIDEREAL+EQUATORIAL: the reference ephemeris outputs mean ecliptic (no nutation)
    # Nutation handling for ecliptic-direct bodies.
    # Mean bodies are stored without nutation; true bodies include it.
    # FLG_NONUT: output on mean ecliptic (no nutation).
    # Velocity is NOT corrected — the Skyfield path also computes
    # velocity from the un-nutated polynomial.
    _sid_eq = bool(iflag & FLG_SIDEREAL) and bool(iflag & FLG_EQUATORIAL)
    _want_nonut_b = bool(iflag & FLG_NONUT)

    if ipl in (MEAN_NODE, MEAN_APOG):
        # Mean bodies stored without nutation. Add dpsi for true ecliptic,
        # UNLESS NONUT or sidereal+equatorial (both want mean ecliptic).
        if not _sid_eq and not _want_nonut_b:
            _, dpsi_rad, _, _ = _frame_data(jd_tt)
            lon = (lon + math.degrees(dpsi_rad)) % 360.0
    elif ipl in (TRUE_NODE, OSCU_APOG, INTP_APOG, INTP_PERG):
        # True/osculating bodies include nutation. Strip dpsi when NONUT
        # or sidereal+equatorial (both want mean ecliptic).
        if _sid_eq or _want_nonut_b:
            _, dpsi_rad, _, _ = _frame_data(jd_tt)
            lon = (lon - math.degrees(dpsi_rad)) % 360.0

    # Coordinate transforms for ecliptic-direct bodies.
    # Input coords are always ecliptic of date.
    #
    # FLG_J2000 is honored for ALL bodies, including TrueNode, OscuApog,
    # IntpApog, and IntpPerg.  the reference ephemeris silently ignores J2000 for these
    # four bodies when FLG_SIDEREAL is set — this is a behavioral bug
    # (ayanamsha and J2000 ecliptic precession are geometrically distinct
    # and composable operations).  LibEphemeris intentionally fixes this.
    # See docs/reference/se-bug-sidereal-j2000-nodes.md
    _effective_j2000 = bool(iflag & FLG_J2000)

    if (iflag & FLG_EQUATORIAL) and _effective_j2000:
        # J2000 equatorial: precess ecliptic-of-date -> J2000 ecliptic,
        # then rotate J2000 ecliptic -> J2000 equatorial.
        eps = OBLIQUITY_J2000_DEG

        def _ecl_date_to_eq_j2000(lo: float, la: float) -> tuple[float, float]:
            lo_j, la_j = _precess_ecliptic(lo, la, jd_tt, J2000)
            return _cotrans(lo_j, la_j, -eps)

        # Velocity via finite difference on original ecliptic coords
        dt_step = 0.001  # days
        eq_now_lon, eq_now_lat = _ecl_date_to_eq_j2000(lon, lat)
        eq_fwd_lon, eq_fwd_lat = _ecl_date_to_eq_j2000(
            lon + dlon * dt_step, lat + dlat * dt_step
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
        # J2000 ecliptic: precess from ecliptic of date to J2000
        # Velocity must also be transformed via finite difference
        dt_step = 0.001  # days
        j_now_lon, j_now_lat = _precess_ecliptic(lon, lat, jd_tt, J2000)
        j_fwd_lon, j_fwd_lat = _precess_ecliptic(
            lon + dlon * dt_step, lat + dlat * dt_step, jd_tt, J2000
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
    """Pipeline C: evaluate heliocentric ecliptic bodies (Uranians, Transpluto).

    LEB stores these as heliocentric J2000 ecliptic (lon, lat, dist).
    Default output is geocentric ecliptic of date, matching the Skyfield
    reference path in planets.py.

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
            # Apply light-time correction to match Skyfield's .observe():
            # the reference path uses sun.at(t).observe(earth) which evaluates
            # Earth at the retarded time (t - light_time_earth_sun).
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
    # The Skyfield reference path precesses J2000 → ecliptic of date ONLY when
    # the J2000 flag is NOT set.  For EQ+J2000, it converts J2000 ecliptic to
    # equatorial using obliquity of date (not J2000).  We must match this.
    is_j2000 = bool(iflag & FLG_J2000)
    is_equatorial = bool(iflag & FLG_EQUATORIAL)

    if is_equatorial and is_j2000:
        # EQ+J2000: Skyfield strips J2000 before _maybe_equatorial_convert,
        # so it uses true obliquity of date on J2000 ecliptic coords.
        # Sidereal mode uses mean obliquity (no nutation), matching the reference ephemeris.
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

    elif is_equatorial:
        # Equatorial of date: precess J2000 → date, then ecliptic → equatorial.
        # Sidereal mode uses mean obliquity (no nutation), matching the reference ephemeris.
        lon, lat = _precess_ecliptic(lon, lat, J2000, jd_tt)
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

    elif is_j2000:
        # J2000 ecliptic: already in J2000, no transform needed
        pass

    else:
        # Default: ecliptic of date
        lon, lat = _precess_ecliptic(lon, lat, J2000, jd_tt)

    return lon, lat, dist, dlon, dlat, ddist


# =============================================================================
# ENTRY POINTS
# =============================================================================


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
        KeyError: If body is not in the .leb file (caller should fall back).
        ValueError: If JD is outside the .leb file's range.
    """
    # FLG_ICRS: not yet implemented in LEB — fall back to Skyfield
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
        from .state import _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0

        sid_mode = _SIDEREAL_MODE
        sid_t0 = _SIDEREAL_T0
        sid_ayan_t0 = _SIDEREAL_AYAN_T0

    # Delta-T conversion: UT -> TT
    from .time_utils import deltat

    delta_t = deltat(tjd_ut)
    jd_tt = tjd_ut + delta_t

    if iflag & FLG_TOPOCTR:
        topo_offset = _topocentric_offset(topo_geopos, jd_tt, tjd_ut, reader)  # type: ignore[possibly-undefined]

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
    )


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
        from .state import _SIDEREAL_MODE, _SIDEREAL_T0, _SIDEREAL_AYAN_T0

        sid_mode = _SIDEREAL_MODE
        sid_t0 = _SIDEREAL_T0
        sid_ayan_t0 = _SIDEREAL_AYAN_T0

    tjd_ut = tjd_tt - reader.delta_t(tjd_tt)

    if iflag & FLG_TOPOCTR:
        topo_offset = _topocentric_offset(topo_geopos, tjd_tt, tjd_ut, reader)  # type: ignore[possibly-undefined]

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
    )


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
    topo_offset: Optional[Tuple[Tuple[float, float, float], Tuple[float, float, float]]] = None,
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

    # Check if body is in the .leb file
    if not reader.has_body(ipl):
        raise KeyError(f"Body {ipl} not in LEB file")

    body = reader._bodies[ipl]

    # Flag for deferred J2000 precession — set True by Pipeline B for any
    # ecliptic-direct body with SID+J2K (no EQ).  All Pipeline B bodies
    # (MeanNode, MeanApog, TrueNode, OscuApog, IntpApog, IntpPerg) use
    # the same deferred precession pattern: subtract ayanamsha first in
    # ecliptic-of-date coords, then precess to J2000.
    _deferred_sid_j2k = False

    _pipeline_a = False

    # Dispatch to appropriate pipeline based on coordinate type.
    # When XYZ+SIDEREAL, force Pipeline A to return spherical so the
    # shared sidereal correction (longitude subtraction) can run first;
    # XYZ conversion then happens in the post-processing stage.
    _xyz_sid = bool(iflag & FLG_XYZ) and bool(iflag & FLG_SIDEREAL)
    _pipe_iflag = (iflag & ~FLG_XYZ) if _xyz_sid else iflag

    if body.coord_type == COORD_ICRS_BARY:
        _pipeline_a = not _xyz_sid
        if iflag & FLG_SPEED:
            lon, lat, dist, dlon, dlat, ddist = _pipeline_icrs(
                reader, jd_tt, ipl, _pipe_iflag, want_velocity=True,
                topo_offset=topo_offset,
            )
        else:
            lon, lat, dist = _pipeline_icrs(
                reader, jd_tt, ipl, _pipe_iflag, topo_offset=topo_offset,
            )
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    elif body.coord_type == COORD_ICRS_BARY_SYSTEM:
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

    elif body.coord_type == COORD_ECLIPTIC:
        if iflag & FLG_TOPOCTR:
            raise KeyError("FLG_TOPOCTR not supported for ecliptic-direct bodies")
        # Pipeline B: ecliptic direct
        #
        # For all ecliptic-direct bodies with SID+J2K (no EQ), defer J2000
        # precession so ayanamsha can be subtracted first in ecliptic-of-date
        # coordinates.  The correct order is:
        #   ecl_date → −aya → precess_to_J2000
        # NOT:
        #   ecl_date → precess_to_J2000 → −aya  (wrong — non-commutative)
        #
        # This applies to ALL Pipeline B bodies uniformly: MeanNode, MeanApog,
        # TrueNode, OscuApog, IntpApog, IntpPerg.  the reference ephemeris only does this
        # for mean bodies (silently ignoring J2000 for the others when sidereal
        # is set) — LibEphemeris intentionally corrects this behavioral bug.
        # See docs/reference/se-bug-sidereal-j2000-nodes.md
        _deferred_sid_j2k = (
            bool(iflag & FLG_SIDEREAL)
            and bool(iflag & FLG_J2000)
            and not bool(iflag & FLG_EQUATORIAL)
        )
        _pipe_flags = (iflag & ~FLG_J2000) if _deferred_sid_j2k else iflag
        lon, lat, dist, dlon, dlat, ddist = _pipeline_ecliptic(
            reader, jd_tt, ipl, _pipe_flags
        )
        if not (iflag & FLG_SPEED):
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    elif body.coord_type == COORD_HELIO_ECL:
        if iflag & FLG_TOPOCTR:
            raise KeyError("FLG_TOPOCTR not supported for heliocentric bodies")
        # Pipeline C: heliocentric
        lon, lat, dist, dlon, dlat, ddist = _pipeline_helio(reader, jd_tt, ipl, iflag)
        if not (iflag & FLG_SPEED):
            dlon, dlat, ddist = 0.0, 0.0, 0.0

    else:
        raise ValueError(f"Unknown coord_type {body.coord_type}")

    # Sidereal correction
    #
    # When EQUATORIAL is set, the reference ephemeris does NOT subtract ayanamsha from RA
    # for any body type:
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
            # This applies uniformly to ALL bodies.  the reference ephemeris uses true
            # ayanamsha for TrueNode/OscuApog/IntpApog/IntpPerg even when
            # FLG_J2000 is set (because it skips J2000 precession for them).
            # LibEphemeris intentionally fixes this: when J2000 is requested,
            # mean ayanamsha is always used, for all bodies.
            # See docs/reference/se-bug-sidereal-j2000-nodes.md
            _eff_mean_aya = bool(iflag & (FLG_J2000 | FLG_NONUT))
            if _eff_mean_aya:
                aya = mean_aya
            else:
                try:
                    _, dpsi_sid, _, _ = _frame_data(jd_tt)
                    nutation_deg = math.degrees(dpsi_sid)
                    aya = mean_aya + nutation_deg
                except (KeyError, ValueError) as _nut_exc:
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
            # Fires for every sidereal SPEED request — including FLG_J2000, where
            # the drift is still present (the ayanamsa offset is constant in the
            # J2000 frame, but the of-date longitude this dlon was built from
            # still carries it). It is the FLG_SPEED-absent case that must be
            # excluded: dlon was zeroed above and must stay 0.0.
            # The deferred ecliptic-direct bodies (nodes/apogees) keep this
            # of-date drift; their _deferred_sid_j2k rebuild below re-precesses
            # the POSITION (not the rate) to J2000. Pipeline-A bodies get the
            # J2000 frame term applied separately, just after this block.
            if (iflag & FLG_SPEED) and (
                _deferred_sid_j2k or not (iflag & FLG_J2000) or _pipeline_a
            ):
                dlon -= _general_precession_rate_deg_day(jd_tt)

        except KeyError:
            # Star-based sidereal mode, fall back
            raise

    # J2000 longitude-speed frame conversion (Pipeline A, spherical output only).
    # _pipeline_icrs rotates the velocity into the J2000 ecliptic with a fixed
    # (instantaneous) matrix, which omits the time-derivative of that rotation —
    # the of-date→J2000 equinox motion. That motion removes the general-
    # precession rate from the longitude speed, so apply the subtraction here.
    # It runs for BOTH tropical and sidereal J2000 ecliptic output (for sidereal
    # it is the second ~p term, on top of the ayanamsa drift above). Excluded:
    # the deferred ecliptic-direct bodies (handled by the rebuild below;
    # _pipeline_a is False for them), equatorial output (its own frame), and XYZ
    # output (dlon then holds a Cartesian velocity component, not a longitude
    # rate — Pipeline-A XYZ returns the vectors straight from _pipeline_icrs).
    if (
        _pipeline_a
        and (iflag & FLG_SPEED)
        and (iflag & FLG_J2000)
        and not (iflag & FLG_EQUATORIAL)
        and not (iflag & FLG_XYZ)
    ):
        dlon -= _general_precession_rate_deg_day(jd_tt)

    # Deferred J2000 precession for Pipeline B bodies with SID+J2K.
    # The pipeline was run without FLG_J2000 so ayanamsha could be
    # subtracted in ecliptic-of-date first.  Now precess the sidereal-
    # corrected coords to J2000.  This applies to ALL ecliptic-direct
    # bodies uniformly (MeanNode, MeanApog, TrueNode, OscuApog, IntpApog,
    # IntpPerg).  See docs/reference/se-bug-sidereal-j2000-nodes.md
    if _deferred_sid_j2k:
        dt_step = 0.001  # days
        j_now_lon, j_now_lat = _precess_ecliptic(lon, lat, jd_tt, J2000)
        j_fwd_lon, j_fwd_lat = _precess_ecliptic(
            lon + dlon * dt_step, lat + dlat * dt_step, jd_tt, J2000
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

    # FLG_XYZ post-processing for Pipeline B/C (spherical → Cartesian)
    # Pipeline A handles XYZ internally via _pipeline_icrs; skip here.
    if (iflag & FLG_XYZ) and not _pipeline_a:
        return _spherical_to_cartesian_with_velocity(
            lon, lat, dist, dlon, dlat, ddist
        ), iflag

    # FLG_RADIANS: convert angular output from degrees to radians
    # Skip when XYZ is active — values are Cartesian AU, not angles
    if (iflag & FLG_RADIANS) and not (iflag & FLG_XYZ):
        _d2r = math.pi / 180.0
        lon = lon * _d2r
        lat = lat * _d2r
        dlon = dlon * _d2r
        dlat = dlat * _d2r

    return (lon, lat, dist, dlon, dlat, ddist), iflag
