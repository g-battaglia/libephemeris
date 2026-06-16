# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
IAU 2006 Precession, IAU 2000B Nutation, and stellar aberration utilities.

Extracted from the former analytical computation package for use by the JPL/SPK calculation
pipeline (spk.py).  This module has NO dependencies on Skyfield or SPK files.
Only numpy and (optionally) pyerfa are used.

Functions provided:
- nutation_angles: IAU 2000B nutation (delta_psi, delta_epsilon)
- precess_from_j2000: Ecliptic precession from J2000 to date
- apply_aberration_to_position: Stellar aberration correction

References:
- Capitaine, Wallace, Chapront (2003), A&A 412, 567-586
- McCarthy & Luzum (2003), Celestial Mechanics 85, 37-49
- IAU SOFA software collection
"""

from __future__ import annotations

import math
from typing import Tuple, Any

import numpy as np

# NOTE: the ``precession_vondrak`` import is deliberately deferred to inside
# ``_mean_obliquity`` (not hoisted to module scope). This module keeps an
# import-time-erfa-optional contract -- the pyerfa import below is guarded, and
# ``test_cov100_astrometry`` execs this file as a standalone module (no package
# parent) to verify that guard. A top-level relative import would break that
# standalone load. ``precession_vondrak`` requires erfa, so importing it lazily
# preserves both the contract and the module's standalone-loadability.

# =============================================================================
# CONSTANTS
# =============================================================================

# Speed of light in AU/day
C_LIGHT_AU_DAY: float = 173.1446326846693

# Julian day constants
J2000: float = 2451545.0
JD_PER_CENTURY: float = 36525.0

# Angular conversion
TWO_PI: float = 2.0 * math.pi
DEG_TO_RAD: float = math.pi / 180.0
RAD_TO_DEG: float = 180.0 / math.pi
ARCSEC_TO_RAD: float = math.pi / (180.0 * 3600.0)

# Mean obliquity of the ecliptic at J2000.0 (arcseconds) – IAU 2006
OBLIQUITY_J2000_ARCSEC: float = 84381.406

# IAU 2006 polynomial for mean obliquity
OBLIQUITY_COEFFS: Tuple[float, ...] = (
    84381.406,
    -46.836769,
    -0.0001831,
    0.00200340,
    -0.000000576,
    -0.0000000434,
)

# =============================================================================
# HELPERS
# =============================================================================


def _jd_to_julian_centuries(jd_tt: float) -> float:
    """Julian centuries from J2000.0."""
    return (jd_tt - J2000) / JD_PER_CENTURY


def _normalize_angle(angle: float) -> float:
    """Normalize angle to [0, 360) degrees."""
    angle = angle % 360.0
    if angle < 0:  # pragma: no cover - modulo by positive 360.0 is always >= 0
        angle += 360.0
    return angle


def _cartesian_to_spherical(x: float, y: float, z: float) -> Tuple[float, float, float]:
    """Convert Cartesian to spherical (lon_deg, lat_deg, dist)."""
    dist = math.sqrt(x * x + y * y + z * z)
    if dist < 1e-30:
        return 0.0, 0.0, 0.0
    lon = math.atan2(y, x) * RAD_TO_DEG
    lat = math.asin(z / dist) * RAD_TO_DEG
    lon = _normalize_angle(lon)
    return lon, lat, dist


# =============================================================================
# PYERFA ACCELERATION (optional)
# =============================================================================

_HAS_ERFA = False
_erfa: Any = None

try:
    import erfa as _erfa_module

    _erfa = _erfa_module
    _HAS_ERFA = True
except ImportError:
    pass

# =============================================================================
# MEAN OBLIQUITY
# =============================================================================


def _julian_epoch(jd_tt: float) -> float:
    """Julian epoch (e.g. 2000.0 at J2000.0) for the ERFA long-term routines."""
    return 2000.0 + (jd_tt - J2000) / 365.25


def _mean_obliquity(jd_tt: float) -> float:
    """Of-date mean obliquity of the ecliptic in degrees.

    With erfa available (the normal runtime path — pyerfa is a required
    dependency) this is the Vondrák 2011 long-term obliquity (the angle between
    the Vondrák equator and ecliptic poles), long-term valid and consistent with
    the precession used everywhere else. The ``_HAS_ERFA`` guard keeps this
    function symmetric with the other reductions in this module: when erfa is
    monkeypatched off (the no-erfa unit tests) it falls back to the IAU 2006
    obliquity polynomial rather than crashing on the erfa-backed Vondrák path.
    """
    if _HAS_ERFA and _erfa is not None:
        # Deferred import: keeps this module loadable without erfa at import time
        # (see the module-scope NOTE and test_cov100_astrometry).
        from .precession_vondrak import vondrak_mean_obliquity_deg

        return vondrak_mean_obliquity_deg(jd_tt)

    # No-erfa fallback: IAU 2006 obliquity polynomial (valid near J2000).
    t = _jd_to_julian_centuries(jd_tt)
    eps_arcsec = sum(c * (t**i) for i, c in enumerate(OBLIQUITY_COEFFS))
    return eps_arcsec / 3600.0


# =============================================================================
# NUTATION (IAU 2000B – 77 terms)
# =============================================================================


def _fundamental_arguments(t: float) -> Tuple[float, float, float, float, float]:
    """Delaunay arguments for nutation (el, elp, F, D, Omega) in radians."""
    el = 485868.249036 + t * (
        1717915923.2178 + t * (31.8792 + t * (0.051635 + t * (-0.00024470)))
    )
    el = (el * ARCSEC_TO_RAD) % TWO_PI

    elp = 1287104.79305 + t * (
        129596581.0481 + t * (-0.5532 + t * (0.000136 + t * (-0.00001149)))
    )
    elp = (elp * ARCSEC_TO_RAD) % TWO_PI

    F = 335779.526232 + t * (
        1739527262.8478 + t * (-12.7512 + t * (-0.001037 + t * (0.00000417)))
    )
    F = (F * ARCSEC_TO_RAD) % TWO_PI

    D = 1072260.70369 + t * (
        1602961601.2090 + t * (-6.3706 + t * (0.006593 + t * (-0.00003169)))
    )
    D = (D * ARCSEC_TO_RAD) % TWO_PI

    Omega = 450160.398036 + t * (
        -6962890.5431 + t * (7.4722 + t * (0.007702 + t * (-0.00005939)))
    )
    Omega = (Omega * ARCSEC_TO_RAD) % TWO_PI

    return el, elp, F, D, Omega


# IAU 2000B nutation series – 77 terms (0.1 microarcseconds)
_NUTATION_TERMS_IAU2000B: Tuple[Tuple[int, ...], ...] = (
    (0, 0, 0, 0, 1, -172064161, -174666, 92052331, 9086),
    (0, 0, 2, -2, 2, -13170906, -1675, 5730336, -3015),
    (0, 0, 2, 0, 2, -2276413, -234, 978459, -485),
    (0, 0, 0, 0, 2, 2074554, 207, -897492, 470),
    (0, 1, 0, 0, 0, 1475877, -3633, 73871, -184),
    (0, 1, 2, -2, 2, -516821, 1226, 224386, -677),
    (1, 0, 0, 0, 0, 711159, 73, -6750, 0),
    (0, 0, 2, 0, 1, -387298, -367, 200728, 18),
    (1, 0, 2, 0, 2, -301461, -36, 129025, -63),
    (0, -1, 2, -2, 2, 215829, -494, -95929, 299),
    (0, 0, 2, -2, 1, 128227, 137, -68982, -9),
    (-1, 0, 2, 0, 2, 123457, 11, -53311, 32),
    (-1, 0, 0, 2, 0, 156994, 10, -1235, 0),
    (1, 0, 0, 0, 1, 63110, 63, -33228, 0),
    (-1, 0, 0, 0, 1, -57976, -63, 31429, 0),
    (-1, 0, 2, 2, 2, -59641, -11, 25543, -11),
    (1, 0, 2, 0, 1, -51613, -42, 26366, 0),
    (-2, 0, 2, 0, 1, 45893, 50, -24236, -10),
    (0, 0, 0, 2, 0, 63384, 11, -1220, 0),
    (0, 0, 2, 2, 2, -38571, -1, 16452, -11),
    (0, -1, 2, 0, 2, 32481, 0, -13870, 0),
    (2, 0, 0, -2, 0, -47722, 0, 477, 0),
    (2, 0, 2, 0, 2, -31046, -1, 13238, -11),
    (1, 0, 2, -2, 2, 28593, 0, -12338, 10),
    (-1, 0, 2, 0, 1, 20441, 21, -10758, 0),
    (2, 0, 0, 0, 0, 29243, 0, -609, 0),
    (0, 0, 2, 0, 0, 25887, 0, -550, 0),
    (-1, 0, 0, 2, 1, -14053, -25, 8551, -2),
    (0, 2, 0, 0, 0, 15164, 10, -167, 0),
    (0, 2, 2, -2, 2, -15794, 72, 6850, -42),
    (-1, 0, 0, 2, -1, 21783, 0, -167, 0),
    (0, 1, 0, 0, 1, -12873, -10, 6953, 0),
    (1, 0, 0, -2, 1, -12654, 11, 6415, 0),
    (0, -1, 0, 0, 1, -10204, 0, 5222, 0),
    (0, 0, 2, -2, 0, 16707, -85, 168, -1),
    (2, 0, 2, -2, 2, -7691, 0, 3268, 0),
    (1, 0, 0, 2, 0, -11024, 0, 104, 0),
    (1, 0, 2, -2, 1, 7566, -21, -4100, 0),
    (0, 0, 0, 2, 1, -6637, -11, 3614, 0),
    (-1, 0, 2, 2, 1, -7141, 21, 2991, 0),
    (0, 2, 0, 0, -1, -6302, -11, 3334, 0),
    (1, 0, 0, -2, -1, 5800, 10, -3040, 0),
    (0, -1, 2, 0, 1, 6443, 0, -2768, 0),
    (-1, 0, 2, 2, 2, -5774, -11, 2499, 0),
    (1, 1, 0, -2, 0, -5350, 0, 111, 0),
    (-2, 0, 2, 0, 0, -4752, -11, -3, 0),
    (0, 1, 2, 0, 2, -4940, -11, 2107, 0),
    (0, -1, 2, 2, 2, -3987, 0, 1681, 0),
    (-1, 0, 0, 0, 2, -3673, 0, 1600, 0),
    (1, 1, 0, 0, 0, -2964, 0, -78, 0),
    (0, 1, 2, -2, 1, 4183, 0, -2056, 0),
    (-1, 0, 2, 0, 0, 3596, 0, -77, 0),
    (0, -1, 0, 2, 0, 4155, 0, -41, 0),
    (0, 0, 0, 1, 0, 3319, 0, -33, 0),
    (-1, 1, 0, 0, 0, 2885, 0, -29, 0),
    (-1, 0, 0, -1, 0, -2904, 0, 29, 0),
    (0, 0, 2, 1, 2, -2972, 0, 1261, 0),
    (1, 0, 0, 0, 2, -2827, 0, 1221, 0),
    (2, 0, 2, 0, 1, -2732, 0, 1167, 0),
    (-1, 1, 0, 2, 0, -2489, 0, 25, 0),
    (0, 0, 2, -1, 2, -2412, 0, 1022, 0),
    (2, 0, 0, 0, 1, -2499, 0, 1064, 0),
    (0, 0, 0, 0, 3, 2571, 0, -1103, 0),
    (-1, -1, 0, 2, 0, -2344, 0, 24, 0),
    (1, -1, 0, 0, 0, 2305, 0, -23, 0),
    (0, 0, 2, 2, 1, -2215, 0, 942, 0),
    (1, 0, 2, 2, 2, 2043, 0, -864, 0),
    (-2, 0, 0, 2, 1, -1861, 0, 799, 0),
    (0, -1, 2, -2, 1, 2204, 0, -945, 0),
    (2, 0, 2, 0, 0, -1955, 0, 42, 0),
    (-1, 1, 2, 0, 2, 1788, 0, -760, 0),
    (0, 1, 0, -2, 0, 1710, 0, -17, 0),
    (-1, -1, 2, 2, 2, 1645, 0, -693, 0),
    (0, -1, 0, 0, 2, 1620, 0, -697, 0),
    (0, 1, -2, 0, 0, 1541, 0, -15, 0),
    (1, 0, 2, -2, 0, -1487, 0, 15, 0),
)


def _nutation_angles_numpy(jd_tt: float) -> Tuple[float, float]:
    """IAU 2000B nutation in longitude/obliquity (pure numpy)."""
    t = _jd_to_julian_centuries(jd_tt)
    el, elp, F, D, Omega = _fundamental_arguments(t)

    dpsi = 0.0
    deps = 0.0

    for term in _NUTATION_TERMS_IAU2000B:
        arg = term[0] * el + term[1] * elp + term[2] * F + term[3] * D + term[4] * Omega
        sin_arg = math.sin(arg)
        cos_arg = math.cos(arg)
        dpsi += (term[5] + term[6] * t) * sin_arg
        deps += (term[7] + term[8] * t) * cos_arg

    dpsi_deg = dpsi * 0.1e-6 * ARCSEC_TO_RAD * RAD_TO_DEG
    deps_deg = deps * 0.1e-6 * ARCSEC_TO_RAD * RAD_TO_DEG
    return dpsi_deg, deps_deg


def nutation_angles(jd_tt: float) -> Tuple[float, float]:
    """Calculate nutation in longitude and obliquity (IAU 2000B).

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        Tuple of (delta_psi, delta_epsilon) in degrees.
    """
    if _HAS_ERFA and _erfa is not None:
        dpsi, deps = _erfa.nut00b(J2000, jd_tt - J2000)
        return float(dpsi) * RAD_TO_DEG, float(deps) * RAD_TO_DEG
    return _nutation_angles_numpy(jd_tt)


# =============================================================================
# PRECESSION
# =============================================================================

# Equatorial precession angles from J2000 (Lieske 1977), polynomial in
# Julian centuries TT since J2000: c0 + c1*t + c2*t^2 + c3*t^3 (arcsec)
ZETA_A_COEFFS: Tuple[float, ...] = (0.0, 2306.2181, 0.30188, 0.017998)
Z_A_COEFFS: Tuple[float, ...] = (0.0, 2306.2181, 1.09468, 0.018203)
THETA_A_COEFFS: Tuple[float, ...] = (0.0, 2004.3109, -0.42665, -0.041833)


def _precession_matrix_j2000_to_date(jd_tt: float) -> np.ndarray:
    """Precession rotation matrix from J2000 to date (Vondrák 2011 via erfa).

    Uses the long-term ``erfa.ltp`` (J2000 mean equator/equinox -> date), valid
    over ±200,000 years. The IAU 2006 / Lieske polynomial path below is a
    fallback only when pyerfa is unavailable.
    """
    if _HAS_ERFA and _erfa is not None:
        return np.array(_erfa.ltp(_julian_epoch(jd_tt)))

    t = _jd_to_julian_centuries(jd_tt)

    zeta = sum(c * (t**i) for i, c in enumerate(ZETA_A_COEFFS))
    z = sum(c * (t**i) for i, c in enumerate(Z_A_COEFFS))
    theta = sum(c * (t**i) for i, c in enumerate(THETA_A_COEFFS))

    zeta_rad = zeta * ARCSEC_TO_RAD
    z_rad = z * ARCSEC_TO_RAD
    theta_rad = theta * ARCSEC_TO_RAD

    cos_zeta = math.cos(zeta_rad)
    sin_zeta = math.sin(zeta_rad)
    cos_z = math.cos(z_rad)
    sin_z = math.sin(z_rad)
    cos_theta = math.cos(theta_rad)
    sin_theta = math.sin(theta_rad)

    r11 = cos_zeta * cos_theta * cos_z - sin_zeta * sin_z
    r12 = -sin_zeta * cos_theta * cos_z - cos_zeta * sin_z
    r13 = -sin_theta * cos_z
    r21 = cos_zeta * cos_theta * sin_z + sin_zeta * cos_z
    r22 = -sin_zeta * cos_theta * sin_z + cos_zeta * cos_z
    r23 = -sin_theta * sin_z
    r31 = cos_zeta * sin_theta
    r32 = -sin_zeta * sin_theta
    r33 = cos_theta

    return np.array([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])


def _nutation_matrix(jd_tt: float) -> np.ndarray:
    """Nutation rotation matrix."""
    if _HAS_ERFA and _erfa is not None:
        dpsi, deps = _erfa.nut00b(J2000, jd_tt - J2000)
        eps = _erfa.obl06(J2000, jd_tt - J2000)
        return np.array(_erfa.numat(eps, dpsi, deps))

    dpsi_deg, deps_deg = nutation_angles(jd_tt)
    eps_mean_deg = _mean_obliquity(jd_tt)

    dpsi = dpsi_deg * DEG_TO_RAD
    deps = deps_deg * DEG_TO_RAD
    eps = eps_mean_deg * DEG_TO_RAD
    eps_true = eps + deps

    cos_eps = math.cos(eps)
    sin_eps = math.sin(eps)
    cos_eps_true = math.cos(eps_true)
    sin_eps_true = math.sin(eps_true)
    cos_dpsi = math.cos(dpsi)
    sin_dpsi = math.sin(dpsi)

    r11 = cos_dpsi
    r12 = -sin_dpsi * cos_eps
    r13 = -sin_dpsi * sin_eps
    r21 = sin_dpsi * cos_eps_true
    r22 = cos_dpsi * cos_eps_true * cos_eps + sin_eps_true * sin_eps
    r23 = cos_dpsi * cos_eps_true * sin_eps - sin_eps_true * cos_eps
    r31 = sin_dpsi * sin_eps_true
    r32 = cos_dpsi * sin_eps_true * cos_eps - cos_eps_true * sin_eps
    r33 = cos_dpsi * sin_eps_true * sin_eps + cos_eps_true * cos_eps

    return np.array([[r11, r12, r13], [r21, r22, r23], [r31, r32, r33]])


def _precession_nutation_matrix(jd_tt: float) -> np.ndarray:
    """Combined precession-nutation matrix from J2000 to true equator of date."""
    if _HAS_ERFA and _erfa is not None:
        return np.array(_erfa.pnm06a(J2000, jd_tt - J2000))
    P = _precession_matrix_j2000_to_date(jd_tt)
    N = _nutation_matrix(jd_tt)
    return N @ P


def _true_obliquity(jd_tt: float) -> float:
    """True obliquity = mean obliquity + nutation in obliquity (degrees)."""
    eps_mean = _mean_obliquity(jd_tt)
    _, d_eps = nutation_angles(jd_tt)
    return eps_mean + d_eps


def _precess_ecliptic(
    lon: float, lat: float, from_jd: float, to_jd: float
) -> Tuple[float, float]:
    """Precess ecliptic coordinates from one epoch to another."""
    if abs(from_jd - to_jd) < 1e-6:
        return lon, lat

    lon_rad = lon * DEG_TO_RAD
    lat_rad = lat * DEG_TO_RAD

    eps1 = _mean_obliquity(from_jd) * DEG_TO_RAD
    eps2 = _mean_obliquity(to_jd) * DEG_TO_RAD

    cos_lat = math.cos(lat_rad)
    sin_lat = math.sin(lat_rad)
    cos_lon = math.cos(lon_rad)
    sin_lon = math.sin(lon_rad)
    cos_eps1 = math.cos(eps1)
    sin_eps1 = math.sin(eps1)

    x1 = cos_lat * cos_lon
    y1 = cos_lat * sin_lon * cos_eps1 - sin_lat * sin_eps1
    z1 = cos_lat * sin_lon * sin_eps1 + sin_lat * cos_eps1

    # Equatorial precession from_jd -> to_jd via Vondrák 2011 (long-term valid):
    # erfa.ltp(to) @ erfa.ltp(from)^T composes the two J2000-referred matrices.
    # The Lieske (1977) two-epoch polynomial is a fallback when pyerfa is absent.
    if _HAS_ERFA and _erfa is not None:
        prec = _erfa.ltp(_julian_epoch(to_jd)) @ _erfa.ltp(_julian_epoch(from_jd)).T
    else:
        t1 = _jd_to_julian_centuries(from_jd)
        dt = _jd_to_julian_centuries(to_jd) - t1
        zeta = (
            (2306.2181 + 1.39656 * t1 - 0.000139 * t1**2) * dt
            + (0.30188 - 0.000344 * t1) * dt**2
            + 0.017998 * dt**3
        ) * ARCSEC_TO_RAD
        z = (
            (2306.2181 + 1.39656 * t1 - 0.000139 * t1**2) * dt
            + (1.09468 + 0.000066 * t1) * dt**2
            + 0.018203 * dt**3
        ) * ARCSEC_TO_RAD
        theta = (
            (2004.3109 - 0.85330 * t1 - 0.000217 * t1**2) * dt
            - (0.42665 + 0.000217 * t1) * dt**2
            - 0.041833 * dt**3
        ) * ARCSEC_TO_RAD
        cz, sz = math.cos(zeta), math.sin(zeta)
        cZ, sZ = math.cos(z), math.sin(z)
        ct, st = math.cos(theta), math.sin(theta)
        prec = np.array(
            [
                [cz * ct * cZ - sz * sZ, -sz * ct * cZ - cz * sZ, -st * cZ],
                [cz * ct * sZ + sz * cZ, -sz * ct * sZ + cz * cZ, -st * sZ],
                [cz * st, -sz * st, ct],
            ]
        )

    x2 = prec[0, 0] * x1 + prec[0, 1] * y1 + prec[0, 2] * z1
    y2 = prec[1, 0] * x1 + prec[1, 1] * y1 + prec[1, 2] * z1
    z2 = prec[2, 0] * x1 + prec[2, 1] * y1 + prec[2, 2] * z1

    cos_eps2 = math.cos(eps2)
    sin_eps2 = math.sin(eps2)

    y_ecl = y2 * cos_eps2 + z2 * sin_eps2
    z_ecl = -y2 * sin_eps2 + z2 * cos_eps2

    new_lon = math.atan2(y_ecl, x2) * RAD_TO_DEG
    new_lat = math.asin(max(-1.0, min(1.0, z_ecl))) * RAD_TO_DEG

    return _normalize_angle(new_lon), new_lat


def precess_from_j2000(lon: float, lat: float, jd_tt: float) -> Tuple[float, float]:
    """Precess ecliptic coordinates from J2000.0 to date.

    Args:
        lon: Ecliptic longitude in degrees (J2000.0).
        lat: Ecliptic latitude in degrees (J2000.0).
        jd_tt: Julian Day (TT) of the target epoch.

    Returns:
        Tuple of (lon_date, lat_date) in degrees.
    """
    return _precess_ecliptic(lon, lat, J2000, jd_tt)


# =============================================================================
# ABERRATION
# =============================================================================


def _annual_aberration_cartesian(
    target_direction: Tuple[float, float, float],
    earth_velocity: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """Stellar aberration correction (first-order Bradley formula, Cartesian)."""
    rx, ry, rz = target_direction
    r_mag = math.sqrt(rx * rx + ry * ry + rz * rz)
    if r_mag < 1e-30:
        return (0.0, 0.0, 0.0)
    rx /= r_mag
    ry /= r_mag
    rz /= r_mag

    vx = earth_velocity[0] / C_LIGHT_AU_DAY
    vy = earth_velocity[1] / C_LIGHT_AU_DAY
    vz = earth_velocity[2] / C_LIGHT_AU_DAY

    r_dot_v = rx * vx + ry * vy + rz * vz

    dx = vx - r_dot_v * rx
    dy = vy - r_dot_v * ry
    dz = vz - r_dot_v * rz

    return (dx, dy, dz)


def apply_aberration_to_position(
    position: Tuple[float, float, float],
    earth_velocity: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """Apply stellar aberration correction to a geocentric position.

    Args:
        position: Geocentric position (x, y, z) in AU.
        earth_velocity: Earth heliocentric velocity (vx, vy, vz) in AU/day.

    Returns:
        Apparent (aberrated) position (x, y, z) in AU.
    """
    x, y, z = position
    distance = math.sqrt(x * x + y * y + z * z)
    if distance < 1e-30:
        return position

    ux, uy, uz = x / distance, y / distance, z / distance
    dx, dy, dz = _annual_aberration_cartesian((ux, uy, uz), earth_velocity)

    ax = ux + dx
    ay = uy + dy
    az = uz + dz

    a_mag = math.sqrt(ax * ax + ay * ay + az * az)
    if a_mag < 1e-30:
        return position
    ax /= a_mag
    ay /= a_mag
    az /= a_mag

    return (ax * distance, ay * distance, az * distance)
