# SPDX-License-Identifier: MIT
# Copyright (c) 2005 Johannes Gajdosik
# Third-party/adapted code — see file docstring and THIRD_PARTY_NOTICES.md
"""
TASS 1.7 - Theorie Analytique des Satellites de Saturne.

Provenance:
    This file is not project-owned code. It is the identified
    Johannes-Gajdosik Stellarium implementation retained under MIT, evaluating
    the Vienne-Duriez TASS theory. The upstream copyright and license remain at
    the top of the file and in ``THIRD_PARTY_NOTICES.md``; LibEphemeris only
    supplies the surrounding adapter boundary.

This is a Python port of the TASS 1.7 theory by Alain Vienne and Luc Duriez,
as implemented by Johannes Gajdosik for Stellarium (MIT license).

The theory provides positions for Saturn's 8 major satellites:
- 0: Mimas
- 1: Enceladus
- 2: Tethys
- 3: Dione
- 4: Rhea
- 5: Titan (most important for barycenter correction)
- 6: Iapetus
- 7: Hyperion

Note: In the original TASS 1.7, Hyperion was index 7 and Iapetus was index 6.
We follow the same convention.

Reference:
- Vienne, A. & Duriez, L. (1995) "TASS1.6", A&A 297, 588-605
- Original Fortran: ftp://ftp.imcce.fr/pub/ephem/satel/tass17/
- Stellarium implementation: Johannes Gajdosik (MIT license)

License: MIT (Stellarium implementation by Johannes Gajdosik)

Precision: ~50-100 km for all satellites (sufficient for ~0.05 arcsec at Saturn distance)
"""

import math
from typing import Tuple, List

from .tass17_data import (
    TASS17_BODIES,
    TASS17_TO_VSOP87,
    MIMAS_0,
    MIMAS_1,
    MIMAS_2,
    MIMAS_3,
    ENCELADUS_0,
    ENCELADUS_1,
    ENCELADUS_2,
    ENCELADUS_3,
    TETHYS_0,
    TETHYS_1,
    TETHYS_2,
    TETHYS_3,
    DIONE_0,
    DIONE_1,
    DIONE_2,
    DIONE_3,
    RHEA_0,
    RHEA_1,
    RHEA_2,
    RHEA_3,
    TITAN_0,
    TITAN_1,
    TITAN_2,
    TITAN_3,
    IAPETUS_0,
    IAPETUS_1,
    IAPETUS_2,
    IAPETUS_3,
    HYPERION_0,
    HYPERION_1,
    HYPERION_2,
    HYPERION_3,
)

# =============================================================================
# CONSTANTS
# =============================================================================

# TASS epoch: JD 2444240.0 = 1980-01-01 00:00 TT
TASS_EPOCH_JD: float = 2444240.0

# AU in km (IAU 2012)
AU_KM: float = 149597870.7

# Satellite names (matching body indices)
SATELLITE_NAMES = [
    "Mimas",
    "Enceladus",
    "Tethys",
    "Dione",
    "Rhea",
    "Titan",
    "Iapetus",
    "Hyperion",
]

# Body series lookup - maps body index to series arrays
_BODY_SERIES = [
    [MIMAS_0, MIMAS_1, MIMAS_2, MIMAS_3],
    [ENCELADUS_0, ENCELADUS_1, ENCELADUS_2, ENCELADUS_3],
    [TETHYS_0, TETHYS_1, TETHYS_2, TETHYS_3],
    [DIONE_0, DIONE_1, DIONE_2, DIONE_3],
    [RHEA_0, RHEA_1, RHEA_2, RHEA_3],
    [TITAN_0, TITAN_1, TITAN_2, TITAN_3],
    [IAPETUS_0, IAPETUS_1, IAPETUS_2, IAPETUS_3],
    [HYPERION_0, HYPERION_1, HYPERION_2, HYPERION_3],
]


def _calc_lon(t: float) -> List[float]:
    """
    Calculate mean longitudes for the first 7 satellites (excluding Hyperion).

    Args:
        t: Time since TASS epoch in days

    Returns:
        List of 7 mean longitudes (radians)
    """
    lon = [0.0] * 7

    for i in range(7):
        # Get the series[1] (longitude series) for this body
        series_1 = _BODY_SERIES[i][1]
        # First multiterm (indices all zero) contains the base longitude terms
        if len(series_1) > 0:
            indices, terms = series_1[0]
            for amplitude, phase, freq in terms:
                lon[i] += amplitude * math.sin(phase + freq * t)

    return lon


def _calc_tass17_elem(t: float, lon: List[float], body: int) -> List[float]:
    """
    Calculate the 6 orbital elements for a satellite.

    Args:
        t: Time since TASS epoch in days
        lon: Mean longitudes for all 7 satellites
        body: Body index (0-7)

    Returns:
        List of 6 orbital elements:
            [0] n: mean motion (rad/day)
            [1] lambda: mean longitude (rad)
            [2] e*cos(varpi): eccentricity vector x
            [3] e*sin(varpi): eccentricity vector y
            [4] sin(i/2)*cos(Omega): inclination vector x
            [5] sin(i/2)*sin(Omega): inclination vector y
    """
    # Get body parameters
    name, mu, aam, s0, series_names = TASS17_BODIES[body]

    # Initialize elements with constant terms (s0)
    elem = list(s0)

    # Series 0: perturbations to mean motion (n)
    series_0 = _BODY_SERIES[body][0]
    for indices, terms in series_0:
        arg = sum(indices[i] * lon[i] for i in range(7))
        for amplitude, phase, freq in terms:
            elem[0] += amplitude * math.cos(phase + freq * t + arg)
    elem[0] = aam * (1.0 + elem[0])

    # Series 1: perturbations to mean longitude (lambda)
    series_1 = _BODY_SERIES[body][1]
    if body != 7:  # Not Hyperion - first multiterm already used for lon
        start_idx = 1
        elem[1] += lon[body]
    else:
        start_idx = 0

    for multiterm_idx in range(start_idx, len(series_1)):
        indices, terms = series_1[multiterm_idx]
        arg = sum(indices[i] * lon[i] for i in range(7))
        for amplitude, phase, freq in terms:
            elem[1] += amplitude * math.sin(phase + freq * t + arg)
    elem[1] += aam * t

    # Series 2: perturbations to eccentricity vector (k = e*cos(varpi), h = e*sin(varpi))
    series_2 = _BODY_SERIES[body][2]
    for indices, terms in series_2:
        arg = sum(indices[i] * lon[i] for i in range(7))
        for amplitude, phase, freq in terms:
            x = phase + freq * t + arg
            elem[2] += amplitude * math.cos(x)
            elem[3] += amplitude * math.sin(x)

    # Series 3: perturbations to inclination vector (q = sin(i/2)*cos(Omega), p = sin(i/2)*sin(Omega))
    series_3 = _BODY_SERIES[body][3]
    for indices, terms in series_3:
        arg = sum(indices[i] * lon[i] for i in range(7))
        for amplitude, phase, freq in terms:
            x = phase + freq * t + arg
            elem[4] += amplitude * math.cos(x)
            elem[5] += amplitude * math.sin(x)

    return elem


def _elliptic_to_rectangular(
    mu: float, elem: List[float], dt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Convert orbital elements to rectangular coordinates.

    This implements the EllipticToRectangularN function from Stellarium.

    Args:
        mu: Gravitational parameter (Saturn GM / Saturn GM + satellite GM) in AU^3/day^2 units
        elem: 6 orbital elements [n, lambda, k, h, q, p]
        dt: Time since reference epoch (days)

    Returns:
        Tuple of (x, y, z, vx, vy, vz) in AU and AU/day
    """
    n = elem[0]  # Mean motion (rad/day)
    lam = elem[1]  # Mean longitude (rad)
    k = elem[2]  # e * cos(varpi)
    h = elem[3]  # e * sin(varpi)
    q = elem[4]  # sin(i/2) * cos(Omega)
    p = elem[5]  # sin(i/2) * sin(Omega)

    # Calculate semi-major axis from mean motion and mu
    # n = sqrt(GM / a^3), so a = (GM / n^2)^(1/3)
    # Here mu is already in appropriate units
    a_cubed = mu / (n * n)
    a = a_cubed ** (1.0 / 3.0)

    # Eccentricity
    e = math.sqrt(k * k + h * h)

    # Mean anomaly (adjusted for dt)
    M = lam + n * dt

    # Longitude of perihelion
    if e > 1e-12:
        varpi = math.atan2(h, k)
    else:
        varpi = 0.0

    # Mean anomaly relative to perihelion
    M_anom = M - varpi
    M_anom = M_anom % (2.0 * math.pi)

    # Solve Kepler's equation: E - e*sin(E) = M
    E = M_anom
    for _ in range(15):
        E_new = M_anom + e * math.sin(E)
        if abs(E_new - E) < 1e-14:
            break
        E = E_new

    cos_E = math.cos(E)
    sin_E = math.sin(E)

    # True anomaly
    beta = math.sqrt(1.0 - e * e)
    cos_nu = (cos_E - e) / (1.0 - e * cos_E)
    sin_nu = beta * sin_E / (1.0 - e * cos_E)

    # Distance
    r = a * (1.0 - e * cos_E)

    # Position in orbital plane (x towards perihelion, y perpendicular)
    x_orb = r * cos_nu
    y_orb = r * sin_nu

    # Velocity in orbital plane
    v_factor = math.sqrt(mu / a) / (1.0 - e * cos_E)
    vx_orb = -v_factor * sin_E
    vy_orb = v_factor * beta * cos_E

    # Inclination and node from q, p
    # q = sin(i/2) * cos(Omega), p = sin(i/2) * sin(Omega)
    sin_half_i_sq = q * q + p * p
    cos_half_i = math.sqrt(max(0.0, 1.0 - sin_half_i_sq))
    sin_half_i = math.sqrt(sin_half_i_sq)

    cos_i = 1.0 - 2.0 * sin_half_i_sq
    sin_i = 2.0 * sin_half_i * cos_half_i

    if sin_half_i > 1e-12:
        cos_Omega = q / sin_half_i
        sin_Omega = p / sin_half_i
    else:
        cos_Omega = 1.0
        sin_Omega = 0.0

    # Argument of perihelion: omega = varpi - Omega
    Omega = math.atan2(sin_Omega, cos_Omega)
    omega = varpi - Omega

    cos_omega = math.cos(omega)
    sin_omega = math.sin(omega)

    # Rotation from orbital plane to reference plane
    # First rotate by omega (argument of perihelion) in orbital plane
    x_rot = x_orb * cos_omega - y_orb * sin_omega
    y_rot = x_orb * sin_omega + y_orb * cos_omega

    vx_rot = vx_orb * cos_omega - vy_orb * sin_omega
    vy_rot = vx_orb * sin_omega + vy_orb * cos_omega

    # Then rotate by inclination around x-axis
    y_incl = y_rot * cos_i
    z_incl = y_rot * sin_i
    vy_incl = vy_rot * cos_i
    vz_incl = vy_rot * sin_i

    # Then rotate by Omega (longitude of node) around z-axis
    x_final = x_rot * cos_Omega - y_incl * sin_Omega
    y_final = x_rot * sin_Omega + y_incl * cos_Omega
    z_final = z_incl

    vx_final = vx_rot * cos_Omega - vy_incl * sin_Omega
    vy_final = vx_rot * sin_Omega + vy_incl * cos_Omega
    vz_final = vz_incl

    return (x_final, y_final, z_final, vx_final, vy_final, vz_final)


def _apply_rotation(x: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """
    Apply TASS17 to VSOP87 rotation matrix to transform coordinates.

    Args:
        x: Position vector in TASS17 reference frame

    Returns:
        Position vector in J2000 ecliptic (VSOP87) frame
    """
    R = TASS17_TO_VSOP87
    return (
        R[0] * x[0] + R[1] * x[1] + R[2] * x[2],
        R[3] * x[0] + R[4] * x[1] + R[5] * x[2],
        R[6] * x[0] + R[7] * x[1] + R[8] * x[2],
    )


def saturn_moon_position(jd: float, body: int = 5) -> Tuple[float, float, float]:
    """
    Calculate position of a Saturn satellite using TASS 1.7 theory.

    Args:
        jd: Julian Date (TT)
        body: Satellite index (0-7):
              0=Mimas, 1=Enceladus, 2=Tethys, 3=Dione,
              4=Rhea, 5=Titan, 6=Iapetus, 7=Hyperion

    Returns:
        Tuple (x, y, z) in AU, satellite position relative to Saturn center
        in the J2000 ecliptic frame (ICRF).
    """
    if body < 0 or body > 7:
        raise ValueError(f"Invalid body index: {body}. Must be 0-7.")

    # Time since TASS epoch in days
    t = jd - TASS_EPOCH_JD

    # Get body parameters
    name, mu, aam, s0, series_names = TASS17_BODIES[body]

    # Calculate mean longitudes for all satellites
    lon = _calc_lon(t)

    # Calculate orbital elements for this body
    elem = _calc_tass17_elem(t, lon, body)

    # Convert to rectangular coordinates (in Saturn-centric frame)
    xyz = _elliptic_to_rectangular(mu, elem, 0.0)

    # Apply rotation to J2000 ecliptic
    x_ecl, y_ecl, z_ecl = _apply_rotation((xyz[0], xyz[1], xyz[2]))

    return (x_ecl, y_ecl, z_ecl)


def saturn_moon_position_velocity(
    jd: float, body: int = 5
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate position and velocity of a Saturn satellite using TASS 1.7 theory.

    Args:
        jd: Julian Date (TT)
        body: Satellite index (0-7):
              0=Mimas, 1=Enceladus, 2=Tethys, 3=Dione,
              4=Rhea, 5=Titan, 6=Iapetus, 7=Hyperion

    Returns:
        Tuple (x, y, z, vx, vy, vz) in AU and AU/day, satellite position and velocity
        relative to Saturn center in the J2000 ecliptic frame (ICRF).
    """
    if body < 0 or body > 7:
        raise ValueError(f"Invalid body index: {body}. Must be 0-7.")

    # Time since TASS epoch in days
    t = jd - TASS_EPOCH_JD

    # Get body parameters
    name, mu, aam, s0, series_names = TASS17_BODIES[body]

    # Calculate mean longitudes for all satellites
    lon = _calc_lon(t)

    # Calculate orbital elements for this body
    elem = _calc_tass17_elem(t, lon, body)

    # Convert to rectangular coordinates (in Saturn-centric frame)
    xyz = _elliptic_to_rectangular(mu, elem, 0.0)

    # Apply rotation to J2000 ecliptic
    x_ecl, y_ecl, z_ecl = _apply_rotation((xyz[0], xyz[1], xyz[2]))
    vx_ecl, vy_ecl, vz_ecl = _apply_rotation((xyz[3], xyz[4], xyz[5]))

    return (x_ecl, y_ecl, z_ecl, vx_ecl, vy_ecl, vz_ecl)


# Backwards compatibility - simplified Titan function
def _titan_position(jd: float) -> Tuple[float, float, float]:
    """
    Calculate Titan's position using full TASS 1.7 theory.

    This is a wrapper for backwards compatibility.
    """
    return saturn_moon_position(jd, body=5)
