"""
Uranian moon orbital theory for Uranus barycenter correction.

Uranus has 5 major moons that dominate its satellite system mass:
- Titania (largest, ~45% of moon mass)
- Oberon (~32%)
- Ariel (~15%)
- Umbriel (~14%)
- Miranda (~1%)

The barycenter offset is computed from the weighted sum of moon positions.

References:
- Laskar, J. (1987) "GUST86 - an analytical ephemeris of the Uranian satellites"
  A&A 188, 212-218
- JPL SSD: Planetary Satellite Mean Elements
  https://ssd.jpl.nasa.gov/sats/elem
- Jacobson, R.A. (2014) "The orbits of the Uranian satellites and rings"
  AJ 148, 76

Precision: mean circular elements without the GUST86 periodic terms —
moon positions are good to a few thousand km (phases and rounded nodes
dominate).  Weighted by the ~1e-5 moon/planet mass ratios, the resulting
COB offset error is tens of km, i.e. sub-milliarcsecond at Uranus
opposition distance, far below the planet_centers SPK alternative.
"""

from __future__ import annotations

import math
from typing import Tuple

# =============================================================================
# PHYSICAL CONSTANTS
# =============================================================================

# GM values (km³/s²) from URA111 / JPL
GM_URANUS: float = 5793951.3  # URA111

# Major moon GM values (approximate, from literature)
GM_ARIEL: float = 14.0
GM_UMBRIEL: float = 13.0
GM_TITANIA: float = 41.0
GM_OBERON: float = 29.0
GM_MIRANDA: float = 0.7

# Total moon GM
GM_URANUS_MOONS_TOTAL: float = (
    GM_ARIEL + GM_UMBRIEL + GM_TITANIA + GM_OBERON + GM_MIRANDA
)

# AU in km
AU_KM: float = 149597870.7

# J2000 epoch
EPOCH_JD: float = 2451545.0

# =============================================================================
# ORBITAL ELEMENTS (J2000, from JPL SSD / GUST86)
# =============================================================================
# All elements are with respect to Uranus equator, then transformed to ICRF

# Semi-major axes (km)
ARIEL_A_KM: float = 190945.0
UMBRIEL_A_KM: float = 266035.0
TITANIA_A_KM: float = 436300.0
OBERON_A_KM: float = 583520.0
MIRANDA_A_KM: float = 129390.0

# Orbital periods (days)
ARIEL_PERIOD: float = 2.520379
UMBRIEL_PERIOD: float = 4.144177
TITANIA_PERIOD: float = 8.705872
OBERON_PERIOD: float = 13.463238
MIRANDA_PERIOD: float = 1.413479

# Mean motions (rad/day)
ARIEL_N: float = 2.0 * math.pi / ARIEL_PERIOD
UMBRIEL_N: float = 2.0 * math.pi / UMBRIEL_PERIOD
TITANIA_N: float = 2.0 * math.pi / TITANIA_PERIOD
OBERON_N: float = 2.0 * math.pi / OBERON_PERIOD
MIRANDA_N: float = 2.0 * math.pi / MIRANDA_PERIOD

# Eccentricities (very small for major moons)
ARIEL_E: float = 0.0012
UMBRIEL_E: float = 0.0040
TITANIA_E: float = 0.0014
OBERON_E: float = 0.0016
MIRANDA_E: float = 0.0013

# Mean anomaly at J2000 (radians, approximate)
ARIEL_M0: float = math.radians(123.0)
UMBRIEL_M0: float = math.radians(285.0)
TITANIA_M0: float = math.radians(107.0)
OBERON_M0: float = math.radians(222.0)
MIRANDA_M0: float = math.radians(33.0)

# Inclinations to Uranus equator (radians) - very small
ARIEL_I: float = math.radians(0.04)
UMBRIEL_I: float = math.radians(0.13)
TITANIA_I: float = math.radians(0.08)
OBERON_I: float = math.radians(0.07)
MIRANDA_I: float = math.radians(4.34)

# Longitudes of ascending node on Uranus equator (radians)
# All moons orbit in Uranus equatorial plane with small nodes
ARIEL_NODE: float = math.radians(167.0)
UMBRIEL_NODE: float = math.radians(197.0)
TITANIA_NODE: float = math.radians(0.0)
OBERON_NODE: float = math.radians(0.0)
MIRANDA_NODE: float = math.radians(97.0)

# Arguments of pericenter (radians)
ARIEL_OMEGA: float = math.radians(92.0)
UMBRIEL_OMEGA: float = math.radians(144.0)
TITANIA_OMEGA: float = math.radians(100.0)
OBERON_OMEGA: float = math.radians(282.0)
MIRANDA_OMEGA: float = math.radians(343.0)

# =============================================================================
# URANUS POLE ORIENTATION (ICRF)
# =============================================================================
# Uranus' north pole direction in ICRF/J2000
# From IAU 2006/2015 report

URANUS_POLE_RA: float = math.radians(257.311)  # Right ascension of pole
URANUS_POLE_DEC: float = math.radians(-15.175)  # Declination of pole

# Pre-matrix: rotation from Uranus equatorial frame to ICRF
# This is computed once for efficiency


def _compute_uranus_to_icrf_matrix() -> Tuple[Tuple[float, ...], ...]:
    """Compute rotation matrix from Uranus equatorial frame to ICRF.

    The Uranus-equator frame is the IAU non-rotating body frame:
    +z along the IAU pole (alpha0, delta0), +x along the ascending node
    of the body equator on the ICRF equator (the IAU node Q, at right
    ascension alpha0 + 90deg).  The moon node longitudes above are
    measured from that origin.

    The matrix columns are the frame's basis vectors expressed in ICRF,
    so for a vector u in the body-equator frame, v_icrf = M @ u.
    """
    cos_ra = math.cos(URANUS_POLE_RA)
    sin_ra = math.sin(URANUS_POLE_RA)
    cos_dec = math.cos(URANUS_POLE_DEC)
    sin_dec = math.sin(URANUS_POLE_DEC)

    # x-axis: IAU node direction (RA = alpha0 + 90deg, Dec = 0)
    x_axis = (-sin_ra, cos_ra, 0.0)
    # z-axis: pole direction
    z_axis = (cos_dec * cos_ra, cos_dec * sin_ra, sin_dec)
    # y-axis: z cross x (right-handed)
    y_axis = (-sin_dec * cos_ra, -sin_dec * sin_ra, cos_dec)

    return (
        (x_axis[0], y_axis[0], z_axis[0]),
        (x_axis[1], y_axis[1], z_axis[1]),
        (x_axis[2], y_axis[2], z_axis[2]),
    )


# Precompute rotation matrix
_URANUS_TO_ICRF = _compute_uranus_to_icrf_matrix()


# =============================================================================
# KEPLERIAN POSITION CALCULATION
# =============================================================================


def _solve_kepler(M: float, e: float, tolerance: float = 1e-12) -> float:
    """Solve Kepler's equation M = E - e*sin(E) for E.

    Uses Newton-Raphson iteration.

    Args:
        M: Mean anomaly (radians)
        e: Eccentricity
        tolerance: Convergence tolerance

    Returns:
        Eccentric anomaly E (radians)
    """
    E = M  # Initial guess
    for _ in range(20):
        delta = (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))
        E -= delta
        if abs(delta) < tolerance:
            break
    return E


def _moon_position_kepler(
    jd: float,
    a_km: float,
    e: float,
    n: float,
    M0: float,
    omega: float,
    i: float,
    node: float,
) -> Tuple[float, float, float]:
    """Calculate moon position in ICRF coordinates using Keplerian elements.

    Args:
        jd: Julian Date (TT)
        a_km: Semi-major axis in km
        e: Eccentricity
        n: Mean motion in rad/day
        M0: Mean anomaly at J2000 (radians)
        omega: Argument of pericenter (radians)
        i: Inclination (radians)
        node: Longitude of ascending node (radians)

    Returns:
        Tuple (x, y, z) position in km relative to Uranus center, in ICRF frame
    """
    # Time since J2000
    dt = jd - EPOCH_JD

    # Mean anomaly at time t
    M = M0 + n * dt
    M = M % (2.0 * math.pi)

    # Solve Kepler's equation
    E = _solve_kepler(M, e)

    # True anomaly
    cos_nu = (math.cos(E) - e) / (1.0 - e * math.cos(E))
    sin_nu = math.sqrt(1.0 - e * e) * math.sin(E) / (1.0 - e * math.cos(E))
    nu = math.atan2(sin_nu, cos_nu)

    # Radius
    r = a_km * (1.0 - e * e) / (1.0 + e * math.cos(nu))

    # Position in orbital plane (perifocal frame)
    x_peri = r * math.cos(nu)
    y_peri = r * math.sin(nu)
    z_peri = 0.0

    # Rotation matrices
    cos_omega = math.cos(omega)
    sin_omega = math.sin(omega)
    cos_i = math.cos(i)
    sin_i = math.sin(i)
    cos_node = math.cos(node)
    sin_node = math.sin(node)

    # Rotate by argument of pericenter (in orbital plane)
    x1 = x_peri * cos_omega - y_peri * sin_omega
    y1 = x_peri * sin_omega + y_peri * cos_omega
    z1 = 0.0

    # Rotate by inclination (around node axis)
    x2 = x1
    y2 = y1 * cos_i - z1 * sin_i
    z2 = y1 * sin_i + z1 * cos_i

    # Rotate by node (around z-axis)
    x_eq = x2 * cos_node - y2 * sin_node
    y_eq = x2 * sin_node + y2 * cos_node
    z_eq = z2

    # Now rotate from Uranus equatorial frame to ICRF
    # Uranus' equatorial plane is tilted ~98 degrees from its orbital plane
    # and its pole points in a specific direction in ICRF

    # Simplified rotation: Uranus equator -> ICRF
    # Uranus has an axial tilt of ~97.77 degrees
    # The equatorial plane is nearly perpendicular to the ecliptic

    # Apply rotation matrix
    m = _URANUS_TO_ICRF
    x_icrf = m[0][0] * x_eq + m[0][1] * y_eq + m[0][2] * z_eq
    y_icrf = m[1][0] * x_eq + m[1][1] * y_eq + m[1][2] * z_eq
    z_icrf = m[2][0] * x_eq + m[2][1] * y_eq + m[2][2] * z_eq

    return (x_icrf, y_icrf, z_icrf)


# =============================================================================
# PUBLIC API
# =============================================================================


def uranian_moon_positions(
    jd: float,
) -> Tuple[
    Tuple[float, float, float],  # Ariel
    Tuple[float, float, float],  # Umbriel
    Tuple[float, float, float],  # Titania
    Tuple[float, float, float],  # Oberon
    Tuple[float, float, float],  # Miranda
]:
    """Calculate positions of the 5 major Uranian moons.

    Args:
        jd: Julian Date (TT)

    Returns:
        Tuple of 5 (x, y, z) positions in km relative to Uranus center,
        in ICRF/J2000 frame. Order: Ariel, Umbriel, Titania, Oberon, Miranda.
    """
    ariel = _moon_position_kepler(
        jd, ARIEL_A_KM, ARIEL_E, ARIEL_N, ARIEL_M0, ARIEL_OMEGA, ARIEL_I, ARIEL_NODE
    )
    umbriel = _moon_position_kepler(
        jd,
        UMBRIEL_A_KM,
        UMBRIEL_E,
        UMBRIEL_N,
        UMBRIEL_M0,
        UMBRIEL_OMEGA,
        UMBRIEL_I,
        UMBRIEL_NODE,
    )
    titania = _moon_position_kepler(
        jd,
        TITANIA_A_KM,
        TITANIA_E,
        TITANIA_N,
        TITANIA_M0,
        TITANIA_OMEGA,
        TITANIA_I,
        TITANIA_NODE,
    )
    oberon = _moon_position_kepler(
        jd,
        OBERON_A_KM,
        OBERON_E,
        OBERON_N,
        OBERON_M0,
        OBERON_OMEGA,
        OBERON_I,
        OBERON_NODE,
    )
    miranda = _moon_position_kepler(
        jd,
        MIRANDA_A_KM,
        MIRANDA_E,
        MIRANDA_N,
        MIRANDA_M0,
        MIRANDA_OMEGA,
        MIRANDA_I,
        MIRANDA_NODE,
    )

    return (ariel, umbriel, titania, oberon, miranda)


def uranus_cob_offset(jd: float) -> Tuple[float, float, float]:
    """Calculate Uranus center-of-body offset from barycenter.

    The barycenter is offset from Uranus' center by:
    offset = -Σ(m_i * r_i) / M_total

    where m_i are moon masses and r_i are moon positions relative to Uranus.

    Args:
        jd: Julian Date (TT)

    Returns:
        Tuple (dx, dy, dz) offset in AU, in ICRF/J2000 frame.
        Add this to the barycenter position to get the planet center.
    """
    ariel, umbriel, titania, oberon, miranda = uranian_moon_positions(jd)

    # Weighted sum: barycenter offset from Uranus center
    total_gm = GM_URANUS + GM_URANUS_MOONS_TOTAL

    bary_x = (
        GM_ARIEL * ariel[0]
        + GM_UMBRIEL * umbriel[0]
        + GM_TITANIA * titania[0]
        + GM_OBERON * oberon[0]
        + GM_MIRANDA * miranda[0]
    ) / total_gm

    bary_y = (
        GM_ARIEL * ariel[1]
        + GM_UMBRIEL * umbriel[1]
        + GM_TITANIA * titania[1]
        + GM_OBERON * oberon[1]
        + GM_MIRANDA * miranda[1]
    ) / total_gm

    bary_z = (
        GM_ARIEL * ariel[2]
        + GM_UMBRIEL * umbriel[2]
        + GM_TITANIA * titania[2]
        + GM_OBERON * oberon[2]
        + GM_MIRANDA * miranda[2]
    ) / total_gm

    # Convert from km to AU
    # The offset from barycenter to COB is the negative of barycenter offset
    offset_x = -bary_x / AU_KM
    offset_y = -bary_y / AU_KM
    offset_z = -bary_z / AU_KM

    return (offset_x, offset_y, offset_z)
