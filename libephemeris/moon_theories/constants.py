# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Physical constants for planetary moon theories.

All GM values are from JPL satellite ephemerides (2024).
Units: km³/s² for GM, km for distances.

References:
- JUP365: Jacobson (2021) - Jupiter system
- SAT441: Jacobson (2022) - Saturn system
- NEP097: Brozović (2020) - Neptune system
- PLU060: Brozović & Jacobson (2024) - Pluto system
"""

from __future__ import annotations

import math
from typing import Tuple

# =============================================================================
# GRAVITATIONAL PARAMETERS (GM) in km³/s²
# =============================================================================

# Planet GM values
GM_JUPITER: float = 126686531.9  # JUP365
GM_SATURN: float = 37931206.23  # SAT441
GM_URANUS: float = 5793951.3  # URA111
GM_NEPTUNE: float = 6835099.97  # NEP097
GM_PLUTO: float = 869.3  # PLU060

# Jupiter's Galilean moons (JUP365)
GM_IO: float = 5959.91547
GM_EUROPA: float = 3202.71210
GM_GANYMEDE: float = 9887.83275
GM_CALLISTO: float = 7179.28340

# Total Galilean moons GM
GM_GALILEAN_TOTAL: float = GM_IO + GM_EUROPA + GM_GANYMEDE + GM_CALLISTO

# Saturn's major moons (SAT441)
GM_MIMAS: float = 2.50349
GM_ENCELADUS: float = 7.21037
GM_TETHYS: float = 41.21353
GM_DIONE: float = 73.11607
GM_RHEA: float = 153.94175
GM_TITAN: float = 8978.13710
GM_HYPERION: float = 0.37049
GM_IAPETUS: float = 120.51511

# Total Saturn major moons GM (dominated by Titan)
GM_SATURN_MOONS_TOTAL: float = (
    GM_MIMAS
    + GM_ENCELADUS
    + GM_TETHYS
    + GM_DIONE
    + GM_RHEA
    + GM_TITAN
    + GM_HYPERION
    + GM_IAPETUS
)

# Neptune's moons (NEP097)
GM_TRITON: float = 1428.49546

# Uranus' major moons (URA111 / GUST86)
GM_ARIEL: float = 14.0
GM_UMBRIEL: float = 13.0
GM_TITANIA: float = 41.0
GM_OBERON: float = 29.0
GM_MIRANDA: float = 0.7

# Total Uranus major moons GM
GM_URANUS_MOONS_TOTAL: float = (
    GM_ARIEL + GM_UMBRIEL + GM_TITANIA + GM_OBERON + GM_MIRANDA
)

# Pluto's moons (PLU060)
GM_CHARON: float = 106.1

# =============================================================================
# ORBITAL CONSTANTS
# =============================================================================

# Jupiter equatorial radius (km) - for Galilean moon theory output
JUPITER_RADIUS_KM: float = 71492.0

# Astronomical Unit in km
AU_KM: float = 149597870.7

# Days per Julian century
DAYS_PER_CENTURY: float = 36525.0

# J2000.0 epoch
J2000_JD: float = 2451545.0

# TASS 1.7 epoch: JD 2444240.0 = 1980-01-01
TASS_EPOCH_JD: float = 2444240.0

# =============================================================================
# CHARON ORBITAL ELEMENTS (PLU060)
# =============================================================================
# Pluto-Charon system - binary with barycenter outside Pluto

CHARON_SEMIMAJOR_KM: float = 19591.0  # Charon distance from Pluto
CHARON_PERIOD_DAYS: float = 6.387230  # Orbital period
CHARON_ECCENTRICITY: float = 0.0002  # Nearly circular
CHARON_INCLINATION_DEG: float = 96.145  # To ecliptic

# Barycenter offset: Pluto center is offset from barycenter by:
# offset = a_charon * (M_charon / (M_pluto + M_charon))
# M_charon/M_total = GM_CHARON / (GM_PLUTO + GM_CHARON) ≈ 0.1088
PLUTO_BARYCENTER_OFFSET_KM: float = CHARON_SEMIMAJOR_KM * (
    GM_CHARON / (GM_PLUTO + GM_CHARON)
)  # ≈ 2130 km

# =============================================================================
# TRITON ORBITAL ELEMENTS (NEP097)
# =============================================================================

TRITON_SEMIMAJOR_KM: float = 354759.0
TRITON_PERIOD_DAYS: float = 5.876854
TRITON_ECCENTRICITY: float = 0.000016  # Nearly circular
TRITON_INCLINATION_DEG: float = 156.865  # Retrograde orbit!

# Neptune's J2 causes nodal precession with ~688 year period
TRITON_NODE_RATE_DEG_PER_DAY: float = -360.0 / (688.0 * 365.25)

# =============================================================================
# MASS RATIOS FOR BARYCENTER CORRECTION
# =============================================================================

# Jupiter: sum of Galilean moons / Jupiter
JUPITER_MOON_MASS_RATIO: float = GM_GALILEAN_TOTAL / GM_JUPITER  # ≈ 0.000207

# Saturn: Titan dominates (96% of moon mass)
SATURN_MOON_MASS_RATIO: float = GM_SATURN_MOONS_TOTAL / GM_SATURN  # ≈ 0.000247
SATURN_TITAN_RATIO: float = GM_TITAN / GM_SATURN  # ≈ 0.000237

# Neptune: Triton is 99.5% of moon mass
NEPTUNE_MOON_MASS_RATIO: float = GM_TRITON / GM_NEPTUNE  # ≈ 0.000209

# Pluto: Charon is 12.2% of system mass (binary!)
PLUTO_CHARON_MASS_RATIO: float = GM_CHARON / (GM_PLUTO + GM_CHARON)  # ≈ 0.109


# =============================================================================
# MAIN API: get_cob_offset
# =============================================================================


def get_cob_offset(
    planet_name: str,
    t,  # Skyfield Time object
) -> Tuple[float, float, float]:
    """
    Calculate the offset from planet barycenter to center of body (COB).

    This function computes the position of the planet's center relative to
    the system barycenter, based on analytical theories for the major moons.

    The offset is: COB = barycenter + offset

    Args:
        planet_name: Planet identifier (e.g., "jupiter barycenter", "saturn barycenter")
        t: Skyfield Time object

    Returns:
        Tuple (dx, dy, dz) offset in AU, in the equatorial ICRF/J2000
        frame (the frame Skyfield SSB vectors use).  Add this to the
        barycenter position to get the planet center.

    Notes:
        - For Jupiter: uses Galilean moon positions (E5/Meeus theory)
        - For Saturn: uses TASS 1.7 theory (primarily Titan)
        - For Uranus: uses Keplerian elements for 5 major moons
        - For Neptune: uses Triton Keplerian elements
        - For Pluto: uses Charon 2-body solution
        - For other planets: returns (0, 0, 0)

        The moon theories output the J2000 ecliptic frame (the Uranus
        theory outputs equatorial ICRF directly), so this function
        applies the fixed ecliptic -> equatorial rotation before
        returning.
    """
    # Get Julian Date from Skyfield Time
    jd = t.tt

    # Normalize planet name
    name_lower = planet_name.lower().replace(" barycenter", "").strip()

    if name_lower == "jupiter":
        return _ecliptic_to_equatorial_j2000(_jupiter_cob_offset(jd))
    elif name_lower == "saturn":
        return _ecliptic_to_equatorial_j2000(_saturn_cob_offset(jd))
    elif name_lower == "uranus":
        return _uranus_cob_offset(jd)  # already equatorial ICRF
    elif name_lower == "neptune":
        return _ecliptic_to_equatorial_j2000(_neptune_cob_offset(jd))
    elif name_lower in ("pluto", "pluto barycenter"):
        return _ecliptic_to_equatorial_j2000(_pluto_cob_offset(jd))
    else:
        # No correction needed for other planets
        return (0.0, 0.0, 0.0)


# Mean obliquity at J2000.0: 84381.406 arcsec (IAU 2006)
_EPS_J2000_RAD: float = math.radians(84381.406 / 3600.0)
_COS_EPS: float = math.cos(_EPS_J2000_RAD)
_SIN_EPS: float = math.sin(_EPS_J2000_RAD)


def _ecliptic_to_equatorial_j2000(
    xyz: Tuple[float, float, float],
) -> Tuple[float, float, float]:
    """Rotate a J2000-ecliptic vector to the equatorial ICRF frame."""
    x, y, z = xyz
    return (
        x,
        y * _COS_EPS - z * _SIN_EPS,
        y * _SIN_EPS + z * _COS_EPS,
    )


def _jupiter_cob_offset(jd: float) -> Tuple[float, float, float]:
    """
    Calculate Jupiter COB offset using Galilean moon positions.

    The barycenter is offset from Jupiter's center by:
    offset = -Σ(m_i * r_i) / M_total

    where m_i are moon masses and r_i are moon positions relative to Jupiter.
    """
    from .galilean import galilean_moon_positions

    # Get moon positions relative to Jupiter in km (J2000 ecliptic)
    io_xyz, europa_xyz, ganymede_xyz, callisto_xyz = galilean_moon_positions(jd)

    # Weighted sum: barycenter offset from Jupiter center
    total_gm = GM_JUPITER + GM_GALILEAN_TOTAL

    bary_x = (
        GM_IO * io_xyz[0]
        + GM_EUROPA * europa_xyz[0]
        + GM_GANYMEDE * ganymede_xyz[0]
        + GM_CALLISTO * callisto_xyz[0]
    ) / total_gm

    bary_y = (
        GM_IO * io_xyz[1]
        + GM_EUROPA * europa_xyz[1]
        + GM_GANYMEDE * ganymede_xyz[1]
        + GM_CALLISTO * callisto_xyz[1]
    ) / total_gm

    bary_z = (
        GM_IO * io_xyz[2]
        + GM_EUROPA * europa_xyz[2]
        + GM_GANYMEDE * ganymede_xyz[2]
        + GM_CALLISTO * callisto_xyz[2]
    ) / total_gm

    # Convert from km to AU
    # The offset from barycenter to COB is the negative of barycenter offset
    offset_x = -bary_x / AU_KM
    offset_y = -bary_y / AU_KM
    offset_z = -bary_z / AU_KM

    return (offset_x, offset_y, offset_z)


def _saturn_cob_offset(jd: float) -> Tuple[float, float, float]:
    """
    Calculate Saturn COB offset using TASS 1.7 theory.

    Titan dominates (96% of moon mass), so we primarily use Titan's position.
    """
    from .tass17 import saturn_moon_position

    # Get Titan position (index 5 in TASS 1.7)
    titan_xyz = saturn_moon_position(jd, 5)  # Returns (x, y, z) in AU

    # For more precision, we could add other moons, but Titan dominates
    # Mass ratio: Titan / (Saturn + all moons)
    total_gm = GM_SATURN + GM_SATURN_MOONS_TOTAL

    # Barycenter offset from Saturn center (in AU)
    bary_x = GM_TITAN * titan_xyz[0] / total_gm
    bary_y = GM_TITAN * titan_xyz[1] / total_gm
    bary_z = GM_TITAN * titan_xyz[2] / total_gm

    # COB offset = -barycenter offset
    return (-bary_x, -bary_y, -bary_z)


def _uranus_cob_offset(jd: float) -> Tuple[float, float, float]:
    """
    Calculate Uranus COB offset using Keplerian elements for major moons.

    Titania and Oberon dominate (~77% of moon mass), but we include all 5
    major moons for better accuracy.

    Precision: ~100-200 km (~0.005-0.01 arcsec at Uranus opposition)
    """
    from .uranian import uranus_cob_offset

    return uranus_cob_offset(jd)


def _neptune_cob_offset(jd: float) -> Tuple[float, float, float]:
    """
    Calculate Neptune COB offset using Triton's Keplerian elements.

    Triton has a retrograde orbit and dominates Neptune's moon mass.
    """
    from .triton import triton_position

    # Get Triton position relative to Neptune in AU (ICRF)
    triton_xyz = triton_position(jd)

    # Mass ratio
    total_gm = GM_NEPTUNE + GM_TRITON

    # Barycenter offset from Neptune center
    bary_x = GM_TRITON * triton_xyz[0] / total_gm
    bary_y = GM_TRITON * triton_xyz[1] / total_gm
    bary_z = GM_TRITON * triton_xyz[2] / total_gm

    # COB offset = -barycenter offset
    return (-bary_x, -bary_y, -bary_z)


def _pluto_cob_offset(jd: float) -> Tuple[float, float, float]:
    """
    Calculate Pluto COB offset using Charon's 2-body solution.

    Pluto-Charon is a binary system with the barycenter outside Pluto.
    """
    from .charon import charon_position

    # Get Charon position relative to Pluto in AU (ICRF)
    charon_xyz = charon_position(jd)

    # Mass ratio (Charon is ~12% of system mass)
    total_gm = GM_PLUTO + GM_CHARON

    # Barycenter offset from Pluto center
    bary_x = GM_CHARON * charon_xyz[0] / total_gm
    bary_y = GM_CHARON * charon_xyz[1] / total_gm
    bary_z = GM_CHARON * charon_xyz[2] / total_gm

    # COB offset = -barycenter offset
    return (-bary_x, -bary_y, -bary_z)
