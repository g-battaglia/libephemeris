# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Charon orbital theory for Pluto-Charon barycenter correction.

This implements a simple two-body Keplerian solution for Charon's orbit
around Pluto. The Pluto-Charon system is essentially a binary, with the
barycenter located ~2130 km from Pluto's center (outside the planet!).

Reference:
- Brozović, M. & Jacobson, R.A. (2024) "Post-New Horizons orbits and masses
  for the satellites of Pluto", AJ 167:256
- Orbital elements from PLU060 ephemeris

Precision: ~1-5 km (sufficient for ~0.01 arcsec at Pluto's distance)

Provenance:
    Orbital period, scale, orientation, phase, and mass inputs are attributed to
    the cited Brozovic/Jacobson Pluto-system solution and its public JPL
    ephemeris record. The runtime uses an independently written two-body Kepler
    solver. Treating unavailable secular rates as zero and the quoted accuracy
    are explicit approximation limits; no compatibility output is fitted.
"""

import math
from typing import Tuple

# =============================================================================
# ORBITAL ELEMENTS (PLU060, epoch J2000)
# =============================================================================

# Semi-major axis in km
CHARON_A_KM: float = 19591.4

# Orbital period in days
CHARON_PERIOD_DAYS: float = 6.3872304

# Mean motion in radians per day
CHARON_N_RAD_DAY: float = 2.0 * math.pi / CHARON_PERIOD_DAYS

# Eccentricity (nearly circular)
CHARON_E: float = 0.0002

# Reference epoch (J2000.0)
EPOCH_JD: float = 2451545.0

# Mean anomaly at epoch (radians) - from PLU060
M0_RAD: float = math.radians(71.255)

# Argument of perihelion at epoch (radians)
OMEGA_RAD: float = math.radians(180.0)

# Inclination to J2000 ecliptic (radians)
# Pluto's orbital plane is highly inclined
I_RAD: float = math.radians(96.145)

# Longitude of ascending node at epoch (radians)
NODE_RAD: float = math.radians(223.046)

# Precession rate of the node (very slow for Pluto system)
NODE_RATE_RAD_DAY: float = 0.0  # Negligible over centuries

# AU in km for conversion
AU_KM: float = 149597870.7


def charon_position(jd: float) -> Tuple[float, float, float]:
    """
    Calculate Charon's position relative to Pluto in ICRF/J2000 coordinates.

    Args:
        jd: Julian Date (TT)

    Returns:
        Tuple (x, y, z) in AU, Charon position relative to Pluto center
        in the J2000 ecliptic/equatorial frame (ICRF).
    """
    # Time since epoch in days
    dt = jd - EPOCH_JD

    # Mean anomaly at time t
    M = M0_RAD + CHARON_N_RAD_DAY * dt

    # Normalize to 0-2π
    M = M % (2.0 * math.pi)

    # Solve Kepler's equation for eccentric anomaly E
    # For small eccentricity, use simple iteration
    E = M
    for _ in range(5):  # 5 iterations is plenty for e=0.0002
        E = M + CHARON_E * math.sin(E)

    # True anomaly
    cos_nu = (math.cos(E) - CHARON_E) / (1.0 - CHARON_E * math.cos(E))
    sin_nu = (math.sqrt(1.0 - CHARON_E**2) * math.sin(E)) / (
        1.0 - CHARON_E * math.cos(E)
    )
    nu = math.atan2(sin_nu, cos_nu)

    # Radius (distance from Pluto)
    r_km = CHARON_A_KM * (1.0 - CHARON_E**2) / (1.0 + CHARON_E * math.cos(nu))

    # Position in orbital plane (x towards perihelion, y perpendicular)
    x_orb = r_km * math.cos(nu)
    y_orb = r_km * math.sin(nu)
    z_orb = 0.0

    # Current node position (with slow precession)
    node = NODE_RAD + NODE_RATE_RAD_DAY * dt

    # Rotation matrices: orbital plane → J2000 ecliptic
    # 1. Rotate by argument of perihelion (ω)
    # 2. Rotate by inclination (i)
    # 3. Rotate by longitude of node (Ω)

    cos_omega = math.cos(OMEGA_RAD)
    sin_omega = math.sin(OMEGA_RAD)
    cos_i = math.cos(I_RAD)
    sin_i = math.sin(I_RAD)
    cos_node = math.cos(node)
    sin_node = math.sin(node)

    # Combined rotation (standard orbital → ecliptic transformation)
    # Using the standard aerospace rotation sequence

    # Position after rotating by ω in orbital plane
    x1 = x_orb * cos_omega - y_orb * sin_omega
    y1 = x_orb * sin_omega + y_orb * cos_omega
    z1 = z_orb

    # Rotate by inclination around x-axis
    x2 = x1
    y2 = y1 * cos_i - z1 * sin_i
    z2 = y1 * sin_i + z1 * cos_i

    # Rotate by node around z-axis
    x_ecl = x2 * cos_node - y2 * sin_node
    y_ecl = x2 * sin_node + y2 * cos_node
    z_ecl = z2

    # Convert to AU
    x_au = x_ecl / AU_KM
    y_au = y_ecl / AU_KM
    z_au = z_ecl / AU_KM

    return (x_au, y_au, z_au)
