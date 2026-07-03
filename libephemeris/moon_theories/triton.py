# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Triton orbital theory for Neptune barycenter correction.

Triton has a retrograde orbit (i ≈ 157°) and dominates Neptune's moon mass
(99.5% of total satellite mass). Its orbit is nearly circular but has
significant nodal precession due to Neptune's oblateness (J2).

Reference:
- Jacobson, R.A. (2009) "The Orbits of the Neptunian Satellites and the
  Orientation of the Pole of Neptune", AJ 137, 4322-4329
- Orbital elements from NEP097 ephemeris

Precision: ~20-50 km (sufficient for ~0.003 arcsec at Neptune's distance)
"""

from __future__ import annotations

import math
from typing import Tuple

# =============================================================================
# ORBITAL ELEMENTS (NEP097)
# =============================================================================

# Semi-major axis in km
TRITON_A_KM: float = 354759.0

# Orbital period in days (retrograde)
TRITON_PERIOD_DAYS: float = 5.8768541

# Mean motion in radians per day
TRITON_N_RAD_DAY: float = 2.0 * math.pi / TRITON_PERIOD_DAYS

# Eccentricity (nearly circular)
TRITON_E: float = 0.000016

# Reference epoch (J2000.0)
EPOCH_JD: float = 2451545.0

# Mean anomaly at epoch (radians) - from NEP097
M0_RAD: float = math.radians(358.656)

# Argument of perihelion (not well defined for e≈0)
OMEGA_RAD: float = math.radians(0.0)

# Inclination to Neptune's equator (retrograde!)
# ~156.865° to J2000 ecliptic
I_RAD: float = math.radians(156.865)

# Longitude of ascending node at J2000 (radians)
NODE_J2000_RAD: float = math.radians(177.608)

# Node precession rate (radians per day)
# Triton's node precesses with ~688 year period due to Neptune's J2
# This is retrograde precession for a retrograde orbit
NODE_RATE_RAD_DAY: float = -2.0 * math.pi / (688.0 * 365.25)

# AU in km for conversion
AU_KM: float = 149597870.7


def triton_position(jd: float) -> Tuple[float, float, float]:
    """
    Calculate Triton's position relative to Neptune in ICRF/J2000 coordinates.

    Uses Keplerian elements with secular nodal precession from Neptune's J2.

    Args:
        jd: Julian Date (TT)

    Returns:
        Tuple (x, y, z) in AU, Triton position relative to Neptune center
        in the J2000 ecliptic/equatorial frame (ICRF).
    """
    # Time since epoch in days
    dt = jd - EPOCH_JD

    # Mean anomaly at time t
    M = M0_RAD + TRITON_N_RAD_DAY * dt

    # Normalize to 0-2π
    M = M % (2.0 * math.pi)

    # Solve Kepler's equation for eccentric anomaly E
    # For e ≈ 0.000016, one iteration is sufficient
    E = M + TRITON_E * math.sin(M)

    # True anomaly
    cos_nu = (math.cos(E) - TRITON_E) / (1.0 - TRITON_E * math.cos(E))
    sin_nu = (math.sqrt(1.0 - TRITON_E**2) * math.sin(E)) / (
        1.0 - TRITON_E * math.cos(E)
    )
    nu = math.atan2(sin_nu, cos_nu)

    # Radius (distance from Neptune)
    r_km = TRITON_A_KM * (1.0 - TRITON_E**2) / (1.0 + TRITON_E * math.cos(nu))

    # Position in orbital plane
    x_orb = r_km * math.cos(nu)
    y_orb = r_km * math.sin(nu)

    # Current node position (with J2 precession)
    node = NODE_J2000_RAD + NODE_RATE_RAD_DAY * dt

    # Rotation matrices: orbital plane → J2000 ecliptic
    cos_omega = math.cos(OMEGA_RAD)
    sin_omega = math.sin(OMEGA_RAD)
    cos_i = math.cos(I_RAD)
    sin_i = math.sin(I_RAD)
    cos_node = math.cos(node)
    sin_node = math.sin(node)

    # Position after rotating by ω in orbital plane
    x1 = x_orb * cos_omega - y_orb * sin_omega
    y1 = x_orb * sin_omega + y_orb * cos_omega
    z1 = 0.0

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
