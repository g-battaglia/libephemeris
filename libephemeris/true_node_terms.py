# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
ELP2000-82B True Node Perturbation Term Table for libephemeris.

This module contains the perturbation coefficients for computing the True Lunar Node
based on the ELP2000-82B theory by Chapront-Touze & Chapront.

The table includes 99 terms organized into the categories returned by
get_terms_by_category():
1. Main Solar Perturbation Terms (main_solar)
2. Second-Order Solar Terms (second_order_solar)
3. Third-Order Terms (third_order)
4. F-Related (Inclination) Terms (f_inclination)
5. Venus Perturbation Terms (venus)
6. Mars Perturbation Terms (mars)
7. Jupiter Perturbation Terms (jupiter)
8. Saturn Perturbation Terms (saturn)
9. Evection Terms (evection)
10. Variation Terms (variation)
11. Annual-Equation Terms (annual_equation)
12. Second-Order Coupling Terms (second_order_coupling)
13. Parallactic Terms (parallactic)

Each term is represented as a NamedTuple containing:
- amplitude: Coefficient in degrees
- d_mult: Multiplier for mean elongation D
- m_mult: Multiplier for Sun's mean anomaly M
- mp_mult: Multiplier for Moon's mean anomaly M'
- f_mult: Multiplier for argument of latitude F
- e_power: Power of Earth's eccentricity correction (0, 1, or 2)
- trig_func: Trigonometric function ('sin' or 'cos')
- planet: Optional planetary longitude argument ('venus', 'mars', 'jupiter', 'saturn', None)
- planet_mult: Multiplier for planetary longitude (default 1)
- description: Human-readable description of the term

References:
    - Chapront-Touze, M. & Chapront, J. "ELP 2000-82B: A semi-analytical
      lunar ephemeris adequate for historical times" (1988)
    - Chapront-Touze, M. & Chapront, J. "Lunar Tables and Programs from
      4000 B.C. to A.D. 8000" (1991)
    - Meeus, J. "Astronomical Algorithms" (2nd ed., 1998), Chapter 47
    - Brown, E.W. "Tables of the Moon" (1919)
"""

from typing import NamedTuple, Optional, Tuple


class TrueNodeTerm(NamedTuple):
    """
    A single perturbation term for the True Lunar Node calculation.

    The term contributes to the perturbation as:
        amplitude * E^e_power * trig_func(d_mult*D + m_mult*M + mp_mult*M' + f_mult*F
                                          + planet_mult*L_planet)

    where:
        - D is the mean elongation of Moon from Sun
        - M is the mean anomaly of the Sun
        - M' is the mean anomaly of the Moon
        - F is the mean argument of latitude of the Moon
        - E is the eccentricity of Earth's orbit
        - L_planet is the mean longitude of the specified planet (if any)
    """

    amplitude: float  # Coefficient in degrees
    d_mult: int  # Multiplier for D (mean elongation)
    m_mult: int  # Multiplier for M (Sun's mean anomaly)
    mp_mult: int  # Multiplier for M' (Moon's mean anomaly)
    f_mult: int  # Multiplier for F (argument of latitude)
    e_power: int  # Power of E (0, 1, or 2)
    trig_func: str  # 'sin' or 'cos'
    planet: Optional[str]  # 'venus', 'mars', 'jupiter', 'saturn', or None
    planet_mult: int  # Multiplier for planetary longitude
    description: str  # Human-readable description


# =============================================================================
# MAIN SOLAR PERTURBATION TERMS (First Order)
# =============================================================================
# These are the dominant terms arising from solar gravitational influence.
# The largest term (-1.5233 sin(2D)) represents the fortnightly oscillation.
# Amplitudes are in degrees.

MAIN_SOLAR_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=-1.5233,
        d_mult=2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Dominant fortnightly term (2D)",
    ),
    TrueNodeTerm(
        amplitude=0.0595,
        d_mult=2,
        m_mult=-1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Solar-elongation coupling (2D-M)",
    ),
    TrueNodeTerm(
        amplitude=-0.0145,
        d_mult=2,
        m_mult=1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Solar-elongation coupling (2D+M)",
    ),
    TrueNodeTerm(
        amplitude=-0.1226,
        d_mult=2,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Evection primary term (2D-M')",
    ),
    TrueNodeTerm(
        amplitude=0.0490,
        d_mult=2,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Lunar anomaly term (2D+M')",
    ),
    TrueNodeTerm(
        amplitude=0.0316,
        d_mult=0,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Direct Moon anomaly term (M')",
    ),
    TrueNodeTerm(
        amplitude=0.1176,
        d_mult=0,
        m_mult=0,
        mp_mult=0,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude argument term (2F)",
    ),
    TrueNodeTerm(
        amplitude=-0.0801,
        d_mult=2,
        m_mult=0,
        mp_mult=0,
        f_mult=-2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Elongation-latitude term (2D-2F)",
    ),
    TrueNodeTerm(
        amplitude=0.0122,
        d_mult=2,
        m_mult=0,
        mp_mult=1,
        f_mult=-2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Combined term (2D-2F+M')",
    ),
)

# =============================================================================
# SECOND-ORDER SOLAR PERTURBATION TERMS
# =============================================================================

SECOND_ORDER_SOLAR_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0187,
        d_mult=2,
        m_mult=-1,
        mp_mult=-1,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Sun-Moon combined term (2D-M-M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0095,
        d_mult=2,
        m_mult=-1,
        mp_mult=1,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Sun-Moon combined term (2D-M+M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0154,
        d_mult=1,
        m_mult=1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Synodic-solar term (D+M)",
    ),
    TrueNodeTerm(
        amplitude=-0.0144,
        d_mult=1,
        m_mult=-1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Synodic-solar term (D-M)",
    ),
    TrueNodeTerm(
        amplitude=-0.0392,
        d_mult=2,
        m_mult=0,
        mp_mult=-2,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Double Moon anomaly term (2D-2M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0232,
        d_mult=0,
        m_mult=0,
        mp_mult=2,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Second harmonic Moon anomaly (2M')",
    ),
    TrueNodeTerm(
        amplitude=0.0109,
        d_mult=2,
        m_mult=0,
        mp_mult=2,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Combined term (2D+2M')",
    ),
    TrueNodeTerm(
        amplitude=0.0093,
        d_mult=0,
        m_mult=0,
        mp_mult=-1,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude-anomaly term (2F-M')",
    ),
    TrueNodeTerm(
        amplitude=0.0072,
        d_mult=0,
        m_mult=0,
        mp_mult=1,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude-anomaly term (2F+M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0279,
        d_mult=1,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Parallactic inequality (D)",
    ),
    TrueNodeTerm(
        amplitude=0.0054,
        d_mult=3,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Third harmonic elongation (3D)",
    ),
    TrueNodeTerm(
        amplitude=-0.0038,
        d_mult=4,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Fourth harmonic elongation (4D)",
    ),
)

# =============================================================================
# THIRD-ORDER AND HIGHER SOLAR TERMS
# =============================================================================

THIRD_ORDER_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0058,
        d_mult=0,
        m_mult=1,
        mp_mult=1,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Sun-Moon direct term (M+M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0053,
        d_mult=0,
        m_mult=1,
        mp_mult=-1,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Sun-Moon difference term (M-M')",
    ),
    TrueNodeTerm(
        amplitude=0.0038,
        d_mult=2,
        m_mult=-2,
        mp_mult=0,
        f_mult=0,
        e_power=2,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="E-squared term (2D-2M)",
    ),
    TrueNodeTerm(
        amplitude=0.0031,
        d_mult=2,
        m_mult=0,
        mp_mult=-3,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Third harmonic Moon anomaly (2D-3M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0025,
        d_mult=2,
        m_mult=0,
        mp_mult=3,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Third harmonic Moon anomaly (2D+3M')",
    ),
    TrueNodeTerm(
        amplitude=0.0023,
        d_mult=2,
        m_mult=1,
        mp_mult=-1,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Combined term (2D+M-M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0021,
        d_mult=2,
        m_mult=-1,
        mp_mult=2,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Combined term (2D-M+2M')",
    ),
    TrueNodeTerm(
        amplitude=0.0019,
        d_mult=0,
        m_mult=0,
        mp_mult=3,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Third harmonic Moon anomaly (3M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0017,
        d_mult=0,
        m_mult=1,
        mp_mult=2,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Sun-Moon combined (M+2M')",
    ),
    TrueNodeTerm(
        amplitude=0.0015,
        d_mult=0,
        m_mult=1,
        mp_mult=-2,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Sun-Moon combined (M-2M')",
    ),
)

# =============================================================================
# F-RELATED TERMS (Inclination Effects)
# =============================================================================

F_INCLINATION_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=-0.0086,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude-elongation term (2F-2D)",
    ),
    TrueNodeTerm(
        amplitude=0.0064,
        d_mult=2,
        m_mult=0,
        mp_mult=0,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude-elongation term (2F+2D)",
    ),
    TrueNodeTerm(
        amplitude=-0.0046,
        d_mult=1,
        m_mult=0,
        mp_mult=0,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude-elongation term (2F+D)",
    ),
    TrueNodeTerm(
        amplitude=0.0039,
        d_mult=-1,
        m_mult=0,
        mp_mult=0,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude-elongation term (2F-D)",
    ),
    TrueNodeTerm(
        amplitude=0.0032,
        d_mult=0,
        m_mult=-1,
        mp_mult=0,
        f_mult=2,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude-solar term (2F-M)",
    ),
    TrueNodeTerm(
        amplitude=-0.0028,
        d_mult=0,
        m_mult=1,
        mp_mult=0,
        f_mult=2,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Latitude-solar term (2F+M)",
    ),
    TrueNodeTerm(
        amplitude=0.0024,
        d_mult=0,
        m_mult=0,
        mp_mult=0,
        f_mult=4,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Fourth harmonic latitude (4F)",
    ),
    TrueNodeTerm(
        amplitude=-0.0018,
        d_mult=-2,
        m_mult=0,
        mp_mult=1,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Combined term (2F-2D+M')",
    ),
    TrueNodeTerm(
        amplitude=0.0016,
        d_mult=2,
        m_mult=0,
        mp_mult=-1,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Combined term (2F+2D-M')",
    ),
)

# =============================================================================
# VENUS PERTURBATION TERMS
# =============================================================================
# Venus gravitational influence on lunar orbit.

VENUS_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0048,
        d_mult=0,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="venus",
        planet_mult=1,
        description="Venus direct term (L_Venus - L_Moon)",
    ),
    TrueNodeTerm(
        amplitude=-0.0037,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="venus",
        planet_mult=1,
        description="Venus-elongation term (L_Venus - 2D)",
    ),
    TrueNodeTerm(
        amplitude=0.0029,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="venus",
        planet_mult=2,
        description="Venus double term (2L_Venus - 2D)",
    ),
    TrueNodeTerm(
        amplitude=-0.0024,
        d_mult=0,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="venus",
        planet_mult=1,
        description="Venus-anomaly term (L_Venus + M')",
    ),
    TrueNodeTerm(
        amplitude=0.0021,
        d_mult=0,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="venus",
        planet_mult=1,
        description="Venus-anomaly term (L_Venus - M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0018,
        d_mult=-2,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="venus",
        planet_mult=1,
        description="Venus combined (L_Venus - 2D + M')",
    ),
    TrueNodeTerm(
        amplitude=0.0015,
        d_mult=-2,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="venus",
        planet_mult=1,
        description="Venus combined (L_Venus - 2D - M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0012,
        d_mult=0,
        m_mult=1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet="venus",
        planet_mult=1,
        description="Venus-solar term (L_Venus + M)",
    ),
    TrueNodeTerm(
        amplitude=0.0010,
        d_mult=0,
        m_mult=-1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet="venus",
        planet_mult=1,
        description="Venus-solar term (L_Venus - M)",
    ),
)

# =============================================================================
# MARS PERTURBATION TERMS
# =============================================================================
# Mars gravitational influence on lunar orbit.

MARS_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0036,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="mars",
        planet_mult=1,
        description="Mars-elongation term (L_Mars - 2D)",
    ),
    TrueNodeTerm(
        amplitude=-0.0027,
        d_mult=0,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="mars",
        planet_mult=1,
        description="Mars direct term (L_Mars)",
    ),
    TrueNodeTerm(
        amplitude=0.0022,
        d_mult=0,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="mars",
        planet_mult=1,
        description="Mars-anomaly term (L_Mars - M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0018,
        d_mult=0,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="mars",
        planet_mult=1,
        description="Mars-anomaly term (L_Mars + M')",
    ),
    TrueNodeTerm(
        amplitude=0.0017,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="mars",
        planet_mult=2,
        description="Mars double term (2L_Mars - 2D)",
    ),
    TrueNodeTerm(
        amplitude=0.0014,
        d_mult=-2,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="mars",
        planet_mult=1,
        description="Mars combined (L_Mars - 2D + M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0011,
        d_mult=-2,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="mars",
        planet_mult=1,
        description="Mars combined (L_Mars - 2D - M')",
    ),
    TrueNodeTerm(
        amplitude=0.0009,
        d_mult=0,
        m_mult=-1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet="mars",
        planet_mult=1,
        description="Mars-solar term (L_Mars - M)",
    ),
    TrueNodeTerm(
        amplitude=-0.0008,
        d_mult=0,
        m_mult=1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet="mars",
        planet_mult=1,
        description="Mars-solar term (L_Mars + M)",
    ),
)

# =============================================================================
# JUPITER PERTURBATION TERMS
# =============================================================================
# Jupiter gravitational influence on lunar orbit.

JUPITER_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0033,
        d_mult=0,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="jupiter",
        planet_mult=1,
        description="Jupiter direct term (L_Jupiter)",
    ),
    TrueNodeTerm(
        amplitude=-0.0028,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="jupiter",
        planet_mult=1,
        description="Jupiter-elongation term (L_Jupiter - 2D)",
    ),
    TrueNodeTerm(
        amplitude=0.0021,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="jupiter",
        planet_mult=2,
        description="Jupiter double term (2L_Jupiter - 2D)",
    ),
    TrueNodeTerm(
        amplitude=-0.0017,
        d_mult=0,
        m_mult=1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet="jupiter",
        planet_mult=1,
        description="Jupiter-solar term (L_Jupiter + M)",
    ),
    TrueNodeTerm(
        amplitude=0.0014,
        d_mult=0,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="jupiter",
        planet_mult=1,
        description="Jupiter-anomaly term (L_Jupiter - M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0012,
        d_mult=0,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="jupiter",
        planet_mult=1,
        description="Jupiter-anomaly term (L_Jupiter + M')",
    ),
    TrueNodeTerm(
        amplitude=0.0009,
        d_mult=-2,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="jupiter",
        planet_mult=1,
        description="Jupiter combined (L_Jupiter - 2D + M')",
    ),
)

# =============================================================================
# SATURN PERTURBATION TERMS
# =============================================================================
# Saturn gravitational influence on lunar orbit.

SATURN_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0026,
        d_mult=0,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="saturn",
        planet_mult=1,
        description="Saturn direct term (L_Saturn)",
    ),
    TrueNodeTerm(
        amplitude=-0.0022,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="saturn",
        planet_mult=1,
        description="Saturn-elongation term (L_Saturn - 2D)",
    ),
    TrueNodeTerm(
        amplitude=0.0018,
        d_mult=-2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="saturn",
        planet_mult=2,
        description="Saturn double term (2L_Saturn - 2D)",
    ),
    TrueNodeTerm(
        amplitude=-0.0014,
        d_mult=0,
        m_mult=1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet="saturn",
        planet_mult=1,
        description="Saturn-solar term (L_Saturn + M)",
    ),
    TrueNodeTerm(
        amplitude=0.0012,
        d_mult=0,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="saturn",
        planet_mult=1,
        description="Saturn-anomaly term (L_Saturn - M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0010,
        d_mult=0,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="saturn",
        planet_mult=1,
        description="Saturn-anomaly term (L_Saturn + M')",
    ),
    TrueNodeTerm(
        amplitude=0.0008,
        d_mult=-2,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet="saturn",
        planet_mult=1,
        description="Saturn combined (L_Saturn - 2D + M')",
    ),
)

# =============================================================================
# LONG-PERIOD TERMS (Evection, Variation, Annual Equation)
# =============================================================================

EVECTION_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0467,
        d_mult=2,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Primary evection term (2D-M')",
    ),
    TrueNodeTerm(
        amplitude=0.0156,
        d_mult=2,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Evection-anomaly coupling (evection_arg + M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0134,
        d_mult=2,
        m_mult=0,
        mp_mult=-2,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Evection-anomaly coupling (evection_arg - M')",
    ),
    TrueNodeTerm(
        amplitude=0.0089,
        d_mult=2,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Evection second harmonic (evection_arg + 2M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0072,
        d_mult=2,
        m_mult=0,
        mp_mult=-3,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Evection second harmonic (evection_arg - 2M')",
    ),
)

VARIATION_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0523,
        d_mult=2,
        m_mult=0,
        mp_mult=0,
        f_mult=1,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Primary variation-inclination term (2D+F)",
    ),
    TrueNodeTerm(
        amplitude=-0.0478,
        d_mult=2,
        m_mult=0,
        mp_mult=0,
        f_mult=-1,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Variation-inclination term (2D-F)",
    ),
    TrueNodeTerm(
        amplitude=0.0067,
        d_mult=4,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Double variation term (4D)",
    ),
    TrueNodeTerm(
        amplitude=0.0043,
        d_mult=4,
        m_mult=0,
        mp_mult=0,
        f_mult=1,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Double variation-inclination term (4D+F)",
    ),
    TrueNodeTerm(
        amplitude=-0.0039,
        d_mult=4,
        m_mult=0,
        mp_mult=0,
        f_mult=-1,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Double variation-inclination term (4D-F)",
    ),
)

ANNUAL_EQUATION_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=-0.1860,
        d_mult=0,
        m_mult=1,
        mp_mult=0,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Primary annual equation term (M)",
    ),
    TrueNodeTerm(
        amplitude=0.0098,
        d_mult=0,
        m_mult=1,
        mp_mult=1,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Annual-anomaly coupling (M+M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0082,
        d_mult=0,
        m_mult=1,
        mp_mult=-1,
        f_mult=0,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Annual-anomaly coupling (M-M')",
    ),
    TrueNodeTerm(
        amplitude=0.0037,
        d_mult=0,
        m_mult=1,
        mp_mult=0,
        f_mult=2,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Annual-inclination coupling (M+2F)",
    ),
    TrueNodeTerm(
        amplitude=-0.0032,
        d_mult=0,
        m_mult=1,
        mp_mult=0,
        f_mult=-2,
        e_power=1,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Annual-inclination coupling (M-2F)",
    ),
    TrueNodeTerm(
        amplitude=0.0024,
        d_mult=0,
        m_mult=2,
        mp_mult=0,
        f_mult=0,
        e_power=2,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Second harmonic annual equation (2M)",
    ),
)

# =============================================================================
# SECOND-ORDER COUPLING TERMS (Products of first-order perturbations)
# =============================================================================

SECOND_ORDER_COUPLING_TERMS: Tuple[TrueNodeTerm, ...] = (
    # Evection x Variation coupling
    TrueNodeTerm(
        amplitude=0.0024,
        d_mult=0,
        m_mult=0,
        mp_mult=1,
        f_mult=0,
        e_power=0,
        trig_func="cos",
        planet=None,
        planet_mult=0,
        description="Evection-Variation coupling: cos(M')",
    ),
    TrueNodeTerm(
        amplitude=-0.0018,
        d_mult=4,
        m_mult=0,
        mp_mult=-1,
        f_mult=0,
        e_power=0,
        trig_func="cos",
        planet=None,
        planet_mult=0,
        description="Evection-Variation coupling: cos(4D-M')",
    ),
    # Evection x Annual Equation coupling
    TrueNodeTerm(
        amplitude=0.0019,
        d_mult=2,
        m_mult=-1,
        mp_mult=-1,
        f_mult=0,
        e_power=1,
        trig_func="cos",
        planet=None,
        planet_mult=0,
        description="Evection-Annual coupling: cos(2D-M'-M)",
    ),
    TrueNodeTerm(
        amplitude=-0.0016,
        d_mult=2,
        m_mult=1,
        mp_mult=-1,
        f_mult=0,
        e_power=1,
        trig_func="cos",
        planet=None,
        planet_mult=0,
        description="Evection-Annual coupling: cos(2D-M'+M)",
    ),
    # Self-coupling terms
    TrueNodeTerm(
        amplitude=0.0014,
        d_mult=4,
        m_mult=0,
        mp_mult=-2,
        f_mult=0,
        e_power=0,
        trig_func="cos",
        planet=None,
        planet_mult=0,
        description="Evection self-coupling: cos(4D-2M')",
    ),
    TrueNodeTerm(
        amplitude=0.0011,
        d_mult=4,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="cos",
        planet=None,
        planet_mult=0,
        description="Variation self-coupling: cos(4D)",
    ),
    TrueNodeTerm(
        amplitude=0.0008,
        d_mult=0,
        m_mult=0,
        mp_mult=0,
        f_mult=2,
        e_power=0,
        trig_func="cos",
        planet=None,
        planet_mult=0,
        description="Latitude self-coupling: cos(2F)",
    ),
    TrueNodeTerm(
        amplitude=0.0006,
        d_mult=0,
        m_mult=2,
        mp_mult=0,
        f_mult=0,
        e_power=2,
        trig_func="cos",
        planet=None,
        planet_mult=0,
        description="Annual self-coupling: cos(2M)",
    ),
)

# =============================================================================
# PARALLACTIC TERMS
# =============================================================================

PARALLACTIC_TERMS: Tuple[TrueNodeTerm, ...] = (
    TrueNodeTerm(
        amplitude=0.0350,
        d_mult=1,
        m_mult=0,
        mp_mult=0,
        f_mult=0,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Primary parallactic inequality (D)",
    ),
    TrueNodeTerm(
        amplitude=0.0026,
        d_mult=2,
        m_mult=0,
        mp_mult=-1,
        f_mult=2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Secondary parallactic term (2D-M'+2F)",
    ),
    TrueNodeTerm(
        amplitude=-0.0022,
        d_mult=2,
        m_mult=0,
        mp_mult=1,
        f_mult=-2,
        e_power=0,
        trig_func="sin",
        planet=None,
        planet_mult=0,
        description="Secondary parallactic term (2D+M'-2F)",
    ),
)


# =============================================================================
# COMPLETE TERM TABLE
# =============================================================================
# All terms combined for convenience.

ALL_TRUE_NODE_TERMS: Tuple[TrueNodeTerm, ...] = (
    MAIN_SOLAR_TERMS
    + SECOND_ORDER_SOLAR_TERMS
    + THIRD_ORDER_TERMS
    + F_INCLINATION_TERMS
    + VENUS_TERMS
    + MARS_TERMS
    + JUPITER_TERMS
    + SATURN_TERMS
    + EVECTION_TERMS
    + VARIATION_TERMS
    + ANNUAL_EQUATION_TERMS
    + SECOND_ORDER_COUPLING_TERMS
    + PARALLACTIC_TERMS
)


def get_term_count() -> int:
    """Return the total number of terms in the table."""
    return len(ALL_TRUE_NODE_TERMS)


def get_terms_by_category() -> dict[str, Tuple[TrueNodeTerm, ...]]:
    """
    Return all terms organized by category.

    Returns:
        Dictionary mapping category names to tuples of TrueNodeTerm objects.
    """
    return {
        "main_solar": MAIN_SOLAR_TERMS,
        "second_order_solar": SECOND_ORDER_SOLAR_TERMS,
        "third_order": THIRD_ORDER_TERMS,
        "f_inclination": F_INCLINATION_TERMS,
        "venus": VENUS_TERMS,
        "mars": MARS_TERMS,
        "jupiter": JUPITER_TERMS,
        "saturn": SATURN_TERMS,
        "evection": EVECTION_TERMS,
        "variation": VARIATION_TERMS,
        "annual_equation": ANNUAL_EQUATION_TERMS,
        "second_order_coupling": SECOND_ORDER_COUPLING_TERMS,
        "parallactic": PARALLACTIC_TERMS,
    }


def get_sorted_terms_by_amplitude() -> list[TrueNodeTerm]:
    """
    Return all terms sorted by absolute amplitude (largest first).

    Returns:
        List of TrueNodeTerm objects sorted by |amplitude|.
    """
    return sorted(ALL_TRUE_NODE_TERMS, key=lambda t: abs(t.amplitude), reverse=True)
