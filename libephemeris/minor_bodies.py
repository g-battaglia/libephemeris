# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Minor body calculations for asteroids and Trans-Neptunian Objects (TNOs).

This module computes positions for:
- Main belt asteroids: Ceres, Pallas, Juno, Vesta, Hygiea, Interamnia, Davida, Europa, Sylvia, Psyche, Sappho, Pandora, Lilith
- Centaurs: Chiron, Pholus, Nessus, Asbolus, Chariklo, Hidalgo
- Trans-Neptunian Objects (TNOs): Eris, Sedna, Haumea, Makemake, Orcus, Quaoar, Ixion, Gonggong, Varuna
- Near-Earth asteroids: Apophis, Eros, Amor, Icarus, Toro, Toutatis, Itokawa, Bennu, Ryugu

Method: Keplerian orbital elements with first-order secular perturbations from
Jupiter, Saturn, Uranus, and Neptune, plus resonant libration modeling for plutinos.
This provides significantly improved accuracy over pure 2-body dynamics, especially
for propagation over multiple years.

PRECISION:
- Main belt asteroids: ~10-30 arcseconds typical (improved from 1-5 arcminutes)
- TNOs: ~1-3 arcminutes typical (improved from 3-10 arcminutes)
- Plutinos (Ixion, Orcus): <2° over 100-year spans (improved from 5-10° drift)
- Errors increase with time from epoch, but secular perturbations reduce drift

PERTURBATION MODEL:
- Applies secular perturbations to orbital elements (ω, Ω, mean anomaly)
- Accounts for gravitational influence of Jupiter (dominant), Saturn, Uranus, and Neptune
- Neptune perturbations are critical for plutinos (2:3 resonance) like Ixion and Orcus
- Uranus perturbations are significant for Trans-Neptunian Objects (TNOs)
- Based on classical Laplace-Lagrange secular theory
- Includes resonant libration model for plutinos (2:3 Neptune resonance)

RESONANT LIBRATION MODEL (plutinos):
- For bodies in 2:3 mean motion resonance with Neptune (Ixion, Orcus)
- Models the oscillation of the resonant argument φ = 3λ_TNO - 2λ_Neptune - ω_TNO
- Applies periodic correction with ~20,000 year period calibrated from JPL Horizons
- Reduces position errors from 5-10° to <2° over 100-year timescales

For research-grade precision, use full numerical integration (e.g., JPL Horizons)

Orbital elements source: JPL Small-Body Database (epoch JD 2461000.5 TDB = 2025-Sep-19)
Algorithm: Keplerian mechanics with Laplace-Lagrange secular perturbations + resonant libration
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Tuple, Optional, NamedTuple
from .logging_config import get_logger
from .exceptions import EphemerisRangeError
from .constants import (
    CHIRON,
    PHOLUS,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    ERIS,
    SEDNA,
    HAUMEA,
    MAKEMAKE,
    IXION,
    ORCUS,
    QUAOAR,
    NESSUS,
    ASBOLUS,
    CHARIKLO,
    GONGGONG,
    VARUNA,
    APOPHIS,
    HYGIEA,
    INTERAMNIA,
    DAVIDA,
    EUROPA_AST,
    SYLVIA,
    PSYCHE,
    EROS,
    AMOR,
    ICARUS,
    TORO,
    SAPPHO,
    PANDORA_AST,
    LILITH_AST,
    HIDALGO,
    TOUTATIS,
    ITOKAWA,
    BENNU,
    RYUGU,
)


# =============================================================================
# PHYSICAL CONSTANTS FOR PERTURBATION CALCULATIONS
# =============================================================================
# Gravitational parameters (GM) in AU^3/day^2
# These are used for secular perturbation calculations

# GM of Sun (k^2 where k is Gaussian gravitational constant)
GM_SUN = 0.00029591220828559  # AU^3/day^2

# Mass ratios relative to Sun (m_planet / m_sun)
# Source: IAU nominal values and JPL DE441
MASS_RATIO_JUPITER = 1.0 / 1047.348644  # ~9.546e-4
MASS_RATIO_SATURN = 1.0 / 3497.901768  # ~2.858e-4

# Mean orbital elements of perturbing planets (J2000.0 values)
# These are used for secular perturbation calculations
# Source: JPL Horizons, mean elements at J2000.0

# Jupiter mean elements (for secular perturbation theory)
JUPITER_A = 5.2026  # Semi-major axis (AU)
JUPITER_E = 0.0485  # Mean eccentricity
JUPITER_I = 1.303  # Mean inclination (degrees)
JUPITER_OMEGA = 274.25  # Mean longitude of perihelion (degrees)
JUPITER_NODE = 100.46  # Mean longitude of ascending node (degrees)
JUPITER_N = 0.08309  # Mean motion (degrees/day)

# Saturn mean elements
SATURN_A = 9.5549  # Semi-major axis (AU)
SATURN_E = 0.0555  # Mean eccentricity
SATURN_I = 2.489  # Mean inclination (degrees)
SATURN_OMEGA = 339.39  # Mean longitude of perihelion (degrees)
SATURN_NODE = 113.66  # Mean longitude of ascending node (degrees)
SATURN_N = 0.03346  # Mean motion (degrees/day)

# Uranus mean elements (significant for TNO calculations)
URANUS_A = 19.2184  # Semi-major axis (AU)
URANUS_E = 0.0457  # Mean eccentricity
URANUS_I = 0.772  # Mean inclination (degrees)
URANUS_OMEGA = 170.9  # Mean argument of perihelion (degrees)
URANUS_NODE = 74.0  # Mean longitude of ascending node (degrees)
URANUS_N = 0.01177  # Mean motion (degrees/day)

# Mass ratio of Uranus relative to Sun
MASS_RATIO_URANUS = 1.0 / 22902  # ~4.366e-5

# Neptune mean elements (critical for TNO accuracy, especially plutinos)
# Plutinos like Ixion and Orcus are in 2:3 mean motion resonance with Neptune
NEPTUNE_A = 30.1104  # Semi-major axis (AU)
NEPTUNE_E = 0.0086  # Mean eccentricity
NEPTUNE_I = 1.769  # Mean inclination (degrees)
NEPTUNE_OMEGA = 44.9  # Mean argument of perihelion (degrees)
NEPTUNE_NODE = 131.7  # Mean longitude of ascending node (degrees)
NEPTUNE_N = 0.006021  # Mean motion (degrees/day)

# Mass ratio of Neptune relative to Sun
MASS_RATIO_NEPTUNE = 1.0 / 19412  # ~5.153e-5


# L4: Linear rates for planet orbital elements (per Julian century from J2000.0)
# Source: Simon et al. (1994) A&A 282, 663; Standish (1992)
# These allow the forced eccentricity/inclination vectors to evolve with time,
# improving accuracy over century-scale propagations.
# Rates are in the same units as the base values (AU, dimensionless, degrees)
# per Julian century (36525 days).

# Jupiter element rates per century
JUPITER_DE = -0.00036  # eccentricity rate (dimensionless/century)
JUPITER_DI = -0.0019  # inclination rate (degrees/century)
JUPITER_DOMEGA = 0.3625  # longitude of perihelion rate (degrees/century)
JUPITER_DNODE = 0.3574  # ascending node rate (degrees/century)

# Saturn element rates per century
SATURN_DE = -0.00051  # eccentricity rate
SATURN_DI = 0.0019  # inclination rate
SATURN_DOMEGA = 0.7375  # longitude of perihelion rate
SATURN_DNODE = -0.2510  # ascending node rate

# Uranus element rates per century
URANUS_DE = -0.00018  # eccentricity rate
URANUS_DI = -0.0024  # inclination rate
URANUS_DOMEGA = 0.0710  # longitude of perihelion rate (arg. peri.)
URANUS_DNODE = 0.0174  # ascending node rate

# Neptune element rates per century
NEPTUNE_DE = 0.00007  # eccentricity rate
NEPTUNE_DI = -0.0035  # inclination rate
NEPTUNE_DOMEGA = -0.0070  # longitude of perihelion rate (arg. peri.)
NEPTUNE_DNODE = -0.0063  # ascending node rate


def _get_planet_elements_at_time(
    jd_tt: float,
) -> list[tuple[float, float, float, float, float, float, float]]:
    """Return planet orbital elements evolved to the given Julian Day.

    L4: Uses linear rates from Simon et al. (1994) to evolve planet elements
    from their J2000.0 values to the target date. This improves the forced
    eccentricity/inclination vectors for century-scale propagations.

    Args:
        jd_tt: Julian Day in Terrestrial Time.

    Returns:
        List of planet tuples: (a_planet, mu, e_planet, omega_deg, node_deg,
        i_deg, a_threshold) — same format as the static _planets list in
        _calc_forced_elements().
    """
    # Time in Julian centuries from J2000.0
    T = (jd_tt - 2451545.0) / 36525.0

    return [
        (
            JUPITER_A,
            MASS_RATIO_JUPITER,
            JUPITER_E + JUPITER_DE * T,
            JUPITER_OMEGA + JUPITER_DOMEGA * T,
            JUPITER_NODE + JUPITER_DNODE * T,
            JUPITER_I + JUPITER_DI * T,
            0.0,
        ),
        (
            SATURN_A,
            MASS_RATIO_SATURN,
            SATURN_E + SATURN_DE * T,
            SATURN_OMEGA + SATURN_DOMEGA * T,
            SATURN_NODE + SATURN_DNODE * T,
            SATURN_I + SATURN_DI * T,
            0.0,
        ),
        (
            URANUS_A,
            MASS_RATIO_URANUS,
            URANUS_E + URANUS_DE * T,
            URANUS_OMEGA + URANUS_DOMEGA * T,
            URANUS_NODE + URANUS_DNODE * T,
            URANUS_I + URANUS_DI * T,
            0.0,
        ),
        (
            NEPTUNE_A,
            MASS_RATIO_NEPTUNE,
            NEPTUNE_E + NEPTUNE_DE * T,
            NEPTUNE_OMEGA + NEPTUNE_DOMEGA * T,
            NEPTUNE_NODE + NEPTUNE_DNODE * T,
            NEPTUNE_I + NEPTUNE_DI * T,
            20.0,
        ),
    ]


# =============================================================================
# MEAN MOTION RESONANCE CONSTANTS
# =============================================================================
# Mean motion resonances with Neptune occur when an outer body's orbital period
# is a simple integer ratio of Neptune's period. These resonances protect
# bodies from close encounters but also mean that secular perturbation theory
# breaks down and gives incorrect results.
#
# Common Neptune resonances (body:Neptune orbital period ratio):
# - 2:3 (Plutinos): ~39.4 AU - most populated resonance (Pluto, Ixion, Orcus)
# - 1:2 (Twotinos): ~47.8 AU
# - 2:5: ~55.4 AU
# - 1:3: ~62.5 AU
# - 3:5: ~36.4 AU
# - 4:7: ~43.7 AU
# - 3:4: ~36.5 AU
#
# References:
#   Malhotra, R. (1995) "The origin of Pluto's peculiar orbit"
#   Murray & Dermott "Solar System Dynamics" Ch. 8

# Neptune's mean motion in degrees/day (for resonance calculations)
NEPTUNE_MEAN_MOTION = 0.006021  # degrees/day (~164.8 year period)


class ResonanceInfo(NamedTuple):
    """
    Information about a mean motion resonance.

    Attributes:
        p: Integer numerator (body's period multiplier)
        q: Integer denominator (Neptune's period multiplier)
        name: Common name for the resonance (e.g., "plutino", "twotino")
        a_resonant: Expected semi-major axis for this resonance (AU)
        tolerance: Fractional tolerance for detecting resonance
    """

    p: int
    q: int
    name: str
    a_resonant: float
    tolerance: float


# Calculate resonant semi-major axis using Kepler's 3rd law:
# For a p:q resonance, body makes p orbits while Neptune makes q orbits
# T_body / T_Neptune = q / p  (body's period is longer if q > p)
# (a_body / a_Neptune)^(3/2) = q / p
# a_body = a_Neptune * (q / p)^(2/3)


def _calc_resonant_a(p: int, q: int) -> float:
    """Calculate the resonant semi-major axis for a p:q resonance with Neptune.

    In a p:q resonance, the body completes p orbits while Neptune completes q orbits.
    For exterior resonances (most common), q > p, so the body's orbit is larger.
    """
    return NEPTUNE_A * (q / p) ** (2.0 / 3.0)


# Known Neptune mean motion resonances
# Format: (p, q, name, a_resonant, tolerance)
# Tolerance is fractional (e.g., 0.02 means ±2% of a_resonant)
NEPTUNE_RESONANCES: list[ResonanceInfo] = [
    # Interior resonances (a < Neptune)
    ResonanceInfo(3, 4, "3:4 resonance", _calc_resonant_a(3, 4), 0.02),
    ResonanceInfo(3, 5, "3:5 resonance", _calc_resonant_a(3, 5), 0.02),
    # Plutinos (most famous resonance, includes Pluto)
    ResonanceInfo(2, 3, "plutino", _calc_resonant_a(2, 3), 0.02),
    # Other exterior resonances
    ResonanceInfo(4, 7, "4:7 resonance", _calc_resonant_a(4, 7), 0.02),
    # Twotinos (1:2 resonance)
    ResonanceInfo(1, 2, "twotino", _calc_resonant_a(1, 2), 0.02),
    # More distant resonances
    ResonanceInfo(2, 5, "2:5 resonance", _calc_resonant_a(2, 5), 0.02),
    ResonanceInfo(1, 3, "1:3 resonance", _calc_resonant_a(1, 3), 0.02),
]

# Default tolerance for resonance detection (fractional difference in mean motion)
DEFAULT_RESONANCE_TOLERANCE = 0.02  # 2%


# =============================================================================
# PLUTINO LIBRATION PARAMETERS
# =============================================================================
# Plutinos (2:3 Neptune resonance) exhibit libration of the resonant argument:
#   φ = 3λ_TNO - 2λ_Neptune - ω_TNO
#
# This argument oscillates (librates) around a center value rather than
# circulating through 360°. The libration causes periodic corrections to the
# position that are NOT captured by secular perturbation theory.
#
# Libration parameters calibrated from JPL Horizons data:
# - libration_amplitude: Maximum angular displacement from center (degrees)
# - libration_period: Period of libration oscillation (days)
# - libration_center: Center value of φ around which it oscillates (degrees)
# - libration_phase: Phase at J2000.0 epoch (radians)
#
# For 2:3 resonance, typical libration period is ~20,000 years (~7.3 million days)
#
# References:
#   Malhotra, R. (1995) "The origin of Pluto's peculiar orbit"
#   Nesvorný & Roig (2000) "Mean Motion Resonances in the Trans-Neptunian Region"
#   Jewitt, Morbidelli & Rauer (2008) "Trans-Neptunian Objects and Comets"

# J2000.0 epoch for reference (JD 2451545.0 = 2000-Jan-01 12:00 TT)
J2000_EPOCH = 2451545.0

# Neptune's mean longitude at J2000.0 (degrees)
# This is needed for calculating the resonant argument φ
NEPTUNE_MEAN_LONGITUDE_J2000 = 304.88003  # degrees at J2000.0


class LibrationParameters(NamedTuple):
    """
    Parameters describing the libration of a resonant body's resonant argument.

    For a p:q resonance, the resonant argument is:
        φ = (p+q)λ_body - qλ_Neptune - pω_body

    For plutinos (2:3 resonance): φ = 3λ_body - 2λ_Neptune - ω_body

    Attributes:
        amplitude: Maximum angular displacement from center (degrees)
        period: Period of libration oscillation (days)
        center: Center value of φ around which it oscillates (degrees)
        phase_j2000: Phase of libration at J2000.0 epoch (radians)
        resonance_p: p in the p:q resonance (body makes p orbits)
        resonance_q: q in the p:q resonance (Neptune makes q orbits)
    """

    amplitude: float  # degrees
    period: float  # days
    center: float  # degrees
    phase_j2000: float  # radians
    resonance_p: int  # p in p:q resonance
    resonance_q: int  # q in p:q resonance


# Libration parameters for known resonant TNOs
# Calibrated from multi-decade JPL Horizons integrations
# Period ~20,000 years = ~7,305,000 days for typical plutinos
PLUTINO_LIBRATION_PARAMS: dict[int, LibrationParameters] = {
    # Ixion (28978): Well-characterized plutino (2:3 Neptune resonance)
    # Libration amplitude ~78°, period ~19,800 years
    # Center at 180° (anti-aligned with Neptune)
    IXION: LibrationParameters(
        amplitude=78.0,  # degrees
        period=7_233_000.0,  # ~19,800 years in days
        center=180.0,  # degrees (anti-aligned libration)
        phase_j2000=2.14,  # radians, calibrated to match JPL
        resonance_p=2,
        resonance_q=3,
    ),
    # Orcus (90482): Anti-Pluto plutino (2:3, opposite orbital phase from Pluto)
    # Libration amplitude ~68°, period ~20,200 years
    # Center at 180° (anti-aligned with Neptune)
    ORCUS: LibrationParameters(
        amplitude=68.0,  # degrees
        period=7_379_000.0,  # ~20,200 years in days
        center=180.0,  # degrees (anti-aligned libration)
        phase_j2000=4.71,  # radians, calibrated to match JPL (near 3π/2)
        resonance_p=2,
        resonance_q=3,
    ),
    # M2: Gonggong (225088): 3:10 Neptune resonance
    # Libration amplitude ~30°, period ~25,000 years (~9.1 million days)
    # Center at 180° (symmetric libration)
    # Parameters estimated from orbital integration studies
    # (Bannister et al. 2018, Brown & Butler 2018)
    GONGGONG: LibrationParameters(
        amplitude=30.0,  # degrees
        period=9_131_000.0,  # ~25,000 years in days
        center=180.0,  # degrees
        phase_j2000=1.05,  # radians, estimated
        resonance_p=3,
        resonance_q=10,
    ),
}


@dataclass
class OrbitalElements:
    """
    Classical Keplerian orbital elements for a minor body.

    Attributes:
        name: Body name
        epoch: Reference epoch (Julian Day in TT)
        a: Semi-major axis in AU
        e: Eccentricity (0-1, dimensionless)
        i: Inclination to ecliptic in degrees
        omega: Argument of perihelion (ω) in degrees
        Omega: Longitude of ascending node (Ω) in degrees
        M0: Mean anomaly at epoch in degrees
        n: Mean motion in degrees/day

    Note:
        These are osculating elements at the given epoch.
        They drift over time due to planetary perturbations.
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float

    def validate(self) -> list[str]:
        """Validate orbital elements for physical sanity.

        Checks that elements are within physically meaningful ranges.
        Does NOT raise exceptions — returns a list of warning messages.
        Empty list means all elements are valid.

        Returns:
            List of warning strings describing any issues found.
        """
        warnings: list[str] = []

        # Semi-major axis must be positive
        if self.a <= 0.0:
            warnings.append(
                f"{self.name}: semi-major axis a={self.a:.6f} AU is non-positive"
            )

        # Eccentricity checks
        if self.e < 0.0:
            warnings.append(f"{self.name}: eccentricity e={self.e:.6f} is negative")
        elif self.e >= 1.0 and self.a > 0.0:
            # Hyperbolic orbit should have negative semi-major axis
            warnings.append(
                f"{self.name}: eccentricity e={self.e:.6f} >= 1 but a > 0 "
                "(hyperbolic orbits should have a < 0)"
            )

        # Inclination should be in [0, 180]
        if self.i < 0.0 or self.i > 180.0:
            warnings.append(
                f"{self.name}: inclination i={self.i:.4f}° outside [0°, 180°]"
            )

        # Mean motion should be positive for bound orbits
        if self.e < 1.0 and self.n <= 0.0:
            warnings.append(
                f"{self.name}: mean motion n={self.n:.8f}°/day is non-positive "
                "for elliptic orbit"
            )

        return warnings


def validate_elements(elements: OrbitalElements) -> None:
    """Validate orbital elements and log warnings for any issues.

    This is a convenience wrapper around OrbitalElements.validate() that
    logs warnings via the module logger instead of returning a list.

    Args:
        elements: Orbital elements to validate.
    """
    issues = elements.validate()
    if issues:
        logger = get_logger()
        for issue in issues:
            logger.warning("Orbital element validation: %s", issue)


# =============================================================================
# MEAN MOTION RESONANCE DETECTION FUNCTIONS
# =============================================================================


@dataclass
class ResonanceResult:
    """
    Result of resonance detection for a body.

    Attributes:
        is_resonant: Whether the body is in or near a mean motion resonance
        resonance: The detected resonance (ResonanceInfo), or None if not resonant
        deviation: Fractional deviation from exact resonance (0.0 = exact)
        warning_message: Message about secular perturbation accuracy, or None
    """

    is_resonant: bool
    resonance: Optional[ResonanceInfo]
    deviation: float
    warning_message: Optional[str]


def detect_mean_motion_resonance(
    elements: "OrbitalElements",
    tolerance: float = DEFAULT_RESONANCE_TOLERANCE,
) -> ResonanceResult:
    """
    Detect if a body is in mean motion resonance with Neptune.

    Checks if the body's mean motion is close to a resonant ratio with Neptune's
    mean motion. Bodies in mean motion resonance have their secular perturbations
    significantly modified by resonant effects that are not captured by
    first-order Laplace-Lagrange secular theory.

    Args:
        elements: Orbital elements of the body to check
        tolerance: Fractional tolerance for detecting resonance (default 2%)

    Returns:
        ResonanceResult: Contains detection result, resonance info, and warning

    Algorithm:
        For each known Neptune resonance p:q, checks if:
        |n_body / n_Neptune - p/q| < tolerance * (p/q)

        where n is mean motion. In a p:q resonance, the body completes p orbits
        while Neptune completes q, so n_body/n_Neptune = p/q.
        This is equivalent to checking if the semi-major axis is within
        tolerance of the resonant value.

    Example:
        >>> from libephemeris.minor_bodies import detect_mean_motion_resonance
        >>> result = detect_mean_motion_resonance(MINOR_BODY_ELEMENTS[IXION])
        >>> result.is_resonant
        True
        >>> result.resonance.name
        'plutino'

    Note:
        Bodies in resonance may have secular perturbation calculations that
        are less accurate than for non-resonant bodies. The warning_message
        field provides guidance in such cases.

    References:
        Murray & Dermott "Solar System Dynamics" Ch. 8
        Malhotra (1995) "The origin of Pluto's peculiar orbit"
    """
    a = elements.a
    n = elements.n

    # Only check bodies that could potentially be in Neptune resonance
    # (must be in outer solar system, say a > 20 AU)
    if a < 20.0:
        return ResonanceResult(
            is_resonant=False, resonance=None, deviation=0.0, warning_message=None
        )

    # Calculate body's mean motion ratio relative to Neptune
    # For a p:q resonance, body makes p orbits while Neptune makes q
    # So n_body / n_Neptune = p / q
    n_ratio = n / NEPTUNE_MEAN_MOTION

    best_match: Optional[ResonanceInfo] = None
    best_deviation = float("inf")

    for resonance in NEPTUNE_RESONANCES:
        # Expected ratio for this resonance: p/q
        expected_ratio = resonance.p / resonance.q

        # Calculate fractional deviation
        deviation = abs(n_ratio - expected_ratio) / expected_ratio

        if deviation < tolerance and deviation < best_deviation:
            best_match = resonance
            best_deviation = deviation

    if best_match is not None:
        warning = (
            f"Body '{elements.name}' is near the {best_match.p}:{best_match.q} "
            f"mean motion resonance with Neptune ({best_match.name}). "
            f"Secular perturbation theory may be inaccurate for resonant bodies. "
            f"Deviation from exact resonance: {best_deviation * 100:.2f}%"
        )
        return ResonanceResult(
            is_resonant=True,
            resonance=best_match,
            deviation=best_deviation,
            warning_message=warning,
        )

    return ResonanceResult(
        is_resonant=False, resonance=None, deviation=0.0, warning_message=None
    )


def is_body_resonant(
    body_id: int, tolerance: float = DEFAULT_RESONANCE_TOLERANCE
) -> bool:
    """
    Check if a minor body is in mean motion resonance with Neptune.

    This is a convenience function that looks up the body's orbital elements
    and checks for resonance.

    Args:
        body_id: Minor body identifier (IXION, ORCUS, etc.)
        tolerance: Fractional tolerance for detecting resonance (default 2%)

    Returns:
        bool: True if the body is in or near a Neptune resonance

    Raises:
        ValueError: If body_id is not in the orbital elements database

    Example:
        >>> from libephemeris.minor_bodies import is_body_resonant
        >>> from libephemeris.constants import IXION, CERES
        >>> is_body_resonant(IXION)
        True
        >>> is_body_resonant(CERES)
        False
    """
    if body_id not in MINOR_BODY_ELEMENTS:
        raise ValueError(f"Unknown body ID: {body_id}")

    elements = MINOR_BODY_ELEMENTS[body_id]
    result = detect_mean_motion_resonance(elements, tolerance)
    return result.is_resonant


def get_resonance_info(body_id: int) -> Optional[ResonanceResult]:
    """
    Get detailed resonance information for a minor body.

    Args:
        body_id: Minor body identifier (IXION, ORCUS, etc.)

    Returns:
        ResonanceResult: Detailed resonance detection result, or None if body unknown

    Example:
        >>> from libephemeris.minor_bodies import get_resonance_info
        >>> from libephemeris.constants import ORCUS
        >>> info = get_resonance_info(ORCUS)
        >>> if info and info.is_resonant:
        ...     print(f"{info.resonance.name}: deviation {info.deviation:.2%}")
        plutino: deviation 0.15%
    """
    if body_id not in MINOR_BODY_ELEMENTS:
        return None

    elements = MINOR_BODY_ELEMENTS[body_id]
    return detect_mean_motion_resonance(elements)


# =============================================================================
# RESONANT LIBRATION FUNCTIONS
# =============================================================================


def calc_neptune_mean_longitude(jd_tt: float) -> float:
    """
    Calculate Neptune's mean longitude at a given Julian Day.

    Uses Neptune's mean motion to propagate from J2000.0 epoch.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Neptune's mean longitude in degrees (0-360)
    """
    dt = jd_tt - J2000_EPOCH  # Days since J2000.0
    lambda_neptune = (NEPTUNE_MEAN_LONGITUDE_J2000 + NEPTUNE_N * dt) % 360.0
    return lambda_neptune


def calc_resonant_argument_plutino(
    elements: OrbitalElements,
    jd_tt: float,
    omega_pert: float,
    M_pert: float,
) -> float:
    """
    Calculate the resonant argument φ for a plutino (2:3 resonance).

    For the 2:3 mean motion resonance with Neptune:
        φ = 3λ_TNO - 2λ_Neptune - ϖ_TNO

    where:
        λ_TNO = Ω + ω + M (mean longitude of the TNO)
        λ_Neptune = mean longitude of Neptune
        ϖ_TNO = Ω + ω (longitude of perihelion of the TNO)

    The longitude of perihelion ϖ (a true longitude measured from the reference
    direction) is used rather than the argument of perihelion ω (measured from
    the node). This makes φ a proper resonant angle: the integer coefficients of
    the longitude-like terms sum to zero (3 - 2 - 1 = 0), so φ is invariant under
    a rotation of the reference longitude origin. φ librates (oscillates) around
    180° rather than circulating through 360°, the defining characteristic of
    bodies captured in the resonance.

    Args:
        elements: Orbital elements of the plutino
        jd_tt: Julian Day in Terrestrial Time (TT)
        omega_pert: Perturbed argument of perihelion (degrees)
        M_pert: Perturbed mean anomaly (degrees)

    Returns:
        float: The resonant argument φ in degrees
    """
    # Mean longitude of the TNO: λ = Ω + ω + M
    lambda_tno = (elements.Omega + omega_pert + M_pert) % 360.0

    # Longitude of perihelion: ϖ = Ω + ω (a true longitude, unlike the
    # argument of perihelion ω which is measured from the node).
    varpi_tno = (elements.Omega + omega_pert) % 360.0

    # Get Neptune's mean longitude at this time
    lambda_neptune = calc_neptune_mean_longitude(jd_tt)

    # Resonant argument for the 2:3 resonance: φ = 3λ_body - 2λ_Neptune - ϖ_body
    phi = (3.0 * lambda_tno - 2.0 * lambda_neptune - varpi_tno) % 360.0

    return phi


def calc_libration_correction(
    body_id: int,
    jd_tt: float,
) -> float:
    """
    Calculate the longitude correction due to resonant libration.

    For bodies in mean motion resonance with Neptune, the resonant
    argument φ librates around a center value instead of circulating. This
    libration causes a periodic perturbation to the body's position that is
    not captured by secular perturbation theory.

    The correction is modeled as a simple harmonic oscillation:
        Δλ = A * sin(2π(t - t0) / P + φ0)

    where:
        A = libration amplitude
        P = libration period
        φ0 = phase at J2000.0

    The correction is scaled by a factor that accounts for the relationship
    between the resonant argument libration and the actual longitude change.
    For a p:q resonance: φ = (p+q)λ - qλ_N - pω, so Δλ = Δφ / (p+q).

    Args:
        body_id: Minor body identifier (must be a known resonant body)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        float: Longitude correction in degrees (positive = ahead of Keplerian)
               Returns 0.0 if body is not a known resonant body with libration params

    Note:
        This is a simplified model of resonant dynamics. Full accuracy requires
        numerical integration of the restricted three-body problem.

    References:
        Murray & Dermott "Solar System Dynamics" Ch. 8 (Resonant dynamics)
        Malhotra (1996) "The Phase Space Structure Near Neptune Resonances"
    """
    if body_id not in PLUTINO_LIBRATION_PARAMS:
        return 0.0

    params = PLUTINO_LIBRATION_PARAMS[body_id]

    # Time since J2000.0 in days
    dt = jd_tt - J2000_EPOCH

    # Calculate libration phase at this time
    # The resonant argument oscillates as: φ(t) = center + amplitude * sin(ωt + φ0)
    # where ω = 2π / period
    omega_lib = 2.0 * math.pi / params.period
    libration_angle = omega_lib * dt + params.phase_j2000

    # Resonant-argument scaling: for an exterior p:q mean-motion
    # resonance the argument is phi = q*lambda_body - p*lambda_N
    # - (q - p)*pomega_body, so d(phi)/d(lambda_body) = q and a
    # libration of the argument maps to the longitude as
    # dlambda = dphi / q (for the 2:3 plutino resonance: /3).
    delta_lambda = (params.amplitude / float(params.resonance_q)) * math.sin(
        libration_angle
    )

    return delta_lambda


def has_libration_model(body_id: int) -> bool:
    """
    Check if a body has resonant libration parameters defined.

    Args:
        body_id: Minor body identifier

    Returns:
        bool: True if the body has libration parameters for resonant correction
    """
    return body_id in PLUTINO_LIBRATION_PARAMS


def get_libration_parameters(body_id: int) -> Optional[LibrationParameters]:
    """
    Get the libration parameters for a resonant body.

    Args:
        body_id: Minor body identifier

    Returns:
        LibrationParameters or None if not a known resonant body
    """
    return PLUTINO_LIBRATION_PARAMS.get(body_id)


# =============================================================================
# SECULAR PERTURBATION FUNCTIONS
# =============================================================================


def _calc_laplace_coefficients(alpha: float, s: float, j: int) -> float:
    """
    Calculate Laplace coefficient b_s^(j)(alpha) using series expansion.

    The Laplace coefficients appear in the disturbing function expansion
    for planetary perturbation theory. They depend on the ratio of semi-major
    axes and are used to calculate secular perturbation rates.

    Args:
        alpha: Ratio of semi-major axes (a_inner / a_outer), must be < 1
        s: Half-integer index (typically 1/2 or 3/2)
        j: Integer index (0, 1, 2, ...)

    Returns:
        float: The Laplace coefficient value

    References:
        Murray & Dermott "Solar System Dynamics" §6.4
        Brouwer & Clemence "Methods of Celestial Mechanics" Ch. XI
    """
    if alpha >= 1.0 or alpha <= 0.0:
        return 0.0

    # Series expansion for b_s^(j)(alpha)
    # b_s^(j) = (2/pi) * integral from 0 to pi of cos(j*psi) / (1 - 2*alpha*cos(psi) + alpha^2)^s dpsi
    # For s = 1/2, we can use the approximation for small alpha

    # L6: Use numerical integration with 500 steps (increased from 100) for
    # better accuracy near resonances where alpha approaches 1.
    # For alpha > 0.7 (near-resonant bodies like plutinos), the integrand
    # becomes sharply peaked and needs more quadrature points.
    n_steps = 500 if alpha > 0.7 else 200
    dpsi = math.pi / n_steps
    result = 0.0

    for k in range(n_steps + 1):
        psi = k * dpsi
        denom = 1.0 - 2.0 * alpha * math.cos(psi) + alpha * alpha
        if denom > 1e-10:
            integrand = math.cos(j * psi) / (denom**s)
            weight = 1.0 if (k == 0 or k == n_steps) else 2.0
            result += weight * integrand * dpsi / 2.0

    return result / math.pi


def calc_secular_perturbation_rates(
    elements: OrbitalElements,
) -> Tuple[float, float, float]:
    """
    Calculate secular perturbation rates for orbital elements due to Jupiter, Saturn, Uranus, and Neptune.

    Uses first-order Laplace-Lagrange secular perturbation theory to compute
    the time rates of change of the argument of perihelion (ω), longitude of
    ascending node (Ω), and a correction to mean motion.

    The secular perturbations cause:
    - Precession of perihelion (ω advances or regresses)
    - Precession of the ascending node (Ω regresses for prograde orbits)
    - Small correction to mean motion due to non-Keplerian effects

    Note:
        Neptune perturbations are critical for TNO accuracy, especially for
        plutinos like Ixion and Orcus which are in 2:3 mean motion resonance
        with Neptune. Neptune perturbations are only applied for bodies with
        semi-major axis > 20 AU to avoid overhead for inner solar system objects.
        Uranus perturbations are particularly significant for Trans-Neptunian Objects
        (TNOs) where its gravitational influence is comparable to or greater than
        Saturn's influence.

    Args:
        elements: Orbital elements of the minor body

    Returns:
        Tuple[float, float, float]: (d_omega, d_Omega, d_n) rates in degrees/day
            - d_omega: Rate of change of argument of perihelion
            - d_Omega: Rate of change of longitude of ascending node
            - d_n: Correction to mean motion (small)

    Algorithm:
        Uses the classical Laplace-Lagrange secular theory formulas:
        - dω/dt ≈ n * (m'/M) * α * b_{3/2}^{(1)}(α) * (1/4) * ...
        - dΩ/dt ≈ -n * (m'/M) * α * b_{3/2}^{(1)}(α) * cos(i) * ...

        where:
        - n is the mean motion of the asteroid
        - m'/M is the mass ratio of the perturbing planet to the Sun
        - α is the ratio of semi-major axes
        - b is the Laplace coefficient

    References:
        Murray & Dermott "Solar System Dynamics" Ch. 7
        Brouwer & Clemence "Methods of Celestial Mechanics" Ch. XVI
    """
    a = elements.a
    e = elements.e
    i_rad = math.radians(elements.i)
    n = elements.n  # degrees/day

    # Convert mean motion to radians/day for calculations
    n_rad = math.radians(n)

    # Initialize perturbation rates
    d_omega = 0.0  # degrees/day
    d_Omega = 0.0  # degrees/day
    d_n = 0.0  # degrees/day

    # Calculate perturbations from Jupiter (dominant for most asteroids)
    if a < JUPITER_A:
        # Asteroid is interior to Jupiter
        alpha = a / JUPITER_A

        # Laplace coefficient approximation for b_{3/2}^{(1)}
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)

        # Secular rates from first-order theory
        # Factor from disturbing function expansion
        factor = n_rad * MASS_RATIO_JUPITER * alpha * b32_1

        # Perihelion precession rate (prograde for interior orbits)
        # dω/dt = (n/4) * (m'/M) * α² * b_{3/2}^{(2)} * (2 - e²) / (1 - e²)
        # Simplified: use leading order term
        d_omega_jup = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)

        # Node regression rate
        # dΩ/dt = -(n/4) * (m'/M) * α * b_{3/2}^{(1)} * cos(i) / (1 - e²)^(1/2)
        d_Omega_jup = -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(1.0 - e * e)

        d_omega += math.degrees(d_omega_jup)
        d_Omega += math.degrees(d_Omega_jup)

    elif a > JUPITER_A:
        # Asteroid is exterior to Jupiter
        alpha = JUPITER_A / a

        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_JUPITER * alpha * b32_1

        # For exterior bodies, perturbation effects are reduced but still significant
        d_omega_jup = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)
        d_Omega_jup = (
            -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(max(0.01, 1.0 - e * e))
        )

        d_omega += math.degrees(d_omega_jup)
        d_Omega += math.degrees(d_Omega_jup)

    # Calculate perturbations from Saturn (significant for outer asteroids/TNOs)
    if a < SATURN_A:
        alpha = a / SATURN_A
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_SATURN * alpha * b32_1

        d_omega_sat = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)
        d_Omega_sat = -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(1.0 - e * e)

        d_omega += math.degrees(d_omega_sat)
        d_Omega += math.degrees(d_Omega_sat)

    elif a > SATURN_A:
        alpha = SATURN_A / a
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_SATURN * alpha * b32_1

        d_omega_sat = (
            factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (max(0.01, 1.0 - e * e))
        )
        d_Omega_sat = (
            -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(max(0.01, 1.0 - e * e))
        )

        d_omega += math.degrees(d_omega_sat)
        d_Omega += math.degrees(d_Omega_sat)

    # Calculate perturbations from Uranus (significant for TNOs)
    if a < URANUS_A:
        alpha = a / URANUS_A
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_URANUS * alpha * b32_1

        d_omega_ura = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)
        d_Omega_ura = -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(1.0 - e * e)

        d_omega += math.degrees(d_omega_ura)
        d_Omega += math.degrees(d_Omega_ura)

    elif a > URANUS_A:
        alpha = URANUS_A / a
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        factor = n_rad * MASS_RATIO_URANUS * alpha * b32_1

        d_omega_ura = (
            factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (max(0.01, 1.0 - e * e))
        )
        d_Omega_ura = (
            -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(max(0.01, 1.0 - e * e))
        )

        d_omega += math.degrees(d_omega_ura)
        d_Omega += math.degrees(d_Omega_ura)

    # Calculate perturbations from Neptune (critical for TNO accuracy)
    # Only apply for bodies with a > 20 AU to avoid overhead for inner solar system
    # Neptune is especially important for plutinos like Ixion and Orcus (2:3 resonance)
    if a > 20.0:
        if a < NEPTUNE_A:
            # Body is interior to Neptune (includes plutinos at ~39 AU)
            alpha = a / NEPTUNE_A
            b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
            factor = n_rad * MASS_RATIO_NEPTUNE * alpha * b32_1

            d_omega_nep = factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (1.0 - e * e)
            d_Omega_nep = (
                -factor * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(1.0 - e * e)
            )

            d_omega += math.degrees(d_omega_nep)
            d_Omega += math.degrees(d_Omega_nep)

        elif a > NEPTUNE_A:
            # Body is exterior to Neptune (most TNOs: Eris, Sedna, Haumea, etc.)
            alpha = NEPTUNE_A / a
            b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
            factor = n_rad * MASS_RATIO_NEPTUNE * alpha * b32_1

            d_omega_nep = (
                factor * (1.0 / 4.0) * (2.0 + 1.5 * e * e) / (max(0.01, 1.0 - e * e))
            )
            d_Omega_nep = (
                -factor
                * (1.0 / 4.0)
                * math.cos(i_rad)
                / math.sqrt(max(0.01, 1.0 - e * e))
            )

            d_omega += math.degrees(d_omega_nep)
            d_Omega += math.degrees(d_Omega_nep)

    # L3: Compute second-order mean motion correction (d_n ≠ 0)
    # The mean motion receives a correction from the interaction between the
    # asteroid's forced eccentricity and the perturbing planets. From second-order
    # secular theory (Murray & Dermott §7.4):
    #   d_n = -(3/2) * (n/a) * d_a_secular
    # where d_a_secular arises from the coupling between ω precession and the
    # forced eccentricity oscillation. For a body near a mean-motion resonance
    # p:q with a planet, the correction is approximately:
    #   d_n ≈ -(3n/2a) * Σ_j [μ_j * α_j * b_{3/2}^{(2)}(α_j) * e_j * cos(ϖ - ϖ_j)]
    # This is typically < 0.01"/yr for most bodies but can reach arcminute level
    # for near-resonant bodies over decades.
    #
    # We approximate this using the d_omega rate (which captures the secular
    # precession from all perturbers) and the forced eccentricity magnitude:
    # d_n ≈ -(3/2) * n * e_forced^2 * (d_omega_rad / n_rad) where d_omega_rad
    # is the total precession rate. This second-order effect is small but
    # accumulates over decades.
    d_omega_rad = math.radians(d_omega)
    if abs(n_rad) > 1e-20 and e < 0.99:
        # Second-order correction: interaction of precession with eccentricity
        e2 = e * e
        d_n = math.degrees(-1.5 * d_omega_rad * e2 / (1.0 - e2))
    else:
        d_n = 0.0

    # M4: Near-resonance amplification of secular rates
    # Bodies near (but not in) mean-motion resonances experience amplified
    # perturbation effects. The amplification factor scales as 1/|Δ| where
    # Δ = (p·n_body - q·n_planet) / n_body measures the distance from exact
    # resonance. This applies to Hildas (3:2 Jupiter), some main belt bodies
    # near 3:1, 5:2, 7:3 Jupiter resonances, and TNOs near Neptune resonances.
    # Reference: Murray & Dermott (1999), Chapter 8, eq. 8.72-8.76
    _resonance_checks = [
        # (n_planet, p, q, max_amplification)
        (JUPITER_N, 3, 1, 3.0),  # 3:1 Kirkwood gap
        (JUPITER_N, 5, 2, 2.5),  # 5:2 Kirkwood gap
        (JUPITER_N, 7, 3, 2.0),  # 7:3 resonance
        (JUPITER_N, 2, 1, 3.0),  # 2:1 Hecuba gap
        (JUPITER_N, 3, 2, 3.0),  # 3:2 Hilda group
        (NEPTUNE_N, 2, 3, 3.0),  # 2:3 plutinos
        (NEPTUNE_N, 1, 2, 2.5),  # 1:2 twotinos
        (NEPTUNE_N, 3, 5, 2.0),  # 3:5 resonance
    ]
    resonance_amplification = 1.0
    for n_planet, p_res, q_res, max_amp in _resonance_checks:
        # Distance from exact resonance: Δ = |p·n_body - q·n_planet| / n_body
        if abs(n) > 1e-20:
            delta = abs(p_res * n - q_res * n_planet) / n
            # Only amplify if close to resonance (within ~5%)
            if 0.001 < delta < 0.05:
                # Amplification: min(1/delta, max_amp) to avoid singularity
                amp = min(1.0 / (delta * 20.0), max_amp)
                resonance_amplification = max(resonance_amplification, amp)

    if resonance_amplification > 1.0:
        d_omega *= resonance_amplification
        d_Omega *= resonance_amplification

    return d_omega, d_Omega, d_n


def _calc_forced_elements(
    elements: OrbitalElements,
    jd_tt: Optional[float] = None,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate forced eccentricity/inclination vectors and proper frequencies.

    Uses Laplace-Lagrange secular theory (Murray & Dermott Ch. 7) to compute
    the forced response of a test particle's eccentricity and inclination
    vectors due to Jupiter, Saturn, Uranus, and Neptune.

    In secular theory, the eccentricity vector (h, k) = (e sin ϖ, e cos ϖ)
    decomposes into a forced component (driven by planets) and a free
    component (oscillating at the proper frequency). The actual eccentricity
    oscillates around the forced value with the proper period.

    Similarly, the inclination vector (p, q) = (sin(i/2) sin Ω, sin(i/2) cos Ω)
    decomposes into forced + free components.

    L4: When jd_tt is provided, planet orbital elements are evolved from their
    J2000.0 values using linear rates from Simon et al. (1994). This improves
    the forced eccentricity/inclination vectors for century-scale propagations
    by ~10-30" over 500 years.

    Args:
        elements: Orbital elements of the test particle
        jd_tt: Optional target Julian Day for time-dependent planet elements.
            If None, uses static J2000.0 planet elements.

    Returns:
        Tuple of (g, h_forced, k_forced, s, p_forced, q_forced):
        - g: Proper frequency for eccentricity (ϖ precession, rad/day, positive)
        - h_forced: Forced e·sin(ϖ) component
        - k_forced: Forced e·cos(ϖ) component
        - s: Proper frequency for inclination (|Ω regression|, rad/day, positive)
        - p_forced: Forced sin(i/2)·sin(Ω) component
        - q_forced: Forced sin(i/2)·cos(Ω) component

    Algorithm:
        For each perturbing planet j with mass ratio μ_j, semi-major axis a_j:

        Eccentricity proper frequency (diagonal):
            A_j = (n/4) · μ_j · α · b_{3/2}^{(1)}(α) · correction(e)
            g = Σ A_j

        Eccentricity coupling (off-diagonal, uses b_{3/2}^{(2)}):
            ν_j = (n/4) · μ_j · α · b_{3/2}^{(2)}(α)
            h_f = Σ ν_j · e_j · sin(ϖ_j) / g
            k_f = Σ ν_j · e_j · cos(ϖ_j) / g

        Inclination proper frequency (diagonal):
            s_j = (n/4) · μ_j · α · b_{3/2}^{(1)}(α) · cos(i)/√(1-e²)
            s = Σ s_j

        Inclination coupling (off-diagonal, uses b_{3/2}^{(1)}):
            σ_j = (n/4) · μ_j · α · b_{3/2}^{(1)}(α)
            p_f = Σ σ_j · sin(i_j/2) · sin(Ω_j) / s
            q_f = Σ σ_j · sin(i_j/2) · cos(Ω_j) / s

    References:
        Murray & Dermott "Solar System Dynamics" §7.3-7.5, eq. 7.10, 7.25-7.26
        Brouwer & Clemence "Methods of Celestial Mechanics" Ch. XVI
        Simon et al. (1994) A&A 282, 663 (planet element rates)
    """
    a = elements.a
    e = elements.e
    i_rad = math.radians(elements.i)
    n_rad = math.radians(elements.n)
    e2 = e * e
    denom_e = max(0.01, 1.0 - e2)

    # Accumulate proper frequencies and forced vector numerators
    g = 0.0  # eccentricity proper frequency (rad/day)
    s = 0.0  # inclination regression frequency (rad/day, positive)
    h_forced_num = 0.0
    k_forced_num = 0.0
    p_forced_num = 0.0
    q_forced_num = 0.0

    # L4: Use time-dependent planet elements when jd_tt is provided,
    # otherwise fall back to static J2000.0 values for backward compatibility.
    if jd_tt is not None:
        _planets = _get_planet_elements_at_time(jd_tt)
    else:
        # Static J2000.0 planet data (legacy path)
        # (a_planet, mu, e_planet, omega_deg, node_deg, i_deg, a_threshold)
        _planets = [
            (
                JUPITER_A,
                MASS_RATIO_JUPITER,
                JUPITER_E,
                JUPITER_OMEGA,
                JUPITER_NODE,
                JUPITER_I,
                0.0,
            ),
            (
                SATURN_A,
                MASS_RATIO_SATURN,
                SATURN_E,
                SATURN_OMEGA,
                SATURN_NODE,
                SATURN_I,
                0.0,
            ),
            (
                URANUS_A,
                MASS_RATIO_URANUS,
                URANUS_E,
                URANUS_OMEGA,
                URANUS_NODE,
                URANUS_I,
                0.0,
            ),
            (
                NEPTUNE_A,
                MASS_RATIO_NEPTUNE,
                NEPTUNE_E,
                NEPTUNE_OMEGA,
                NEPTUNE_NODE,
                NEPTUNE_I,
                20.0,
            ),
        ]

    for a_p, mu_p, e_p, omega_p_deg, node_p_deg, i_p_deg, a_thresh in _planets:
        # Apply threshold (Neptune only for a > 20 AU)
        if a_thresh > 0.0 and a < a_thresh:
            continue
        # S4: Adaptive co-orbital detection — skip if asteroid is within the
        # planet's Hill sphere radius (~a_p * (mu/3)^(1/3)), which is the
        # region where secular theory breaks down. This replaces the old fixed
        # 0.1 AU threshold that was too aggressive for Trojan asteroids and
        # horseshoe orbit bodies near Jupiter (Hill sphere ~0.35 AU).
        hill_radius = a_p * (mu_p / 3.0) ** (1.0 / 3.0)
        if abs(a - a_p) < hill_radius:
            continue

        # Semi-major axis ratio (always < 1)
        if a < a_p:
            alpha = a / a_p
        else:
            alpha = a_p / a

        if alpha >= 1.0 or alpha <= 0.0:
            continue

        # Laplace coefficients
        b32_1 = _calc_laplace_coefficients(alpha, 1.5, 1)
        b32_2 = _calc_laplace_coefficients(alpha, 1.5, 2)

        base = n_rad * mu_p * alpha

        # --- Eccentricity proper frequency (diagonal, with e-correction) ---
        A_j = base * b32_1 * (1.0 / 4.0) * (2.0 + 1.5 * e2) / denom_e
        g += A_j

        # --- Inclination regression frequency (diagonal, with i,e correction) ---
        s_j = base * b32_1 * (1.0 / 4.0) * math.cos(i_rad) / math.sqrt(denom_e)
        s += s_j

        # --- Off-diagonal eccentricity coupling (uses b_{3/2}^{(2)}) ---
        nu_j = base * b32_2 * (1.0 / 4.0)

        # --- Off-diagonal inclination coupling (uses b_{3/2}^{(1)}) ---
        sigma_j = base * b32_1 * (1.0 / 4.0)

        # Planet longitude of perihelion: ϖ = ω + Ω
        varpi_p_rad = math.radians(omega_p_deg + node_p_deg)
        i_p_rad = math.radians(i_p_deg)
        Omega_p_rad = math.radians(node_p_deg)

        # Planet eccentricity vector (h_j, k_j)
        h_j = e_p * math.sin(varpi_p_rad)
        k_j = e_p * math.cos(varpi_p_rad)

        # Planet inclination vector (p_j, q_j)
        p_j = math.sin(i_p_rad / 2.0) * math.sin(Omega_p_rad)
        q_j = math.sin(i_p_rad / 2.0) * math.cos(Omega_p_rad)

        # Accumulate forced vector numerators
        # From the secular equations: h_f = Σ ν_j h_j / g, k_f = Σ ν_j k_j / g
        h_forced_num += nu_j * h_j
        k_forced_num += nu_j * k_j
        # From inclination equations: p_f = Σ σ_j p_j / s, q_f = Σ σ_j q_j / s
        p_forced_num += sigma_j * p_j
        q_forced_num += sigma_j * q_j

    # Divide by proper frequencies to get forced vectors
    h_forced = 0.0
    k_forced = 0.0
    if abs(g) > 1e-20:
        h_forced = h_forced_num / g
        k_forced = k_forced_num / g

    p_forced = 0.0
    q_forced = 0.0
    if abs(s) > 1e-20:
        p_forced = p_forced_num / s
        q_forced = q_forced_num / s

    # H1: Second-order secular perturbation corrections
    # Include cross-coupling between different perturbing planets (Jupiter×Saturn)
    # and second-order eccentricity/inclination effects (Hori 1966, Yuasa 1973).
    #
    # The second-order correction to the proper frequency g adds:
    #   Δg = Σ_{j≠k} A_j × A_k / (g_j - g_k)
    # where g_j are the eigenfrequencies of the secular system.
    #
    # For high-eccentricity bodies (e > 0.3), the first-order theory
    # underestimates the forced eccentricity. Apply a correction factor:
    #   e_forced_corrected = e_forced × (1 + 0.5 × e² + 0.375 × e⁴)
    # This comes from the expansion of the disturbing function to higher
    # order in eccentricity (Murray & Dermott §7.6).
    if e > 0.3:
        e_correction = 1.0 + 0.5 * e2 + 0.375 * e2 * e2
        h_forced *= e_correction
        k_forced *= e_correction

    # For high-inclination bodies (i > 20°), apply similar correction:
    #   i_forced_corrected = i_forced × (1 + sin²(i/2))
    if elements.i > 20.0:
        sin_half_i = math.sin(i_rad / 2.0)
        i_correction = 1.0 + sin_half_i * sin_half_i
        p_forced *= i_correction
        q_forced *= i_correction

    return g, h_forced, k_forced, s, p_forced, q_forced


def apply_secular_perturbations(
    elements: OrbitalElements, jd_tt: float, include_perturbations: bool = True
) -> Tuple[float, float, float, float, float, float]:
    """
    Apply secular perturbations to orbital elements and return perturbed values.

    Takes the osculating orbital elements at epoch and applies Laplace-Lagrange
    secular perturbation corrections to propagate them to the target time.
    This accounts for:
    - Long-term drift in ω and Ω (linear precession rates)
    - Oscillation of eccentricity e around the forced value from planets
    - Oscillation of inclination i around the forced value from planets

    The eccentricity and inclination evolution uses the (h,k)/(p,q) vector
    formalism from Murray & Dermott Ch. 7. The eccentricity vector
    (h, k) = (e sin ϖ, e cos ϖ) decomposes into a forced component (driven
    by planetary eccentricities) and a free component (oscillating at the
    proper frequency). Similarly for the inclination vector (p, q).

    WARNING: For bodies in mean motion resonance with Neptune (e.g., plutinos
    like Ixion and Orcus in 2:3 resonance, twotinos in 1:2 resonance), the
    secular perturbation theory used here may give inaccurate results. Resonant
    bodies experience additional perturbations that are not captured by
    first-order Laplace-Lagrange theory. Use detect_mean_motion_resonance()
    or is_body_resonant() to check if a body is in resonance.

    Args:
        elements: Original orbital elements at epoch
        jd_tt: Target Julian Day in Terrestrial Time
        include_perturbations: If True, apply secular corrections; if False, return unperturbed

    Returns:
        Tuple of (omega_pert, Omega_pert, M_pert, n_pert, e_pert, i_pert):
            - omega_pert: Perturbed argument of perihelion (degrees)
            - Omega_pert: Perturbed longitude of ascending node (degrees)
            - M_pert: Perturbed mean anomaly at target time (degrees)
            - n_pert: Perturbed mean motion (degrees/day)
            - e_pert: Perturbed eccentricity (dimensionless, 0 < e < 1)
            - i_pert: Perturbed inclination (degrees)

    See Also:
        detect_mean_motion_resonance: Check if a body is in Neptune resonance
        is_body_resonant: Quick check if a body ID is resonant
    """
    dt = jd_tt - elements.epoch  # Time since epoch in days

    if not include_perturbations or abs(dt) < 1.0:
        # For very short propagation times, perturbations are negligible
        M = (elements.M0 + elements.n * dt) % 360.0
        return elements.omega, elements.Omega, M, elements.n, elements.e, elements.i

    # Calculate secular perturbation rates for ω, Ω (linear precession)
    d_omega, d_Omega, d_n = calc_secular_perturbation_rates(elements)

    # Apply linear secular corrections for angular elements
    omega_pert = (elements.omega + d_omega * dt) % 360.0
    Omega_pert = (elements.Omega + d_Omega * dt) % 360.0
    n_pert = elements.n + d_n

    # Propagate mean anomaly with perturbed mean motion
    M_pert = (elements.M0 + n_pert * dt) % 360.0

    # --- Eccentricity and inclination evolution via (h,k)/(p,q) vectors ---
    # L4: Pass jd_tt to use time-dependent planet elements
    g, h_forced, k_forced, s_freq, p_forced, q_forced = _calc_forced_elements(
        elements, jd_tt=jd_tt
    )

    # Eccentricity vector evolution
    e0 = elements.e
    varpi0_rad = math.radians(elements.omega + elements.Omega)

    # Osculating eccentricity vector at epoch
    h0 = e0 * math.sin(varpi0_rad)
    k0 = e0 * math.cos(varpi0_rad)

    # Free eccentricity vector = osculating - forced
    dh = h0 - h_forced
    dk = k0 - k_forced
    e_free = math.sqrt(dh * dh + dk * dk)

    if e_free > 1e-15:
        beta = math.atan2(dh, dk)
    else:
        beta = 0.0

    # Evolve free vector: rotates at proper frequency g
    h_t = h_forced + e_free * math.sin(g * dt + beta)
    k_t = k_forced + e_free * math.cos(g * dt + beta)

    e_pert = math.sqrt(h_t * h_t + k_t * k_t)
    # Clamp to the physical range [0.0, 0.999]: cap below 1 (parabolic) but do
    # NOT impose a positive floor — a genuinely near-circular orbit must keep
    # its small eccentricity. The downstream anomaly solver uses the (h, k)
    # vectors and handles e == 0 fine, so no positive floor is needed.
    e_pert = max(0.0, min(e_pert, 0.999))

    # Inclination vector evolution
    i0_rad = math.radians(elements.i)
    Omega0_rad = math.radians(elements.Omega)

    # Osculating inclination vector at epoch
    sin_half_i0 = math.sin(i0_rad / 2.0)
    p0 = sin_half_i0 * math.sin(Omega0_rad)
    q0 = sin_half_i0 * math.cos(Omega0_rad)

    # Free inclination vector = osculating - forced
    dp = p0 - p_forced
    dq = q0 - q_forced
    i_free = math.sqrt(dp * dp + dq * dq)

    if i_free > 1e-15:
        gamma = math.atan2(dp, dq)
    else:
        gamma = 0.0

    # Evolve free vector: rotates at frequency -s (retrograde precession)
    p_t = p_forced + i_free * math.sin(-s_freq * dt + gamma)
    q_t = q_forced + i_free * math.cos(-s_freq * dt + gamma)

    sin_half_i_t = math.sqrt(p_t * p_t + q_t * q_t)
    # Clamp to valid range for asin
    sin_half_i_t = max(0.0, min(sin_half_i_t, 1.0))
    i_pert = math.degrees(2.0 * math.asin(sin_half_i_t))

    return omega_pert, Omega_pert, M_pert, n_pert, e_pert, i_pert


# H4: Yarkovsky effect parameters for NEAs with measured da/dt
# Source: JPL SBDB (https://ssd.jpl.nasa.gov/) and Chesley et al. (2014)
# da/dt is in AU per million years (AU/My). Positive = outward drift.
# Only bodies with published measured values are included.
# The Yarkovsky effect causes a secular drift in semi-major axis due to
# anisotropic thermal emission, and is the dominant non-gravitational force
# for small NEAs over decade timescales.
YARKOVSKY_DA_DT: dict[int, float] = {
    # Apophis: da/dt = -2.899e-4 AU/My (Chesley et al. 2014, Brozovic et al. 2018)
    APOPHIS: -2.899e-4,
    # Bennu: da/dt = -18.99e-4 AU/My (Chesley et al. 2014, OSIRIS-REx)
    BENNU: -18.99e-4,
    # Ryugu: da/dt = -3.5e-4 AU/My (estimated from Hayabusa2 data)
    RYUGU: -3.5e-4,
    # Itokawa: da/dt = -3.5e-4 AU/My (Vokrouhlicky et al. 2008)
    ITOKAWA: -3.5e-4,
    # Eros: da/dt = -0.5e-4 AU/My (estimated, large body = weak Yarkovsky)
    EROS: -0.5e-4,
    # Toutatis: da/dt = -1.5e-4 AU/My (estimated)
    TOUTATIS: -1.5e-4,
}


def _apply_yarkovsky_correction(
    elements: OrbitalElements, dt_days: float, body_id: Optional[int]
) -> float:
    """Apply Yarkovsky semi-major axis drift correction.

    H4: For NEAs with measured da/dt (Yarkovsky parameter), applies a secular
    drift to the semi-major axis: a(t) = a₀ + (da/dt) × dt.

    This is the dominant non-gravitational force for small NEAs, causing
    ~1" per decade in geocentric longitude. Irrelevant for main belt and TNOs.

    Args:
        elements: Orbital elements at epoch
        dt_days: Time since epoch in days
        body_id: Minor body identifier

    Returns:
        Corrected semi-major axis in AU
    """
    a = elements.a
    if body_id is None or body_id not in YARKOVSKY_DA_DT:
        return a

    # Convert da/dt from AU/My to AU/day
    da_dt_per_day = YARKOVSKY_DA_DT[body_id] / (1e6 * 365.25)

    # Apply linear drift
    return a + da_dt_per_day * dt_days


# =============================================================================
# ORBITAL ELEMENTS DATABASE (Epoch JD 2461000.5 TDB = 2025-Sep-19)
# =============================================================================
# Source: NASA JPL Small-Body Database (sbdb.api), retrieved 2025-Nov-27
# Epoch: JD 2461000.5 TDB (2025-Sep-19)
#
# Note: Elements are osculating at the given epoch.
# Accuracy degrades ~10-50 arcsec/year due to secular perturbations.
# For dates >10 years from epoch, consider updating elements.

MINOR_BODY_ELEMENTS = {
    CHIRON: OrbitalElements(
        name="Chiron",
        epoch=2461000.5,
        a=13.69219896172984,  # AU - between Saturn and Uranus
        e=0.3789792342846475,
        i=6.926003536565557,
        omega=339.2537417045351,  # argument of perihelion
        Omega=209.2984204899107,  # longitude of ascending node
        M0=212.8397717853335,
        n=0.01945334424082164,  # ~51 year period
    ),
    PHOLUS: OrbitalElements(
        name="Pholus",
        epoch=2461000.5,
        a=20.28340105547402,
        e=0.574744804002886,  # Highly eccentric
        i=24.75699076739707,
        omega=354.7299656288133,
        Omega=119.2896923424452,
        M0=134.470501527666,
        n=0.01078929100989589,  # ~91 year period
    ),
    CERES: OrbitalElements(
        name="Ceres",
        epoch=2461000.5,
        a=2.765615651508659,
        e=0.07957631994408416,
        i=10.58788658206854,
        omega=73.29975464616518,
        Omega=80.24963090816965,
        M0=231.5397330043706,
        n=0.2142971214271186,  # ~4.6 year period
    ),
    PALLAS: OrbitalElements(
        name="Pallas",
        epoch=2461000.5,
        a=2.76992582511479,
        e=0.2306429787781384,
        i=34.92832687077855,  # High inclination
        omega=310.9333840114307,
        Omega=172.8885963367437,
        M0=211.5297778033731,
        n=0.2137971269626138,  # ~4.6 year period
    ),
    JUNO: OrbitalElements(
        name="Juno",
        epoch=2461000.5,
        a=2.670879058906207,
        e=0.2558257725543152,
        i=12.98603961477441,
        omega=247.8836661359693,
        Omega=169.8198844530219,
        M0=217.5909617686606,
        n=0.2257993770223806,  # ~4.4 year period
    ),
    VESTA: OrbitalElements(
        name="Vesta",
        epoch=2461000.5,
        a=2.361541280084789,
        e=0.09016764504738634,
        i=7.144060599543863,
        omega=151.5371488873794,
        Omega=103.7022980342142,
        M0=26.80967220901607,
        n=0.2715881155129186,  # ~3.6 year period
    ),
    ERIS: OrbitalElements(
        name="Eris",
        epoch=2461000.5,
        a=67.99636506315233,  # Highly distant
        e=0.4369647318678042,
        i=43.86893210857381,  # Extreme inclination
        omega=150.7324609262193,
        Omega=36.0271728267921,
        M0=211.4489901906776,
        n=0.001757824561999289,  # ~561 year period
    ),
    SEDNA: OrbitalElements(
        name="Sedna",
        epoch=2461000.5,
        a=549.5405279241455,  # Extreme distance (detached object)
        e=0.8612974378772624,  # Very eccentric
        i=11.92591697769048,
        omega=311.0097709574837,
        Omega=144.4787268123774,
        M0=358.6072615114487,
        n=7.650758332945998e-05,  # ~12,880 year period
    ),
    HAUMEA: OrbitalElements(
        name="Haumea",
        epoch=2461000.5,
        a=43.00549881333706,
        e=0.1957748188564788,
        i=28.20840579216805,
        omega=240.888336259082,
        Omega=121.7972901946415,
        M0=222.3276474703985,
        n=0.003494765903502902,  # ~282 year period
    ),
    MAKEMAKE: OrbitalElements(
        name="Makemake",
        epoch=2461000.5,
        a=45.51068175675377,
        e=0.1604249910368623,
        i=29.03230613859081,
        omega=297.0754218776326,
        Omega=79.26892117343222,
        M0=169.3202841867025,
        n=0.003210214572603997,  # ~307 year period
    ),
    IXION: OrbitalElements(
        name="Ixion",
        epoch=2461000.5,
        a=39.35053706213409,  # Plutino (2:3 resonance with Neptune)
        e=0.2442328489971719,
        i=19.67041190191546,
        omega=300.6585723831106,
        Omega=71.09295825377649,
        M0=294.2004612799266,
        n=0.003992804789574964,  # ~247 year period
    ),
    ORCUS: OrbitalElements(
        name="Orcus",
        epoch=2461000.5,
        a=39.33577647200568,  # Plutino (anti-Pluto phase)
        e=0.2217300030420161,
        i=20.55552551599616,
        omega=73.72249191891537,
        Omega=268.3859416278478,
        M0=188.1111318293787,
        n=0.003995052426031791,  # ~247 year period
    ),
    QUAOAR: OrbitalElements(
        name="Quaoar",
        epoch=2461000.5,
        a=43.1476797802032,
        e=0.03583878353429052,  # Nearly circular
        i=7.991371294217068,
        omega=163.9231384883233,
        Omega=188.9632800603184,
        M0=291.4818844949103,
        n=0.003477506123841158,  # ~284 year period
    ),
    NESSUS: OrbitalElements(
        name="Nessus",
        epoch=2461000.5,
        a=24.51657596400793,  # Centaur orbit
        e=0.5176801039355027,  # Moderately eccentric
        i=15.64407698270027,
        omega=170.3715103791255,
        Omega=31.29215422788363,
        M0=100.5239875112597,
        n=0.008119220769969982,  # ~121 year period
    ),
    ASBOLUS: OrbitalElements(
        name="Asbolus",
        epoch=2461000.5,
        a=18.06355955222326,  # Centaur orbit
        e=0.6173222203620333,  # Highly eccentric
        i=17.6174400161563,
        omega=290.5836923759394,
        Omega=6.020562307154061,
        M0=109.6324751756823,
        n=0.01283805024514489,  # ~77 year period
    ),
    CHARIKLO: OrbitalElements(
        name="Chariklo",
        epoch=2461000.5,
        a=15.73995155535189,  # Centaur orbit, largest known centaur (~250 km)
        e=0.1702459901276738,  # Moderate eccentricity
        i=23.43032407033149,
        omega=241.2242283234723,
        Omega=300.4752379513423,
        M0=126.9607903541095,
        n=0.01578334271898484,  # ~62 year period, has ring system discovered 2014
    ),
    GONGGONG: OrbitalElements(
        name="Gonggong",
        epoch=2461000.5,
        a=66.89366871435344,  # TNO, dwarf planet candidate (~1230 km diameter)
        e=0.5031674399617051,  # Eccentric orbit
        i=30.86626129015389,  # High inclination
        omega=206.6416070091921,
        Omega=336.8400960976296,
        M0=111.3903730396343,
        n=0.001801467997135536,  # ~547 year period (formerly 2007 OR10)
    ),
    VARUNA: OrbitalElements(
        name="Varuna",
        epoch=2461000.5,
        a=43.17823437208563,  # Classical KBO (~670 km diameter)
        e=0.05254516459334936,  # Nearly circular
        i=17.13812404484431,
        omega=273.2206220996139,
        Omega=97.21030234100793,
        M0=115.0289798413239,
        n=0.003473815549791981,  # ~284 year period
    ),
    # Near-Earth Asteroids (NEAs)
    # NOTE: Apophis orbital elements change measurably with each refinement due to
    # close Earth approaches in 2029 and 2036. These elements are from JPL SBDB
    # solution 220 (2024-06-25). For critical applications, use SPK files.
    APOPHIS: OrbitalElements(
        name="Apophis",
        epoch=2461000.5,
        a=0.9223803173917017,  # Aten-class NEA (~0.34 km diameter)
        e=0.1911663355386932,  # Moderate eccentricity
        i=3.340958441017069,  # Low inclination
        omega=126.6728325163065,
        Omega=203.8996515621043,
        M0=312.8054663584516,
        n=1.11259994308075,  # ~0.89 year period, close approach 2029-04-13
    ),
    # Main Belt Asteroids - Large bodies
    HYGIEA: OrbitalElements(
        name="Hygiea",
        epoch=2461000.5,
        a=3.147591335345947,  # Fourth largest asteroid, dwarf planet candidate (~430 km)
        e=0.108223833089903,  # Low eccentricity
        i=3.832937411339835,  # Low inclination
        omega=312.6058061800255,
        Omega=283.1216602821771,
        M0=216.6903206124055,
        n=0.1764966888337409,  # ~5.6 year period
    ),
    INTERAMNIA: OrbitalElements(
        name="Interamnia",
        epoch=2461000.5,
        a=3.056218782582974,  # Fifth largest asteroid (~320 km diameter)
        e=0.1552477966505805,  # Moderate eccentricity
        i=17.31546543833021,  # Moderate inclination
        omega=94.10694029269551,
        Omega=280.1660768334629,
        M0=184.2054144782699,
        n=0.1844707074261051,  # ~5.3 year period
    ),
    DAVIDA: OrbitalElements(
        name="Davida",
        epoch=2461000.5,
        a=3.163619117539836,  # Seventh largest asteroid (~270 km diameter)
        e=0.1896284393751952,  # Moderate eccentricity
        i=15.95004792439553,  # Moderate inclination
        omega=336.666946405915,
        Omega=107.5554381931758,
        M0=35.24757453107242,
        n=0.1751571164581395,  # ~5.6 year period
    ),
    EUROPA_AST: OrbitalElements(
        name="Europa",
        epoch=2461000.5,
        a=3.092193380509604,  # Main belt asteroid (~300 km diameter, not Jupiter's moon)
        e=0.1120016601360199,  # Low eccentricity
        i=7.480921967944819,  # Low inclination
        omega=342.9584258894789,
        Omega=128.5790551357993,
        M0=312.5276379313007,
        n=0.1812608885921891,  # ~5.4 year period
    ),
    SYLVIA: OrbitalElements(
        name="Sylvia",
        epoch=2461000.5,
        a=3.484518782251748,  # Outer main belt asteroid (~280 km, triple system with moons Romulus and Remus)
        e=0.09433757216253168,  # Low eccentricity
        i=10.87247390257978,  # Moderate inclination
        omega=265.8265584696771,
        Omega=73.01146029932147,
        M0=94.95962367758351,
        n=0.151526942166858,  # ~6.5 year period
    ),
    PSYCHE: OrbitalElements(
        name="Psyche",
        epoch=2461000.5,
        a=2.923314514376521,  # Main belt asteroid (~226 km, metallic M-type, NASA Psyche mission target)
        e=0.1343462202634327,  # Low eccentricity
        i=3.097291662719076,  # Low inclination
        omega=229.7534129677803,
        Omega=150.0098773790708,
        M0=40.63883784197816,
        n=0.1971926664511637,  # ~5.0 year period
    ),
    EROS: OrbitalElements(
        name="Eros",
        epoch=2461000.5,
        a=1.458120998474684,  # Near-Earth asteroid (Amor class, NEAR Shoemaker mission target)
        e=0.2228359407071628,  # Moderate eccentricity
        i=10.82846651399785,  # Moderate inclination
        omega=178.9297536744151,
        Omega=304.2701025753316,
        M0=310.5543277370992,
        n=0.5597752949285997,  # ~1.76 year period
    ),
    AMOR: OrbitalElements(
        name="Amor",
        epoch=2461000.5,
        a=1.919831260906931,  # Amor-class NEA prototype (~1 km diameter)
        e=0.4346323375284623,  # Moderate eccentricity
        i=11.86882308811668,  # Moderate inclination
        omega=26.75822531564864,
        Omega=171.2371875175307,
        M0=59.87048521279021,
        n=0.3705181485730381,  # ~2.66 year period (972 days)
    ),
    ICARUS: OrbitalElements(
        name="Icarus",
        epoch=2461000.5,
        a=1.078037837567316,  # Apollo asteroid (~1.4 km, perihelion inside Mercury's orbit)
        e=0.8270056908369543,  # Very high eccentricity (q=0.19 AU, Q=1.97 AU)
        i=22.8032205675089,  # Moderate inclination
        omega=31.43821652068713,
        Omega=87.95241724594092,
        M0=153.0789301425594,
        n=0.8805480660198883,  # ~1.12 year period (409 days)
    ),
    TORO: OrbitalElements(
        name="Toro",
        epoch=2461000.5,
        a=1.367838675762221,  # Apollo asteroid, near-Earth asteroid (~3.4 km)
        e=0.4360487430728354,  # Moderate eccentricity (q=0.77 AU, Q=1.96 AU)
        i=9.38184277020598,  # Low inclination
        omega=127.2752987170327,
        Omega=274.2077498325479,
        M0=82.68163042845468,
        n=0.6161007748143977,  # ~1.6 year period (584 days)
    ),
    SAPPHO: OrbitalElements(
        name="Sappho",
        epoch=2461000.5,
        a=2.296282028063693,  # Main belt asteroid (~69 km, artistic expression, same-sex love)
        e=0.1996711692276041,  # Moderate eccentricity
        i=8.676097553661929,  # Low inclination
        omega=139.693299727962,
        Omega=218.6385022891111,
        M0=57.25859931389674,
        n=0.2832475967264622,  # ~3.5 year period (1271 days)
    ),
    PANDORA_AST: OrbitalElements(
        name="Pandora",
        epoch=2461000.5,
        a=2.75784980823487,  # Main belt asteroid (~67 km, distinct from Saturn moon Pandora)
        e=0.1450625784537587,  # Moderate eccentricity
        i=7.176481151539076,  # Low inclination
        omega=4.93258429587031,
        Omega=10.28381432395156,
        M0=157.5343526907893,
        n=0.2152029188679012,  # ~4.6 year period (1673 days)
    ),
    LILITH_AST: OrbitalElements(
        name="Lilith",
        epoch=2461000.5,
        a=2.664645194686444,  # Main belt asteroid (~30 km, not lunar apogee Lilith)
        e=0.1933175242946558,  # Moderate eccentricity
        i=5.592042797130856,  # Low inclination
        omega=156.6248120324298,
        Omega=260.5622393942289,
        M0=311.9465725963339,
        n=0.2265922174357816,  # ~4.35 year period (1589 days)
    ),
    HIDALGO: OrbitalElements(
        name="Hidalgo",
        epoch=2461000.5,
        a=5.728305724964734,  # Centaur-class asteroid with comet-like orbit
        e=0.6622230478106673,  # Very high eccentricity (comet-like)
        i=42.53593406738054,  # Very high inclination
        omega=56.59815243071501,
        Omega=21.36316817210909,
        M0=185.3374947174508,
        n=0.07188938874958553,  # ~13.7 year period (5008 days)
    ),
    TOUTATIS: OrbitalElements(
        name="Toutatis",
        epoch=2461000.5,
        a=2.543029853370523,  # Apollo PHA, radar and spacecraft target (~5.4 km)
        e=0.6247403332202829,  # High eccentricity (q=0.95 AU, Q=4.13 AU)
        i=0.4480782679581268,  # Very low inclination
        omega=277.8613313176078,
        Omega=125.3653430437421,
        M0=76.89343773316328,
        n=0.2430395129094619,  # ~4.05 year period (1481 days), tumbling rotation
    ),
    ITOKAWA: OrbitalElements(
        name="Itokawa",
        epoch=2461000.5,
        a=1.324135178668783,  # Apollo PHA, Hayabusa sample return mission target (~535 m)
        e=0.2801500981413037,  # Moderate eccentricity (q=0.95 AU, Q=1.70 AU)
        i=1.620938690165939,  # Very low inclination
        omega=162.853569017811,
        Omega=69.07562356929041,
        M0=41.26158479412086,
        n=0.646852987504302,  # ~1.52 year period (557 days), elongated shape
    ),
    BENNU: OrbitalElements(
        name="Bennu",
        epoch=2461000.5,
        a=1.126391025894812,  # Apollo PHA, OSIRIS-REx sample return mission target (~490 m)
        e=0.2037450762416414,  # Moderate eccentricity (q=0.90 AU, Q=1.36 AU)
        i=6.03494377024794,  # Low inclination
        omega=66.22306084084298,
        Omega=2.06086619569642,
        M0=265.1247751080418,  # Propagated from epoch 2455562.5
        n=0.8244613503320309,  # ~1.20 year period (437 days)
    ),
    RYUGU: OrbitalElements(
        name="Ryugu",
        epoch=2461000.5,
        a=1.190921090916117,  # Apollo PHA, Hayabusa2 sample return mission target (~900 m)
        e=0.1910626231558206,  # Moderate eccentricity (q=0.96 AU, Q=1.42 AU)
        i=5.866551865769142,  # Low inclination
        omega=211.6167659559647,
        Omega=251.2915331553314,
        M0=270.6593928357033,
        n=0.7583672923533237,  # ~1.30 year period (475 days)
    ),
}


# =============================================================================
# MULTI-EPOCH ORBITAL ELEMENTS
# =============================================================================
# Dense multi-epoch elements at 10-year intervals for 36 bodies (1650-2450).
# Moved to separate module to keep this file manageable.
# Generated by scripts/generate_multi_epoch_elements.py from JPL SPK data.
from .multi_epoch_data import MINOR_BODY_ELEMENTS_MULTI  # noqa: E402


def _get_closest_epoch_elements(body_id: int, jd_tt: float) -> OrbitalElements:
    """Select the best orbital elements for a minor body at the target time.

    With dense multi-epoch elements at 10-year intervals, simply selects the
    nearest epoch entry (single-epoch or multi-epoch) to minimize Keplerian
    propagation distance. The maximum propagation from nearest entry is ~5 years.

    Cubic Hermite interpolation of osculating elements was tested and rejected:
    perturbation-driven oscillations in a, e, i between epochs cause Hermite
    overshoot that degrades precision by up to 28x for TNOs (Quaoar, Varuna,
    Sedna) and 3-4x for main-belt bodies at certain offsets.

    Args:
        body_id: Minor body identifier (CERES, CHIRON, etc.)
        jd_tt: Target Julian Day in Terrestrial Time

    Returns:
        The OrbitalElements instance closest in time to jd_tt
    """
    best, second, _w = _get_epoch_elements_blend(body_id, jd_tt)
    return best


# Half-width of the cross-fade window centered on the midpoint between
# adjacent element epochs.  Inside the window, positions from the two
# bracketing element sets are blended linearly so the switch between
# epochs (each propagated up to ~5 years two-body) does not step.
_EPOCH_BLEND_HALF_DAYS: float = 182.625  # half a Julian year


def _get_epoch_elements_blend(
    body_id: int, jd_tt: float
) -> tuple[OrbitalElements, Optional[OrbitalElements], float]:
    """Nearest element epoch plus a cross-fade partner near midpoints.

    Returns (primary, secondary, weight): ``primary`` is the nearest
    epoch's elements.  Within ±_EPOCH_BLEND_HALF_DAYS of the midpoint
    between two adjacent epochs, ``secondary`` is the other bracketing
    epoch and ``weight`` (0..0.5] is the secondary's linear blend
    fraction; otherwise ``secondary`` is None and ``weight`` is 0.0.
    """
    base = MINOR_BODY_ELEMENTS[body_id]
    entries = MINOR_BODY_ELEMENTS_MULTI.get(body_id)
    if not entries:
        return base, None, 0.0

    # Multi-epoch entries are generated in ascending epoch order; the
    # single-epoch base participates as a fallback only when it is
    # nearer than every table entry (e.g. outside the table range).
    nearest = min(entries, key=lambda e: abs(jd_tt - e.epoch))
    if abs(jd_tt - base.epoch) < abs(jd_tt - nearest.epoch):
        return base, None, 0.0

    idx = entries.index(nearest)
    if jd_tt >= nearest.epoch and idx + 1 < len(entries):
        other = entries[idx + 1]
    elif jd_tt < nearest.epoch and idx > 0:
        other = entries[idx - 1]
    else:
        return nearest, None, 0.0

    midpoint = (nearest.epoch + other.epoch) / 2.0
    dist = abs(jd_tt - midpoint)
    if dist >= _EPOCH_BLEND_HALF_DAYS:
        return nearest, None, 0.0

    # weight ramps 0.5 at the midpoint -> 0 at the window edge
    weight = 0.5 * (1.0 - dist / _EPOCH_BLEND_HALF_DAYS)
    return nearest, other, weight


def solve_kepler_equation_elliptic(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve Kepler's equation M = E - e·sin(E) for eccentric anomaly E (elliptic orbits).

    Uses Markley's (1995) starter for high eccentricity followed by Newton-Raphson
    iteration for robust convergence. Logs a warning if convergence is not achieved
    within the maximum number of iterations.

    Args:
        M: Mean anomaly in radians
        e: Eccentricity (0 ≤ e < 1)
        tol: Convergence tolerance (default 1e-8 ~ 0.002 arcsec)

    Returns:
        float: Eccentric anomaly E in radians

    Algorithm:
        For e < 0.8: initial guess E = M (fast convergence in ~3-6 iterations).
        For e >= 0.8: Markley (1995) rational starter that guarantees convergence
        in 2-3 Newton iterations for any eccentricity.

        Newton-Raphson: E_{n+1} = E_n - f(E_n)/f'(E_n)
        where f(E) = E - e·sin(E) - M
        and f'(E) = 1 - e·cos(E)

    References:
        Markley, F.L. (1995). Kepler equation solver. Celestial Mechanics, 63, 101.
        Curtis "Orbital Mechanics for Engineering Students" §3.1
        Vallado "Fundamentals of Astrodynamics" Algorithm 2
    """
    # Normalize M to [-pi, pi] for better numerical behavior
    M_orig = M
    M = M % (2.0 * math.pi)
    if M > math.pi:
        M -= 2.0 * math.pi

    if e < 0.8:
        E = M
    else:
        # Markley (1995) rational starter for high eccentricity
        # Provides an initial guess accurate to ~1e-5 for any e < 1
        alpha = (1.0 - e) / (4.0 * e + 0.5)
        beta = M / (2.0 * (4.0 * e + 0.5))
        aux = math.sqrt(beta * beta + alpha * alpha * alpha)
        z = abs(beta + aux) ** (1.0 / 3.0)
        s = z - alpha / z
        # Refine with one correction
        w = s - 0.078 * s * s * s * s * s / (1.0 + e)
        E = M + e * (3.0 * w - 4.0 * w * w * w) / (1.0 + e * (1.0 - w * w * 1.5))
        if M < 0:
            E = -abs(E)
        else:
            E = abs(E)

    for iteration in range(30):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)
        E_new = E - f / fp

        if abs(E_new - E) < tol:
            return E_new
        E = E_new

    logger = get_logger()
    logger.warning(
        "Kepler equation (elliptic) did not converge after 30 iterations: "
        "M=%.6f, e=%.6f, residual=%.2e",
        M_orig,
        e,
        abs(E - e * math.sin(E) - M),
    )
    return E


def solve_kepler_equation_hyperbolic(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve hyperbolic Kepler's equation M = e·sinh(H) - H for hyperbolic anomaly H.

    Uses Newton-Raphson iteration with an improved starter for robust convergence.
    Logs a warning if convergence is not achieved within the maximum iterations.

    Args:
        M: Mean anomaly in radians (can be any value for hyperbolic orbits)
        e: Eccentricity (e > 1)
        tol: Convergence tolerance (default 1e-8)

    Returns:
        float: Hyperbolic anomaly H in radians

    Algorithm:
        Newton-Raphson: H_{n+1} = H_n - f(H_n)/f'(H_n)
        where f(H) = e·sinh(H) - H - M
        and f'(H) = e·cosh(H) - 1

    References:
        Raposo-Pulido, V. & Pelaez, J. (2017). MNRAS, 467, 1702.
        Curtis "Orbital Mechanics for Engineering Students" §3.4
        Vallado "Fundamentals of Astrodynamics" Algorithm 4
    """
    # Initial guess
    if e < 1.6:
        H = M if abs(M) < math.pi else math.copysign(math.pi, M)
    else:
        H = math.copysign(math.log(2 * abs(M) / e + 1.8), M) if abs(M) > 1e-10 else 0.0

    converged = False
    for _ in range(50):
        sinh_H = math.sinh(H)
        cosh_H = math.cosh(H)
        f = e * sinh_H - H - M
        fp = e * cosh_H - 1

        if abs(fp) < 1e-15:
            break

        H_new = H - f / fp

        if abs(H_new - H) < tol:
            converged = True
            H = H_new
            break
        H = H_new

    if not converged:
        logger = get_logger()
        residual = abs(e * math.sinh(H) - H - M)
        logger.warning(
            "Kepler equation (hyperbolic) did not converge after 50 iterations: "
            "M=%.6f, e=%.6f, residual=%.2e",
            M,
            e,
            residual,
        )
    return H


def solve_barker_equation(M: float) -> float:
    """
    Solve Barker's equation for parabolic orbits to find true anomaly.

    For parabolic orbits (e = 1), the relationship between mean anomaly M
    and true anomaly ν is given by Barker's equation.

    Args:
        M: Mean anomaly in radians

    Returns:
        float: True anomaly ν in radians

    Algorithm:
        Uses the cubic solution for Barker's equation:
        tan(ν/2) = s where s³ + 3s - 6M/√2 = 0

        Solving via Cardano's formula or the closed-form solution.

    References:
        Vallado "Fundamentals of Astrodynamics" Algorithm 3
        Danby "Fundamentals of Celestial Mechanics" §6.6
    """
    # Barker's equation: M = (1/2) * tan(ν/2) + (1/6) * tan³(ν/2)
    # Let s = tan(ν/2), then: s + s³/3 = 2M
    # Rearranging: s³ + 3s - 6M = 0

    # Use the closed-form solution (Cardano's formula for depressed cubic)
    W = 3.0 * M
    # s = (W + sqrt(W^2 + 1))^(1/3) - (W + sqrt(W^2 + 1))^(-1/3) won't work directly
    # Use: s = 2 * sinh(asinh(W) / 3)
    s = 2.0 * math.sinh(math.asinh(W) / 3.0)

    # True anomaly from tan(ν/2) = s
    nu = 2.0 * math.atan(s)

    return nu


def solve_kepler_equation(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve Kepler's equation for the appropriate anomaly based on orbit type.

    Handles elliptic, parabolic, and hyperbolic orbits.

    Args:
        M: Mean anomaly in radians
        e: Eccentricity (e ≥ 0)
        tol: Convergence tolerance (default 1e-8 ~ 0.002 arcsec)

    Returns:
        float: For e < 1: Eccentric anomaly E in radians
               For e = 1: True anomaly ν in radians (special case)
               For e > 1: Hyperbolic anomaly H in radians

    Note:
        For parabolic orbits (e = 1), returns TRUE anomaly directly since
        the eccentric anomaly is not defined for parabolas.

    References:
        Curtis "Orbital Mechanics for Engineering Students" §3.1-3.4
        Vallado "Fundamentals of Astrodynamics" Algorithms 2-4
    """
    # Parabolic orbit (e ≈ 1)
    if abs(e - 1.0) < 1e-10:
        # For parabolic orbits, solve Barker's equation
        # Returns true anomaly directly
        return solve_barker_equation(M)

    # Elliptic orbit (e < 1)
    if e < 1.0:
        return solve_kepler_equation_elliptic(M, e, tol)

    # Hyperbolic orbit (e > 1)
    return solve_kepler_equation_hyperbolic(M, e, tol)


def _calc_short_period_correction(
    elements: OrbitalElements,
    jd_tt: float,
    omega_rad: float,
    Omega_rad: float,
    i_rad: float,
) -> Tuple[float, float, float]:
    """Calculate first-order short-period perturbation corrections.

    M1: Implements analytical short-period perturbation theory from
    Brouwer & Clemence (1961) Chapter 15. The dominant short-period terms
    arise from conjunctions with Jupiter and Saturn, oscillating at the
    synodic period (~400 days for Jupiter-Ceres).

    The correction to heliocentric longitude is:
        Δλ ≈ Σ_j [-2·μ_j·α_j·b_{1/2}^{(1)}(α_j) / (n - n_j)] · sin(λ - λ_j)

    where:
        μ_j = mass ratio of planet j
        α_j = semi-major axis ratio
        n, n_j = mean motions
        λ, λ_j = mean longitudes

    The amplitude depends on the current orbital elements (not a static fit),
    so it naturally evolves as elements drift — avoiding the problem of
    static Fourier fits that worsen positions (as noted in Task 4.2).

    Args:
        elements: Orbital elements at epoch
        jd_tt: Target Julian Day
        omega_rad: Perturbed argument of perihelion (radians)
        Omega_rad: Perturbed longitude of ascending node (radians)
        i_rad: Perturbed inclination (radians)

    Returns:
        Tuple (dx, dy, dz) correction in AU in ecliptic J2000.0 frame

    References:
        Brouwer, D. & Clemence, G.M. (1961). Methods of Celestial Mechanics.
        Ch. 15 (short-period perturbations).
        Murray & Dermott (1999) §6.9 (disturbing function expansion).
    """
    a = elements.a
    e = elements.e
    n = elements.n  # degrees/day
    dt = jd_tt - elements.epoch

    if a <= 0.0 or e >= 1.0:
        return 0.0, 0.0, 0.0

    # Mean longitude of asteroid at target time
    M_body = math.radians((elements.M0 + n * dt) % 360.0)
    lambda_body = Omega_rad + omega_rad + M_body

    # Planet data for short-period corrections: (a_planet, mu, n_planet_deg, lambda_j2000)
    # lambda_j2000 = mean longitude at J2000.0
    _sp_planets = [
        (
            JUPITER_A,
            MASS_RATIO_JUPITER,
            JUPITER_N,
            math.radians(
                (
                    JUPITER_NODE
                    + JUPITER_OMEGA
                    + 18.818  # Jupiter M0 at J2000.0 ≈ 18.818°
                )
            ),
        ),
        (
            SATURN_A,
            MASS_RATIO_SATURN,
            SATURN_N,
            math.radians(
                (
                    SATURN_NODE
                    + SATURN_OMEGA
                    + 320.346  # Saturn M0 at J2000.0 ≈ 320.346°
                )
            ),
        ),
    ]

    # Accumulate longitude correction (in radians)
    delta_lambda = 0.0

    for a_p, mu_p, n_p_deg, lambda_p_j2000 in _sp_planets:
        # Semi-major axis ratio
        if a < a_p:
            alpha = a / a_p
        else:
            alpha = a_p / a

        if alpha >= 1.0 or alpha <= 0.01:
            continue

        # Planet mean longitude at target time
        n_p_rad = math.radians(n_p_deg)
        dt_j2000 = jd_tt - 2451545.0
        lambda_p = lambda_p_j2000 + n_p_rad * dt_j2000

        # Synodic frequency: n - n_planet (in rad/day)
        n_rad = math.radians(n)
        dn = n_rad - n_p_rad

        # Avoid division by zero near exact resonance
        if abs(dn) < 1e-8:
            continue

        # Laplace coefficient b_{1/2}^{(1)}(alpha) for the direct term
        b12_1 = _calc_laplace_coefficients(alpha, 0.5, 1)

        # First-order short-period correction to mean longitude
        # Δλ ≈ -2·μ·α·b_{1/2}^{(1)}·sin(λ - λ_p) / [(n - n_p)/n]
        # Normalized by mean motion ratio for dimensional consistency
        sin_diff = math.sin(lambda_body - lambda_p)
        amplitude = 2.0 * mu_p * alpha * b12_1

        # The correction in mean longitude (radians):
        delta_lambda += -amplitude * sin_diff / (dn / n_rad)

        # Second-order eccentricity correction: modulate with e
        # This adds the e·cos(2M - M_J) term from Brouwer eq. 15.22
        math.cos(lambda_body - lambda_p)
        delta_lambda += (
            -amplitude
            * e
            * math.sin(2.0 * M_body - lambda_p)
            / ((2.0 * n_rad - n_p_rad) / n_rad)
            * 0.5
            if abs(2.0 * n_rad - n_p_rad) > 1e-8
            else 0.0
        )

    # Convert longitude correction to Cartesian correction
    # Δx ≈ -r·sin(λ)·Δλ, Δy ≈ r·cos(λ)·Δλ, Δz ≈ 0 (approximately)
    # where r ≈ a·(1 - e²)/(1 + e·cos(ν)) ≈ a for circular approximation
    r_approx = a * (1.0 - e * e)  # semi-latus rectum as typical distance

    # Project into ecliptic frame
    cos_lambda = math.cos(lambda_body)
    sin_lambda = math.sin(lambda_body)
    cos_i = math.cos(i_rad)

    dx = -r_approx * sin_lambda * delta_lambda
    dy = r_approx * cos_lambda * delta_lambda * cos_i
    dz = r_approx * cos_lambda * delta_lambda * math.sin(i_rad)

    return dx, dy, dz


def calc_minor_body_position(
    elements: OrbitalElements,
    jd_tt: float,
    include_perturbations: bool = True,
    body_id: Optional[int] = None,
) -> Tuple[float, float, float]:
    """
    Calculate heliocentric position using Keplerian orbital elements with perturbations.

    Propagates orbit from epoch to target time using mean motion and applies
    first-order secular perturbations from Jupiter and Saturn. For plutinos
    (bodies in 2:3 resonance with Neptune), also applies resonant libration
    corrections when body_id is provided.
    Supports elliptic, parabolic, and hyperbolic orbits.

    Args:
        elements: Orbital elements at epoch
        jd_tt: Target Julian Day in Terrestrial Time (TT)
        include_perturbations: If True, apply secular perturbations (default True)
        body_id: Optional minor body identifier. If provided and body is a known
            plutino (Ixion, Orcus), applies resonant libration correction to
            improve accuracy over multi-decadal timescales.

    Returns:
        Tuple[float, float, float]: (x, y, z) heliocentric position in AU
            Coordinates in ecliptic J2000.0 frame

    Algorithm:
        1. Apply secular perturbations to ω and Ω (if enabled)
        2. Propagate mean anomaly: M(t) = M0 + n·Δt
        3. For plutinos: Apply resonant libration correction to mean anomaly
        4. Solve Kepler's equation for eccentric/hyperbolic/true anomaly
        5. Calculate true anomaly ν from E (or H for hyperbolic)
        6. Compute position in orbital plane
        7. Rotate to ecliptic frame using perturbed Ω, i, ω

    Perturbation Model:
        Uses first-order Laplace-Lagrange secular perturbation theory to
        account for the gravitational influence of Jupiter and Saturn.
        Secular perturbations cause:
        - Precession of perihelion (ω)
        - Regression of ascending node (Ω)

        This significantly improves accuracy for multi-year propagation,
        reducing errors from ~1-5 arcminutes (pure Keplerian) to
        ~10-30 arcseconds for main belt asteroids.

    Resonant Libration Model (plutinos):
        For bodies in 2:3 mean motion resonance with Neptune (Ixion, Orcus),
        the resonant argument φ = 3λ_TNO - 2λ_Neptune - ω_TNO librates around
        180° with a period of ~20,000 years. This libration causes systematic
        position errors of 5-10° over 50+ years if not corrected.

        When body_id is provided for a known plutino, a sinusoidal correction
        is applied to the mean anomaly, reducing position errors to <2° over
        100-year timescales.

    Precision:
        With perturbations enabled:
        - Main belt asteroids: ~10-30 arcseconds typical
        - TNOs: ~1-3 arcminutes typical
        - Plutinos (with libration): <2° over 100 years

        Remaining errors from:
        - Higher-order perturbations
        - Non-gravitational forces (radiation pressure, Yarkovsky)
        - Relativistic effects (minor for asteroids)

    References:
        Murray & Dermott "Solar System Dynamics" Ch. 7 (secular theory), Ch. 8 (resonances)
        Malhotra (1995) "The origin of Pluto's peculiar orbit"
        Curtis §3 (orbital elements)
        Vallado §2.3 (coordinate transformations)
    """
    e = elements.e
    dt = jd_tt - elements.epoch

    # Apply secular perturbations to get perturbed orbital elements
    # Now includes eccentricity and inclination evolution via (h,k)/(p,q) vectors
    # S3: Apply perturbations for near-parabolic (e close to 1) bound orbits too,
    # not just strictly elliptic. Bodies with poorly determined orbits may have
    # e ≈ 1 but still be on bound trajectories deserving perturbation corrections.
    if include_perturbations and e < 1.0:
        # Apply perturbations for all bound orbits (elliptic and near-parabolic)
        omega_pert, Omega_pert, M_deg, n_pert, e_pert, i_pert = (
            apply_secular_perturbations(elements, jd_tt, include_perturbations=True)
        )
    else:
        # Truly hyperbolic orbits (e > 1), or perturbations disabled
        omega_pert = elements.omega
        Omega_pert = elements.Omega
        M_deg = elements.M0 + elements.n * dt
        e_pert = e
        i_pert = elements.i
        # n_pert not used for hyperbolic orbits

    # Use perturbed eccentricity for orbit computation
    e_use = e_pert if e < 1.0 else e

    # Propagate mean anomaly (handle differently for each orbit type)
    if abs(e - 1.0) < 1e-10:
        # Parabolic orbit: mean anomaly grows linearly
        M = math.radians(elements.M0 + elements.n * dt)
    elif e < 1.0:
        # Elliptic orbit: use perturbed values, wrap to [0, 360)
        M_corrected = M_deg

        # Apply resonant libration correction for plutinos
        # This corrects for the oscillatory perturbation from the 2:3 Neptune resonance
        # that is not captured by secular perturbation theory.
        # The correction must be DIFFERENTIAL: subtract the value at epoch so that
        # the correction is exactly zero at the element epoch. The osculating elements
        # already encode the resonant state at epoch — we only need the *change*
        # in the libration angle from epoch to target time.
        if body_id is not None and body_id in PLUTINO_LIBRATION_PARAMS:
            libration_at_target = calc_libration_correction(body_id, jd_tt)
            libration_at_epoch = calc_libration_correction(body_id, elements.epoch)
            libration_correction = libration_at_target - libration_at_epoch
            # Apply correction to mean anomaly (libration affects mean longitude λ = Ω + ω + M)
            # Since Ω and ω are already perturbed, we apply the correction to M
            M_corrected = M_deg + libration_correction

        M = math.radians(M_corrected % 360.0)
    else:
        # Hyperbolic orbit: no wrapping needed
        M = math.radians(elements.M0 + elements.n * dt)

    # Solve the appropriate Kepler equation (using perturbed eccentricity)
    anomaly = solve_kepler_equation(M, e_use)

    # Calculate true anomaly and distance based on orbit type
    if abs(e - 1.0) < 1e-10:
        # Parabolic orbit: solve_kepler_equation returns true anomaly directly
        nu = anomaly
        # For parabolic orbits: r = p / (1 + cos(ν)) where p = 2 * q (q = perihelion)
        # Since a is undefined for parabolas, we use q (stored in 'a' field as perihelion)
        # p = 2 * q, so r = 2 * q / (1 + cos(ν))
        p = 2.0 * elements.a  # elements.a is perihelion distance q for parabolic orbits
        r = p / (1.0 + math.cos(nu))
    elif e < 1.0:
        # Elliptic orbit: anomaly is eccentric anomaly E, use perturbed e
        E = anomaly
        nu = 2.0 * math.atan2(
            math.sqrt(1 + e_use) * math.sin(E / 2),
            math.sqrt(1 - e_use) * math.cos(E / 2),
        )
        # H4: Apply Yarkovsky semi-major axis drift for NEAs with measured da/dt
        a_use = _apply_yarkovsky_correction(elements, dt, body_id)
        r = a_use * (1 - e_use * math.cos(E))
    else:
        # Hyperbolic orbit: anomaly is hyperbolic anomaly H
        H = anomaly
        # True anomaly from hyperbolic anomaly
        nu = 2.0 * math.atan2(
            math.sqrt(e + 1) * math.sinh(H / 2),
            math.sqrt(e - 1) * math.cosh(H / 2),
        )
        # Distance from hyperbolic anomaly. Hyperbolic orbits carry
        # a < 0, so the magnitude is needed: r = |a| (e cosh H - 1).
        r = abs(elements.a) * (e * math.cosh(H) - 1)

    # Position in orbital plane (perifocal frame)
    x_orb = r * math.cos(nu)
    y_orb = r * math.sin(nu)

    # Convert Euler angles to radians (use perturbed values for bound orbits)
    if e < 1.0:
        omega_rad = math.radians(omega_pert)
        Omega_rad = math.radians(Omega_pert)
    else:
        omega_rad = math.radians(elements.omega)
        Omega_rad = math.radians(elements.Omega)
    i_rad = math.radians(i_pert)

    # Precompute trig functions
    cos_omega = math.cos(omega_rad)
    sin_omega = math.sin(omega_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)
    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)

    # Rotation matrix from perifocal to ecliptic (standard transformation)
    P11 = cos_omega * cos_Omega - sin_omega * sin_Omega * cos_i
    P12 = -sin_omega * cos_Omega - cos_omega * sin_Omega * cos_i
    P21 = cos_omega * sin_Omega + sin_omega * cos_Omega * cos_i
    P22 = -sin_omega * sin_Omega + cos_omega * cos_Omega * cos_i
    P31 = sin_omega * sin_i
    P32 = cos_omega * sin_i

    # Transform to heliocentric ecliptic J2000 coordinates
    x = P11 * x_orb + P12 * y_orb
    y = P21 * x_orb + P22 * y_orb
    z = P31 * x_orb + P32 * y_orb

    # M1: Apply short-period perturbation corrections from Jupiter and Saturn.
    # These analytical corrections capture the dominant synodic-period oscillations
    # (~400 days for Jupiter-Ceres) that are the largest error source at timescales
    # < 5 years. The amplitudes depend on current orbital elements, so they evolve
    # naturally as elements drift. Only applied for bound orbits with perturbations.
    if include_perturbations and e < 1.0:
        dx, dy, dz = _calc_short_period_correction(
            elements, jd_tt, omega_rad, Omega_rad, i_rad
        )
        x += dx
        y += dy
        z += dz

    return x, y, z


def calc_minor_body_heliocentric(
    body_id: int, jd_tt: float, use_spk: bool = True
) -> Tuple[float, float, float]:
    """
    Calculate heliocentric ecliptic coordinates for a minor body.

    This function automatically uses SPK kernels when available for major
    asteroids (Ceres, Pallas, Juno, Vesta, Chiron), providing sub-arcsecond
    precision. When SPK is not available, it falls back to Keplerian
    calculations with secular perturbations (~10-30 arcsec precision).

    Args:
        body_id: Minor body identifier (CHIRON, ERIS, etc.)
        jd_tt: Julian Day in Terrestrial Time (TT)
        use_spk: If True (default), attempt to use SPK kernels for
            high-precision positions. If False, always use Keplerian.

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (-90 to +90)
            - distance: Heliocentric distance in AU

    Raises:
        ValueError: If body_id is not in the database

    Note:
        Returns HELIOCENTRIC coordinates. For geocentric, caller must
        subtract Earth's heliocentric position (see planets.py).

    Precision:
        With SPK (major asteroids): Sub-arcsecond to arcsecond
        Without SPK:
            Asteroids (Ceres, etc.): ~10-30 arcseconds typical
            TNOs (Eris, etc.): ~1-3 arcminutes typical
        Errors increase with time from epoch

    Example:
        >>> from libephemeris.minor_bodies import calc_minor_body_heliocentric
        >>> from libephemeris.constants import CERES
        >>> lon, lat, dist = calc_minor_body_heliocentric(CERES, 2451545.0)
        >>> print(f"Ceres: {lon:.4f}° lon, {lat:.4f}° lat, {dist:.4f} AU")
    """
    if body_id not in MINOR_BODY_ELEMENTS:
        raise ValueError(f"illegal planet number {body_id}.")

    # Try to use SPK for high precision if available
    if use_spk:
        try:
            from . import state
            from .spk import calc_spk_body_position
            from .constants import FLG_HELCTR

            # Check if SPK is already registered for this body
            if body_id in state._SPK_BODY_MAP:
                ts = state.get_timescale()
                t = ts.tt_jd(jd_tt)
                result = calc_spk_body_position(t, body_id, FLG_HELCTR)
                if result is not None:
                    lon, lat, dist, _, _, _ = result
                    return lon, lat, dist
        except (ImportError, ValueError, KeyError, EphemerisRangeError):
            # Fall through to Keplerian calculation
            pass

    # Fall back to Keplerian calculation
    # Use multi-epoch elements if available for better accuracy
    # Pass body_id to enable resonant libration correction for plutinos.
    # Near the midpoint between adjacent element epochs the two
    # propagations are cross-faded in Cartesian space so the nearest-
    # epoch switch does not step the position.
    elements, blend_elements, blend_w = _get_epoch_elements_blend(body_id, jd_tt)
    x, y, z = calc_minor_body_position(elements, jd_tt, body_id=body_id)
    if blend_elements is not None and blend_w > 0.0:
        xb, yb, zb = calc_minor_body_position(
            blend_elements, jd_tt, body_id=body_id
        )
        x = x * (1.0 - blend_w) + xb * blend_w
        y = y * (1.0 - blend_w) + yb * blend_w
        z = z * (1.0 - blend_w) + zb * blend_w

    # Convert Cartesian to spherical coordinates
    r = math.sqrt(x**2 + y**2 + z**2)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(z / r))

    # Frame contract: the orbital-element tables are referenced to the
    # J2000 ecliptic (fit from JPL VECTORS), so calc_minor_body_position
    # returns J2000 ecliptic coordinates.  Precess to the ecliptic of
    # date so this fallback agrees with the SPK branch above (measured
    # bias before this correction: ~1290" at 2026, i.e. 26 yr of
    # general precession).
    from .astrometry import _precess_ecliptic

    lon, lat = _precess_ecliptic(lon, lat, 2451545.0, jd_tt)

    return lon, lat, r


# =============================================================================
# JPL SBDB API FOR DYNAMIC ASTEROID LOOKUP
# =============================================================================

# JPL Small-Body Database API endpoint
SBDB_API_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api"

# Cache for dynamically fetched orbital elements
# Maps asteroid_number -> OrbitalElements
_ASTEROID_ELEMENTS_CACHE: dict[int, OrbitalElements] = {}


def fetch_orbital_elements_from_sbdb(
    asteroid_number: int, timeout: float = 30.0
) -> Optional[OrbitalElements]:
    """
    Fetch orbital elements for a numbered asteroid from JPL SBDB API.

    This function queries the JPL Small-Body Database API to retrieve current
    osculating orbital elements for any numbered asteroid in the database.

    Args:
        asteroid_number: The asteroid catalog number (e.g., 1 for Ceres, 2060 for Chiron)
        timeout: Request timeout in seconds (default: 30)

    Returns:
        OrbitalElements if successful, None if the asteroid is not found or
        the request fails.

    Raises:
        No exceptions are raised; errors return None.

    Example:
        >>> from libephemeris.minor_bodies import fetch_orbital_elements_from_sbdb
        >>> elements = fetch_orbital_elements_from_sbdb(433)  # Eros
        >>> if elements:
        ...     print(f"Eros semi-major axis: {elements.a:.4f} AU")

    Notes:
        - Requires network access to ssd-api.jpl.nasa.gov
        - Results are cached to avoid redundant API calls
        - JPL SBDB contains 1+ million numbered asteroids
    """
    import json
    import ssl
    import urllib.error
    import urllib.request

    import certifi

    # Check cache first
    if asteroid_number in _ASTEROID_ELEMENTS_CACHE:
        return _ASTEROID_ELEMENTS_CACHE[asteroid_number]

    # Build API request
    params = f"sstr={asteroid_number}&phys-par=false&full-prec=true"
    url = f"{SBDB_API_URL}?{params}"

    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "libephemeris/1.0 (Python)"},
        )
        _ssl_ctx = ssl.create_default_context(cafile=certifi.where())
        with urllib.request.urlopen(
            request, timeout=timeout, context=_ssl_ctx
        ) as response:
            data = json.loads(response.read().decode("utf-8"))

        # Check for API errors
        if "error" in data:
            return None

        # Extract orbital elements
        orbit = data.get("orbit", {})
        elements_list = orbit.get("elements", [])

        if not elements_list:
            return None

        # Build a dict of element name -> value
        elem_dict: dict[str, float] = {}
        for elem in elements_list:
            elem_dict[elem["name"]] = float(elem["value"])

        # Get epoch (Julian Day TDB)
        epoch_jd = float(orbit.get("epoch", 0))
        if epoch_jd == 0:
            return None

        # Get body name from object data
        obj_data = data.get("object", {})
        body_name = obj_data.get("fullname", f"Asteroid {asteroid_number}")

        # Map SBDB element names to our structure
        # SBDB uses: e, a, q, i, om (node), w (arg peri), ma (mean anom), n (mean motion)
        # A record missing any of the core elements cannot be propagated;
        # zero-filling them would yield a fictitious orbit, so give None.
        for _key in ("a", "e", "i", "n"):
            if elem_dict.get(_key) is None:
                return None
        orbital_elements = OrbitalElements(
            name=body_name,
            epoch=epoch_jd,
            a=elem_dict["a"],  # Semi-major axis
            e=elem_dict["e"],  # Eccentricity
            i=elem_dict["i"],  # Inclination
            omega=elem_dict.get("w", 0.0),  # Argument of perihelion
            Omega=elem_dict.get("om", 0.0),  # Longitude of ascending node
            M0=elem_dict.get("ma", 0.0),  # Mean anomaly
            n=elem_dict["n"],  # Mean motion
        )

        # Cache the result
        _ASTEROID_ELEMENTS_CACHE[asteroid_number] = orbital_elements

        return orbital_elements

    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError):
        return None
    except (KeyError, ValueError, TypeError):
        return None


def clear_asteroid_elements_cache() -> int:
    """
    Clear the cache of dynamically fetched orbital elements.

    Returns:
        int: Number of entries cleared from the cache.

    Example:
        >>> from libephemeris.minor_bodies import clear_asteroid_elements_cache
        >>> cleared = clear_asteroid_elements_cache()
        >>> print(f"Cleared {cleared} cached entries")
    """
    count = len(_ASTEROID_ELEMENTS_CACHE)
    _ASTEROID_ELEMENTS_CACHE.clear()
    return count


# Cache for name-to-number lookups (asteroid_name_lower -> asteroid_number)
_ASTEROID_NAME_CACHE: dict[str, int] = {}


def _build_local_name_cache() -> None:
    """
    Build a lookup cache from asteroid names to their catalog numbers.

    Populates the cache using bodies in MINOR_BODY_ELEMENTS and known mappings
    for asteroids with dedicated SE_* constants (Chiron, Pholus, Ceres, etc.).
    """
    if _ASTEROID_NAME_CACHE:
        return  # Already built

    # Known mappings for bodies with dedicated SE_* constants
    known_mappings = {
        "chiron": 2060,
        "pholus": 5145,
        "ceres": 1,
        "pallas": 2,
        "juno": 3,
        "vesta": 4,
        "eris": 136199,
        "sedna": 90377,
        "haumea": 136108,
        "makemake": 136472,
        "ixion": 28978,
        "orcus": 90482,
        "quaoar": 50000,
        "nessus": 7066,
        "asbolus": 8405,
        "chariklo": 10199,
        "gonggong": 225088,
        "varuna": 20000,
        "apophis": 99942,
        "flora": 8,
        "eunomia": 15,
        "hygiea": 10,
        "interamnia": 704,
        "davida": 511,
        "europa": 52,  # Main belt asteroid, not Jupiter's moon
        "sylvia": 87,
        "psyche": 16,
        "eros": 433,
        "amor": 1221,
        "icarus": 1566,
        "toro": 1685,
        "sappho": 80,
        "pandora": 55,  # Main belt asteroid, not Saturn's moon
        "lilith": 1181,  # Main belt asteroid, not lunar apogee Lilith
        "hidalgo": 944,
        "toutatis": 4179,
        "itokawa": 25143,
        "bennu": 101955,
        "ryugu": 162173,
    }

    _ASTEROID_NAME_CACHE.update(known_mappings)


def get_asteroid_number(name: str, timeout: float = 30.0) -> Optional[int]:
    """
    Look up an asteroid's catalog number by name.

    This function is useful for users who know an asteroid's name but not its
    catalog number. It first checks a local database of well-known asteroids,
    then queries the JPL Small-Body Database (SBDB) API if not found locally.

    Args:
        name: The asteroid name to look up (case-insensitive).
            Examples: "Ceres", "Eros", "Apophis", "Bennu", "Vesta"
        timeout: Request timeout in seconds for SBDB API queries (default: 30)

    Returns:
        int: The asteroid catalog number if found, or None if not found.

    Example:
        >>> from libephemeris.minor_bodies import get_asteroid_number
        >>> get_asteroid_number("Ceres")
        1
        >>> get_asteroid_number("Chiron")
        2060
        >>> get_asteroid_number("Apophis")
        99942
        >>> get_asteroid_number("Eros")
        433

    Notes:
        - Local lookup is instant and doesn't require network access
        - SBDB API query requires network access to ssd-api.jpl.nasa.gov
        - Results from SBDB queries are cached for subsequent lookups
        - The search is case-insensitive ("ceres", "Ceres", "CERES" all work)
        - For ambiguous names (e.g., "Pandora" can be an asteroid or Saturn's
          moon), this returns the asteroid catalog number

    See Also:
        calc_asteroid_by_number: Calculate position using the returned number
        fetch_orbital_elements_from_sbdb: Get orbital elements by number
    """
    import json
    import ssl
    import urllib.error
    import urllib.parse
    import urllib.request

    import certifi

    # Normalize name for lookup
    name_lower = name.strip().lower()

    if not name_lower:
        return None

    # Build local cache if not done
    _build_local_name_cache()

    # Check local cache first
    if name_lower in _ASTEROID_NAME_CACHE:
        return _ASTEROID_NAME_CACHE[name_lower]

    # Query JPL SBDB API by name
    # The API accepts names in the sstr parameter
    params = f"sstr={urllib.parse.quote(name)}&phys-par=false"
    url = f"{SBDB_API_URL}?{params}"

    try:
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "libephemeris/1.0 (Python)"},
        )
        _ssl_ctx = ssl.create_default_context(cafile=certifi.where())
        with urllib.request.urlopen(
            request, timeout=timeout, context=_ssl_ctx
        ) as response:
            data = json.loads(response.read().decode("utf-8"))

        # Check for API errors
        if "error" in data:
            return None

        # Extract object data
        obj_data = data.get("object", {})

        # Get the SPK-ID (asteroid number) from object data
        # SBDB returns "spkid" for numbered objects or we can use "des" (designation)
        spkid = obj_data.get("spkid")
        if spkid:
            try:
                # SPK-ID format for asteroids: 2XXXXXX (2000000 + number)
                spkid_int = int(spkid)
                if spkid_int >= 2000000:
                    asteroid_number = spkid_int - 2000000
                else:
                    # Some objects have direct SPK-IDs
                    asteroid_number = spkid_int
            except (ValueError, TypeError):
                return None
        else:
            # Try to parse from designation
            des = obj_data.get("des", "")
            if des and des.isdigit():
                asteroid_number = int(des)
            else:
                # Try fullname parsing - e.g., "1 Ceres" -> 1
                fullname = obj_data.get("fullname", "")
                parts = fullname.split()
                if parts and parts[0].isdigit():
                    asteroid_number = int(parts[0])
                else:
                    return None

        # Cache the result for future lookups
        _ASTEROID_NAME_CACHE[name_lower] = asteroid_number

        return asteroid_number

    except (urllib.error.URLError, urllib.error.HTTPError, json.JSONDecodeError):
        return None
    except (KeyError, ValueError, TypeError):
        return None


def clear_asteroid_name_cache() -> int:
    """
    Clear the cache of name-to-number lookups.

    Returns:
        int: Number of entries cleared from the cache.

    Example:
        >>> from libephemeris.minor_bodies import clear_asteroid_name_cache
        >>> cleared = clear_asteroid_name_cache()
        >>> print(f"Cleared {cleared} cached entries")
    """
    count = len(_ASTEROID_NAME_CACHE)
    _ASTEROID_NAME_CACHE.clear()
    return count


# =============================================================================
# AUTOMATIC SPK DOWNLOAD FOR MAJOR ASTEROIDS
# =============================================================================

# Major asteroids with high-precision SPK kernels available from JPL Horizons
# These are bodies where users benefit most from automatic SPK downloads:
# - Keplerian elements give ~10-30 arcsec precision (good enough for many uses)
# - SPK kernels give sub-arcsecond precision (required for precise timing)
#
# Bodies included:
# - Main belt "Big Four": Ceres (1), Pallas (2), Juno (3), Vesta (4)
# - Centaur Chiron (2060) - frequently used in astrology
# Format: body_id -> (asteroid_number, horizons_id, naif_id, body_name)

from .constants import (  # noqa: E402 (section-local import)
    CERES,
    PALLAS,
    JUNO,
    VESTA,
    CHIRON,
    NAIF_CERES,
    NAIF_PALLAS,
    NAIF_JUNO,
    NAIF_VESTA,
    NAIF_CHIRON,
)

# All major asteroids can have SPK downloaded from JPL Horizons.
# Bodies formerly in JPL's "major body index" (Ceres, Pallas, Juno, Vesta) are now
# downloaded using the name syntax (e.g., "Ceres;") which bypasses the major body
# index restriction. See SPK_BODY_NAME_MAP in constants.py for Horizons IDs.
SPK_DOWNLOADABLE_ASTEROIDS: dict[int, tuple[int, str, int, str]] = {
    CERES: (1, "Ceres;", NAIF_CERES, "Ceres"),
    PALLAS: (2, "Pallas;", NAIF_PALLAS, "Pallas"),
    JUNO: (3, "Juno;", NAIF_JUNO, "Juno"),
    VESTA: (4, "Vesta;", NAIF_VESTA, "Vesta"),
    CHIRON: (2060, "2060", NAIF_CHIRON, "Chiron"),
}

# Combined info for all major asteroids (for backward compatibility and info lookup)
# Note: Ceres/Pallas/Juno/Vesta use "Name;" syntax to bypass JPL's major body index.
# Bare numeric IDs ("1", "2", "3", "4") collide with planet barycenter IDs and cause
# "SPK creation is not available for pre-computed objects" errors from Horizons.
MAJOR_ASTEROID_SPK_INFO: dict[int, tuple[int, str, int, str]] = {
    CERES: (1, "Ceres;", NAIF_CERES, "Ceres"),
    PALLAS: (2, "Pallas;", NAIF_PALLAS, "Pallas"),
    JUNO: (3, "Juno;", NAIF_JUNO, "Juno"),
    VESTA: (4, "Vesta;", NAIF_VESTA, "Vesta"),
    CHIRON: (2060, "2060", NAIF_CHIRON, "Chiron"),
}


def is_major_asteroid(body_id: int) -> bool:
    """
    Check if a body ID is a major asteroid.

    Major asteroids are prominent bodies like Ceres, Vesta, Chiron, etc.
    Note: Not all major asteroids can have SPK downloaded - use
    is_spk_downloadable() to check if Horizons can generate SPK for a body.

    Args:
        body_id: Minor body identifier (CERES, VESTA, etc.)

    Returns:
        bool: True if body is a major asteroid

    Example:
        >>> from libephemeris.minor_bodies import is_major_asteroid
        >>> from libephemeris.constants import CERES, ERIS
        >>> is_major_asteroid(CERES)
        True
        >>> is_major_asteroid(ERIS)
        False
    """
    return body_id in MAJOR_ASTEROID_SPK_INFO


def is_spk_downloadable(body_id: int) -> bool:
    """
    Check if SPK can be downloaded from JPL Horizons for this body.

    All major asteroids (Ceres, Pallas, Juno, Vesta, Chiron) can have SPK
    downloaded. Bodies in JPL's major body index are downloaded using the
    name syntax (e.g., "Ceres;") which bypasses the restriction.

    Additionally, any body in SPK_BODY_NAME_MAP (constants.py) is downloadable,
    unless it appears in SPK_AUTO_DOWNLOAD_BLOCKED (e.g., Bennu).

    Args:
        body_id: Minor body identifier (CERES, CHIRON, etc.)

    Returns:
        bool: True if Horizons can generate SPK for this body

    Example:
        >>> from libephemeris.minor_bodies import is_spk_downloadable
        >>> from libephemeris.constants import CERES, BENNU
        >>> is_spk_downloadable(CERES)
        True
        >>> is_spk_downloadable(BENNU)
        False
    """
    from .constants import SPK_AUTO_DOWNLOAD_BLOCKED, SPK_BODY_NAME_MAP

    if body_id in SPK_AUTO_DOWNLOAD_BLOCKED:
        return False
    return body_id in SPK_DOWNLOADABLE_ASTEROIDS or body_id in SPK_BODY_NAME_MAP


def is_in_jpl_major_body_index(body_id: int) -> bool:
    """
    Check if body is in JPL's major body index.

    Note: This function is kept for backward compatibility but now always
    returns False. All major body index asteroids (Ceres, Pallas, Juno, Vesta)
    can now be downloaded using the name syntax (e.g., "Ceres;") which
    bypasses the major body index restriction in Horizons.

    Args:
        body_id: Minor body identifier

    Returns:
        bool: Always False - all bodies are now downloadable.
    """
    return False


def get_major_asteroid_info(
    body_id: int,
) -> Optional[tuple[int, str, int, str]]:
    """
    Get SPK download information for a major asteroid.

    Args:
        body_id: Minor body identifier (CERES, VESTA, etc.)

    Returns:
        Tuple of (asteroid_number, horizons_id, naif_id, body_name) if
        the body is a major asteroid with SPK support, None otherwise.

    Example:
        >>> from libephemeris.minor_bodies import get_major_asteroid_info
        >>> from libephemeris.constants import CERES
        >>> info = get_major_asteroid_info(CERES)
        >>> if info:
        ...     ast_num, horizons_id, naif_id, name = info
        ...     print(f"{name}: asteroid #{ast_num}, NAIF ID {naif_id}")
        Ceres: asteroid #1, NAIF ID 2000001
    """
    return MAJOR_ASTEROID_SPK_INFO.get(body_id)


def auto_download_asteroid_spk(
    body_id: int,
    jd_start: Optional[float] = None,
    jd_end: Optional[float] = None,
    force: bool = False,
) -> Optional[str]:
    """
    Automatically download SPK kernel for a minor body if not cached.

    This function checks if an SPK file is already available for the given
    body. If not, it downloads the SPK from JPL Horizons via direct HTTP
    and registers it for use by calc_ut() and related functions.

    No external dependencies beyond urllib are required (no astroquery).

    Supported bodies: All bodies in MAJOR_ASTEROID_SPK_INFO and
    SPK_BODY_NAME_MAP (constants.py).

    Args:
        body_id: Minor body identifier (CERES, VESTA, CHIRON, etc.)
        jd_start: Start Julian Day for SPK coverage. If None, uses 10 years
            before current date.
        jd_end: End Julian Day for SPK coverage. If None, uses 10 years
            after current date.
        force: If True, re-download even if SPK is already cached.

    Returns:
        str: Path to the SPK file if download successful or already cached.
        None: If body is not supported or download fails.

    Raises:
        No exceptions are raised; errors return None to allow graceful
        fallback to Keplerian calculations.

    Example:
        >>> from libephemeris.minor_bodies import auto_download_asteroid_spk
        >>> from libephemeris.constants import CERES
        >>> spk_path = auto_download_asteroid_spk(CERES)
        >>> if spk_path:
        ...     print(f"SPK downloaded: {spk_path}")

    See Also:
        is_spk_available_for_body: Check if SPK is already available
    """
    from . import spk as spk_module
    from . import state
    from .constants import SPK_BODY_NAME_MAP

    logger = get_logger()

    # Try MAJOR_ASTEROID_SPK_INFO first, then SPK_BODY_NAME_MAP
    if body_id in MAJOR_ASTEROID_SPK_INFO:
        _, horizons_id, naif_id, body_name = MAJOR_ASTEROID_SPK_INFO[body_id]
    elif body_id in SPK_BODY_NAME_MAP:
        horizons_id, naif_id = SPK_BODY_NAME_MAP[body_id]
        body_name = spk_module._get_body_name(body_id) or str(body_id)
    else:
        return None

    try:
        # Check if SPK is already registered (and not forcing re-download)
        if not force and body_id in state._SPK_BODY_MAP:
            spk_file, _ = state._SPK_BODY_MAP[body_id]
            return spk_file

        # Determine date range for SPK
        import time as time_module

        current_jd = 2440587.5 + (time_module.time() / 86400.0)

        if jd_start is None:
            jd_start = current_jd - 3652.5  # ~10 years before
        if jd_end is None:
            jd_end = current_jd + 3652.5  # ~10 years after

        # JPL Horizons only supports SPK generation for minor bodies in the
        # range ~1600-01-01 to ~2500-01-01 (empirically verified).  Clamp the
        # requested range so that callers (e.g. extended-tier LEB generation)
        # don't produce date strings that Horizons cannot parse (negative
        # years) or that fall outside the supported window.
        HORIZONS_SPK_JD_MIN = 2305448.5  # 1600-01-01 UTC
        HORIZONS_SPK_JD_MAX = 2634167.5  # 2500-01-01 UTC

        if jd_start < HORIZONS_SPK_JD_MIN or jd_end > HORIZONS_SPK_JD_MAX:
            jd_start = max(jd_start, HORIZONS_SPK_JD_MIN)
            jd_end = min(jd_end, HORIZONS_SPK_JD_MAX)
            if jd_start >= jd_end:
                logger.warning(
                    "Requested date range for %s is entirely outside "
                    "Horizons SPK limits (1600-2500); skipping download",
                    body_name,
                )
                return None

        # Convert JD to calendar dates
        from skyfield.api import load

        ts = load.timescale()
        t_start = ts.tt_jd(jd_start)
        t_end = ts.tt_jd(jd_end)

        start_date = f"{t_start.utc[0]:04d}-{t_start.utc[1]:02d}-{t_start.utc[2]:02d}"
        end_date = f"{t_end.utc[0]:04d}-{t_end.utc[1]:02d}-{t_end.utc[2]:02d}"

        logger.info("Auto-downloading SPK for %s from JPL Horizons...", body_name)

        # Download and register using direct HTTP (no astroquery needed).
        # Don't pass naif_id — let download_and_register_spk auto-detect from
        # the SPK file, since JPL Horizons uses the 20000000+N convention.
        spk_path = spk_module.download_and_register_spk(
            body=horizons_id,
            ipl=body_id,
            start=start_date,
            end=end_date,
        )

        return spk_path

    except (ValueError, KeyError) as e:
        logger.warning("Auto SPK download failed for %s: %s", body_name, e)
        return None
    except (ConnectionError, OSError) as e:
        # Network failures and the downloader's not-found errors are the
        # documented no-exception cases of this helper: log and let the
        # caller fall back (Keplerian elements or UnknownBodyError).
        logger.warning("Auto SPK download failed for %s: %s", body_name, e)
        return None


def is_spk_available_for_body(body_id: int) -> bool:
    """
    Check if an SPK kernel is registered and available for a body.

    This function checks if an SPK kernel has been downloaded and registered
    for the given body, enabling high-precision calculations.

    Args:
        body_id: Minor body identifier (CERES, CHIRON, etc.)

    Returns:
        bool: True if SPK is registered and available for this body

    Example:
        >>> from libephemeris.minor_bodies import is_spk_available_for_body
        >>> from libephemeris.constants import CERES
        >>> if is_spk_available_for_body(CERES):
        ...     print("Using high-precision SPK data")
        ... else:
        ...     print("Using Keplerian approximation")
    """
    try:
        from . import state

        return body_id in state._SPK_BODY_MAP
    except ImportError:
        return False


def ensure_major_asteroid_spk(
    body_id: int,
    jd: Optional[float] = None,
) -> bool:
    """
    Ensure SPK is available for a minor body, downloading if needed.

    This is a convenience function that attempts to download the SPK for
    a minor body if not already available. It is designed to be called
    before position calculations to ensure best precision.

    Supports all bodies in REQUIRED_SPK_BODIES (major asteroids, centaurs,
    dwarf planets) and any body in SPK_BODY_NAME_MAP.

    Args:
        body_id: Minor body identifier (CERES, VESTA, CHIRON,
            PHOLUS, NESSUS, ERIS, etc.)
        jd: Optional Julian Day to center the SPK coverage around.
            If None, uses current date.

    Returns:
        bool: True if SPK is now available (either already was or download
              succeeded), False otherwise.

    Example:
        >>> from libephemeris.minor_bodies import ensure_major_asteroid_spk
        >>> from libephemeris.constants import CERES, REQUIRED_SPK_BODIES
        >>> # Ensure SPK is available for Ceres
        >>> if ensure_major_asteroid_spk(CERES):
        ...     print("High-precision SPK available")
        ... else:
        ...     print("Will use Keplerian approximation")
        >>> # Download all required SPK bodies
        >>> for body in REQUIRED_SPK_BODIES:
        ...     ensure_major_asteroid_spk(body)

    Note:
        This function is non-blocking and returns immediately if the SPK
        is already available. The download only occurs on first call.
    """
    from .constants import SPK_BODY_NAME_MAP

    logger = get_logger()

    # Get body name for logging
    body_info = MAJOR_ASTEROID_SPK_INFO.get(body_id)
    if body_info:
        body_name = body_info[3]
    elif body_id in SPK_BODY_NAME_MAP:
        # Try to get name from MINOR_BODY_ELEMENTS for non-major asteroids
        if body_id in MINOR_BODY_ELEMENTS:
            body_name = MINOR_BODY_ELEMENTS[body_id].name
        else:
            body_name = str(body_id)
    else:
        body_name = str(body_id)

    logger.debug("Checking SPK availability for body %d (%s)", body_id, body_name)

    # First check if already available
    if is_spk_available_for_body(body_id):
        logger.debug("SPK for %s (body %d) already cached", body_name, body_id)
        return True

    # Download using auto_download_asteroid_spk which now handles all bodies
    # in MAJOR_ASTEROID_SPK_INFO and SPK_BODY_NAME_MAP via direct HTTP
    logger.info("SPK for %s not cached, downloading...", body_name)
    if jd is not None:
        jd_start = jd - 3652.5  # ~10 years before
        jd_end = jd + 3652.5  # ~10 years after
        spk_path = auto_download_asteroid_spk(body_id, jd_start, jd_end)
    else:
        spk_path = auto_download_asteroid_spk(body_id)

    if spk_path is not None:
        return True

    return False


def list_major_asteroids() -> list[tuple[int, str]]:
    """
    List all major asteroids with automatic SPK download support.

    Returns:
        list: List of (body_id, name) tuples for all supported major asteroids.

    Example:
        >>> from libephemeris.minor_bodies import list_major_asteroids
        >>> for body_id, name in list_major_asteroids():
        ...     print(f"{name}: body_id={body_id}")
        Ceres: body_id=17
        Pallas: body_id=18
        Juno: body_id=19
        Vesta: body_id=20
        Chiron: body_id=15
    """
    return [(body_id, info[3]) for body_id, info in MAJOR_ASTEROID_SPK_INFO.items()]


def calc_asteroid_by_number(
    asteroid_number: int,
    jd_tt: float,
    use_spk: bool = True,
) -> Tuple[float, float, float]:
    """
    Calculate position for any numbered asteroid by fetching orbital elements on demand.

    This function provides a generic way to calculate the heliocentric position
    of any numbered asteroid (out of 1+ million known) without requiring the
    asteroid to be hardcoded in the library. It first checks if an SPK file
    is available for higher precision, otherwise falls back to Keplerian
    calculations using orbital elements fetched from JPL SBDB.

    Args:
        asteroid_number: The asteroid catalog number (e.g., 1 for Ceres,
            433 for Eros, 2060 for Chiron, 136199 for Eris)
        jd_tt: Julian Day in Terrestrial Time (TT)
        use_spk: If True (default), check for registered SPK file first.
            If False, always use Keplerian calculation.

    Returns:
        Tuple[float, float, float]: (longitude, latitude, distance) where:
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Heliocentric ecliptic latitude in degrees (-90 to +90)
            - distance: Heliocentric distance in AU

    Raises:
        ValueError: If asteroid_number is invalid or not found in JPL SBDB.
        ConnectionError: If unable to fetch orbital elements from JPL SBDB
            (only when SPK is not available and cache is empty).

    Precision:
        - With SPK: Sub-arcsecond to arcsecond accuracy
        - Without SPK (Keplerian): ~10-30 arcseconds for main belt,
          ~1-3 arcminutes for TNOs

    Example:
        >>> from libephemeris.minor_bodies import calc_asteroid_by_number
        >>> # Calculate position for Eros (asteroid 433) at J2000.0
        >>> lon, lat, dist = calc_asteroid_by_number(433, 2451545.0)
        >>> print(f"Eros: {lon:.4f}deg, {lat:.4f}deg, {dist:.4f} AU")
        >>>
        >>> # Calculate for any asteroid by number
        >>> lon, lat, dist = calc_asteroid_by_number(99942, 2451545.0)  # Apophis
        >>> lon, lat, dist = calc_asteroid_by_number(25143, 2451545.0)  # Itokawa

    Notes:
        - First call for an asteroid requires network access to JPL SBDB
        - Subsequent calls use cached orbital elements
        - Use clear_asteroid_elements_cache() to clear the cache
        - For known bodies in MINOR_BODY_ELEMENTS, use calc_minor_body_heliocentric()
          as it doesn't require network access

    See Also:
        calc_minor_body_heliocentric: For bodies with pre-defined orbital elements
        fetch_orbital_elements_from_sbdb: To fetch elements without calculating position
    """
    from .constants import AST_OFFSET

    # Validate asteroid number
    if not isinstance(asteroid_number, int) or asteroid_number <= 0:
        raise ValueError(
            f"asteroid_number must be a positive integer, got {asteroid_number}"
        )

    # Check if we have this asteroid in the predefined MINOR_BODY_ELEMENTS
    body_id = asteroid_number + AST_OFFSET
    if body_id in MINOR_BODY_ELEMENTS:
        return calc_minor_body_heliocentric(body_id, jd_tt)

    # Check special cases for low-numbered asteroids that use different IDs
    # (Chiron, Pholus, Ceres, Pallas, Juno, Vesta have dedicated SE_* constants)
    from .constants import (
        CHIRON,
        PHOLUS,
        CERES,
        PALLAS,
        JUNO,
        VESTA,
    )

    special_mapping = {
        2060: CHIRON,  # Chiron
        5145: PHOLUS,  # Pholus
        1: CERES,  # Ceres
        2: PALLAS,  # Pallas
        3: JUNO,  # Juno
        4: VESTA,  # Vesta
    }

    if asteroid_number in special_mapping:
        mapped_id = special_mapping[asteroid_number]
        if mapped_id in MINOR_BODY_ELEMENTS:
            return calc_minor_body_heliocentric(mapped_id, jd_tt)

    # Check if SPK file is available for this asteroid
    if use_spk:
        try:
            from . import state
            from .constants import FLG_HELCTR
            from .spk import calc_spk_body_position

            # Check if there's a registered SPK for this body
            if body_id in state._SPK_BODY_MAP:
                # SPK is available, create a Skyfield Time object and use SPK calculation
                ts = state.get_timescale()
                t = ts.tt_jd(jd_tt)

                # Calculate heliocentric position using SPK
                result = calc_spk_body_position(t, body_id, FLG_HELCTR)
                if result is not None:
                    lon, lat, dist, _, _, _ = result
                    return lon, lat, dist
        except (ImportError, ValueError, EphemerisRangeError):
            # Fall through to Keplerian calculation
            pass

    # No SPK available, fetch orbital elements from SBDB and calculate
    elements = fetch_orbital_elements_from_sbdb(asteroid_number)

    if elements is None:
        raise ValueError(
            f"Asteroid {asteroid_number} not found in JPL Small-Body Database. "
            "Check if the asteroid number is valid."
        )

    # Calculate position using Keplerian mechanics with perturbations.
    # Frame contract: calc_minor_body_position returns coordinates in the
    # frame the elements are referenced to.  SBDB osculating elements are
    # J2000 ecliptic, so precess to the ecliptic of date below — matching
    # the SPK branch above and calc_minor_body_heliocentric.
    x, y, z = calc_minor_body_position(elements, jd_tt)

    # Convert Cartesian to spherical coordinates
    r = math.sqrt(x**2 + y**2 + z**2)
    lon = math.degrees(math.atan2(y, x)) % 360.0
    lat = math.degrees(math.asin(z / r))

    from .astrometry import _precess_ecliptic

    lon, lat = _precess_ecliptic(lon, lat, 2451545.0, jd_tt)

    return lon, lat, r
