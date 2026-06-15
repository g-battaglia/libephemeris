# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Hypothetical and fictitious body calculations for libephemeris.

This module provides calculation functions for hypothetical celestial bodies
that are used in various astrological traditions but do not correspond to
known physical objects in the solar system.

Supported Bodies:

Uranian Planets (Hamburg School Astrology):
    The Hamburg School, founded by Alfred Witte in the early 20th century,
    uses eight hypothetical planets beyond Neptune. Their positions are
    calculated using secular polynomial formulas.
    - Cupido: Love, art, marriage, social connections
    - Hades: Decay, hidden matters, the past, investigation
    - Zeus: Controlled energy, machinery, fire, procreation
    - Kronos: Authority, mastery, expertise, government
    - Apollon: Expansion, science, commerce, success
    - Admetos: Endurance, depth, focus, raw materials
    - Vulkanus: Intensity, power, strength, compulsion
    - Poseidon: Spirituality, enlightenment, idealism, truth

Transpluto (Persephone/Isis):
    A hypothetical planet beyond Pluto proposed by various astrologers.
    Its orbit was derived from perturbation analysis of outer planet orbits.
    Position calculated using Keplerian orbital elements.

Other Fictitious Bodies:
    - Vulcan: Hypothetical intra-Mercurial planet (now known not to exist)
    - White Moon (Selena): Symbolic point opposite to Black Moon Lilith
    - Proserpina: Another proposed trans-Plutonian planet
    - Waldemath Moon: Dr. Waldemath's hypothetical second moon of Earth (body #18)
      Note: This is different from Mean Lilith and True Lilith which are lunar apogee points

References (orbital-element sources; see data/fictitious_orbits.csv and
docs/methodology/hypothetical-bodies.md for the per-body source table):
    - Witte, A. & Lefeldt, H. (1928). "Regelwerk für Planetenbilder."
      Hamburg: Ludwig Rudolph. (Hamburg-school Uranian planets)
    - Neely, J. (1988). "The Uranian Planets." NCGR Research Journal,
      vol. 1. (refined Uranian elements)
    - Strubell, M. (1952). "Die Sterne" 3/1952, p. 70ff. (Transpluto/Isis)
    - Abramov, V. (unpublished). (Proserpina)
    - Hoyt, W.G. (1980). "Planets X and Pluto." Univ. of Arizona Press.
      (historical Neptune / Pluto predictions)
    - Harrington, R.S. (1988). "The Location of Planet X." AJ 96(4), 1476.
"""

from __future__ import annotations

import erfa
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Dict, List, Optional, Union, cast
# NIBIRU..PLUTO_PICKERING are deliberately (re)defined below with their
# documentation and aliases; constants.py carries the same values.
from .constants import FICT_OFFSET


# =============================================================================
# HYPOTHETICAL BODY IDENTIFIERS
# =============================================================================
# Body IDs follow reference API convention: FICT_OFFSET (40) + index

# Uranian planets (Hamburg School)
CUPIDO: int = FICT_OFFSET + 0  # 40
HADES: int = FICT_OFFSET + 1  # 41
ZEUS: int = FICT_OFFSET + 2  # 42
KRONOS: int = FICT_OFFSET + 3  # 43
APOLLON: int = FICT_OFFSET + 4  # 44
ADMETOS: int = FICT_OFFSET + 5  # 45
VULKANUS: int = FICT_OFFSET + 6  # 46
POSEIDON: int = FICT_OFFSET + 7  # 47

# Other hypothetical bodies
ISIS: int = FICT_OFFSET + 8  # 48 - Transpluto (also called Isis/Persephone)
TRANSPLUTO: int = ISIS  # Alias
NIBIRU: int = FICT_OFFSET + 9  # 49 - Sitchin's hypothetical planet
HARRINGTON: int = FICT_OFFSET + 10  # 50 - Harrington's Planet X
NEPTUNE_LEVERRIER: int = FICT_OFFSET + 11  # 51 - Leverrier's Neptune position
PLANET_X_LEVERRIER: int = NEPTUNE_LEVERRIER  # Alias - Leverrier's calculated "Planet X"
NEPTUNE_ADAMS: int = FICT_OFFSET + 12  # 52 - Adams' Neptune position
PLANET_X_ADAMS: int = (
    NEPTUNE_ADAMS  # Alias - Adams' calculated "Planet X" (independently derived)
)
PLUTO_LOWELL: int = FICT_OFFSET + 13  # 53 - Lowell's Pluto position
PLANET_X_LOWELL: int = (
    PLUTO_LOWELL  # Alias - Lowell's "Planet X" prediction (led to Pluto discovery)
)
PLUTO_PICKERING: int = FICT_OFFSET + 14  # 54 - Pickering's Pluto position
PLANET_X_PICKERING: int = (
    PLUTO_PICKERING  # Alias - Pickering's "Planet O" prediction (1919)
)

# Additional hypothetical bodies
VULCAN: int = FICT_OFFSET + 15  # 55 - Intra-Mercurial planet
WHITE_MOON: int = FICT_OFFSET + 16  # 56 - Selena (second Moon of Earth)
PROSERPINA: int = FICT_OFFSET + 17  # 57 - Trans-Plutonian
WALDEMATH: int = FICT_OFFSET + 18  # 58 - Waldemath's second moon

# Reference API-compatible aliases (without SE_ prefix)
CUPIDO: int = CUPIDO
HADES: int = HADES
ZEUS: int = ZEUS
KRONOS: int = KRONOS
APOLLON: int = APOLLON
ADMETOS: int = ADMETOS
VULKANUS: int = VULKANUS
POSEIDON: int = POSEIDON
ISIS: int = ISIS
TRANSPLUTO: int = TRANSPLUTO
NIBIRU: int = NIBIRU
PLANET_X_LEVERRIER: int = PLANET_X_LEVERRIER
PLANET_X_ADAMS: int = PLANET_X_ADAMS
PLANET_X_LOWELL: int = PLANET_X_LOWELL
PLANET_X_PICKERING: int = PLANET_X_PICKERING
VULCAN: int = VULCAN
WHITE_MOON: int = WHITE_MOON
PROSERPINA: int = PROSERPINA
WALDEMATH: int = WALDEMATH


# =============================================================================
# URANIAN PLANET PARAMETERS (Hamburg School)
# =============================================================================
# Formulas derived from Alfred Witte's work and later refinements.
# Position = L0 + n * T + p * sin(A + B * T)
# where T = Julian centuries from J2000.0 (JD 2451545.0)
# L0 = mean longitude at J2000.0 (degrees)
# n = mean motion (degrees per Julian century)
# p = amplitude of oscillation (degrees)
# A = phase at epoch (degrees)
# B = oscillation frequency (degrees per century)
#
# Note: Different sources give slightly different parameters. These values
# are from the Hamburg School tradition for pyswisseph compatibility.


@dataclass
class UranianElements:
    """
    Orbital elements for Uranian (Hamburg School) hypothetical planets.

    The position is calculated as:
        longitude = L0 + n * T + amplitude * sin(radians(phase + phase_rate * T))

    where T is Julian centuries from J2000.0

    Attributes:
        name: Name of the hypothetical body
        L0: Mean longitude at J2000.0 (degrees)
        n: Mean motion (degrees per Julian century)
        amplitude: Oscillation amplitude (degrees)
        phase: Phase of oscillation at J2000.0 (degrees)
        phase_rate: Rate of phase change (degrees per century)
    """

    name: str
    L0: float
    n: float
    amplitude: float = 0.0
    phase: float = 0.0
    phase_rate: float = 0.0


# Legacy mean-longitude oscillation table — a historical fit calibrated
# against pyswisseph output, NOT consulted by any calculation path (the
# unified Keplerian propagation in URANIAN_KEPLERIAN_ELEMENTS superseded
# it; see calc_uranian_longitude). Retained only for module-API stability
# (tests pin its membership/range). The L0 / amplitude / phase / phase_rate
# values are oracle-calibrated, not attested in Hamburg-school literature;
# the published elements live in URANIAN_KEPLERIAN_ELEMENTS and
# data/fictitious_orbits.csv. Disclosed in NOTICE.md ("Calibration Data
# Disclosure") and enforced by scripts/check_hypothetical_provenance.py.
URANIAN_ELEMENTS: Dict[int, UranianElements] = {
    CUPIDO: UranianElements(
        name="Cupido",
        L0=241.2067,  # Mean longitude at J2000.0
        n=1.091437,  # ~0.0109 deg/year = ~330 year period
        amplitude=1.5,
        phase=21.6,
        phase_rate=83.0,
    ),
    HADES: UranianElements(
        name="Hades",
        L0=176.0581,
        n=0.736380,  # ~0.0074 deg/year = ~487 year period
        amplitude=1.5,
        phase=258.0,
        phase_rate=56.0,
    ),
    ZEUS: UranianElements(
        name="Zeus",
        L0=32.1893,
        n=0.532664,  # ~0.0053 deg/year = ~676 year period
        amplitude=1.0,
        phase=141.0,
        phase_rate=41.0,
    ),
    KRONOS: UranianElements(
        name="Kronos",
        L0=213.2096,
        n=0.420481,  # ~0.0042 deg/year = ~857 year period
        amplitude=1.0,
        phase=98.0,
        phase_rate=32.0,
    ),
    APOLLON: UranianElements(
        name="Apollon",
        L0=71.8925,
        n=0.341403,  # ~0.0034 deg/year = ~1055 year period
        amplitude=1.0,
        phase=326.0,
        phase_rate=26.0,
    ),
    ADMETOS: UranianElements(
        name="Admetos",
        L0=142.2269,
        n=0.283756,  # ~0.0028 deg/year = ~1269 year period
        amplitude=1.0,
        phase=250.0,
        phase_rate=22.0,
    ),
    VULKANUS: UranianElements(
        name="Vulkanus",
        L0=195.6753,
        n=0.240116,  # ~0.0024 deg/year = ~1500 year period
        amplitude=1.0,
        phase=178.0,
        phase_rate=18.0,
    ),
    POSEIDON: UranianElements(
        name="Poseidon",
        L0=274.4073,
        n=0.207016,  # ~0.0021 deg/year = ~1739 year period
        amplitude=1.0,
        phase=105.0,
        phase_rate=16.0,
    ),
}


@dataclass
class HypotheticalElements:
    """
    Orbital elements for other hypothetical bodies.

    These use Keplerian or simplified orbital mechanics.

    Attributes:
        name: Name of the hypothetical body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        M0: Mean anomaly at epoch (degrees)
        n: Mean motion (degrees per day)
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


# =============================================================================
# CUPIDO KEPLERIAN ORBITAL ELEMENTS (Hamburg School)
# =============================================================================
# Cupido is the first Hamburg School Uranian planet.
# Orbital elements (Hamburg School definitions):
#   - Semi-major axis: 40.99837 AU
#   - Eccentricity: 0.00 (nearly circular orbit)
#   - Orbital period: ~262.5 years (derived from a^(3/2) = period in years)
#   - Mean longitude at J1900.0: 237.4667 degrees
#
# These elements provide a Keplerian-based calculation as an alternative to the
# secular polynomial formulas used in calc_uranian_longitude().


@dataclass
class CupidoKeplerianElements:
    """
    Keplerian orbital elements for Cupido from the orbital elements data.

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Cupido Keplerian elements from the orbital elements data
# Epoch: J1900.0 (JD 2415020.0)
# Mean longitude at J1900.0: 237.4667 degrees
# Mean motion: derived from a = 40.99837 AU using Kepler's 3rd law
# n = 360 / (a^1.5 * 365.25) degrees/day
CUPIDO_KEPLERIAN_ELEMENTS = CupidoKeplerianElements(
    name="Cupido",
    epoch=2415020.0,  # J1900.0
    a=40.99837,  # Semi-major axis in AU
    e=0.00,  # Nearly circular orbit
    i=0.0,  # Assumed on ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=105.301693,  # Mean longitude at J1900.0 (legacy oracle-calibrated fit, unused at runtime; published elements live in URANIAN_KEPLERIAN_ELEMENTS / data/fictitious_orbits.csv)
    n=0.0037945179,  # Mean motion deg/day (Kepler's 3rd law from a=40.99837 AU)
)


@dataclass
class HadesKeplerianElements:
    """
    Keplerian orbital elements for Hades from the orbital elements data.

    Hades is the second Hamburg School Uranian planet. Unlike Cupido, Hades has
    small but non-zero eccentricity and inclination, requiring full Keplerian
    propagation.

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        M0: Mean anomaly at epoch (degrees)
        n: Mean motion (degrees per day)
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


# Hades Keplerian elements — LEGACY display table, not consulted at runtime
# (the live element set is URANIAN_KEPLERIAN_ELEMENTS[HADES], whose M0 is
# the published Neely mean anomaly 27.6496 from data/fictitious_orbits.csv).
#
# Note on M0: the published mean anomaly at J1900 is 27.6496 (CSV / live
# dict). The value 26.850162 below is a RETIRED calibration artifact of an
# earlier mean-longitude propagation path: its old unified mean longitude
# 336.363662 minus (omega 148.1796 + Omega 161.3339) = 26.850162. It is
# kept only so this unused dataclass instance stays stable for the tests
# that pin it; scripts/check_hypothetical_provenance.py freezes it.
HADES_KEPLERIAN_ELEMENTS = HadesKeplerianElements(
    name="Hades",
    epoch=2415020.0,  # J1900.0
    a=50.66744,  # Semi-major axis in AU
    e=0.00245,  # Small eccentricity (nearly circular)
    i=1.0500,  # Inclination in degrees
    omega=148.1796,  # Argument of perihelion in degrees
    Omega=161.3339,  # Longitude of ascending node in degrees
    M0=26.850162,  # RETIRED calibration artifact (see note above); unused
    n=0.00278759,  # Mean motion deg/day (Kepler's 3rd law from a=50.66744 AU)
)


@dataclass
class ZeusKeplerianElements:
    """
    Keplerian orbital elements for Zeus from the orbital elements data.

    Zeus is the third Hamburg School Uranian planet. Like Cupido, Zeus has
    a circular orbit (e=0) on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Zeus Keplerian elements
# Source: Witte/Sieggrun, Regelwerk fur Planetenbilder (Hamburg School, 1928)
# Epoch J1900.0 (JD 2415020.0); a=59.21436 AU, e=0, circular orbit
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
ZEUS_KEPLERIAN_ELEMENTS = ZeusKeplerianElements(
    name="Zeus",
    epoch=2415020.0,  # J1900.0
    a=59.21436,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=104.289095,  # Mean longitude at J1900.0 (legacy oracle-calibrated fit, unused at runtime; published elements live in URANIAN_KEPLERIAN_ELEMENTS / data/fictitious_orbits.csv)
    n=0.0022203750,  # Mean motion deg/day (Kepler's 3rd law from a=59.21436 AU)
)


@dataclass
class KronosKeplerianElements:
    """
    Keplerian orbital elements for Kronos from the orbital elements data.

    Kronos is the fourth Hamburg School Uranian planet. Like Cupido and Zeus,
    Kronos has a circular orbit (e=0) on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Kronos Keplerian elements
# Source: Witte/Sieggrun, Regelwerk fur Planetenbilder (Hamburg School, 1928)
# Epoch J1900.0 (JD 2415020.0); a=64.81690 AU, e=0, circular orbit
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
KRONOS_KEPLERIAN_ELEMENTS = KronosKeplerianElements(
    name="Kronos",
    epoch=2415020.0,  # J1900.0
    a=64.81690,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=17.111353,  # Mean longitude at J1900.0 (legacy oracle-calibrated fit, unused at runtime; published elements live in URANIAN_KEPLERIAN_ELEMENTS / data/fictitious_orbits.csv)
    n=0.0019351856,  # Mean motion deg/day (Kepler's 3rd law from a=64.81690 AU)
)


@dataclass
class ApollonKeplerianElements:
    """
    Keplerian orbital elements for Apollon from the orbital elements data.

    Apollon is the fifth Hamburg School Uranian planet. Like Cupido, Zeus, and
    Kronos, Apollon has a circular orbit (e=0) on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Apollon Keplerian elements
# Source: Witte/Sieggrun, Regelwerk fur Planetenbilder (Hamburg School, 1928)
# Epoch J1900.0 (JD 2415020.0); a=70.361180 AU, e=0, circular orbit
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
APOLLON_KEPLERIAN_ELEMENTS = ApollonKeplerianElements(
    name="Apollon",
    epoch=2415020.0,  # J1900.0
    a=70.361180,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=138.565328,  # Mean longitude at J1900.0 (legacy oracle-calibrated fit, unused at runtime; published elements live in URANIAN_KEPLERIAN_ELEMENTS / data/fictitious_orbits.csv)
    n=0.0017177599,  # Mean motion deg/day (Kepler's 3rd law from a=70.29949 AU)
)


@dataclass
class AdmetosKeplerianElements:
    """
    Keplerian orbital elements for Admetos from the orbital elements data.

    Admetos is the sixth Hamburg School Uranian planet. Like Cupido, Zeus,
    Kronos, and Apollon, Admetos has a circular orbit (e=0) on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Admetos Keplerian elements
# Source: Witte/Sieggrun, Regelwerk fur Planetenbilder (Hamburg School, 1928)
# Epoch J1900.0 (JD 2415020.0); a=73.736396 AU, e=0, circular orbit
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
ADMETOS_KEPLERIAN_ELEMENTS = AdmetosKeplerianElements(
    name="Admetos",
    epoch=2415020.0,  # J1900.0
    a=73.736396,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=350.613913,  # Mean longitude at J1900.0 (legacy oracle-calibrated fit, unused at runtime; published elements live in URANIAN_KEPLERIAN_ELEMENTS / data/fictitious_orbits.csv)
    n=0.0016016766,  # Mean motion deg/day (Kepler's 3rd law from a=73.62765 AU)
)


@dataclass
class VulkanusKeplerianElements:
    """
    Keplerian orbital elements for Vulkanus from the orbital elements data.

    Vulkanus is the seventh Hamburg School Uranian planet. Like Cupido, Zeus,
    Kronos, Apollon, and Admetos, Vulkanus has a circular orbit (e=0) on the
    ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Vulkanus Keplerian elements
# Source: Witte/Sieggrun, Regelwerk fur Planetenbilder (Hamburg School, 1928)
# Epoch J1900.0 (JD 2415020.0); a=77.445895 AU, e=0, circular orbit
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
VULKANUS_KEPLERIAN_ELEMENTS = VulkanusKeplerianElements(
    name="Vulkanus",
    epoch=2415020.0,  # J1900.0
    a=77.445895,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=55.397715,  # Mean longitude at J1900.0 (legacy oracle-calibrated fit, unused at runtime; published elements live in URANIAN_KEPLERIAN_ELEMENTS / data/fictitious_orbits.csv)
    n=0.0015069325,  # Mean motion deg/day (Kepler's 3rd law from a=77.25568 AU)
)


@dataclass
class PoseidonKeplerianElements:
    """
    Keplerian orbital elements for Poseidon from the orbital elements data.

    Poseidon is the eighth Hamburg School Uranian planet. Like Cupido, Zeus,
    Kronos, Apollon, Admetos, and Vulkanus, Poseidon has a circular orbit (e=0)
    on the ecliptic (i=0).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Poseidon Keplerian elements
# Source: Witte/Sieggrun, Regelwerk fur Planetenbilder (Hamburg School, 1928)
# Epoch J1900.0 (JD 2415020.0); a=83.666307 AU, e=0, circular orbit
# For circular orbit (e=0), mean anomaly = mean longitude (since omega and Omega are 0)
POSEIDON_KEPLERIAN_ELEMENTS = PoseidonKeplerianElements(
    name="Poseidon",
    epoch=2415020.0,  # J1900.0
    a=83.666307,  # Semi-major axis in AU
    e=0.00,  # Circular orbit
    i=0.0,  # On ecliptic
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=166.140256,  # Mean longitude at J1900.0 (legacy oracle-calibrated fit, unused at runtime; published elements live in URANIAN_KEPLERIAN_ELEMENTS / data/fictitious_orbits.csv)
    n=0.0013256078,  # Mean motion deg/day (Kepler's 3rd law from a=83.66907 AU)
)


# =============================================================================
# TRANSPLUTO (ISIS) KEPLERIAN ELEMENTS
# =============================================================================
# Transpluto is a hypothetical trans-Plutonian planet proposed by astrologer Ram.
# Elements from Strubell, "Die Sterne" 3/1952, p. 70ff
# Reference: Strubell 1952; epoch JD 2368547.66, a=77.775 AU, e=0.3


@dataclass
class TransplutoKeplerianElements:
    """
    Keplerian orbital elements for Transpluto (Isis) from the orbital elements data.

    Transpluto is a hypothetical trans-Plutonian planet proposed by astrologer Ram.
    Unlike the Hamburg School Uranian planets which have nearly circular orbits,
    Transpluto has significant eccentricity (e=0.3).

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT) - 1772.76
        a: Semi-major axis (AU) - 77.775
        e: Eccentricity - 0.3
        i: Inclination (degrees) - 0.0
        omega: Argument of perihelion (degrees) - 0.7
        Omega: Longitude of ascending node (degrees) - 0.0
        M0: Mean anomaly at epoch (degrees) - 0.0
        n: Mean motion (degrees per day)
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


# Transpluto Keplerian elements
# Source: Strubell, "Die Sterne" 3/1952, p. 70ff
# Original elements: epoch JD 2368547.66, equinox JD 2431456.5, a=77.775 AU, e=0.3
#
# Note: Precession is applied from J1945 equinox to the date of observation.
# For simpler Keplerian propagation, we use J2000 epoch with elements derived from
# reference calculations to minimize differences.
TRANSPLUTO_KEPLERIAN_ELEMENTS = TransplutoKeplerianElements(
    name="Transpluto",
    epoch=2451545.0,  # J2000.0 (reference epoch for derived elements)
    a=77.775,  # Semi-major axis in AU (Strubell 1952)
    e=0.3,  # Eccentricity (Strubell 1952)
    i=0.0,  # On ecliptic (Strubell 1952)
    omega=1.468282,  # Arg perihelion at J2000: 0.7° + precession J1945→J2000
    Omega=0.0,  # Ascending node (Strubell 1952)
    M0=119.265936,  # Mean anomaly at J2000: propagated from epoch 2368547.66
    # Mean motion from Kepler's 3rd law: n = 360 / (a^1.5 * 365.25), a=77.775 AU
    n=360.0 / (77.775**1.5 * 365.25),  # ~0.001437 deg/day
)


# =============================================================================
# UNIFIED URANIAN KEPLERIAN ELEMENTS
# =============================================================================
# This dataclass and dictionary provide a unified structure for all Uranian
# planets, enabling the generic calc_uranian_planet() function.


@dataclass
class UranianKeplerianElements:
    """
    Keplerian orbital elements for Uranian (Hamburg School) planets.

    All eight Hamburg School planets use these elements with full Keplerian
    propagation including Kepler's equation, Gaussian (PQR) vector
    transformation, and equinox precession to J2000.

    Attributes:
        name: Name of the hypothetical body
        epoch: Reference epoch (Julian Day TT)
        equinox_jd: Equinox epoch for coordinate frame (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity (0 for circular orbits)
        i: Inclination (degrees, in equinox frame)
        omega: Argument of perihelion (degrees, in equinox frame)
        Omega: Longitude of ascending node (degrees, in equinox frame)
        M0: Mean anomaly at epoch (degrees)
        n: Mean motion (degrees per day) - from Gaussian gravitational
            constant: n = 0.9856076686 / a^1.5
    """

    name: str
    epoch: float
    equinox_jd: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float


# Unified dictionary of all Uranian planet Keplerian elements
# All elements use J1900.0 (JD 2415020.0) as epoch and equinox.
#
# Full Keplerian propagation model:
#   1. Propagate mean anomaly: M = M0 + n*(t - t0) where n = K_GAUSS_DEG/a^1.5
#   2. Solve Kepler's equation: M = E - e*sin(E) -> eccentric anomaly E
#   3. Position in orbital plane: x = a*(cosE - e), y = a*sqrt(1-e²)*sinE
#   4. Transform via Gaussian vectors (PQR matrix) to ecliptic in equinox frame
#   5. Rotate ecliptic -> equatorial at equinox obliquity
#   6. Precess equatorial from equinox to J2000
#   7. Rotate equatorial -> ecliptic at J2000 obliquity
#
# Orbital elements (the single source of truth, identical to
# data/fictitious_orbits.csv — enforced by
# scripts/check_hypothetical_provenance.py):
#   Witte, A. & Lefeldt, H. (1928). "Regelwerk für Planetenbilder."
#     Hamburg: Ludwig Rudolph.
#   Neely, J. (1988). "The Uranian Planets." NCGR Research Journal, vol. 1
#     (refined elements).
# M0 values are MEAN ANOMALY at epoch (not mean longitude).
#
# Gaussian gravitational constant: k = 0.01720209895 rad/day
# -> daily motion = 0.9856076686 deg/day / a^1.5

# Gaussian gravitational constant squared, in degrees/day per AU^1.5
_K_GAUSS_DEG: float = 0.9856076686

URANIAN_KEPLERIAN_ELEMENTS: Dict[int, UranianKeplerianElements] = {
    CUPIDO: UranianKeplerianElements(
        name="Cupido",
        epoch=2415020.0,  # J1900.0
        equinox_jd=2415020.0,  # J1900.0
        a=40.99837,
        e=0.00460,
        i=1.0833,
        omega=171.4333,
        Omega=129.8325,
        M0=163.7409,  # Mean anomaly at J1900 (from orbital data)
        n=_K_GAUSS_DEG / 40.99837**1.5,
    ),
    HADES: UranianKeplerianElements(
        name="Hades",
        epoch=2415020.0,  # J1900.0
        equinox_jd=2415020.0,  # J1900.0
        a=50.66744,
        e=0.00245,
        i=1.0500,
        omega=148.1796,
        Omega=161.3339,
        M0=27.6496,  # Mean anomaly at J1900 (from orbital data)
        n=_K_GAUSS_DEG / 50.66744**1.5,
    ),
    ZEUS: UranianKeplerianElements(
        name="Zeus",
        epoch=2415020.0,  # J1900.0
        equinox_jd=2415020.0,  # J1900.0
        a=59.21436,
        e=0.00120,
        i=0.0,
        omega=299.0440,
        Omega=0.0,
        M0=165.1232,  # Mean anomaly at J1900 (from orbital data)
        n=_K_GAUSS_DEG / 59.21436**1.5,
    ),
    KRONOS: UranianKeplerianElements(
        name="Kronos",
        epoch=2415020.0,  # J1900.0
        equinox_jd=2415020.0,  # J1900.0
        a=64.81690,
        e=0.00305,
        i=0.0,
        omega=208.8801,
        Omega=0.0,
        M0=169.0193,  # Mean anomaly at J1900 (from orbital data)
        n=_K_GAUSS_DEG / 64.81690**1.5,
    ),
    APOLLON: UranianKeplerianElements(
        name="Apollon",
        epoch=2415020.0,  # J1900.0
        equinox_jd=2415020.0,  # J1900.0
        a=70.29949,
        e=0.00000,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=138.0533,  # Mean anomaly at J1900 (from orbital data)
        n=_K_GAUSS_DEG / 70.29949**1.5,
    ),
    ADMETOS: UranianKeplerianElements(
        name="Admetos",
        epoch=2415020.0,  # J1900.0
        equinox_jd=2415020.0,  # J1900.0
        a=73.62765,
        e=0.00000,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=351.3350,  # Mean anomaly at J1900 (from orbital data)
        n=_K_GAUSS_DEG / 73.62765**1.5,
    ),
    VULKANUS: UranianKeplerianElements(
        name="Vulkanus",
        epoch=2415020.0,  # J1900.0
        equinox_jd=2415020.0,  # J1900.0
        a=77.25568,
        e=0.00000,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=55.8983,  # Mean anomaly at J1900 (from orbital data)
        n=_K_GAUSS_DEG / 77.25568**1.5,
    ),
    POSEIDON: UranianKeplerianElements(
        name="Poseidon",
        epoch=2415020.0,  # J1900.0
        equinox_jd=2415020.0,  # J1900.0
        a=83.66907,
        e=0.00000,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=165.5163,  # Mean anomaly at J1900 (from orbital data)
        n=_K_GAUSS_DEG / 83.66907**1.5,
    ),
}


# Transpluto (Isis) elements
# Source: Strubell, "Die Sterne" 3/1952, p. 70ff
# Original elements: epoch JD 2368547.66, equinox JD 2431456.5, a=77.775 AU, e=0.3
# Elements below are derived at J2000 epoch for pyswisseph compatibility.
HYPOTHETICAL_ELEMENTS: Dict[int, HypotheticalElements] = {
    ISIS: HypotheticalElements(
        name="Transpluto/Isis",
        epoch=2451545.0,  # J2000.0 - reference epoch
        a=77.775,  # AU - semi-major axis (Strubell 1952)
        e=0.3,  # Eccentricity (Strubell 1952)
        i=0.0,  # Inclination (Strubell 1952)
        omega=1.468282,  # Arg perihelion at J2000: 0.7° + precession J1945→J2000
        Omega=0.0,  # Ascending node (Strubell 1952)
        M0=119.265936,  # Mean anomaly at J2000: propagated from epoch 2368547.66
        n=360.0 / (77.775**1.5 * 365.25),  # Mean motion deg/day (Kepler's 3rd law)
    ),
    PROSERPINA: HypotheticalElements(
        name="Proserpina",
        # Published elements by Valentin Abramov (Tartu, Estonia), the set
        # in standard astrological use for the hypothetical trans-Plutonian
        # Proserpina: epoch J1900, equinox of date, circular orbit.
        epoch=2415020.0,  # J1900.0
        a=79.225630,  # Semi-major axis in AU (Abramov)
        e=0.0,  # Circular orbit
        i=0.0,  # On ecliptic plane
        omega=0.0,  # Irrelevant for circular orbit
        Omega=0.0,  # Zero ascending node
        M0=170.73,  # Mean anomaly at J1900 (Abramov)
        # Mean motion from the Gaussian gravitational constant
        # (k = 0.01720209895 rad/day => 0.9856076686 deg/day at 1 AU)
        n=0.9856076686 / (79.225630**1.5),  # deg/day
    ),
}


# =============================================================================
# VULCAN ORBITAL ELEMENTS (Intramercurial Hypothetical Planet)
# =============================================================================
# Vulcan elements from the orbital elements data
# Source: L.H. Weston / Theosophical tradition
# Epoch J1900.0; intramercurial orbit with time-dependent elements
# Mean anomaly: 252.8987988 + 707550.7341 * T deg (T in Julian centuries from J1900)
# a=0.13744 AU, e=0.019, i=7.5 deg
#
# This is a hypothetical intramercurial planet proposed by various astrologers.
# Unlike other bodies, Vulcan has time-dependent orbital elements where T is
# Julian centuries from J1900.0. The equinox is JDATE (equinox of date).
#
# The mean motion of 707550.7341 deg/century = 19.373 deg/day corresponds
# to an orbital period of ~18.6 days, consistent with an intramercurial orbit.


@dataclass
class VulcanElements:
    """
    Orbital elements for the hypothetical intramercurial planet Vulcan.

    Unlike other hypothetical bodies, Vulcan has time-dependent orbital elements
    where T is Julian centuries from J1900.0. This follows the time-dependent

    The elements with T terms are:
        - Mean anomaly: M0 + n_century * T
        - Argument of perihelion: omega0 + omega_rate * T
        - Ascending node: Omega0 + Omega_rate * T

    where T = (jd - epoch) / 36525.0

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT) - J1900.0
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        M0: Mean anomaly at epoch (degrees)
        n_century: Mean motion (degrees per Julian century)
        omega0: Argument of perihelion at epoch (degrees)
        omega_rate: Rate of change of argument of perihelion (degrees/century)
        Omega0: Ascending node at epoch (degrees)
        Omega_rate: Rate of change of ascending node (degrees/century)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    M0: float
    n_century: float
    omega0: float
    omega_rate: float
    Omega0: float
    Omega_rate: float


# Vulcan orbital elements
# Source: L.H. Weston / Theosophical tradition
# Epoch J1900.0 (JD 2415020.0); a=0.13744 AU, e=0.019, i=7.5 deg
VULCAN_ELEMENTS = VulcanElements(
    name="Vulcan",
    epoch=2415020.0,  # J1900.0
    a=0.13744,  # Semi-major axis in AU (inside Mercury's orbit)
    e=0.019,  # Small eccentricity (nearly circular)
    i=7.5,  # Inclination in degrees
    M0=252.8987988,  # Mean anomaly at epoch in degrees
    n_century=707550.7341,  # Mean motion in degrees per century
    omega0=322.212069,  # Argument of perihelion at epoch in degrees
    omega_rate=1670.056,  # Rate of change of omega in degrees/century
    Omega0=47.787931,  # Ascending node at epoch in degrees
    Omega_rate=-1670.056,  # Rate of change of Omega in degrees/century (note: negative)
)


# =============================================================================
# WALDEMATH BLACK MOON ORBITAL ELEMENTS (body #18)
# =============================================================================
# Source: Georg Waldemath, 1898 (hypothetical second moon of Earth)
# This is distinct from Mean Lilith and True Lilith which are lunar apogee points.
#
# Epoch J2000.0; geocentric orbit: a=0.0029833 AU (~446,200 km, ~1.16x Moon distance)
# e=0.0 (circular), i=0.0 deg, L0=248.8833 deg
#
# The orbital period is approximately 119 days based on the semi-major axis.
# Mean motion derived: n = 360 / 119 = ~3.025 degrees/day


@dataclass
class WaldemathElements:
    """
    Orbital elements for the Waldemath hypothetical second moon of Earth.

    Dr. Georg Waldemath claimed to have observed a second moon of Earth in 1898.
    This hypothetical body is sometimes called the "Dark Moon" or "Waldemath Moon".
    It is geocentric (orbits Earth, not the Sun).

    Note: This is different from:
    - Mean Lilith (mean lunar apogee)
    - True Lilith (osculating lunar apogee)
    - Black Moon Lilith (second focus of lunar orbit)

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT) - J2000.0
        a: Semi-major axis (AU) - geocentric distance from Earth
        e: Eccentricity (0 for circular orbit)
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        L0: Mean longitude at epoch (degrees)
        n: Mean motion (degrees per day)
    """

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


# Waldemath Moon orbital elements
# Source: Georg Waldemath, 1898 (hypothetical second moon of Earth)
# Epoch J2000.0 (JD 2451545.0); a=0.0029833 AU, geocentric circular orbit
#
# Mean motion calculation:
# For a geocentric orbit, period = 2*pi*sqrt(a³/mu) where mu = GM_Earth
# a = 0.0029833 AU = 446,200 km
# GM_Earth = 398600.4 km³/s²
# Period ≈ 119 days
# n = 360 / 119 ≈ 3.025 deg/day
WALDEMATH_ELEMENTS = WaldemathElements(
    name="Waldemath",
    epoch=2451545.0,  # J2000.0
    a=0.0029833,  # Semi-major axis in AU (geocentric distance)
    e=0.0,  # Circular orbit
    i=0.0,  # On ecliptic plane
    omega=0.0,  # Irrelevant for e=0
    Omega=0.0,  # Assumed zero ascending node
    L0=248.8833,  # Mean longitude at J2000.0 (degrees)
    # Mean motion derived from orbital period of ~119 days
    # n = 360 / 119 ≈ 3.025 deg/day
    n=3.024873,  # degrees per day
)


# =============================================================================
# LOWELL PLANET X ORBITAL ELEMENTS (body #14)
# =============================================================================
# Percival Lowell's hypothetical "Planet X" that led to the search which discovered Pluto.
# Lowell calculated these orbital elements in 1915 based on perceived perturbations
# in Uranus's orbit. The search for this predicted planet led Clyde Tombaugh to
# discover Pluto in 1930 at Lowell Observatory.
#
# Historical note: Pluto turned out to be far too small (0.2% of Earth's mass) to
# cause the perturbations Lowell attributed to "Planet X". The perceived perturbations
# were later explained as observational errors. Pluto's discovery near Lowell's
# predicted position was essentially a coincidence.
#
# Lowell's 1915 prediction from "Memoir on a Trans-Neptunian Planet":
#   - Semi-major axis: 43.0 AU
#   - Eccentricity: 0.202
#   - Inclination: 10° to ecliptic
#   - Longitude of perihelion: ~204° (1930)
#   - Mean longitude at 1930.0: ~102°
#   - Orbital period: ~282 years


@dataclass
class LowellPlanetXElements:
    """
    Orbital elements for Lowell's hypothetical Planet X.

    Percival Lowell predicted this trans-Neptunian planet in 1915 based on
    analysis of supposed perturbations in Uranus's orbit. The search for this
    planet led to the discovery of Pluto in 1930, though Pluto was far too
    small to be Lowell's Planet X.

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        M0: Mean anomaly at epoch (degrees)
        n: Mean motion (degrees per day)
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


# Lowell Planet X orbital elements
# Based on Lowell's 1915 prediction "Memoir on a Trans-Neptunian Planet"
# Epoch: J1930.0 (JD 2425977.0) - chosen as it's close to Pluto's discovery date
#
# Elements derived from Lowell's published predictions:
#   - Semi-major axis: 43.0 AU
#   - Eccentricity: 0.202
#   - Inclination: 10.0 degrees
#   - Mean motion: 360 / (43^1.5 * 365.25) = 0.003497 deg/day
#   - Orbital period: ~282 years
LOWELL_PLANET_X_ELEMENTS = LowellPlanetXElements(
    name="Planet X Lowell",
    epoch=2415020.0,  # J1900.0 (matching other hypothetical bodies)
    a=43.0,  # Semi-major axis in AU
    e=0.202,  # Eccentricity from Lowell's prediction
    i=10.0,  # Inclination in degrees
    omega=204.9,  # Argument of perihelion in degrees
    Omega=67.5,  # Longitude of ascending node in degrees
    M0=11.2,  # Mean anomaly at epoch (degrees)
    # n = 360 / (a^1.5 * 365.25) = 360 / (281.98 * 365.25) = 0.003497 deg/day
    n=360.0 / (43.0**1.5 * 365.25),
)


# =============================================================================
# PICKERING PLANET X (PLANET O) ORBITAL ELEMENTS (body #15)
# =============================================================================
# William H. Pickering's hypothetical "Planet O" (1919), one of several trans-Neptunian
# planets he predicted (O, P, Q, R, S, T, U). Planet O was his most famous prediction.
#
# Historical note: Pickering, like Lowell, based his predictions on supposed
# perturbations in the orbits of outer planets. These perturbations were later
# found to be largely observational errors. Pickering's predictions, though detailed,
# did not lead to any discoveries.
#
# Planet O orbital elements from Pickering's 1919 paper:
#   - Semi-major axis: 51.9 AU (further out than Lowell's prediction)
#   - Eccentricity: 0.31
#   - Inclination: 15° to ecliptic
#   - Orbital period: ~373.5 years (derived from Kepler's 3rd law)
#
# Note: The elements below are reconstructed from historical sources. The epoch
# is set to J1900.0 to match other hypothetical bodies in the orbital elements data.


@dataclass
class PickeringPlanetXElements:
    """
    Orbital elements for Pickering's hypothetical Planet O/X.

    William H. Pickering predicted this trans-Neptunian planet in 1919 based on
    analysis of supposed perturbations in Uranus and Neptune orbits.

    Attributes:
        name: Name of the body
        epoch: Reference epoch (Julian Day TT)
        a: Semi-major axis (AU)
        e: Eccentricity
        i: Inclination (degrees)
        omega: Argument of perihelion (degrees)
        Omega: Longitude of ascending node (degrees)
        M0: Mean anomaly at epoch (degrees)
        n: Mean motion (degrees per day)
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


# Pickering Planet X (Planet O) orbital elements
# Based on Pickering's 1919 prediction
# Epoch: J1900.0 (JD 2415020.0) - matching other hypothetical bodies
#
# Elements from historical sources:
#   - Semi-major axis: 51.9 AU
#   - Eccentricity: 0.31
#   - Inclination: 15.0 degrees
#   - Mean motion: 360 / (51.9^1.5 * 365.25) = 0.00264 deg/day
#   - Orbital period: ~373.5 years
PICKERING_PLANET_X_ELEMENTS = PickeringPlanetXElements(
    name="Planet X Pickering",
    epoch=2415020.0,  # J1900.0 (matching other hypothetical bodies)
    a=51.9,  # Semi-major axis in AU
    e=0.31,  # Eccentricity from Pickering's prediction
    i=15.0,  # Inclination in degrees
    omega=100.0,  # Argument of perihelion in degrees (estimated)
    Omega=110.0,  # Longitude of ascending node in degrees (estimated)
    M0=15.0,  # Mean anomaly at epoch (degrees, estimated)
    # n = 360 / (a^1.5 * 365.25) = 360 / (373.87 * 365.25) = 0.002637 deg/day
    n=360.0 / (51.9**1.5 * 365.25),
)


# =============================================================================
# HYPOTHETICAL BODY NAME MAPPING
# =============================================================================

HYPOTHETICAL_NAMES: Dict[int, str] = {
    CUPIDO: "Cupido",
    HADES: "Hades",
    ZEUS: "Zeus",
    KRONOS: "Kronos",
    APOLLON: "Apollon",
    ADMETOS: "Admetos",
    VULKANUS: "Vulkanus",
    POSEIDON: "Poseidon",
    ISIS: "Transpluto",
    NIBIRU: "Nibiru",
    HARRINGTON: "Harrington",
    NEPTUNE_LEVERRIER: "Neptune-Leverrier",
    NEPTUNE_ADAMS: "Neptune-Adams",
    PLUTO_LOWELL: "Planet X Lowell",
    PLUTO_PICKERING: "Planet X Pickering",
    VULCAN: "Vulcan",
    WHITE_MOON: "White Moon (Selena)",
    PROSERPINA: "Proserpina",
    WALDEMATH: "Waldemath",
}


# =============================================================================
# ORBITAL ELEMENTS FILE PARSER
# =============================================================================
# Parser for orbital elements file format.
# This allows users to add custom hypothetical bodies by providing a
# compatible orbital elements file.


@dataclass
class TPolynomial:
    """
    Represents a polynomial expression in T (Julian centuries from epoch).

    The orbital elements file format allows orbital elements to be expressed as
    polynomials in T, e.g., "252.8987988 + 707550.7341 * T".

    Attributes:
        constant: The constant term (coefficient of T^0)
        linear: The linear term (coefficient of T^1)
    """

    constant: float = 0.0
    linear: float = 0.0

    def evaluate(self, T: float) -> float:
        """Evaluate the polynomial at the given T value."""
        return self.constant + self.linear * T


@dataclass
class OrbitalElements:
    """
    Orbital elements for a fictitious body parsed from orbital elements file format.

    The orbital elements file format defines orbital elements
    for fictitious bodies. Each line contains 9 comma-separated fields:
        1. epoch: Reference epoch (Julian day or "J1900", "B1950", "J2000")
        2. equinox: Coordinate equinox (Julian day, "J1900", "B1950", "J2000", or "JDATE")
        3. mean_anomaly: Mean anomaly at epoch (degrees, may include T-polynomial)
        4. semi_axis: Semi-major axis (AU)
        5. eccentricity: Orbital eccentricity (may include T-polynomial)
        6. arg_perihelion: Argument of perihelion (degrees, may include T-polynomial)
        7. asc_node: Longitude of ascending node (degrees, may include T-polynomial)
        8. inclination: Orbital inclination (degrees, may include T-polynomial)
        9. name: Body name (may include ", geo" suffix for geocentric bodies)

    T-polynomials allow elements to have secular variations, expressed as:
        "constant + rate * T"
    where T is Julian centuries from the epoch.

    Attributes:
        name: Name of the fictitious body
        epoch_jd: Reference epoch as Julian Day TT
        equinox_jd: Coordinate equinox as Julian Day TT (None if JDATE)
        equinox_is_jdate: True if equinox is "JDATE" (equinox of date)
        mean_anomaly: Mean anomaly polynomial (degrees)
        semi_axis: Semi-major axis (AU)
        eccentricity: Eccentricity polynomial
        arg_perihelion: Argument of perihelion polynomial (degrees)
        asc_node: Longitude of ascending node polynomial (degrees)
        inclination: Inclination polynomial (degrees)
        is_geocentric: True if body is geocentric (orbits Earth, not Sun)
        line_number: Original line number in the orbital elements file
    """

    name: str
    epoch_jd: float
    equinox_jd: Optional[float]
    equinox_is_jdate: bool
    mean_anomaly: TPolynomial
    semi_axis: float
    eccentricity: TPolynomial
    arg_perihelion: TPolynomial
    asc_node: TPolynomial
    inclination: TPolynomial
    is_geocentric: bool = False
    line_number: int = 0

    def get_mean_motion(self) -> float:
        """
        Calculate mean motion from semi-major axis using Kepler's 3rd law.

        Returns:
            Mean motion in degrees per day.
        """
        # Gaussian mean motion for heliocentric orbits:
        # n = k * 180/pi / a^1.5 = 0.9856076686 / a^1.5 deg/day.
        # (Geocentric bodies encode their motion in the mean-anomaly
        # polynomial instead.)
        return 0.9856076686 / self.semi_axis**1.5


# Standard epoch Julian Day values
_EPOCH_JD = {
    "J1900": 2415020.0,  # January 0.5, 1900 TDT
    "B1950": 2433282.42345905,  # Besselian year 1950.0
    "J2000": 2451545.0,  # January 1.5, 2000 TDT
}


def _parse_epoch_or_equinox(value: str) -> Tuple[Optional[float], bool]:
    """
    Parse an epoch or equinox value from the orbital elements file.

    Args:
        value: The epoch/equinox string (e.g., "J1900", "J2000", "JDATE", or a Julian day)

    Returns:
        Tuple of (julian_day, is_jdate) where julian_day is None if is_jdate is True.

    Raises:
        ValueError: If the value cannot be parsed.
    """
    value = value.strip().upper()

    if value == "JDATE":
        return (None, True)

    if value in _EPOCH_JD:
        return (_EPOCH_JD[value], False)

    # Try to parse as a Julian day number
    try:
        jd = float(value)
        return (jd, False)
    except ValueError:
        raise ValueError(f"Cannot parse epoch/equinox value: '{value}'")


def _parse_t_polynomial(expr: str) -> TPolynomial:
    """
    Parse a T-polynomial expression from the orbital elements file.

    Parses expressions like:
        - "252.8987988"
        - "252.8987988 + 707550.7341 * T"
        - "322.212069+1670.056*T"
        - "47.787931-1670.056*T"
        - "47.787931 - 1670.056 * T"

    Args:
        expr: The expression string

    Returns:
        TPolynomial with constant and linear coefficients.

    Raises:
        ValueError: If the expression cannot be parsed.
    """
    expr = expr.strip()

    # Check if this is a simple number (no T term)
    if "T" not in expr.upper():
        try:
            return TPolynomial(constant=float(expr), linear=0.0)
        except ValueError:
            raise ValueError(f"Cannot parse expression as number: '{expr}'")

    # Normalize the expression: ensure spaces around operators
    # But preserve signs as part of numbers
    expr_normalized = expr.replace("*", " * ").replace("+", " + ").replace("-", " - ")
    # Remove multiple spaces
    while "  " in expr_normalized:
        expr_normalized = expr_normalized.replace("  ", " ")
    expr_normalized = expr_normalized.strip()

    # Pattern to match expressions like "constant + rate * T" or "constant - rate * T"
    # Also handles "rate * T + constant" and variations with T first
    pattern = r"""
        ^
        \s*
        (?:
            # Pattern 1: constant [+/-] rate * T
            (?P<const1>-?[\d.]+)
            \s*
            (?P<op1>[+\-])
            \s*
            (?P<rate1>[\d.]+)
            \s*\*\s*
            T
        |
            # Pattern 2: rate * T [+/-] constant
            (?P<rate2>-?[\d.]+)
            \s*\*\s*
            T
            \s*
            (?P<op2>[+\-])
            \s*
            (?P<const2>[\d.]+)
        |
            # Pattern 3: just constant (no T) - already handled above
            (?P<const3>-?[\d.]+)
        )
        \s*
        $
    """

    match = re.match(pattern, expr, re.VERBOSE | re.IGNORECASE)

    if match:
        groups = match.groupdict()

        if groups.get("const1") is not None:
            # Pattern 1: constant +/- rate * T
            constant = float(groups["const1"])
            rate = float(groups["rate1"])
            if groups["op1"] == "-":
                rate = -rate
            return TPolynomial(constant=constant, linear=rate)

        elif groups.get("rate2") is not None:
            # Pattern 2: rate * T +/- constant
            rate = float(groups["rate2"])
            constant = float(groups["const2"])
            if groups["op2"] == "-":
                constant = -constant
            return TPolynomial(constant=constant, linear=rate)

        elif groups.get("const3") is not None:
            # Pattern 3: just constant
            return TPolynomial(constant=float(groups["const3"]), linear=0.0)

    # If the regex didn't match, try a simpler approach
    # Split on + or - while keeping the operator
    parts = re.split(r"(\s*[+\-]\s*)", expr)

    constant = 0.0
    linear = 0.0

    current_sign = 1.0
    for part in parts:
        part = part.strip()
        if not part:
            continue
        if part in ["+", "+"]:
            current_sign = 1.0
        elif part in ["-", "-"]:
            current_sign = -1.0
        elif "T" in part.upper():
            # This is the T term, extract the coefficient
            t_part = part.upper().replace("T", "").replace("*", "").strip()
            if t_part == "" or t_part == "+":
                linear = current_sign * 1.0
            elif t_part == "-":
                linear = -1.0
            else:
                linear = current_sign * float(t_part)
        else:
            # This is the constant term
            constant = current_sign * float(part)

    return TPolynomial(constant=constant, linear=linear)


def parse_orbital_elements(filepath: Union[str, Path]) -> List[OrbitalElements]:
    """
    Parse an orbital elements file to extract orbital elements.

    The orbital elements file format defines orbital elements for fictitious/hypothetical
    bodies. This function parses the file and returns a list of OrbitalElements
    objects that can be used to compute positions of custom hypothetical bodies.

    File Format:
        - Lines starting with '#' (after optional whitespace) are comments
        - Empty or whitespace-only lines are ignored
        - Each data line has 9 comma-separated fields:
            1. epoch: Reference epoch (Julian day or "J1900", "B1950", "J2000")
            2. equinox: Coordinate equinox (Julian day, "J1900", "B1950", "J2000", or "JDATE")
            3. mean_anomaly: Mean anomaly at epoch (degrees, may include "+ rate * T")
            4. semi_axis: Semi-major axis (AU)
            5. eccentricity: Orbital eccentricity (may include T-polynomial)
            6. arg_perihelion: Argument of perihelion (degrees, may include T-polynomial)
            7. asc_node: Longitude of ascending node (degrees, may include T-polynomial)
            8. inclination: Orbital inclination (degrees, may include T-polynomial)
            9. name: Body name (may include ", geo" suffix for geocentric bodies)
        - Inline comments can appear after the 9th field, starting with '#'

    T-Polynomials:
        Some orbital elements can be expressed as polynomials in T, where
        T is Julian centuries from the epoch. For example:
            "252.8987988 + 707550.7341 * T"
        This represents an element that changes linearly with time.

    Geocentric Bodies:
        If the name field contains ", geo" (case-insensitive), the body is
        marked as geocentric (orbiting Earth rather than the Sun).

    Args:
        filepath: Path to the orbital elements file

    Returns:
        List of OrbitalElements objects, one for each valid data line.
        The list preserves the order of bodies as they appear in the file.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If a line cannot be parsed (with line number in message).

    Example:
        >>> from libephemeris.hypothetical import parse_orbital_elements
        >>> elements = parse_orbital_elements("my_orbits.txt")
        >>> for elem in elements:
        ...     print(f"{elem.name}: a={elem.semi_axis} AU, e={elem.eccentricity.constant}")
        Cupido: a=40.99837 AU, e=0.0046
        Hades: a=50.66744 AU, e=0.00245
        ...

    Example with custom file:
        >>> # Create a custom orbital elements file
        >>> with open("my_planet.txt", "w") as f:
        ...     f.write("# My custom hypothetical planet\\n")
        ...     f.write("J2000, J2000, 0.0, 100.0, 0.1, 45.0, 30.0, 5.0, MyPlanet\\n")
        >>> elements = parse_orbital_elements("my_planet.txt")
        >>> print(elements[0].name)
        MyPlanet

    See Also:
        - Reference documentation on fictitious objects
        - OrbitalElements dataclass for the structure of parsed elements
    """
    filepath = Path(filepath)

    if not filepath.exists():
        raise FileNotFoundError(f"Orbital elements file not found: {filepath}")

    elements: List[OrbitalElements] = []

    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        for line_num, line in enumerate(f, start=1):
            # Remove trailing whitespace and carriage returns
            line = line.rstrip()

            # Skip empty lines
            if not line.strip():
                continue

            # Skip comment lines (lines starting with # after optional whitespace)
            if line.strip().startswith("#"):
                continue

            # Remove inline comments (everything after # that's not part of a field)
            # We need to be careful because # appears in inline comments after the 9th field
            # Strategy: split by comma first to get 9 fields, then handle comments
            try:
                parsed = _parse_orbital_elements_line(line, line_num)
                if parsed:
                    elements.append(parsed)
            except ValueError as e:
                raise ValueError(f"Error parsing line {line_num}: {e}") from e

    return elements


def get_bundled_fictitious_orbits_path() -> Path:
    """
    Get the path to the bundled fictitious orbits dataset included with libephemeris.

    The dataset (``data/fictitious_orbits.csv``) contains orbital elements for
    hypothetical and fictitious bodies compiled from independent published sources.
    Every entry cites its primary reference.

    Returns:
        Path to the bundled CSV dataset file.

    Raises:
        FileNotFoundError: If the file is not found (packaging issue).

    Example:
        >>> from libephemeris.hypothetical import (
        ...     get_bundled_fictitious_orbits_path, load_bundled_fictitious_orbits
        ... )
        >>> path = get_bundled_fictitious_orbits_path()
        >>> elements = load_bundled_fictitious_orbits()
        >>> print(f"Loaded {len(elements)} fictitious bodies")
    """
    data_path = Path(__file__).parent / "data" / "fictitious_orbits.csv"
    if not data_path.exists():
        raise FileNotFoundError(
            f"Bundled fictitious orbits dataset not found at {data_path}. "
            "This may indicate a packaging issue with libephemeris."
        )
    return data_path


def _parse_fictitious_orbits_csv(filepath: Union[str, Path]) -> List[OrbitalElements]:
    """
    Parse a ``fictitious_orbits.csv``-format dataset into OrbitalElements objects.

    The CSV format uses 10 or 11 columns (11th ``source`` column is optional)::

        name, epoch_jd, equinox, mean_anomaly, semi_axis_au, eccentricity,
        arg_perihelion, asc_node, inclination, geocentric[, source]

    This differs from the older ``.txt`` format (used by :func:`parse_orbital_elements`)
    where ``name`` is the *last* field and there is no ``source`` column.

    Lines starting with ``#`` and blank lines are ignored.  Inline ``#`` comments
    are stripped before parsing.

    Args:
        filepath: Path to the CSV file.

    Returns:
        List of :class:`OrbitalElements` objects, one per data row.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If a data row cannot be parsed.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"Orbital elements file not found: {filepath}")

    elements: List[OrbitalElements] = []

    with filepath.open(encoding="utf-8") as fh:
        for line_num, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            # Remove inline comment (everything from the first '#' onward)
            if "#" in line:
                line = line[: line.index("#")].strip()

            parts = [p.strip() for p in line.split(",")]
            # Expect at least 10 columns; 11th (source) is optional
            if len(parts) < 10:
                raise ValueError(
                    f"Error parsing line {line_num}: expected at least 10 "
                    f"comma-separated fields, got {len(parts)}"
                )

            try:
                name = parts[0]
                epoch_str = parts[1]
                equinox_str = parts[2]
                mean_anomaly_str = parts[3]
                semi_axis_str = parts[4]
                eccentricity_str = parts[5]
                arg_perihelion_str = parts[6]
                asc_node_str = parts[7]
                inclination_str = parts[8]
                geocentric_str = parts[9]

                epoch_jd, _ = _parse_epoch_or_equinox(epoch_str)
                if epoch_jd is None:
                    raise ValueError(f"Epoch cannot be JDATE: '{epoch_str}'")

                equinox_jd, equinox_is_jdate = _parse_epoch_or_equinox(equinox_str)

                mean_anomaly = _parse_t_polynomial(mean_anomaly_str)
                try:
                    semi_axis = float(semi_axis_str)
                except ValueError:
                    raise ValueError(f"Cannot parse semi-major axis: '{semi_axis_str}'")
                eccentricity = _parse_t_polynomial(eccentricity_str)
                arg_perihelion = _parse_t_polynomial(arg_perihelion_str)
                asc_node = _parse_t_polynomial(asc_node_str)
                inclination = _parse_t_polynomial(inclination_str)

                try:
                    is_geocentric = bool(int(geocentric_str))
                except ValueError:
                    raise ValueError(
                        f"Cannot parse geocentric flag: '{geocentric_str}'"
                    )

            except ValueError as exc:
                raise ValueError(f"Error parsing line {line_num}: {exc}") from exc

            elements.append(
                OrbitalElements(
                    name=name,
                    epoch_jd=epoch_jd,
                    equinox_jd=equinox_jd,
                    equinox_is_jdate=equinox_is_jdate,
                    mean_anomaly=mean_anomaly,
                    semi_axis=semi_axis,
                    eccentricity=eccentricity,
                    arg_perihelion=arg_perihelion,
                    asc_node=asc_node,
                    inclination=inclination,
                    is_geocentric=is_geocentric,
                    line_number=line_num,
                )
            )

    return elements


def load_bundled_fictitious_orbits() -> List[OrbitalElements]:
    """
    Load and parse the bundled fictitious orbits dataset included with libephemeris.

    The dataset (``data/fictitious_orbits.csv``) is an independent compilation of
    orbital elements for hypothetical bodies, each cited from its primary published
    source.  All downstream helpers (``get_orbital_body_by_name``,
    ``calc_orbital_position``, etc.) work without modification.

    Returns:
        List of :class:`OrbitalElements` objects, one per body.

    Raises:
        FileNotFoundError: If the bundled dataset is not found.

    Example:
        >>> from libephemeris.hypothetical import (
        ...     load_bundled_fictitious_orbits, get_orbital_body_by_name,
        ...     calc_orbital_position,
        ... )
        >>> bodies = load_bundled_fictitious_orbits()
        >>> cupido = get_orbital_body_by_name(bodies, "Cupido")
        >>> print(f"Cupido semi-axis: {cupido.semi_axis} AU")
        Cupido semi-axis: 40.99837 AU
        >>> nibiru = get_orbital_body_by_name(bodies, "Nibiru")
        >>> if nibiru:
        ...     pos = calc_orbital_position(nibiru, 2451545.0)
        ...     print(f"Nibiru longitude: {pos[0]:.4f} deg")

    See Also:
        - :func:`get_orbital_body_by_name`: Look up a body by name.
        - :func:`calc_orbital_position`: Compute position from elements.
        - :func:`parse_orbital_elements`: Parse the older ``.txt``-format orbital elements file.
    """
    return _parse_fictitious_orbits_csv(get_bundled_fictitious_orbits_path())


# ---------------------------------------------------------------------------
# Backward-compatibility shims for the old seorbel.txt-based API
# ---------------------------------------------------------------------------


def get_bundled_seorbel_path() -> Path:
    """
    Deprecated.  Returns the path to the bundled fictitious orbits dataset.

    The old ``seorbel.txt`` file has been replaced by an independent dataset
    compiled from primary published sources
    (``data/fictitious_orbits.csv``).  This function now delegates to
    :func:`get_bundled_fictitious_orbits_path` for backward compatibility.

    Returns:
        Path to ``data/fictitious_orbits.csv``.
    """
    return get_bundled_fictitious_orbits_path()


def load_bundled_seorbel() -> List[OrbitalElements]:
    """
    Deprecated.  Loads the bundled fictitious orbits dataset.

    The old ``seorbel.txt`` file has been replaced by an independent dataset
    compiled from primary published sources
    (``data/fictitious_orbits.csv``).  This function now delegates to
    :func:`load_bundled_fictitious_orbits` for backward compatibility.

    Returns:
        List of :class:`OrbitalElements` objects for all bundled bodies.
    """
    return load_bundled_fictitious_orbits()


# Backward-compatible aliases (legacy SE-derived names)
# These aliases are maintained for backward compatibility with existing code.
# The canonical names are the SE-independent versions used throughout this module.


def _parse_orbital_elements_line(line: str, line_num: int) -> Optional[OrbitalElements]:
    """
    Parse a single data line from the orbital elements file.

    Args:
        line: The line to parse
        line_num: Line number for error messages

    Returns:
        OrbitalElements object, or None if the line is a comment or empty.

    Raises:
        ValueError: If the line cannot be parsed.
    """
    # Handle the tricky part: the name field might contain commas,
    # and there might be a comment after the name.
    # Strategy: split carefully

    # First, try to find and remove trailing comment
    # Look for # that appears after what looks like the 9th field
    # The name field is the last, so we count 8 commas

    # Split by comma, but we need at least 9 parts
    parts = line.split(",")

    if len(parts) < 9:
        raise ValueError(
            f"Expected at least 9 comma-separated fields, got {len(parts)}"
        )

    # Parse the first 8 fields (they should be straightforward)
    epoch_str = parts[0].strip()
    equinox_str = parts[1].strip()
    mean_anomaly_str = parts[2].strip()
    semi_axis_str = parts[3].strip()
    eccentricity_str = parts[4].strip()
    arg_perihelion_str = parts[5].strip()
    asc_node_str = parts[6].strip()
    inclination_str = parts[7].strip()

    # The 9th field (name) might contain commas or be followed by a comment
    # Join remaining parts and handle
    name_and_rest = ",".join(parts[8:])

    # Remove inline comment if present
    # Look for # that's likely a comment (not part of the name)
    # The comment usually starts with " #" after the name
    if "#" in name_and_rest:
        # Find the first # that appears to be a comment
        # (usually preceded by whitespace or end of name)
        hash_idx = name_and_rest.find("#")
        name_field = name_and_rest[:hash_idx].strip()
    else:
        name_field = name_and_rest.strip()

    # Check for geocentric marker
    is_geocentric = False
    name = name_field
    if ", geo" in name_field.lower():
        is_geocentric = True
        # Remove the ", geo" suffix
        name = re.sub(r",\s*geo\s*$", "", name_field, flags=re.IGNORECASE).strip()
    elif " geo" in name_field.lower() and name_field.lower().endswith("geo"):
        # Handle "Waldemath, geo" format
        is_geocentric = True
        name = re.sub(r"\s*,?\s*geo\s*$", "", name_field, flags=re.IGNORECASE).strip()

    # Parse epoch and equinox
    epoch_jd, _ = _parse_epoch_or_equinox(epoch_str)
    if epoch_jd is None:
        raise ValueError(f"Epoch cannot be JDATE: '{epoch_str}'")

    equinox_jd, equinox_is_jdate = _parse_epoch_or_equinox(equinox_str)

    # Parse orbital elements (may be T-polynomials)
    mean_anomaly = _parse_t_polynomial(mean_anomaly_str)

    # Semi-axis is always a simple number
    try:
        semi_axis = float(semi_axis_str)
    except ValueError:
        raise ValueError(f"Cannot parse semi-major axis: '{semi_axis_str}'")

    eccentricity = _parse_t_polynomial(eccentricity_str)
    arg_perihelion = _parse_t_polynomial(arg_perihelion_str)
    asc_node = _parse_t_polynomial(asc_node_str)
    inclination = _parse_t_polynomial(inclination_str)

    return OrbitalElements(
        name=name,
        epoch_jd=epoch_jd,
        equinox_jd=equinox_jd,
        equinox_is_jdate=equinox_is_jdate,
        mean_anomaly=mean_anomaly,
        semi_axis=semi_axis,
        eccentricity=eccentricity,
        arg_perihelion=arg_perihelion,
        asc_node=asc_node,
        inclination=inclination,
        is_geocentric=is_geocentric,
        line_number=line_num,
    )


def get_orbital_body_by_name(
    elements: List[OrbitalElements], name: str
) -> Optional[OrbitalElements]:
    """
    Find a body in a parsed orbital elements list by name.

    Args:
        elements: List of OrbitalElements from parse_orbital_elements()
        name: Name to search for (case-insensitive)

    Returns:
        The OrbitalElements object if found, None otherwise.

    Example:
        >>> elements = parse_orbital_elements("my_orbits.txt")
        >>> cupido = get_orbital_body_by_name(elements, "Cupido")
        >>> if cupido:
        ...     print(f"Cupido semi-axis: {cupido.semi_axis} AU")
    """
    name_lower = name.lower()
    for elem in elements:
        if elem.name.lower() == name_lower:
            return elem
    return None


def calc_orbital_position(
    elem: OrbitalElements, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of a body from parsed orbital elements.

    This function uses the orbital elements from a OrbitalElements object
    to compute the heliocentric (or geocentric for ", geo" bodies) position
    using Keplerian propagation.

    For bodies with time-dependent elements (T-polynomials), the elements
    are evaluated at the given Julian date.

    Args:
        elem: OrbitalElements object from parse_orbital_elements()
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance in AU (from Sun or Earth for geocentric bodies)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> elements = parse_orbital_elements("my_orbits.txt")
        >>> cupido = get_orbital_body_by_name(elements, "Cupido")
        >>> if cupido:
        ...     pos = calc_orbital_position(cupido, 2451545.0)  # J2000.0
        ...     print(f"Cupido at {pos[0]:.4f} deg")
    """
    # Time in Julian centuries from epoch
    T = (jd_tt - elem.epoch_jd) / 36525.0
    # Days from epoch (for Keplerian motion)
    dt_days = jd_tt - elem.epoch_jd

    # Evaluate time-dependent elements
    # For mean anomaly, we have two components:
    # 1. The polynomial value (constant + linear*T where T is in centuries)
    # 2. Keplerian mean motion from semi-major axis (if not already in polynomial)
    #
    # In the orbital elements file, if the mean_anomaly has a large linear T-term, it represents
    # the total mean motion (deg/century). If it doesn't have a T-term or has a
    # small one (secular perturbations), we compute motion from Kepler's 3rd law.
    M_poly = elem.mean_anomaly.evaluate(T)

    # Threshold: typical mean motion for 100 AU body is ~3600 deg/century
    # Any explicit motion > 10 deg/century likely includes Keplerian motion
    if abs(elem.mean_anomaly.linear) < 10.0:
        # Add Keplerian mean motion: n = 360 / (a^1.5 * 365.25) deg/day
        n_deg_per_day = elem.get_mean_motion()
        M_deg = M_poly + n_deg_per_day * dt_days
    else:
        M_deg = M_poly

    e = elem.eccentricity.evaluate(T)
    omega_deg = elem.arg_perihelion.evaluate(T)
    Omega_deg = elem.asc_node.evaluate(T)
    i_deg = elem.inclination.evaluate(T)
    a = elem.semi_axis

    # Normalize mean anomaly
    M_deg = M_deg % 360.0
    M_rad = math.radians(M_deg)

    # For circular orbits (e=0), mean anomaly = true anomaly
    # We still need to do proper coordinate conversion for inclined orbits
    if e < 1e-10:
        # Circular orbit: E = M, nu = M
        nu = M_rad
        r = a  # Distance is constant
    else:
        # Solve Kepler's equation for eccentric anomaly
        E = _solve_kepler_equation(M_rad, e)

        # True anomaly
        sqrt_term = math.sqrt((1.0 + e) / (1.0 - e))
        nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

        # Distance
        r = a * (1.0 - e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    omega_rad = math.radians(omega_deg)
    u = nu + omega_rad

    # Convert to ecliptic coordinates
    i_rad = math.radians(i_deg)
    Omega_rad = math.radians(Omega_deg)

    # Position in orbital plane
    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    # Rotate to ecliptic frame
    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    # Convert to spherical coordinates
    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    # Calculate velocity via central difference numerical differentiation
    dt_step = 1.0  # 1 day step for daily velocity
    pos_prev = _calc_orbital_position_raw(elem, jd_tt - dt_step)
    pos_next = _calc_orbital_position_raw(elem, jd_tt + dt_step)

    # Elements referred to a fixed equinox are precessed to J2000 (the
    # frame this function documents); equinox-of-date elements (JDATE)
    # stay in their of-date frame. The same rotation applies to the
    # speed samples so the rates stay frame-consistent.
    tequ = elem.equinox_jd
    if tequ is not None and abs(tequ - 2451545.0) > 1e-6:
        from .astrometry import _precess_ecliptic

        longitude, latitude = _precess_ecliptic(longitude, latitude, tequ, 2451545.0)
        p_lon, p_lat = _precess_ecliptic(pos_prev[0], pos_prev[1], tequ, 2451545.0)
        n_lon, n_lat = _precess_ecliptic(pos_next[0], pos_next[1], tequ, 2451545.0)
        pos_prev = (p_lon, p_lat, pos_prev[2])
        pos_next = (n_lon, n_lat, pos_next[2])

    dlon = (pos_next[0] - pos_prev[0]) / (2.0 * dt_step)
    # Handle wrap-around
    if dlon > 180.0 / (2.0 * dt_step):
        dlon -= 360.0 / (2.0 * dt_step)
    elif dlon < -180.0 / (2.0 * dt_step):
        dlon += 360.0 / (2.0 * dt_step)

    dlat = (pos_next[1] - pos_prev[1]) / (2.0 * dt_step)
    ddist = (pos_next[2] - pos_prev[2]) / (2.0 * dt_step)

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_orbital_position_raw(
    elem: OrbitalElements, jd_tt: float
) -> Tuple[float, float, float]:
    """
    Calculate raw position without velocity (helper for differentiation).
    """
    # Time in Julian centuries from epoch
    T = (jd_tt - elem.epoch_jd) / 36525.0
    # Days from epoch (for Keplerian motion)
    dt_days = jd_tt - elem.epoch_jd

    # Evaluate time-dependent elements
    # Same logic as calc_orbital_position for mean anomaly
    M_poly = elem.mean_anomaly.evaluate(T)
    if abs(elem.mean_anomaly.linear) < 10.0:
        n_deg_per_day = elem.get_mean_motion()
        M_deg = M_poly + n_deg_per_day * dt_days
    else:
        M_deg = M_poly

    e = elem.eccentricity.evaluate(T)
    omega_deg = elem.arg_perihelion.evaluate(T)
    Omega_deg = elem.asc_node.evaluate(T)
    i_deg = elem.inclination.evaluate(T)
    a = elem.semi_axis

    # Normalize mean anomaly
    M_deg = M_deg % 360.0
    M_rad = math.radians(M_deg)

    # For circular orbits (e=0), mean anomaly = true anomaly
    if e < 1e-10:
        nu = M_rad
        r = a
    else:
        # Solve Kepler's equation
        E = _solve_kepler_equation(M_rad, e)

        # True anomaly
        sqrt_term = math.sqrt((1.0 + e) / (1.0 - e))
        nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

        # Distance
        r = a * (1.0 - e * math.cos(E))

    # Argument of latitude
    omega_rad = math.radians(omega_deg)
    u = nu + omega_rad

    # Convert to ecliptic coordinates
    i_rad = math.radians(i_deg)
    Omega_rad = math.radians(Omega_deg)

    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0

    return (longitude, latitude, r)


# =============================================================================
# CALCULATION FUNCTIONS
# =============================================================================


def is_hypothetical_body(ipl: int) -> bool:
    """
    Check if a body ID corresponds to a hypothetical body.

    Args:
        ipl: Planet/body ID

    Returns:
        True if the body is a hypothetical body, False otherwise.

    Example:
        >>> from libephemeris.hypothetical import is_hypothetical_body, CUPIDO
        >>> is_hypothetical_body(CUPIDO)
        True
        >>> is_hypothetical_body(0)  # SUN
        False
    """
    return FICT_OFFSET <= ipl < FICT_OFFSET + 30


def get_hypothetical_name(ipl: int) -> str:
    """
    Get the name of a hypothetical body.

    Args:
        ipl: Planet/body ID

    Returns:
        Name of the hypothetical body, or "Unknown" if not found.

    Example:
        >>> from libephemeris.hypothetical import get_hypothetical_name, CUPIDO
        >>> get_hypothetical_name(CUPIDO)
        'Cupido'
    """
    return HYPOTHETICAL_NAMES.get(ipl, f"Unknown ({ipl})")


# Orbital elements of the classical hypothetical/predicted planets the
# reference exposes as bodies 49-54. The values are historical data from
# the published predictions (collected, like every fictitious-body
# element set, in the standard orbital-elements file format):
# - Nibiru: elements as circulated in the astrological community for
#   Z. Sitchin's popular "12th planet" lore.
# - Harrington: R.S. Harrington's Planet X search orbit (AJ 96, 1476,
#   1988).
# - Leverrier / Adams: the 1846 theoretical orbits of Neptune by
#   U. Le Verrier (CRAS 23) and J.C. Adams.
# - Lowell: P. Lowell, "Memoir on a Trans-Neptunian Planet" (1915).
# - Pickering: W.H. Pickering's Planet P prediction (1928).
FICTITIOUS_ORBITAL_ELEMENTS = {
    NIBIRU: OrbitalElements(
        name="Nibiru",
        epoch_jd=1856113.380954,
        equinox_jd=1856113.380954,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(0.0),
        semi_axis=234.8921,
        eccentricity=TPolynomial(0.981092),
        arg_perihelion=TPolynomial(103.966),
        asc_node=TPolynomial(-44.567),
        inclination=TPolynomial(158.708),
    ),
    HARRINGTON: OrbitalElements(
        name="Harrington",
        epoch_jd=2374696.5,
        equinox_jd=2451545.0,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(0.0),
        semi_axis=101.2,
        eccentricity=TPolynomial(0.411),
        arg_perihelion=TPolynomial(208.5),
        asc_node=TPolynomial(275.4),
        inclination=TPolynomial(32.4),
    ),
    NEPTUNE_LEVERRIER: OrbitalElements(
        name="Leverrier",
        epoch_jd=2395662.5,
        equinox_jd=2395662.5,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(34.05),
        semi_axis=36.15,
        eccentricity=TPolynomial(0.10761),
        arg_perihelion=TPolynomial(284.75),
        asc_node=TPolynomial(0.0),
        inclination=TPolynomial(0.0),
    ),
    NEPTUNE_ADAMS: OrbitalElements(
        name="Adams",
        epoch_jd=2395662.5,
        equinox_jd=2395662.5,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(24.28),
        semi_axis=37.25,
        eccentricity=TPolynomial(0.12062),
        arg_perihelion=TPolynomial(299.11),
        asc_node=TPolynomial(0.0),
        inclination=TPolynomial(0.0),
    ),
    PLUTO_LOWELL: OrbitalElements(
        name="Lowell",
        epoch_jd=2425977.5,
        equinox_jd=2425977.5,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(281.0),
        semi_axis=43.0,
        eccentricity=TPolynomial(0.202),
        arg_perihelion=TPolynomial(204.9),
        asc_node=TPolynomial(0.0),
        inclination=TPolynomial(0.0),
    ),
    PLUTO_PICKERING: OrbitalElements(
        name="Pickering",
        epoch_jd=2425977.5,
        equinox_jd=2425977.5,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(48.95),
        semi_axis=55.1,
        eccentricity=TPolynomial(0.31),
        arg_perihelion=TPolynomial(280.1),
        asc_node=TPolynomial(100.0),
        inclination=TPolynomial(15.0),
    ),
}


def calc_fictitious_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """Heliocentric position of a classical predicted planet (49-54).

    Full Keplerian propagation of the published prediction orbits with
    equinox precession; see FICTITIOUS_ORBITAL_ELEMENTS for the sources.

    Returns (lon, lat, dist, dlon, dlat, ddist), J2000 ecliptic.
    """
    elem = FICTITIOUS_ORBITAL_ELEMENTS.get(ipl)
    if elem is None:
        from .exceptions import UnknownBodyError

        raise UnknownBodyError(message=f"no fictitious elements for body {ipl}")

    # Propagate through _keplerian_to_ecliptic_j2000: it applies the
    # element-equinox precession to J2000 correctly and uses the
    # Gaussian mean motion (calc_orbital_position does neither).
    from types import SimpleNamespace

    ad = SimpleNamespace(
        epoch=elem.epoch_jd,
        M0=elem.mean_anomaly.evaluate(0.0),
        n=0.9856076686 / elem.semi_axis**1.5,
        e=elem.eccentricity.evaluate(0.0),
        a=elem.semi_axis,
        omega=elem.arg_perihelion.evaluate(0.0),
        Omega=elem.asc_node.evaluate(0.0),
        i=elem.inclination.evaluate(0.0),
        equinox_jd=elem.equinox_jd,
    )
    ad_typed = cast(UranianKeplerianElements, ad)
    lon, lat, dist = _keplerian_to_ecliptic_j2000(ad_typed, jd_tt)
    h = 0.5
    lon_p, lat_p, dist_p = _keplerian_to_ecliptic_j2000(ad_typed, jd_tt - h)
    lon_n, lat_n, dist_n = _keplerian_to_ecliptic_j2000(ad_typed, jd_tt + h)
    dlon = ((lon_n - lon_p + 180.0) % 360.0) - 180.0
    return (lon, lat, dist, dlon, lat_n - lat_p, dist_n - dist_p)


def calc_uranian_longitude(ipl: int, jd_tt: float) -> float:
    """
    Calculate the ecliptic longitude of a Uranian (Hamburg School) planet.

    Uses the secular polynomial formula:
        longitude = L0 + n * T + amplitude * sin(radians(phase + phase_rate * T))

    where T is Julian centuries from J2000.0

    Args:
        ipl: Uranian planet ID (CUPIDO through POSEIDON)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Ecliptic longitude in degrees (0-360)

    Raises:
        ValueError: If ipl is not a valid Uranian planet ID.

    Example:
        >>> from libephemeris.hypothetical import calc_uranian_longitude, CUPIDO
        >>> lon = calc_uranian_longitude(CUPIDO, 2451545.0)  # J2000.0
        >>> print(f"Cupido longitude: {lon:.4f}")
    """
    if ipl not in URANIAN_KEPLERIAN_ELEMENTS:
        raise ValueError(f"Body ID {ipl} is not a valid Uranian planet")

    # Delegate to the verified Keplerian element set; the legacy
    # mean-longitude oscillation table (URANIAN_ELEMENTS) carries
    # contradictory values and is no longer consulted.
    return calc_uranian_planet(ipl, jd_tt)[0]


def calc_uranian_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the full position of a Uranian (Hamburg School) planet.

    Returns ecliptic longitude, latitude, distance, and their rates of change.
    Uranian planets are assumed to be on the ecliptic (latitude = 0) and at
    a fixed distance (based on their orbital period).

    Args:
        ipl: Uranian planet ID (CUPIDO through POSEIDON)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees
            - latitude: Ecliptic latitude in degrees (always 0)
            - distance: Distance in AU (estimated from mean motion)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0)

    Raises:
        ValueError: If ipl is not a valid Uranian planet ID.

    Example:
        >>> from libephemeris.hypothetical import calc_uranian_position, KRONOS
        >>> pos = calc_uranian_position(KRONOS, 2451545.0)
        >>> print(f"Kronos at {pos[0]:.4f} deg, dist {pos[2]:.2f} AU")
    """
    if ipl not in URANIAN_KEPLERIAN_ELEMENTS:
        raise ValueError(f"Body ID {ipl} is not a valid Uranian planet")

    # Delegate to the verified Keplerian element set (see
    # calc_uranian_longitude).
    return calc_uranian_planet(ipl, jd_tt)


def calc_cupido(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Cupido using Keplerian propagation.

    Cupido is the first Hamburg School Uranian planet. This function uses
    orbital elements derived from Witte/Sieggrun 1928 ephemeris tables.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Cupido)
            - distance: Distance from Sun in AU (40.99837 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_cupido
        >>> pos = calc_cupido(2451545.0)  # J2000.0
        >>> print(f"Cupido at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(CUPIDO, jd_tt)


def calc_hades(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Hades using Keplerian propagation.

    Hades is the second Hamburg School Uranian planet. This function uses
    orbital elements derived from Witte/Sieggrun 1928 ephemeris tables.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_hades
        >>> pos = calc_hades(2451545.0)  # J2000.0
        >>> print(f"Hades at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(HADES, jd_tt)


def calc_zeus(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Zeus using Keplerian propagation.

    Zeus is the third Hamburg School Uranian planet. This function uses
    orbital elements derived from Witte/Sieggrun 1928 ephemeris tables.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Zeus)
            - distance: Distance from Sun in AU (59.21436 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_zeus
        >>> pos = calc_zeus(2451545.0)  # J2000.0
        >>> print(f"Zeus at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(ZEUS, jd_tt)


def calc_kronos(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Kronos using Keplerian propagation.

    Kronos is the fourth Hamburg School Uranian planet. This function uses
    orbital elements derived from Witte/Sieggrun 1928 ephemeris tables.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Kronos)
            - distance: Distance from Sun in AU (64.81690 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_kronos
        >>> pos = calc_kronos(2451545.0)  # J2000.0
        >>> print(f"Kronos at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(KRONOS, jd_tt)


def calc_apollon(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Apollon using Keplerian propagation.

    Apollon is the fifth Hamburg School Uranian planet. This function uses
    orbital elements derived from Witte/Sieggrun 1928 ephemeris tables.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Apollon)
            - distance: Distance from Sun in AU (70.361180 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_apollon
        >>> pos = calc_apollon(2451545.0)  # J2000.0
        >>> print(f"Apollon at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(APOLLON, jd_tt)


def calc_admetos(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Admetos using Keplerian propagation.

    Admetos is the sixth Hamburg School Uranian planet. This function uses
    orbital elements derived from Witte/Sieggrun 1928 ephemeris tables.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Admetos)
            - distance: Distance from Sun in AU (73.736396 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_admetos
        >>> pos = calc_admetos(2451545.0)  # J2000.0
        >>> print(f"Admetos at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(ADMETOS, jd_tt)


def calc_vulkanus(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Vulkanus using Keplerian propagation.

    Vulkanus is the seventh Hamburg School Uranian planet. This function uses
    orbital elements derived from Witte/Sieggrun 1928 ephemeris tables.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Vulkanus)
            - distance: Distance from Sun in AU (77.445895 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_vulkanus
        >>> pos = calc_vulkanus(2451545.0)  # J2000.0
        >>> print(f"Vulkanus at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(VULKANUS, jd_tt)


def calc_poseidon(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Poseidon using Keplerian propagation.

    Poseidon is the eighth Hamburg School Uranian planet. This function uses
    orbital elements derived from Witte/Sieggrun 1928 ephemeris tables.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (always 0 for Poseidon)
            - distance: Distance from Sun in AU (83.666307 AU, constant)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (always 0)
            - ddist: Daily distance change in AU/day (always 0 for e=0)

    Example:
        >>> from libephemeris.hypothetical import calc_poseidon
        >>> pos = calc_poseidon(2451545.0)  # J2000.0
        >>> print(f"Poseidon at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    return calc_uranian_planet(POSEIDON, jd_tt)


def calc_transpluto(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Transpluto (Isis) using Keplerian propagation.

    Transpluto is a hypothetical trans-Plutonian planet proposed by astrologer Ram,
    with orbital elements from Strubell, "Die Sterne" 3/1952. Unlike the Hamburg
    School Uranian planets which have nearly circular orbits, Transpluto has
    significant eccentricity (e=0.3).

    Orbital elements:
        - Epoch: 2368547.66 (1772.76)
        - Semi-major axis: 77.775 AU
        - Eccentricity: 0.3
        - Inclination: 0 degrees
        - Argument of perihelion: 0.7 degrees
        - Ascending node: 0 degrees
        - Mean anomaly at epoch: 0 degrees (at perihelion)

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (0 for Transpluto)
            - distance: Distance from Sun in AU
            - dlon: Daily heliocentric longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (0)
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_transpluto
        >>> pos = calc_transpluto(2451545.0)  # J2000.0
        >>> print(f"Transpluto at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    elements = TRANSPLUTO_KEPLERIAN_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    u = nu + math.radians(elements.omega)

    # Convert to ecliptic coordinates
    # For zero inclination and ascending node, this simplifies considerably
    i_rad = math.radians(elements.i)
    Omega_rad = math.radians(elements.Omega)

    # Position in orbital plane
    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    # Rotate to ecliptic frame
    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    # Convert to spherical coordinates
    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    # Calculate velocity via central difference numerical differentiation
    dt_step = 1.0  # 1 day step for daily velocity
    pos_prev = _calc_transpluto_raw(jd_tt - dt_step)
    pos_next = _calc_transpluto_raw(jd_tt + dt_step)

    dlon = (pos_next[0] - pos_prev[0]) / (2.0 * dt_step)
    # Handle wrap-around
    if dlon > 180.0 / (2.0 * dt_step):
        dlon -= 360.0 / (2.0 * dt_step)
    elif dlon < -180.0 / (2.0 * dt_step):
        dlon += 360.0 / (2.0 * dt_step)

    dlat = (pos_next[1] - pos_prev[1]) / (2.0 * dt_step)
    ddist = (pos_next[2] - pos_prev[2]) / (2.0 * dt_step)

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_transpluto_raw(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Transpluto position without velocity (helper for differentiation).
    """
    elements = TRANSPLUTO_KEPLERIAN_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    u = nu + math.radians(elements.omega)

    # Convert to ecliptic coordinates
    i_rad = math.radians(elements.i)
    Omega_rad = math.radians(elements.Omega)

    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    return (longitude, latitude, distance)


# J2000.0 epoch Julian Day for precession calculations
_J2000_JD: float = 2451545.0


def _keplerian_to_ecliptic_j2000(
    elements: UranianKeplerianElements, jd_tt: float
) -> Tuple[float, float, float]:
    """
    Full Keplerian propagation with equinox precession to J2000 ecliptic frame.

    Implements the standard celestial mechanics algorithm for propagating
    orbital elements:

    1. Propagate mean anomaly from epoch using Gaussian gravitational constant
    2. Solve Kepler's equation (Newton-Raphson) for eccentric anomaly E
    3. Compute position in orbital plane: x = a*(cosE - e), y = a*sqrt(1-e^2)*sinE
    4. Apply Gaussian vectors (PQR matrix) to transform to ecliptic Cartesian
       coordinates in the equinox reference frame
    5. If equinox != J2000: rotate ecliptic -> equatorial at equinox obliquity,
       precess equatorial from equinox epoch to J2000, then rotate back to
       ecliptic J2000

    This algorithm matches standard Keplerian orbit propagation as described
    in Meeus "Astronomical Algorithms" Ch. 33 and other celestial mechanics
    references.

    Args:
        elements: Orbital elements including equinox information
        jd_tt: Julian Day in Terrestrial Time

    Returns:
        Tuple of (longitude_deg, latitude_deg, distance_AU) in J2000 ecliptic
    """
    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly at date (radians)
    M_deg = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M_deg)

    # Solve Kepler's equation: M = E - e*sin(E)
    e = elements.e
    if e < 1e-10:
        # Circular orbit: E = M, no need to solve
        E = M_rad
    else:
        E = _solve_kepler_equation(M_rad, e)

    # Position in orbital plane (Cartesian, origin at focus)
    cos_E = math.cos(E)
    sin_E = math.sin(E)
    fac = math.sqrt((1.0 - e) * (1.0 + e))  # sqrt(1 - e^2), numerically stable

    x_orb = elements.a * (cos_E - e)
    y_orb = elements.a * fac * sin_E

    # Gaussian vectors (PQR matrix) - standard orbital mechanics
    # Transforms from orbital plane to ecliptic Cartesian in the equinox frame
    # Reference: Meeus "Astronomical Algorithms", Ch. 33
    omega_rad = math.radians(elements.omega)
    Omega_rad = math.radians(elements.Omega)
    incl_rad = math.radians(elements.i)

    cos_omega = math.cos(omega_rad)
    sin_omega = math.sin(omega_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)
    cos_incl = math.cos(incl_rad)
    sin_incl = math.sin(incl_rad)

    # PQR matrix elements (Gaussian vectors)
    # P = direction of perihelion, Q = 90 degrees ahead in orbit
    Px = cos_omega * cos_Omega - sin_omega * cos_incl * sin_Omega
    Qx = -sin_omega * cos_Omega - cos_omega * cos_incl * sin_Omega
    Py = cos_omega * sin_Omega + sin_omega * cos_incl * cos_Omega
    Qy = -sin_omega * sin_Omega + cos_omega * cos_incl * cos_Omega
    Pz = sin_omega * sin_incl
    Qz = cos_omega * sin_incl

    # Ecliptic Cartesian in equinox frame
    x_ecl = Px * x_orb + Qx * y_orb
    y_ecl = Py * x_orb + Qy * y_orb
    z_ecl = Pz * x_orb + Qz * y_orb

    # If equinox != J2000, transform to J2000 frame via precession
    tequ = elements.equinox_jd
    if abs(tequ - _J2000_JD) > 1e-6:
        # Step 1: Ecliptic -> equatorial at equinox epoch obliquity
        # Rotation around x-axis by -eps (ecliptic to equator)
        eps = float(erfa.obl06(_J2000_JD, tequ - _J2000_JD))
        cos_eps = math.cos(eps)
        sin_eps = math.sin(eps)

        x_eq = x_ecl
        y_eq = y_ecl * cos_eps - z_ecl * sin_eps
        z_eq = y_ecl * sin_eps + z_ecl * cos_eps

        # Step 2: Precess equatorial from equinox epoch to J2000
        # erfa.pmat06 gives the rotation matrix FROM J2000 TO date.
        # To go from date TO J2000, apply the transpose (inverse of rotation).
        P = erfa.pmat06(_J2000_JD, tequ - _J2000_JD)
        # P^T maps date -> J2000
        x_j2000 = P[0][0] * x_eq + P[1][0] * y_eq + P[2][0] * z_eq
        y_j2000 = P[0][1] * x_eq + P[1][1] * y_eq + P[2][1] * z_eq
        z_j2000 = P[0][2] * x_eq + P[1][2] * y_eq + P[2][2] * z_eq

        # Step 3: Equatorial J2000 -> ecliptic J2000
        # Rotation around x-axis by +eps_j2000 (equator to ecliptic)
        eps_j2000 = float(erfa.obl06(_J2000_JD, 0.0))
        cos_eps0 = math.cos(eps_j2000)
        sin_eps0 = math.sin(eps_j2000)

        x_ecl_final = x_j2000
        y_ecl_final = y_j2000 * cos_eps0 + z_j2000 * sin_eps0
        z_ecl_final = -y_j2000 * sin_eps0 + z_j2000 * cos_eps0
    else:
        x_ecl_final = x_ecl
        y_ecl_final = y_ecl
        z_ecl_final = z_ecl

    # Cartesian -> spherical coordinates
    r = math.sqrt(x_ecl_final**2 + y_ecl_final**2 + z_ecl_final**2)
    longitude = math.degrees(math.atan2(y_ecl_final, x_ecl_final)) % 360.0
    if r > 0:
        sin_lat = max(-1.0, min(1.0, z_ecl_final / r))
        latitude = math.degrees(math.asin(sin_lat))
    else:
        latitude = 0.0

    return (longitude, latitude, r)


def calc_uranian_planet(
    body_id: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of any Uranian planet using full Keplerian propagation.

    Handles all eight Hamburg School Uranian planets (Cupido, Hades, Zeus,
    Kronos, Apollon, Admetos, Vulkanus, Poseidon) using proper Keplerian
    orbit propagation with:

    - Mean anomaly propagation using Gaussian gravitational constant
    - Kepler's equation solving (Newton-Raphson iteration)
    - Gaussian vector (PQR) transformation to ecliptic coordinates
    - Equinox precession from J1900 coordinate frame to J2000

    Velocity is computed via central-difference numerical differentiation.

    Args:
        body_id: Uranian planet ID (CUPIDO through POSEIDON, i.e., 40-47)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees (0-360), J2000 frame
            - latitude: Ecliptic latitude in degrees, J2000 frame
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Raises:
        ValueError: If body_id is not a valid Uranian planet ID.

    Example:
        >>> from libephemeris.hypothetical import calc_uranian_planet, CUPIDO
        >>> pos = calc_uranian_planet(CUPIDO, 2451545.0)  # J2000.0
        >>> print(f"Cupido at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    if body_id not in URANIAN_KEPLERIAN_ELEMENTS:
        raise ValueError(
            f"Body ID {body_id} is not a valid Uranian planet. "
            f"Valid IDs: {list(URANIAN_KEPLERIAN_ELEMENTS.keys())}"
        )

    elements = URANIAN_KEPLERIAN_ELEMENTS[body_id]

    # Position at requested time
    pos = _keplerian_to_ecliptic_j2000(elements, jd_tt)
    longitude, latitude, distance = pos

    # Velocity via central-difference numerical differentiation
    dt_step = 1.0  # 1 day step for daily velocity
    pos_prev = _keplerian_to_ecliptic_j2000(elements, jd_tt - dt_step)
    pos_next = _keplerian_to_ecliptic_j2000(elements, jd_tt + dt_step)

    dlon = (pos_next[0] - pos_prev[0]) / (2.0 * dt_step)
    # Handle wrap-around at the 0/360 boundary. The positions are normalised to
    # [0, 360), so a boundary crossing shows up as a ~360 deg jump in the raw
    # difference; after dividing by 2*dt_step the artefact is ~360/(2*dt_step),
    # so the detection threshold and correction must carry the same factor.
    # Valid only because these bodies are slow (<<90 deg/day at dt_step=1): a
    # real speed above 180/(2*dt_step) cannot occur, so it can only be a wrap.
    # Do NOT reuse this for fast bodies -- it would misread real motion.
    if dlon > 180.0 / (2.0 * dt_step):
        dlon -= 360.0 / (2.0 * dt_step)
    elif dlon < -180.0 / (2.0 * dt_step):
        dlon += 360.0 / (2.0 * dt_step)

    dlat = (pos_next[1] - pos_prev[1]) / (2.0 * dt_step)
    ddist = (pos_next[2] - pos_prev[2]) / (2.0 * dt_step)

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_uranian_planet_raw(body_id: int, jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Uranian planet position without velocity (helper for differentiation).

    Uses full Keplerian propagation with Gaussian vectors and equinox precession
    to J2000.

    Args:
        body_id: Uranian planet ID
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance)
    """
    elements = URANIAN_KEPLERIAN_ELEMENTS[body_id]
    return _keplerian_to_ecliptic_j2000(elements, jd_tt)


def calc_vulcan(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Vulcan using the standard orbital elements.

    Vulcan is a hypothetical intramercurial planet. Various astrologers proposed
    different orbital elements for Vulcan; this implementation uses the standard
    orbital elements (body #16).

    Unlike other hypothetical bodies, Vulcan has time-dependent orbital elements:
        - Mean anomaly: 252.8987988 + 707550.7341 * T degrees
        - Argument of perihelion: 322.212069 + 1670.056 * T degrees
        - Ascending node: 47.787931 - 1670.056 * T degrees

    where T = Julian centuries from J1900.0 (JD 2415020.0)

    Orbital elements:
        - Epoch: J1900.0 (JD 2415020.0)
        - Equinox: JDATE (equinox of date)
        - Semi-major axis: 0.13744 AU (inside Mercury's orbit)
        - Eccentricity: 0.019 (nearly circular)
        - Inclination: 7.5 degrees
        - Orbital period: ~18.6 days (derived from mean motion)

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily heliocentric longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_vulcan
        >>> pos = calc_vulcan(2451545.0)  # J2000.0
        >>> print(f"Vulcan at {pos[0]:.4f} deg, distance {pos[2]:.4f} AU")
    """
    pos = _calc_vulcan_raw(jd_tt)
    longitude, latitude, distance = pos

    # Calculate velocity via central difference numerical differentiation
    dt_step = 1.0  # 1 day step for daily velocity
    pos_prev = _calc_vulcan_raw(jd_tt - dt_step)
    pos_next = _calc_vulcan_raw(jd_tt + dt_step)

    dlon = (pos_next[0] - pos_prev[0]) / (2.0 * dt_step)
    # Handle wrap-around (Vulcan moves ~19 deg/day, so this is important)
    if dlon > 180.0 / (2.0 * dt_step):
        dlon -= 360.0 / (2.0 * dt_step)
    elif dlon < -180.0 / (2.0 * dt_step):
        dlon += 360.0 / (2.0 * dt_step)

    dlat = (pos_next[1] - pos_prev[1]) / (2.0 * dt_step)
    ddist = (pos_next[2] - pos_prev[2]) / (2.0 * dt_step)

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_vulcan_raw(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Vulcan position without velocity (helper for differentiation).

    Implements the time-dependent orbital elements from the orbital elements data.
    """
    elements = VULCAN_ELEMENTS

    # Time in Julian centuries from epoch (J1900.0)
    T = (jd_tt - elements.epoch) / 36525.0

    # Compute time-dependent elements
    # Mean anomaly: M0 + n_century * T
    M = (elements.M0 + elements.n_century * T) % 360.0
    M_rad = math.radians(M)

    # Argument of perihelion: omega0 + omega_rate * T
    omega = elements.omega0 + elements.omega_rate * T
    omega_rad = math.radians(omega)

    # Ascending node: Omega0 + Omega_rate * T
    Omega = elements.Omega0 + elements.Omega_rate * T
    Omega_rad = math.radians(Omega)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    u = nu + omega_rad

    # Convert to ecliptic coordinates
    i_rad = math.radians(elements.i)

    # Position in orbital plane
    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    # Rotate to ecliptic frame
    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    # Convert to spherical coordinates
    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    return (longitude, latitude, distance)


def calc_waldemath(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the geocentric position of the Waldemath hypothetical second moon.

    Dr. Georg Waldemath's hypothetical second moon of Earth (1898). This is a different body from:
    - Mean Lilith (MEAN_APOG) - the mean lunar apogee
    - True Lilith (OSCU_APOG) - the osculating lunar apogee
    - The second focus of the lunar orbit

    Waldemath claimed to have observed this body in 1898, describing it as a dark
    moon with an orbital period of approximately 119 days.

    Orbital elements:
        - Epoch: J2000.0 (JD 2451545.0)
        - Semi-major axis: 0.0029833 AU (~446,200 km, ~1.16x Moon distance)
        - Eccentricity: 0.0 (circular orbit)
        - Inclination: 0.0 degrees (on ecliptic)
        - Mean longitude at epoch: 248.8833 degrees
        - Orbital period: ~119 days

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Geocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (0 for circular ecliptic orbit)
            - distance: Distance from Earth in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (0)
            - ddist: Daily distance change in AU/day (0 for circular orbit)

    Example:
        >>> from libephemeris.hypothetical import calc_waldemath
        >>> pos = calc_waldemath(2451545.0)  # J2000.0
        >>> print(f"Waldemath at {pos[0]:.4f} deg, distance {pos[2]:.6f} AU")
    """
    elements = WALDEMATH_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # For a circular orbit (e = 0), mean longitude = true longitude
    # Simply propagate the mean longitude
    longitude = (elements.L0 + elements.n * dt) % 360.0

    # Waldemath is assumed to be on the ecliptic (zero inclination)
    latitude = 0.0

    # Distance is constant for circular orbit (equal to semi-major axis)
    distance = elements.a

    # Daily motion is simply the mean motion for circular orbit
    dlon = elements.n

    # No latitude or distance change for circular orbit on ecliptic
    dlat = 0.0
    ddist = 0.0

    return (longitude, latitude, distance, dlon, dlat, ddist)


def calc_transpluto_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of Transpluto (Isis/Persephone).

    Uses Keplerian orbital mechanics with the elements stored in
    HYPOTHETICAL_ELEMENTS[ISIS].

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_transpluto_position
        >>> pos = calc_transpluto_position(2451545.0)  # J2000.0
        >>> print(f"Transpluto at {pos[0]:.4f} deg")
    """
    return _calc_keplerian_position(ISIS, jd_tt)


def _solve_kepler_equation(M: float, e: float, tol: float = 1e-8) -> float:
    """
    Solve Kepler's equation M = E - e * sin(E) for eccentric anomaly E.

    Uses Newton-Raphson iteration.

    Args:
        M: Mean anomaly in radians
        e: Eccentricity (0 <= e < 1)
        tol: Convergence tolerance

    Returns:
        Eccentric anomaly E in radians
    """
    # Initial guess
    E = M if e < 0.8 else math.pi

    for _ in range(30):
        f = E - e * math.sin(E) - M
        fp = 1.0 - e * math.cos(E)
        E_new = E - f / fp

        if abs(E_new - E) < tol:
            return E_new
        E = E_new

    return E


def _calc_keplerian_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate position using Keplerian orbital elements.

    Internal function for hypothetical bodies with Keplerian orbits.

    Args:
        ipl: Body ID
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
    """
    if ipl not in HYPOTHETICAL_ELEMENTS:
        raise ValueError(f"Body ID {ipl} does not have Keplerian elements")

    elements = HYPOTHETICAL_ELEMENTS[ipl]

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = elements.M0 + elements.n * dt
    M_rad = math.radians(M % 360.0)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    u = nu + math.radians(elements.omega)

    # Convert to ecliptic coordinates
    # For zero inclination and ascending node, this simplifies considerably
    i_rad = math.radians(elements.i)
    Omega_rad = math.radians(elements.Omega)

    # Position in orbital plane
    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    # Rotate to ecliptic frame
    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    # Convert to spherical coordinates
    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    # Calculate velocity via central difference numerical differentiation
    dt_step = 1.0 / 86400.0  # 1 second
    pos_prev = _calc_keplerian_position_raw(ipl, jd_tt - dt_step)
    pos_next = _calc_keplerian_position_raw(ipl, jd_tt + dt_step)

    dlon = (pos_next[0] - pos_prev[0]) / (2.0 * dt_step)
    # Handle wrap-around
    if dlon > 180.0 / (2.0 * dt_step):
        dlon -= 360.0 / (2.0 * dt_step)
    elif dlon < -180.0 / (2.0 * dt_step):
        dlon += 360.0 / (2.0 * dt_step)

    dlat = (pos_next[1] - pos_prev[1]) / (2.0 * dt_step)
    ddist = (pos_next[2] - pos_prev[2]) / (2.0 * dt_step)

    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_keplerian_position_raw(ipl: int, jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Keplerian position without velocity (helper for differentiation).
    """
    if ipl not in HYPOTHETICAL_ELEMENTS:
        raise ValueError(f"Body ID {ipl} does not have Keplerian elements")

    elements = HYPOTHETICAL_ELEMENTS[ipl]

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = elements.M0 + elements.n * dt
    M_rad = math.radians(M % 360.0)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude
    u = nu + math.radians(elements.omega)

    # Convert to ecliptic coordinates
    i_rad = math.radians(elements.i)
    Omega_rad = math.radians(elements.Omega)

    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    return (longitude, latitude, distance)


def calc_white_moon_position(
    jd_tt: float,
    use_true_lilith: bool = False,
) -> Tuple[float, float, float, float, float, float]:
    """Position of the White Moon (Selena).

    Uses the published circular geocentric orbit (the standard
    orbital-elements convention for this body): mean longitude
    242.2205555 + 5143.5418158*T degrees (T in Julian centuries from
    J2000, equinox of date), radius 0.05280098949 AU, zero eccentricity
    and inclination. The longitude advances 0.140822 degrees/day
    (a 7-year cycle).

    The ``use_true_lilith`` parameter is kept for backward compatibility
    and is ignored: the reference defines Selena by these elements, not
    as the anti-apogee.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)
        use_true_lilith: Ignored (kept for API compatibility).

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist) -
        geocentric ecliptic of date, distance in AU.
    """
    t_cent = (jd_tt - 2451545.0) / 36525.0
    lon = (242.2205555 + 5143.5418158 * t_cent) % 360.0
    dlon = 5143.5418158 / 36525.0
    return (lon, 0.0, 0.05280098949, dlon, 0.0, 0.0)


def calc_waldemath_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of the Waldemath Moon (hypothetical second moon of Earth).

    Dr. Georg Waldemath's hypothetical second moon of Earth (1898). This is a different body from:
    - Mean Lilith (MEAN_APOG) - the mean lunar apogee
    - True Lilith (OSCU_APOG) - the osculating lunar apogee
    - The second focus of the lunar orbit

    Waldemath claimed to have observed this body in 1898, describing it as a dark
    moon with an orbital period of approximately 119 days.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Geocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Earth in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Note:
        This is a hypothetical body. Waldemath's observations were never confirmed
        and no such second moon has been found to exist.

    Example:
        >>> from libephemeris.hypothetical import calc_waldemath_position
        >>> pos = calc_waldemath_position(2451545.0)  # J2000.0
        >>> print(f"Waldemath at {pos[0]:.4f} deg")
    """
    return calc_waldemath(jd_tt)


def calc_proserpina(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Proserpina using Keplerian propagation.

    Proserpina is a hypothetical trans-Plutonian planet used by some astrologers.
    This is NOT the same as the asteroid 26 Proserpina. Unlike other hypothetical
    bodies documented in the orbital elements data, Proserpina is not part of
    the standard fictitious bodies set.

    The name "Proserpina" refers to the Roman goddess of the underworld (Greek:
    Persephone), wife of Pluto. In astrological usage, it represents transformation,
    cycles of death and rebirth, and the shadow self.

    Orbital elements used (traditional astrological):
        - Epoch: J2000.0 (JD 2451545.0)
        - Semi-major axis: 81.0 AU (beyond Neptune and Pluto)
        - Eccentricity: 0.0 (circular orbit)
        - Inclination: 0.0 degrees (on ecliptic)
        - Orbital period: ~729 years (derived from Kepler's 3rd law)

    Note: Different astrologers may use different orbital elements for Proserpina.
    This implementation uses a simple circular orbit model.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees (0 for circular orbit on ecliptic)
            - distance: Distance from Sun in AU (81.0 AU, constant for circular orbit)
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day (0)
            - ddist: Daily distance change in AU/day (0 for circular orbit)

    Example:
        >>> from libephemeris.hypothetical import calc_proserpina
        >>> pos = calc_proserpina(2451545.0)  # J2000.0
        >>> print(f"Proserpina at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    elements = HYPOTHETICAL_ELEMENTS[PROSERPINA]

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # For a circular orbit (e = 0), mean longitude = true longitude
    # Simply propagate the mean longitude
    longitude = (elements.M0 + elements.n * dt) % 360.0

    # Proserpina is assumed to be on the ecliptic (zero inclination)
    latitude = 0.0

    # Distance is constant for circular orbit (equal to semi-major axis)
    distance = elements.a

    # Daily motion is simply the mean motion for circular orbit
    dlon = elements.n

    # No latitude or distance change for circular orbit on ecliptic
    dlat = 0.0
    ddist = 0.0

    return (longitude, latitude, distance, dlon, dlat, ddist)


def calc_planet_x_lowell(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Lowell's Planet X using Keplerian propagation.

    Percival Lowell's hypothetical "Planet X" was a trans-Neptunian planet he predicted
    in 1915 based on perceived perturbations in Uranus's orbit. The search for this
    planet at Lowell Observatory led to Clyde Tombaugh's discovery of Pluto in 1930.

    Historical note: Pluto was far too small (0.2% of Earth's mass) to cause the
    perturbations Lowell attributed to "Planet X". The perceived perturbations were
    later explained as observational errors. Pluto's discovery near Lowell's predicted
    position was essentially a fortunate coincidence.

    Orbital elements from Lowell's 1915 prediction "Memoir on a Trans-Neptunian Planet":
        - Semi-major axis: 43.0 AU
        - Eccentricity: 0.202
        - Inclination: 10 degrees to ecliptic
        - Orbital period: ~282 years

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_planet_x_lowell
        >>> pos = calc_planet_x_lowell(2451545.0)  # J2000.0
        >>> print(f"Planet X Lowell at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    # Delegate to the canonical published-element path (Hoyt 1980, via
    # data/fictitious_orbits.csv) that calc_hypothetical_position and the
    # public calc() API already use for this body, so every route agrees.
    # The legacy LOWELL_PLANET_X_ELEMENTS set and _calc_planet_x_lowell_raw
    # are no longer consulted at runtime (see CLEAN.md).
    return calc_fictitious_position(PLUTO_LOWELL, jd_tt)


def _calc_planet_x_lowell_raw(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Planet X Lowell position without velocity (helper for differentiation).
    """
    elements = LOWELL_PLANET_X_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    u = nu + math.radians(elements.omega)

    # Convert to ecliptic coordinates
    i_rad = math.radians(elements.i)
    Omega_rad = math.radians(elements.Omega)

    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    return (longitude, latitude, distance)


def calc_planet_x_pickering(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the heliocentric position of Pickering's Planet X using Keplerian propagation.

    William H. Pickering's hypothetical "Planet O" was a trans-Neptunian planet he predicted
    in 1919 based on supposed perturbations in Uranus and Neptune orbits. Pickering proposed
    several hypothetical planets (O, P, Q, R, S, T, U), with Planet O being the most famous.

    Historical note: Like Lowell's Planet X, Pickering's predictions were based on supposed
    perturbations that later proved to be observational errors. No planet was ever found
    at Pickering's predicted positions.

    Orbital elements from Pickering's 1919 prediction:
        - Semi-major axis: 51.9 AU
        - Eccentricity: 0.31
        - Inclination: 15 degrees to ecliptic
        - Orbital period: ~373.5 years (derived from Kepler's 3rd law)

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Heliocentric ecliptic longitude in degrees (0-360)
            - latitude: Ecliptic latitude in degrees
            - distance: Distance from Sun in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Example:
        >>> from libephemeris.hypothetical import calc_planet_x_pickering
        >>> pos = calc_planet_x_pickering(2451545.0)  # J2000.0
        >>> print(f"Planet X Pickering at {pos[0]:.4f} deg, distance {pos[2]:.2f} AU")
    """
    # Delegate to the canonical published-element path (Hoyt 1980, via
    # data/fictitious_orbits.csv) that calc_hypothetical_position and the
    # public calc() API already use for this body, so every route agrees.
    # The legacy PICKERING_PLANET_X_ELEMENTS set and
    # _calc_planet_x_pickering_raw are no longer consulted at runtime
    # (see CLEAN.md).
    return calc_fictitious_position(PLUTO_PICKERING, jd_tt)


def _calc_planet_x_pickering_raw(jd_tt: float) -> Tuple[float, float, float]:
    """
    Calculate raw Planet X Pickering position without velocity (helper for differentiation).
    """
    elements = PICKERING_PLANET_X_ELEMENTS

    # Time since epoch in days
    dt = jd_tt - elements.epoch

    # Mean anomaly
    M = (elements.M0 + elements.n * dt) % 360.0
    M_rad = math.radians(M)

    # Solve Kepler's equation for eccentric anomaly
    E = _solve_kepler_equation(M_rad, elements.e)

    # True anomaly
    sqrt_term = math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
    nu = 2.0 * math.atan(sqrt_term * math.tan(E / 2.0))

    # Distance from Sun (heliocentric)
    r = elements.a * (1.0 - elements.e * math.cos(E))

    # Argument of latitude (measured from ascending node)
    u = nu + math.radians(elements.omega)

    # Convert to ecliptic coordinates
    i_rad = math.radians(elements.i)
    Omega_rad = math.radians(elements.Omega)

    x_orb = r * math.cos(u)
    y_orb = r * math.sin(u)

    cos_i = math.cos(i_rad)
    sin_i = math.sin(i_rad)
    cos_Omega = math.cos(Omega_rad)
    sin_Omega = math.sin(Omega_rad)

    x_ecl = cos_Omega * x_orb - sin_Omega * cos_i * y_orb
    y_ecl = sin_Omega * x_orb + cos_Omega * cos_i * y_orb
    z_ecl = sin_i * y_orb

    longitude = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
    latitude = math.degrees(math.asin(z_ecl / r)) if r > 0 else 0.0
    distance = r

    return (longitude, latitude, distance)


def calc_hypothetical_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of any hypothetical body.

    This is the main entry point for hypothetical body calculations.
    Routes the request to the appropriate calculation function based on body ID.

    Args:
        ipl: Hypothetical body ID (CUPIDO through WALDEMATH)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist)
            - longitude: Ecliptic longitude in degrees
            - latitude: Ecliptic latitude in degrees
            - distance: Distance in AU
            - dlon: Daily longitude change in degrees/day
            - dlat: Daily latitude change in degrees/day
            - ddist: Daily distance change in AU/day

    Raises:
        ValueError: If ipl is not a valid hypothetical body ID.

    Example:
        >>> from libephemeris.hypothetical import calc_hypothetical_position, CUPIDO
        >>> pos = calc_hypothetical_position(CUPIDO, 2451545.0)
        >>> print(f"Cupido: {pos[0]:.4f} deg")
    """
    if not is_hypothetical_body(ipl):
        raise ValueError(f"Body ID {ipl} is not a hypothetical body")

    # Uranian planets (Hamburg School) - use Keplerian propagation
    # (Witte/Sieggruen elements refined by Neely 1988)
    if ipl in URANIAN_KEPLERIAN_ELEMENTS:
        return calc_uranian_planet(ipl, jd_tt)

    if ipl in FICTITIOUS_ORBITAL_ELEMENTS:
        return calc_fictitious_position(ipl, jd_tt)

    # Transpluto / Isis
    if ipl == ISIS:
        return calc_transpluto_position(jd_tt)

    # Vulcan (intramercurial hypothetical planet)
    if ipl == VULCAN:
        return calc_vulcan(jd_tt)

    # White Moon (Selena)
    if ipl == WHITE_MOON:
        return calc_white_moon_position(jd_tt)

    # Waldemath Black Moon
    if ipl == WALDEMATH:
        return calc_waldemath_position(jd_tt)

    # Proserpina (hypothetical trans-Plutonian planet)
    if ipl == PROSERPINA:
        return calc_proserpina(jd_tt)

    # Planet X Lowell (Lowell's predicted trans-Neptunian planet)
    if ipl == PLANET_X_LOWELL:
        return calc_planet_x_lowell(jd_tt)

    # Planet X Pickering (Pickering's predicted trans-Neptunian planet)
    if ipl == PLANET_X_PICKERING:
        return calc_planet_x_pickering(jd_tt)

    # Other Keplerian bodies
    if ipl in HYPOTHETICAL_ELEMENTS:
        return _calc_keplerian_position(ipl, jd_tt)

    # Unknown hypothetical body - return zeros
    return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def list_hypothetical_bodies() -> Dict[int, str]:
    """
    Get a dictionary of all supported hypothetical body IDs and names.

    Returns:
        Dictionary mapping body ID to body name.

    Example:
        >>> from libephemeris.hypothetical import list_hypothetical_bodies
        >>> bodies = list_hypothetical_bodies()
        >>> for id, name in bodies.items():
        ...     print(f"{id}: {name}")
    """
    return HYPOTHETICAL_NAMES.copy()


# ---------------------------------------------------------------------------
# Backward-compatible aliases (legacy SE-derived names)
# ---------------------------------------------------------------------------
# These aliases preserve backward compatibility for code using the old
# SE-derived naming convention. New code should use the canonical names:
#   OrbitalElements, parse_orbital_elements, get_orbital_body_by_name,
#   calc_orbital_position, _parse_orbital_elements_line
# ---------------------------------------------------------------------------

SeorbelElements = OrbitalElements
"""Backward-compatible alias for :class:`OrbitalElements`."""

parse_seorbel = parse_orbital_elements
"""Backward-compatible alias for :func:`parse_orbital_elements`."""

get_seorbel_body_by_name = get_orbital_body_by_name
"""Backward-compatible alias for :func:`get_orbital_body_by_name`."""

calc_seorbel_position = calc_orbital_position
"""Backward-compatible alias for :func:`calc_orbital_position`."""

_parse_seorbel_line = _parse_orbital_elements_line
"""Backward-compatible alias for :func:`_parse_orbital_elements_line`."""
