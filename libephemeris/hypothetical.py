# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Hypothetical and fictitious body calculations for LibEphemeris.

The historical bodies exposed as IDs 40--58 are mathematical conventions, not
claims that the objects exist.  Their source classification is recorded in
``HYPOTHETICAL_PROVENANCE`` and in
``docs/methodology/hypothetical-bodies.md``.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Dict, List, Optional, Union, cast, NoReturn

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

# Additional reference-compatible hypothetical-body identifiers
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
@dataclass
class UranianElements:
    """Mean-longitude element container for Uranian hypothetical bodies."""

    name: str
    L0: float
    n: float
    amplitude: float = 0.0
    phase: float = 0.0
    phase_rate: float = 0.0


URANIAN_ELEMENTS: Dict[int, UranianElements] = {
    CUPIDO: UranianElements("Cupido", 241.2067, 1.091437, 1.5, 21.6, 83.0),
    HADES: UranianElements("Hades", 176.0581, 0.736380, 1.5, 258.0, 56.0),
    ZEUS: UranianElements("Zeus", 32.1893, 0.532664, 1.0, 141.0, 41.0),
    KRONOS: UranianElements("Kronos", 213.2096, 0.420481, 1.0, 98.0, 32.0),
    APOLLON: UranianElements("Apollon", 71.8925, 0.341403, 1.0, 326.0, 26.0),
    ADMETOS: UranianElements("Admetos", 142.2269, 0.283756, 1.0, 250.0, 22.0),
    VULKANUS: UranianElements("Vulkanus", 195.6753, 0.240116, 1.0, 178.0, 18.0),
    POSEIDON: UranianElements("Poseidon", 274.4073, 0.207016, 1.0, 105.0, 16.0),
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


@dataclass
class TransplutoKeplerianElements(HypotheticalElements):
    """Public rc7-compatible element type for Transpluto (Isis)."""


@dataclass
class CupidoKeplerianElements:
    """Legacy rc7 public element type for Cupido."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


CUPIDO_KEPLERIAN_ELEMENTS = CupidoKeplerianElements(
    name="Cupido",
    epoch=2415020.0,
    a=40.99837,
    e=0.0,
    i=0.0,
    omega=0.0,
    Omega=0.0,
    L0=105.301693,
    n=0.0037945179,
)


@dataclass
class HadesKeplerianElements:
    """Legacy rc7 public element type for Hades."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float


HADES_KEPLERIAN_ELEMENTS = HadesKeplerianElements(
    name="Hades",
    epoch=2415020.0,
    a=50.66744,
    e=0.00245,
    i=1.0500,
    omega=148.1796,
    Omega=161.3339,
    M0=26.850162,
    n=0.00278759,
)


@dataclass
class ZeusKeplerianElements:
    """Legacy rc7 public element type for Zeus."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


ZEUS_KEPLERIAN_ELEMENTS = ZeusKeplerianElements(
    name="Zeus",
    epoch=2415020.0,
    a=59.21436,
    e=0.0,
    i=0.0,
    omega=0.0,
    Omega=0.0,
    L0=104.289095,
    n=0.0022203750,
)


@dataclass
class KronosKeplerianElements:
    """Legacy rc7 public element type for Kronos."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


KRONOS_KEPLERIAN_ELEMENTS = KronosKeplerianElements(
    name="Kronos",
    epoch=2415020.0,
    a=64.81690,
    e=0.0,
    i=0.0,
    omega=0.0,
    Omega=0.0,
    L0=17.111353,
    n=0.0019351856,
)


@dataclass
class ApollonKeplerianElements:
    """Legacy rc7 public element type for Apollon."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


APOLLON_KEPLERIAN_ELEMENTS = ApollonKeplerianElements(
    name="Apollon",
    epoch=2415020.0,
    a=70.361180,
    e=0.0,
    i=0.0,
    omega=0.0,
    Omega=0.0,
    L0=138.565328,
    n=0.0017177599,
)


@dataclass
class AdmetosKeplerianElements:
    """Legacy rc7 public element type for Admetos."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


ADMETOS_KEPLERIAN_ELEMENTS = AdmetosKeplerianElements(
    name="Admetos",
    epoch=2415020.0,
    a=73.736396,
    e=0.0,
    i=0.0,
    omega=0.0,
    Omega=0.0,
    L0=350.613913,
    n=0.0016016766,
)


@dataclass
class VulkanusKeplerianElements:
    """Legacy rc7 public element type for Vulkanus."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


VULKANUS_KEPLERIAN_ELEMENTS = VulkanusKeplerianElements(
    name="Vulkanus",
    epoch=2415020.0,
    a=77.445895,
    e=0.0,
    i=0.0,
    omega=0.0,
    Omega=0.0,
    L0=55.397715,
    n=0.0015069325,
)


@dataclass
class PoseidonKeplerianElements:
    """Legacy rc7 public element type for Poseidon."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float


POSEIDON_KEPLERIAN_ELEMENTS = PoseidonKeplerianElements(
    name="Poseidon",
    epoch=2415020.0,
    a=83.666307,
    e=0.0,
    i=0.0,
    omega=0.0,
    Omega=0.0,
    L0=166.140256,
    n=0.0013256078,
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

    Container reserved for independently reviewed Hamburg-school elements.
    Propagation uses Kepler's equation, a Gaussian (PQR) vector
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


@dataclass
class VulcanElements:
    """Time-dependent heliocentric elements for historical Vulcan ID 55."""

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


@dataclass
class WaldemathElements:
    """Time-dependent geocentric elements for Waldemath's second moon."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    L0: float
    n: float
    M0: float = 0.0
    M_rate: float = 0.0
    omega_rate: float = 0.0
    Omega_rate: float = 0.0
    equinox: float = 0.0


@dataclass
class LowellPlanetXElements:
    """Legacy public element container for Lowell's Planet X."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float


@dataclass
class PickeringPlanetXElements:
    """Legacy public element container for Pickering's Planet O/X."""

    name: str
    epoch: float
    a: float
    e: float
    i: float
    omega: float
    Omega: float
    M0: float
    n: float


# Full Keplerian propagation model:
#   1. Propagate mean anomaly: M = M0 + n*(t - t0) where n = K_GAUSS_DEG/a^1.5
#   2. Solve Kepler's equation: M = E - e*sin(E) -> eccentric anomaly E
#   3. Position in orbital plane: x = a*(cosE - e), y = a*sqrt(1-e²)*sinE
#   4. Transform via Gaussian vectors (PQR matrix) to ecliptic in equinox frame
#   5. Precess the ecliptic direction from the element equinox to J2000 via
#      astrometry._precess_ecliptic (Vondrák 2011 long-term precession, the
#      same realization used downstream for the J2000 -> of-date step)
#
# Gaussian gravitational constant: k = 0.01720209895 rad/day
# -> daily motion = 0.9856076686 deg/day / a^1.5

# Gaussian gravitational constant squared, in degrees/day per AU^1.5
_K_GAUSS_DEG: float = 0.9856076686

URANIAN_KEPLERIAN_ELEMENTS: Dict[int, UranianKeplerianElements] = {
    CUPIDO: UranianKeplerianElements(
        "Cupido",
        2415020.0,
        2415020.0,
        40.99837,
        0.00460,
        1.0833,
        171.4333,
        129.8325,
        163.7409,
        _K_GAUSS_DEG / 40.99837**1.5,
    ),
    HADES: UranianKeplerianElements(
        "Hades",
        2415020.0,
        2415020.0,
        50.66744,
        0.00245,
        1.0500,
        148.1796,
        161.3339,
        27.6496,
        _K_GAUSS_DEG / 50.66744**1.5,
    ),
    ZEUS: UranianKeplerianElements(
        "Zeus",
        2415020.0,
        2415020.0,
        59.21436,
        0.00120,
        0.0,
        299.0440,
        0.0,
        165.1232,
        _K_GAUSS_DEG / 59.21436**1.5,
    ),
    KRONOS: UranianKeplerianElements(
        "Kronos",
        2415020.0,
        2415020.0,
        64.81690,
        0.00305,
        0.0,
        208.8801,
        0.0,
        169.0193,
        _K_GAUSS_DEG / 64.81690**1.5,
    ),
    APOLLON: UranianKeplerianElements(
        "Apollon",
        2415020.0,
        2415020.0,
        70.29949,
        0.0,
        0.0,
        0.0,
        0.0,
        138.0533,
        _K_GAUSS_DEG / 70.29949**1.5,
    ),
    ADMETOS: UranianKeplerianElements(
        "Admetos",
        2415020.0,
        2415020.0,
        73.62765,
        0.0,
        0.0,
        0.0,
        0.0,
        351.3350,
        _K_GAUSS_DEG / 73.62765**1.5,
    ),
    VULKANUS: UranianKeplerianElements(
        "Vulkanus",
        2415020.0,
        2415020.0,
        77.25568,
        0.0,
        0.0,
        0.0,
        0.0,
        55.8983,
        _K_GAUSS_DEG / 77.25568**1.5,
    ),
    POSEIDON: UranianKeplerianElements(
        "Poseidon",
        2415020.0,
        2415020.0,
        83.66907,
        0.0,
        0.0,
        0.0,
        0.0,
        165.5163,
        _K_GAUSS_DEG / 83.66907**1.5,
    ),
}

# Strubell, "Die Sterne" 3/1952, p. 70ff: epoch JD 2368547.66,
# equinox JD 2431456.5, a=77.775 AU, e=0.3, i=0, perihelion argument 0.7.
# The epoch and angular elements below are that orbit propagated/precessed to
# J2000 by LibEphemeris.
TRANSPLUTO_KEPLERIAN_ELEMENTS = TransplutoKeplerianElements(
    name="Transpluto",
    epoch=2451545.0,
    a=77.775,
    e=0.3,
    i=0.0,
    omega=1.468282,
    Omega=0.0,
    M0=119.265936,
    n=360.0 / (77.775**1.5 * 365.25),
)

HYPOTHETICAL_ELEMENTS: Dict[int, HypotheticalElements] = {
    ISIS: TRANSPLUTO_KEPLERIAN_ELEMENTS,
    PROSERPINA: HypotheticalElements(
        name="Proserpina",
        epoch=2415020.0,
        a=79.225630,
        e=0.0,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        M0=170.73,
        n=_K_GAUSS_DEG / 79.225630**1.5,
    ),
}

# These public constants are restored for source compatibility.  Runtime
# Lowell/Pickering calculations use the canonical historical table below.
LOWELL_PLANET_X_ELEMENTS = LowellPlanetXElements(
    "Planet X Lowell",
    2415020.0,
    43.0,
    0.202,
    10.0,
    204.9,
    67.5,
    11.2,
    360.0 / (43.0**1.5 * 365.25),
)
PICKERING_PLANET_X_ELEMENTS = PickeringPlanetXElements(
    "Planet X Pickering",
    2415020.0,
    51.9,
    0.31,
    15.0,
    100.0,
    110.0,
    15.0,
    360.0 / (51.9**1.5 * 365.25),
)

# Historical convention (L. H. Weston); the documented element source is
# data/fictitious_orbits.csv (row Vulcan).
VULCAN_ELEMENTS = VulcanElements(
    "Vulcan",
    2415020.0,
    0.13744,
    0.019,
    7.5,
    252.8987988,
    707550.7341,
    322.212069,
    1670.056,
    47.787931,
    -1670.056,
)

# Historical convention (Waldemath 1898 claim); the documented element
# source is data/fictitious_orbits.csv (row Waldemath).
WALDEMATH_ELEMENTS = WaldemathElements(
    "Waldemath",
    2414290.95827875,
    0.0068400705250028,
    0.1587,
    2.5,
    8.14049594,
    136.24878256,
    (70.3407215 + 8.14049594 + 136.24878256) % 360.0,
    109023.2634989 / 36525.0,
    70.3407215,
    109023.2634989,
    2393.47417444,
    -1131.71719709,
    2414290.95827875,
)

HYPOTHETICAL_PROVENANCE: Dict[int, Tuple[str, str]] = {
    CUPIDO: ("built-in", "Witte 1923; Neely 1988 refinement"),
    HADES: ("built-in", "Witte 1924; Neely 1988 refinement"),
    ZEUS: ("built-in", "Witte 1924; Neely 1988 refinement"),
    KRONOS: ("built-in", "Witte 1924; Neely 1988 refinement"),
    APOLLON: ("built-in", "Sieggruen 1929/1937"),
    ADMETOS: ("built-in", "Sieggruen 1929/1937"),
    VULKANUS: ("built-in", "Sieggruen 1929/1937"),
    POSEIDON: ("built-in", "Sieggruen 1929/1937"),
    ISIS: ("built-in", 'Strubell, "Die Sterne" 3/1952 p.70ff'),
    NIBIRU: ("built-in", "conventional element set"),
    HARRINGTON: ("built-in", "Harrington, AJ 96 (1988), p.1478"),
    NEPTUNE_LEVERRIER: ("built-in", "Le Verrier, CRAS 23 (1846)"),
    NEPTUNE_ADAMS: ("built-in", "Adams, 1846 prediction"),
    PLUTO_LOWELL: ("built-in", "Lowell, 1915 memoir"),
    PLUTO_PICKERING: ("built-in", "Pickering, 1919/1928 predictions"),
    VULCAN: ("built-in", "L. H. Weston historical element convention"),
    WHITE_MOON: (
        "built-in",
        "Russian-school seven-year Selena convention",
    ),
    PROSERPINA: ("built-in", "V. Abramov convention"),
    WALDEMATH: ("built-in", "Waldemath 1898 claim"),
}


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
    WHITE_MOON: "White Moon",
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
    polynomials in T, e.g., "12.5 + 3456.75 * T".

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
        - "12.5"
        - "12.5 + 3456.75 * T"
        - "98.25+432.5*T"
        - "47.5-84.5*T"
        - "47.5 - 84.5 * T"

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
    """Parse a custom hypothetical-body orbital-elements file.

    Data lines contain nine comma-separated fields: epoch, equinox, mean
    anomaly, semi-major axis, eccentricity, argument of perihelion, ascending
    node, inclination, and body name. Blank lines and lines beginning with
    ``#`` are ignored. Supported elements may use a linear polynomial in
    Julian centuries, for example ``12.5 + 3456.75 * T``.

    Args:
        filepath: Path to the orbital-elements file.

    Returns:
        Parsed ``OrbitalElements`` objects in file order.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If a data line cannot be parsed.

    Example:
        >>> elements = parse_orbital_elements("my_orbits.txt")
        >>> [element.name for element in elements]
        ['MyPlanet']
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

    The dataset contains the built-in element sets for IDs 40--58.

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

    The dataset (``data/fictitious_orbits.csv``) contains one independently
    transcribed element set. All downstream helpers
    (``get_orbital_body_by_name``, ``calc_orbital_position``, etc.) work
    without modification.

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
        >>> harrington = get_orbital_body_by_name(bodies, "Harrington")
        >>> print(f"Harrington semi-axis: {harrington.semi_axis} AU")
        Harrington semi-axis: 101.2 AU

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

    The old dataset has been retired. The clean-room Harrington dataset
    (``data/fictitious_orbits.csv``) is returned instead. This function delegates to
    :func:`get_bundled_fictitious_orbits_path` for backward compatibility.

    Returns:
        Path to ``data/fictitious_orbits.csv``.
    """
    return get_bundled_fictitious_orbits_path()


def load_bundled_seorbel() -> List[OrbitalElements]:
    """
    Deprecated.  Loads the bundled fictitious orbits dataset.

    The old dataset has been retired. The clean-room Harrington dataset
    (``data/fictitious_orbits.csv``) is loaded instead. This function delegates to
    :func:`load_bundled_fictitious_orbits` for backward compatibility.

    Returns:
        List of :class:`OrbitalElements` objects for all bundled bodies.
    """
    return load_bundled_fictitious_orbits()


# Backward-compatible aliases for legacy public names
# These aliases are maintained for backward compatibility with existing code.
# The canonical names are the reference-independent versions used throughout this module.


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

    # Calculate velocity via central difference numerical differentiation.
    # The bodies routed here (the predicted trans-Neptunian planets) are slow
    # (<<1 deg/day), so a 1-day step gives the true derivative of the reported
    # position; the fast short-period bodies (Vulcan, Waldemath) have their own
    # samplers with tighter steps.
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


def _keplerian_to_ecliptic_j2000(
    elements: UranianKeplerianElements, jd_tt: float
) -> Tuple[float, float, float]:
    """Propagate fixed Keplerian elements into the J2000 ecliptic frame.

    This is a direct implementation of the standard two-body perifocal-to-
    ecliptic rotation. It is shared by reviewed or caller-supplied element
    sets; no external ephemeris implementation or coefficient table is
    involved.
    """
    mean_anomaly = math.radians(
        (elements.M0 + elements.n * (jd_tt - elements.epoch)) % 360.0
    )
    eccentric_anomaly = _solve_kepler_equation(mean_anomaly, elements.e)

    # Cartesian position in the orbital plane, with the x-axis at periapsis.
    x_orbit = elements.a * (math.cos(eccentric_anomaly) - elements.e)
    y_orbit = (
        elements.a
        * math.sqrt(max(0.0, 1.0 - elements.e * elements.e))
        * math.sin(eccentric_anomaly)
    )

    argument = math.radians(elements.omega)
    node = math.radians(elements.Omega)
    inclination = math.radians(elements.i)
    cos_arg, sin_arg = math.cos(argument), math.sin(argument)
    cos_node, sin_node = math.cos(node), math.sin(node)
    cos_inc, sin_inc = math.cos(inclination), math.sin(inclination)

    # R3(node) R1(inclination) R3(argument) applied to (x_orbit, y_orbit, 0).
    x_ecliptic = (cos_node * cos_arg - sin_node * sin_arg * cos_inc) * x_orbit + (
        -cos_node * sin_arg - sin_node * cos_arg * cos_inc
    ) * y_orbit
    y_ecliptic = (sin_node * cos_arg + cos_node * sin_arg * cos_inc) * x_orbit + (
        -sin_node * sin_arg + cos_node * cos_arg * cos_inc
    ) * y_orbit
    z_ecliptic = sin_arg * sin_inc * x_orbit + cos_arg * sin_inc * y_orbit

    distance = math.sqrt(
        x_ecliptic * x_ecliptic + y_ecliptic * y_ecliptic + z_ecliptic * z_ecliptic
    )
    longitude = math.degrees(math.atan2(y_ecliptic, x_ecliptic)) % 360.0
    latitude = math.degrees(math.atan2(z_ecliptic, math.hypot(x_ecliptic, y_ecliptic)))

    if abs(elements.equinox_jd - 2451545.0) > 1e-6:
        from .astrometry import _precess_ecliptic

        longitude, latitude = _precess_ecliptic(
            longitude, latitude, elements.equinox_jd, 2451545.0
        )
    return (float(longitude), float(latitude), float(distance))


# Built-in prediction elements for bodies 49--54.
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
    """Return the heliocentric position of a classical prediction.

    Propagates the elements in ``FICTITIOUS_ORBITAL_ELEMENTS``. Per-body source
    classifications are recorded in ``HYPOTHETICAL_PROVENANCE``.

    Returns (lon, lat, dist, dlon, dlat, ddist), J2000 ecliptic.
    """
    elem = FICTITIOUS_ORBITAL_ELEMENTS.get(ipl)
    if elem is None:
        from .exceptions import UnknownBodyError

        raise UnknownBodyError(
            message=f"no fictitious elements for body {ipl}", body_id=ipl
        )

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
    """Calculate the longitude of a restored Hamburg-school body.

    Args:
        ipl: Uranian planet ID (CUPIDO through POSEIDON)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Raises:
        UnknownBodyError: The body has no registered Uranian element set.
    """
    if ipl not in URANIAN_KEPLERIAN_ELEMENTS:
        _raise_unsupported_builtin_fictitious(ipl)

    # Use the single consolidated Keplerian element set shared by all Uranian
    # calculation entry points.
    return calc_uranian_planet(ipl, jd_tt)[0]


def calc_uranian_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """Calculate the state of a restored Hamburg-school body."""
    if ipl not in URANIAN_KEPLERIAN_ELEMENTS:
        _raise_unsupported_builtin_fictitious(ipl)
    return calc_uranian_planet(ipl, jd_tt)


def calc_cupido(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored Cupido state."""
    return calc_uranian_planet(CUPIDO, jd_tt)


def calc_hades(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored Hades state."""
    return calc_uranian_planet(HADES, jd_tt)


def calc_zeus(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored Zeus state."""
    return calc_uranian_planet(ZEUS, jd_tt)


def calc_kronos(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored Kronos state."""
    return calc_uranian_planet(KRONOS, jd_tt)


def calc_apollon(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored Apollon state."""
    return calc_uranian_planet(APOLLON, jd_tt)


def calc_admetos(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored Admetos state."""
    return calc_uranian_planet(ADMETOS, jd_tt)


def calc_vulkanus(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored Vulkanus state."""
    return calc_uranian_planet(VULKANUS, jd_tt)


def calc_poseidon(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored Poseidon state."""
    return calc_uranian_planet(POSEIDON, jd_tt)


def calc_transpluto(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Return the heliocentric J2000 state of Strubell's Transpluto."""
    longitude, latitude, distance = _calc_transpluto_raw(jd_tt)
    step = 1.0
    previous = _calc_transpluto_raw(jd_tt - step)
    following = _calc_transpluto_raw(jd_tt + step)
    dlon = ((following[0] - previous[0] + 180.0) % 360.0 - 180.0) / (2.0 * step)
    dlat = (following[1] - previous[1]) / (2.0 * step)
    ddist = (following[2] - previous[2]) / (2.0 * step)
    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_transpluto_raw(jd_tt: float) -> Tuple[float, float, float]:
    """Return Transpluto position without numerical velocity."""
    return _calc_keplerian_position_raw(ISIS, jd_tt)


def calc_uranian_planet(
    body_id: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """Calculate a restored Hamburg-school body.

    The independent propagation implementation uses Kepler's equation,
    Gaussian-vector rotation, and central-difference velocity. The compatibility
    element registry covers Cupido through Poseidon; its per-body source status
    is recorded in ``HYPOTHETICAL_PROVENANCE``.

    Args:
        body_id: Uranian planet ID (CUPIDO through POSEIDON, i.e., 40-47)
        jd_tt: Julian Day in Terrestrial Time (TT)

    Raises:
        UnknownBodyError: The body has no registered Uranian element set.
    """
    if body_id not in URANIAN_KEPLERIAN_ELEMENTS:
        _raise_unsupported_builtin_fictitious(body_id)

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


def _raise_unsupported_builtin_fictitious(body_id: int) -> NoReturn:
    """Raise for an ID absent from the restored element registries."""
    from .exceptions import UnknownBodyError

    raise UnknownBodyError(
        message=f"no built-in fictitious elements for body {body_id}",
        body_id=body_id,
    )


def calc_vulcan(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the historical time-dependent Vulcan convention."""
    longitude, latitude, distance = _calc_vulcan_raw(jd_tt)
    step = 0.01
    previous = _calc_vulcan_raw(jd_tt - step)
    following = _calc_vulcan_raw(jd_tt + step)
    dlon = ((following[0] - previous[0] + 180.0) % 360.0 - 180.0) / (2.0 * step)
    dlat = (following[1] - previous[1]) / (2.0 * step)
    ddist = (following[2] - previous[2]) / (2.0 * step)
    return (longitude, latitude, distance, dlon, dlat, ddist)


def _calc_vulcan_raw(jd_tt: float) -> Tuple[float, float, float]:
    """Return Vulcan's heliocentric ecliptic position without velocity."""
    elements = VULCAN_ELEMENTS
    centuries = (jd_tt - elements.epoch) / 36525.0
    mean_anomaly = math.radians((elements.M0 + elements.n_century * centuries) % 360.0)
    argument = math.radians(elements.omega0 + elements.omega_rate * centuries)
    node = math.radians(elements.Omega0 + elements.Omega_rate * centuries)
    eccentric_anomaly = _solve_kepler_equation(mean_anomaly, elements.e)
    true_anomaly = 2.0 * math.atan(
        math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
        * math.tan(eccentric_anomaly / 2.0)
    )
    distance = elements.a * (1.0 - elements.e * math.cos(eccentric_anomaly))
    argument_of_latitude = true_anomaly + argument
    inclination = math.radians(elements.i)
    x_orbit = distance * math.cos(argument_of_latitude)
    y_orbit = distance * math.sin(argument_of_latitude)
    cos_i, sin_i = math.cos(inclination), math.sin(inclination)
    cos_node, sin_node = math.cos(node), math.sin(node)
    x = cos_node * x_orbit - sin_node * cos_i * y_orbit
    y = sin_node * x_orbit + cos_node * cos_i * y_orbit
    z = sin_i * y_orbit
    longitude = math.degrees(math.atan2(y, x)) % 360.0
    latitude = math.degrees(math.asin(z / distance))
    return (longitude, latitude, distance)


def calc_waldemath(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the historical geocentric Waldemath convention."""
    from .astrometry import _precess_ecliptic

    elements = WALDEMATH_ELEMENTS
    equinox = elements.equinox or elements.epoch

    def raw(jd: float) -> Tuple[float, float, float]:
        centuries = (jd - elements.epoch) / 36525.0
        mean_anomaly = math.radians((elements.M0 + elements.M_rate * centuries) % 360.0)
        argument = math.radians(
            (elements.omega + elements.omega_rate * centuries) % 360.0
        )
        node = math.radians((elements.Omega + elements.Omega_rate * centuries) % 360.0)
        eccentric_anomaly = _solve_kepler_equation(mean_anomaly, elements.e)
        true_anomaly = 2.0 * math.atan(
            math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
            * math.tan(eccentric_anomaly / 2.0)
        )
        distance = elements.a * (1.0 - elements.e * math.cos(eccentric_anomaly))
        argument_of_latitude = true_anomaly + argument
        inclination = math.radians(elements.i)
        x_orbit = distance * math.cos(argument_of_latitude)
        y_orbit = distance * math.sin(argument_of_latitude)
        cos_i, sin_i = math.cos(inclination), math.sin(inclination)
        cos_node, sin_node = math.cos(node), math.sin(node)
        x = cos_node * x_orbit - sin_node * cos_i * y_orbit
        y = sin_node * x_orbit + cos_node * cos_i * y_orbit
        z = sin_i * y_orbit
        longitude = math.degrees(math.atan2(y, x)) % 360.0
        latitude = math.degrees(math.asin(max(-1.0, min(1.0, z / distance))))
        longitude, latitude = _precess_ecliptic(longitude, latitude, equinox, jd)
        return (longitude, latitude, distance)

    longitude, latitude, distance = raw(jd_tt)
    step = 0.01
    previous = raw(jd_tt - step)
    following = raw(jd_tt + step)
    dlon = ((following[0] - previous[0] + 180.0) % 360.0 - 180.0) / (2.0 * step)
    dlat = (following[1] - previous[1]) / (2.0 * step)
    ddist = (following[2] - previous[2]) / (2.0 * step)
    return (longitude, latitude, distance, dlon, dlat, ddist)


def calc_transpluto_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate the position of Transpluto (Isis/Persephone)."""
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
    """Return the seven-year White Moon/Selena convention.

    ``use_true_lilith`` is retained for API compatibility and ignored.
    """
    del use_true_lilith
    centuries = (jd_tt - 2451545.0) / 36525.0
    longitude = (242.2205555 + 5143.5418158 * centuries) % 360.0
    longitude_speed = 5143.5418158 / 36525.0
    return (longitude, 0.0, 0.05280098949, longitude_speed, 0.0, 0.0)


def calc_waldemath_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate Waldemath's historical second-moon convention."""
    return calc_waldemath(jd_tt)


def calc_proserpina(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the restored circular Proserpina compatibility convention."""
    elements = HYPOTHETICAL_ELEMENTS[PROSERPINA]
    longitude = (elements.M0 + elements.n * (jd_tt - elements.epoch)) % 360.0
    return (longitude, 0.0, elements.a, elements.n, 0.0, 0.0)


def calc_planet_x_lowell(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate Lowell's historical Planet X prediction."""
    return calc_fictitious_position(PLUTO_LOWELL, jd_tt)


def calc_planet_x_pickering(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate Pickering's historical Planet O/X prediction."""
    return calc_fictitious_position(PLUTO_PICKERING, jd_tt)


def calc_hypothetical_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of any hypothetical body.

    This is the main entry point for hypothetical body calculations.
    Routes the request to the appropriate calculation function based on body ID.

    Args:
        ipl: Hypothetical body ID (CUPIDO through PLUTO_PICKERING)
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

    if ipl in URANIAN_KEPLERIAN_ELEMENTS:
        return calc_uranian_planet(ipl, jd_tt)

    if ipl in FICTITIOUS_ORBITAL_ELEMENTS:
        return calc_fictitious_position(ipl, jd_tt)

    if ipl == ISIS:
        return calc_transpluto_position(jd_tt)
    if ipl == VULCAN:
        return calc_vulcan(jd_tt)
    if ipl == WHITE_MOON:
        return calc_white_moon_position(jd_tt)
    if ipl == PROSERPINA:
        return calc_proserpina(jd_tt)
    if ipl == WALDEMATH:
        return calc_waldemath_position(jd_tt)
    if ipl in HYPOTHETICAL_ELEMENTS:
        return _calc_keplerian_position(ipl, jd_tt)

    _raise_unsupported_builtin_fictitious(ipl)


def list_hypothetical_bodies() -> Dict[int, str]:
    """
    Get a dictionary of recognised hypothetical body IDs and names.

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
# Backward-compatible aliases (legacy public names)
# ---------------------------------------------------------------------------
# These aliases preserve backward compatibility for code using the old naming
# convention. New code should use the canonical names:
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
