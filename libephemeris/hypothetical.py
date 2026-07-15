# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Independently sourced hypothetical-body calculations.

The historical bodies exposed as IDs 40--58 are mathematical conventions, not
claims that the objects exist.  A built-in calculation is enabled only when
every numerical field can be transcribed or derived reproducibly from an
identified public source.  IDs whose complete element set could not be
recovered remain importable and discoverable for API compatibility, but fail
closed with :class:`~libephemeris.exceptions.UnknownBodyError`.

The supported set is deliberately explicit:

* IDs 40--47 use James Neely's published J1900 Hamburg-school elements;
* ID 50 uses Harrington's 1988 nominal Planet-X orbit;
* IDs 51--52 use the primary 1846 predictions of Le Verrier and Adams; and
* ID 53 uses Lowell's preferred 1915 Planet-X solution; and
* ID 56 is the published uniform seven-year Selena convention.

The source transcription, frame choices and arithmetic transformations are
documented beside the constants below, in ``fictitious_orbits.csv`` and in
``docs/methodology/hypothetical-bodies.md``.  Missing fields for disabled IDs
are inventoried in ``docs/methodology/missing-hypothetical-models.md``.
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


# Public element-container objects for *unsupported* bodies use NaN fields.
# This preserves their import names and type identity without silently
# distributing, approximating, or calculating from an unverified data table.
_UNAVAILABLE_FLOAT = float("nan")


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
# NEELY 1980 PRIMARY TRANSCRIPTION (IDS 40--47)
# =============================================================================
@dataclass(frozen=True)
class _NeelySourceRow:
    """One literal row from Neely's published Table I.

    ``mean_anomaly_rate_century`` is retained even though the runtime derives
    mean motion from ``a``.  Keeping the printed rate makes the transcription
    independently reviewable and lets the provenance gate verify that the
    semimajor-axis propagation agrees with the author's rounded rate.
    """

    name: str
    a: float
    e: float
    i: float
    omega: float
    node: float
    mean_anomaly: float
    mean_anomaly_rate_century: float

    @property
    def mean_longitude(self) -> float:
        """Return M + omega + node at the printed epoch, in degrees."""
        return (self.mean_anomaly + self.omega + self.node) % 360.0


# Primary source and transcription policy
# ---------------------------------------
# James Neely, "Orbital Elements for the Transneptunian Planets", Matrix
# Magazine, issue VII (1980), pp. 6--12, Table I on p. 10.  A public scan is
# catalogued as item 18 by Alexandria iBase:
#
#   https://special-collections.alexandriaibase.org/items/show/18
#
# SHA-256 of the reviewed 63-page scan:
#   f3d2e834a4a684736383ec387c9c2786064ca4668256724865da2f3d8fa03020
#
# Neely defines T as Julian centuries from JD 2415020.0.  Table I prints the
# seven fields represented below plus a common node term +1.3960*T.  That term
# expresses the moving equinox used by the source.  LibEphemeris instead keeps
# the T=0 angles in the fixed J1900 ecliptic frame and applies its documented
# IAU/Vondrak precession exactly once in ``_keplerian_to_ecliptic_j2000``.
# Adding Neely's moving-equinox term as well would double-count precession.
#
# Apollon, Admetos, Vulkanus and Poseidon are circular and coplanar in Table I.
# For e=i=0, perihelion and node are geometrically degenerate: only the sum
# M+omega+node affects position.  We retain Neely's literal convention M=0 and
# phase in omega.  Earlier representations that moved the phase into M were
# numerically equivalent but were not literal transcriptions.
_NEELY_EPOCH_JD = 2415020.0
_NEELY_SOURCE_ROWS: Dict[int, _NeelySourceRow] = {
    CUPIDO: _NeelySourceRow(
        "Cupido", 40.99837, 0.00460, 1.0833, 171.4333, 129.8325, 163.7409, 137.13640
    ),
    HADES: _NeelySourceRow(
        "Hades", 50.66744, 0.00245, 1.0500, 148.1796, 161.3339, 27.6496, 99.81803
    ),
    ZEUS: _NeelySourceRow(
        "Zeus", 59.21436, 0.00120, 0.0, 299.0440, 0.0, 165.1232, 79.00633
    ),
    KRONOS: _NeelySourceRow(
        "Kronos", 64.81690, 0.00305, 0.0, 208.8801, 0.0, 169.0193, 68.98746
    ),
    APOLLON: _NeelySourceRow(
        "Apollon", 70.29949, 0.0, 0.0, 138.0533, 0.0, 0.0, 61.0765
    ),
    ADMETOS: _NeelySourceRow(
        "Admetos", 73.62765, 0.0, 0.0, 351.3350, 0.0, 0.0, 56.98245
    ),
    VULKANUS: _NeelySourceRow(
        "Vulkanus", 77.25568, 0.0, 0.0, 55.8983, 0.0, 0.0, 53.015987
    ),
    POSEIDON: _NeelySourceRow(
        "Poseidon", 83.66907, 0.0, 0.0, 165.5163, 0.0, 0.0, 47.03868
    ),
}


# =============================================================================
# URANIAN PLANET PARAMETERS (Hamburg School)
# =============================================================================
@dataclass
class UranianElements:
    """Legacy mean-longitude view of a reviewed Neely element row.

    ``L0`` is the source row's mean longitude at JD 2415020.0 and ``n`` is
    Neely's printed mean-anomaly rate converted from degrees per Julian century
    to degrees per day.  The three sinusoid fields are zero because the retired
    empirical correction table had no independently established provenance;
    no such correction is present in Neely's primary element table.
    """

    name: str
    L0: float
    n: float
    amplitude: float = 0.0
    phase: float = 0.0
    phase_rate: float = 0.0


URANIAN_ELEMENTS: Dict[int, UranianElements] = {
    body_id: UranianElements(
        name=row.name,
        L0=row.mean_longitude,
        n=row.mean_anomaly_rate_century / 36525.0,
        amplitude=0.0,
        phase=0.0,
        phase_rate=0.0,
    )
    for body_id, row in _NEELY_SOURCE_ROWS.items()
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
    epoch=_NEELY_EPOCH_JD,
    a=_NEELY_SOURCE_ROWS[CUPIDO].a,
    e=_NEELY_SOURCE_ROWS[CUPIDO].e,
    i=_NEELY_SOURCE_ROWS[CUPIDO].i,
    omega=_NEELY_SOURCE_ROWS[CUPIDO].omega,
    Omega=_NEELY_SOURCE_ROWS[CUPIDO].node,
    L0=_NEELY_SOURCE_ROWS[CUPIDO].mean_longitude,
    n=_NEELY_SOURCE_ROWS[CUPIDO].mean_anomaly_rate_century / 36525.0,
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
    epoch=_NEELY_EPOCH_JD,
    a=_NEELY_SOURCE_ROWS[HADES].a,
    e=_NEELY_SOURCE_ROWS[HADES].e,
    i=_NEELY_SOURCE_ROWS[HADES].i,
    omega=_NEELY_SOURCE_ROWS[HADES].omega,
    Omega=_NEELY_SOURCE_ROWS[HADES].node,
    M0=_NEELY_SOURCE_ROWS[HADES].mean_anomaly,
    n=_NEELY_SOURCE_ROWS[HADES].mean_anomaly_rate_century / 36525.0,
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
    epoch=_NEELY_EPOCH_JD,
    a=_NEELY_SOURCE_ROWS[ZEUS].a,
    e=_NEELY_SOURCE_ROWS[ZEUS].e,
    i=_NEELY_SOURCE_ROWS[ZEUS].i,
    omega=_NEELY_SOURCE_ROWS[ZEUS].omega,
    Omega=_NEELY_SOURCE_ROWS[ZEUS].node,
    L0=_NEELY_SOURCE_ROWS[ZEUS].mean_longitude,
    n=_NEELY_SOURCE_ROWS[ZEUS].mean_anomaly_rate_century / 36525.0,
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
    epoch=_NEELY_EPOCH_JD,
    a=_NEELY_SOURCE_ROWS[KRONOS].a,
    e=_NEELY_SOURCE_ROWS[KRONOS].e,
    i=_NEELY_SOURCE_ROWS[KRONOS].i,
    omega=_NEELY_SOURCE_ROWS[KRONOS].omega,
    Omega=_NEELY_SOURCE_ROWS[KRONOS].node,
    L0=_NEELY_SOURCE_ROWS[KRONOS].mean_longitude,
    n=_NEELY_SOURCE_ROWS[KRONOS].mean_anomaly_rate_century / 36525.0,
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
    epoch=_NEELY_EPOCH_JD,
    a=_NEELY_SOURCE_ROWS[APOLLON].a,
    e=_NEELY_SOURCE_ROWS[APOLLON].e,
    i=_NEELY_SOURCE_ROWS[APOLLON].i,
    omega=_NEELY_SOURCE_ROWS[APOLLON].omega,
    Omega=_NEELY_SOURCE_ROWS[APOLLON].node,
    L0=_NEELY_SOURCE_ROWS[APOLLON].mean_longitude,
    n=_NEELY_SOURCE_ROWS[APOLLON].mean_anomaly_rate_century / 36525.0,
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
    epoch=_NEELY_EPOCH_JD,
    a=_NEELY_SOURCE_ROWS[ADMETOS].a,
    e=_NEELY_SOURCE_ROWS[ADMETOS].e,
    i=_NEELY_SOURCE_ROWS[ADMETOS].i,
    omega=_NEELY_SOURCE_ROWS[ADMETOS].omega,
    Omega=_NEELY_SOURCE_ROWS[ADMETOS].node,
    L0=_NEELY_SOURCE_ROWS[ADMETOS].mean_longitude,
    n=_NEELY_SOURCE_ROWS[ADMETOS].mean_anomaly_rate_century / 36525.0,
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
    epoch=_NEELY_EPOCH_JD,
    a=_NEELY_SOURCE_ROWS[VULKANUS].a,
    e=_NEELY_SOURCE_ROWS[VULKANUS].e,
    i=_NEELY_SOURCE_ROWS[VULKANUS].i,
    omega=_NEELY_SOURCE_ROWS[VULKANUS].omega,
    Omega=_NEELY_SOURCE_ROWS[VULKANUS].node,
    L0=_NEELY_SOURCE_ROWS[VULKANUS].mean_longitude,
    n=_NEELY_SOURCE_ROWS[VULKANUS].mean_anomaly_rate_century / 36525.0,
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
    epoch=_NEELY_EPOCH_JD,
    a=_NEELY_SOURCE_ROWS[POSEIDON].a,
    e=_NEELY_SOURCE_ROWS[POSEIDON].e,
    i=_NEELY_SOURCE_ROWS[POSEIDON].i,
    omega=_NEELY_SOURCE_ROWS[POSEIDON].omega,
    Omega=_NEELY_SOURCE_ROWS[POSEIDON].node,
    L0=_NEELY_SOURCE_ROWS[POSEIDON].mean_longitude,
    n=_NEELY_SOURCE_ROWS[POSEIDON].mean_anomaly_rate_century / 36525.0,
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
    """Legacy public container for historical Vulcan ID 55.

    The exported instance is deliberately populated with ``NaN`` sentinels:
    no source-complete Weston convention has been recovered, so these fields
    must never be treated as a calculable orbit.  The type remains available
    solely to avoid removing an established import name.
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


@dataclass
class WaldemathElements:
    """Legacy public container for the unsupported Waldemath ID 58.

    Its exported instance contains only ``NaN`` sentinels.  This preserves
    import and attribute compatibility while making accidental numerical use
    fail visibly instead of resurrecting an unverified historical table.
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


# =============================================================================
# WHITE MOON / SELENA PUBLISHED SEVEN-YEAR CONVENTION (ID 56)
# =============================================================================
# Primary source:
#
#   Felix Velichko and Max Larin, "Lilith, Selena, Proserpina: Articles and
#   Ephemerides 1800--2050" (Mir Uranii, Moscow, 2007), ISBN
#   5-900191-12-5, pp. 17, 18, 20, 29, and 45.
#
#   https://djvu.online/file/lBJKQiKf9n4Gd
#
# The implementation below is a fresh reconstruction from the publication,
# not a retained compatibility coefficient:
#
# * p. 17 says that Selena's cycle is seven years and that the point may be
#   treated as moving uniformly through the zodiac;
# * p. 18 discusses a different "seven years minus fifteen days" cycle and
#   explicitly says that this alternative is *not* used in the printed
#   ephemerides, so it must not be substituted here;
# * p. 20 prints January 1800 as 6 degrees 08 arcminutes Taurus, i.e.
#   30 + 6 + 8/60 = 36.133333... degrees ecliptic longitude;
# * p. 45 prints January 2000 as 2 degrees 18 arcminutes Sagittarius, i.e.
#   240 + 2 + 18/60 = 242.3 degrees ecliptic longitude.
#
# A monthly table does not identify a sub-day instant or a dynamical time
# scale.  LibEphemeris therefore makes the missing convention explicit rather
# than hiding it in a fitted phase: a month row is interpreted at 00:00 TT on
# its first civil day.  The January 1800 and January 2000 rows are consequently
# JD 2378496.5 TT and JD 2451544.5 TT.  Any fixed alternative sample instant
# within January would cancel from their 200-year time difference, so this
# declared convention affects the absolute phase but not the derived rate.
#
# The two source longitudes and the stated approximately seven-year cycle
# determine the unwrapped angular displacement without any reference-engine
# output.  Two hundred years contain 28 complete seven-year turns and the
# printed residual arc from 36 deg 08 arcmin to 242 deg 18 arcmin:
#
#   elapsed days = 2451544.5 - 2378496.5 = 73048
#   displacement = 28*360 + (242 deg 18' - 36 deg 08')
#                = 10286 deg 10'
#   rate = displacement / elapsed days
#
# Choosing 29 complete turns would imply about 6.76 years per turn and is
# excluded by p. 17.  The 28-turn result has a period of 2556.567558... days,
# or about 7.000 years, and is therefore the unique unwrap consistent with the
# publication.  It also independently predicts three rows that are *not* used
# in the derivation: March 1879 (p. 29), December 2000 (p. 45), and January
# 2007 (p. 45), each within 0.36 arcminute of the printed value.  That is inside
# the table's one-arcminute resolution and is a materially stronger check than
# assigning a modern year length that the authors never specified.
#
# No public reference-API output is used, persisted, or fitted here.  In
# particular, neither phase nor rate is inherited from a legacy ephemeris
# table: the source checkpoints, full-turn count, calendar-day difference, and
# every arithmetic transformation are written out above.
#
# The source treats Selena as a symbolic point, so it supplies no physical
# distance.  The six-component compatibility API nevertheless requires a
# finite positive radial coordinate.  We use a clearly labelled *display
# convention*: the circular-orbit radius whose period is the source-derived
# 2556.567558... days under the IAU 2015 Resolution B3 nominal terrestrial mass
# parameter.  The SI-to-AU conversion uses the exact astronomical unit adopted
# by IAU 2012 Resolution B2.  This radius is reproducible from public standards
# and is not a claim that Selena is a material Earth satellite.
_WHITE_MOON_EARLY_SOURCE_EPOCH_JD: float = 2378496.5
_WHITE_MOON_EARLY_SOURCE_LONGITUDE_DEG: float = 30.0 + 6.0 + 8.0 / 60.0
_WHITE_MOON_SOURCE_EPOCH_JD: float = 2451544.5
_WHITE_MOON_SOURCE_LONGITUDE_DEG: float = 240.0 + 2.0 + 18.0 / 60.0
_WHITE_MOON_COMPLETE_TURNS: int = 28
_WHITE_MOON_SOURCE_ELAPSED_DAYS: float = (
    _WHITE_MOON_SOURCE_EPOCH_JD - _WHITE_MOON_EARLY_SOURCE_EPOCH_JD
)
_WHITE_MOON_SOURCE_DISPLACEMENT_DEG: float = (
    _WHITE_MOON_COMPLETE_TURNS * 360.0
    + _WHITE_MOON_SOURCE_LONGITUDE_DEG
    - _WHITE_MOON_EARLY_SOURCE_LONGITUDE_DEG
)
_WHITE_MOON_RATE_DEG_PER_DAY: float = (
    _WHITE_MOON_SOURCE_DISPLACEMENT_DEG / _WHITE_MOON_SOURCE_ELAPSED_DAYS
)
_WHITE_MOON_PERIOD_DAYS: float = 360.0 / _WHITE_MOON_RATE_DEG_PER_DAY

_NOMINAL_EARTH_GM_M3_S2: float = 3.986004e14
_ASTRONOMICAL_UNIT_M: float = 149_597_870_700.0
_SECONDS_PER_DAY: float = 86_400.0
_WHITE_MOON_DISTANCE_AU: float = (
    _NOMINAL_EARTH_GM_M3_S2
    * (_WHITE_MOON_PERIOD_DAYS * _SECONDS_PER_DAY / (2.0 * math.pi)) ** 2
) ** (1.0 / 3.0) / _ASTRONOMICAL_UNIT_M

# Runtime projection of the literal Neely rows.
#
# Mean motion is derived from the *transcribed semimajor axis* and the Gaussian
# gravitational constant instead of copying the rounded rate column.  This is
# the standard two-body relation n=k/a^(3/2), gives one internally consistent
# orbit, and differs from every printed Neely rate by less than
# 0.0027 deg/century (the expected effect of the table's rounded ``a`` and
# rate).  The primary rate remains in ``_NEELY_SOURCE_ROWS`` for audit checks.
URANIAN_KEPLERIAN_ELEMENTS: Dict[int, UranianKeplerianElements] = {
    body_id: UranianKeplerianElements(
        name=row.name,
        epoch=_NEELY_EPOCH_JD,
        equinox_jd=_NEELY_EPOCH_JD,
        a=row.a,
        e=row.e,
        i=row.i,
        omega=row.omega,
        Omega=row.node,
        M0=row.mean_anomaly,
        n=_K_GAUSS_DEG / row.a**1.5,
    )
    for body_id, row in _NEELY_SOURCE_ROWS.items()
}

TRANSPLUTO_KEPLERIAN_ELEMENTS = TransplutoKeplerianElements(
    name="Transpluto",
    epoch=_UNAVAILABLE_FLOAT,
    a=_UNAVAILABLE_FLOAT,
    e=_UNAVAILABLE_FLOAT,
    i=_UNAVAILABLE_FLOAT,
    omega=_UNAVAILABLE_FLOAT,
    Omega=_UNAVAILABLE_FLOAT,
    M0=_UNAVAILABLE_FLOAT,
    n=_UNAVAILABLE_FLOAT,
)

# This legacy registry intentionally remains empty.  Complete historical
# predictions live in ``FICTITIOUS_ORBITAL_ELEMENTS`` with explicit equinox
# metadata; Hamburg-school rows live in ``URANIAN_KEPLERIAN_ELEMENTS``; and
# Selena has its own source-defined uniform model.  Adding an unsupported ID
# here would bypass the provenance-aware dispatch and is therefore forbidden.
HYPOTHETICAL_ELEMENTS: Dict[int, HypotheticalElements] = {}

# Lowell's public compatibility container is populated from the same literal
# J1850 source elements used by the runtime registry below.  Pickering's
# container stays importable but unavailable: Circular 215 announces a
# separate memoir without printing the complete elements, and that memoir was
# not recovered during this review (see
# ``docs/methodology/missing-hypothetical-models.md``).
LOWELL_PLANET_X_ELEMENTS = LowellPlanetXElements(
    "Planet X Lowell",
    2396757.5,
    43.0,
    0.202,
    0.0,
    203.8,
    0.0,
    178.3,
    _K_GAUSS_DEG / 43.0**1.5,
)
PICKERING_PLANET_X_ELEMENTS = PickeringPlanetXElements(
    "Planet X Pickering",
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
)

VULCAN_ELEMENTS = VulcanElements(
    "Vulcan",
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
)

WALDEMATH_ELEMENTS = WaldemathElements(
    "Waldemath",
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
    _UNAVAILABLE_FLOAT,
)

HYPOTHETICAL_PROVENANCE: Dict[int, Tuple[str, str]] = {
    CUPIDO: ("primary-transcription", "Neely 1980, Matrix VII, p. 10, Table I"),
    HADES: ("primary-transcription", "Neely 1980, Matrix VII, p. 10, Table I"),
    ZEUS: ("primary-transcription", "Neely 1980, Matrix VII, p. 10, Table I"),
    KRONOS: ("primary-transcription", "Neely 1980, Matrix VII, p. 10, Table I"),
    APOLLON: ("primary-transcription", "Neely 1980, Matrix VII, p. 10, Table I"),
    ADMETOS: ("primary-transcription", "Neely 1980, Matrix VII, p. 10, Table I"),
    VULKANUS: ("primary-transcription", "Neely 1980, Matrix VII, p. 10, Table I"),
    POSEIDON: ("primary-transcription", "Neely 1980, Matrix VII, p. 10, Table I"),
    ISIS: (
        "unsupported",
        "Strubell-attributed primary element publication not recovered; see missing inventory",
    ),
    NIBIRU: (
        "unsupported",
        "no credible source-complete primary orbit identified; see missing inventory",
    ),
    HARRINGTON: ("primary-transcription", "Harrington, AJ 96 (1988), p. 1478"),
    NEPTUNE_LEVERRIER: (
        "primary-transcription",
        "Le Verrier, CRAS 23 (1846), p. 432; Recherches pp. 234-236",
    ),
    NEPTUNE_ADAMS: (
        "primary-transcription",
        "Adams 1846, sections 47-48, printed pp. 25-26",
    ),
    PLUTO_LOWELL: (
        "primary-transcription",
        "Lowell 1915, pp. 5, 9 and 105",
    ),
    PLUTO_PICKERING: (
        "unsupported",
        "Circular 215 omits complete elements; announced memoir not recovered",
    ),
    VULCAN: (
        "unsupported",
        "exact Weston-attributed epoch, frame, and elements not recovered",
    ),
    WHITE_MOON: (
        "published-model",
        "Velichko and Larin 2007, pp. 17, 18, 20, 29, and 45; "
        "source-unwrapped uniform v3 realization",
    ),
    PROSERPINA: (
        "unsupported",
        "complete primary Abramov publication not recovered",
    ),
    WALDEMATH: (
        "unsupported",
        "historical claim omits the complete computational convention",
    ),
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

    The dataset contains the 12 independently transcribed primary-source
    orbital models documented in ``docs/methodology/hypothetical-bodies.md``.
    The six recognised IDs listed in
    ``docs/methodology/missing-hypothetical-models.md`` are intentionally absent;
    Selena is a separately implemented symbolic point, not an orbital CSV row.

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


# Classical predictions with complete, independently reviewed source chains.
#
# These rows intentionally use the epochs printed by the authors.  Keeping the
# source epoch avoids an opaque propagation to a later compatibility epoch and
# makes every constant traceable by elementary arithmetic.  The generic
# propagator below converts from each stated ecliptic/equinox to J2000.
FICTITIOUS_ORBITAL_ELEMENTS = {
    # R. S. Harrington, "The Location of Planet X", AJ 96 (1988), p. 1478:
    # https://articles.adsabs.harvard.edu/pdf/1988AJ.....96.1476H
    #
    # The printed perihelion date 1789-08-06 is JD 2374696.5, hence M=0.
    # Harrington directly prints a, e, omega, node and inclination.  The paper
    # does not specify a modern machine frame, so J2000 is an explicit project
    # convention rather than a claim about wording in the source.
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
    # U.-J. Le Verrier, CRAS 23 (31 August 1846), pp. 428--438, p. 432,
    # together with the more precise table in Recherches, pp. 234--236:
    # https://fr.wikisource.org/wiki/Livre:Comptes_rendus_hebdomadaires_des_s%C3%A9ances_de_l%E2%80%99Acad%C3%A9mie_des_sciences,_tome_023,_1846.djvu
    #
    # Source quantities at 1847-01-01:
    #   a=36.1539 AU, e=0.10761, varpi=284 deg 45 arcmin,
    #   M=34 deg 1 arcmin 56 arcsec.
    # The source explicitly begins with an ecliptic-plane solution.  Therefore
    # node=i=0 are model assumptions printed by the author, not missing values.
    NEPTUNE_LEVERRIER: OrbitalElements(
        name="Leverrier",
        epoch_jd=2395662.5,
        equinox_jd=2395662.5,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(34.0 + 1.0 / 60.0 + 56.0 / 3600.0),
        semi_axis=36.1539,
        eccentricity=TPolynomial(0.10761),
        arg_perihelion=TPolynomial(284.0 + 45.0 / 60.0),
        asc_node=TPolynomial(0.0),
        inclination=TPolynomial(0.0),
    ),
    # J. C. Adams, "An Explanation of the Observed Irregularities in the
    # Motion of Uranus ..." (1846), sections 47--48, printed pp. 25--26:
    # https://archive.org/details/explanationofobs00adam
    #
    # Adams prints mean longitude 323 deg 02 arcmin on 1846-10-06 and
    # longitude of perihelion 299 deg 11 arcmin.  Thus the mean anomaly is the
    # direct difference 23 deg 51 arcmin.  He prints e=0.120615 and his second
    # distance hypothesis a(Uranus)/a(body)=0.515.  Benjamin A. Gould's 1850
    # Smithsonian Report on the History of the Discovery of Neptune, p. 30,
    # explicitly renders that second mean-distance hypothesis as 37.25 AU and
    # cites these same Adams sections.  This near-contemporary independent
    # statement closes the unit conversion without borrowing a later element
    # table.  Section 60 says latitude observations did not determine a
    # satisfactory node/inclination, so the published prediction is planar.
    #
    # The explicit divisions by 60 are intentional audit evidence: the former
    # value 299.11 confused arcminutes with decimal hundredths.  Correcting it
    # is a scientific accuracy fix, not a compatibility fit.
    NEPTUNE_ADAMS: OrbitalElements(
        name="Adams",
        epoch_jd=2395575.5,
        equinox_jd=2395575.5,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(23.0 + 51.0 / 60.0),
        semi_axis=37.25,
        eccentricity=TPolynomial(0.120615),
        arg_perihelion=TPolynomial(299.0 + 11.0 / 60.0),
        asc_node=TPolynomial(0.0),
        inclination=TPolynomial(0.0),
    ),
    # Percival Lowell, Memoir on a Trans-Neptunian Planet (1915), pp. 5, 9,
    # and 105; HathiTrust record 009026348:
    # https://catalog.hathitrust.org/Record/009026348
    #
    # Page 5 defines epsilon as mean longitude at the epoch and varpi as the
    # longitude of perihelion.  Page 9 identifies the epoch as 1850.0.  The
    # preferred solution on p. 105 gives epsilon=22.1, a=43.0, e=0.202 and
    # varpi=203.8, hence M=(22.1-203.8) mod 360=178.3 degrees.  The same page
    # says latitude perturbations did not provide a trustworthy orbital plane;
    # an analogy-based expectation of about 10 degrees is not an element.
    PLUTO_LOWELL: OrbitalElements(
        name="Lowell",
        epoch_jd=2396757.5,
        equinox_jd=2396757.5,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial((22.1 - 203.8) % 360.0),
        semi_axis=43.0,
        eccentricity=TPolynomial(0.202),
        arg_perihelion=TPolynomial(203.8),
        asc_node=TPolynomial(0.0),
        inclination=TPolynomial(0.0),
    ),
}


def calc_fictitious_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """Return a reviewed classical prediction as a heliocentric J2000 state.

    The registry contains Harrington, Le Verrier, Adams and Lowell.  Other
    classical-prediction compatibility IDs raise :class:`UnknownBodyError`.

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
    """Calculate the longitude of a registered Hamburg-school body.

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
    """Calculate the state of a registered Hamburg-school body."""
    if ipl not in URANIAN_KEPLERIAN_ELEMENTS:
        _raise_unsupported_builtin_fictitious(ipl)
    return calc_uranian_planet(ipl, jd_tt)


def calc_cupido(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate Cupido from Neely's reviewed Table-I elements."""
    return calc_uranian_planet(CUPIDO, jd_tt)


def calc_hades(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate Hades from Neely's reviewed Table-I elements."""
    return calc_uranian_planet(HADES, jd_tt)


def calc_zeus(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate Zeus from Neely's reviewed Table-I elements."""
    return calc_uranian_planet(ZEUS, jd_tt)


def calc_kronos(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate Kronos from Neely's reviewed Table-I elements."""
    return calc_uranian_planet(KRONOS, jd_tt)


def calc_apollon(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate Apollon from Neely's reviewed Table-I elements."""
    return calc_uranian_planet(APOLLON, jd_tt)


def calc_admetos(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate Admetos from Neely's reviewed Table-I elements."""
    return calc_uranian_planet(ADMETOS, jd_tt)


def calc_vulkanus(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate Vulkanus from Neely's reviewed Table-I elements."""
    return calc_uranian_planet(VULKANUS, jd_tt)


def calc_poseidon(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate Poseidon from Neely's reviewed Table-I elements."""
    return calc_uranian_planet(POSEIDON, jd_tt)


def calc_transpluto(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Raise until Transpluto has an independently reviewed element set."""
    del jd_tt
    _raise_unsupported_builtin_fictitious(ISIS)


def _calc_transpluto_raw(jd_tt: float) -> Tuple[float, float, float]:
    """Raise until Transpluto has an independently reviewed element set."""
    del jd_tt
    _raise_unsupported_builtin_fictitious(ISIS)


def calc_uranian_planet(
    body_id: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """Calculate a registered Hamburg-school body.

    The independent propagation implementation uses Kepler's equation,
    Gaussian-vector rotation, fixed-J1900-to-J2000 precession, and a centered
    velocity sample.  Every registered row comes from Neely (1980), Table I;
    no empirical correction table is applied.

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
    if body_id not in URANIAN_KEPLERIAN_ELEMENTS:
        _raise_unsupported_builtin_fictitious(body_id)
    elements = URANIAN_KEPLERIAN_ELEMENTS[body_id]
    return _keplerian_to_ecliptic_j2000(elements, jd_tt)


def _raise_unsupported_builtin_fictitious(body_id: int) -> NoReturn:
    """Raise for an ID absent from independently reviewed registries."""
    from .exceptions import UnknownBodyError

    raise UnknownBodyError(
        message=f"no built-in fictitious elements for body {body_id}",
        body_id=body_id,
    )


def calc_vulcan(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Raise until Vulcan has an independently reviewed element set."""
    del jd_tt
    _raise_unsupported_builtin_fictitious(VULCAN)


def _calc_vulcan_raw(jd_tt: float) -> Tuple[float, float, float]:
    """Raise until Vulcan has an independently reviewed element set."""
    del jd_tt
    _raise_unsupported_builtin_fictitious(VULCAN)


def calc_waldemath(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Raise until Waldemath has an independently reviewed element set."""
    del jd_tt
    _raise_unsupported_builtin_fictitious(WALDEMATH)


def calc_transpluto_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Raise until Transpluto has an independently reviewed element set."""
    del jd_tt
    _raise_unsupported_builtin_fictitious(ISIS)


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
    """Return the published uniform seven-year White Moon convention.

    The point is a symbolic geocentric construction, not a physical satellite.
    ``use_true_lilith`` is retained for API compatibility and does not alter
    this independently defined Selena model.
    """
    del use_true_lilith
    elapsed_days = jd_tt - _WHITE_MOON_SOURCE_EPOCH_JD
    longitude = (
        _WHITE_MOON_SOURCE_LONGITUDE_DEG + _WHITE_MOON_RATE_DEG_PER_DAY * elapsed_days
    ) % 360.0
    return (
        float(longitude),
        0.0,
        float(_WHITE_MOON_DISTANCE_AU),
        float(_WHITE_MOON_RATE_DEG_PER_DAY),
        0.0,
        0.0,
    )


def calc_waldemath_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Raise until Waldemath has an independently reviewed element set."""
    return calc_waldemath(jd_tt)


def calc_proserpina(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Raise until Proserpina has an independently reviewed element set."""
    del jd_tt
    _raise_unsupported_builtin_fictitious(PROSERPINA)


def calc_planet_x_lowell(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate Lowell's preferred 1915 Planet-X solution."""
    return calc_fictitious_position(PLUTO_LOWELL, jd_tt)


def calc_planet_x_pickering(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Raise until Pickering's orbit is independently transcribed."""
    del jd_tt
    _raise_unsupported_builtin_fictitious(PLUTO_PICKERING)


def calc_hypothetical_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate the position of a supported hypothetical body.

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
        UnknownBodyError: If the recognised ID has no independently reviewed
            built-in element set.  Supported IDs are 40--47, 50--53 and 56.
    """
    if not is_hypothetical_body(ipl):
        raise ValueError(f"Body ID {ipl} is not a hypothetical body")

    if ipl in URANIAN_KEPLERIAN_ELEMENTS:
        return calc_uranian_planet(ipl, jd_tt)

    if ipl in FICTITIOUS_ORBITAL_ELEMENTS:
        return calc_fictitious_position(ipl, jd_tt)

    if ipl == WHITE_MOON:
        return calc_white_moon_position(jd_tt)

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
