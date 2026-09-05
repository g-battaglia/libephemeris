# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Independently sourced hypothetical-body calculations.

The historical bodies exposed as IDs 40--58 are mathematical conventions, not
claims that the objects exist.  A built-in calculation is enabled only when
its numerical model can be transcribed or reconstructed from identified
public publications.  IDs whose model could not be recovered remain
importable and discoverable for API compatibility, but fail closed with
:class:`~libephemeris.exceptions.UnknownBodyError`.

The supported set is deliberately explicit:

* IDs 40--47 use James Neely's published J1900 Hamburg-school elements;
* ID 48 uses the published Landscheidt-tradition Transpluto orbit;
* ID 50 uses Harrington's 1988 nominal Planet-X orbit;
* IDs 51--52 use the primary 1846 predictions of Le Verrier and Adams;
* ID 53 uses Lowell's preferred 1915 Planet-X solution;
* ID 54 uses Pickering's published Planet-O prediction;
* ID 55 uses the intramercurial Vulcan convention (compatibility values);
* ID 56 is the published uniform seven-year Selena convention; and
* ID 57 uses the circular Proserpina convention (compatibility values).

* ID 58 (Waldemath) follows Sepharial's published dark-moon convention
  (Sepharial, "The Science of Foreknowledge", 1918: 177-day synodic figure
  and the predicted 1898-02-02 transit anchor).

ID 49 (Nibiru) stays fail-closed: no publicly published element set has
been recovered for it.

The source transcription, frame choices and arithmetic transformations are
documented beside the constants below, in ``fictitious_orbits.csv`` and in
``docs/methodology/hypothetical-bodies.md``.  Missing fields for disabled IDs
are inventoried in ``docs/methodology/missing-hypothetical-models.md``.

Provenance:
    Each enabled ID's record names its public source and the project
    conversions applied to it.  Runtime propagation uses only those
    transcriptions plus published two-body geometry.  The remaining
    source-incomplete IDs have no hidden replacement constants and raise
    ``UnknownBodyError``.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import (
    Callable,
    Dict,
    List,
    NoReturn,
    Optional,
    Tuple,
    Type,
    TypeVar,
    Union,
)

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
# SHARED CONSTANTS
# =============================================================================
# Epoch and equinox J1900.0 = JD 2415020.0 TT.  The Hamburg-school, Transpluto,
# Vulcan and Proserpina element sets are all published for this epoch.
_J1900_EPOCH_JD: float = 2415020.0

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


# =============================================================================
# WALDEMATH / SEPHARIAL DARK MOON (ID 58) — PUBLISHED UNIFORM MODEL
# =============================================================================
# The astrological "Waldemath Dark Moon" is, in the primary literature, a
# mean-longitude-only construct, not a Keplerian orbit:
#
# * Sepharial [Walter Gorn Old], "The Science of Foreknowledge", W. Foulsham,
#   London, 1918, chapter "The New Satellite — Lilith" (pp. 39-45 of the
#   public-domain scan), prints the operative rule: "Multiply the number of
#   days from the day of last conjunction to the day of birth by three, and
#   add that number of degrees to the longitude of the last conjunction" —
#   i.e. uniform 3 deg/day geocentric tropical motion — together with the
#   177-day synodic period ("synodical revolution at 177 days") and a dated
#   table of Lilith-Sun conjunctions 1854-1906.
# * Georg Waltemath's 1898 Hamburg announcements (reported in Science, New
#   Series, Vol. 8, No. 189, 1898-08-12, p. 185, and summarized by Joseph
#   Ashbrook, "The Many Moons of Dr. Waltemath", Sky & Telescope, Vol. 28,
#   p. 218, October 1964) give the physical scale: mean distance about
#   1.03 million km, sidereal period 119 days, synodic 177 days, and the
#   predicted transit of the Sun on 1898 February 2-4.
#
# Realization: the anchor is Waltemath's predicted 1898 February 2 transit
# (a Lilith-Sun conjunction in Sepharial's table), taken at 00:00 GMT — the
# midnight-GMT convention of the Delphine Jay AFA ephemerides (Jay, "Lilith
# Ephemeris 1900-2000 A.D.", AFA, 1983). At a conjunction the Lilith
# longitude equals the Sun's apparent geocentric longitude; at that instant
# (JD 2414322.5 UT) DE440 gives 313.2150483608214 deg (reproducible with
# this library's own Sun).
#
# The mean motion comes from Sepharial's printed synodic period: the point
# returns to conjunction with the Sun every 177 days, so
#   n = 360/177 + 360/365.2422 deg/day = 3.0195456... deg/day,
# where 365.2422 days is the mean tropical year (Meeus, Astronomical
# Algorithms, 2nd ed., 1998). Sepharial's own consistency statements pin
# this reading of the primary sources: his "the satellite Lilith returns
# to the same longitude on the same day in 126 years" is exactly
# 126 x 365.2422 / 177 = 260.00 synodic revolutions, and the implied
# sidereal period 360/n = 119.2 days matches Waltemath's published
# 119-day sidereal period. (Sepharial's hand-computation rule rounds the
# rate to "3 deg per day"; used globally that rounding cannot reproduce
# his own 1854-1906 conjunction table, so the printed synodic period is
# the defining quantity and the 3 deg/day wording is its practical
# approximation.)
#
# The model is longitude-only: uniform motion on the ecliptic (zero
# latitude), circular by construction, with the published mean distance as
# the fixed radial coordinate. Eccentricity, inclination and node values
# that circulate in legacy tabulations (e = 0.1587, i = 2.5) have no
# traceable primary source — the 2.5 figure is most plausibly a conflation
# with Waltemath's ~2.5 arcminute apparent DIAMETER — and are deliberately
# not used.
_WALDEMATH_ANCHOR_JD_UT: float = 2414322.5
_WALDEMATH_ANCHOR_APPARENT_LON_DEG: float = 313.2150483608214
_WALDEMATH_SYNODIC_DAYS: float = 177.0
_WALDEMATH_TROPICAL_YEAR_DAYS: float = 365.2422
_WALDEMATH_RATE_DEG_PER_DAY: float = (
    360.0 / _WALDEMATH_SYNODIC_DAYS + 360.0 / _WALDEMATH_TROPICAL_YEAR_DAYS
)
# Waltemath's published mean distance, converted with the exact astronomical
# unit adopted by IAU 2012 Resolution B2 (149 597 870 700 m).
_WALDEMATH_DISTANCE_AU: float = 1.03e6 / 149_597_870.7


# =============================================================================
# THE HYPOTHETICAL-BODY TABLE
# =============================================================================
# Every identifier from 40 to 58 has exactly one row.  A row carries the body's
# identity, the provenance of its numbers and the elements the runtime
# propagates.  ``model`` names which propagator reads the row, because these
# historical constructions are not all Keplerian orbits and the two families
# that are do not share a velocity stencil.  The public element containers and
# the ``calc_*`` entry points further down are views on this table.

#: Neely's Hamburg-school orbits: Keplerian propagation from the fixed J1900
#: ecliptic, velocity from a one-day centred difference.
MODEL_HAMBURG = "hamburg"
#: The classical trans-Neptunian predictions (Harrington, Le Verrier, Adams,
#: Lowell, Pickering): the same propagation from each author's own equinox,
#: with a half-day centred difference.
MODEL_CLASSICAL = "classical"
#: Transpluto: the orbit lies in the J1900 ecliptic, so the longitude is the
#: perihelion longitude plus the true anomaly and the latitude is exactly zero.
MODEL_TRANSPLUTO = "transpluto"
#: Vulcan: Keplerian, but the mean anomaly, the perihelion and the node move
#: linearly in Julian centuries and the angles are referred to the equinox of
#: date.
MODEL_VULCAN = "vulcan"
#: Selena: uniform tropical longitude at a fixed radius, not an orbit.
MODEL_SELENA = "selena"
#: The Waldemath Dark Moon: uniform tropical longitude like Selena, but its
#: anchor is an apparent longitude, so the propagator removes the nutation of
#: the anchor epoch before advancing.
MODEL_DARK_MOON = "dark-moon"
#: Proserpina: a uniform circular heliocentric orbit, so the longitude is the
#: mean longitude and both the latitude and the radial speed are zero.
MODEL_CIRCULAR = "circular"
#: An identifier the public API recognises but cannot compute: no primary
#: element set has been recovered for it, so it fails closed.
MODEL_UNSUPPORTED = "unsupported"

# Provenance categories: a literal transcription of one identified primary
# publication, or a convention realized from documented public sources.  A row
# with neither is categorised by the same word as its model.
_PRIMARY_TRANSCRIPTION = "primary-transcription"
_PUBLISHED_MODEL = "published-model"
_NO_PRIMARY_SOURCE = MODEL_UNSUPPORTED

_NEELY_SOURCE = "Neely 1980, Matrix VII, p. 8, Table I"


@dataclass(frozen=True)
class HypotheticalBody:
    """One row of the hypothetical-body table.

    Angles are in degrees, distances in astronomical units and rates per day
    unless the field name says otherwise.  A field the row's model does not use
    keeps its zero default.

    Attributes:
        body_id: Public identifier, 40 to 58.
        name: Name the public API answers with.
        model: Which propagator reads this row; one of the ``MODEL_*``
            constants above.
        category: Provenance class: ``primary-transcription`` for a literal
            transcription of one identified publication, ``published-model``
            for a convention realized from documented public sources, and
            ``unsupported`` when no primary element set was recovered.
        source: The publication behind the numbers, or the note recording that
            the primary source is not verified.
        epoch: Julian Day TT the elements belong to.  The Dark Moon anchor is
            the one exception and is a Julian Day UT, as its source prints it.
        equinox_jd: Julian Day TT of the ecliptic and equinox the angles are
            referred to; 0.0 means the equinox of date, which has no fixed
            frame to precess from.
        a: Semi-major axis, or the fixed radius of the longitude-only models.
        e: Eccentricity.
        i: Inclination.
        omega: Argument of perihelion at ``epoch``.
        Omega: Longitude of the ascending node at ``epoch``.
        M0: Mean anomaly at ``epoch``; for the circular models it is the mean
            longitude, and for the Dark Moon the anchor's apparent longitude.
        n: Mean motion in degrees per day.
        n_century: Mean motion in degrees per Julian century, read instead of
            ``n`` by the Vulcan model.
        omega_rate: Secular motion of ``omega``, degrees per Julian century.
        Omega_rate: Secular motion of ``Omega``, degrees per Julian century.
        printed_rate_century: The mean-anomaly rate the author printed, in
            degrees per Julian century.  Only Neely's rows have one; the
            runtime derives ``n`` from ``a`` instead, and the provenance gate
            checks that the two agree to the table's print precision.
    """

    body_id: int
    name: str
    model: str
    category: str
    source: str
    epoch: float = 0.0
    equinox_jd: float = 0.0
    a: float = 0.0
    e: float = 0.0
    i: float = 0.0
    omega: float = 0.0
    Omega: float = 0.0
    M0: float = 0.0
    n: float = 0.0
    n_century: float = 0.0
    omega_rate: float = 0.0
    Omega_rate: float = 0.0
    printed_rate_century: float = 0.0

    @property
    def mean_longitude(self) -> float:
        """Return M0 + omega + Omega at ``epoch``, in degrees."""
        return (self.M0 + self.omega + self.Omega) % 360.0


def _gaussian_row(
    body_id: int,
    name: str,
    *,
    model: str,
    category: str,
    source: str,
    epoch: float,
    equinox_jd: float,
    a: float,
    e: float = 0.0,
    i: float = 0.0,
    omega: float = 0.0,
    Omega: float = 0.0,
    M0: float = 0.0,
) -> HypotheticalBody:
    """Build a row whose mean motion is the two-body Gaussian value.

    Every orbit built here is a small body around the Sun, so the semimajor
    axis fixes the mean motion, n = k / a^(3/2).  Deriving it gives one
    internally consistent orbit per row and keeps the semimajor axis written
    down exactly once.
    """
    return HypotheticalBody(
        body_id=body_id,
        name=name,
        model=model,
        category=category,
        source=source,
        epoch=epoch,
        equinox_jd=equinox_jd,
        a=a,
        e=e,
        i=i,
        omega=omega,
        Omega=Omega,
        M0=M0,
        n=_K_GAUSS_DEG / a**1.5,
    )


def _neely_row(
    body_id: int,
    name: str,
    a: float,
    e: float,
    i: float,
    argp: float,
    node: float,
    mean_anomaly: float,
    printed_rate_century: float,
) -> HypotheticalBody:
    """Build the runtime row for one printed Neely Table-I orbit.

    The arguments after ``body_id`` are the eight published columns in the
    order Table I prints them.  Two documented transformations turn them into
    a propagable row.

    Degenerate-angle normalization: for the four circular coplanar rows
    (e=i=node=0: Apollon, Admetos, Vulkanus, Poseidon) the perihelion
    direction is geometrically undefined and only the sum M+omega+node is
    observable.  The row stores that observable phase in ``M0`` with
    omega=Omega=0, so ``M0`` *is* the mean longitude at the epoch, matching the
    historical public contract of this table.  For e=0 the propagation depends
    on the angles only through their sum, so this representation is
    mathematically identical to the literal column placement kept in the
    reviewed CSV (Kepler's equation is the identity E=M at e=0).

    Mean motion: derived from the transcribed semimajor axis rather than copied
    from the rounded rate column.  The two differ by less than
    0.0027 deg/century, the expected effect of the table's rounded ``a`` and
    rate; the printed rate stays in the row for the audit checks and for the
    legacy mean-longitude view.
    """
    degenerate = e == 0.0 and i == 0.0
    return HypotheticalBody(
        body_id=body_id,
        name=name,
        model=MODEL_HAMBURG,
        category=_PRIMARY_TRANSCRIPTION,
        source=_NEELY_SOURCE,
        epoch=_J1900_EPOCH_JD,
        equinox_jd=_J1900_EPOCH_JD,
        a=a,
        e=e,
        i=i,
        omega=0.0 if degenerate else argp,
        Omega=0.0 if degenerate else node,
        M0=(mean_anomaly + argp + node) % 360.0 if degenerate else mean_anomaly,
        n=_K_GAUSS_DEG / a**1.5,
        printed_rate_century=printed_rate_century,
    )


# Primary source and transcription policy for IDs 40--47
# ------------------------------------------------------
# James Neely, "Orbital Elements for the Transneptunian Planets", Matrix
# Magazine, issue VII (1980), pp. 6--12, Table I on p. 8.  A public scan is
# catalogued as item 18 by Alexandria iBase:
#
#   https://special-collections.alexandriaibase.org/items/show/18
#
# SHA-256 of the reviewed 63-page scan:
#   f3d2e834a4a684736383ec387c9c2786064ca4668256724865da2f3d8fa03020
#
# Neely defines T as Julian centuries from JD 2415020.0.  Table I prints the
# columns of ``_NEELY_TABLE_I`` below plus a common node term +1.3960*T.  That term expresses the moving equinox used by the source.
# LibEphemeris instead keeps the T=0 angles in the fixed J1900 ecliptic frame
# and applies its documented IAU/Vondrak precession exactly once in
# ``_keplerian_to_ecliptic_j2000``.  Adding Neely's moving-equinox term as well
# would double-count precession.
#
# Apollon, Admetos, Vulkanus and Poseidon are circular and coplanar in Table I.
# Neely's literal convention for them is M=0 with the phase in omega; that
# placement is preserved in the reviewed ``fictitious_orbits.csv``, while the
# rows below carry the equivalent propagable form ``_neely_row`` describes.

#: Table I in identifier order, IDs 40 (Cupido) to 47 (Poseidon), in the
#: source's own column order: name, semimajor axis, eccentricity, inclination,
#: argument of perihelion, node, mean anomaly, printed mean-anomaly rate.
_NEELY_TABLE_I: Tuple[
    Tuple[str, float, float, float, float, float, float, float], ...
] = (
    ("Cupido", 40.99837, 0.00460, 1.0833, 171.4333, 129.8325, 163.7409, 137.13640),
    ("Hades", 50.66744, 0.00245, 1.0500, 148.1796, 161.3339, 27.6496, 99.81803),
    ("Zeus", 59.21436, 0.00120, 0.0, 299.0440, 0.0, 165.1232, 79.00633),
    ("Kronos", 64.81690, 0.00305, 0.0, 208.8801, 0.0, 169.0193, 68.98746),
    ("Apollon", 70.29949, 0.0, 0.0, 138.0533, 0.0, 0.0, 61.0765),
    ("Admetos", 73.62765, 0.0, 0.0, 351.3350, 0.0, 0.0, 56.98245),
    ("Vulkanus", 77.25568, 0.0, 0.0, 55.8983, 0.0, 0.0, 53.015987),
    ("Poseidon", 83.66907, 0.0, 0.0, 165.5163, 0.0, 0.0, 47.03868),
)

_TABLE_ROWS: Tuple[HypotheticalBody, ...] = (
    *(
        _neely_row(body_id, *columns)
        for body_id, columns in zip(
            range(CUPIDO, POSEIDON + 1), _NEELY_TABLE_I, strict=True
        )
    ),
    # Transpluto (Isis/Bacchus), ID 48.
    #
    # The astrological Transpluto orbit family originates with Francis M. E.
    # Sevin (1946), was realized as the Landscheidt 1972 graphical ephemeris
    # (Kosmobiologische Akademie Aalen), and is printed as a complete element
    # set in John Robert Hawkins, "Transpluto, or Should We Call Him Bacchus,
    # the Ruler of Taurus?" (Hawkins Enterprising Publications, Dallas, 1976;
    # 2nd ed. 1978, p. 79; AFA reprint ISBN 0-86690-386-0).  Hawkins prints,
    # for the J1900.0 epoch and equinox:
    #
    #   a = 77.755 AU, e = 0.3, i = 0, node = 0,
    #   longitude of perihelion = 0.0438748 deg,
    #   mean anomaly = 66.806096 deg,
    #   perihelion passage 1772.76, sidereal period 685.65 Julian years.
    #
    # Internal consistency: (1900.0 - 1772.76) * (360/685.65) = 66.807 deg,
    # reproducing the printed mean anomaly to ~0.001 deg.  The mean motion is
    # the printed sidereal period rather than the Gaussian value, so this row
    # is written out in full.  With i = 0 the orbit lies in the J1900 ecliptic
    # plane and the runtime reports the point on the ecliptic (latitude exactly
    # zero, the historical output convention for this body).
    #
    # A variant semimajor axis 77.775 AU circulates in software-derived
    # material; no printed publication for it was recovered, so the book value
    # 77.755 is used.
    HypotheticalBody(
        body_id=ISIS,
        name="Transpluto",
        model=MODEL_TRANSPLUTO,
        category=_PUBLISHED_MODEL,
        source=(
            "Sevin 1946 / Landscheidt 1972 lineage as printed in Hawkins 1978, "
            "p. 79: a=77.755, e=0.3, planar, M(J1900)=66.806096, P=685.65 yr"
        ),
        epoch=_J1900_EPOCH_JD,
        equinox_jd=_J1900_EPOCH_JD,
        a=77.755,
        e=0.3,
        omega=0.0438748,
        M0=66.806096,
        n=360.0 / (685.65 * 365.25),
    ),
    # Nibiru, ID 49: recognised by name and number, computed by nothing.  See
    # docs/methodology/missing-hypothetical-models.md.
    HypotheticalBody(
        body_id=NIBIRU,
        name="Nibiru",
        model=MODEL_UNSUPPORTED,
        category=_NO_PRIMARY_SOURCE,
        source=(
            "no credible source-complete primary orbit identified; "
            "see missing inventory"
        ),
    ),
    # The classical predictions, IDs 50--54, intentionally use the epochs
    # printed by their authors.  Keeping the source epoch avoids an opaque
    # propagation to a later compatibility epoch and makes every constant
    # traceable by elementary arithmetic; the propagator converts from each
    # stated ecliptic and equinox to J2000.
    #
    # R. S. Harrington, "The Location of Planet X", AJ 96 (1988), p. 1478:
    # https://articles.adsabs.harvard.edu/pdf/1988AJ.....96.1476H
    #
    # The printed perihelion date 1789-08-06 is JD 2374696.5, hence M=0.
    # Harrington directly prints a, e, omega, node and inclination.  The paper
    # does not state the equinox of its angles; referring them to the J2000
    # equinox is a convention of this library, not a statement about the
    # source.
    _gaussian_row(
        HARRINGTON,
        "Harrington",
        model=MODEL_CLASSICAL,
        category=_PRIMARY_TRANSCRIPTION,
        source="Harrington, AJ 96 (1988), p. 1478",
        epoch=2374696.5,
        equinox_jd=2451545.0,
        a=101.2,
        e=0.411,
        i=32.4,
        omega=208.5,
        Omega=275.4,
        M0=0.0,
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
    _gaussian_row(
        NEPTUNE_LEVERRIER,
        "Neptune-Leverrier",
        model=MODEL_CLASSICAL,
        category=_PRIMARY_TRANSCRIPTION,
        source="Le Verrier, CRAS 23 (1846), p. 432; Recherches pp. 234-236",
        epoch=2395662.5,
        equinox_jd=2395662.5,
        a=36.1539,
        e=0.10761,
        i=0.0,
        omega=284.0 + 45.0 / 60.0,
        Omega=0.0,
        M0=34.0 + 1.0 / 60.0 + 56.0 / 3600.0,
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
    _gaussian_row(
        NEPTUNE_ADAMS,
        "Neptune-Adams",
        model=MODEL_CLASSICAL,
        category=_PRIMARY_TRANSCRIPTION,
        source="Adams 1846, sections 47-48, printed pp. 25-26",
        epoch=2395575.5,
        equinox_jd=2395575.5,
        a=37.25,
        e=0.120615,
        i=0.0,
        omega=299.0 + 11.0 / 60.0,
        Omega=0.0,
        M0=23.0 + 51.0 / 60.0,
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
    # an analogy-based expectation of about 10 degrees is not an element, so
    # this propagable row stays planar.
    _gaussian_row(
        PLUTO_LOWELL,
        "Planet X Lowell",
        model=MODEL_CLASSICAL,
        category=_PRIMARY_TRANSCRIPTION,
        source="Lowell 1915, pp. 5, 9 and 105",
        epoch=2396757.5,
        equinox_jd=2396757.5,
        a=43.0,
        e=0.202,
        i=0.0,
        omega=203.8,
        Omega=0.0,
        M0=(22.1 - 203.8) % 360.0,
    ),
    # W. H. Pickering, "The Transneptunian Planet", Annals of Harvard College
    # Observatory 82, No. 3 (1919), pp. 49-59; complete element list on p. 59:
    # https://articles.adsabs.harvard.edu/pdf/1937AnHar..82...49P
    #
    # Pickering prints: epoch 1920.0, mean distance 55.1, period 409, mean
    # annual motion 0.880 deg, longitude of the perihelion 280 deg, date of
    # the perihelion 1720.0, e = 0.31, perihelion distance 38.0, aphelion
    # distance 72.2, and closes with "The longitude of the Ascending Node and
    # the Inclination may perhaps be estimated at 100 +/- 5 deg and
    # 15 +/- 5 deg, respectively."  The printed perihelion date is the
    # element anchor: M = 0 at Julian epoch 1720.0, i.e.
    # JD 2451545 - 280*365.25 = 2349275.0.  The argument of perihelion is the
    # direct difference 280 - 100 = 180 deg.  The paper's angles belong to its
    # own epoch's equinox; Julian epoch 1920.0 = JD 2422325.0 is used as the
    # frame, converted to J2000 by the shared propagator.
    _gaussian_row(
        PLUTO_PICKERING,
        "Planet X Pickering",
        model=MODEL_CLASSICAL,
        category=_PRIMARY_TRANSCRIPTION,
        source="Pickering 1919, Annals of Harvard College Observatory 82(3), p. 59",
        epoch=2349275.0,
        equinox_jd=2422325.0,
        a=55.1,
        e=0.31,
        i=15.0,
        omega=280.0 - 100.0,
        Omega=100.0,
        M0=0.0,
    ),
    # Vulcan, ID 55: the intramercurial Vulcan of the astrological tradition,
    # as a J1900 linear parameterization.  The mean anomaly, the argument of
    # perihelion and the ascending node each move linearly in Julian centuries
    # TT from JD 2415020, and the angles are referred to the ecliptic and
    # equinox of date.  The perihelion and node rates are equal and opposite,
    # so the longitude of perihelion omega + Omega stays fixed at
    # 370 deg = 10 deg.  All of these numbers are compatibility values,
    # primary source not verified.
    HypotheticalBody(
        body_id=VULCAN,
        name="Vulcan",
        model=MODEL_VULCAN,
        category=_PUBLISHED_MODEL,
        source=(
            "intramercurial Vulcan convention, J1900 linear parameterization: "
            "compatibility values, primary source not verified"
        ),
        epoch=_J1900_EPOCH_JD,
        a=0.13744,
        e=0.019,
        i=7.5,
        omega=322.212069,
        Omega=47.787931,
        M0=252.8987988,
        n_century=707550.7341,
        omega_rate=1670.056,
        Omega_rate=-1670.056,
    ),
    # White Moon / Selena, ID 56: the uniform seven-year convention derived in
    # the block above.  ``a`` is the declared display radius, not a claim that
    # Selena is a material satellite.
    HypotheticalBody(
        body_id=WHITE_MOON,
        name="White Moon",
        model=MODEL_SELENA,
        category=_PUBLISHED_MODEL,
        source=(
            "Globa-school uniform seven-year Selena convention (Kutalev "
            "encyclopedia; Velichko and Larin 2007); fixed compatibility "
            "realization inside the published arcminute-level scatter"
        ),
        epoch=_WHITE_MOON_SOURCE_EPOCH_JD,
        a=_WHITE_MOON_DISTANCE_AU,
        M0=_WHITE_MOON_SOURCE_LONGITUDE_DEG,
        n=_WHITE_MOON_RATE_DEG_PER_DAY,
    ),
    # Proserpina, ID 57: the circular trans-Plutonian convention of the
    # astrological tradition, referred to the ecliptic and mean equinox of
    # date.  With e = i = 0 the perihelion and the node are degenerate, so the
    # observable phase is the mean longitude carried in ``M0`` with
    # omega = Omega = 0, a convention of this library.  Compatibility values,
    # primary source not verified.
    _gaussian_row(
        PROSERPINA,
        "Proserpina",
        model=MODEL_CIRCULAR,
        category=_PUBLISHED_MODEL,
        source=(
            "circular trans-Plutonian convention: a=79.22563 AU, mean longitude "
            "170.73 deg at J1900, Gaussian mean motion; compatibility values, "
            "primary source not verified"
        ),
        epoch=_J1900_EPOCH_JD,
        equinox_jd=0.0,
        a=79.22563,
        M0=170.73,
    ),
    # Waldemath Dark Moon, ID 58: the uniform model derived in the block above.
    # ``epoch`` is the anchor instant in UT, as the source prints it, and
    # ``M0`` is the Sun's apparent longitude there; the propagator removes the
    # anchor's nutation before advancing at ``n``.
    HypotheticalBody(
        body_id=WALDEMATH,
        name="Waldemath",
        model=MODEL_DARK_MOON,
        category=_PUBLISHED_MODEL,
        source=(
            "Sepharial 1918 (Science of Foreknowledge, ch. 'The New Satellite — "
            "Lilith'): uniform 3 deg/day tropical longitude anchored at the "
            "printed 1898 Feb 2 Lilith-Sun conjunction; distance 1.03e6 km from "
            "Waltemath 1898 (Science 8/189 p. 185; Ashbrook, S&T 28, p. 218)"
        ),
        epoch=_WALDEMATH_ANCHOR_JD_UT,
        a=_WALDEMATH_DISTANCE_AU,
        M0=_WALDEMATH_ANCHOR_APPARENT_LON_DEG,
        n=_WALDEMATH_RATE_DEG_PER_DAY,
    ),
)

#: Every recognised hypothetical body, keyed by public identifier.
HYPOTHETICAL_BODIES: Dict[int, HypotheticalBody] = {
    body.body_id: body for body in _TABLE_ROWS
}

HYPOTHETICAL_PROVENANCE: Dict[int, Tuple[str, str]] = {
    body_id: (body.category, body.source)
    for body_id, body in HYPOTHETICAL_BODIES.items()
}

HYPOTHETICAL_NAMES: Dict[int, str] = {
    body_id: body.name for body_id, body in HYPOTHETICAL_BODIES.items()
}


# =============================================================================
# LEGACY PUBLIC ELEMENT CONTAINERS
# =============================================================================
# Read-only views on the table above, kept because they are part of the public
# surface.  Every number in them is a table field except the four marked at
# their construction: quantities Lowell and Pickering published which the
# runtime deliberately does not propagate.


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
class HadesKeplerianElements(HypotheticalElements):
    """Legacy rc7 public element type for Hades."""


@dataclass
class LowellPlanetXElements(HypotheticalElements):
    """Legacy public element container for Lowell's Planet X."""


@dataclass
class PickeringPlanetXElements(HypotheticalElements):
    """Legacy public element container for Pickering's Planet O/X."""


@dataclass
class _CircularElements:
    """Element type of the legacy circular Hamburg containers.

    The same nine fields as :class:`HypotheticalElements`, except that the
    epoch phase is a mean longitude ``L0`` rather than a mean anomaly, which is
    what a circular orbit actually determines.
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


@dataclass
class CupidoKeplerianElements(_CircularElements):
    """Legacy rc7 public element type for Cupido."""


@dataclass
class ZeusKeplerianElements(_CircularElements):
    """Legacy rc7 public element type for Zeus."""


@dataclass
class KronosKeplerianElements(_CircularElements):
    """Legacy rc7 public element type for Kronos."""


@dataclass
class ApollonKeplerianElements(_CircularElements):
    """Legacy rc7 public element type for Apollon."""


@dataclass
class AdmetosKeplerianElements(_CircularElements):
    """Legacy rc7 public element type for Admetos."""


@dataclass
class VulkanusKeplerianElements(_CircularElements):
    """Legacy rc7 public element type for Vulkanus."""


@dataclass
class PoseidonKeplerianElements(_CircularElements):
    """Legacy rc7 public element type for Poseidon."""


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
    """Public element container for the historical intramercurial Vulcan (55).

    The fields describe the astrological Vulcan convention as this library
    realizes it: a fast intramercurial Keplerian orbit referred to the
    ecliptic and equinox of date, with linear secular motion of the mean
    anomaly, argument of perihelion and ascending node in Julian centuries
    from the J1900 epoch.  The perihelion and node rates are equal and
    opposite, so the longitude of perihelion ``omega + Omega`` is constant.
    The numbers are compatibility values, primary source not verified.
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
    """Legacy public container for the Waldemath Dark Moon (ID 58).

    The exported instance describes the published uniform model (Sepharial
    1918: 3 deg/day geocentric tropical longitude anchored at a printed
    Lilith-Sun conjunction; Waltemath 1898: mean distance ~1.03 million km)
    — see the WALDEMATH model block for the full citations. ``a`` carries
    the published mean distance in AU, ``L0``/``M0`` the anchor longitude,
    ``n`` the 3 deg/day rate; eccentricity, inclination and the orientation
    angles are structural zeros of the circular longitude-only model (the
    e/i values in legacy tabulations have no traceable primary source).
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


_KeplerianView = TypeVar("_KeplerianView", bound=HypotheticalElements)
_CircularView = TypeVar("_CircularView", bound=_CircularElements)


def _printed_rate_per_day(body_id: int) -> float:
    """Return Neely's printed mean-anomaly rate in degrees per day."""
    return HYPOTHETICAL_BODIES[body_id].printed_rate_century / 36525.0


def _keplerian_view(
    container: Type[_KeplerianView],
    body_id: int,
    *,
    i: Optional[float] = None,
    M0: Optional[float] = None,
    n: Optional[float] = None,
) -> _KeplerianView:
    """Project one table row onto a legacy nine-field element container.

    ``i``, ``M0`` and ``n`` may be overridden where a legacy container reports
    a published quantity the runtime row deliberately does not propagate;
    every other field comes from the row.
    """
    body = HYPOTHETICAL_BODIES[body_id]
    return container(
        body.name,
        body.epoch,
        body.a,
        body.e,
        body.i if i is None else i,
        body.omega,
        body.Omega,
        body.M0 if M0 is None else M0,
        body.n if n is None else n,
    )


def _circular_view(container: Type[_CircularView], body_id: int) -> _CircularView:
    """Project a Hamburg row onto one of the legacy circular containers.

    The Hamburg school's classical description is a uniform circle in the
    ecliptic: a radius, a mean longitude at epoch and the author's printed rate
    (see, e.g., Neely 1980, Matrix VII, which presents its Table I as a
    refinement of that earlier circular convention).  The eccentric and
    inclined parts of the reviewed row are therefore dropped here and the epoch
    phase is the row's mean longitude.  These containers are display tables:
    the runtime propagation is unaffected by them.
    """
    body = HYPOTHETICAL_BODIES[body_id]
    return container(
        name=body.name,
        epoch=body.epoch,
        a=body.a,
        e=0.0,
        i=0.0,
        omega=0.0,
        Omega=0.0,
        L0=body.mean_longitude,
        n=_printed_rate_per_day(body_id),
    )


CUPIDO_KEPLERIAN_ELEMENTS = _circular_view(CupidoKeplerianElements, CUPIDO)
ZEUS_KEPLERIAN_ELEMENTS = _circular_view(ZeusKeplerianElements, ZEUS)
KRONOS_KEPLERIAN_ELEMENTS = _circular_view(KronosKeplerianElements, KRONOS)
APOLLON_KEPLERIAN_ELEMENTS = _circular_view(ApollonKeplerianElements, APOLLON)
ADMETOS_KEPLERIAN_ELEMENTS = _circular_view(AdmetosKeplerianElements, ADMETOS)
VULKANUS_KEPLERIAN_ELEMENTS = _circular_view(VulkanusKeplerianElements, VULKANUS)
POSEIDON_KEPLERIAN_ELEMENTS = _circular_view(PoseidonKeplerianElements, POSEIDON)

# Hades is the one Hamburg body whose legacy container always carried Neely's
# full eccentric, inclined row, so it keeps doing so; only its rate is the
# printed one rather than the Gaussian value the runtime uses.
HADES_KEPLERIAN_ELEMENTS = _keplerian_view(
    HadesKeplerianElements, HADES, n=_printed_rate_per_day(HADES)
)

URANIAN_ELEMENTS: Dict[int, UranianElements] = {
    body_id: UranianElements(
        name=body.name,
        L0=body.mean_longitude,
        n=_printed_rate_per_day(body_id),
        amplitude=0.0,
        phase=0.0,
        phase_rate=0.0,
    )
    for body_id, body in HYPOTHETICAL_BODIES.items()
    if body.model == MODEL_HAMBURG
}

URANIAN_KEPLERIAN_ELEMENTS: Dict[int, UranianKeplerianElements] = {
    body_id: UranianKeplerianElements(
        name=body.name,
        epoch=body.epoch,
        equinox_jd=body.equinox_jd,
        a=body.a,
        e=body.e,
        i=body.i,
        omega=body.omega,
        Omega=body.Omega,
        M0=body.M0,
        n=body.n,
    )
    for body_id, body in HYPOTHETICAL_BODIES.items()
    if body.model == MODEL_HAMBURG
}

TRANSPLUTO_KEPLERIAN_ELEMENTS = _keplerian_view(TransplutoKeplerianElements, ISIS)

HYPOTHETICAL_ELEMENTS: Dict[int, HypotheticalElements] = {
    # Transpluto's exported rc7 container subclasses HypotheticalElements, so
    # the same object serves both public views.
    ISIS: TRANSPLUTO_KEPLERIAN_ELEMENTS,
    PROSERPINA: _keplerian_view(HypotheticalElements, PROSERPINA),
}

# Lowell's public container reports the memoir's stated inclination
# expectation of about 10 degrees (Lowell 1915; the latitude discussion
# accompanying the preferred p. 105 solution), because this table documents
# what Lowell published rather than the reduced planar propagation model.  The
# reproducible runtime row keeps node=i=0: Lowell did not determine a node, so
# the 10-degree figure cannot orient an orbital plane.  Its mean anomaly is the
# printed 178.3 degrees rather than the subtraction the runtime evaluates; the
# two differ in the last bits of the double.
LOWELL_PLANET_X_ELEMENTS = _keplerian_view(
    LowellPlanetXElements, PLUTO_LOWELL, i=10.0, M0=178.3
)

# Pickering's legacy display container preserves the rc7 mixed historical
# view: the mean distance and Gaussian rate of his first Planet-O solution
# ("A Search for a Planet beyond Neptune", Annals of the Astronomical
# Observatory of Harvard College 61, Part II (1909), p. 162: epoch 1900.0,
# longitude 105.8 deg, mean distance 51.9, period 373.5, an effectively
# circular two-perturbation solution) combined with the eccentricity,
# inclination and node estimates of his 1919 revision (Annals 82, p. 59).
# M0 is the 1909 epoch longitude measured from the 1919 perihelion,
# (105.8 - 280) mod 360.  The runtime row uses the complete 1919 set instead,
# so the 1909 epoch and mean distance below belong to this container alone.
_PICKERING_1909_MEAN_DISTANCE_AU: float = 51.9
_PICKERING_ROW = HYPOTHETICAL_BODIES[PLUTO_PICKERING]
PICKERING_PLANET_X_ELEMENTS = PickeringPlanetXElements(
    _PICKERING_ROW.name,
    _J1900_EPOCH_JD,
    _PICKERING_1909_MEAN_DISTANCE_AU,
    _PICKERING_ROW.e,
    _PICKERING_ROW.i,
    _PICKERING_ROW.omega,
    _PICKERING_ROW.Omega,
    (105.8 - 280.0) % 360.0,
    _K_GAUSS_DEG / _PICKERING_1909_MEAN_DISTANCE_AU**1.5,
)

_VULCAN_ROW = HYPOTHETICAL_BODIES[VULCAN]
VULCAN_ELEMENTS = VulcanElements(
    _VULCAN_ROW.name,
    _VULCAN_ROW.epoch,
    _VULCAN_ROW.a,
    _VULCAN_ROW.e,
    _VULCAN_ROW.i,
    _VULCAN_ROW.M0,
    _VULCAN_ROW.n_century,
    _VULCAN_ROW.omega,
    _VULCAN_ROW.omega_rate,
    _VULCAN_ROW.Omega,
    _VULCAN_ROW.Omega_rate,
)

_WALDEMATH_ROW = HYPOTHETICAL_BODIES[WALDEMATH]
WALDEMATH_ELEMENTS = WaldemathElements(
    _WALDEMATH_ROW.name,
    _WALDEMATH_ROW.epoch,
    _WALDEMATH_ROW.a,
    0.0,  # circular by construction (uniform longitude model)
    0.0,  # longitude-only model: the point rides the ecliptic
    0.0,  # degenerate for a circular planar orbit (stated-convention zero)
    0.0,  # degenerate for a circular planar orbit (stated-convention zero)
    _WALDEMATH_ROW.M0,
    _WALDEMATH_ROW.n,
    _WALDEMATH_ROW.M0,
)


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
        raise ValueError(f"cannot parse the epoch or equinox value {value!r}")


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
            raise ValueError(f"cannot parse {expr!r} as a number")

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
    given = filepath
    filepath = Path(filepath)

    # exists() plus an explicit directory check, not is_file(): Path("") is
    # Path("."), and a directory exists, so the guard was bypassed and open()
    # raised IsADirectoryError instead of the documented FileNotFoundError.
    # is_file() would close that gap but also reject the special files a
    # caller may legitimately pass (/dev/null, a FIFO, /dev/stdin) and, on
    # some Python versions, report an unreadable path as missing. The message
    # names what the caller actually passed rather than the resolved ".".
    if not filepath.exists() or filepath.is_dir():
        raise FileNotFoundError(f"orbital elements file not found: {given!r}")

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
                raise ValueError(f"cannot parse line {line_num}: {e}") from e

    return elements


def get_bundled_fictitious_orbits_path() -> Path:
    """
    Get the path to the bundled fictitious orbits dataset included with libephemeris.

    The dataset contains the 12 independently transcribed primary-source
    orbital models documented in ``docs/methodology/hypothetical-bodies.md``.
    The other supported bodies (Transpluto, Pickering, Vulcan, Selena,
    Proserpina and Waldemath) are defined as module-level element sets
    rather than CSV rows; the one ID without a recoverable public model
    (Nibiru) is inventoried in
    ``docs/methodology/missing-hypothetical-models.md``.

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
            f"the bundled fictitious-orbits dataset is missing at {data_path}; "
            "the package may be incomplete"
        )
    return data_path


# Columns of the bundled fictitious-orbits dataset, in file order.  The file
# repeats them as its header line so that a stale layout is refused instead of
# being read with the columns silently permuted.
_CSV_COLUMNS = (
    "name",
    "epoch_jd",
    "equinox_token",
    "equinox_jd",
    "a_au",
    "e",
    "i_deg",
    "node_deg",
    "argp_deg",
    "mean_anomaly_deg",
    "source",
)


def _parse_column(text: str, column: str) -> float:
    """Read one numeric column, naming it when the text is not a number.

    Args:
        text: The column's text as written in the file.
        column: Human-readable column name for the error message.

    Returns:
        The column value.

    Raises:
        ValueError: If the text is not a number.
    """
    try:
        return float(text)
    except ValueError:
        raise ValueError(f"cannot parse the {column} {text!r}") from None


def _resolve_equinox(token: str, jd_text: str) -> float:
    """Resolve a row's disjoint equinox pair into a Julian Day TT.

    A row states its equinox either as a standard name from ``_EPOCH_JD``, such
    as ``J1900``, or, for sources whose elements are referred to their own
    date, as an explicit Julian Day.  Exactly one of the two columns carries a
    value.

    Args:
        token: Content of the ``equinox_token`` column, possibly empty.
        jd_text: Content of the ``equinox_jd`` column, possibly empty.

    Returns:
        The equinox as a Julian Day TT.

    Raises:
        ValueError: If both columns are filled, both are empty, or the token
            is not a standard epoch name.
    """
    if bool(token) == bool(jd_text):
        raise ValueError(
            "exactly one of the equinox_token and equinox_jd columns must be filled"
        )
    if not token:
        return _parse_column(jd_text, "equinox Julian Day")
    equinox_jd = _EPOCH_JD.get(token.upper())
    if equinox_jd is None:
        raise ValueError(f"cannot parse the equinox token {token!r}")
    return equinox_jd


def _parse_row(columns: List[str], line_num: int) -> OrbitalElements:
    """Turn one data row of the fictitious-orbits dataset into elements.

    Args:
        columns: The row's fields, already stripped, in ``_CSV_COLUMNS`` order.
        line_num: Line number of the row, quoted in error messages.

    Returns:
        The parsed :class:`OrbitalElements`.

    Raises:
        ValueError: If the row has the wrong number of columns, or a column
            cannot be read.
    """
    if len(columns) != len(_CSV_COLUMNS):
        raise ValueError(
            f"cannot parse line {line_num}: expected {len(_CSV_COLUMNS)} "
            f"comma-separated fields, got {len(columns)}"
        )

    # The trailing source column is prose for reviewers, not runtime data.
    try:
        name, epoch_text, equinox_token, equinox_text = columns[:4]
        epoch_jd = _parse_column(epoch_text, "epoch")
        equinox_jd = _resolve_equinox(equinox_token, equinox_text)
        semi_axis = _parse_column(columns[4], "semi-major axis")
        eccentricity = _parse_column(columns[5], "eccentricity")
        inclination = _parse_column(columns[6], "inclination")
        asc_node = _parse_column(columns[7], "ascending node")
        arg_perihelion = _parse_column(columns[8], "argument of perihelion")
        mean_anomaly = _parse_column(columns[9], "mean anomaly")
    except ValueError as exc:
        raise ValueError(f"cannot parse line {line_num}: {exc}") from exc

    # Every row is a fixed-equinox heliocentric orbit stated as plain numbers,
    # so the dataset carries neither an equinox-of-date marker, nor a
    # geocentric flag, nor secular rates.
    return OrbitalElements(
        name=name,
        epoch_jd=epoch_jd,
        equinox_jd=equinox_jd,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(constant=mean_anomaly),
        semi_axis=semi_axis,
        eccentricity=TPolynomial(constant=eccentricity),
        arg_perihelion=TPolynomial(constant=arg_perihelion),
        asc_node=TPolynomial(constant=asc_node),
        inclination=TPolynomial(constant=inclination),
        line_number=line_num,
    )


def _parse_fictitious_orbits_csv(filepath: Union[str, Path]) -> List[OrbitalElements]:
    """
    Parse a ``fictitious_orbits.csv`` dataset into OrbitalElements objects.

    The file has eleven columns, named on a header line that precedes the data
    rows::

        name, epoch_jd, equinox_token, equinox_jd, a_au, e, i_deg, node_deg,
        argp_deg, mean_anomaly_deg, source

    Blank lines and lines starting with ``#`` are ignored, so the file can
    document each transcription beside its row.  The final ``source`` column
    runs to the end of the line and may contain commas.

    Args:
        filepath: Path to the CSV file.

    Returns:
        List of :class:`OrbitalElements` objects, one per data row.

    Raises:
        FileNotFoundError: If the file does not exist.
        ValueError: If the header line is missing or a data row cannot be
            parsed.
    """
    filepath = Path(filepath)
    if not filepath.exists():
        raise FileNotFoundError(f"orbital elements file not found: {filepath}")

    elements: List[OrbitalElements] = []
    header_seen = False

    with filepath.open(encoding="utf-8") as fh:
        for line_num, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue

            columns = [part.strip() for part in line.split(",", len(_CSV_COLUMNS) - 1)]
            if not header_seen:
                if tuple(columns) != _CSV_COLUMNS:
                    raise ValueError(
                        f"{filepath} does not start with the fictitious-orbits "
                        "column header"
                    )
                header_seen = True
                continue

            elements.append(_parse_row(columns, line_num))

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
            f"expected at least 9 comma-separated fields, got {len(parts)}"
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
        raise ValueError(
            f"the epoch {epoch_str!r} is JDATE, but a fixed epoch is required"
        )

    equinox_jd, equinox_is_jdate = _parse_epoch_or_equinox(equinox_str)

    # Parse orbital elements (may be T-polynomials)
    mean_anomaly = _parse_t_polynomial(mean_anomaly_str)

    # Semi-axis is always a simple number
    try:
        semi_axis = float(semi_axis_str)
    except ValueError:
        raise ValueError(f"cannot parse the semi-major axis {semi_axis_str!r}")

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
    elements: HypotheticalBody, jd_tt: float
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


# Public view of the classical-prediction rows in the shape the orbital
# elements parser produces, so a caller reads the historical predictions the
# same way it reads a custom file.  The labels are the short names these rows
# have always carried; every number comes from the table.
_CLASSICAL_ELEMENT_LABELS: Dict[int, str] = {
    HARRINGTON: "Harrington",
    NEPTUNE_LEVERRIER: "Leverrier",
    NEPTUNE_ADAMS: "Adams",
    PLUTO_LOWELL: "Lowell",
    PLUTO_PICKERING: "Pickering",
}


def _classical_orbital_elements(body_id: int, label: str) -> OrbitalElements:
    """Project one classical-prediction row onto an ``OrbitalElements``."""
    body = HYPOTHETICAL_BODIES[body_id]
    return OrbitalElements(
        name=label,
        epoch_jd=body.epoch,
        equinox_jd=body.equinox_jd,
        equinox_is_jdate=False,
        mean_anomaly=TPolynomial(body.M0),
        semi_axis=body.a,
        eccentricity=TPolynomial(body.e),
        arg_perihelion=TPolynomial(body.omega),
        asc_node=TPolynomial(body.Omega),
        inclination=TPolynomial(body.i),
    )


FICTITIOUS_ORBITAL_ELEMENTS: Dict[int, OrbitalElements] = {
    body_id: _classical_orbital_elements(body_id, label)
    for body_id, label in _CLASSICAL_ELEMENT_LABELS.items()
}


def calc_fictitious_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """Return a reviewed classical prediction as a heliocentric J2000 state.

    The registry contains Harrington, Le Verrier, Adams, Lowell and
    Pickering.  Other classical-prediction compatibility IDs raise
    :class:`UnknownBodyError`.

    Returns (lon, lat, dist, dlon, dlat, ddist), J2000 ecliptic.
    """
    body = HYPOTHETICAL_BODIES.get(ipl)
    if body is None or body.model != MODEL_CLASSICAL:
        from .exceptions import UnknownBodyError

        raise UnknownBodyError(
            message=f"no fictitious elements for body {ipl}", body_id=ipl
        )

    # Half-day centred difference: the outer samples are one day apart, so
    # their difference is already the daily rate and needs no division.
    lon, lat, dist = _keplerian_to_ecliptic_j2000(body, jd_tt)
    h = 0.5
    lon_p, lat_p, dist_p = _keplerian_to_ecliptic_j2000(body, jd_tt - h)
    lon_n, lat_n, dist_n = _keplerian_to_ecliptic_j2000(body, jd_tt + h)
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
    return calc_uranian_planet(ipl, jd_tt)[0]


def calc_uranian_position(
    ipl: int, jd_tt: float
) -> Tuple[float, float, float, float, float, float]:
    """Calculate the state of a registered Hamburg-school body."""
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
    """Calculate Transpluto (Isis) from the published Hawkins element set.

    Propagates the Sevin/Landscheidt/Hawkins J1900 orbit (a = 77.755 AU,
    e = 0.3, planar) with Kepler's equation and converts the longitude from
    the J1900 equinox to J2000.  The point is reported on the ecliptic
    (latitude exactly zero), the historical output convention for this body.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT).

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist),
        heliocentric J2000 ecliptic.
    """
    lon, lat, dist = _calc_transpluto_raw(jd_tt)
    dt_step = 1.0
    lon_p, _lat_p, dist_p = _calc_transpluto_raw(jd_tt - dt_step)
    lon_n, _lat_n, dist_n = _calc_transpluto_raw(jd_tt + dt_step)
    dlon = ((lon_n - lon_p + 180.0) % 360.0 - 180.0) / (2.0 * dt_step)
    ddist = (dist_n - dist_p) / (2.0 * dt_step)
    return (lon, lat, dist, dlon, 0.0, ddist)


def _calc_transpluto_raw(jd_tt: float) -> Tuple[float, float, float]:
    """Geometric heliocentric J2000 position of the published Transpluto."""
    elements = HYPOTHETICAL_BODIES[ISIS]
    mean_anomaly = math.radians(
        (elements.M0 + elements.n * (jd_tt - elements.epoch)) % 360.0
    )
    eccentric_anomaly = _solve_kepler_equation(mean_anomaly, elements.e)
    true_anomaly = 2.0 * math.atan(
        math.sqrt((1.0 + elements.e) / (1.0 - elements.e))
        * math.tan(eccentric_anomaly / 2.0)
    )
    distance = elements.a * (1.0 - elements.e * math.cos(eccentric_anomaly))
    # i = node = 0: the orbit lies in the J1900 ecliptic and the longitude is
    # perihelion longitude + true anomaly in that frame.
    lon_j1900 = (elements.omega + math.degrees(true_anomaly)) % 360.0
    from .astrometry import _precess_ecliptic

    lon_j2000, _ = _precess_ecliptic(lon_j1900, 0.0, elements.equinox_jd, 2451545.0)
    return (float(lon_j2000), 0.0, float(distance))


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
    elements = _table_row(body_id, MODEL_HAMBURG)

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
    return _keplerian_to_ecliptic_j2000(_table_row(body_id, MODEL_HAMBURG), jd_tt)


def _raise_unsupported_builtin_fictitious(body_id: int) -> NoReturn:
    """Raise for an ID absent from independently reviewed registries."""
    from .exceptions import UnknownBodyError

    raise UnknownBodyError(
        message=f"no built-in fictitious elements for body {body_id}",
        body_id=body_id,
    )


def _table_row(body_id: int, model: str) -> HypotheticalBody:
    """Return the table row of ``body_id``, or fail closed.

    Args:
        body_id: Hypothetical body ID.
        model: The model the caller knows how to propagate.

    Raises:
        UnknownBodyError: The ID has no row, or its row has another model.
    """
    body = HYPOTHETICAL_BODIES.get(body_id)
    if body is None or body.model != model:
        _raise_unsupported_builtin_fictitious(body_id)
    return body


def calc_vulcan(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the intramercurial Vulcan (heliocentric, of-date).

    Propagates the Vulcan row of ``HYPOTHETICAL_BODIES`` (compatibility
    values, primary source not verified): mean anomaly, argument of
    perihelion and node move linearly in Julian
    centuries from the J1900 epoch, and the orbit is referred to the ecliptic
    and mean equinox of date.  The returned state is geometric (no
    light-time); apparent reductions are applied by the ``calc_ut`` layer.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT).

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist) with the
        angles in degrees, the distance in AU and speeds per day.
    """
    lon, lat, dist = _calc_vulcan_raw(jd_tt)
    # Central difference with a tight step: Vulcan's sidereal period is only
    # ~18.6 days, so a one-day stencil would materially smooth the speed.
    dt_step = 0.01
    lon_p, lat_p, dist_p = _calc_vulcan_raw(jd_tt - dt_step)
    lon_n, lat_n, dist_n = _calc_vulcan_raw(jd_tt + dt_step)
    dlon = ((lon_n - lon_p + 180.0) % 360.0 - 180.0) / (2.0 * dt_step)
    dlat = (lat_n - lat_p) / (2.0 * dt_step)
    ddist = (dist_n - dist_p) / (2.0 * dt_step)
    return (lon, lat, dist, dlon, dlat, ddist)


def _calc_vulcan_raw(jd_tt: float) -> Tuple[float, float, float]:
    """Geometric heliocentric of-date position of Vulcan (ID 55)."""
    elements = HYPOTHETICAL_BODIES[VULCAN]
    T = (jd_tt - elements.epoch) / 36525.0
    mean_anomaly = math.radians((elements.M0 + elements.n_century * T) % 360.0)
    argument = math.radians((elements.omega + elements.omega_rate * T) % 360.0)
    node = math.radians((elements.Omega + elements.Omega_rate * T) % 360.0)
    inclination = math.radians(elements.i)

    eccentric_anomaly = _solve_kepler_equation(mean_anomaly, elements.e)
    x_orbit = elements.a * (math.cos(eccentric_anomaly) - elements.e)
    y_orbit = (
        elements.a
        * math.sqrt(1.0 - elements.e * elements.e)
        * math.sin(eccentric_anomaly)
    )

    cos_arg, sin_arg = math.cos(argument), math.sin(argument)
    cos_node, sin_node = math.cos(node), math.sin(node)
    cos_inc, sin_inc = math.cos(inclination), math.sin(inclination)
    x = (cos_node * cos_arg - sin_node * sin_arg * cos_inc) * x_orbit + (
        -cos_node * sin_arg - sin_node * cos_arg * cos_inc
    ) * y_orbit
    y = (sin_node * cos_arg + cos_node * sin_arg * cos_inc) * x_orbit + (
        -sin_node * sin_arg + cos_node * cos_arg * cos_inc
    ) * y_orbit
    z = sin_arg * sin_inc * x_orbit + cos_arg * sin_inc * y_orbit

    distance = math.sqrt(x * x + y * y + z * z)
    longitude = math.degrees(math.atan2(y, x)) % 360.0
    latitude = math.degrees(math.asin(z / distance))
    return (float(longitude), float(latitude), float(distance))


def calc_waldemath(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Return the published uniform Waldemath/Sepharial Dark Moon model.

    Sepharial (The Science of Foreknowledge, 1918, ch. "The New Satellite —
    Lilith") defines the point by uniform 3 deg/day geocentric tropical
    motion from the Sun's longitude at a printed Lilith-Sun conjunction;
    the anchor here is Waltemath's predicted 1898 February 2 transit at
    00:00 GMT (see the WALDEMATH model block above for the full sources).
    The stored anchor is the Sun's APPARENT longitude at that instant, so
    the nutation at the anchor epoch is removed to yield the mean-of-date
    longitude this layer must return (the apparent-place layer in calc_ut
    adds the nutation of the request date, exactly as for Selena).
    """
    import math as _math

    from .cache import get_cached_nutation
    from .time_utils import deltat

    elements = HYPOTHETICAL_BODIES[WALDEMATH]
    anchor_tt = elements.epoch + deltat(elements.epoch)
    dpsi_rad, _ = get_cached_nutation(anchor_tt)
    anchor_mean = elements.M0 - _math.degrees(dpsi_rad)
    longitude = (anchor_mean + elements.n * (jd_tt - anchor_tt)) % 360.0
    return (
        float(longitude),
        0.0,
        float(elements.a),
        float(elements.n),
        0.0,
        0.0,
    )


def calc_transpluto_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate Transpluto (Isis); alias entry point for ID 48.

    See :func:`calc_transpluto` for the published element set and frame.
    """
    return calc_transpluto(jd_tt)


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
    longitude, latitude, distance = _calc_keplerian_position_raw(ipl, jd_tt)

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
        raise ValueError(f"body id {ipl} has no Keplerian elements")

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
    elements = HYPOTHETICAL_BODIES[WHITE_MOON]
    elapsed_days = jd_tt - elements.epoch
    longitude = (elements.M0 + elements.n * elapsed_days) % 360.0
    return (
        float(longitude),
        0.0,
        float(elements.a),
        float(elements.n),
        0.0,
        0.0,
    )


def calc_waldemath_position(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Return the Waldemath Dark Moon; alias entry point for ID 58.

    See :func:`calc_waldemath` for the published model and its sources.
    """
    return calc_waldemath(jd_tt)


def calc_proserpina(jd_tt: float) -> Tuple[float, float, float, float, float, float]:
    """Calculate the trans-Plutonian Proserpina (heliocentric, of-date).

    The convention realized here is a uniform circular ecliptic orbit:
    a = 79.22563 AU, mean longitude 170.73 degrees at the J1900 epoch,
    Gaussian mean motion (~51.05 deg/Julian century, period ~705 years),
    referred to the ecliptic and mean equinox of date.  These are
    compatibility values, primary source not verified.  Being exactly
    circular and planar, the latitude and the radial speed are zero.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT).

    Returns:
        Tuple of (longitude, latitude, distance, dlon, dlat, ddist).
    """
    elements = HYPOTHETICAL_BODIES[PROSERPINA]
    longitude = (elements.M0 + elements.n * (jd_tt - elements.epoch)) % 360.0
    return (float(longitude), 0.0, float(elements.a), float(elements.n), 0.0, 0.0)


def calc_planet_x_lowell(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate Lowell's preferred 1915 Planet-X solution."""
    return calc_fictitious_position(PLUTO_LOWELL, jd_tt)


def calc_planet_x_pickering(
    jd_tt: float,
) -> Tuple[float, float, float, float, float, float]:
    """Calculate Pickering's 1919 Planet O prediction.

    Uses the complete element list of Annals of Harvard College Observatory
    82, No. 3 (1919), p. 59 (a = 55.1 AU, e = 0.31, perihelion longitude
    280 deg, perihelion date 1720.0, node 100 deg, inclination 15 deg).
    """
    return calc_fictitious_position(PLUTO_PICKERING, jd_tt)


#: Which function propagates each model.  Every propagator takes the body ID
#: and the Julian Day; the single-body models ignore the ID.
_MODEL_PROPAGATORS: Dict[
    str, Callable[[int, float], Tuple[float, float, float, float, float, float]]
] = {
    MODEL_HAMBURG: calc_uranian_planet,
    MODEL_CLASSICAL: calc_fictitious_position,
    MODEL_TRANSPLUTO: lambda _body_id, jd_tt: calc_transpluto_position(jd_tt),
    MODEL_VULCAN: lambda _body_id, jd_tt: calc_vulcan(jd_tt),
    MODEL_SELENA: lambda _body_id, jd_tt: calc_white_moon_position(jd_tt),
    MODEL_DARK_MOON: lambda _body_id, jd_tt: calc_waldemath_position(jd_tt),
    MODEL_CIRCULAR: lambda _body_id, jd_tt: calc_proserpina(jd_tt),
}


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
        UnknownBodyError: If the recognised ID has no publicly sourced
            built-in model.  Supported IDs are 40--48 and 50--58; ID 49
            (Nibiru, see
            ``docs/methodology/missing-hypothetical-models.md``) and the
            unnamed tail 59--69 of the recognised fictitious range fail
            closed.
    """
    if not is_hypothetical_body(ipl):
        raise ValueError(f"body id {ipl} is not a hypothetical body")

    body = HYPOTHETICAL_BODIES.get(ipl)
    propagate = None if body is None else _MODEL_PROPAGATORS.get(body.model)
    if propagate is None:
        _raise_unsupported_builtin_fictitious(ipl)
    return propagate(ipl, jd_tt)


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
