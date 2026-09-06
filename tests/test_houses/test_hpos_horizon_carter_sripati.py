# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Structural identities of the horizontal, Carter and Sripati positions.

Each of the three systems is fixed by naming one circle, one origin on it and
one direction of counting: the horizon from its east point for ``H``, the
celestial equator from the ascendant's right ascension for ``F``, the ecliptic
from half a house before the ascendant for ``S``. The identities that follow
from those three statements need no reference values at all, and they are what
this file asserts.

The obliquity is the axis the rest of the suite barely moves: a chart carries
the true obliquity of date, so every recorded frame sits within a twentieth of
a degree of 23.44. Each system is therefore also exercised at an obliquity far
from that band, where a construction that had quietly folded the obliquity
into the wrong step would part company with the geometry.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as ephem
from libephemeris.houses import (
    _armc_to_mc,
    _eastern_ascendant,
    _ecliptic_to_ra_simple,
    _hpos_carter,
    _hpos_horizon,
)

#: Obliquities far below and far above the band the charts occupy. The polar
#: circle of the larger one reaches 30 degrees of latitude, so the frames it
#: is used on stay inside the tropics.
SMALL_EPS = 5.0
LARGE_EPS = 60.0
#: The band every recorded chart sits in.
CHART_EPS = 23.4392911

ARMCS = (0.0, 47.0, 123.0, 211.0, 305.0)


def _horizontal_to_equatorial(
    azimuth: float, altitude: float, geolat: float
) -> tuple[float, float]:
    """Place a body by its horizon coordinates.

    The inverse of the rotation the horizontal system reads its answer off:
    the azimuth runs from the north point of the horizon towards the east
    point, and the arc returned is the body's distance from the upper
    meridian, whose negative is the hour angle.
    """
    north = math.cos(math.radians(altitude)) * math.cos(math.radians(azimuth))
    east = math.cos(math.radians(altitude)) * math.sin(math.radians(azimuth))
    up = math.sin(math.radians(altitude))
    dec = math.degrees(
        math.asin(
            math.cos(math.radians(geolat)) * north + math.sin(math.radians(geolat)) * up
        )
    )
    hour = math.degrees(
        math.atan2(
            -east,
            -math.sin(math.radians(geolat)) * north
            + math.cos(math.radians(geolat)) * up,
        )
    )
    return -hour, dec


def _equatorial_to_ecliptic(ra: float, dec: float, eps: float) -> tuple[float, float]:
    """Carry a body back from the equator to the ecliptic."""
    lat = math.degrees(
        math.asin(
            math.sin(math.radians(dec)) * math.cos(math.radians(eps))
            - math.cos(math.radians(dec))
            * math.sin(math.radians(eps))
            * math.sin(math.radians(ra))
        )
    )
    lon = (
        math.degrees(
            math.atan2(
                math.sin(math.radians(ra)) * math.cos(math.radians(eps))
                + math.tan(math.radians(dec)) * math.sin(math.radians(eps)),
                math.cos(math.radians(ra)),
            )
        )
        % 360.0
    )
    return lon, lat


@pytest.mark.unit
class TestHorizontalPosition:
    """``H``: the horizon divided from its east point."""

    @pytest.mark.parametrize(
        "azimuth,cusp", [(90.0, 1.0), (0.0, 4.0), (270.0, 7.0), (180.0, 10.0)]
    )
    @pytest.mark.parametrize("geolat", [-89.0, -45.0, -10.0, 0.0, 10.0, 45.0, 89.0])
    def test_cardinal_points_of_the_horizon_are_cusps(self, azimuth, cusp, geolat):
        """The east point opens the first house, the north point is the
        fourth cusp, the west point the seventh and the south point the
        tenth -- at every latitude, in both hemispheres alike, because the
        ring is anchored to the horizon and not to the meridian."""
        md_upper, dec = _horizontal_to_equatorial(azimuth, 17.0, geolat)
        assert _hpos_horizon(md_upper, dec, geolat) == pytest.approx(cusp, abs=1e-9)

    @pytest.mark.parametrize("geolat", [-70.0, -23.44, 0.0, 23.44, 70.0])
    @pytest.mark.parametrize("azimuth", [12.0, 97.0, 188.0, 264.0, 341.0])
    def test_the_altitude_does_not_enter(self, geolat, azimuth):
        """Every house of this system is bounded by vertical circles, so
        sliding a body along its own vertical circle cannot move it."""
        answers = [
            _hpos_horizon(*_horizontal_to_equatorial(azimuth, altitude, geolat), geolat)
            for altitude in (-80.0, -35.0, 0.0, 35.0, 80.0)
        ]
        for answer in answers[1:]:
            assert answer == pytest.approx(answers[0], abs=1e-11)

    @pytest.mark.parametrize("geolat", [-62.0, -20.0, 0.0, 20.0, 62.0])
    @pytest.mark.parametrize("armc", ARMCS)
    def test_the_obliquity_does_not_enter(self, geolat, armc):
        """The obliquity belongs to the conversion that precedes this system,
        not to the system: one equatorial place read through three very
        different obliquities is three different ecliptic places and one
        house position."""
        answers = []
        for eps in (SMALL_EPS, CHART_EPS, LARGE_EPS):
            lon, lat = _equatorial_to_ecliptic(armc + 143.0, 31.0, eps)
            answers.append(ephem.house_pos(armc, geolat, eps, ord("H"), lon, lat))
        for answer in answers[1:]:
            assert answer == pytest.approx(answers[0], abs=1e-11)


@pytest.mark.unit
class TestCarterPosition:
    """``F``: the equator divided from the ascendant's right ascension."""

    @pytest.mark.parametrize("eps,geolat", [(SMALL_EPS, 62.0), (LARGE_EPS, 20.0)])
    @pytest.mark.parametrize("armc", ARMCS)
    def test_the_ascendant_and_the_descendant_are_cusps(self, eps, geolat, armc):
        """The two anchors of the system, at an obliquity far from the band
        the charts occupy. The meridian is not a cusp here, so these two are
        the whole of what is fixed."""
        for latitude in (geolat, -geolat):
            asc = _eastern_ascendant(armc, eps, latitude)
            assert ephem.house_pos(
                armc, latitude, eps, ord("F"), asc, 0.0
            ) == pytest.approx(1.0, abs=1e-9)
            assert ephem.house_pos(
                armc, latitude, eps, ord("F"), (asc + 180.0) % 360.0, 0.0
            ) == pytest.approx(7.0, abs=1e-9)

    @pytest.mark.parametrize("eps,geolat", [(SMALL_EPS, 62.0), (LARGE_EPS, 20.0)])
    @pytest.mark.parametrize("ra", [0.0, 61.0, 137.0, 249.0, 318.0])
    def test_ninety_degrees_of_right_ascension_is_three_houses(self, eps, geolat, ra):
        """The equator is cut into twelve equal arcs of right ascension, so a
        quarter of it is three houses however the ring is placed."""
        for armc in ARMCS:
            here = _hpos_carter(ra, armc, eps, geolat)
            there = _hpos_carter((ra + 90.0) % 360.0, armc, eps, geolat)
            assert (there - here) % 12.0 == pytest.approx(3.0, abs=1e-12)

    @pytest.mark.parametrize("eps,geolat", [(SMALL_EPS, 62.0), (LARGE_EPS, 20.0)])
    @pytest.mark.parametrize("armc", ARMCS)
    def test_the_declination_does_not_enter(self, eps, geolat, armc):
        """The projection runs along hour circles: two bodies on one hour
        circle share a house position whatever their distance from the
        equator."""
        answers = []
        for dec in (-55.0, 0.0, 55.0):
            lon, lat = _equatorial_to_ecliptic(armc + 76.0, dec, eps)
            answers.append(ephem.house_pos(armc, geolat, eps, ord("F"), lon, lat))
        for answer in answers[1:]:
            assert answer == pytest.approx(answers[0], abs=1e-11)


@pytest.mark.unit
class TestSripatiPosition:
    """``S``: the Porphyry ring pulled back by half a house."""

    @pytest.mark.parametrize("eps,geolat", [(SMALL_EPS, 62.0), (LARGE_EPS, 20.0)])
    @pytest.mark.parametrize("armc", ARMCS)
    def test_the_four_angles_fall_in_the_middle_of_their_houses(
        self, eps, geolat, armc
    ):
        """Every Sripati cusp sits at the midpoint of a Porphyry house, so
        the ascendant, the IC, the descendant and the MC land half a house
        after the cusps they would open in Porphyry."""
        for latitude in (geolat, -geolat):
            asc = _eastern_ascendant(armc, eps, latitude)
            mc = _armc_to_mc(armc, eps)
            for angle, expected in (
                (asc, 1.5),
                ((mc + 180.0) % 360.0, 4.5),
                ((asc + 180.0) % 360.0, 7.5),
                (mc, 10.5),
            ):
                assert ephem.house_pos(
                    armc, latitude, eps, ord("S"), angle, 0.0
                ) == pytest.approx(expected, abs=1e-9)

    @pytest.mark.parametrize("eps,geolat", [(SMALL_EPS, 62.0), (LARGE_EPS, 20.0)])
    @pytest.mark.parametrize("armc", ARMCS)
    def test_it_is_the_porphyry_position_advanced_by_half_a_house(
        self, eps, geolat, armc
    ):
        """The construction stated as an identity against the system it is
        built on, taken round the circle so that the last half of the twelfth
        Porphyry house opens the first Sripati one."""
        for lon in range(0, 360, 23):
            porphyry = ephem.house_pos(armc, geolat, eps, ord("O"), float(lon), 0.0)
            assert (
                ephem.house_pos(armc, geolat, eps, ord("S"), float(lon), 0.0)
                == ((porphyry - 1.0 + 0.5) % 12.0) + 1.0
            )

    @pytest.mark.parametrize("eps,geolat", [(SMALL_EPS, 62.0), (LARGE_EPS, 20.0)])
    @pytest.mark.parametrize("armc", ARMCS)
    def test_the_ecliptic_latitude_does_not_enter(self, eps, geolat, armc):
        """The system divides the ecliptic itself, so nothing off the
        ecliptic can reach the answer: the equality is bit for bit, not
        approximate."""
        for lon in range(0, 360, 29):
            base = ephem.house_pos(armc, geolat, eps, ord("S"), float(lon), 0.0)
            for lat_body in (-17.0, -4.0, 4.0, 17.0):
                assert (
                    ephem.house_pos(armc, geolat, eps, ord("S"), float(lon), lat_body)
                    == base
                )


@pytest.mark.unit
class TestCarterAscendantIsShared:
    """``F`` reads the ascendant the module already owns."""

    @pytest.mark.parametrize("geolat", [-78.0, -40.0, 0.0, 40.0, 78.0])
    def test_the_origin_is_the_ascendant_carried_to_the_equator(self, geolat):
        """The eastern test in the ascendant is what keeps the ring from
        landing six houses out inside the polar circles, so the origin has to
        be that ascendant and not the closed form nearer the meridian."""
        for armc in ARMCS:
            asc_ra = _ecliptic_to_ra_simple(
                _eastern_ascendant(armc, CHART_EPS, geolat), CHART_EPS
            )
            for ra in (11.0, 98.0, 203.0, 299.0):
                assert _hpos_carter(ra, armc, CHART_EPS, geolat) == pytest.approx(
                    ((ra - asc_ra) % 360.0) / 30.0 + 1.0, abs=1e-12
                )
