# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Guards on the house-position dispatcher behind :func:`house_pos`.

Three groups of checks live here.

*Degenerate polar cells.* At a geographic pole the Regiomontanus house
circles are no longer distinguishable and the construction collapses onto the
fourth and the tenth cusp; a body on the celestial equator is the exception,
because there the unbounded tangent of the latitude multiplies a vanishing one
and the closed form stays finite. Those cells, and the Koch cells at the same
configuration, are pinned one by one so that no registered tolerance on the
polar rows can hide a change in them.

*The ecliptic poles.* No chart puts a body at an ecliptic pole, so nothing
else exercises the unbounded tangent in the ecliptic-to-equatorial conversion.
The geometry is unambiguous and is asserted directly.

*Structural identities.* The relations that hold by construction and need no
reference values: the Gauquelin sector is the Placidus position re-expressed,
the aliases answer identically, the two calling conventions agree bit for bit,
a body on a cusp is answered with that cusp's number, and a Placidus position
reaches the upper half of the circle exactly while the body is above the
horizon.
"""

from __future__ import annotations

import itertools
import math

import pytest

import libephemeris as ephem

#: The frame the recorded polar grid uses.
POLAR_ARMC = 197.5428262395338
POLAR_EPS = 23.438088989372396


@pytest.mark.unit
class TestPolarRegiomontanus:
    """The two-cusp collapse at a geographic pole, and its one exception."""

    @pytest.mark.parametrize("geolat", [-90.0, 90.0])
    def test_equinoctial_body_keeps_the_finite_closed_form(self, geolat):
        """A body of zero declination does not collapse: the product of the
        unbounded tangent of the latitude and the vanishing tangent of the
        declination is zero, and the position is the body's own meridian
        distance carried a quarter turn."""
        assert (
            ephem.house_pos(POLAR_ARMC, geolat, POLAR_EPS, ord("R"), 0.0, 0.0)
            == 3.415239125348873
        )
        assert (
            ephem.house_pos(POLAR_ARMC, geolat, POLAR_EPS, (0.0, 0.0), b"R")
            == 3.415239125348873
        )

    @pytest.mark.parametrize(
        "geolat,body_lon,body_lat,expected",
        [
            # Off the celestial equator the circles collapse: the tenth cusp
            # above the horizon, the fourth below it. At the north pole the
            # visible half is the northern one.
            (90.0, 90.0, 10.0, 10.0),
            (-90.0, 90.0, 10.0, 4.0),
            (90.0, 270.0, -10.0, 4.0),
            (-90.0, 270.0, -10.0, 10.0),
            # The autumn equinox: its declination is a few units in the last
            # place rather than zero, so the collapse applies there too.
            (90.0, 180.0, 0.0, 10.0),
            (-90.0, 180.0, 0.0, 4.0),
        ],
    )
    def test_collapse_onto_the_meridian_pair(
        self, geolat, body_lon, body_lat, expected
    ):
        got = ephem.house_pos(
            POLAR_ARMC, geolat, POLAR_EPS, ord("R"), body_lon, body_lat
        )
        assert got == expected

    @pytest.mark.parametrize("geolat,expected", [(90.0, 10.0), (-90.0, 4.0)])
    def test_the_autumn_equinox_collapses_like_any_other_body(self, geolat, expected):
        """The conversion of the ecliptic place ``(180, 0)`` leaves a
        declination of a few units in the last place rather than zero, so the
        body is off the celestial equator and the collapse applies to it. A
        construction that reaches the pole through a clamped latitude answers
        an ordinary position there instead, and the two are far apart."""
        assert (
            ephem.house_pos(POLAR_ARMC, geolat, POLAR_EPS, ord("R"), 180.0, 0.0)
            == expected
        )
        assert (
            ephem.house_pos(POLAR_ARMC, geolat, POLAR_EPS, (180.0, 0.0), b"R")
            == expected
        )

    def test_every_polar_body_off_the_equator_lands_on_a_cusp(self):
        """The collapse is the whole polar row, not a handful of cells."""
        for geolat, body_lat in itertools.product(
            (90.0, -90.0), (-10.0, -5.0, 5.0, 10.0)
        ):
            for body_lon in range(0, 360, 5):
                got = ephem.house_pos(
                    POLAR_ARMC, geolat, POLAR_EPS, ord("R"), float(body_lon), body_lat
                )
                assert got in (4.0, 10.0)


@pytest.mark.unit
class TestPolarKoch:
    """The Koch cells at a pole with the body at an equinoctial point.

    The body's ascensional difference is the product of the unbounded tangent
    of the latitude and the vanishing tangent of the declination. The autumn
    equinox carries a declination of a few units in the last place rather than
    zero, so it is circumpolar and its ascensional difference is a full
    quarter turn; the spring equinox carries an exact zero and keeps an
    ascensional difference of zero. The two answers differ by six houses and
    are pinned here so that no tolerance on the polar rows can hide a change.
    """

    @pytest.mark.parametrize(
        "eps",
        [
            23.869822743656616,
            24.224964526003276,
            23.81159319154222,
            23.437679741216634,
            22.704643125829467,
            22.806677013875014,
            23.438088989372396,
        ],
    )
    def test_autumn_equinox_at_the_south_pole(self, eps):
        assert (
            ephem.house_pos(POLAR_ARMC, -90.0, eps, ord("K"), 180.0, 0.0)
            == 6.707619562674436
        )

    def test_spring_equinox_at_the_south_pole(self):
        assert (
            ephem.house_pos(POLAR_ARMC, -90.0, POLAR_EPS, ord("K"), 0.0, 0.0)
            == 2.2076195626744366
        )

    def test_both_calling_conventions_agree_there(self):
        six = ephem.house_pos(POLAR_ARMC, -90.0, POLAR_EPS, ord("K"), 180.0, 0.0)
        five = ephem.house_pos(POLAR_ARMC, -90.0, POLAR_EPS, (180.0, 0.0), b"K")
        assert six == five == 6.707619562674436


@pytest.mark.unit
class TestEclipticPoles:
    """A body at an ecliptic pole is a legitimate input.

    It sits at right ascension 270 or 90 and declination plus or minus
    ``90 - eps``; nothing in the conversion is undefined, although the tangent
    of the ecliptic latitude is unbounded there.
    """

    @pytest.mark.parametrize("eps", [23.438088989372396, 23.45654434517177])
    @pytest.mark.parametrize("body_lon", [0.0, 71.3, 180.0, 302.6])
    def test_north_ecliptic_pole_is_on_the_meridian_system_scale(self, eps, body_lon):
        # The Meridian system places a body by its own right ascension, so it
        # reads the pole's right ascension back directly.
        armc = 33.75
        got = ephem.house_pos(armc, 40.0, eps, ord("X"), body_lon, 90.0)
        expected = ((270.0 - armc - 90.0) % 360.0) / 30.0 + 1.0
        assert got == pytest.approx(expected, abs=1e-12)

    @pytest.mark.parametrize("eps", [23.438088989372396, 23.45654434517177])
    def test_south_ecliptic_pole(self, eps):
        armc = 33.75
        got = ephem.house_pos(armc, 40.0, eps, ord("X"), 123.4, -90.0)
        expected = ((90.0 - armc - 90.0) % 360.0) / 30.0 + 1.0
        assert got == pytest.approx(expected, abs=1e-12)

    def test_the_two_poles_are_half_a_circle_apart(self):
        eps = 23.438088989372396
        north = ephem.house_pos(120.0, 51.5, eps, ord("P"), 200.0, 90.0)
        south = ephem.house_pos(120.0, 51.5, eps, ord("P"), 200.0, -90.0)
        assert north != south
        for pos in (north, south):
            assert 1.0 <= pos < 13.0


# A grid deliberately free of cardinal coincidences: a body exactly on an
# angle sits on a house boundary, where a half-open interval and a closed
# above-the-horizon test necessarily disagree about which side it is on.
GRID_ARMC = [7.3, 63.7, 128.9, 211.4, 296.1, 341.8]
GRID_LAT = [1.7, 19.3, 38.6, -27.4, 55.2, -61.9]
GRID_EPS = [23.438088989372396, 23.45654434517177]
GRID_LON = [3.7, 51.2, 118.6, 176.4, 233.9, 301.5, 349.1]
GRID_BLAT = [-9.3, -2.1, 0.0, 4.7, 16.2]

ALL_SELECTORS = "PKORCAEWBMXHTGJVUYFLQSDNIi"


@pytest.mark.unit
class TestStructuralIdentities:
    """Relations that hold by construction, checked without reference values."""

    def test_gauquelin_is_the_placidus_position_re_expressed(self):
        for armc, geolat, eps, lon, blat in itertools.product(
            GRID_ARMC, GRID_LAT + [90.0, -90.0], GRID_EPS[:1], GRID_LON, GRID_BLAT
        ):
            placidus = ephem.house_pos(armc, geolat, eps, ord("P"), lon, blat)
            sector = ephem.house_pos(armc, geolat, eps, ord("G"), lon, blat)
            expected = ((360.0 - 30.0 * (placidus - 1.0)) % 360.0) / 10.0 + 1.0
            assert sector == expected
            assert 1.0 <= sector < 37.0

    @pytest.mark.parametrize("left,right", [("A", "E"), ("G", "g"), ("I", "i")])
    def test_aliases_answer_identically(self, left, right):
        for armc, geolat, eps, lon, blat in itertools.product(
            GRID_ARMC, GRID_LAT, GRID_EPS, GRID_LON[:4], GRID_BLAT[:3]
        ):
            a = ephem.house_pos(armc, geolat, eps, ord(left), lon, blat)
            b = ephem.house_pos(armc, geolat, eps, ord(right), lon, blat)
            assert a == b

    def test_the_two_calling_conventions_agree_bitwise(self):
        for armc, geolat, eps, lon, blat in itertools.product(
            GRID_ARMC[:3], GRID_LAT, GRID_EPS[:1], GRID_LON[:4], GRID_BLAT
        ):
            for selector in ALL_SELECTORS:
                six = ephem.house_pos(armc, geolat, eps, ord(selector), lon, blat)
                five = ephem.house_pos(
                    armc, geolat, eps, (lon, blat), selector.encode("latin-1")
                )
                assert six == five

    def test_a_one_element_place_means_zero_ecliptic_latitude(self):
        for selector in ALL_SELECTORS:
            short = ephem.house_pos(
                63.7, 19.3, GRID_EPS[0], (118.6,), selector.encode()
            )
            full = ephem.house_pos(
                63.7, 19.3, GRID_EPS[0], (118.6, 0.0), selector.encode()
            )
            assert short == full

    def test_an_omitted_selector_selects_placidus(self):
        assert ephem.house_pos(
            63.7, 19.3, GRID_EPS[0], (118.6, 2.0)
        ) == ephem.house_pos(63.7, 19.3, GRID_EPS[0], ord("P"), 118.6, 2.0)

    @pytest.mark.parametrize("selector", list("PKRCBMXAEVWDNOLQTJUIYFS") + ["i"])
    def test_a_body_on_a_cusp_is_answered_with_that_cusp(self, selector):
        for armc, geolat, eps in itertools.product(GRID_ARMC, GRID_LAT, GRID_EPS[:1]):
            cusps, _ = ephem.houses_armc(armc, geolat, eps, selector.encode("latin-1"))
            for number, cusp in enumerate(cusps[:12], start=1):
                got = ephem.house_pos(armc, geolat, eps, ord(selector), cusp, 0.0)
                if selector == "K" and got == 0.0 and number in (4, 10):
                    # The fourth and the tenth cusp are the two edges of the
                    # pair of quadrants the Koch construction spans, so a body
                    # exactly on one of them sits on the boundary of the
                    # placeable band and rounding decides which side it takes.
                    continue
                assert abs(((got - number + 6.0) % 12.0) - 6.0) < 1e-6

    def test_the_answer_stays_in_the_half_open_interval(self):
        for armc, geolat, eps, lon, blat in itertools.product(
            GRID_ARMC[:3], GRID_LAT, GRID_EPS[:1], GRID_LON, GRID_BLAT[:3]
        ):
            for selector in ALL_SELECTORS:
                pos = ephem.house_pos(armc, geolat, eps, ord(selector), lon, blat)
                top = 37.0 if selector in ("G", "g") else 13.0
                assert pos == 0.0 or 1.0 <= pos < top

    def test_the_first_cusp_never_comes_back_as_thirteen(self):
        # The interval is half-open at both ends. A body fed back on the first
        # cusp lands on the seam between the twelfth house and the first, and
        # the answer sits on one side of it or the other -- never at 13.0,
        # which is the same point named from outside the interval.
        # 'H' is left out: near the equator its cusp construction and its
        # position helper anchor the azimuth scale half a circle apart, which
        # belongs to the horizontal system and not to the dispatcher.
        for selector in ALL_SELECTORS.replace("H", ""):
            for armc, geolat in itertools.product(GRID_ARMC, GRID_LAT):
                cusps, _ = ephem.houses_armc(
                    armc, geolat, GRID_EPS[0], selector.encode("latin-1")
                )
                top = 37.0 if selector in ("G", "g") else 13.0
                pos = ephem.house_pos(
                    armc, geolat, GRID_EPS[0], ord(selector), cusps[0], 0.0
                )
                if pos == 0.0:
                    continue
                assert 1.0 <= pos < top
                assert min(pos - 1.0, top - pos) < 1e-6

    def test_placidus_reaches_the_upper_half_exactly_above_the_horizon(self):
        for armc, geolat, eps, lon, blat in itertools.product(
            GRID_ARMC, GRID_LAT + [90.0, -90.0, 89.999], GRID_EPS, GRID_LON, GRID_BLAT
        ):
            dec = _declination(lon, blat, eps)
            if abs(geolat) == 90.0:
                # At a pole the horizon is the celestial equator.
                above = (
                    math.sin(math.radians(geolat)) * math.sin(math.radians(dec)) >= 0.0
                )
            else:
                hour_angle = (armc - _right_ascension(lon, blat, eps)) % 360.0
                above = (
                    math.tan(math.radians(dec)) * math.tan(math.radians(geolat))
                    + math.cos(math.radians(hour_angle))
                ) >= 0.0
            pos = ephem.house_pos(armc, geolat, eps, ord("P"), lon, blat)
            assert (pos >= 7.0) is above


def _declination(lon: float, lat: float, eps: float) -> float:
    """Declination of an ecliptic place, in degrees."""
    sin_dec = math.sin(math.radians(lat)) * math.cos(math.radians(eps)) + math.cos(
        math.radians(lat)
    ) * math.sin(math.radians(lon)) * math.sin(math.radians(eps))
    return math.degrees(math.asin(max(-1.0, min(1.0, sin_dec))))


def _right_ascension(lon: float, lat: float, eps: float) -> float:
    """Right ascension of an ecliptic place, in degrees."""
    return (
        math.degrees(
            math.atan2(
                math.sin(math.radians(lon)) * math.cos(math.radians(eps))
                - math.tan(math.radians(lat)) * math.sin(math.radians(eps)),
                math.cos(math.radians(lon)),
            )
        )
        % 360.0
    )
