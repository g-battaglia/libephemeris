# SPDX-License-Identifier: AGPL-3.0-only
# Copyright (c) 2025-2026 Giacomo Battaglia
"""Guards on the Sunshine and APC house positions (``I``, ``i`` and ``Y``).

The three selectors divide one pencil of house circles -- the great circles
through the north and south points of the horizon -- and differ only in the
parallel of declination whose semi-arcs they cut. :func:`house_pos` is handed
no Sun, so the Sunshine reference parallel is the celestial equator and the
two Sunshine selectors answer the Regiomontanus question; the APC keeps its
own reference, the Ascendant's parallel.

Two things are pinned here that nothing else pins.

*The refusal at a geographic pole.* There the pencil has no anchors of its
own: the horizon is the celestial equator and every house circle crosses it
in the same pair of points, so the construction has nothing left to divide.
All three selectors refuse, with the same exception, and they refuse only
there -- a thousandth of a degree inside the pole they answer.

*The obliquity axis.* The recorded frames span a degree and a half of
obliquity around the obliquity of date, which is no test at all for the APC,
whose reference parallel is a function of the obliquity through the
Ascendant. The identities below are asserted far away from it.
"""

from __future__ import annotations

import itertools
import math

import pytest

import libephemeris as ephem

#: The frame the recorded polar grid uses.
POLAR_ARMC = 285.2170818087561
POLAR_EPS = 23.45654434517177

#: A body the recorded family carries at every latitude.
BODY_LON = 285.38339935017507
BODY_LAT = 0.00014133248984594553

#: Obliquities far from the 23.4 degrees the recorded family spans.
FAR_EPS = [5.0, 60.0]

SELECTORS = ["I", "i", "Y"]


def _circular(left: float, right: float) -> float:
    """Distance between two house positions across the wrap at twelve."""
    gap = abs(left - right)
    return min(gap, abs(gap - 12.0))


@pytest.mark.unit
class TestPolarRefusal:
    """The pencil has no anchors at a geographic pole, and none of the three
    selectors answers there."""

    @pytest.mark.parametrize("selector", SELECTORS)
    @pytest.mark.parametrize("geolat", [90.0, -90.0])
    def test_a_geographic_pole_is_refused(self, selector, geolat):
        with pytest.raises(ZeroDivisionError) as raised:
            ephem.house_pos(POLAR_ARMC, geolat, POLAR_EPS, selector, BODY_LON, BODY_LAT)
        assert str(raised.value) == "float division by zero"

    @pytest.mark.parametrize("selector", SELECTORS)
    def test_the_five_argument_form_refuses_alike(self, selector):
        with pytest.raises(ZeroDivisionError):
            ephem.house_pos(
                POLAR_ARMC,
                90.0,
                POLAR_EPS,
                (BODY_LON, BODY_LAT),
                selector.encode("ascii"),
            )

    @pytest.mark.parametrize("selector", SELECTORS)
    @pytest.mark.parametrize("geolat", [89.999, -89.999, 89.0, 0.0, -66.56])
    def test_every_other_latitude_answers(self, selector, geolat):
        pos = ephem.house_pos(
            POLAR_ARMC, geolat, POLAR_EPS, selector, BODY_LON, BODY_LAT
        )
        assert math.isfinite(pos)
        assert 1.0 <= pos < 13.0

    @pytest.mark.parametrize("selector", SELECTORS)
    def test_the_answer_approaches_a_meridian_cusp_at_the_pole(self, selector):
        """What the refusal stands in for: the limit taken from inside."""
        previous = None
        for geolat in (89.9, 89.99, 89.999, 89.9999):
            pos = ephem.house_pos(
                POLAR_ARMC, geolat, POLAR_EPS, selector, BODY_LON, BODY_LAT
            )
            gap = abs(pos - 4.0)
            assert previous is None or gap < previous
            previous = gap
        assert previous < 1e-4


@pytest.mark.unit
class TestAtAnObliquityFarFromTheObliquityOfDate:
    """The identities that hold by construction, asserted where the recorded
    family carries no frames."""

    @pytest.mark.parametrize("eps", FAR_EPS)
    def test_sunshine_is_the_regiomontanus_position(self, eps):
        """No Sun reaches this entry point, so the reference parallel is the
        equator and the correction vanishes identically."""
        for armc, geolat, lon, blat in itertools.product(
            (0.0, 73.0, 197.5428262395338, 291.0),
            (-80.0, -45.0, -0.001, 0.0, 23.44, 66.56, 89.999),
            (0.0, 47.5, 130.0, 285.38339935017507),
            (-17.5, 0.0, 12.0),
        ):
            regiomontanus = ephem.house_pos(armc, geolat, eps, "R", lon, blat)
            for selector in ("I", "i"):
                assert (
                    ephem.house_pos(armc, geolat, eps, selector, lon, blat)
                    == regiomontanus
                )

    @pytest.mark.parametrize("eps", FAR_EPS)
    def test_the_alternative_solution_answers_like_the_first(self, eps):
        """``i`` is a system of its own on the cusp side; on the position side
        both solutions are handed the same reference parallel."""
        for armc, geolat in itertools.product(
            (0.0, 73.0, 197.5428262395338, 291.0), (-70.0, -10.0, 45.0, 80.0)
        ):
            first = ephem.house_pos(armc, geolat, eps, "I", BODY_LON, BODY_LAT)
            second = ephem.house_pos(armc, geolat, eps, "i", BODY_LON, BODY_LAT)
            assert first == second

    @pytest.mark.parametrize("eps", FAR_EPS)
    def test_apc_meets_sunshine_where_the_ascendant_is_on_the_equator(self, eps):
        """At three quarters of a turn of the ARMC the Ascendant is an
        equinox, so the APC divides the equator as the Sunshine systems do."""
        for geolat, lon, blat in itertools.product(
            (-80.0, -45.0, 0.0, 23.44, 70.0), (0.0, 47.5, 130.0, 285.0), (-17.5, 0.0)
        ):
            apc = ephem.house_pos(270.0, geolat, eps, "Y", lon, blat)
            sunshine = ephem.house_pos(270.0, geolat, eps, "I", lon, blat)
            assert _circular(apc, sunshine) < 1e-12

    @pytest.mark.parametrize("eps", FAR_EPS)
    def test_apc_divides_the_equator_at_the_equator(self, eps):
        """At latitude zero the house circles are the hour circles, so every
        parallel is met at the same meridian distance and the reference
        parallel cannot matter."""
        for armc, lon, blat in itertools.product(
            (0.0, 73.0, 197.5428262395338, 291.0),
            (0.0, 47.5, 130.0, 285.0),
            (-17.5, 0.0, 12.0),
        ):
            apc = ephem.house_pos(armc, 0.0, eps, "Y", lon, blat)
            meridian = ephem.house_pos(armc, 0.0, eps, "X", lon, blat)
            assert _circular(apc, meridian) < 1e-12

    @pytest.mark.parametrize("eps", FAR_EPS)
    @pytest.mark.parametrize("selector", SELECTORS)
    def test_the_cusps_of_a_frame_come_back_as_their_own_numbers(self, selector, eps):
        """The horizon and the meridian are two members of the pencil, so the
        first cusp answers 1 and the tenth answers 10; the eight others are
        the division points the system defines.

        The latitudes stay outside the circle of perpetual apparition of the
        obliquity in play (``|latitude| < 90 - obliquity``), where the cusps
        of the alternative Sunshine solution take a branch of their own: that
        branch belongs to the cusp side and is not the inverse of any house
        position.
        """
        for armc, geolat in itertools.product((37.0, 197.5428262395338), (-25.0, 12.0)):
            cusps, _ = ephem.houses_armc(armc, geolat, eps, selector)
            for number, cusp in enumerate(cusps, start=1):
                pos = ephem.house_pos(armc, geolat, eps, selector, cusp, 0.0)
                assert _circular(pos, float(number)) < 1e-9

    @pytest.mark.parametrize("eps", FAR_EPS)
    @pytest.mark.parametrize("selector", SELECTORS)
    def test_the_upper_half_of_the_wheel_is_the_visible_sky(self, selector, eps):
        """Seven and above is above the horizon, for every reference
        parallel: the horizon is a member of the pencil and carries the whole
        of the parallel's ascensional difference."""
        for armc, geolat, lon, blat in itertools.product(
            (0.0, 73.0, 291.0),
            (-70.0, -23.44, 10.0, 66.56),
            (0.0, 47.5, 130.0, 285.0),
            (-17.5, 0.0, 12.0),
        ):
            dec = math.degrees(
                math.asin(
                    math.sin(math.radians(blat)) * math.cos(math.radians(eps))
                    + math.cos(math.radians(blat))
                    * math.sin(math.radians(lon))
                    * math.sin(math.radians(eps))
                )
            )
            ra = math.degrees(
                math.atan2(
                    math.sin(math.radians(lon)) * math.cos(math.radians(eps))
                    - math.tan(math.radians(blat)) * math.sin(math.radians(eps)),
                    math.cos(math.radians(lon)),
                )
            )
            hour_angle = (armc - ra) % 360.0
            above = (
                math.tan(math.radians(dec)) * math.tan(math.radians(geolat))
                + math.cos(math.radians(hour_angle))
            ) >= 0.0
            pos = ephem.house_pos(armc, geolat, eps, selector, lon, blat)
            assert (pos >= 7.0) is above
