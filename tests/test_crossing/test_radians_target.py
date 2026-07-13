"""Unit-invariance tests for crossing targets expressed with FLG_RADIANS."""

from __future__ import annotations

import math

import pytest

from libephemeris import crossing
from libephemeris.constants import FLG_RADIANS, FLG_SWIEPH, MARS

JD = 2451545.0


@pytest.mark.unit
class TestCrossingRadiansTarget:
    def test_solcross_radians_equals_degrees(self):
        deg = crossing.solcross_ut(90.0, JD, flags=FLG_SWIEPH)
        rad = crossing.solcross_ut(
            math.radians(90.0), JD, flags=FLG_SWIEPH | FLG_RADIANS
        )
        assert abs(rad - deg) < 1e-6

    def test_solcross_tt_radians_equals_degrees(self):
        deg = crossing.solcross(90.0, JD, flags=FLG_SWIEPH)
        rad = crossing.solcross(math.radians(90.0), JD, flags=FLG_SWIEPH | FLG_RADIANS)
        assert abs(rad - deg) < 1e-6

    def test_mooncross_radians_equals_degrees(self):
        deg = crossing.mooncross_ut(200.0, JD, flags=FLG_SWIEPH)
        rad = crossing.mooncross_ut(
            math.radians(200.0), JD, flags=FLG_SWIEPH | FLG_RADIANS
        )
        assert abs(rad - deg) < 1e-6

    def test_mooncross_tt_radians_equals_degrees(self):
        deg = crossing.mooncross(200.0, JD, flags=FLG_SWIEPH)
        rad = crossing.mooncross(
            math.radians(200.0), JD, flags=FLG_SWIEPH | FLG_RADIANS
        )
        assert abs(rad - deg) < 1e-6

    def test_helio_cross_radians_equals_degrees(self):
        # Coherent unit conversion applies to the heliocentric family too.
        deg = crossing.helio_cross_ut(MARS, 120.0, JD, FLG_SWIEPH)
        rad = crossing.helio_cross_ut(
            MARS, math.radians(120.0), JD, FLG_SWIEPH | FLG_RADIANS
        )
        assert abs(rad - deg) < 1e-6

    def test_cross_ut_radians_equals_degrees(self):
        deg = crossing.cross_ut(MARS, 45.0, JD, FLG_SWIEPH)
        rad = crossing.cross_ut(MARS, math.radians(45.0), JD, FLG_SWIEPH | FLG_RADIANS)
        assert abs(rad - deg) < 1e-6
