"""Crossing functions honor FLG_RADIANS on the target longitude.

Measured against the reference oracle: solcross_ut(radians(90), jd,
FLG_RADIANS) returns the same instant as the plain 90-degree call
(+171.5748 d from J2000), and mooncross the same (+25.4152 d). An earlier
revision left the target in degrees while calc_ut returned radians,
corrupting the Newton/wrap math (solcross drifted to +19894 d, a crossing
~54 years out). The reference's helio family mishandles the flag on its own
side (divergent search); libephemeris keeps the coherent radians reading
there too — see known-differences.md §11.
"""

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
        # Oracle: both +171.5748 d from J2000
        assert abs(rad - deg) < 1e-6
        assert abs((deg - JD) - 171.5748) < 0.01

    def test_solcross_tt_radians_equals_degrees(self):
        deg = crossing.solcross(90.0, JD, flags=FLG_SWIEPH)
        rad = crossing.solcross(math.radians(90.0), JD, flags=FLG_SWIEPH | FLG_RADIANS)
        assert abs(rad - deg) < 1e-6

    def test_mooncross_radians_equals_degrees(self):
        deg = crossing.mooncross_ut(200.0, JD, flags=FLG_SWIEPH)
        rad = crossing.mooncross_ut(
            math.radians(200.0), JD, flags=FLG_SWIEPH | FLG_RADIANS
        )
        # Oracle: both +25.4152 d from J2000
        assert abs(rad - deg) < 1e-6
        assert abs((deg - JD) - 25.4152) < 0.01

    def test_mooncross_tt_radians_equals_degrees(self):
        deg = crossing.mooncross(200.0, JD, flags=FLG_SWIEPH)
        rad = crossing.mooncross(
            math.radians(200.0), JD, flags=FLG_SWIEPH | FLG_RADIANS
        )
        assert abs(rad - deg) < 1e-6

    def test_helio_cross_radians_equals_degrees(self):
        # Coherent radians reading (the reference's own helio family
        # mishandles the flag — documented divergence, kd §11).
        deg = crossing.helio_cross_ut(MARS, 120.0, JD, FLG_SWIEPH)
        rad = crossing.helio_cross_ut(
            MARS, math.radians(120.0), JD, FLG_SWIEPH | FLG_RADIANS
        )
        assert abs(rad - deg) < 1e-6

    def test_cross_ut_radians_equals_degrees(self):
        deg = crossing.cross_ut(MARS, 45.0, JD, FLG_SWIEPH)
        rad = crossing.cross_ut(MARS, math.radians(45.0), JD, FLG_SWIEPH | FLG_RADIANS)
        assert abs(rad - deg) < 1e-6
