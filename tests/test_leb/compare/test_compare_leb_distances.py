from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    FLG_SPEED,
    FLG_HELCTR,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
)

from .conftest import (
    TOLS,
    ICRS_PLANETS,
    ASTEROID_BODIES,
    CompareHelper,
    filter_asteroid_dates,
)

DISTANCE_BODIES = ICRS_PLANETS + [(15, "Chiron"), (17, "Ceres")]
HELIO_BODIES = [
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
]


class TestGeocentricDistance:
    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", DISTANCE_BODIES)
    def test_geocentric_distance(
        self,
        compare: CompareHelper,
        test_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        max_err = 0.0
        worst_jd = 0.0

        dates = filter_asteroid_dates(test_dates_100, body_id)
        for jd in dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


class TestHeliocentricDistance:
    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", HELIO_BODIES)
    def test_heliocentric_distance(
        self,
        compare: CompareHelper,
        test_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        flags = FLG_SPEED | FLG_HELCTR
        max_err = 0.0
        worst_jd = 0.0

        for jd in test_dates_100:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS.DISTANCE_AU, (
            f"{body_name} helio: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
