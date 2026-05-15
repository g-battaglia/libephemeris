from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.exceptions import EphemerisRangeError
from libephemeris.constants import (
    FLG_SPEED,
    FLG_HELCTR,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
)

from tests.test_leb.compare.conftest import (
    ICRS_PLANETS,
    CompareHelper,
    filter_asteroid_dates,
)

from .conftest import TOLS_BASE

DISTANCE_BODIES = ICRS_PLANETS + [(15, "Chiron"), (17, "Ceres")]
HELIO_BODIES = [
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
]


class TestBaseGeocentricDistance:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", DISTANCE_BODIES)
    def test_geocentric_distance(
        self,
        compare: CompareHelper,
        base_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        dates = filter_asteroid_dates(base_dates_150, body_id)
        max_err = 0.0
        worst_jd = 0.0

        for jd in dates:
            try:
                ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
                leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)
            except (KeyError, ValueError, EphemerisRangeError):
                continue

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )


class TestBaseHeliocentricDistance:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", HELIO_BODIES)
    def test_heliocentric_distance(
        self,
        compare: CompareHelper,
        base_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        flags = FLG_SPEED | FLG_HELCTR
        max_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_150:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[2] - leb[2])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.DISTANCE_AU, (
            f"{body_name} helio: max dist error = {max_err:.2e} AU at JD {worst_jd:.1f}"
        )
