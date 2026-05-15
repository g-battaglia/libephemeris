from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MOON,
    MARS,
    JUPITER,
    FLG_SPEED,
    FLG_EQUATORIAL,
    FLG_J2000,
    FLG_HELCTR,
    FLG_BARYCTR,
    FLG_TRUEPOS,
    FLG_NOABERR,
)

from tests.test_leb.compare.conftest import (
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_BASE

FLAG_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]

FLAG_COMBINATIONS = [
    (FLG_SPEED, "default"),
    (FLG_SPEED | FLG_EQUATORIAL, "equatorial"),
    (FLG_SPEED | FLG_J2000, "J2000"),
    (FLG_SPEED | FLG_EQUATORIAL | FLG_J2000, "equatorial_J2000"),
    (FLG_SPEED | FLG_HELCTR, "heliocentric"),
    (FLG_SPEED | FLG_BARYCTR, "barycentric"),
    (FLG_SPEED | FLG_TRUEPOS, "truepos"),
    (FLG_SPEED | FLG_NOABERR, "noaberr"),
]


class TestBaseFlagCombinations:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", FLAG_BODIES)
    @pytest.mark.parametrize("flags,flag_name", FLAG_COMBINATIONS)
    def test_flags(
        self,
        compare: CompareHelper,
        base_dates_100: list[float],
        body_id: int,
        body_name: str,
        flags: int,
        flag_name: str,
    ):
        max_err = 0.0
        max_speed_err = 0.0
        worst_jd = 0.0

        for jd in base_dates_100:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0
            err = max(lon_err, lat_err)
            speed_err = abs(ref[3] - leb[3])

            if err > max_err:
                max_err = err
                worst_jd = jd
            max_speed_err = max(max_speed_err, speed_err)

        assert max_err < TOLS_BASE.EQUATORIAL_ARCSEC, (
            f'{body_name} {flag_name}: max error = {max_err:.4f}" at JD {worst_jd:.1f}'
        )
        assert max_speed_err < TOLS_BASE.SPEED_LON_DEG_DAY, (
            f"{body_name} {flag_name}: max speed error = {max_speed_err:.6f} deg/day"
        )
