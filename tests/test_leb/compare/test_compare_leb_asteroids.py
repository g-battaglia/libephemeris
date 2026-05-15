from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED

from .conftest import (
    TOLS,
    ASTEROID_BODIES,
    CompareHelper,
    filter_asteroid_dates,
    lon_error_arcsec,
    _ASTEROID_BODY_IDS,
)


class TestAsteroidPrecision:
    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_all_components(
        self,
        compare: CompareHelper,
        test_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        max_lon_err = 0.0
        max_lat_err = 0.0
        max_dist_err = 0.0
        max_speed_lon_err = 0.0

        dates = filter_asteroid_dates(test_dates_200, body_id)
        for jd in dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0
            dist_err = abs(ref[2] - leb[2])
            speed_lon_err = abs(ref[3] - leb[3])

            max_lon_err = max(max_lon_err, lon_err)
            max_lat_err = max(max_lat_err, lat_err)
            max_dist_err = max(max_dist_err, dist_err)
            max_speed_lon_err = max(max_speed_lon_err, speed_lon_err)

        assert max_lon_err < TOLS.ASTEROID_ARCSEC, (
            f'{body_name}: lon error {max_lon_err:.4f}"'
        )
        assert max_lat_err < TOLS.ASTEROID_ARCSEC, (
            f'{body_name}: lat error {max_lat_err:.4f}"'
        )
        assert max_dist_err < TOLS.DISTANCE_AU, (
            f"{body_name}: dist error {max_dist_err:.2e} AU"
        )
        tol_speed = (
            TOLS.ASTEROID_SPEED_LON_DEG_DAY
            if body_id in _ASTEROID_BODY_IDS
            else TOLS.SPEED_LON_DEG_DAY
        )
        assert max_speed_lon_err < tol_speed, (
            f"{body_name}: speed error {max_speed_lon_err:.6f} deg/day"
        )
