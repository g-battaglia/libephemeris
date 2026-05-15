from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED
from libephemeris.exceptions import EphemerisRangeError

from tests.test_leb.compare.conftest import (
    ASTEROID_BODIES,
    CompareHelper,
    filter_asteroid_dates,
    lon_error_arcsec,
)

from .conftest import TOLS_BASE


class TestBaseAsteroidPrecision:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_all_components(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        dates = filter_asteroid_dates(base_dates_300, body_id)
        assert len(dates) > 0, f"{body_name}: no dates in SPK coverage range"

        max_lon_err = 0.0
        max_lat_err = 0.0
        max_dist_err = 0.0
        max_speed_err = 0.0
        worst_lon_jd = 0.0
        worst_lat_jd = 0.0
        worst_dist_jd = 0.0
        worst_speed_jd = 0.0

        for jd in dates:
            try:
                ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
                leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)
            except (KeyError, ValueError, EphemerisRangeError):
                continue

            lon_err = lon_error_arcsec(ref[0], leb[0])
            lat_err = abs(ref[1] - leb[1]) * 3600.0
            dist_err = abs(ref[2] - leb[2])
            speed_err = abs(ref[3] - leb[3])

            if lon_err > max_lon_err:
                max_lon_err = lon_err
                worst_lon_jd = jd
            if lat_err > max_lat_err:
                max_lat_err = lat_err
                worst_lat_jd = jd
            if dist_err > max_dist_err:
                max_dist_err = dist_err
                worst_dist_jd = jd
            if speed_err > max_speed_err:
                max_speed_err = speed_err
                worst_speed_jd = jd

        assert max_lon_err < TOLS_BASE.ASTEROID_ARCSEC, (
            f'{body_name}: max lon error = {max_lon_err:.4f}" at JD {worst_lon_jd:.1f}'
        )
        assert max_lat_err < TOLS_BASE.ASTEROID_ARCSEC, (
            f'{body_name}: max lat error = {max_lat_err:.4f}"'
        )
        assert max_dist_err < TOLS_BASE.DISTANCE_AU, (
            f"{body_name}: max dist error = {max_dist_err:.2e} AU at JD {worst_dist_jd:.1f}"
        )
        assert max_speed_err < TOLS_BASE.SPEED_LON_DEG_DAY, (
            f"{body_name}: max speed error = {max_speed_err:.6f} deg/day at JD {worst_speed_jd:.1f}"
        )
