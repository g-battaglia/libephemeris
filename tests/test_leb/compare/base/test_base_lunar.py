from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED

from tests.test_leb.compare.conftest import (
    ECLIPTIC_BODIES,
    ECLIPTIC_TOLERANCES,
    CompareHelper,
    lon_error_arcsec,
)

from .conftest import TOLS_BASE


class TestBaseLunarPrecision:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_all_components(
        self,
        compare: CompareHelper,
        base_dates_300: list[float],
        body_id: int,
        body_name: str,
    ):
        tol_lon = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "lon", TOLS_BASE.ECLIPTIC_ARCSEC
        )
        tol_lat = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "lat", TOLS_BASE.ECLIPTIC_ARCSEC
        )
        tol_dist = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "dist", TOLS_BASE.DISTANCE_AU
        )
        tol_speed = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "speed", TOLS_BASE.SPEED_LON_DEG_DAY
        )

        max_lon_err = 0.0
        max_lat_err = 0.0
        max_dist_err = 0.0
        max_speed_err = 0.0
        worst_lon_jd = 0.0
        worst_lat_jd = 0.0
        worst_dist_jd = 0.0
        worst_speed_jd = 0.0

        for jd in base_dates_300:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

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

        assert max_lon_err < tol_lon, (
            f'{body_name}: max lon error = {max_lon_err:.4f}" at JD {worst_lon_jd:.1f} '
            f'(tol={tol_lon}")'
        )
        assert max_lat_err < tol_lat, (
            f'{body_name}: max lat error = {max_lat_err:.4f}" at JD {worst_lat_jd:.1f}'
        )
        assert max_dist_err < tol_dist, (
            f"{body_name}: max dist error = {max_dist_err:.2e} AU at JD {worst_dist_jd:.1f}"
        )
        assert max_speed_err < tol_speed, (
            f"{body_name}: max speed error = {max_speed_err:.6f} deg/day at JD {worst_speed_jd:.1f} "
            f"(tol={tol_speed})"
        )
