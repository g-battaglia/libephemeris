from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED, FLG_EQUATORIAL
from libephemeris.exceptions import EphemerisRangeError

from tests.test_leb.compare.conftest import (
    ALL_LEB_BODIES,
    ASTEROID_BODIES,
    ECLIPTIC_TOLERANCES,
    ICRS_PLANETS,
    CompareHelper,
    filter_asteroid_dates,
)

from .conftest import TOLS_BASE


class TestBaseVelocityEcliptic:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_LEB_BODIES)
    def test_ecliptic_velocity(
        self,
        compare: CompareHelper,
        base_dates_150: list[float],
        body_id: int,
        body_name: str,
    ):
        dates = filter_asteroid_dates(base_dates_150, body_id)

        ecl_tol = ECLIPTIC_TOLERANCES.get(body_id, {}).get("speed")
        tol_lon = ecl_tol if ecl_tol is not None else TOLS_BASE.SPEED_LON_DEG_DAY
        if ecl_tol is not None:
            tol_lat = ecl_tol
        else:
            asteroid_ids = {b[0] for b in ASTEROID_BODIES}
            tol_lat = (
                TOLS_BASE.ASTEROID_SPEED_LAT_DEG_DAY
                if body_id in asteroid_ids
                else TOLS_BASE.SPEED_LAT_DEG_DAY
            )

        max_lon_err = 0.0
        max_lat_err = 0.0
        worst_lon_jd = 0.0
        worst_lat_jd = 0.0

        for jd in dates:
            try:
                ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
                leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)
            except (KeyError, ValueError, EphemerisRangeError):
                continue

            lon_err = abs(ref[3] - leb[3])
            lat_err = abs(ref[4] - leb[4])

            if lon_err > max_lon_err:
                max_lon_err = lon_err
                worst_lon_jd = jd
            if lat_err > max_lat_err:
                max_lat_err = lat_err
                worst_lat_jd = jd

        assert max_lon_err < tol_lon, (
            f"{body_name}: max lon speed error = {max_lon_err:.6f} deg/day at JD {worst_lon_jd:.1f}"
        )
        assert max_lat_err < tol_lat, (
            f"{body_name}: max lat speed error = {max_lat_err:.6f} deg/day at JD {worst_lat_jd:.1f}"
        )


class TestBaseDistanceVelocity:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS + ASTEROID_BODIES)
    def test_speed_distance(
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

            err = abs(ref[5] - leb[5])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_BASE.SPEED_DIST_AU_DAY, (
            f"{body_name}: max dist speed error = {max_err:.2e} AU/day at JD {worst_jd:.1f}"
        )


class TestBaseVelocityEquatorial:
    @pytest.mark.leb_compare_base
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_equatorial_velocity(
        self,
        compare: CompareHelper,
        base_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        flags = FLG_SPEED | FLG_EQUATORIAL
        max_err = 0.0

        for jd in base_dates_100:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS_BASE.SPEED_LON_DEG_DAY, (
            f"{body_name}: max equatorial speed error = {max_err:.6f} deg/day"
        )
