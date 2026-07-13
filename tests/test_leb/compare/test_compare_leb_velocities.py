from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED

from .conftest import (
    TOLS,
    ICRS_PLANETS,
    ECLIPTIC_BODIES,
    ECLIPTIC_TOLERANCES,
    ASTEROID_BODIES,
    CompareHelper,
    filter_asteroid_dates,
    _ASTEROID_BODY_IDS,
)

ALL_LEB_BODIES = ICRS_PLANETS + ECLIPTIC_BODIES + ASTEROID_BODIES


class TestVelocityEcliptic:
    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ALL_LEB_BODIES)
    def test_ecliptic_velocity(
        self,
        compare: CompareHelper,
        test_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        max_lon_err = 0.0
        max_lat_err = 0.0

        dates = filter_asteroid_dates(test_dates_100, body_id)
        for jd in dates:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            lon_err = abs(ref[3] - leb[3])
            lat_err = abs(ref[4] - leb[4])

            max_lon_err = max(max_lon_err, lon_err)
            max_lat_err = max(max_lat_err, lat_err)

        ecl_tol = ECLIPTIC_TOLERANCES.get(body_id, {}).get("speed")
        if ecl_tol is not None:
            tol_lon = ecl_tol
        elif body_id in _ASTEROID_BODY_IDS:
            tol_lon = TOLS.ASTEROID_SPEED_LON_DEG_DAY
        else:
            tol_lon = TOLS.SPEED_LON_DEG_DAY

        if ecl_tol is not None:
            tol_lat = ecl_tol
        elif body_id in _ASTEROID_BODY_IDS:
            tol_lat = TOLS.ASTEROID_SPEED_LAT_DEG_DAY
        else:
            tol_lat = TOLS.SPEED_LAT_DEG_DAY

        assert max_lon_err < tol_lon, (
            f"{body_name}: lon speed error {max_lon_err:.6f} deg/day (tol={tol_lon})"
        )
        assert max_lat_err < tol_lat, (
            f"{body_name}: lat speed error {max_lat_err:.6f} deg/day (tol={tol_lat})"
        )


class TestDistanceVelocity:
    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS + ASTEROID_BODIES)
    def test_speed_distance(
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

            err = abs(ref[5] - leb[5])
            if err > max_err:
                max_err = err
                worst_jd = jd

        tol = (
            TOLS.ASTEROID_SPEED_DIST_AU_DAY
            if body_id in _ASTEROID_BODY_IDS
            else TOLS.SPEED_DIST_AU_DAY
        )
        assert max_err < tol, (
            f"{body_name}: max dist speed error = {max_err:.2e} AU/day "
            f"at JD {worst_jd:.1f}"
        )


class TestVelocityEquatorial:
    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_equatorial_velocity(
        self,
        compare: CompareHelper,
        test_dates_50: list[float],
        body_id: int,
        body_name: str,
    ):
        from libephemeris.constants import FLG_EQUATORIAL

        flags = FLG_SPEED | FLG_EQUATORIAL
        max_err = 0.0

        for jd in test_dates_50:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS.SPEED_LON_DEG_DAY, (
            f"{body_name}: max equatorial speed error = {max_err:.6f} deg/day"
        )
