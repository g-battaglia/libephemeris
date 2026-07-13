"""
LEB vs Skyfield Comparison: Velocity Precision (Extended Tier).

Exhaustive velocity validation for all LEB bodies (excluding asteroids
whose SPK coverage does not span the extended range).
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED, FLG_EQUATORIAL

from tests.test_leb.compare.conftest import (
    ECLIPTIC_BODIES,
    ECLIPTIC_TOLERANCES,
    ICRS_PLANETS,
    CompareHelper,
)

from .conftest import TOLS_EXT

EXT_BODIES = ICRS_PLANETS + ECLIPTIC_BODIES


class TestExtVelocityEcliptic:
    """Longitude and latitude velocity for planets and ecliptic bodies."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", EXT_BODIES)
    def test_ecliptic_velocity(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        ecl_tol = ECLIPTIC_TOLERANCES.get(body_id, {}).get("speed")
        tol_lon = ecl_tol if ecl_tol is not None else TOLS_EXT.SPEED_LON_DEG_DAY
        tol_lat = ecl_tol if ecl_tol is not None else TOLS_EXT.SPEED_LAT_DEG_DAY

        worst_lon = (0.0, 0.0)
        worst_lat = (0.0, 0.0)

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err_lon = abs(ref[3] - leb[3])
            if err_lon > worst_lon[0]:
                worst_lon = (err_lon, jd)

            err_lat = abs(ref[4] - leb[4])
            if err_lat > worst_lat[0]:
                worst_lat = (err_lat, jd)

        e_lon, jd_lon = worst_lon
        assert e_lon < tol_lon, (
            f"{body_name}: max lon speed error = {e_lon:.6f} deg/day at JD {jd_lon:.1f}"
        )

        e_lat, jd_lat = worst_lat
        assert e_lat < tol_lat, (
            f"{body_name}: max lat speed error = {e_lat:.6f} deg/day at JD {jd_lat:.1f}"
        )


class TestExtDistanceVelocity:
    """Distance velocity (result[5]) for ICRS planets."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_speed_distance(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        max_err = 0.0
        worst_jd = 0.0

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err = abs(ref[5] - leb[5])
            if err > max_err:
                max_err = err
                worst_jd = jd

        assert max_err < TOLS_EXT.SPEED_DIST_AU_DAY, (
            f"{body_name}: max dist speed error = {max_err:.2e} AU/day "
            f"at JD {worst_jd:.1f}"
        )


class TestExtVelocityEquatorial:
    """Velocity in equatorial coordinates for ICRS planets."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ICRS_PLANETS)
    def test_equatorial_velocity(
        self,
        compare: CompareHelper,
        ext_dates_100: list[float],
        body_id: int,
        body_name: str,
    ):
        flags = FLG_SPEED | FLG_EQUATORIAL
        max_err = 0.0

        for jd in ext_dates_100:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, flags)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, flags)

            err = abs(ref[3] - leb[3])
            max_err = max(max_err, err)

        assert max_err < TOLS_EXT.SPEED_LON_DEG_DAY, (
            f"{body_name}: max equatorial speed error = {max_err:.6f} deg/day"
        )
