"""
LEB vs Skyfield Comparison: Lunar/Ecliptic Bodies (Extended Tier).

Validates all 6 Pipeline B (ecliptic-direct) bodies with per-body tolerances
across the extended tier range (-5000 to 5000).
"""

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

from .conftest import TOLS_EXT


class TestExtLunarPrecision:
    """All-in-one precision check for lunar/ecliptic bodies."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ECLIPTIC_BODIES)
    def test_all_components(
        self,
        compare: CompareHelper,
        ext_dates_200: list[float],
        body_id: int,
        body_name: str,
    ):
        tol_lon = max(
            ECLIPTIC_TOLERANCES.get(body_id, {}).get("lon", TOLS_EXT.ECLIPTIC_ARCSEC),
            TOLS_EXT.ECLIPTIC_ARCSEC,
        )
        tol_lat = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "lat", TOLS_EXT.ECLIPTIC_ARCSEC
        )
        tol_dist = 0.01
        tol_speed_lon = max(
            ECLIPTIC_TOLERANCES.get(body_id, {}).get(
                "speed", TOLS_EXT.SPEED_LON_DEG_DAY
            ),
            TOLS_EXT.SPEED_LON_DEG_DAY,
        )
        tol_speed_lat = ECLIPTIC_TOLERANCES.get(body_id, {}).get(
            "speed", TOLS_EXT.SPEED_LAT_DEG_DAY
        )
        tol_speed_dist = TOLS_EXT.SPEED_DIST_AU_DAY

        worst = {
            "lon": (0.0, 0.0),
            "lat": (0.0, 0.0),
            "dist": (0.0, 0.0),
            "speed_lon": (0.0, 0.0),
            "speed_lat": (0.0, 0.0),
            "speed_dist": (0.0, 0.0),
        }

        for jd in ext_dates_200:
            ref, _ = compare.skyfield(ephem.calc_ut, jd, body_id, FLG_SPEED)
            leb, _ = compare.leb(ephem.calc_ut, jd, body_id, FLG_SPEED)

            err_lon = lon_error_arcsec(ref[0], leb[0])
            if err_lon > worst["lon"][0]:
                worst["lon"] = (err_lon, jd)

            err_lat = abs(ref[1] - leb[1]) * 3600.0
            if err_lat > worst["lat"][0]:
                worst["lat"] = (err_lat, jd)

            err_dist = abs(ref[2] - leb[2])
            if err_dist > worst["dist"][0]:
                worst["dist"] = (err_dist, jd)

            err_speed_lon = abs(ref[3] - leb[3])
            if err_speed_lon > worst["speed_lon"][0]:
                worst["speed_lon"] = (err_speed_lon, jd)

            err_speed_lat = abs(ref[4] - leb[4])
            if err_speed_lat > worst["speed_lat"][0]:
                worst["speed_lat"] = (err_speed_lat, jd)

            err_speed_dist = abs(ref[5] - leb[5])
            if err_speed_dist > worst["speed_dist"][0]:
                worst["speed_dist"] = (err_speed_dist, jd)

        e_lon, jd_lon = worst["lon"]
        assert e_lon < tol_lon, (
            f'{body_name}: max lon error = {e_lon:.4f}" at JD {jd_lon:.1f} '
            f'(tol={tol_lon}")'
        )

        e_lat, jd_lat = worst["lat"]
        assert e_lat < tol_lat, (
            f'{body_name}: max lat error = {e_lat:.4f}" at JD {jd_lat:.1f}'
        )

        e_dist, jd_dist = worst["dist"]
        assert e_dist < tol_dist, (
            f"{body_name}: max distance error = {e_dist:.2e} AU at JD {jd_dist:.1f}"
        )

        e_slon, jd_slon = worst["speed_lon"]
        assert e_slon < tol_speed_lon, (
            f"{body_name}: max speed error = {e_slon:.6f} deg/day at JD {jd_slon:.1f} "
            f"(tol={tol_speed_lon})"
        )

        e_slat, jd_slat = worst["speed_lat"]
        assert e_slat < tol_speed_lat, (
            f"{body_name}: max lat speed error = {e_slat:.6f} deg/day "
            f"at JD {jd_slat:.1f}"
        )

        e_sdist, jd_sdist = worst["speed_dist"]
        assert e_sdist < tol_speed_dist, (
            f"{body_name}: max dist speed error = {e_sdist:.2e} AU/day "
            f"at JD {jd_sdist:.1f}"
        )
