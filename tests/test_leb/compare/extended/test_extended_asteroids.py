"""
LEB vs Skyfield Comparison: Asteroid Positions (Extended Tier).

Validates ecliptic longitude, latitude, distance, and speed for all 5
asteroid bodies: Chiron, Ceres, Pallas, Juno, Vesta.

Asteroid SPK coverage is limited to ~1900-2100 CE.  Test dates are filtered
to this range to avoid catastrophically wrong Keplerian fallback data.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED

from tests.test_leb.compare.conftest import (
    ASTEROID_BODIES,
    CompareHelper,
    filter_asteroid_dates,
    lon_error_arcsec,
    generate_test_dates,
    year_to_jd,
)

from .conftest import TOLS_EXT

_SPK_START = year_to_jd(1910)
_SPK_END = year_to_jd(2090)


@pytest.fixture(scope="session")
def asteroid_dates() -> list[float]:
    return generate_test_dates(100, _SPK_START, _SPK_END)


class TestExtAsteroidPrecision:
    """All-in-one precision check for asteroids within SPK coverage."""

    @pytest.mark.leb_compare_extended
    @pytest.mark.slow
    @pytest.mark.parametrize("body_id,body_name", ASTEROID_BODIES)
    def test_all_components(
        self,
        compare: CompareHelper,
        asteroid_dates: list[float],
        body_id: int,
        body_name: str,
    ):
        dates = filter_asteroid_dates(asteroid_dates, body_id)
        if not dates:
            pytest.skip(f"No valid dates for {body_name}")

        tol_lon = TOLS_EXT.ASTEROID_ARCSEC
        tol_lat = TOLS_EXT.ASTEROID_ARCSEC
        tol_dist = TOLS_EXT.DISTANCE_AU
        tol_speed_lon = TOLS_EXT.ASTEROID_SPEED_LON_DEG_DAY
        tol_speed_lat = TOLS_EXT.ASTEROID_SPEED_LAT_DEG_DAY

        worst = {
            "lon": (0.0, 0.0),
            "lat": (0.0, 0.0),
            "dist": (0.0, 0.0),
            "speed_lon": (0.0, 0.0),
            "speed_lat": (0.0, 0.0),
        }

        for jd in dates:
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

        e_lon, jd_lon = worst["lon"]
        assert e_lon < tol_lon, (
            f'{body_name}: max lon error = {e_lon:.4f}" at JD {jd_lon:.1f}'
        )

        e_lat, jd_lat = worst["lat"]
        assert e_lat < tol_lat, (
            f'{body_name}: max lat error = {e_lat:.4f}" at JD {jd_lat:.1f}'
        )

        e_dist, jd_dist = worst["dist"]
        assert e_dist < tol_dist, (
            f"{body_name}: max dist error = {e_dist:.2e} AU at JD {jd_dist:.1f}"
        )

        e_slon, jd_slon = worst["speed_lon"]
        assert e_slon < tol_speed_lon, (
            f"{body_name}: max speed error = {e_slon:.6f} deg/day at JD {jd_slon:.1f}"
        )

        e_slat, jd_slat = worst["speed_lat"]
        assert e_slat < tol_speed_lat, (
            f"{body_name}: max lat speed error = {e_slat:.6f} deg/day"
        )
