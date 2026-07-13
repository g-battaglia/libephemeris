"""
LEB vs Skyfield Comparison: Solar Eclipses.

Validates solar eclipse search and circumstances.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem

from .conftest import TOLS, CompareHelper, year_to_jd


TEST_LOCATIONS = [
    ("Rome", 41.9, 12.5, 0),
    ("NewYork", 40.7, -74.0, 0),
    ("Sydney", -33.9, 151.2, 0),
    ("Equator", 0.0, 0.0, 0),
]


# geopos for sol_eclipse_when_loc is (lon, lat, alt), not (lat, lon, alt)
def _geopos(lat: float, lon: float, alt: float) -> tuple[float, float, float]:
    """Convert (lat, lon, alt) parametrize order to (lon, lat, alt) geopos."""
    return (lon, lat, alt)


class TestSolarEclipseGlobal:
    """Global solar eclipse search."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_eclipse_search_2024_2028(self, compare: CompareHelper):
        """Solar eclipse global timing matches within tolerance."""
        jd_start = year_to_jd(2024)

        ref_result = compare.skyfield(
            ephem.sol_eclipse_when_glob, jd_start, 2, 0, False
        )
        leb_result = compare.leb(ephem.sol_eclipse_when_glob, jd_start, 2, 0, False)

        # Check first eclipse found
        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_jd = ref_result[1][0]
            leb_jd = leb_result[1][0]
            diff_sec = abs(ref_jd - leb_jd) * 86400.0

            assert diff_sec < TOLS.ECLIPSE_TIMING_SEC, (
                f"Global eclipse timing diff = {diff_sec:.1f}s"
            )


class TestSolarEclipseLocal:
    """Local solar eclipse circumstances."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    def test_local_eclipse(
        self, compare: CompareHelper, name: str, lat: float, lon: float, alt: float
    ):
        """Local eclipse timing matches within tolerance."""
        jd_start = year_to_jd(2024)
        geopos = _geopos(lat, lon, alt)

        ref_result = compare.skyfield(
            ephem.sol_eclipse_when_loc, jd_start, geopos, 0, False
        )
        leb_result = compare.leb(ephem.sol_eclipse_when_loc, jd_start, geopos, 0, False)

        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_jd = ref_result[1][0]
            leb_jd = leb_result[1][0]
            diff_sec = abs(ref_jd - leb_jd) * 86400.0

            assert diff_sec < TOLS.ECLIPSE_TIMING_SEC, (
                f"Local eclipse at {name}: timing diff = {diff_sec:.1f}s"
            )


class TestSolarEclipseWhere:
    """Solar eclipse central line."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_eclipse_where(self, compare: CompareHelper):
        """Eclipse central line position matches."""
        jd = year_to_jd(2024) + 60  # Around April 2024 eclipse

        ref_result = compare.skyfield(ephem.sol_eclipse_where, jd, 0)
        leb_result = compare.leb(ephem.sol_eclipse_where, jd, 0)

        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_lat = ref_result[1][0]
            ref_lon = ref_result[1][1]
            leb_lat = leb_result[1][0]
            leb_lon = leb_result[1][1]

            lat_diff = abs(ref_lat - leb_lat)
            lon_diff = abs(ref_lon - leb_lon)

            assert lat_diff < TOLS.ECLIPSE_POSITION_DEG, f"Lat diff = {lat_diff:.4f}°"
            assert lon_diff < TOLS.ECLIPSE_POSITION_DEG, f"Lon diff = {lon_diff:.4f}°"


class TestSolarEclipseHow:
    """Solar eclipse circumstances at location."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS[:2])
    def test_eclipse_how(
        self, compare: CompareHelper, name: str, lat: float, lon: float, alt: float
    ):
        """Eclipse magnitude matches within tolerance."""
        jd = year_to_jd(2024) + 60  # Around April 2024 eclipse
        geopos = _geopos(lat, lon, alt)

        ref_result = compare.skyfield(ephem.sol_eclipse_how, jd, geopos)
        leb_result = compare.leb(ephem.sol_eclipse_how, jd, geopos)

        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_mag = ref_result[1][0]
            leb_mag = leb_result[1][0]

            mag_diff = abs(ref_mag - leb_mag)
            assert mag_diff < TOLS.ECLIPSE_MAGNITUDE, (
                f"Eclipse magnitude at {name}: diff = {mag_diff:.4f}"
            )
