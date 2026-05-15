"""
LEB vs Skyfield Comparison: Lunar Eclipses.

Validates lunar eclipse search and circumstances.
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


class TestLunarEclipseGlobal:
    """Global lunar eclipse search."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    def test_eclipse_search_2024_2028(self, compare: CompareHelper):
        """Lunar eclipse global timing matches within tolerance."""
        jd_start = year_to_jd(2024)

        ref_result = compare.skyfield(ephem.lun_eclipse_when, jd_start, 2, 0)
        leb_result = compare.leb(ephem.lun_eclipse_when, jd_start, 2, 0)

        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_jd = ref_result[1][0]
            leb_jd = leb_result[1][0]
            diff_sec = abs(ref_jd - leb_jd) * 86400.0

            assert diff_sec < TOLS.ECLIPSE_TIMING_SEC, (
                f"Global lunar eclipse timing diff = {diff_sec:.1f}s"
            )


class TestLunarEclipseLocal:
    """Local lunar eclipse circumstances."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    def test_local_eclipse(
        self, compare: CompareHelper, name: str, lat: float, lon: float, alt: float
    ):
        """Local lunar eclipse timing matches within tolerance."""
        jd_start = year_to_jd(2024)

        geopos = (lon, lat, alt)
        ref_result = compare.skyfield(
            ephem.lun_eclipse_when_loc, jd_start, geopos, 2
        )
        leb_result = compare.leb(ephem.lun_eclipse_when_loc, jd_start, geopos, 2)

        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_jd = ref_result[1][0]
            leb_jd = leb_result[1][0]
            diff_sec = abs(ref_jd - leb_jd) * 86400.0

            assert diff_sec < TOLS.ECLIPSE_TIMING_SEC, (
                f"Local lunar eclipse at {name}: timing diff = {diff_sec:.1f}s"
            )


class TestLunarEclipseHow:
    """Lunar eclipse circumstances."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS[:2])
    def test_eclipse_how(
        self, compare: CompareHelper, name: str, lat: float, lon: float, alt: float
    ):
        """Lunar eclipse magnitude matches within tolerance."""
        jd = year_to_jd(2025) + 90  # March 2025 lunar eclipse
        geopos = (lat, lon, alt)

        ref_result = compare.skyfield(ephem.lun_eclipse_how, jd, geopos)
        leb_result = compare.leb(ephem.lun_eclipse_how, jd, geopos)

        if ref_result[0] != 0 and leb_result[0] != 0:
            ref_mag = ref_result[1][0]
            leb_mag = leb_result[1][0]

            mag_diff = abs(ref_mag - leb_mag)
            assert mag_diff < TOLS.ECLIPSE_MAGNITUDE, (
                f"Lunar eclipse magnitude at {name}: diff = {mag_diff:.4f}"
            )
