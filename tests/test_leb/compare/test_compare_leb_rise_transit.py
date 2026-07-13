"""
LEB vs Skyfield Comparison: Rise/Transit/Set.

Validates rise/transit/set timing. These functions use calc_ut
in iterative position refinement loops.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import SUN, MOON, VENUS

from .conftest import TOLS, CompareHelper, year_to_jd


RISE_BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (VENUS, "Venus"),
]

TEST_LOCATIONS = [
    ("Rome", 41.9, 12.5, 0),
    ("NewYork", 40.7, -74.0, 0),
    ("Sydney", -33.9, 151.2, 0),
    ("Equator", 0.0, 0.0, 0),
]

TEST_DATES = [
    (year_to_jd(2024) + 0, "equinox"),
    (year_to_jd(2024) + 91, "solstice"),
    (year_to_jd(2024) + 200, "random"),
]


class TestSunRiseSet:
    """Sun rise/set/transit timing."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    @pytest.mark.parametrize("jd,date_name", TEST_DATES)
    def test_sunrise(
        self,
        compare: CompareHelper,
        name: str,
        lat: float,
        lon: float,
        alt: float,
        jd: float,
        date_name: str,
    ):
        """Sun rise timing matches within tolerance."""
        ref_result = compare.skyfield(ephem.rise_trans, jd, SUN, 1, (lon, lat, alt))
        leb_result = compare.leb(ephem.rise_trans, jd, SUN, 1, (lon, lat, alt))

        if ref_result[0] != 0 and leb_result[0] != 0:
            diff_sec = abs(ref_result[0] - leb_result[0]) * 86400.0
            assert diff_sec < TOLS.RISE_TRANSIT_SEC, (
                f"Sunrise at {name} ({date_name}): diff = {diff_sec:.1f}s"
            )

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    @pytest.mark.parametrize("jd,date_name", TEST_DATES)
    def test_sunset(
        self,
        compare: CompareHelper,
        name: str,
        lat: float,
        lon: float,
        alt: float,
        jd: float,
        date_name: str,
    ):
        """Sun set timing matches within tolerance."""
        ref_result = compare.skyfield(ephem.rise_trans, jd, SUN, 2, (lon, lat, alt))
        leb_result = compare.leb(ephem.rise_trans, jd, SUN, 2, (lon, lat, alt))

        if ref_result[0] != 0 and leb_result[0] != 0:
            diff_sec = abs(ref_result[0] - leb_result[0]) * 86400.0
            assert diff_sec < TOLS.RISE_TRANSIT_SEC, (
                f"Sunset at {name} ({date_name}): diff = {diff_sec:.1f}s"
            )

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    @pytest.mark.parametrize("jd,date_name", TEST_DATES)
    def test_sun_transit(
        self,
        compare: CompareHelper,
        name: str,
        lat: float,
        lon: float,
        alt: float,
        jd: float,
        date_name: str,
    ):
        """Sun transit timing matches within tolerance."""
        ref_result = compare.skyfield(ephem.rise_trans, jd, SUN, 4, (lon, lat, alt))
        leb_result = compare.leb(ephem.rise_trans, jd, SUN, 4, (lon, lat, alt))

        if ref_result[0] != 0 and leb_result[0] != 0:
            diff_sec = abs(ref_result[0] - leb_result[0]) * 86400.0
            assert diff_sec < TOLS.RISE_TRANSIT_SEC, (
                f"Sun transit at {name} ({date_name}): diff = {diff_sec:.1f}s"
            )


class TestMoonRiseSet:
    """Moon rise/set timing."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    @pytest.mark.parametrize("jd,date_name", TEST_DATES)
    def test_moonrise(
        self,
        compare: CompareHelper,
        name: str,
        lat: float,
        lon: float,
        alt: float,
        jd: float,
        date_name: str,
    ):
        """Moon rise timing matches within tolerance."""
        ref_result = compare.skyfield(ephem.rise_trans, jd, MOON, 1, (lon, lat, alt))
        leb_result = compare.leb(ephem.rise_trans, jd, MOON, 1, (lon, lat, alt))

        if ref_result[0] != 0 and leb_result[0] != 0:
            diff_sec = abs(ref_result[0] - leb_result[0]) * 86400.0
            # Moon can have slightly larger tolerance due to parallax
            assert diff_sec < TOLS.RISE_TRANSIT_SEC * 2, (
                f"Moonrise at {name} ({date_name}): diff = {diff_sec:.1f}s"
            )

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS)
    @pytest.mark.parametrize("jd,date_name", TEST_DATES)
    def test_moonset(
        self,
        compare: CompareHelper,
        name: str,
        lat: float,
        lon: float,
        alt: float,
        jd: float,
        date_name: str,
    ):
        """Moon set timing matches within tolerance."""
        ref_result = compare.skyfield(ephem.rise_trans, jd, MOON, 2, (lon, lat, alt))
        leb_result = compare.leb(ephem.rise_trans, jd, MOON, 2, (lon, lat, alt))

        if ref_result[0] != 0 and leb_result[0] != 0:
            diff_sec = abs(ref_result[0] - leb_result[0]) * 86400.0
            assert diff_sec < TOLS.RISE_TRANSIT_SEC * 2, (
                f"Moonset at {name} ({date_name}): diff = {diff_sec:.1f}s"
            )


class TestVenusRise:
    """Venus rise timing."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon,alt", TEST_LOCATIONS[:2])
    @pytest.mark.parametrize("jd,date_name", TEST_DATES)
    def test_venus_rise(
        self,
        compare: CompareHelper,
        name: str,
        lat: float,
        lon: float,
        alt: float,
        jd: float,
        date_name: str,
    ):
        """Venus rise timing matches within tolerance."""
        ref_result = compare.skyfield(ephem.rise_trans, jd, VENUS, 1, (lon, lat, alt))
        leb_result = compare.leb(ephem.rise_trans, jd, VENUS, 1, (lon, lat, alt))

        if ref_result[0] != 0 and leb_result[0] != 0:
            diff_sec = abs(ref_result[0] - leb_result[0]) * 86400.0
            assert diff_sec < TOLS.RISE_TRANSIT_SEC, (
                f"Venus rise at {name} ({date_name}): diff = {diff_sec:.1f}s"
            )
