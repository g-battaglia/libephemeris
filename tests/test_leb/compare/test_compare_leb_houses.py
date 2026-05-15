"""
LEB vs Skyfield Comparison: House Calculations.

Validates that house calculations are unaffected by LEB mode for
non-Sunshine systems, and correct for Sunshine ('I') system.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import FLG_SIDEREAL

from .conftest import TOLS, CompareHelper, year_to_jd, angular_diff


HOUSE_SYSTEMS_NON_SUNSHINE = [
    "P",
    "K",
    "R",
    "C",
    "E",
    "A",
    "W",
    "O",
    "B",
    "T",
    "M",
    "X",
    "V",
    "H",
    "F",
    "U",
    "N",
    "Y",
    "D",
    "L",
    "S",
    "G",
    "Q",
]

TEST_LOCATIONS = [
    ("Rome", 41.9, 12.5),
    ("NewYork", 40.7, -74.0),
    ("Sydney", -33.9, 151.2),
    ("Equator", 0.0, 0.0),
]

TEST_DATES = [year_to_jd(2000), year_to_jd(2024), year_to_jd(2050)]


class TestNonSunshineIdentical:
    """Non-Sunshine systems must be bit-for-bit identical."""

    @pytest.mark.leb_compare
    @pytest.mark.slow
    @pytest.mark.parametrize("hsys", HOUSE_SYSTEMS_NON_SUNSHINE)
    @pytest.mark.parametrize("name,lat,lon", TEST_LOCATIONS)
    @pytest.mark.parametrize("jd", TEST_DATES)
    def test_cusps_identical(
        self,
        compare: CompareHelper,
        hsys: str,
        name: str,
        lat: float,
        lon: float,
        jd: float,
    ):
        """Non-Sunshine house cusps must be identical in both modes."""
        ref_cusps, ref_ascmc = compare.skyfield(ephem.houses, jd, lat, lon, hsys)
        leb_cusps, leb_ascmc = compare.leb(ephem.houses, jd, lat, lon, hsys)

        # All cusps must be exactly identical
        for i in range(min(len(ref_cusps), len(leb_cusps))):
            assert ref_cusps[i] == pytest.approx(leb_cusps[i], rel=1e-12), (
                f"{hsys} at {name}: cusp {i + 1} differs"
            )

        # ASC/MC must be exactly identical
        for i in range(min(len(ref_ascmc), len(leb_ascmc))):
            assert ref_ascmc[i] == pytest.approx(leb_ascmc[i], rel=1e-12), (
                f"{hsys} at {name}: ASCMC[{i}] differs"
            )


class TestSunshinePrecision:
    """Sunshine system precision (uses calc_ut for Sun declination)."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("name,lat,lon", TEST_LOCATIONS)
    @pytest.mark.parametrize("jd", TEST_DATES)
    def test_sunshine_cusps(
        self, compare: CompareHelper, name: str, lat: float, lon: float, jd: float
    ):
        """Sunshine cusps match within tolerance."""
        ref_cusps, _ = compare.skyfield(ephem.houses, jd, lat, lon, "I")
        leb_cusps, _ = compare.leb(ephem.houses, jd, lat, lon, "I")

        for i in range(min(len(ref_cusps), len(leb_cusps))):
            diff = angular_diff(ref_cusps[i], leb_cusps[i]) * 3600.0
            assert diff < TOLS.HOUSE_SUNSHINE_ARCSEC, (
                f'Sunshine at {name}: cusp {i + 1} diff = {diff:.4f}"'
            )


class TestHousesExSidereal:
    """houses_ex with sidereal flag."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("sid_mode", [0, 1, 2, 3])
    @pytest.mark.parametrize("name,lat,lon", TEST_LOCATIONS[:2])
    def test_sidereal_houses(
        self, compare: CompareHelper, sid_mode: int, name: str, lat: float, lon: float
    ):
        """Sidereal houses must be identical (both use same calc_ut path)."""
        jd = year_to_jd(2024)
        flags = FLG_SIDEREAL
        ephem.set_sid_mode(sid_mode, 2451545.0, 0.0)

        hsys = ord("P")  # Placidus
        ref_cusps, _ = compare.skyfield(ephem.houses_ex, jd, lat, lon, hsys, flags)
        leb_cusps, _ = compare.leb(ephem.houses_ex, jd, lat, lon, hsys, flags)

        for i in range(12):
            assert ref_cusps[i] == pytest.approx(leb_cusps[i], rel=1e-10), (
                f"Sidereal mode {sid_mode} at {name}: cusp {i + 1} differs"
            )
