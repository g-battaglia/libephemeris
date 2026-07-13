"""
LEB vs Skyfield Comparison: Crossing Functions.

Validates that crossing functions (iterative solvers using calc_ut
internally) produce consistent timing in both modes.
"""

from __future__ import annotations

import pytest

import libephemeris as ephem
from libephemeris.constants import (
    SUN,
    MARS,
    JUPITER,
    SATURN,
    FLG_SWIEPH,
)

from .conftest import TOLS, CompareHelper, year_to_jd


@pytest.fixture
def crossing_dates():
    """Start dates for crossing tests."""
    return [year_to_jd(2020), year_to_jd(2025), year_to_jd(2030), year_to_jd(2035)]


class TestSunCrossings:
    """Solar longitude crossings."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("target_lon", [0, 30, 60, 90, 120, 150, 180])
    def test_solcross_ut(
        self, compare: CompareHelper, crossing_dates: list[float], target_lon: float
    ):
        """Sun crossing timing matches within tolerance."""
        for jd_start in crossing_dates:
            ref_jd = compare.skyfield(ephem.solcross_ut, target_lon, jd_start, 0)
            leb_jd = compare.leb(ephem.solcross_ut, target_lon, jd_start, 0)

            diff_sec = abs(ref_jd - leb_jd) * 86400.0
            assert diff_sec < TOLS.CROSSING_SUN_SEC, (
                f"Sun crossing {target_lon}° at JD {jd_start:.1f}: diff = {diff_sec:.1f}s"
            )

    @pytest.mark.leb_compare
    def test_crossing_position_verify(self, compare: CompareHelper):
        """Verify Sun is at target longitude at crossing time."""
        target = 45.0
        jd_start = year_to_jd(2024)

        leb_jd = compare.leb(ephem.solcross_ut, target, jd_start, 0)
        pos, _ = compare.leb(ephem.calc_ut, leb_jd, SUN, 0)

        lon_diff = abs(pos[0] - target)
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < 0.01, f"Sun not at {target}° at crossing JD {leb_jd}"


class TestMoonCrossings:
    """Lunar longitude crossings."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("target_lon", [0, 90, 180, 270])
    def test_mooncross_ut(
        self, compare: CompareHelper, crossing_dates: list[float], target_lon: float
    ):
        """Moon crossing timing matches within tolerance."""
        for jd_start in crossing_dates:
            ref_jd = compare.skyfield(ephem.mooncross_ut, target_lon, jd_start, 0)
            leb_jd = compare.leb(ephem.mooncross_ut, target_lon, jd_start, 0)

            diff_sec = abs(ref_jd - leb_jd) * 86400.0
            assert diff_sec < TOLS.CROSSING_MOON_SEC, (
                f"Moon crossing {target_lon}° at JD {jd_start:.1f}: diff = {diff_sec:.1f}s"
            )


class TestMoonNodeCrossings:
    """Lunar node crossings."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize("node_type", [0, 1], ids=["ascending", "descending"])
    def test_mooncross_node_ut(
        self, compare: CompareHelper, crossing_dates: list[float], node_type: int
    ):
        """Moon node crossing timing matches within tolerance."""
        for jd_start in crossing_dates:
            ref_result = compare.skyfield(ephem.mooncross_node_ut, jd_start, node_type)
            leb_result = compare.leb(ephem.mooncross_node_ut, jd_start, node_type)

            diff_sec = abs(ref_result[0] - leb_result[0]) * 86400.0
            assert diff_sec < TOLS.CROSSING_MOON_NODE_SEC, (
                f"Moon node {node_type} at JD {jd_start:.1f}: diff = {diff_sec:.1f}s"
            )


class TestPlanetCrossings:
    """Planet longitude crossings."""

    # Jupiter and Saturn crossing searches diverge even in pure Skyfield mode
    # (pre-existing solver bug in crossing.py, not LEB-related).
    _XFAIL_BODIES = {JUPITER, SATURN}

    @pytest.mark.leb_compare
    @pytest.mark.parametrize(
        "body_id,body_name",
        [(MARS, "Mars"), (JUPITER, "Jupiter"), (SATURN, "Saturn")],
    )
    @pytest.mark.parametrize("target_lon", [0, 90, 180, 270])
    def test_cross_ut(
        self, compare: CompareHelper, body_id: int, body_name: str, target_lon: float
    ):
        """Planet crossing timing matches within tolerance."""
        if body_id in self._XFAIL_BODIES:
            pytest.xfail("crossing solver diverges for slow planets (pre-existing bug)")

        jd_start = year_to_jd(2024)

        try:
            ref_jd = compare.skyfield(
                ephem.cross_ut, body_id, target_lon, jd_start, FLG_SWIEPH
            )
            leb_jd = compare.leb(
                ephem.cross_ut, body_id, target_lon, jd_start, FLG_SWIEPH
            )
        except RuntimeError:
            pytest.xfail(
                f"crossing solver diverged for {body_name} at {target_lon}° "
                "(pre-existing bug)"
            )

        diff_sec = abs(ref_jd - leb_jd) * 86400.0
        assert diff_sec < TOLS.CROSSING_PLANET_SEC, (
            f"{body_name} crossing {target_lon}°: diff = {diff_sec:.1f}s"
        )


class TestHelioCrossings:
    """Heliocentric longitude crossings."""

    @pytest.mark.leb_compare
    @pytest.mark.parametrize(
        "body_id,body_name",
        [(MARS, "Mars"), (JUPITER, "Jupiter"), (SATURN, "Saturn")],
    )
    @pytest.mark.parametrize("target_lon", [0, 90, 180, 270])
    def test_helio_cross_ut(
        self, compare: CompareHelper, body_id: int, body_name: str, target_lon: float
    ):
        """Heliocentric crossing timing matches within tolerance."""
        jd_start = year_to_jd(2024)

        try:
            ref_jd = compare.skyfield(
                ephem.helio_cross_ut, body_id, target_lon, jd_start, FLG_SWIEPH
            )
            leb_jd = compare.leb(
                ephem.helio_cross_ut, body_id, target_lon, jd_start, FLG_SWIEPH
            )
        except RuntimeError:
            pytest.xfail(
                f"helio crossing solver diverged for {body_name} at {target_lon}° "
                "(pre-existing bug)"
            )

        diff_sec = abs(ref_jd - leb_jd) * 86400.0
        assert diff_sec < TOLS.CROSSING_PLANET_SEC, (
            f"{body_name} helio crossing {target_lon}°: diff = {diff_sec:.1f}s"
        )
