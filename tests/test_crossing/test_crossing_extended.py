"""
Tests for crossing and station functions.

Verifies solcross_ut, mooncross_ut, mooncross_node_ut,
and find_station_ut produce valid results.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    MERCURY,
    MARS,
    FLG_SWIEPH,
)


class TestSolcross:
    """Test Sun crossing calculations."""

    @pytest.mark.unit
    def test_solcross_vernal_equinox(self):
        """Sun crosses 0 Aries near March 20."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        jd_cross = swe.solcross_ut(0.0, jd_start, FLG_SWIEPH)
        assert isinstance(jd_cross, float)
        assert math.isfinite(jd_cross)
        # Should be around March 20, 2024
        y, m, d, h = swe.revjul(jd_cross)
        assert y == 2024
        assert m == 3
        assert 19 <= d <= 21, f"Vernal equinox: {y}-{m}-{d}"

    @pytest.mark.unit
    def test_solcross_summer_solstice(self):
        """Sun crosses 90 deg near June 21."""
        jd_start = swe.julday(2024, 4, 1, 0.0)
        jd_cross = swe.solcross_ut(90.0, jd_start, FLG_SWIEPH)
        y, m, d, h = swe.revjul(jd_cross)
        assert y == 2024
        assert m == 6
        assert 20 <= d <= 22

    @pytest.mark.unit
    def test_solcross_autumnal_equinox(self):
        """Sun crosses 180 deg near September 22."""
        jd_start = swe.julday(2024, 7, 1, 0.0)
        jd_cross = swe.solcross_ut(180.0, jd_start, FLG_SWIEPH)
        y, m, d, h = swe.revjul(jd_cross)
        assert y == 2024
        assert m == 9
        assert 21 <= d <= 23

    @pytest.mark.unit
    def test_solcross_winter_solstice(self):
        """Sun crosses 270 deg near December 21."""
        jd_start = swe.julday(2024, 10, 1, 0.0)
        jd_cross = swe.solcross_ut(270.0, jd_start, FLG_SWIEPH)
        y, m, d, h = swe.revjul(jd_cross)
        assert y == 2024
        assert m == 12
        assert 20 <= d <= 22

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "degree",
        [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0],
    )
    def test_solcross_all_signs(self, degree: float):
        """Sun should cross every 30 degree mark within a year."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        jd_cross = swe.solcross_ut(degree, jd_start, FLG_SWIEPH)
        assert math.isfinite(jd_cross)
        assert jd_cross > jd_start
        # Should be within ~366 days
        assert jd_cross - jd_start < 367


class TestMooncross:
    """Test Moon crossing calculations."""

    @pytest.mark.unit
    def test_mooncross_0_aries(self):
        """Moon crosses 0 Aries within ~28 days."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        jd_cross = swe.mooncross_ut(0.0, jd_start, FLG_SWIEPH)
        assert isinstance(jd_cross, float)
        assert math.isfinite(jd_cross)
        assert jd_cross > jd_start
        # Moon crosses each point roughly every 27.3 days
        assert jd_cross - jd_start < 29

    @pytest.mark.unit
    @pytest.mark.parametrize("degree", [0.0, 90.0, 180.0, 270.0])
    def test_mooncross_quadrants(self, degree: float):
        """Moon crosses quadrant points within ~28 days."""
        jd_start = swe.julday(2024, 6, 1, 0.0)
        jd_cross = swe.mooncross_ut(degree, jd_start, FLG_SWIEPH)
        assert math.isfinite(jd_cross)
        assert jd_cross > jd_start
        assert jd_cross - jd_start < 29

    @pytest.mark.unit
    def test_mooncross_consecutive(self):
        """Two consecutive Moon crossings of same point ~27.3 days apart."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        jd1 = swe.mooncross_ut(90.0, jd_start, FLG_SWIEPH)
        jd2 = swe.mooncross_ut(90.0, jd1 + 1.0, FLG_SWIEPH)
        gap = jd2 - jd1
        assert 26 < gap < 29, f"Gap between Moon crossings: {gap:.2f} days"


class TestMooncrossNode:
    """Test Moon node crossing calculations."""

    @pytest.mark.unit
    def test_mooncross_node_returns_tuple(self):
        """mooncross_node_ut returns (jd, lon, lat)."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        result = swe.mooncross_node_ut(jd_start, FLG_SWIEPH)
        assert len(result) == 3
        jd_cross, lon, lat = result
        assert math.isfinite(jd_cross)
        assert math.isfinite(lon)
        assert math.isfinite(lat)

    @pytest.mark.unit
    def test_mooncross_node_lat_near_zero(self):
        """At node crossing, latitude should be near 0."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        jd_cross, lon, lat = swe.mooncross_node_ut(jd_start, FLG_SWIEPH)
        # At the node, latitude should be very close to 0
        assert abs(lat) < 0.1, f"Moon lat at node: {lat} deg"

    @pytest.mark.unit
    def test_mooncross_node_within_month(self):
        """Moon should cross a node within ~14 days."""
        jd_start = swe.julday(2024, 6, 1, 0.0)
        jd_cross, _, _ = swe.mooncross_node_ut(jd_start, FLG_SWIEPH)
        assert jd_cross > jd_start
        # Moon crosses nodes roughly every 13.6 days
        assert jd_cross - jd_start < 15


class TestFindStation:
    """Test planetary station finding."""

    @pytest.mark.unit
    def test_mercury_station_found(self):
        """Mercury stations should be found within ~4 months."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        result = swe.find_station_ut(MERCURY, jd_start, "any", FLG_SWIEPH)
        # Returns (jd, station_type_str)
        jd_station, stype = result
        assert isinstance(jd_station, float)
        assert math.isfinite(jd_station)
        assert jd_station > jd_start
        assert stype in ("SR", "SD"), f"station type: {stype}"
        # Mercury retrogrades ~3 times/year, so station within ~122 days
        assert jd_station - jd_start < 130

    @pytest.mark.unit
    def test_mars_station_found(self):
        """Mars station should be found."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        result = swe.find_station_ut(MARS, jd_start, "any", FLG_SWIEPH)
        jd_station, stype = result
        assert isinstance(jd_station, float)
        assert math.isfinite(jd_station)
        assert jd_station > jd_start
        assert stype in ("SR", "SD")

    @pytest.mark.unit
    @pytest.mark.parametrize("station_type", ["SR", "SD"])
    def test_mercury_station_types(self, station_type: str):
        """Mercury SR and SD stations should be findable."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        result = swe.find_station_ut(
            MERCURY, jd_start, station_type, FLG_SWIEPH
        )
        jd_station, stype = result
        assert isinstance(jd_station, float)
        assert math.isfinite(jd_station)
        assert stype == station_type, f"Expected {station_type}, got {stype}"
