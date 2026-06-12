"""
Tests for crossing functions and planetary stations.

Verifies solcross_ut, mooncross_ut, mooncross_node_ut,
and find_station_ut return valid results.
"""

from __future__ import annotations


import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    FLG_SPEED,
)


class TestSolcross:
    """Tests for Sun crossing a given longitude."""

    @pytest.mark.unit
    def test_solcross_returns_jd(self):
        """solcross_ut returns a Julian Day number."""
        jd_start = 2451545.0
        result = swe.solcross_ut(0.0, jd_start, 0)
        assert isinstance(result, float)
        assert result > jd_start, "Crossing should be after start"

    @pytest.mark.unit
    @pytest.mark.parametrize("degree", [0, 30, 60, 90, 120, 180, 270, 330])
    def test_solcross_at_various_degrees(self, degree: int):
        """Sun crossing found at various longitudes."""
        jd_start = 2451545.0
        jd_cross = swe.solcross_ut(float(degree), jd_start, 0)
        assert jd_cross > jd_start
        # Sun should be within 1 year
        assert jd_cross - jd_start < 366

        # Verify Sun is actually near the target degree
        result, _ = swe.calc_ut(jd_cross, SUN, 0)
        diff = abs(result[0] - degree)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.01, (
            f"Sun at {result[0]:.4f}° but expected {degree}° (diff {diff:.4f}°)"
        )

    @pytest.mark.unit
    def test_solcross_vernal_equinox(self):
        """Find vernal equinox (Sun at 0° Aries)."""
        # Search from 2024-01-01
        jd_start = swe.julday(2024, 1, 1, 0.0)
        jd_cross = swe.solcross_ut(0.0, jd_start, 0)
        # Should be around March 20, 2024
        y, m, d, h = swe.revjul(jd_cross)
        assert y == 2024
        assert m == 3
        assert 19 <= d <= 21, f"Vernal equinox on March {d}"

    @pytest.mark.unit
    def test_solcross_summer_solstice(self):
        """Find summer solstice (Sun at 90°)."""
        jd_start = swe.julday(2024, 1, 1, 0.0)
        jd_cross = swe.solcross_ut(90.0, jd_start, 0)
        y, m, d, h = swe.revjul(jd_cross)
        assert y == 2024
        assert m == 6
        assert 20 <= d <= 22, f"Summer solstice on June {d}"


class TestMooncross:
    """Tests for Moon crossing a given longitude."""

    @pytest.mark.unit
    def test_mooncross_returns_jd(self):
        """mooncross_ut returns a Julian Day number."""
        jd_start = 2451545.0
        result = swe.mooncross_ut(0.0, jd_start, 0)
        assert isinstance(result, float)
        assert result > jd_start

    @pytest.mark.unit
    @pytest.mark.parametrize("degree", [0, 90, 180, 270])
    def test_mooncross_at_cardinal_degrees(self, degree: int):
        """Moon crossing found at cardinal longitudes."""
        jd_start = 2451545.0
        jd_cross = swe.mooncross_ut(float(degree), jd_start, 0)
        assert jd_cross > jd_start
        # Moon should cross within ~28 days
        assert jd_cross - jd_start < 30

        # Verify Moon is near target
        result, _ = swe.calc_ut(jd_cross, MOON, 0)
        diff = abs(result[0] - degree)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Moon at {result[0]:.4f}° but expected {degree}°"

    @pytest.mark.unit
    def test_mooncross_sequential(self):
        """Find 3 consecutive Moon crossings of 0°."""
        jd = 2451545.0
        crossings = []
        for _ in range(3):
            jd_cross = swe.mooncross_ut(0.0, jd, 0)
            crossings.append(jd_cross)
            jd = jd_cross + 1.0

        # Consecutive crossings should be ~27.3 days apart (sidereal month)
        for i in range(2):
            gap = crossings[i + 1] - crossings[i]
            assert 25 < gap < 30, f"Gap between crossings: {gap:.2f} days"


class TestMooncrossNode:
    """Tests for Moon crossing its own node."""

    @pytest.mark.unit
    def test_mooncross_node_returns_tuple(self):
        """mooncross_node_ut returns (jd, longitude, latitude)."""
        jd_start = 2451545.0
        result = swe.mooncross_node_ut(jd_start, 0)
        assert len(result) >= 2, f"Expected >= 2 elements, got {len(result)}"

    @pytest.mark.unit
    def test_mooncross_node_after_start(self):
        """Node crossing should be after search start."""
        jd_start = 2451545.0
        result = swe.mooncross_node_ut(jd_start, 0)
        jd_cross = result[0]
        assert jd_cross > jd_start

    @pytest.mark.unit
    def test_mooncross_node_within_month(self):
        """Node crossing should be within ~15 days (half nodal period)."""
        jd_start = 2451545.0
        result = swe.mooncross_node_ut(jd_start, 0)
        jd_cross = result[0]
        assert jd_cross - jd_start < 16, (
            f"Node crossing {jd_cross - jd_start:.1f} days after start"
        )

    @pytest.mark.unit
    def test_mooncross_node_sequential(self):
        """Find 4 consecutive node crossings."""
        jd = 2451545.0
        crossings = []
        for _ in range(4):
            result = swe.mooncross_node_ut(jd, 0)
            crossings.append(result[0])
            jd = result[0] + 1.0

        # Node crossings occur ~13.6 days apart
        for i in range(3):
            gap = crossings[i + 1] - crossings[i]
            assert 12 < gap < 16, f"Gap between node crossings: {gap:.2f} days"


class TestFindStation:
    """Tests for finding planetary stations (retrograde/direct)."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_find_station_returns_valid(self, body_id: int, name: str):
        """find_station_ut returns (jd, type_str) for each planet."""
        jd_start = 2451545.0
        result = swe.find_station_ut(body_id, jd_start, "any", 0)
        assert len(result) >= 2, f"{name}: expected >= 2 elements"
        jd_station = result[0]
        assert jd_station > jd_start, (
            f"{name}: station {jd_station} not after start {jd_start}"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,name,max_days",
        [
            (MERCURY, "Mercury", 120),  # Mercury stations ~3-4 times/year
            (JUPITER, "Jupiter", 400),  # Jupiter stations ~once/year
            (SATURN, "Saturn", 400),  # Saturn stations ~once/year
        ],
    )
    def test_station_within_reasonable_time(
        self, body_id: int, name: str, max_days: int
    ):
        """Station should be found within reasonable time."""
        jd_start = 2451545.0
        result = swe.find_station_ut(body_id, jd_start, "any", 0)
        jd_station = result[0]
        gap = jd_station - jd_start
        assert gap < max_days, f"{name}: station {gap:.1f} days after start"

    @pytest.mark.unit
    def test_mercury_station_speed_near_zero(self):
        """At a station, Mercury's speed should be near zero."""
        jd_start = 2451545.0
        result = swe.find_station_ut(MERCURY, jd_start, "any", 0)
        jd_station = result[0]

        # Check speed at station
        pos, _ = swe.calc_ut(jd_station, MERCURY, FLG_SPEED)
        speed = pos[3]
        assert abs(speed) < 0.05, (
            f"Mercury speed at station: {speed}°/day (expected ~0)"
        )
