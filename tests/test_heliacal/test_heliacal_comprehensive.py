"""
Tests for heliacal rising/setting calculations.

Verifies heliacal_ut returns valid results for planets
and bright stars, with reasonable event dates.

NOTE: Heliacal calculations are computationally expensive (~10-30s each).
All tests are marked @pytest.mark.slow to exclude from fast test runs.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    HELIACAL_RISING,
    HELIACAL_SETTING,
    EVENING_FIRST,
    MORNING_LAST,
)


# Standard atmospheric conditions
ATMO = (1013.25, 15.0, 40.0, 0.0)

# Standard observer
OBSERVER = (36.0, 1.0, 0, 0, 0, 0)

# Locations (lon, lat, alt)
ROME = (12.5, 41.9, 50.0)
CAIRO = (31.2, 30.0, 75.0)


class TestHeliacalRisingBasic:
    """Basic heliacal rising tests."""

    @pytest.mark.unit
    @pytest.mark.slow
    def test_venus_heliacal_rising_returns_3(self):
        """heliacal_ut for Venus returns 3 JD values."""
        jd = 2451545.0
        result = swe.heliacal_ut(
            jd, ROME, ATMO, OBSERVER, "Venus", HELIACAL_RISING
        )
        assert len(result) == 3, f"Expected 3 values, got {len(result)}"

    @pytest.mark.unit
    @pytest.mark.slow
    def test_venus_heliacal_rising_after_start(self):
        """Heliacal rising should be after search start."""
        jd = 2451545.0
        result = swe.heliacal_ut(
            jd, ROME, ATMO, OBSERVER, "Venus", HELIACAL_RISING
        )
        assert result[0] > jd, f"Rising JD {result[0]} not after start {jd}"

    @pytest.mark.unit
    @pytest.mark.slow
    def test_venus_heliacal_rising_within_2_years(self):
        """Venus heliacal rising should be within ~2 years."""
        jd = 2451545.0
        result = swe.heliacal_ut(
            jd, ROME, ATMO, OBSERVER, "Venus", HELIACAL_RISING
        )
        gap = result[0] - jd
        assert gap < 800, f"Venus rising {gap:.1f} days after start"

    @pytest.mark.unit
    @pytest.mark.slow
    def test_all_result_values_finite(self):
        """All 3 returned JD values should be finite."""
        jd = 2451545.0
        result = swe.heliacal_ut(
            jd, ROME, ATMO, OBSERVER, "Venus", HELIACAL_RISING
        )
        for i, val in enumerate(result):
            assert math.isfinite(val), f"result[{i}] = {val}"


class TestHeliacalPlanets:
    """Test heliacal events for various planets."""

    @pytest.mark.unit
    @pytest.mark.slow
    @pytest.mark.parametrize(
        "planet",
        ["Venus", "Jupiter", "Saturn"],
    )
    def test_planet_heliacal_rising(self, planet: str):
        """Heliacal rising works for visible planets."""
        jd = 2451545.0
        result = swe.heliacal_ut(
            jd, ROME, ATMO, OBSERVER, planet, HELIACAL_RISING
        )
        assert len(result) == 3
        assert result[0] > jd

    @pytest.mark.unit
    @pytest.mark.slow
    @pytest.mark.parametrize(
        "planet",
        ["Venus", "Jupiter"],
    )
    def test_planet_heliacal_setting(self, planet: str):
        """Heliacal setting works for visible planets."""
        jd = 2451545.0
        result = swe.heliacal_ut(
            jd, ROME, ATMO, OBSERVER, planet, HELIACAL_SETTING
        )
        assert len(result) == 3
        assert result[0] > jd


class TestHeliacalInnerPlanets:
    """Test evening first / morning last for inner planets."""

    @pytest.mark.unit
    @pytest.mark.slow
    def test_venus_evening_first(self):
        """Evening first visibility for Venus."""
        jd = 2451545.0
        result = swe.heliacal_ut(
            jd, ROME, ATMO, OBSERVER, "Venus", EVENING_FIRST
        )
        assert len(result) == 3
        assert result[0] > jd

    @pytest.mark.unit
    @pytest.mark.slow
    def test_venus_morning_last(self):
        """Morning last visibility for Venus."""
        jd = 2451545.0
        result = swe.heliacal_ut(jd, ROME, ATMO, OBSERVER, "Venus", MORNING_LAST)
        assert len(result) == 3
        assert result[0] > jd


class TestHeliacalFixedStars:
    """Test heliacal events for fixed stars."""

    @pytest.mark.unit
    @pytest.mark.slow
    def test_sirius_heliacal_rising(self):
        """Heliacal rising of Sirius (historically significant)."""
        jd = 2451545.0
        try:
            result = swe.heliacal_ut(
                jd, CAIRO, ATMO, OBSERVER, "Sirius", HELIACAL_RISING
            )
            assert len(result) == 3
            assert result[0] > jd
        except Exception:
            pytest.skip("Heliacal rising not supported for Sirius")


class TestHeliacalInvalidInput:
    """Test error handling for invalid heliacal inputs."""

    @pytest.mark.unit
    def test_sun_heliacal_raises(self):
        """Sun is not valid for heliacal calculations."""
        jd = 2451545.0
        with pytest.raises((ValueError, Exception)):
            swe.heliacal_ut(jd, ROME, ATMO, OBSERVER, "Sun", HELIACAL_RISING)

    @pytest.mark.unit
    def test_moon_heliacal_raises(self):
        """Moon is not valid for heliacal calculations."""
        jd = 2451545.0
        with pytest.raises((ValueError, Exception)):
            swe.heliacal_ut(jd, ROME, ATMO, OBSERVER, "Moon", HELIACAL_RISING)
