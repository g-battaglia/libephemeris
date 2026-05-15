"""
Tests for planetary moon calculations.

Planetary moons (9000-series body IDs) currently return placeholder
values (all zeros). These tests verify the API contract: that the
functions accept these body IDs without crashing, and return the
expected 6-element tuple format.

When full planetary moon support is implemented, these tests should
be updated to verify actual positions.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.constants import (
    FLG_SPEED,
    FLG_EQUATORIAL,
)


# Planetary moon body IDs (9000 series)
MOON_IO = 9001
MOON_EUROPA = 9002
MOON_GANYMEDE = 9003
MOON_CALLISTO = 9004

GALILEAN_MOONS = [
    (MOON_IO, "Io"),
    (MOON_EUROPA, "Europa"),
    (MOON_GANYMEDE, "Ganymede"),
    (MOON_CALLISTO, "Callisto"),
]


class TestGalileanMoonsAPI:
    """Test API contract for Galilean moon body IDs."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", GALILEAN_MOONS)
    def test_galilean_moon_returns_6_tuple(self, body_id: int, name: str):
        """Each Galilean moon returns a 6-element tuple without crashing."""
        jd = 2451545.0
        result, retflag = swe.calc_ut(jd, body_id, FLG_SPEED)
        assert len(result) == 6, f"{name}: expected 6 elements, got {len(result)}"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", GALILEAN_MOONS)
    def test_galilean_moon_all_finite(self, body_id: int, name: str):
        """All output values should be finite (even if zero)."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
        for i, val in enumerate(result):
            assert math.isfinite(val), f"{name}: result[{i}] = {val} not finite"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", GALILEAN_MOONS)
    def test_galilean_moon_returns_native_float(self, body_id: int, name: str):
        """All return values should be native Python float."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, 0)
        for i, val in enumerate(result):
            assert type(val) is float, (
                f"{name}: result[{i}] is {type(val).__name__}, expected float"
            )


class TestGalileanMoonsFlagCombos:
    """Test Galilean moons accept various flag combinations without crashing."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "flags,desc",
        [
            (0, "default"),
            (FLG_SPEED, "speed"),
            (FLG_EQUATORIAL, "equatorial"),
            (FLG_SPEED | FLG_EQUATORIAL, "speed+equatorial"),
        ],
    )
    @pytest.mark.parametrize(
        "body_id,name",
        [
            (MOON_IO, "Io"),
            (MOON_CALLISTO, "Callisto"),
        ],
    )
    def test_moon_flag_combo_no_crash(
        self, body_id: int, name: str, flags: int, desc: str
    ):
        """Galilean moons don't crash with various flags."""
        jd = 2451545.0
        result, _ = swe.calc_ut(jd, body_id, flags)
        assert len(result) == 6


class TestGalileanMoonsDateRange:
    """Test Galilean moons across dates — verify no crashes."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", GALILEAN_MOONS)
    def test_galilean_across_years(self, body_id: int, name: str):
        """Galilean moons don't crash at multiple years."""
        for year in [1900, 1950, 2000, 2024, 2050, 2100]:
            jd = swe.julday(year, 1, 1, 12.0)
            result, _ = swe.calc_ut(jd, body_id, FLG_SPEED)
            assert len(result) == 6, f"{name} @ {year}"
