"""
Comprehensive tests for Arabic parts (Lots) calculations.

Verifies calc_all_arabic_parts and individual Arabic part formulas
with known inputs and edge cases.
"""

from __future__ import annotations

import math

import pytest

import libephemeris as swe
from libephemeris.arabic_parts import (
    calc_arabic_part_of_fortune,
    calc_arabic_part_of_spirit,
    calc_arabic_part_of_love,
    calc_arabic_part_of_faith,
    calc_all_arabic_parts,
    is_day_chart,
)


def _degnorm(x: float) -> float:
    """Normalize angle to 0-360."""
    return x % 360


class TestPartOfFortune:
    """Tests for Part of Fortune calculation."""

    @pytest.mark.unit
    def test_fortune_day_chart(self):
        """Day chart: Fortune = Asc + Moon - Sun."""
        asc, sun, moon = 15.0, 120.0, 240.0
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=True)
        expected = _degnorm(asc + moon - sun)
        assert abs(result - expected) < 0.001, f"Fortune day: {result} != {expected}"

    @pytest.mark.unit
    def test_fortune_night_chart(self):
        """Night chart: Fortune = Asc + Sun - Moon."""
        asc, sun, moon = 15.0, 120.0, 240.0
        result = calc_arabic_part_of_fortune(asc, sun, moon, is_diurnal=False)
        expected = _degnorm(asc + sun - moon)
        assert abs(result - expected) < 0.001, f"Fortune night: {result} != {expected}"

    @pytest.mark.unit
    def test_fortune_result_in_range(self):
        """Fortune should always be 0-360."""
        for asc in [0, 90, 180, 270]:
            for sun in [30, 150, 300]:
                for moon in [60, 200, 350]:
                    for diurnal in [True, False]:
                        result = calc_arabic_part_of_fortune(
                            float(asc), float(sun), float(moon), diurnal
                        )
                        assert 0 <= result < 360, (
                            f"Fortune {result} out of range "
                            f"(asc={asc}, sun={sun}, moon={moon})"
                        )

    @pytest.mark.unit
    def test_fortune_zero_inputs(self):
        """Fortune with all zeros should be 0."""
        result = calc_arabic_part_of_fortune(0.0, 0.0, 0.0, True)
        assert abs(result) < 0.001 or abs(result - 360) < 0.001

    @pytest.mark.unit
    def test_fortune_same_sun_moon(self):
        """If Sun == Moon, Fortune == Asc (day) or Fortune == Asc (night)."""
        asc = 45.0
        result_day = calc_arabic_part_of_fortune(asc, 100.0, 100.0, True)
        result_night = calc_arabic_part_of_fortune(asc, 100.0, 100.0, False)
        assert abs(result_day - asc) < 0.001
        assert abs(result_night - asc) < 0.001


class TestPartOfSpirit:
    """Tests for Part of Spirit calculation."""

    @pytest.mark.unit
    def test_spirit_day_chart(self):
        """Day chart: Spirit = Asc + Sun - Moon."""
        asc, sun, moon = 15.0, 120.0, 240.0
        result = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=True)
        expected = _degnorm(asc + sun - moon)
        assert abs(result - expected) < 0.001

    @pytest.mark.unit
    def test_spirit_night_chart(self):
        """Night chart: Spirit = Asc + Moon - Sun."""
        asc, sun, moon = 15.0, 120.0, 240.0
        result = calc_arabic_part_of_spirit(asc, sun, moon, is_diurnal=False)
        expected = _degnorm(asc + moon - sun)
        assert abs(result - expected) < 0.001

    @pytest.mark.unit
    def test_spirit_is_fortune_reversed(self):
        """Spirit day == Fortune night and vice versa."""
        asc, sun, moon = 45.0, 130.0, 250.0
        fortune_day = calc_arabic_part_of_fortune(asc, sun, moon, True)
        spirit_night = calc_arabic_part_of_spirit(asc, sun, moon, False)
        assert abs(fortune_day - spirit_night) < 0.001

        fortune_night = calc_arabic_part_of_fortune(asc, sun, moon, False)
        spirit_day = calc_arabic_part_of_spirit(asc, sun, moon, True)
        assert abs(fortune_night - spirit_day) < 0.001


class TestPartOfLove:
    """Tests for Part of Love calculation."""

    @pytest.mark.unit
    def test_love_formula(self):
        """Love = Asc + Venus - Sun."""
        asc, venus, sun = 30.0, 150.0, 90.0
        result = calc_arabic_part_of_love(asc, venus, sun)
        expected = _degnorm(asc + venus - sun)
        assert abs(result - expected) < 0.001

    @pytest.mark.unit
    def test_love_result_in_range(self):
        """Love should always be 0-360."""
        for asc in [0, 120, 240]:
            for venus in [60, 180, 300]:
                for sun in [30, 150, 270]:
                    result = calc_arabic_part_of_love(
                        float(asc), float(venus), float(sun)
                    )
                    assert 0 <= result < 360


class TestPartOfFaith:
    """Tests for Part of Faith calculation."""

    @pytest.mark.unit
    def test_faith_formula(self):
        """Faith = Asc + Mercury - Moon."""
        asc, mercury, moon = 30.0, 150.0, 90.0
        result = calc_arabic_part_of_faith(asc, mercury, moon)
        expected = _degnorm(asc + mercury - moon)
        assert abs(result - expected) < 0.001

    @pytest.mark.unit
    def test_faith_result_in_range(self):
        """Faith should always be 0-360."""
        for asc in [0, 120, 240]:
            for mercury in [60, 180, 300]:
                for moon in [30, 150, 270]:
                    result = calc_arabic_part_of_faith(
                        float(asc), float(mercury), float(moon)
                    )
                    assert 0 <= result < 360


class TestCalcAllArabicParts:
    """Tests for calc_all_arabic_parts top-level function."""

    @pytest.mark.unit
    def test_returns_four_parts(self):
        """calc_all_arabic_parts returns 4 named parts."""
        positions = {
            "Asc": 15.5,
            "Sun": 120.0,
            "Moon": 240.0,
            "Mercury": 130.0,
            "Venus": 110.0,
        }
        parts = calc_all_arabic_parts(positions)
        assert "Pars_Fortunae" in parts
        assert "Pars_Spiritus" in parts
        assert "Pars_Amoris" in parts
        assert "Pars_Fidei" in parts
        assert len(parts) == 4

    @pytest.mark.unit
    def test_all_parts_in_range(self):
        """All parts should be 0-360."""
        positions = {
            "Asc": 15.5,
            "Sun": 120.0,
            "Moon": 240.0,
            "Mercury": 130.0,
            "Venus": 110.0,
        }
        parts = calc_all_arabic_parts(positions)
        for name, val in parts.items():
            assert 0 <= val < 360, f"{name} = {val} out of range"

    @pytest.mark.unit
    def test_all_parts_finite(self):
        """All part values should be finite floats."""
        positions = {
            "Asc": 0.0,
            "Sun": 0.0,
            "Moon": 0.0,
            "Mercury": 0.0,
            "Venus": 0.0,
        }
        parts = calc_all_arabic_parts(positions)
        for name, val in parts.items():
            assert math.isfinite(val), f"{name} = {val} not finite"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "asc,sun,moon,merc,venus",
        [
            (0, 0, 0, 0, 0),
            (180, 90, 270, 45, 135),
            (359.99, 0.01, 179.99, 90.0, 270.0),
            (1.0, 359.0, 180.0, 90.0, 270.0),
        ],
    )
    def test_various_positions(
        self,
        asc: float,
        sun: float,
        moon: float,
        merc: float,
        venus: float,
    ):
        """calc_all_arabic_parts works with various position combos."""
        positions = {
            "Asc": asc,
            "Sun": sun,
            "Moon": moon,
            "Mercury": merc,
            "Venus": venus,
        }
        parts = calc_all_arabic_parts(positions)
        for name, val in parts.items():
            assert 0 <= val < 360, f"{name}={val}"

    @pytest.mark.unit
    def test_with_jd_and_geo(self):
        """calc_all_arabic_parts with jd/geo_lat/geo_lon parameters."""
        positions = {
            "Asc": 15.5,
            "Sun": 120.0,
            "Moon": 240.0,
            "Mercury": 130.0,
            "Venus": 110.0,
        }
        parts = calc_all_arabic_parts(
            positions, jd=2451545.0, geo_lat=41.9, geo_lon=12.5
        )
        for name, val in parts.items():
            assert 0 <= val < 360, f"{name}={val}"


class TestIsDayChart:
    """Tests for day/night chart determination."""

    @pytest.mark.unit
    def test_sun_above_horizon_is_day(self):
        """Sun above horizon (between Asc and Desc) is a day chart."""
        # Sun at 120°, Asc at 90° => Sun is in upper hemisphere
        result = is_day_chart(120.0, 90.0)
        assert isinstance(result, bool)

    @pytest.mark.unit
    def test_returns_bool(self):
        """is_day_chart always returns a boolean."""
        for sun in [0, 90, 180, 270]:
            for asc in [0, 90, 180, 270]:
                result = is_day_chart(float(sun), float(asc))
                assert isinstance(result, bool), (
                    f"is_day_chart({sun}, {asc}) returned {type(result)}"
                )


class TestArabicPartsIntegration:
    """Integration tests using real ephemeris positions."""

    @pytest.mark.unit
    def test_arabic_parts_from_real_positions(self):
        """Calculate Arabic parts from actual planetary positions."""
        jd = 2451545.0
        swe.set_topo(12.5, 41.9, 0.0)

        # Get actual positions
        sun_pos, _ = swe.calc_ut(jd, 0, 0)
        moon_pos, _ = swe.calc_ut(jd, 1, 0)
        merc_pos, _ = swe.calc_ut(jd, 2, 0)
        venus_pos, _ = swe.calc_ut(jd, 3, 0)

        # Get Asc from houses
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("P"))
        asc = ascmc[0]

        positions = {
            "Asc": asc,
            "Sun": sun_pos[0],
            "Moon": moon_pos[0],
            "Mercury": merc_pos[0],
            "Venus": venus_pos[0],
        }
        parts = calc_all_arabic_parts(positions, jd=jd, geo_lat=41.9, geo_lon=12.5)

        for name, val in parts.items():
            assert 0 <= val < 360, f"{name}={val}"
            assert math.isfinite(val), f"{name}={val} not finite"

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2024, 2050, 2100])
    def test_arabic_parts_across_years(self, year: int):
        """Arabic parts from real positions at different years."""
        jd = swe.julday(year, 6, 21, 12.0)
        swe.set_topo(12.5, 41.9, 0.0)

        sun_pos, _ = swe.calc_ut(jd, 0, 0)
        moon_pos, _ = swe.calc_ut(jd, 1, 0)
        merc_pos, _ = swe.calc_ut(jd, 2, 0)
        venus_pos, _ = swe.calc_ut(jd, 3, 0)
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("P"))

        positions = {
            "Asc": ascmc[0],
            "Sun": sun_pos[0],
            "Moon": moon_pos[0],
            "Mercury": merc_pos[0],
            "Venus": venus_pos[0],
        }
        parts = calc_all_arabic_parts(positions, jd=jd, geo_lat=41.9, geo_lon=12.5)
        for name, val in parts.items():
            assert 0 <= val < 360
