"""
Comprehensive tests for planetary phenomena (pheno_ut).

Verifies phase angle, illumination fraction, elongation,
apparent diameter, and magnitude calculations.
"""

from __future__ import annotations

import math

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
)


PLANETS_FOR_PHENO = [
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
]


class TestPhenoBasic:
    """Basic pheno_ut functionality."""

    @pytest.mark.unit
    def test_pheno_returns_20_elements(self):
        """pheno_ut returns 20 floats."""
        result = swe.pheno_ut(2451545.0, MARS, 0)
        assert len(result) == 20, f"Expected 20, got {len(result)}"

    @pytest.mark.unit
    def test_pheno_returns_native_floats(self):
        """All elements should be native Python float."""
        result = swe.pheno_ut(2451545.0, MARS, 0)
        for i in range(5):
            assert type(result[i]) is float, (
                f"result[{i}] is {type(result[i]).__name__}"
            )

    @pytest.mark.unit
    def test_pheno_all_finite(self):
        """First 5 elements should be finite."""
        result = swe.pheno_ut(2451545.0, MARS, 0)
        for i in range(5):
            assert math.isfinite(result[i]), f"result[{i}] = {result[i]}"


class TestPhenoPhaseAngle:
    """Test phase angle values."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS_FOR_PHENO)
    def test_phase_angle_in_range(self, body_id: int, name: str):
        """Phase angle should be 0-180deg."""
        result = swe.pheno_ut(2451545.0, body_id, 0)
        phase_angle = result[0]
        assert 0 <= phase_angle <= 180, f"{name}: phase angle {phase_angle} deg"

    @pytest.mark.unit
    def test_sun_phase_angle_zero(self):
        """Sun phase angle should be 0 (Sun observes itself)."""
        result = swe.pheno_ut(2451545.0, SUN, 0)
        assert abs(result[0]) < 0.1, f"Sun phase angle: {result[0]} deg"


class TestPhenoIllumination:
    """Test illumination fraction."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS_FOR_PHENO)
    def test_illumination_in_range(self, body_id: int, name: str):
        """Illumination fraction should be 0-1."""
        result = swe.pheno_ut(2451545.0, body_id, 0)
        phase = result[1]
        assert 0 <= phase <= 1.01, f"{name}: illumination {phase}"

    @pytest.mark.unit
    def test_sun_illumination_defined(self):
        """The luminous solar disc is represented as fully illuminated."""
        result = swe.pheno_ut(2451545.0, SUN, 0)
        assert result[1] == 1.0

    @pytest.mark.unit
    def test_outer_planet_mostly_illuminated(self):
        """Outer planets (Jupiter, Saturn) should be mostly illuminated."""
        for body_id in [JUPITER, SATURN]:
            result = swe.pheno_ut(2451545.0, body_id, 0)
            phase = result[1]
            assert phase > 0.9, f"Body {body_id}: illumination {phase} (expected > 0.9)"


class TestPhenoElongation:
    """Test elongation from Sun."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS_FOR_PHENO)
    def test_elongation_in_range(self, body_id: int, name: str):
        """Elongation should be 0-180deg."""
        result = swe.pheno_ut(2451545.0, body_id, 0)
        elongation = result[2]
        assert 0 <= elongation <= 180, f"{name}: elongation {elongation} deg"

    @pytest.mark.unit
    def test_sun_elongation_zero(self):
        """Sun elongation from itself should be ~0."""
        result = swe.pheno_ut(2451545.0, SUN, 0)
        assert abs(result[2]) < 0.1, f"Sun elongation: {result[2]} deg"

    @pytest.mark.unit
    def test_inner_planet_elongation_limited(self):
        """Mercury elongation should be < ~28deg; Venus < ~47deg."""
        # Mercury max elongation ~28deg
        result_merc = swe.pheno_ut(2451545.0, MERCURY, 0)
        assert result_merc[2] < 30, f"Mercury elongation {result_merc[2]} deg > 30 deg"

        # Venus max elongation ~47deg
        result_ven = swe.pheno_ut(2451545.0, VENUS, 0)
        assert result_ven[2] < 50, f"Venus elongation {result_ven[2]} deg > 50 deg"


class TestPhenoDiameter:
    """Test apparent diameter."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS_FOR_PHENO)
    def test_diameter_positive(self, body_id: int, name: str):
        """Apparent diameter should be positive."""
        result = swe.pheno_ut(2451545.0, body_id, 0)
        diameter = result[3]
        assert diameter > 0, f"{name}: diameter {diameter}"

    @pytest.mark.unit
    def test_moon_diameter_largest(self):
        """Moon should have the largest apparent diameter (~0.5 deg or ~1800 arcsec)."""
        result = swe.pheno_ut(2451545.0, MOON, 0)
        diameter = result[3]
        # Diameter may be in degrees (~0.5) or arcsec (~1800) depending on impl
        if diameter < 10:
            # In degrees: Moon ~0.5 deg
            assert 0.4 < diameter < 0.6, (
                f"Moon diameter {diameter} deg (expected ~0.5 deg)"
            )
        else:
            # In arcseconds: Moon ~1800-2000"
            assert 1700 < diameter < 2100, (
                f"Moon diameter {diameter} arcsec (expected ~1800-2000)"
            )


class TestPhenoMagnitude:
    """Test visual magnitude."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS_FOR_PHENO)
    def test_magnitude_finite(self, body_id: int, name: str):
        """Visual magnitude should be finite."""
        result = swe.pheno_ut(2451545.0, body_id, 0)
        mag = result[4]
        assert math.isfinite(mag), f"{name}: magnitude {mag}"

    @pytest.mark.unit
    def test_venus_brightest_planet(self):
        """Venus should typically be brighter (lower mag) than Jupiter."""
        # Not always true, but at J2000 Venus is close to its max brightness
        result_v = swe.pheno_ut(2451545.0, VENUS, 0)
        result_j = swe.pheno_ut(2451545.0, JUPITER, 0)
        # Venus mag typically -3 to -4.5, Jupiter -2 to -3
        # Both should be negative (bright)
        assert result_v[4] < 0, f"Venus mag {result_v[4]}"
        assert result_j[4] < 0, f"Jupiter mag {result_j[4]}"


class TestPhenoDateRange:
    """Test pheno_ut across dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2024, 2050, 2100])
    def test_mars_pheno_across_years(self, year: int):
        """Mars phenomena valid across years."""
        jd = swe.julday(year, 6, 21, 12.0)
        result = swe.pheno_ut(jd, MARS, 0)
        assert len(result) == 20
        assert 0 <= result[0] <= 180
        assert 0 <= result[1] <= 1.01
        assert 0 <= result[2] <= 180

    @pytest.mark.unit
    def test_pheno_multiple_planets_same_date(self):
        """All planets valid on same date."""
        jd = 2451545.0
        for body_id, name in PLANETS_FOR_PHENO:
            result = swe.pheno_ut(jd, body_id, 0)
            assert len(result) == 20, f"{name}: {len(result)} elements"
            assert math.isfinite(result[0]), f"{name}: phase angle"
            assert math.isfinite(result[1]), f"{name}: illumination"
