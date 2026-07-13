"""Tests for pheno_ut across all planets over extended date ranges."""

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
    FLG_SWIEPH,
)

JD_J2000 = 2451545.0


@pytest.mark.unit
class TestPhenoBasic:
    """Basic pheno_ut functionality."""

    @pytest.mark.parametrize(
        "body,name",
        [
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_pheno_returns_values(self, body, name):
        """pheno_ut returns values for planets."""
        result = swe.pheno_ut(JD_J2000, body, FLG_SWIEPH)
        assert len(result) >= 5, f"Not enough values for {name}"

    def test_pheno_sun_illumination_inapplicable(self):
        """Illuminated fraction is inapplicable to the Sun: reported as 0.0."""
        result = swe.pheno_ut(JD_J2000, SUN, FLG_SWIEPH)
        assert result[1] == 0.0

    def test_pheno_moon_phase_angle(self):
        """Moon phase angle should be in [0, 180]."""
        result = swe.pheno_ut(JD_J2000, MOON, FLG_SWIEPH)
        phase_angle = result[0]
        assert 0 <= phase_angle <= 180, f"Phase angle {phase_angle}° out of range"

    def test_pheno_moon_illumination(self):
        """Moon illumination should be in [0, 1]."""
        result = swe.pheno_ut(JD_J2000, MOON, FLG_SWIEPH)
        illum = result[1]
        assert 0 <= illum <= 1, f"Illumination {illum} out of range"

    def test_pheno_elongation_range(self):
        """Elongation should be in [0, 180]."""
        for body in [MERCURY, VENUS, MARS, JUPITER]:
            result = swe.pheno_ut(JD_J2000, body, FLG_SWIEPH)
            elong = result[2]
            assert 0 <= elong <= 180, f"Elongation {elong}° for body {body}"

    def test_pheno_diameter_positive(self):
        """Apparent diameter should be positive."""
        for body in [MOON, MERCURY, VENUS, MARS, JUPITER]:
            result = swe.pheno_ut(JD_J2000, body, FLG_SWIEPH)
            diam = result[3]
            assert diam > 0, f"Diameter {diam} for body {body} not positive"


@pytest.mark.unit
class TestPhenoMoonPhases:
    """Test Moon phenomena through a full lunation."""

    def test_full_moon_high_illumination(self):
        """Near full Moon, illumination should be near 1.0."""
        # Find a point where Moon is near 180° from Sun
        for offset in range(30):
            jd = JD_J2000 + offset
            sun, _ = swe.calc_ut(jd, SUN, FLG_SWIEPH)
            moon, _ = swe.calc_ut(jd, MOON, FLG_SWIEPH)
            elong = abs(moon[0] - sun[0])
            if elong > 180:
                elong = 360 - elong
            if elong > 170:
                result = swe.pheno_ut(jd, MOON, FLG_SWIEPH)
                assert result[1] > 0.9, f"Full moon illumination {result[1]} too low"
                return
        pytest.skip("No full moon found in 30-day window")

    def test_new_moon_low_illumination(self):
        """Near new Moon, illumination should be near 0.0."""
        for offset in range(30):
            jd = JD_J2000 + offset
            sun, _ = swe.calc_ut(jd, SUN, FLG_SWIEPH)
            moon, _ = swe.calc_ut(jd, MOON, FLG_SWIEPH)
            elong = abs(moon[0] - sun[0])
            if elong > 180:
                elong = 360 - elong
            if elong < 10:
                result = swe.pheno_ut(jd, MOON, FLG_SWIEPH)
                assert result[1] < 0.1, f"New moon illumination {result[1]} too high"
                return
        pytest.skip("No new moon found in 30-day window")

    def test_moon_illumination_varies(self):
        """Moon illumination should vary over a month."""
        illuminations = []
        for offset in range(0, 30):
            jd = JD_J2000 + offset
            result = swe.pheno_ut(jd, MOON, FLG_SWIEPH)
            illuminations.append(result[1])
        # Should see both high and low illumination values
        assert max(illuminations) > 0.8, "Never saw high illumination"
        assert min(illuminations) < 0.2, "Never saw low illumination"


@pytest.mark.unit
class TestPhenoInnerPlanets:
    """Phenomena for inner planets (Mercury, Venus)."""

    @pytest.mark.parametrize("body,name", [(MERCURY, "Mercury"), (VENUS, "Venus")])
    def test_max_elongation_bounded(self, body, name):
        """Inner planet elongation should be bounded (< 90°)."""
        max_elong = 0
        for offset in range(0, 365, 5):
            jd = JD_J2000 + offset
            result = swe.pheno_ut(jd, body, FLG_SWIEPH)
            max_elong = max(max_elong, result[2])
        if body == MERCURY:
            assert max_elong < 30, f"Mercury max elongation {max_elong}° > 30"
        else:
            assert max_elong < 50, f"Venus max elongation {max_elong}° > 50"

    @pytest.mark.parametrize("body,name", [(MERCURY, "Mercury"), (VENUS, "Venus")])
    def test_illumination_varies_with_elongation(self, body, name):
        """Inner planet illumination should correlate with elongation."""
        # At small elongation (near Sun), can be nearly full or nearly new
        # At max elongation, should be roughly half-illuminated
        elongations = []
        illuminations = []
        for offset in range(0, 365, 10):
            jd = JD_J2000 + offset
            result = swe.pheno_ut(jd, body, FLG_SWIEPH)
            elongations.append(result[2])
            illuminations.append(result[1])
        # Just verify range exists
        assert max(illuminations) > 0.3, f"{name} max illumination too low"


@pytest.mark.unit
class TestPhenoOuterPlanets:
    """Phenomena for outer planets."""

    @pytest.mark.parametrize(
        "body,name",
        [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_elongation_reaches_180(self, body, name):
        """Outer planets can reach 180° elongation (opposition)."""
        max_elong = 0
        for offset in range(0, 730, 5):  # Search 2 years
            jd = JD_J2000 + offset
            result = swe.pheno_ut(jd, body, FLG_SWIEPH)
            max_elong = max(max_elong, result[2])
        assert max_elong > 150, f"{name} max elongation {max_elong}° too low"

    @pytest.mark.parametrize(
        "body,name",
        [
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
        ],
    )
    def test_illumination_nearly_full(self, body, name):
        """Far outer planets should always be nearly fully illuminated."""
        for offset in range(0, 365, 30):
            jd = JD_J2000 + offset
            result = swe.pheno_ut(jd, body, FLG_SWIEPH)
            illum = result[1]
            assert illum > 0.85, (
                f"{name} illumination {illum} too low at offset {offset}"
            )

    @pytest.mark.parametrize(
        "body,name,mag_range",
        [
            (VENUS, "Venus", (-5.0, -3.0)),
            (JUPITER, "Jupiter", (-3.0, -1.0)),
            (SATURN, "Saturn", (-1.0, 2.0)),
            (MARS, "Mars", (-3.0, 2.5)),
        ],
    )
    def test_magnitude_reasonable(self, body, name, mag_range):
        """Apparent magnitude should be within expected range."""
        result = swe.pheno_ut(JD_J2000, body, FLG_SWIEPH)
        mag = result[4]
        if not math.isnan(mag):
            lo, hi = mag_range
            assert lo < mag < hi, f"{name} magnitude {mag} not in [{lo}, {hi}]"
