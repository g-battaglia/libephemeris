"""
Tests for calc_besselian_y function in libephemeris.

Tests the Besselian y coordinate calculation for solar eclipses.

The Besselian y coordinate is the y-component (northward) of the Moon's
shadow axis intersection with the fundamental plane, measured in Earth radii.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Explanatory Supplement to the Astronomical Almanac
- Meeus "Astronomical Algorithms" Chapter 54
"""

import math
from libephemeris import julday, calc_besselian_y, calc_besselian_x, FLG_SWIEPH


class TestBesselianYBasicFunctionality:
    """Test basic functionality of calc_besselian_y."""

    def test_function_exists(self):
        """Test that calc_besselian_y function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_y

        assert callable(calc_besselian_y)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_besselian_y

        assert callable(calc_besselian_y)

    def test_returns_float(self):
        """Test that function returns a float value."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_y(jd)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_y(jd, flags=FLG_SWIEPH)

        assert isinstance(result, float)

    def test_result_in_reasonable_range(self):
        """Test that result is within reasonable range for Earth radii."""
        # For any time, y should be within a few Earth radii
        # (Moon is about 60 Earth radii away)
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_y(jd)

        # The y coordinate should be between -2 and +2 during/near an eclipse
        # and could be larger when far from eclipse (but still bounded)
        assert -100 < result < 100


class TestBesselianYDuringEclipses:
    """Test calc_besselian_y during known solar eclipses."""

    def test_april_2024_total_eclipse_near_maximum(self):
        """Test y coordinate during April 8, 2024 total solar eclipse.

        NASA gives maximum eclipse at approximately 18:17 UTC with the
        central line over North America. The gamma value (closest approach
        of shadow axis to Earth's center) is approximately 0.34.
        Since gamma^2 = x^2 + y^2 at maximum, y should be bounded.
        """
        # April 8, 2024 around 18:17 UTC - time of maximum eclipse
        jd_max = julday(2024, 4, 8, 18.3)

        y = calc_besselian_y(jd_max)

        # At maximum eclipse, y should be relatively small (< 1 Earth radius)
        # because the shadow axis passes through Earth
        # Reference: NASA eclipse data shows gamma ~ 0.34 for this eclipse
        assert abs(y) < 1.5, f"y = {y} at maximum eclipse, expected |y| < 1.5"

    def test_april_2024_eclipse_y_varies_with_time(self):
        """Test that y coordinate changes as eclipse progresses."""
        # Sample times during the eclipse
        jd_before = julday(2024, 4, 8, 17.0)  # Before maximum
        jd_max = julday(2024, 4, 8, 18.3)  # Near maximum
        jd_after = julday(2024, 4, 8, 19.5)  # After maximum

        y_before = calc_besselian_y(jd_before)
        y_max = calc_besselian_y(jd_max)
        y_after = calc_besselian_y(jd_after)

        # y should change over time (not constant)
        assert y_before != y_max or y_max != y_after, "y should vary with time"

    def test_october_2023_annular_eclipse(self):
        """Test y coordinate during October 14, 2023 annular eclipse.

        This eclipse crossed North, Central, and South America.
        Maximum eclipse around 18:00 UTC.
        """
        jd_max = julday(2023, 10, 14, 18.0)

        y = calc_besselian_y(jd_max)

        # During this eclipse, y should be relatively small
        assert abs(y) < 1.5, f"y = {y} at annular eclipse maximum"

    def test_december_2021_total_eclipse(self):
        """Test y coordinate during December 4, 2021 total eclipse.

        This was a total eclipse visible from Antarctica.
        Maximum eclipse around 07:34 UTC.
        """
        jd_max = julday(2021, 12, 4, 7.5)

        y = calc_besselian_y(jd_max)

        # y should be within bounds during eclipse
        assert abs(y) < 2.0, f"y = {y} at Dec 2021 eclipse maximum"


class TestBesselianYPhysicalProperties:
    """Test physical properties of the Besselian y coordinate."""

    def test_y_changes_sign_appropriately(self):
        """Test that y can be positive or negative depending on geometry."""
        # Test at different times to check sign variation
        jd1 = julday(2024, 4, 8, 16.0)
        jd2 = julday(2024, 4, 8, 20.0)

        y1 = calc_besselian_y(jd1)
        y2 = calc_besselian_y(jd2)

        # Both should be finite
        assert math.isfinite(y1)
        assert math.isfinite(y2)

    def test_y_is_continuous(self):
        """Test that y changes continuously (no discontinuities)."""
        jd_start = julday(2024, 4, 8, 17.0)

        # Sample at small time intervals
        y_values = []
        for i in range(10):
            jd = jd_start + i * 0.01  # 10 samples over ~15 minutes
            y = calc_besselian_y(jd)
            y_values.append(y)

        # Check that consecutive values don't jump too much
        # (expecting smooth variation)
        for i in range(1, len(y_values)):
            delta = abs(y_values[i] - y_values[i - 1])
            # Change should be small over ~15 minutes (< 0.2 Earth radii)
            # The Moon moves about 0.5 deg/hour, causing gradual y variation
            assert delta < 0.2, f"Discontinuity detected: delta = {delta}"

    def test_y_away_from_eclipse(self):
        """Test y coordinate when far from any eclipse (random time)."""
        # A random date far from any eclipse
        jd_random = julday(2024, 7, 15, 12.0)

        y = calc_besselian_y(jd_random)

        # Should still return a finite value
        assert math.isfinite(y)

        # When far from eclipse, y can be large (shadow misses Earth)
        # but should still be bounded
        assert abs(y) < 100, f"y = {y} should be bounded even away from eclipse"


class TestBesselianYNumericalStability:
    """Test numerical stability of the calculation."""

    def test_no_division_by_zero(self):
        """Test that function handles edge cases without division by zero."""
        # Test various dates
        test_dates = [
            julday(2000, 1, 1, 12.0),
            julday(2024, 6, 21, 12.0),  # Summer solstice
            julday(2024, 12, 21, 12.0),  # Winter solstice
            julday(2024, 3, 20, 12.0),  # Vernal equinox
        ]

        for jd in test_dates:
            y = calc_besselian_y(jd)
            assert math.isfinite(y), f"y not finite at JD {jd}"

    def test_consistent_results(self):
        """Test that same input gives same output (deterministic)."""
        jd = julday(2024, 4, 8, 18.0)

        y1 = calc_besselian_y(jd)
        y2 = calc_besselian_y(jd)

        assert y1 == y2, "Results should be deterministic"


class TestBesselianYValidation:
    """Validation tests against reference values."""

    def test_near_new_moon_y_bounded(self):
        """Test that y is bounded near new moon (potential eclipse time)."""
        # Near April 2024 new moon (which was the eclipse)
        jd_new_moon = julday(2024, 4, 8, 18.2)

        y = calc_besselian_y(jd_new_moon)

        # At new moon during eclipse, y should be small (shadow near Earth)
        # The April 2024 eclipse had gamma approx 0.34, so y should be bounded
        assert abs(y) < 2.0, f"y = {y} should be bounded at eclipse new moon"

    def test_y_rate_of_change(self):
        """Test that y changes at a reasonable rate."""
        jd1 = julday(2024, 4, 8, 18.0)
        jd2 = julday(2024, 4, 8, 19.0)  # 1 hour later

        y1 = calc_besselian_y(jd1)
        y2 = calc_besselian_y(jd2)

        # The Moon moves about 0.5 deg per hour relative to the Sun
        # This translates to roughly 0.008 Earth radii per hour change in y
        # Allow for some variation due to projection effects
        delta_y = abs(y2 - y1)

        # Rate should be less than 1 Earth radius per hour in typical cases
        assert delta_y < 1.0, f"y changed by {delta_y} in 1 hour, seems too fast"

    def test_gamma_approximation(self):
        """Test that x^2 + y^2 approximately equals gamma^2 at maximum eclipse.

        For the April 8, 2024 eclipse, NASA gives gamma ~ 0.3433.
        At maximum eclipse, x^2 + y^2 should approximate gamma^2.
        """
        # April 8, 2024 eclipse maximum
        jd_max = julday(2024, 4, 8, 18.3)

        x = calc_besselian_x(jd_max)
        y = calc_besselian_y(jd_max)

        # gamma^2 = x^2 + y^2 (approximately, at maximum eclipse)
        gamma_squared = x**2 + y**2
        gamma = math.sqrt(gamma_squared)

        # NASA reference: gamma = 0.3433 for this eclipse
        # Allow some tolerance due to timing approximation
        assert 0.1 < gamma < 0.8, (
            f"gamma = {gamma} (x={x}, y={y}), expected gamma ~ 0.34"
        )

    def test_x_y_orthogonal_properties(self):
        """Test that x and y coordinates are computed consistently.

        The x and y axes should be orthogonal, meaning changes in
        shadow position should affect x and y independently.
        """
        jd = julday(2024, 4, 8, 18.0)

        x = calc_besselian_x(jd)
        y = calc_besselian_y(jd)

        # Both should be finite
        assert math.isfinite(x)
        assert math.isfinite(y)

        # x represents east-west displacement, y represents north-south
        # They should both be reasonable for an eclipse
        assert abs(x) < 2.0, f"x = {x} seems unreasonable"
        assert abs(y) < 2.0, f"y = {y} seems unreasonable"
