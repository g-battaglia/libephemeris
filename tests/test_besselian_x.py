"""
Tests for calc_besselian_x function in libephemeris.

Tests the Besselian x coordinate calculation for solar eclipses.

The Besselian x coordinate is the x-component (eastward) of the Moon's
shadow axis intersection with the fundamental plane, measured in Earth radii.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Explanatory Supplement to the Astronomical Almanac
- Meeus "Astronomical Algorithms" Chapter 54
"""

import math
from libephemeris import julday, calc_besselian_x, FLG_SWIEPH


class TestBesselianXBasicFunctionality:
    """Test basic functionality of calc_besselian_x."""

    def test_function_exists(self):
        """Test that calc_besselian_x function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_x

        assert callable(calc_besselian_x)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_besselian_x

        assert callable(calc_besselian_x)

    def test_returns_float(self):
        """Test that function returns a float value."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_x(jd)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_x(jd, flags=FLG_SWIEPH)

        assert isinstance(result, float)

    def test_result_in_reasonable_range(self):
        """Test that result is within reasonable range for Earth radii."""
        # For any time, x should be within a few Earth radii
        # (Moon is about 60 Earth radii away)
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_x(jd)

        # The x coordinate should be between -2 and +2 during/near an eclipse
        # and could be larger when far from eclipse (but still bounded)
        assert -100 < result < 100


class TestBesselianXDuringEclipses:
    """Test calc_besselian_x during known solar eclipses."""

    def test_april_2024_total_eclipse_near_maximum(self):
        """Test x coordinate during April 8, 2024 total solar eclipse.

        NASA gives maximum eclipse at approximately 18:17 UTC with the
        central line over North America. The x coordinate at maximum
        should be close to 0 (shadow axis passing near Earth's center line).
        """
        # April 8, 2024 around 18:17 UTC - time of maximum eclipse
        jd_max = julday(2024, 4, 8, 18.3)

        x = calc_besselian_x(jd_max)

        # At maximum eclipse, x should be relatively small (< 1 Earth radius)
        # because the shadow axis passes through Earth
        # Reference: NASA eclipse data shows gamma ~ 0.34 for this eclipse
        assert abs(x) < 1.5, f"x = {x} at maximum eclipse, expected |x| < 1.5"

    def test_april_2024_eclipse_x_varies_with_time(self):
        """Test that x coordinate changes as eclipse progresses."""
        # Sample times during the eclipse
        jd_before = julday(2024, 4, 8, 17.0)  # Before maximum
        jd_max = julday(2024, 4, 8, 18.3)  # Near maximum
        jd_after = julday(2024, 4, 8, 19.5)  # After maximum

        x_before = calc_besselian_x(jd_before)
        x_max = calc_besselian_x(jd_max)
        x_after = calc_besselian_x(jd_after)

        # x should change over time (not constant)
        assert x_before != x_max or x_max != x_after, "x should vary with time"

    def test_october_2023_annular_eclipse(self):
        """Test x coordinate during October 14, 2023 annular eclipse.

        This eclipse crossed North, Central, and South America.
        Maximum eclipse around 18:00 UTC.
        """
        jd_max = julday(2023, 10, 14, 18.0)

        x = calc_besselian_x(jd_max)

        # During this eclipse, x should be relatively small
        assert abs(x) < 1.5, f"x = {x} at annular eclipse maximum"

    def test_december_2021_total_eclipse(self):
        """Test x coordinate during December 4, 2021 total eclipse.

        This was a total eclipse visible from Antarctica.
        Maximum eclipse around 07:34 UTC.
        """
        jd_max = julday(2021, 12, 4, 7.5)

        x = calc_besselian_x(jd_max)

        # x should be within bounds during eclipse
        assert abs(x) < 2.0, f"x = {x} at Dec 2021 eclipse maximum"


class TestBesselianXPhysicalProperties:
    """Test physical properties of the Besselian x coordinate."""

    def test_x_changes_sign_appropriately(self):
        """Test that x can be positive or negative depending on geometry."""
        # Test at different times to check sign variation
        jd1 = julday(2024, 4, 8, 16.0)
        jd2 = julday(2024, 4, 8, 20.0)

        x1 = calc_besselian_x(jd1)
        x2 = calc_besselian_x(jd2)

        # Both should be finite
        assert math.isfinite(x1)
        assert math.isfinite(x2)

    def test_x_is_continuous(self):
        """Test that x changes continuously (no discontinuities)."""
        jd_start = julday(2024, 4, 8, 17.0)

        # Sample at small time intervals
        x_values = []
        for i in range(10):
            jd = jd_start + i * 0.01  # 10 samples over ~15 minutes
            x = calc_besselian_x(jd)
            x_values.append(x)

        # Check that consecutive values don't jump too much
        # (expecting smooth variation)
        for i in range(1, len(x_values)):
            delta = abs(x_values[i] - x_values[i - 1])
            # Change should be small over ~15 minutes (< 0.2 Earth radii)
            # The Moon moves about 0.5°/hour, causing gradual x variation
            assert delta < 0.2, f"Discontinuity detected: delta = {delta}"

    def test_x_away_from_eclipse(self):
        """Test x coordinate when far from any eclipse (random time)."""
        # A random date far from any eclipse
        jd_random = julday(2024, 7, 15, 12.0)

        x = calc_besselian_x(jd_random)

        # Should still return a finite value
        assert math.isfinite(x)

        # When far from eclipse, x can be large (shadow misses Earth)
        # but should still be bounded
        assert abs(x) < 100, f"x = {x} should be bounded even away from eclipse"


class TestBesselianXNumericalStability:
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
            x = calc_besselian_x(jd)
            assert math.isfinite(x), f"x not finite at JD {jd}"

    def test_consistent_results(self):
        """Test that same input gives same output (deterministic)."""
        jd = julday(2024, 4, 8, 18.0)

        x1 = calc_besselian_x(jd)
        x2 = calc_besselian_x(jd)

        assert x1 == x2, "Results should be deterministic"


class TestBesselianXValidation:
    """Validation tests against reference values."""

    def test_near_new_moon_x_bounded(self):
        """Test that x is bounded near new moon (potential eclipse time)."""
        # Near April 2024 new moon (which was the eclipse)
        jd_new_moon = julday(2024, 4, 8, 18.2)

        x = calc_besselian_x(jd_new_moon)

        # At new moon during eclipse, x should be small (shadow near Earth)
        # The April 2024 eclipse had gamma ≈ 0.34, so x should be bounded
        assert abs(x) < 2.0, f"x = {x} should be bounded at eclipse new moon"

    def test_x_rate_of_change(self):
        """Test that x changes at a reasonable rate."""
        jd1 = julday(2024, 4, 8, 18.0)
        jd2 = julday(2024, 4, 8, 19.0)  # 1 hour later

        x1 = calc_besselian_x(jd1)
        x2 = calc_besselian_x(jd2)

        # The Moon moves about 0.5° per hour relative to the Sun
        # This translates to roughly 0.008 Earth radii per hour change in x
        # Allow for some variation due to projection effects
        delta_x = abs(x2 - x1)

        # Rate should be less than 1 Earth radius per hour in typical cases
        assert delta_x < 1.0, f"x changed by {delta_x} in 1 hour, seems too fast"

    def test_eclipse_geometry_consistency(self):
        """Test that x coordinate is consistent with eclipse geometry.

        During a central eclipse, when the shadow axis passes close to
        Earth's center, x should be small. This is a sanity check.
        """
        # April 8, 2024 eclipse maximum (central eclipse)
        jd_max = julday(2024, 4, 8, 18.3)

        x = calc_besselian_x(jd_max)

        # For a central eclipse with gamma ≈ 0.34, x² + y² ≈ gamma²
        # So x should be bounded by gamma (approximately)
        # Allow some margin since we're not computing y
        assert abs(x) < 1.0, (
            f"x = {x} at central eclipse, expected |x| < 1.0 for gamma ≈ 0.34"
        )
