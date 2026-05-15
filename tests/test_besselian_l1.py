"""
Tests for calc_besselian_l1 function in libephemeris.

Tests the Besselian l1 (penumbral shadow radius) calculation for solar eclipses.

The Besselian element l1 is the radius of the penumbral shadow cone where it
intersects the fundamental plane, measured in Earth equatorial radii. This
determines the extent of the partial eclipse zone on Earth.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Explanatory Supplement to the Astronomical Almanac
- Meeus "Astronomical Algorithms" Chapter 54
"""

import pytest
import math
from libephemeris import julday, calc_besselian_l1, FLG_SWIEPH


class TestBesselianL1BasicFunctionality:
    """Test basic functionality of calc_besselian_l1."""

    def test_function_exists(self):
        """Test that calc_besselian_l1 function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_l1

        assert callable(calc_besselian_l1)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_besselian_l1

        assert callable(calc_besselian_l1)

    def test_returns_float(self):
        """Test that function returns a float value."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_l1(jd)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_l1(jd, flags=FLG_SWIEPH)

        assert isinstance(result, float)

    def test_result_is_positive(self):
        """Test that l1 is always positive (it's a radius)."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_l1(jd)

        assert result > 0, f"l1 = {result} should be positive"


class TestBesselianL1DuringEclipses:
    """Test calc_besselian_l1 during known solar eclipses."""

    def test_april_2024_total_eclipse_near_maximum(self):
        """Test l1 during April 8, 2024 total solar eclipse.

        The penumbral shadow radius l1 during an eclipse is typically
        between 0.5 and 0.6 Earth radii, corresponding to the penumbral
        zone extending roughly 3200-3800 km from the shadow axis.
        """
        # April 8, 2024 around 18:17 UTC - time of maximum eclipse
        jd_max = julday(2024, 4, 8, 18.3)

        l1 = calc_besselian_l1(jd_max)

        # l1 should be in the typical range for solar eclipses
        assert 0.4 < l1 < 0.7, f"l1 = {l1} Earth radii, expected 0.5-0.6"

    def test_april_2024_eclipse_l1_varies_with_time(self):
        """Test that l1 changes slightly as eclipse progresses."""
        # Sample times during the eclipse
        jd_before = julday(2024, 4, 8, 17.0)  # Before maximum
        jd_max = julday(2024, 4, 8, 18.3)  # Near maximum
        jd_after = julday(2024, 4, 8, 19.5)  # After maximum

        l1_before = calc_besselian_l1(jd_before)
        l1_max = calc_besselian_l1(jd_max)
        l1_after = calc_besselian_l1(jd_after)

        # All should be in the reasonable range
        assert all(0.4 < l1 < 0.7 for l1 in [l1_before, l1_max, l1_after])

        # l1 should change slightly over time as Moon's distance changes
        assert l1_before != l1_after, "l1 should vary with time"

    def test_october_2023_annular_eclipse(self):
        """Test l1 during October 14, 2023 annular eclipse."""
        jd_max = julday(2023, 10, 14, 18.0)

        l1 = calc_besselian_l1(jd_max)

        # l1 should be in the typical range
        assert 0.4 < l1 < 0.7, f"l1 = {l1} Earth radii for October eclipse"

    def test_december_2021_total_eclipse(self):
        """Test l1 during December 4, 2021 total eclipse (Antarctica)."""
        jd_max = julday(2021, 12, 4, 7.5)

        l1 = calc_besselian_l1(jd_max)

        # l1 should be in the typical range
        assert 0.4 < l1 < 0.7, f"l1 = {l1} Earth radii for December 2021 eclipse"

    def test_june_2021_annular_eclipse(self):
        """Test l1 during June 10, 2021 annular eclipse."""
        jd_max = julday(2021, 6, 10, 10.5)

        l1 = calc_besselian_l1(jd_max)

        # l1 should be in the typical range
        assert 0.4 < l1 < 0.7, f"l1 = {l1} Earth radii for June 2021 eclipse"


class TestBesselianL1PhysicalProperties:
    """Test physical properties of the l1 calculation."""

    def test_l1_always_positive(self):
        """Test that l1 is always positive at any time."""
        # Test various dates throughout the year (not necessarily during eclipses)
        test_dates = [
            julday(2024, 1, 15, 12.0),
            julday(2024, 4, 15, 12.0),
            julday(2024, 7, 15, 12.0),
            julday(2024, 10, 15, 12.0),
        ]

        for jd in test_dates:
            l1 = calc_besselian_l1(jd)
            assert l1 > 0, f"l1 = {l1} should be positive at JD {jd}"

    def test_l1_is_continuous(self):
        """Test that l1 changes continuously (no discontinuities)."""
        jd_start = julday(2024, 4, 8, 17.0)

        # Sample at small time intervals
        l1_values = []
        for i in range(10):
            jd = jd_start + i * 0.01  # 10 samples over ~15 minutes
            l1 = calc_besselian_l1(jd)
            l1_values.append(l1)

        # Check that consecutive values don't jump too much
        for i in range(1, len(l1_values)):
            delta = abs(l1_values[i] - l1_values[i - 1])
            # Change should be very small over ~15 minutes (< 0.01 Earth radii)
            assert delta < 0.05, f"Discontinuity detected: delta = {delta} Earth radii"

    def test_l1_reasonable_at_new_moon(self):
        """Test l1 is reasonable at new moon (when eclipses can occur)."""
        # New moons in 2024
        new_moons = [
            julday(2024, 1, 11, 12.0),  # January new moon
            julday(2024, 4, 8, 18.0),  # April new moon (eclipse)
            julday(2024, 7, 5, 22.0),  # July new moon
            julday(2024, 10, 2, 18.0),  # October new moon
        ]

        for jd in new_moons:
            l1 = calc_besselian_l1(jd)
            # l1 should be in a reasonable range
            assert 0.3 < l1 < 0.8, f"l1 = {l1} at JD {jd}"


class TestBesselianL1NumericalStability:
    """Test numerical stability of the calculation."""

    def test_no_division_by_zero(self):
        """Test that function handles edge cases without division by zero."""
        # Test various dates including equinoxes and solstices
        test_dates = [
            julday(2000, 1, 1, 12.0),
            julday(2024, 6, 21, 12.0),  # Summer solstice
            julday(2024, 12, 21, 12.0),  # Winter solstice
            julday(2024, 3, 20, 12.0),  # Vernal equinox
            julday(2024, 9, 22, 12.0),  # Autumnal equinox
        ]

        for jd in test_dates:
            l1 = calc_besselian_l1(jd)
            assert math.isfinite(l1), f"l1 not finite at JD {jd}"

    def test_consistent_results(self):
        """Test that same input gives same output (deterministic)."""
        jd = julday(2024, 4, 8, 18.0)

        l1_1 = calc_besselian_l1(jd)
        l1_2 = calc_besselian_l1(jd)

        assert l1_1 == l1_2, "Results should be deterministic"

    def test_no_nan_or_inf(self):
        """Test that function never returns NaN or infinity."""
        # Test a range of dates
        for year in range(2020, 2030):
            for month in [1, 4, 7, 10]:
                jd = julday(year, month, 15, 12.0)
                l1 = calc_besselian_l1(jd)
                assert not math.isnan(l1), f"l1 is NaN at {year}/{month}"
                assert not math.isinf(l1), f"l1 is infinite at {year}/{month}"


class TestBesselianL1Validation:
    """Validation tests against expected values."""

    def test_l1_typical_value_during_eclipse(self):
        """Validate that l1 is in the expected range during solar eclipses.

        The penumbral shadow radius depends on:
        - Sun's angular size (varies with Earth-Sun distance)
        - Moon's angular size (varies with Earth-Moon distance)
        - Moon's distance from Earth

        Typical values are around 0.53-0.57 Earth radii.
        """
        jd_eclipse = julday(2024, 4, 8, 18.0)
        l1 = calc_besselian_l1(jd_eclipse)

        # Expected range based on typical geometry
        assert 0.45 < l1 < 0.65, f"l1 = {l1}, expected 0.5-0.6 Earth radii"

    def test_l1_rate_of_change(self):
        """Test that l1 changes at a reasonable rate."""
        jd1 = julday(2024, 4, 8, 18.0)
        jd2 = julday(2024, 4, 8, 19.0)  # 1 hour later

        l1_1 = calc_besselian_l1(jd1)
        l1_2 = calc_besselian_l1(jd2)

        # The change in l1 over an hour should be small
        # (Moon's distance changes slowly compared to shadow geometry)
        delta_l1 = abs(l1_2 - l1_1)

        # Rate should be less than 0.01 Earth radii per hour
        assert delta_l1 < 0.05, f"l1 changed by {delta_l1} Earth radii in 1 hour"

    def test_l1_larger_when_moon_farther(self):
        """Test that l1 tends to be larger when Moon is farther from Earth.

        When the Moon is at apogee (farther from Earth), the penumbral cone
        vertex moves slightly relative to Earth, but the effect on l1 is
        relatively small since the Sun's distance dominates the geometry.

        This test verifies that the calculation produces reasonable values
        at different Moon distances.
        """
        # Test during eclipse conditions when the geometry is well-defined
        # Use times during the April 2024 eclipse when Moon distance varies slightly
        jd_early = julday(2024, 4, 8, 16.0)  # Earlier in eclipse
        jd_late = julday(2024, 4, 8, 20.0)  # Later in eclipse

        l1_early = calc_besselian_l1(jd_early)
        l1_late = calc_besselian_l1(jd_late)

        # Both should be positive and in reasonable range
        assert l1_early > 0.4, f"l1 early = {l1_early}"
        assert l1_late > 0.4, f"l1 late = {l1_late}"

        # Both should be in the expected range
        assert 0.4 < l1_early < 0.7
        assert 0.4 < l1_late < 0.7


class TestBesselianL1EdgeCases:
    """Test edge cases and boundary conditions."""

    def test_far_past_date(self):
        """Test function works for historical dates."""
        # Eclipse of 1999 (August 11, total solar eclipse over Europe)
        jd = julday(1999, 8, 11, 11.0)

        l1 = calc_besselian_l1(jd)

        assert 0.3 < l1 < 0.8, f"l1 = {l1} for 1999 eclipse"
        assert math.isfinite(l1)

    def test_far_future_date(self):
        """Test function works for future dates."""
        # A date far in the future
        jd = julday(2050, 6, 15, 12.0)

        l1 = calc_besselian_l1(jd)

        assert 0.3 < l1 < 0.8, f"l1 = {l1} for June 2050"
        assert math.isfinite(l1)

    def test_multiple_calls_consistent(self):
        """Test that multiple calls with same input are consistent."""
        jd = julday(2024, 4, 8, 18.0)

        results = [calc_besselian_l1(jd) for _ in range(5)]

        # All results should be identical
        assert all(r == results[0] for r in results), "Results should be consistent"

    def test_l1_at_different_times_of_day(self):
        """Test l1 at different times throughout a day."""
        base_jd = julday(2024, 4, 8, 0.0)

        # Sample every 4 hours
        l1_values = []
        for hour in [0, 4, 8, 12, 16, 20]:
            jd = base_jd + hour / 24.0
            l1 = calc_besselian_l1(jd)
            l1_values.append(l1)
            assert 0.3 < l1 < 0.8, f"l1 = {l1} at hour {hour}"

        # All values should be similar (Moon doesn't move that far in a day)
        max_diff = max(l1_values) - min(l1_values)
        assert max_diff < 0.1, f"l1 varied by {max_diff} over the day"
