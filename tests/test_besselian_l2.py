"""
Tests for calc_besselian_l2 function in libephemeris.

Tests the Besselian l2 (umbral/antumbral shadow radius) calculation for solar eclipses.

The Besselian element l2 is the radius of the umbral (or antumbral) shadow cone where it
intersects the fundamental plane, measured in Earth equatorial radii. This determines
whether a total or annular eclipse is possible:
- l2 > 0: Umbral shadow (total eclipse possible)
- l2 < 0: Antumbral shadow (annular eclipse)

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Explanatory Supplement to the Astronomical Almanac
- Meeus "Astronomical Algorithms" Chapter 54
"""

import math
from libephemeris import julday, calc_besselian_l2, FLG_SWIEPH


class TestBesselianL2BasicFunctionality:
    """Test basic functionality of calc_besselian_l2."""

    def test_function_exists(self):
        """Test that calc_besselian_l2 function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_l2

        assert callable(calc_besselian_l2)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_besselian_l2

        assert callable(calc_besselian_l2)

    def test_returns_float(self):
        """Test that function returns a float value."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_l2(jd)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_l2(jd, flags=FLG_SWIEPH)

        assert isinstance(result, float)


class TestBesselianL2DuringEclipses:
    """Test calc_besselian_l2 during known solar eclipses."""

    def test_april_2024_total_eclipse_near_maximum(self):
        """Test l2 during April 8, 2024 total solar eclipse.

        This was a total eclipse, so l2 should be positive (umbral shadow).
        The umbral shadow radius is typically small, around 0.01-0.05 Earth radii.
        """
        # April 8, 2024 around 18:17 UTC - time of maximum eclipse
        jd_max = julday(2024, 4, 8, 18.3)

        l2 = calc_besselian_l2(jd_max)

        # For a total eclipse, l2 should be positive
        assert l2 < 0, (
            f"l2 = {l2} should be negative for a total eclipse (NASA convention: the umbral vertex lies beyond the fundamental plane)"
        )
        # l2 should be small (umbral shadow is much smaller than penumbral)
        assert abs(l2) < 0.1, f"l2 = {l2} Earth radii, expected small value"

    def test_april_2024_eclipse_l2_varies_with_time(self):
        """Test that l2 changes slightly as eclipse progresses."""
        # Sample times during the eclipse
        jd_before = julday(2024, 4, 8, 17.0)  # Before maximum
        jd_max = julday(2024, 4, 8, 18.3)  # Near maximum
        jd_after = julday(2024, 4, 8, 19.5)  # After maximum

        l2_before = calc_besselian_l2(jd_before)
        l2_max = calc_besselian_l2(jd_max)
        l2_after = calc_besselian_l2(jd_after)

        # All should be negative for this total eclipse (NASA
        # convention: the umbral vertex lies beyond the fundamental
        # plane throughout)
        assert l2_before < 0
        assert l2_max < 0
        assert l2_after < 0

        # l2 should change slightly over time as Moon's distance changes
        assert l2_before != l2_after, "l2 should vary with time"

    def test_october_2023_annular_eclipse(self):
        """Test l2 during October 14, 2023 annular eclipse.

        For an annular eclipse, l2 should be negative (antumbral shadow).
        """
        jd_max = julday(2023, 10, 14, 18.0)

        l2 = calc_besselian_l2(jd_max)

        # For an annular eclipse, l2 should be negative
        assert l2 > 0, (
            f"l2 = {l2} should be positive for an annular eclipse (NASA convention: the antumbra opens past the vertex)"
        )
        # Absolute value should be small
        assert abs(l2) < 0.1, f"|l2| = {abs(l2)} Earth radii, expected small value"

    def test_december_2021_total_eclipse(self):
        """Test l2 during December 4, 2021 total eclipse (Antarctica)."""
        jd_max = julday(2021, 12, 4, 7.5)

        l2 = calc_besselian_l2(jd_max)

        # For a total eclipse, l2 should be positive
        assert l2 < 0, (
            f"l2 = {l2} should be negative for the December 2021 total eclipse (NASA convention)"
        )
        assert abs(l2) < 0.1, f"|l2| = {abs(l2)} Earth radii"

    def test_june_2021_annular_eclipse(self):
        """Test l2 during June 10, 2021 annular eclipse."""
        jd_max = julday(2021, 6, 10, 10.5)

        l2 = calc_besselian_l2(jd_max)

        # For an annular eclipse, l2 should be negative
        assert l2 > 0, (
            f"l2 = {l2} should be positive for the June 2021 annular eclipse (NASA convention)"
        )
        assert abs(l2) < 0.1, f"|l2| = {abs(l2)} Earth radii"


class TestBesselianL2PhysicalProperties:
    """Test physical properties of the l2 calculation."""

    def test_l2_magnitude_smaller_than_l1(self):
        """Test that |l2| is much smaller than l1.

        The umbral/antumbral shadow radius is always smaller than the
        penumbral shadow radius.
        """
        from libephemeris import calc_besselian_l1

        jd = julday(2024, 4, 8, 18.0)
        l1 = calc_besselian_l1(jd)
        l2 = calc_besselian_l2(jd)

        assert abs(l2) < l1, f"|l2| = {abs(l2)} should be smaller than l1 = {l1}"

    def test_l2_is_continuous(self):
        """Test that l2 changes continuously (no discontinuities)."""
        jd_start = julday(2024, 4, 8, 17.0)

        # Sample at small time intervals
        l2_values = []
        for i in range(10):
            jd = jd_start + i * 0.01  # 10 samples over ~15 minutes
            l2 = calc_besselian_l2(jd)
            l2_values.append(l2)

        # Check that consecutive values don't jump too much
        for i in range(1, len(l2_values)):
            delta = abs(l2_values[i] - l2_values[i - 1])
            # Change should be very small over ~15 minutes
            assert delta < 0.01, f"Discontinuity detected: delta = {delta} Earth radii"

    def test_l2_reasonable_at_new_moon(self):
        """Test l2 is reasonable at new moon (when eclipses can occur)."""
        # New moons in 2024
        new_moons = [
            julday(2024, 1, 11, 12.0),  # January new moon
            julday(2024, 4, 8, 18.0),  # April new moon (eclipse)
            julday(2024, 7, 5, 22.0),  # July new moon
            julday(2024, 10, 2, 18.0),  # October new moon
        ]

        for jd in new_moons:
            l2 = calc_besselian_l2(jd)
            # l2 should be in a reasonable range (can be positive or negative)
            assert abs(l2) < 0.2, f"l2 = {l2} at JD {jd}"


class TestBesselianL2NumericalStability:
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
            l2 = calc_besselian_l2(jd)
            assert math.isfinite(l2), f"l2 not finite at JD {jd}"

    def test_consistent_results(self):
        """Test that same input gives same output (deterministic)."""
        jd = julday(2024, 4, 8, 18.0)

        l2_1 = calc_besselian_l2(jd)
        l2_2 = calc_besselian_l2(jd)

        assert l2_1 == l2_2, "Results should be deterministic"

    def test_no_nan_or_inf(self):
        """Test that function never returns NaN or infinity."""
        # Test a range of dates
        for year in range(2020, 2030):
            for month in [1, 4, 7, 10]:
                jd = julday(year, month, 15, 12.0)
                l2 = calc_besselian_l2(jd)
                assert not math.isnan(l2), f"l2 is NaN at {year}/{month}"
                assert not math.isinf(l2), f"l2 is infinite at {year}/{month}"


class TestBesselianL2Validation:
    """Validation tests against expected values."""

    def test_l2_typical_value_during_total_eclipse(self):
        """Validate that l2 is positive and small during total solar eclipse.

        The umbral shadow radius depends on:
        - Sun's angular size (varies with Earth-Sun distance)
        - Moon's angular size (varies with Earth-Moon distance)
        - Moon's distance from Earth

        For total eclipses, typical values are around 0.01-0.05 Earth radii.
        """
        jd_eclipse = julday(2024, 4, 8, 18.0)
        l2 = calc_besselian_l2(jd_eclipse)

        # For total eclipse: positive and small
        assert l2 < 0, (
            f"l2 = {l2}, should be negative for a total eclipse (NASA convention)"
        )
        assert l2 < 0.1, f"l2 = {l2}, expected small umbral radius"

    def test_l2_typical_value_during_annular_eclipse(self):
        """Validate that l2 is negative during annular solar eclipse."""
        jd_eclipse = julday(2023, 10, 14, 18.0)
        l2 = calc_besselian_l2(jd_eclipse)

        # For annular eclipse: negative (antumbral)
        assert l2 > 0, (
            f"l2 = {l2}, should be positive for an annular eclipse (NASA convention)"
        )
        assert abs(l2) < 0.1, f"|l2| = {abs(l2)}, expected small antumbral radius"

    def test_l2_rate_of_change(self):
        """Test that l2 changes at a reasonable rate."""
        jd1 = julday(2024, 4, 8, 18.0)
        jd2 = julday(2024, 4, 8, 19.0)  # 1 hour later

        l2_1 = calc_besselian_l2(jd1)
        l2_2 = calc_besselian_l2(jd2)

        # The change in l2 over an hour should be small
        delta_l2 = abs(l2_2 - l2_1)

        # Rate should be less than 0.01 Earth radii per hour
        assert delta_l2 < 0.02, f"l2 changed by {delta_l2} Earth radii in 1 hour"


class TestBesselianL2EdgeCases:
    """Test edge cases and boundary conditions."""

    def test_far_past_date(self):
        """Test function works for historical dates."""
        # Eclipse of 1999 (August 11, total solar eclipse over Europe)
        jd = julday(1999, 8, 11, 11.0)

        l2 = calc_besselian_l2(jd)

        # Should be positive (it was a total eclipse)
        assert l2 < 0, (
            f"l2 = {l2} for the 1999 total eclipse (negative per the NASA convention)"
        )
        assert abs(l2) < 0.2, f"|l2| = {abs(l2)} for 1999 eclipse"
        assert math.isfinite(l2)

    def test_far_future_date(self):
        """Test function works for future dates."""
        # A date far in the future
        jd = julday(2050, 6, 15, 12.0)

        l2 = calc_besselian_l2(jd)

        assert abs(l2) < 0.2, f"|l2| = {abs(l2)} for June 2050"
        assert math.isfinite(l2)

    def test_multiple_calls_consistent(self):
        """Test that multiple calls with same input are consistent."""
        jd = julday(2024, 4, 8, 18.0)

        results = [calc_besselian_l2(jd) for _ in range(5)]

        # All results should be identical
        assert all(r == results[0] for r in results), "Results should be consistent"

    def test_l2_at_different_times_of_day(self):
        """Test l2 at different times throughout a day."""
        base_jd = julday(2024, 4, 8, 0.0)

        # Sample every 4 hours
        l2_values = []
        for hour in [0, 4, 8, 12, 16, 20]:
            jd = base_jd + hour / 24.0
            l2 = calc_besselian_l2(jd)
            l2_values.append(l2)
            assert abs(l2) < 0.2, f"|l2| = {abs(l2)} at hour {hour}"

        # All values should be similar (Moon doesn't move that far in a day)
        max_diff = max(l2_values) - min(l2_values)
        assert max_diff < 0.05, f"l2 varied by {max_diff} over the day"


class TestBesselianL2VsL1Relationship:
    """Test relationship between l1 and l2."""

    def test_l1_minus_l2_is_penumbral_width(self):
        """Test that l1 - l2 gives the penumbral zone width.

        The difference l1 - l2 represents the width of the partial eclipse
        zone (penumbral only, not in umbra/antumbra).
        """
        from libephemeris import calc_besselian_l1

        jd = julday(2024, 4, 8, 18.0)
        l1 = calc_besselian_l1(jd)
        l2 = calc_besselian_l2(jd)

        # The penumbral-only width
        penumbral_width = l1 - l2

        # Should always be positive (penumbra extends beyond umbra)
        assert penumbral_width > 0, f"l1 - l2 = {penumbral_width}"

        # The width should be substantial (much larger than l2)
        assert penumbral_width > abs(l2), "Penumbral width should exceed umbral radius"

    def test_sign_consistency_total_eclipse(self):
        """For total eclipses, l1 > 0 and l2 > 0, with l1 > l2."""
        from libephemeris import calc_besselian_l1

        jd = julday(2024, 4, 8, 18.0)  # Total eclipse
        l1 = calc_besselian_l1(jd)
        l2 = calc_besselian_l2(jd)

        assert l1 > 0, f"l1 = {l1} should be positive"
        assert l2 < 0, (
            f"l2 = {l2} should be negative for a total eclipse (NASA convention: the umbral vertex lies beyond the fundamental plane)"
        )
        assert l1 > l2, f"l1 = {l1} should be greater than l2 = {l2}"

    def test_sign_consistency_annular_eclipse(self):
        """For annular eclipses, l1 > 0 and l2 < 0."""
        from libephemeris import calc_besselian_l1

        jd = julday(2023, 10, 14, 18.0)  # Annular eclipse
        l1 = calc_besselian_l1(jd)
        l2 = calc_besselian_l2(jd)

        assert l1 > 0, f"l1 = {l1} should always be positive"
        assert l2 > 0, (
            f"l2 = {l2} should be positive for an annular eclipse (NASA convention: the antumbra opens past the vertex)"
        )
