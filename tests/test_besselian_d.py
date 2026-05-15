"""
Tests for calc_besselian_d function in libephemeris.

Tests the Besselian d (declination) calculation for solar eclipses.

The Besselian element d is the declination of the Moon's shadow axis,
measured in degrees. This is the angle between the shadow axis direction
and the celestial equatorial plane.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Explanatory Supplement to the Astronomical Almanac
- Meeus "Astronomical Algorithms" Chapter 54
"""

import pytest
import math
from libephemeris import julday, calc_besselian_d, FLG_SWIEPH


class TestBesselianDBasicFunctionality:
    """Test basic functionality of calc_besselian_d."""

    def test_function_exists(self):
        """Test that calc_besselian_d function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_d

        assert callable(calc_besselian_d)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_besselian_d

        assert callable(calc_besselian_d)

    def test_returns_float(self):
        """Test that function returns a float value."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_d(jd)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_d(jd, flags=FLG_SWIEPH)

        assert isinstance(result, float)

    def test_result_in_valid_range(self):
        """Test that result is within valid declination range (-90 to +90)."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_d(jd)

        # Declination must be between -90 and +90 degrees
        assert -90 <= result <= 90, f"d = {result} is outside valid range [-90, 90]"


class TestBesselianDDuringEclipses:
    """Test calc_besselian_d during known solar eclipses."""

    def test_april_2024_total_eclipse_near_maximum(self):
        """Test d during April 8, 2024 total solar eclipse.

        NASA gives d0 (declination at T0) approximately 7.6 degrees for
        this eclipse. The shadow axis points north of the equator
        because the eclipse occurs after the vernal equinox when the
        Sun's declination is positive (around +7 degrees).
        """
        # April 8, 2024 around 18:17 UTC - time of maximum eclipse
        jd_max = julday(2024, 4, 8, 18.3)

        d = calc_besselian_d(jd_max)

        # Sun's declination on April 8 is approximately +7.5 degrees
        # The shadow axis declination should be similar (opposite sign since
        # it points from Sun toward Moon, roughly anti-solar)
        # Actually, for solar eclipses, d tracks the solar declination closely
        # since Moon and Sun are nearly aligned
        assert 5.0 < d < 10.0, f"d = {d} degrees, expected ~7.6 degrees"

    def test_april_2024_eclipse_d_varies_with_time(self):
        """Test that d changes slightly as eclipse progresses."""
        # Sample times during the eclipse
        jd_before = julday(2024, 4, 8, 17.0)  # Before maximum
        jd_max = julday(2024, 4, 8, 18.3)  # Near maximum
        jd_after = julday(2024, 4, 8, 19.5)  # After maximum

        d_before = calc_besselian_d(jd_before)
        d_max = calc_besselian_d(jd_max)
        d_after = calc_besselian_d(jd_after)

        # All should be in the same general range (positive, ~7 degrees)
        assert all(5.0 < d < 10.0 for d in [d_before, d_max, d_after])

        # d should change slightly over time (Sun's declination changes ~1 deg/day)
        # Over 2.5 hours, change should be small but detectable
        assert d_before != d_after, "d should vary with time"

    def test_october_2023_annular_eclipse(self):
        """Test d during October 14, 2023 annular eclipse.

        In mid-October, the Sun's declination is around -8 to -9 degrees.
        """
        jd_max = julday(2023, 10, 14, 18.0)

        d = calc_besselian_d(jd_max)

        # Sun's declination in mid-October is negative (around -9 degrees)
        assert -12.0 < d < -5.0, f"d = {d} degrees for October eclipse"

    def test_december_2021_total_eclipse(self):
        """Test d during December 4, 2021 total eclipse (Antarctica).

        In early December, the Sun's declination is around -22 degrees
        (approaching winter solstice).
        """
        jd_max = julday(2021, 12, 4, 7.5)

        d = calc_besselian_d(jd_max)

        # Sun's declination in early December is around -22 degrees
        assert -25.0 < d < -18.0, f"d = {d} degrees for December 2021 eclipse"

    def test_june_2021_annular_eclipse(self):
        """Test d during June 10, 2021 annular eclipse.

        In early June, the Sun's declination is around +23 degrees
        (approaching summer solstice).
        """
        jd_max = julday(2021, 6, 10, 10.5)

        d = calc_besselian_d(jd_max)

        # Sun's declination in early June is around +23 degrees
        assert 20.0 < d < 25.0, f"d = {d} degrees for June 2021 eclipse"


class TestBesselianDPhysicalProperties:
    """Test physical properties of the declination calculation."""

    def test_d_follows_solar_declination_pattern(self):
        """Test that d roughly follows the Sun's declination through the year.

        Since eclipses occur when Moon and Sun are aligned (new Moon),
        the shadow axis declination should closely match the solar declination.
        Note: Away from eclipses, the Moon can be up to ~5 degrees from the
        ecliptic, which affects the shadow axis direction.
        """
        # Test at different times of year (not necessarily during eclipses)
        # Use wider tolerances since Moon's orbital inclination affects d
        dates = [
            (julday(2024, 3, 20, 12.0), -5.0, 5.0),  # Vernal equinox: d ~ 0
            (julday(2024, 6, 21, 12.0), 18.0, 28.0),  # Summer solstice: d ~ +23
            (julday(2024, 9, 22, 12.0), -5.0, 5.0),  # Autumnal equinox: d ~ 0
            (julday(2024, 12, 21, 12.0), -28.0, -18.0),  # Winter solstice: d ~ -23
        ]

        for jd, d_min, d_max in dates:
            d = calc_besselian_d(jd)
            assert d_min < d < d_max, (
                f"d = {d} at JD {jd}, expected {d_min} < d < {d_max}"
            )

    def test_d_is_continuous(self):
        """Test that d changes continuously (no discontinuities)."""
        jd_start = julday(2024, 4, 8, 17.0)

        # Sample at small time intervals
        d_values = []
        for i in range(10):
            jd = jd_start + i * 0.01  # 10 samples over ~15 minutes
            d = calc_besselian_d(jd)
            d_values.append(d)

        # Check that consecutive values don't jump too much
        for i in range(1, len(d_values)):
            delta = abs(d_values[i] - d_values[i - 1])
            # Change should be very small over ~15 minutes (< 0.01 degrees)
            assert delta < 0.1, f"Discontinuity detected: delta = {delta} degrees"

    def test_d_bounded_by_ecliptic_obliquity(self):
        """Test that d is bounded by the Earth's axial tilt.

        The Sun's declination is bounded by the obliquity of the ecliptic
        (~23.4 degrees). Since d tracks solar declination during eclipses,
        it should also be bounded by this value (with small variations
        due to Moon's orbital inclination of ~5 degrees).
        """
        # Test various dates throughout the year
        test_dates = [
            julday(2024, 1, 15, 12.0),
            julday(2024, 4, 15, 12.0),
            julday(2024, 7, 15, 12.0),
            julday(2024, 10, 15, 12.0),
        ]

        for jd in test_dates:
            d = calc_besselian_d(jd)
            # Max obliquity ~23.4 deg + Moon inclination ~5 deg = ~28.4 deg max
            assert abs(d) < 30.0, f"d = {d} exceeds expected bounds"


class TestBesselianDNumericalStability:
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
            d = calc_besselian_d(jd)
            assert math.isfinite(d), f"d not finite at JD {jd}"

    def test_consistent_results(self):
        """Test that same input gives same output (deterministic)."""
        jd = julday(2024, 4, 8, 18.0)

        d1 = calc_besselian_d(jd)
        d2 = calc_besselian_d(jd)

        assert d1 == d2, "Results should be deterministic"

    def test_no_nan_or_inf(self):
        """Test that function never returns NaN or infinity."""
        # Test a range of dates
        for year in range(2020, 2030):
            for month in [1, 4, 7, 10]:
                jd = julday(year, month, 15, 12.0)
                d = calc_besselian_d(jd)
                assert not math.isnan(d), f"d is NaN at {year}/{month}"
                assert not math.isinf(d), f"d is infinite at {year}/{month}"


class TestBesselianDValidation:
    """Validation tests against reference values."""

    def test_april_2024_eclipse_nasa_reference(self):
        """Validate against NASA eclipse data for April 8, 2024.

        From NASA eclipse data:
        - T0 (reference time): 2024 Apr 8 at 18:00 TDT
        - d0 (declination at T0): approximately 7.58 degrees
        """
        # T0 for this eclipse (approximately 18:00 TDT, close to UT)
        jd_t0 = julday(2024, 4, 8, 18.0)

        d = calc_besselian_d(jd_t0)

        # NASA reference d0 ~ 7.58 degrees
        # Allow tolerance for differences in ephemeris precision
        assert abs(d - 7.58) < 0.5, f"d = {d}, NASA reference ~ 7.58 degrees"

    def test_d_rate_of_change(self):
        """Test that d changes at a reasonable rate."""
        jd1 = julday(2024, 4, 8, 18.0)
        jd2 = julday(2024, 4, 8, 19.0)  # 1 hour later

        d1 = calc_besselian_d(jd1)
        d2 = calc_besselian_d(jd2)

        # The Sun's declination changes about 1 degree per day maximum
        # (around equinoxes). Per hour, this is ~0.04 degrees maximum.
        # The Moon's orbital motion adds a small additional component.
        delta_d = abs(d2 - d1)

        # Rate should be less than 0.1 degrees per hour
        assert delta_d < 0.1, f"d changed by {delta_d} degrees in 1 hour"

    def test_relationship_with_solar_position(self):
        """Test that d correlates with the Sun's declination.

        During a solar eclipse, the Moon and Sun are nearly aligned,
        so the shadow axis direction should match the anti-solar direction,
        and d should be close to the Sun's declination.
        """
        # Use the April 2024 eclipse where we know solar declination
        jd = julday(2024, 4, 8, 18.0)

        d = calc_besselian_d(jd)

        # Sun's geocentric apparent declination on April 8, 2024 is about +7.6 deg
        # The shadow axis points from Sun toward Moon, so its declination
        # should be very close to the solar declination during an eclipse
        sun_dec_approx = 7.6  # degrees

        # Should match within a few degrees (Moon's orbital inclination ~5 deg)
        assert abs(d - sun_dec_approx) < 5.0, (
            f"d = {d}, expected close to solar dec {sun_dec_approx}"
        )


class TestBesselianDEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_far_past_date(self):
        """Test function works for historical dates."""
        # Eclipse of 1999 (August 11, total solar eclipse over Europe)
        jd = julday(1999, 8, 11, 11.0)

        d = calc_besselian_d(jd)

        # In August, Sun's declination is positive, around +15 degrees
        assert 10.0 < d < 20.0, f"d = {d} for 1999 eclipse"
        assert math.isfinite(d)

    def test_far_future_date(self):
        """Test function works for future dates."""
        # A date far in the future
        jd = julday(2050, 6, 15, 12.0)

        d = calc_besselian_d(jd)

        # In June, Sun's declination is around +23 degrees
        assert 20.0 < d < 25.0, f"d = {d} for June 2050"
        assert math.isfinite(d)

    def test_multiple_calls_consistent(self):
        """Test that multiple calls with same input are consistent."""
        jd = julday(2024, 4, 8, 18.0)

        results = [calc_besselian_d(jd) for _ in range(5)]

        # All results should be identical
        assert all(r == results[0] for r in results), "Results should be consistent"
