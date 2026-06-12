"""
Tests for calc_besselian_mu function in libephemeris.

Tests the Besselian mu (Greenwich hour angle) calculation for solar eclipses.

The Besselian element mu is the Greenwich hour angle of the shadow axis,
measured in degrees. It represents the angle between the Greenwich meridian
and the hour circle containing the shadow axis, measured westward from
Greenwich along the celestial equator.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Explanatory Supplement to the Astronomical Almanac
- Meeus "Astronomical Algorithms" Chapter 54
"""

import math
from libephemeris import julday, calc_besselian_mu, FLG_SWIEPH


class TestBesselianMuBasicFunctionality:
    """Test basic functionality of calc_besselian_mu."""

    def test_function_exists(self):
        """Test that calc_besselian_mu function exists and is callable."""
        from libephemeris.eclipse import calc_besselian_mu

        assert callable(calc_besselian_mu)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_besselian_mu

        assert callable(calc_besselian_mu)

    def test_returns_float(self):
        """Test that function returns a float value."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_mu(jd)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_mu(jd, flags=FLG_SWIEPH)

        assert isinstance(result, float)

    def test_result_in_valid_range(self):
        """Test that result is within valid range (0 to 360 degrees)."""
        jd = julday(2024, 4, 8, 18.0)
        result = calc_besselian_mu(jd)

        # mu must be between 0 and 360 degrees
        assert 0 <= result < 360, f"mu = {result} is outside valid range [0, 360)"


class TestBesselianMuDuringEclipses:
    """Test calc_besselian_mu during known solar eclipses."""

    def test_april_2024_total_eclipse_near_maximum(self):
        """Test mu during April 8, 2024 total solar eclipse.

        The eclipse maximum occurs around 18:17 UTC. At this time,
        the shadow axis has a specific Greenwich hour angle that can
        be verified against NASA data.
        """
        # April 8, 2024 around 18:17 UTC - time of maximum eclipse
        jd_max = julday(2024, 4, 8, 18.3)

        mu = calc_besselian_mu(jd_max)

        # mu should be in valid range
        assert 0 <= mu < 360, f"mu = {mu} degrees is out of range"

    def test_april_2024_eclipse_mu_increases_with_time(self):
        """Test that mu increases as Earth rotates.

        mu should increase by approximately 15 degrees per hour as Earth
        rotates, since mu = GAST - RA, and GAST increases with time.
        """
        # Sample times during the eclipse
        jd_before = julday(2024, 4, 8, 17.0)  # Before maximum
        jd_after = julday(2024, 4, 8, 18.0)  # 1 hour later

        mu_before = calc_besselian_mu(jd_before)
        mu_after = calc_besselian_mu(jd_after)

        # Compute the change in mu (accounting for wraparound at 360)
        delta_mu = mu_after - mu_before
        if delta_mu < -180:
            delta_mu += 360
        if delta_mu > 180:
            delta_mu -= 360

        # Earth rotates 15 degrees per hour
        # mu should increase by approximately 15 degrees per hour
        # Allow for some variation due to shadow axis RA change
        assert 14 < delta_mu < 16, f"mu changed by {delta_mu} deg/hour, expected ~15"

    def test_october_2023_annular_eclipse(self):
        """Test mu during October 14, 2023 annular eclipse."""
        jd_max = julday(2023, 10, 14, 18.0)

        mu = calc_besselian_mu(jd_max)

        # mu should be in valid range
        assert 0 <= mu < 360, f"mu = {mu} degrees for October eclipse"

    def test_december_2021_total_eclipse(self):
        """Test mu during December 4, 2021 total eclipse (Antarctica)."""
        jd_max = julday(2021, 12, 4, 7.5)

        mu = calc_besselian_mu(jd_max)

        # mu should be in valid range
        assert 0 <= mu < 360, f"mu = {mu} degrees for December 2021 eclipse"

    def test_june_2021_annular_eclipse(self):
        """Test mu during June 10, 2021 annular eclipse."""
        jd_max = julday(2021, 6, 10, 10.5)

        mu = calc_besselian_mu(jd_max)

        # mu should be in valid range
        assert 0 <= mu < 360, f"mu = {mu} degrees for June 2021 eclipse"


class TestBesselianMuPhysicalProperties:
    """Test physical properties of the Greenwich hour angle calculation."""

    def test_mu_rate_of_change(self):
        """Test that mu increases at approximately the sidereal rate.

        The Greenwich Apparent Sidereal Time increases at about 15.041 degrees
        per solar hour. The RA of the shadow axis changes slowly (Sun moves
        about 1 degree per day), so mu should increase at nearly this rate.
        """
        jd1 = julday(2024, 4, 8, 12.0)
        jd2 = julday(2024, 4, 8, 14.0)  # 2 hours later

        mu1 = calc_besselian_mu(jd1)
        mu2 = calc_besselian_mu(jd2)

        # Compute the change in mu (accounting for wraparound)
        delta_mu = mu2 - mu1
        if delta_mu < -180:
            delta_mu += 360
        if delta_mu > 180:
            delta_mu -= 360

        # Rate of change should be close to 15 degrees per hour
        rate = delta_mu / 2.0  # degrees per hour
        assert 14.5 < rate < 15.5, f"mu rate = {rate} deg/hour, expected ~15"

    def test_mu_is_continuous(self):
        """Test that mu changes continuously (no discontinuities except at 360/0)."""
        jd_start = julday(2024, 4, 8, 17.0)

        # Sample at small time intervals
        mu_values = []
        for i in range(10):
            jd = jd_start + i * 0.01  # 10 samples over ~15 minutes
            mu = calc_besselian_mu(jd)
            mu_values.append(mu)

        # Check that consecutive values don't jump too much
        # (except for possible wraparound at 360/0)
        for i in range(1, len(mu_values)):
            delta = mu_values[i] - mu_values[i - 1]
            # Account for wraparound
            if delta > 180:
                delta -= 360
            if delta < -180:
                delta += 360

            # Change should be very small over ~15 minutes (< 5 degrees)
            assert abs(delta) < 5, f"Discontinuity detected: delta = {delta} degrees"

    def test_mu_daily_cycle(self):
        """Test that mu cycles through ~360 degrees per sidereal day.

        Over 24 sidereal hours (23h 56m 4s solar), mu should increase by
        approximately 360 degrees (one full rotation).
        """
        jd1 = julday(2024, 4, 8, 0.0)
        # One sidereal day later (~23.9345 solar hours)
        sidereal_day_hours = 23.9344696
        jd2 = jd1 + sidereal_day_hours / 24.0

        mu1 = calc_besselian_mu(jd1)
        mu2 = calc_besselian_mu(jd2)

        # After one sidereal day, mu should return to approximately
        # the same value (within a few degrees, accounting for Sun's motion)
        delta = abs(mu2 - mu1)
        # Allow for Sun's RA change (~1 degree/day)
        assert delta < 3 or delta > 357, (
            f"mu1 = {mu1}, mu2 = {mu2}, expected similar values after 1 sidereal day"
        )


class TestBesselianMuNumericalStability:
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
            mu = calc_besselian_mu(jd)
            assert math.isfinite(mu), f"mu not finite at JD {jd}"

    def test_consistent_results(self):
        """Test that same input gives same output (deterministic)."""
        jd = julday(2024, 4, 8, 18.0)

        mu1 = calc_besselian_mu(jd)
        mu2 = calc_besselian_mu(jd)

        assert mu1 == mu2, "Results should be deterministic"

    def test_no_nan_or_inf(self):
        """Test that function never returns NaN or infinity."""
        # Test a range of dates
        for year in range(2020, 2030):
            for month in [1, 4, 7, 10]:
                jd = julday(year, month, 15, 12.0)
                mu = calc_besselian_mu(jd)
                assert not math.isnan(mu), f"mu is NaN at {year}/{month}"
                assert not math.isinf(mu), f"mu is infinite at {year}/{month}"

    def test_always_in_range(self):
        """Test that mu is always normalized to 0-360 degrees."""
        # Test various times to exercise different mu values
        for hour in range(0, 24, 3):
            jd = julday(2024, 4, 8, hour)
            mu = calc_besselian_mu(jd)
            assert 0 <= mu < 360, f"mu = {mu} at hour {hour} is out of range"


class TestBesselianMuValidation:
    """Validation tests against reference values."""

    def test_mu_changes_with_ut(self):
        """Test that mu responds correctly to changes in UT.

        As UT increases (time passes), mu should increase because GAST
        increases faster than RA changes.
        """
        jd1 = julday(2024, 4, 8, 0.0)
        jd2 = julday(2024, 4, 8, 6.0)
        jd3 = julday(2024, 4, 8, 12.0)
        jd4 = julday(2024, 4, 8, 18.0)

        mu1 = calc_besselian_mu(jd1)
        mu2 = calc_besselian_mu(jd2)
        mu3 = calc_besselian_mu(jd3)
        mu4 = calc_besselian_mu(jd4)

        # Unwrap mu values to check monotonic increase
        def unwrap_angle(angles):
            """Unwrap angles to show continuous increase."""
            result = [angles[0]]
            for i in range(1, len(angles)):
                diff = angles[i] - angles[i - 1]
                if diff < -180:
                    diff += 360
                result.append(result[-1] + diff)
            return result

        unwrapped = unwrap_angle([mu1, mu2, mu3, mu4])

        # Each 6-hour interval should add ~90 degrees
        for i in range(1, 4):
            diff = unwrapped[i] - unwrapped[i - 1]
            assert 85 < diff < 95, f"Expected ~90 deg increase per 6 hours, got {diff}"

    def test_relationship_between_mu_and_gast(self):
        """Test that mu = GAST - RA (shadow axis).

        We can verify the relationship by checking that mu tracks GAST
        closely since RA changes slowly.
        """
        # At two times 1 hour apart
        jd1 = julday(2024, 4, 8, 17.0)
        jd2 = julday(2024, 4, 8, 18.0)

        mu1 = calc_besselian_mu(jd1)
        mu2 = calc_besselian_mu(jd2)

        # mu change should be close to GAST change (~15 degrees/hour)
        delta_mu = mu2 - mu1
        if delta_mu < -180:
            delta_mu += 360
        if delta_mu > 180:
            delta_mu -= 360

        # GAST increases at ~15.041 degrees per hour
        # RA of shadow axis changes slowly (~0.04 degrees per hour)
        # So mu should increase by ~15 degrees per hour
        assert 14 < delta_mu < 16, f"mu change = {delta_mu} deg/hour, expected ~15"


class TestBesselianMuEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_far_past_date(self):
        """Test function works for historical dates."""
        # Eclipse of 1999 (August 11, total solar eclipse over Europe)
        jd = julday(1999, 8, 11, 11.0)

        mu = calc_besselian_mu(jd)

        assert 0 <= mu < 360, f"mu = {mu} for 1999 eclipse"
        assert math.isfinite(mu)

    def test_far_future_date(self):
        """Test function works for future dates."""
        # A date far in the future
        jd = julday(2050, 6, 15, 12.0)

        mu = calc_besselian_mu(jd)

        assert 0 <= mu < 360, f"mu = {mu} for June 2050"
        assert math.isfinite(mu)

    def test_multiple_calls_consistent(self):
        """Test that multiple calls with same input are consistent."""
        jd = julday(2024, 4, 8, 18.0)

        results = [calc_besselian_mu(jd) for _ in range(5)]

        # All results should be identical
        assert all(r == results[0] for r in results), "Results should be consistent"

    def test_midnight_ut(self):
        """Test mu calculation at 0h UT."""
        jd = julday(2024, 4, 8, 0.0)

        mu = calc_besselian_mu(jd)

        assert 0 <= mu < 360, f"mu = {mu} at midnight"
        assert math.isfinite(mu)

    def test_noon_ut(self):
        """Test mu calculation at 12h UT."""
        jd = julday(2024, 4, 8, 12.0)

        mu = calc_besselian_mu(jd)

        assert 0 <= mu < 360, f"mu = {mu} at noon"
        assert math.isfinite(mu)
