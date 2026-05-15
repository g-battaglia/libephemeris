"""
Comprehensive tests for utc_to_jd() function.

Tests cover:
- Basic UTC to JD conversion
- Leap second handling
- Gregorian vs Julian calendar
- Comparison with julday() to verify differences
- Edge cases (year boundaries, leap years)
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import GREG_CAL, JUL_CAL


class TestUtcToJdBasic:
    """Test basic UTC to JD conversion."""

    @pytest.mark.unit
    def test_j2000_epoch(self):
        """J2000.0 = 2000-01-01 12:00:00 UTC approximately."""
        jd_tt, jd_ut = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        # J2000.0 is defined at JD(TT) = 2451545.0
        # At this time, Delta T was about 63.8 seconds
        assert jd_tt == pytest.approx(2451545.0, abs=0.001)
        # UT1 should be slightly less than TT (by Delta T)
        assert jd_ut < jd_tt
        # Difference should be approximately Delta T in days
        delta_t_days = jd_tt - jd_ut
        delta_t_seconds = delta_t_days * 86400
        assert 60 < delta_t_seconds < 70  # Delta T was ~64 seconds in 2000

    @pytest.mark.unit
    def test_returns_tuple(self):
        """Function should return a tuple of (jd_tt, jd_ut)."""
        result = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert isinstance(result[0], float)
        assert isinstance(result[1], float)

    @pytest.mark.unit
    def test_midnight(self):
        """Test midnight UTC conversion."""
        jd_tt, jd_ut = ephem.utc_to_jd(2000, 1, 1, 0, 0, 0.0)
        # JD starts at noon, so midnight is .5 fractional day
        assert jd_ut == pytest.approx(2451544.5, abs=0.001)

    @pytest.mark.unit
    def test_noon(self):
        """Test noon UTC conversion."""
        jd_tt, jd_ut = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        # Noon should be close to .0 fractional day
        assert jd_ut == pytest.approx(2451545.0, abs=0.001)

    @pytest.mark.unit
    def test_seconds_precision(self):
        """Test that seconds are properly converted."""
        jd_tt_0, jd_ut_0 = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        jd_tt_30, jd_ut_30 = ephem.utc_to_jd(2000, 1, 1, 12, 0, 30.0)
        # 30 seconds = 30/86400 days
        expected_diff = 30.0 / 86400.0
        assert jd_ut_30 - jd_ut_0 == pytest.approx(expected_diff, rel=1e-6)
        assert jd_tt_30 - jd_tt_0 == pytest.approx(expected_diff, rel=1e-6)


class TestUtcToJdLeapSeconds:
    """Test leap second handling."""

    @pytest.mark.unit
    def test_around_leap_second_2016(self):
        """Test around the 2016-12-31 23:59:60 leap second."""
        # Just before leap second
        jd_tt_before, jd_ut_before = ephem.utc_to_jd(2016, 12, 31, 23, 59, 59.0)
        # Just after leap second (in 2017)
        jd_tt_after, jd_ut_after = ephem.utc_to_jd(2017, 1, 1, 0, 0, 0.0)
        # Should be approximately 2 seconds apart (59->60->0)
        diff_seconds = (jd_ut_after - jd_ut_before) * 86400
        # With leap second, 23:59:59 to 00:00:00 includes the extra second
        assert 1.5 < diff_seconds < 2.5

    @pytest.mark.unit
    def test_tt_ut_difference_increases_after_leap_second(self):
        """TT-UT1 difference should reflect cumulative leap seconds."""
        # Before any leap seconds (1970)
        _, jd_ut_1970 = ephem.utc_to_jd(1970, 1, 1, 0, 0, 0.0)
        # After many leap seconds (2020)
        jd_tt_2020, jd_ut_2020 = ephem.utc_to_jd(2020, 1, 1, 0, 0, 0.0)
        delta_t_2020 = (jd_tt_2020 - jd_ut_2020) * 86400  # seconds
        # Delta T in 2020 is about 69 seconds
        assert 65 < delta_t_2020 < 75

    @pytest.mark.unit
    def test_leap_second_day_2015(self):
        """Test 2015-06-30 leap second day."""
        jd_tt, jd_ut = ephem.utc_to_jd(2015, 6, 30, 23, 59, 59.0)
        # Should produce valid JD values
        assert jd_tt > 2457000  # Around JD for 2015
        assert jd_ut < jd_tt


class TestUtcToJdVsJulday:
    """Compare utc_to_jd() with julday() to verify they differ appropriately."""

    @pytest.mark.unit
    def test_ut_close_to_julday(self):
        """jd_ut from utc_to_jd should be very close to julday() result."""
        year, month, day = 2000, 1, 1
        hour = 12
        minute = 30
        second = 45.5
        decimal_hour = hour + minute / 60.0 + second / 3600.0

        jd_tt, jd_ut = ephem.utc_to_jd(year, month, day, hour, minute, second)
        jd_julday = ephem.julday(year, month, day, decimal_hour)

        # julday assumes UT1 input, utc_to_jd returns proper UT1
        # For modern dates, UTC and UT1 differ by less than 0.9 seconds
        diff_seconds = abs(jd_ut - jd_julday) * 86400
        # The difference should be < 1 second (UTC-UT1 is always < 0.9s)
        assert diff_seconds < 1.0

    @pytest.mark.unit
    def test_tt_differs_from_julday(self):
        """jd_tt from utc_to_jd should differ from julday() by Delta T."""
        year, month, day = 2000, 1, 1
        hour = 12
        minute = 0
        second = 0.0
        decimal_hour = hour + minute / 60.0 + second / 3600.0

        jd_tt, jd_ut = ephem.utc_to_jd(year, month, day, hour, minute, second)
        jd_julday = ephem.julday(year, month, day, decimal_hour)

        # TT differs from UT1 by Delta T (~64 seconds in 2000)
        diff_seconds = (jd_tt - jd_julday) * 86400
        assert 60 < diff_seconds < 70


class TestUtcToJdCalendars:
    """Test Gregorian vs Julian calendar handling."""

    @pytest.mark.unit
    def test_gregorian_default(self):
        """Gregorian calendar should be default."""
        jd_tt_default, jd_ut_default = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        jd_tt_greg, jd_ut_greg = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0, GREG_CAL)
        assert jd_tt_default == jd_tt_greg
        assert jd_ut_default == jd_ut_greg

    @pytest.mark.unit
    def test_julian_calendar_differs(self):
        """Julian calendar should give different result for modern dates."""
        jd_tt_greg, jd_ut_greg = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0, GREG_CAL)
        jd_tt_jul, jd_ut_jul = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0, JUL_CAL)
        # Julian calendar is ~13 days behind Gregorian in 2000
        diff_days = abs(jd_ut_greg - jd_ut_jul)
        assert diff_days > 10

    @pytest.mark.unit
    def test_julian_calendar_before_reform(self):
        """Test Julian calendar for dates before Gregorian reform."""
        # Oct 5, 1582 Julian = Oct 15, 1582 Gregorian
        jd_tt, jd_ut = ephem.utc_to_jd(1582, 10, 5, 12, 0, 0.0, JUL_CAL)
        # Should produce valid JD
        assert jd_ut > 0


class TestUtcToJdEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.unit
    def test_year_boundary(self):
        """Test New Year transition."""
        jd_tt_dec31, jd_ut_dec31 = ephem.utc_to_jd(1999, 12, 31, 23, 59, 59.0)
        jd_tt_jan1, jd_ut_jan1 = ephem.utc_to_jd(2000, 1, 1, 0, 0, 0.0)
        diff_seconds = (jd_ut_jan1 - jd_ut_dec31) * 86400
        # Should be 1 second (plus any leap second)
        assert 0.5 < diff_seconds < 2.5

    @pytest.mark.unit
    def test_leap_year_feb29(self):
        """Test Feb 29 in leap year."""
        jd_tt, jd_ut = ephem.utc_to_jd(2000, 2, 29, 12, 0, 0.0)
        # 2000 is a leap year, so Feb 29 should work
        assert jd_ut > 0
        # Should be 1 day after Feb 28
        jd_tt_28, jd_ut_28 = ephem.utc_to_jd(2000, 2, 28, 12, 0, 0.0)
        assert jd_ut - jd_ut_28 == pytest.approx(1.0, abs=0.001)

    @pytest.mark.unit
    def test_fractional_seconds(self):
        """Test fractional seconds are handled."""
        jd_tt_0, jd_ut_0 = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        jd_tt_half, jd_ut_half = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.5)
        diff_seconds = (jd_ut_half - jd_ut_0) * 86400
        assert diff_seconds == pytest.approx(0.5, rel=1e-3)

    @pytest.mark.unit
    def test_historical_date(self):
        """Test historical date (before UTC was defined)."""
        # 1800 is before UTC (1972), function should still work
        jd_tt, jd_ut = ephem.utc_to_jd(1800, 6, 15, 12, 0, 0.0)
        assert jd_ut > 0
        # Delta T was very different historically
        delta_t_seconds = (jd_tt - jd_ut) * 86400
        # In 1800, Delta T was about 13-14 seconds
        assert 10 < delta_t_seconds < 20

    @pytest.mark.unit
    def test_future_date(self):
        """Test future date."""
        jd_tt, jd_ut = ephem.utc_to_jd(2100, 1, 1, 12, 0, 0.0)
        assert jd_ut > 2451545.0  # After J2000
        assert jd_tt > jd_ut  # TT is always ahead of UT1


class TestUtcToJdConsistency:
    """Test internal consistency."""

    @pytest.mark.unit
    def test_tt_always_greater_than_ut(self):
        """TT (jd_tt) should always be greater than UT1 (jd_ut) for modern dates."""
        # Test several dates
        dates = [
            (1980, 1, 1),
            (2000, 6, 15),
            (2020, 12, 31),
        ]
        for year, month, day in dates:
            jd_tt, jd_ut = ephem.utc_to_jd(year, month, day, 12, 0, 0.0)
            assert jd_tt > jd_ut, f"TT should be > UT1 for {year}-{month}-{day}"

    @pytest.mark.unit
    def test_consecutive_seconds(self):
        """Consecutive seconds should differ by approximately 1/86400 days."""
        jd_tt_0, jd_ut_0 = ephem.utc_to_jd(2020, 6, 15, 12, 30, 0.0)
        jd_tt_1, jd_ut_1 = ephem.utc_to_jd(2020, 6, 15, 12, 30, 1.0)
        expected = 1.0 / 86400.0
        # UT1 has inherent microsecond-level variations due to Earth rotation
        assert jd_ut_1 - jd_ut_0 == pytest.approx(expected, rel=1e-4)

    @pytest.mark.unit
    def test_consecutive_minutes(self):
        """Consecutive minutes should differ by 60/86400 days."""
        jd_tt_0, jd_ut_0 = ephem.utc_to_jd(2020, 6, 15, 12, 30, 0.0)
        jd_tt_1, jd_ut_1 = ephem.utc_to_jd(2020, 6, 15, 12, 31, 0.0)
        expected = 60.0 / 86400.0
        assert jd_ut_1 - jd_ut_0 == pytest.approx(expected, rel=1e-6)

    @pytest.mark.unit
    def test_consecutive_hours(self):
        """Consecutive hours should differ by 1/24 days."""
        jd_tt_0, jd_ut_0 = ephem.utc_to_jd(2020, 6, 15, 12, 0, 0.0)
        jd_tt_1, jd_ut_1 = ephem.utc_to_jd(2020, 6, 15, 13, 0, 0.0)
        expected = 1.0 / 24.0
        assert jd_ut_1 - jd_ut_0 == pytest.approx(expected, rel=1e-6)
