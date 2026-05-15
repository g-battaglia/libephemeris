"""
Comprehensive tests for jdut1_to_utc() function.

Tests cover:
- Basic JD(UT1) to UTC conversion
- Roundtrip consistency with utc_to_jd()
- Difference between UT1 and UTC (DUT1)
- Gregorian vs Julian calendar handling
- Edge cases (historical dates, future dates, year boundaries)
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import GREG_CAL, JUL_CAL


class TestJdut1ToUtcBasic:
    """Test basic JD(UT1) to UTC conversion."""

    @pytest.mark.unit
    def test_j2000_epoch(self):
        """Test conversion near J2000.0 epoch."""
        # Get JD(UT1) for a known date
        jd_tt, jd_ut1 = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        assert year == 2000
        assert month == 1
        assert day == 1
        assert hour == 12
        assert minute == 0
        # UT1 differs from UTC by less than 0.9 seconds
        assert abs(second) < 1.0

    @pytest.mark.unit
    def test_returns_tuple(self):
        """Function should return a tuple of 6 elements."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 6, 15, 12, 0, 0.0)
        result = ephem.jdut1_to_utc(jd_ut1)
        assert isinstance(result, tuple)
        assert len(result) == 6
        # year, month, day, hour, minute should be int
        assert isinstance(result[0], int)
        assert isinstance(result[1], int)
        assert isinstance(result[2], int)
        assert isinstance(result[3], int)
        assert isinstance(result[4], int)
        # second should be float
        assert isinstance(result[5], float)

    @pytest.mark.unit
    def test_noon_ut1(self):
        """Test JD(UT1) at noon."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2023, 2, 24, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        assert year == 2023
        assert month == 2
        assert day == 24
        assert hour == 12
        assert minute == 0
        # UTC should be very close to the original (within DUT1 < 0.9s)
        assert abs(second) < 1.0

    @pytest.mark.unit
    def test_midnight_ut1(self):
        """Test JD(UT1) at midnight."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2023, 2, 25, 0, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        assert year == 2023
        assert month == 2
        assert day == 25
        assert hour == 0
        assert minute == 0
        # UTC should be very close (within DUT1)
        assert abs(second) < 1.0


class TestJdut1ToUtcRoundtrip:
    """Test roundtrip consistency with utc_to_jd()."""

    @pytest.mark.unit
    def test_roundtrip_modern_date(self):
        """Converting UTC -> JD(UT1) -> UTC should give original date."""
        orig_year, orig_month, orig_day = 2020, 6, 15
        orig_hour, orig_minute, orig_second = 14, 30, 45.5

        # Convert to JD (get UT1 component)
        jd_tt, jd_ut1 = ephem.utc_to_jd(
            orig_year, orig_month, orig_day, orig_hour, orig_minute, orig_second
        )

        # Convert back to UTC using JD(UT1)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        # Should match original (within DUT1 tolerance < 0.9 seconds)
        assert year == orig_year
        assert month == orig_month
        assert day == orig_day
        assert hour == orig_hour
        assert minute == orig_minute
        # DUT1 is always < 0.9 seconds
        assert abs(second - orig_second) < 1.0

    @pytest.mark.unit
    def test_roundtrip_j2000(self):
        """Roundtrip for J2000.0 epoch."""
        orig_year, orig_month, orig_day = 2000, 1, 1
        orig_hour, orig_minute, orig_second = 12, 0, 0.0

        jd_tt, jd_ut1 = ephem.utc_to_jd(
            orig_year, orig_month, orig_day, orig_hour, orig_minute, orig_second
        )
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        assert year == orig_year
        assert month == orig_month
        assert day == orig_day
        assert hour == orig_hour
        assert minute == orig_minute
        assert abs(second - orig_second) < 1.0

    @pytest.mark.unit
    def test_roundtrip_multiple_dates(self):
        """Roundtrip for various dates across history."""
        test_dates = [
            (1990, 7, 4, 12, 30, 0.0),  # Summer date
            (2010, 12, 31, 23, 59, 30.0),  # Year end
            (2024, 2, 29, 6, 15, 22.5),  # Leap year
        ]

        for (
            orig_year,
            orig_month,
            orig_day,
            orig_hour,
            orig_minute,
            orig_second,
        ) in test_dates:
            jd_tt, jd_ut1 = ephem.utc_to_jd(
                orig_year, orig_month, orig_day, orig_hour, orig_minute, orig_second
            )
            year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

            assert year == orig_year, (
                f"Year mismatch for {orig_year}-{orig_month}-{orig_day}"
            )
            assert month == orig_month, (
                f"Month mismatch for {orig_year}-{orig_month}-{orig_day}"
            )
            assert day == orig_day, (
                f"Day mismatch for {orig_year}-{orig_month}-{orig_day}"
            )
            # Compare total seconds to handle minute boundary shifts from DUT1
            orig_total_seconds = orig_hour * 3600 + orig_minute * 60 + orig_second
            result_total_seconds = hour * 3600 + minute * 60 + second
            # DUT1 < 0.9 seconds, so time should match within 1 second
            assert abs(result_total_seconds - orig_total_seconds) < 1.0, (
                f"Time mismatch for {orig_year}-{orig_month}-{orig_day}"
            )


class TestJdut1ToUtcDut1:
    """Test that UT1-UTC difference (DUT1) is properly handled."""

    @pytest.mark.unit
    def test_dut1_always_less_than_0_9_seconds(self):
        """DUT1 = UT1 - UTC should always be < 0.9 seconds by definition."""
        # Test several dates
        test_jds = [
            ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)[1],  # JD(UT1) around J2000
            ephem.utc_to_jd(2010, 6, 15, 0, 0, 0.0)[1],  # Mid 2010
            ephem.utc_to_jd(2020, 3, 20, 6, 30, 0.0)[1],  # Equinox 2020
        ]

        for jd_ut1 in test_jds:
            year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

            # Get the UTC time in seconds from midnight
            utc_seconds = hour * 3600 + minute * 60 + second

            # Get the UT1 time by reversing the JD
            y, m, d, decimal_hour = ephem.revjul(jd_ut1)
            ut1_seconds = decimal_hour * 3600

            # The difference should be very small (essentially the DUT1 value)
            # Since we're doing a roundtrip, the difference should be negligible
            # but theoretically |UT1 - UTC| < 0.9 seconds
            # In practice for roundtrip it should be near zero
            diff = abs(ut1_seconds - utc_seconds)
            # Allow for date boundary crossings
            if diff > 43000:  # Near midnight crossing
                diff = 86400 - diff
            assert diff < 1.0, f"DUT1 too large: {diff}s at JD {jd_ut1}"


class TestJdut1ToUtcCalendars:
    """Test Gregorian vs Julian calendar handling."""

    @pytest.mark.unit
    def test_gregorian_default(self):
        """Gregorian calendar should be default."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 6, 15, 12, 0, 0.0)
        result_default = ephem.jdut1_to_utc(jd_ut1)
        result_greg = ephem.jdut1_to_utc(jd_ut1, GREG_CAL)
        assert result_default == result_greg

    @pytest.mark.unit
    def test_julian_calendar_differs(self):
        """Julian calendar should give different date for modern dates."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)
        year_greg, month_greg, day_greg, _, _, _ = ephem.jdut1_to_utc(
            jd_ut1, GREG_CAL
        )
        year_jul, month_jul, day_jul, _, _, _ = ephem.jdut1_to_utc(jd_ut1, JUL_CAL)

        # In 2000, Julian calendar is 13 days behind Gregorian
        # Gregorian 2000-01-01 = Julian 1999-12-19
        assert year_greg == 2000
        assert year_jul == 1999
        assert month_greg == 1
        assert month_jul == 12
        assert day_greg == 1
        assert day_jul == 19

    @pytest.mark.unit
    def test_julian_roundtrip(self):
        """Roundtrip should work for Julian calendar dates."""
        orig_year, orig_month, orig_day = 1999, 12, 19
        orig_hour, orig_minute, orig_second = 12, 0, 0.0

        jd_tt, jd_ut1 = ephem.utc_to_jd(
            orig_year,
            orig_month,
            orig_day,
            orig_hour,
            orig_minute,
            orig_second,
            calendar=JUL_CAL,
        )
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1, JUL_CAL)

        assert year == orig_year
        assert month == orig_month
        assert day == orig_day
        assert hour == orig_hour
        assert minute == orig_minute
        # Within DUT1 tolerance
        assert abs(second - orig_second) < 1.0


class TestJdut1ToUtcVsJdetToUtc:
    """Compare jdut1_to_utc with jdet_to_utc to verify Delta-T handling."""

    @pytest.mark.unit
    def test_ut1_and_tt_give_different_results(self):
        """JD(UT1) and JD(TT) for same instant should give same UTC."""
        # Get both JD(TT) and JD(UT1) for the same UTC instant
        utc_date = (2020, 6, 15, 14, 30, 0.0)
        jd_tt, jd_ut1 = ephem.utc_to_jd(*utc_date)

        # Convert both back to UTC
        result_from_tt = ephem.jdet_to_utc(jd_tt)
        result_from_ut1 = ephem.jdut1_to_utc(jd_ut1)

        # Both should return the same UTC time (within tolerance)
        assert result_from_tt[0] == result_from_ut1[0]  # year
        assert result_from_tt[1] == result_from_ut1[1]  # month
        assert result_from_tt[2] == result_from_ut1[2]  # day
        # Compare total seconds to account for minute boundary shifts from DUT1
        tt_total_sec = (
            result_from_tt[3] * 3600 + result_from_tt[4] * 60 + result_from_tt[5]
        )
        ut1_total_sec = (
            result_from_ut1[3] * 3600 + result_from_ut1[4] * 60 + result_from_ut1[5]
        )
        # Should match within DUT1 tolerance (<0.9 seconds)
        assert abs(tt_total_sec - ut1_total_sec) < 1.0

    @pytest.mark.unit
    def test_delta_t_difference(self):
        """Verify that JD(TT) and JD(UT1) differ by Delta-T."""
        utc_date = (2020, 6, 15, 12, 0, 0.0)
        jd_tt, jd_ut1 = ephem.utc_to_jd(*utc_date)

        # Delta-T = TT - UT1 (in days)
        delta_t_days = jd_tt - jd_ut1
        delta_t_seconds = delta_t_days * 86400

        # Delta-T around 2020 is approximately 69 seconds
        assert 60 < delta_t_seconds < 80


class TestJdut1ToUtcEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.unit
    def test_year_boundary(self):
        """Test New Year transition."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 1, 1, 0, 0, 0.5)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        assert year == 2020
        assert month == 1
        assert day == 1
        assert hour == 0
        assert minute == 0
        assert abs(second - 0.5) < 1.0

    @pytest.mark.unit
    def test_leap_year_feb29(self):
        """Test Feb 29 in leap year."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 2, 29, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        assert year == 2020
        assert month == 2
        assert day == 29
        assert hour == 12
        assert minute == 0
        assert abs(second) < 1.0

    @pytest.mark.unit
    def test_historical_date(self):
        """Test historical date (before UTC was defined)."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(1800, 6, 15, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        assert year == 1800
        assert month == 6
        assert day == 15
        # Hour should be close to 12
        total_seconds = hour * 3600 + minute * 60 + second
        expected_seconds = 12 * 3600
        # Should be within a reasonable tolerance
        assert abs(total_seconds - expected_seconds) < 60

    @pytest.mark.unit
    def test_future_date(self):
        """Test future date."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2100, 1, 1, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        assert year == 2100
        assert month == 1
        assert day == 1
        # Compare total seconds to account for DUT1 uncertainty in future
        total_seconds = hour * 3600 + minute * 60 + second
        expected_seconds = 12 * 3600
        # Future DUT1 is predicted, so allow reasonable tolerance
        assert abs(total_seconds - expected_seconds) < 60

    @pytest.mark.unit
    def test_fractional_seconds_precision(self):
        """Test fractional seconds are preserved."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 6, 15, 14, 30, 45.123)
        year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

        # Within DUT1 tolerance
        assert abs(second - 45.123) < 1.0


class TestJdut1ToUtcConsistency:
    """Test internal consistency."""

    @pytest.mark.unit
    def test_consecutive_seconds(self):
        """Consecutive JD values should give consecutive UTC times."""
        one_second_jd = 1.0 / 86400.0

        jd_tt_0, jd_ut1_0 = ephem.utc_to_jd(2023, 2, 24, 12, 0, 0.0)
        jd_ut1_1 = jd_ut1_0 + one_second_jd

        _, _, _, h0, m0, s0 = ephem.jdut1_to_utc(jd_ut1_0)
        _, _, _, h1, m1, s1 = ephem.jdut1_to_utc(jd_ut1_1)

        utc_seconds_0 = h0 * 3600 + m0 * 60 + s0
        utc_seconds_1 = h1 * 3600 + m1 * 60 + s1

        # Should differ by 1 second
        diff = utc_seconds_1 - utc_seconds_0
        assert diff == pytest.approx(1.0, abs=0.01)

    @pytest.mark.unit
    def test_consecutive_days(self):
        """JD values 1 day apart should give dates 1 day apart."""
        jd_tt_0, jd_ut1_0 = ephem.utc_to_jd(2023, 2, 24, 12, 0, 0.0)
        jd_ut1_1 = jd_ut1_0 + 1.0

        y0, m0, d0, _, _, _ = ephem.jdut1_to_utc(jd_ut1_0)
        y1, m1, d1, _, _, _ = ephem.jdut1_to_utc(jd_ut1_1)

        # Convert to JD for comparison (ignoring time)
        jd_date_0 = ephem.julday(y0, m0, d0, 0.0)
        jd_date_1 = ephem.julday(y1, m1, d1, 0.0)

        # Should differ by 1 day
        assert jd_date_1 - jd_date_0 == pytest.approx(1.0, abs=0.01)

    @pytest.mark.unit
    def test_valid_time_components(self):
        """Returned time components should be in valid ranges."""
        test_dates = [
            (2000, 1, 1, 12, 0, 0.0),
            (2023, 2, 24, 12, 0, 0.0),
            (2023, 2, 25, 0, 0, 0.0),
            (1900, 1, 1, 12, 0, 0.0),
        ]

        for date in test_dates:
            jd_tt, jd_ut1 = ephem.utc_to_jd(*date)
            year, month, day, hour, minute, second = ephem.jdut1_to_utc(jd_ut1)

            assert 1 <= month <= 12, f"Invalid month {month} for {date}"
            assert 1 <= day <= 31, f"Invalid day {day} for {date}"
            assert 0 <= hour <= 23, f"Invalid hour {hour} for {date}"
            assert 0 <= minute <= 59, f"Invalid minute {minute} for {date}"
            # Second can be 0-59.999..., or up to 60.x during leap second
            assert 0 <= second < 61, f"Invalid second {second} for {date}"
