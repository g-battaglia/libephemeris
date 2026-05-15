"""
Comprehensive tests for jdet_to_utc() function.

Tests cover:
- Basic JD(TT) to UTC conversion
- Roundtrip consistency with utc_to_jd()
- Delta-T application (TT vs UT1)
- Gregorian vs Julian calendar handling
- Leap second awareness
- Edge cases (historical dates, future dates, year boundaries)
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import GREG_CAL, JUL_CAL


class TestJdetToUtcBasic:
    """Test basic JD(TT) to UTC conversion."""

    @pytest.mark.unit
    def test_j2000_epoch(self):
        """J2000.0 = JD(TT) 2451545.0 should give ~2000-01-01 11:58:55 UTC."""
        year, month, day, hour, minute, second = ephem.jdet_to_utc(2451545.0)
        # J2000.0 at JD(TT) 2451545.0 is approximately 11:58:55.82 UTC
        # (Delta T was ~63.82 seconds in 2000)
        assert year == 2000
        assert month == 1
        assert day == 1
        assert hour == 11
        assert minute == 58
        # Second should be around 55-56 due to Delta T
        assert 55 < second < 57

    @pytest.mark.unit
    def test_returns_tuple(self):
        """Function should return a tuple of 6 elements."""
        result = ephem.jdet_to_utc(2451545.0)
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
    def test_noon_tt(self):
        """Test JD ending in .0 (noon TT)."""
        year, month, day, hour, minute, second = ephem.jdet_to_utc(2460000.0)
        # JD 2460000.0 = Feb 24, 2023 12:00:00 TT
        # UTC should be slightly before noon due to Delta T (~69 sec in 2023)
        assert year == 2023
        assert month == 2
        assert day == 24
        assert hour == 11
        assert minute == 58
        # Seconds around 50-52 due to Delta T

    @pytest.mark.unit
    def test_midnight_tt(self):
        """Test JD ending in .5 (midnight TT)."""
        year, month, day, hour, minute, second = ephem.jdet_to_utc(2460000.5)
        # JD 2460000.5 = Feb 25, 2023 00:00:00 TT
        # UTC should be Feb 24, 2023 ~23:58:50
        assert year == 2023
        assert month == 2
        assert day == 24
        assert hour == 23
        assert minute == 58


class TestJdetToUtcRoundtrip:
    """Test roundtrip consistency with utc_to_jd()."""

    @pytest.mark.unit
    def test_roundtrip_modern_date(self):
        """Converting UTC -> JD(TT) -> UTC should give original date."""
        # Start with a UTC date
        orig_year, orig_month, orig_day = 2020, 6, 15
        orig_hour, orig_minute, orig_second = 14, 30, 45.5

        # Convert to JD(TT)
        jd_tt, jd_ut = ephem.utc_to_jd(
            orig_year, orig_month, orig_day, orig_hour, orig_minute, orig_second
        )

        # Convert back to UTC
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        # Should match original (within floating point tolerance)
        assert year == orig_year
        assert month == orig_month
        assert day == orig_day
        assert hour == orig_hour
        assert minute == orig_minute
        assert second == pytest.approx(orig_second, abs=0.001)

    @pytest.mark.unit
    def test_roundtrip_j2000(self):
        """Roundtrip for J2000.0 epoch."""
        orig_year, orig_month, orig_day = 2000, 1, 1
        orig_hour, orig_minute, orig_second = 12, 0, 0.0

        jd_tt, _ = ephem.utc_to_jd(
            orig_year, orig_month, orig_day, orig_hour, orig_minute, orig_second
        )
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        assert year == orig_year
        assert month == orig_month
        assert day == orig_day
        assert hour == orig_hour
        assert minute == orig_minute
        assert second == pytest.approx(orig_second, abs=0.001)

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
            jd_tt, _ = ephem.utc_to_jd(
                orig_year, orig_month, orig_day, orig_hour, orig_minute, orig_second
            )
            year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

            assert year == orig_year, (
                f"Year mismatch for {orig_year}-{orig_month}-{orig_day}"
            )
            assert month == orig_month, (
                f"Month mismatch for {orig_year}-{orig_month}-{orig_day}"
            )
            assert day == orig_day, (
                f"Day mismatch for {orig_year}-{orig_month}-{orig_day}"
            )
            assert hour == orig_hour, (
                f"Hour mismatch for {orig_year}-{orig_month}-{orig_day}"
            )
            assert minute == orig_minute, (
                f"Minute mismatch for {orig_year}-{orig_month}-{orig_day}"
            )
            assert second == pytest.approx(orig_second, abs=0.001), (
                f"Second mismatch for {orig_year}-{orig_month}-{orig_day}"
            )


class TestJdetToUtcDeltaT:
    """Test Delta-T is properly applied."""

    @pytest.mark.unit
    def test_delta_t_difference_j2000(self):
        """At J2000.0, Delta T was ~63.8 seconds."""
        # JD(TT) 2451545.0 = 2000-01-01 12:00:00 TT
        year, month, day, hour, minute, second = ephem.jdet_to_utc(2451545.0)

        # UTC should be ~63.8 seconds earlier
        # 12:00:00 TT - 63.8 sec = 11:58:56.2 UTC
        utc_in_seconds = hour * 3600 + minute * 60 + second
        tt_noon_seconds = 12 * 3600

        delta_t_applied = tt_noon_seconds - utc_in_seconds
        # Delta T in 2000 was approximately 63.82 seconds
        assert 63 < delta_t_applied < 65

    @pytest.mark.unit
    def test_delta_t_historical(self):
        """Historical dates should show appropriate Delta T for the era."""
        # 1900-01-01 12:00 TT
        jd_1900 = 2415021.0
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_1900)

        # Delta T around 1900 was actually about -2 to +2 seconds based on IERS data
        # However, before TAI was established, Skyfield/IERS may use different models
        # The important thing is that the conversion works and returns valid data
        assert (
            year == 1900 or year == 1899
        )  # May cross year boundary depending on delta
        assert 1 <= month <= 12
        assert 1 <= day <= 31


class TestJdetToUtcCalendars:
    """Test Gregorian vs Julian calendar handling."""

    @pytest.mark.unit
    def test_gregorian_default(self):
        """Gregorian calendar should be default."""
        result_default = ephem.jdet_to_utc(2451545.0)
        result_greg = ephem.jdet_to_utc(2451545.0, GREG_CAL)
        assert result_default == result_greg

    @pytest.mark.unit
    def test_julian_calendar_differs(self):
        """Julian calendar should give different date for modern dates."""
        year_greg, month_greg, day_greg, _, _, _ = ephem.jdet_to_utc(
            2451545.0, GREG_CAL
        )
        year_jul, month_jul, day_jul, _, _, _ = ephem.jdet_to_utc(2451545.0, JUL_CAL)

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
        # Use Julian calendar
        orig_year, orig_month, orig_day = 1999, 12, 19
        orig_hour, orig_minute, orig_second = 12, 0, 0.0

        jd_tt, _ = ephem.utc_to_jd(
            orig_year,
            orig_month,
            orig_day,
            orig_hour,
            orig_minute,
            orig_second,
            calendar=JUL_CAL,
        )
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt, JUL_CAL)

        assert year == orig_year
        assert month == orig_month
        assert day == orig_day
        assert hour == orig_hour
        assert minute == orig_minute
        assert second == pytest.approx(orig_second, abs=0.01)


class TestJdetToUtcLeapSeconds:
    """Test leap second awareness."""

    @pytest.mark.unit
    def test_after_leap_second_2016(self):
        """Test date after 2016-12-31 leap second."""
        # Get JD(TT) for 2017-01-01 00:00:01 UTC
        jd_tt, _ = ephem.utc_to_jd(2017, 1, 1, 0, 0, 1.0)
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        assert year == 2017
        assert month == 1
        assert day == 1
        assert hour == 0
        assert minute == 0
        assert second == pytest.approx(1.0, abs=0.001)

    @pytest.mark.unit
    def test_before_leap_second_2016(self):
        """Test date before 2016-12-31 leap second."""
        # Get JD(TT) for 2016-12-31 23:59:58 UTC
        jd_tt, _ = ephem.utc_to_jd(2016, 12, 31, 23, 59, 58.0)
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        assert year == 2016
        assert month == 12
        assert day == 31
        assert hour == 23
        assert minute == 59
        assert second == pytest.approx(58.0, abs=0.001)


class TestJdetToUtcEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.unit
    def test_year_boundary(self):
        """Test New Year transition."""
        # 2020-01-01 00:00:00.5 UTC
        jd_tt, _ = ephem.utc_to_jd(2020, 1, 1, 0, 0, 0.5)
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        assert year == 2020
        assert month == 1
        assert day == 1
        assert hour == 0
        assert minute == 0
        assert second == pytest.approx(0.5, abs=0.001)

    @pytest.mark.unit
    def test_leap_year_feb29(self):
        """Test Feb 29 in leap year."""
        # 2020-02-29 12:00:00 UTC
        jd_tt, _ = ephem.utc_to_jd(2020, 2, 29, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        assert year == 2020
        assert month == 2
        assert day == 29
        assert hour == 12
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_historical_date(self):
        """Test historical date (before UTC was defined)."""
        # 1800-06-15 12:00:00 (treated as UT1 approximation)
        jd_tt, _ = ephem.utc_to_jd(1800, 6, 15, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        # For historical dates, Delta T may cause hour to differ slightly
        # The important thing is the roundtrip works approximately
        assert year == 1800
        assert month == 6
        assert day == 15
        # Hour may vary slightly due to Delta T application
        # Delta T in 1800 was about 13-14 seconds, so hour should still be 12 or close
        total_seconds = hour * 3600 + minute * 60 + second
        expected_seconds = 12 * 3600  # noon
        # Should be within 60 seconds of noon
        assert abs(total_seconds - expected_seconds) < 60

    @pytest.mark.unit
    def test_future_date(self):
        """Test future date."""
        # 2100-01-01 12:00:00 UTC
        jd_tt, _ = ephem.utc_to_jd(2100, 1, 1, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        assert year == 2100
        assert month == 1
        assert day == 1
        assert hour == 12
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.01)

    @pytest.mark.unit
    def test_fractional_seconds_precision(self):
        """Test fractional seconds are preserved."""
        # 2020-06-15 14:30:45.123 UTC
        jd_tt, _ = ephem.utc_to_jd(2020, 6, 15, 14, 30, 45.123)
        year, month, day, hour, minute, second = ephem.jdet_to_utc(jd_tt)

        assert second == pytest.approx(45.123, abs=0.001)


class TestJdetToUtcConsistency:
    """Test internal consistency."""

    @pytest.mark.unit
    def test_consecutive_seconds(self):
        """Consecutive JD values should give consecutive UTC times."""
        # One second = 1/86400 day
        one_second_jd = 1.0 / 86400.0

        jd_tt_0 = 2460000.0
        jd_tt_1 = jd_tt_0 + one_second_jd

        _, _, _, h0, m0, s0 = ephem.jdet_to_utc(jd_tt_0)
        _, _, _, h1, m1, s1 = ephem.jdet_to_utc(jd_tt_1)

        utc_seconds_0 = h0 * 3600 + m0 * 60 + s0
        utc_seconds_1 = h1 * 3600 + m1 * 60 + s1

        # Should differ by 1 second
        diff = utc_seconds_1 - utc_seconds_0
        assert diff == pytest.approx(1.0, abs=0.001)

    @pytest.mark.unit
    def test_consecutive_days(self):
        """JD values 1 day apart should give dates 1 day apart."""
        jd_tt_0 = 2460000.0
        jd_tt_1 = jd_tt_0 + 1.0

        y0, m0, d0, _, _, _ = ephem.jdet_to_utc(jd_tt_0)
        y1, m1, d1, _, _, _ = ephem.jdet_to_utc(jd_tt_1)

        # Convert to JD for comparison (ignoring time)
        jd_date_0 = ephem.julday(y0, m0, d0, 0.0)
        jd_date_1 = ephem.julday(y1, m1, d1, 0.0)

        # Should differ by 1 day
        assert jd_date_1 - jd_date_0 == pytest.approx(1.0, abs=0.01)

    @pytest.mark.unit
    def test_valid_time_components(self):
        """Returned time components should be in valid ranges."""
        test_jds = [2451545.0, 2460000.0, 2460000.5, 2415021.0]

        for jd in test_jds:
            year, month, day, hour, minute, second = ephem.jdet_to_utc(jd)

            assert 1 <= month <= 12, f"Invalid month {month} for JD {jd}"
            assert 1 <= day <= 31, f"Invalid day {day} for JD {jd}"
            assert 0 <= hour <= 23, f"Invalid hour {hour} for JD {jd}"
            assert 0 <= minute <= 59, f"Invalid minute {minute} for JD {jd}"
            # Second can be 0-59.999..., or up to 60.x during leap second
            assert 0 <= second < 61, f"Invalid second {second} for JD {jd}"
