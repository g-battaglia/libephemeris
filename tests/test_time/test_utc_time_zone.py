"""
Comprehensive tests for utc_time_zone() function.

Reference ephemeris semantics: utc_time_zone converts LOCAL time to UTC by
SUBTRACTING the timezone offset.  E.g. 10:00 CET (UTC+1) → 09:00 UTC.

Tests cover:
- Basic timezone offset application (local → UTC)
- Positive offsets (east of UTC: CET, JST, etc.)
- Negative offsets (west of UTC: EST, PST, etc.)
- Day boundary crossings
- Month boundary crossings
- Year boundary crossings
- Fractional hour offsets
- Edge cases (leap years, UTC+0)
"""

import pytest
import libephemeris as ephem


class TestUtcTimeZoneBasic:
    """Test basic timezone offset application."""

    @pytest.mark.unit
    def test_returns_tuple(self):
        """Function should return a tuple of 6 elements."""
        result = ephem.utc_time_zone(2024, 1, 15, 10, 30, 0.0, 1)
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
    def test_zero_offset(self):
        """UTC+0 should return the same time."""
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 30, 45.5, 0
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 10
        assert minute == 30
        assert second == pytest.approx(45.5, abs=0.001)

    @pytest.mark.unit
    def test_positive_offset_no_boundary(self):
        """Positive offset subtracts hours (local→UTC)."""
        # Local 11:30 CET (UTC+1) → UTC 10:30
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 11, 30, 0.0, 1
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 10
        assert minute == 30
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_offset_no_boundary(self):
        """Negative offset adds hours (local→UTC)."""
        # Local 10:30 EST (UTC-5) → UTC 15:30
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 10, 30, 0.0, -5
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 15
        assert minute == 30
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneDayBoundary:
    """Test day boundary crossings."""

    @pytest.mark.unit
    def test_positive_offset_previous_day(self):
        """Positive offset crossing into previous day (local→UTC)."""
        # Local 01:00 UTC+2 → UTC 23:00 previous day
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 16, 1, 0, 0.0, 2
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 23
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_offset_next_day(self):
        """Negative offset crossing into next day (local→UTC)."""
        # Local 21:00 EST (UTC-5) → UTC 02:00 next day
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 21, 0, 0.0, -5
        )
        assert year == 2024
        assert month == 1
        assert day == 16
        assert hour == 2
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneMonthBoundary:
    """Test month boundary crossings."""

    @pytest.mark.unit
    def test_positive_offset_previous_month(self):
        """Positive offset crossing into previous month (local→UTC)."""
        # Local 01:00 Feb 1 UTC+2 → UTC 23:00 Jan 31
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 2, 1, 1, 0, 0.0, 2
        )
        assert year == 2024
        assert month == 1
        assert day == 31
        assert hour == 23
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_offset_next_month(self):
        """Negative offset crossing into next month (local→UTC)."""
        # Local 21:00 Jan 31 EST (UTC-5) → UTC 02:00 Feb 1
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 31, 21, 0, 0.0, -5
        )
        assert year == 2024
        assert month == 2
        assert day == 1
        assert hour == 2
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneYearBoundary:
    """Test year boundary crossings."""

    @pytest.mark.unit
    def test_positive_offset_previous_year(self):
        """Positive offset crossing into previous year (local→UTC)."""
        # Local 01:00 Jan 1 2024 UTC+2 → UTC 23:00 Dec 31 2023
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 1, 1, 0, 0.0, 2
        )
        assert year == 2023
        assert month == 12
        assert day == 31
        assert hour == 23
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_offset_next_year(self):
        """Negative offset crossing into next year (local→UTC)."""
        # Local 21:00 Dec 31 2023 EST (UTC-5) → UTC 02:00 Jan 1 2024
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2023, 12, 31, 21, 0, 0.0, -5
        )
        assert year == 2024
        assert month == 1
        assert day == 1
        assert hour == 2
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneFractionalOffsets:
    """Test fractional hour timezone offsets."""

    @pytest.mark.unit
    def test_half_hour_offset(self):
        """Test 30-minute offset (e.g., India UTC+5:30)."""
        # Local 15:30 IST (UTC+5.5) → UTC 10:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 15, 30, 0.0, 5.5
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 10
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_quarter_hour_offset(self):
        """Test 45-minute offset (e.g., Nepal UTC+5:45)."""
        # Local 15:45 Nepal (UTC+5.75) → UTC 10:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 15, 45, 0.0, 5.75
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 10
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_negative_half_hour_offset(self):
        """Test negative 30-minute offset (e.g., Newfoundland UTC-3:30)."""
        # Local 11:30 NST (UTC-3.5) → UTC 15:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 11, 30, 0.0, -3.5
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 15
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)


class TestUtcTimeZoneCommonTimezones:
    """Test common real-world timezone conversions (local→UTC)."""

    @pytest.mark.unit
    def test_cet_to_utc(self):
        """Test CET (UTC+1) to UTC."""
        # Local 13:00 CET → UTC 12:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 6, 15, 13, 0, 0.0, 1
        )
        assert hour == 12
        assert minute == 0

    @pytest.mark.unit
    def test_est_to_utc(self):
        """Test EST (UTC-5) to UTC."""
        # Local 13:00 EST → UTC 18:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 13, 0, 0.0, -5
        )
        assert hour == 18
        assert minute == 0

    @pytest.mark.unit
    def test_pst_to_utc(self):
        """Test PST (UTC-8) to UTC."""
        # Local 12:00 PST → UTC 20:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 12, 0, 0.0, -8
        )
        assert hour == 20
        assert minute == 0

    @pytest.mark.unit
    def test_jst_to_utc(self):
        """Test JST (UTC+9) to UTC."""
        # Local 12:00 JST → UTC 03:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 12, 0, 0.0, 9
        )
        assert hour == 3
        assert minute == 0


class TestUtcTimeZoneEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.mark.unit
    def test_leap_year_feb29(self):
        """Test Feb 29 in leap year."""
        # Local 17:00 Feb 29 UTC+5 → UTC 12:00
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 2, 29, 17, 0, 0.0, 5
        )
        assert year == 2024
        assert month == 2
        assert day == 29
        assert hour == 12
        assert minute == 0

    @pytest.mark.unit
    def test_leap_year_crossing_from_mar1(self):
        """Test crossing from Mar 1 back to Feb 29 in leap year."""
        # Local 01:00 Mar 1 UTC+2 → UTC 23:00 Feb 29
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 3, 1, 1, 0, 0.0, 2
        )
        assert year == 2024
        assert month == 2
        assert day == 29
        assert hour == 23
        assert minute == 0

    @pytest.mark.unit
    def test_fractional_seconds_preserved(self):
        """Test that fractional seconds are preserved."""
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 11, 30, 45.123, 1
        )
        assert second == pytest.approx(45.123, abs=0.01)

    @pytest.mark.unit
    def test_midnight_local_negative_offset(self):
        """Test midnight local with negative offset."""
        # Local 00:00 Jan 15 EST (UTC-5) → UTC 05:00 Jan 15
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 15, 0, 0, 0.0, -5
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 5
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_large_positive_offset(self):
        """Test UTC+14 (Line Islands, Kiribati)."""
        # Local 00:00 Jan 16 UTC+14 → UTC 10:00 Jan 15
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 16, 0, 0, 0.0, 14
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 10
        assert minute == 0

    @pytest.mark.unit
    def test_large_negative_offset(self):
        """Test UTC-12 (Baker Island)."""
        # Local 22:00 Jan 14 UTC-12 → UTC 10:00 Jan 15
        year, month, day, hour, minute, second = ephem.utc_time_zone(
            2024, 1, 14, 22, 0, 0.0, -12
        )
        assert year == 2024
        assert month == 1
        assert day == 15
        assert hour == 10
        assert minute == 0


class TestUtcTimeZoneConsistency:
    """Test internal consistency."""

    @pytest.mark.unit
    def test_roundtrip_offset(self):
        """Applying +N then -N offset should return original time."""
        orig_year, orig_month, orig_day = 2024, 6, 15
        orig_hour, orig_minute, orig_second = 12, 30, 45.5

        # Apply +5 offset (local→UTC: subtract 5)
        y1, m1, d1, h1, min1, s1 = ephem.utc_time_zone(
            orig_year, orig_month, orig_day, orig_hour, orig_minute, orig_second, 5
        )

        # Apply -5 offset (local→UTC: add 5) to get back
        y2, m2, d2, h2, min2, s2 = ephem.utc_time_zone(y1, m1, d1, h1, min1, s1, -5)

        assert y2 == orig_year
        assert m2 == orig_month
        assert d2 == orig_day
        assert h2 == orig_hour
        assert min2 == orig_minute
        assert s2 == pytest.approx(orig_second, abs=0.01)

    @pytest.mark.unit
    def test_consecutive_offsets(self):
        """Test that consecutive integer offsets give consecutive hours."""
        base_hour = 12
        for offset in range(-5, 6):
            year, month, day, hour, minute, second = ephem.utc_time_zone(
                2024, 1, 15, base_hour, 0, 0.0, offset
            )
            # local→UTC: hour = base_hour - offset
            expected_hour = (base_hour - offset) % 24
            assert hour == expected_hour, f"Failed for offset {offset}"
