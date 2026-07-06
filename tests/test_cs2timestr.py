"""
Tests for the cs2timestr centiseconds to time string conversion function.

Tests verify that cs2timestr correctly converts time values in centiseconds to
formatted time strings (e.g., "12:34:56").
"""

import libephemeris as ephem


# Constants for centiseconds calculations (time-based)
CS_HOUR = 360000  # centiseconds per hour (60 * 60 * 100)
CS_MINUTE = 6000  # centiseconds per minute (60 * 100)
CS_SECOND = 100  # centiseconds per second


class TestCs2timestrBasic:
    """Basic functionality tests for cs2timestr."""

    def test_cs2timestr_exported(self):
        """Test that cs2timestr is exported from the package."""
        assert hasattr(ephem, "cs2timestr")
        assert callable(ephem.cs2timestr)

    def test_cs2timestr_zero(self):
        """Test conversion of zero."""
        result = ephem.cs2timestr(0)
        assert "0" in result
        assert ":" in result

    def test_cs2timestr_returns_string(self):
        """Test that cs2timestr returns a string."""
        result = ephem.cs2timestr(CS_HOUR)
        assert isinstance(result, str)

    def test_cs2timestr_one_hour(self):
        """Test conversion of exactly 1 hour."""
        result = ephem.cs2timestr(CS_HOUR)
        # Should be " 1:00:00"
        assert "1:" in result
        assert ":00:00" in result

    def test_cs2timestr_contains_colons(self):
        """Test that result contains the colon separator."""
        result = ephem.cs2timestr(CS_HOUR * 12)
        assert result.count(":") == 2


class TestCs2timestrValues:
    """Tests for various value conversions."""

    def test_cs2timestr_12_hours(self):
        """Test conversion of 12 hours."""
        result = ephem.cs2timestr(CS_HOUR * 12)
        assert "12:" in result

    def test_cs2timestr_23_hours(self):
        """Test conversion of 23 hours."""
        result = ephem.cs2timestr(CS_HOUR * 23)
        assert "23:" in result

    def test_cs2timestr_with_minutes(self):
        """Test conversion with minutes component."""
        # 1 hour 30 minutes = 360000 + 30*6000 = 360000 + 180000 = 540000 cs
        cs_value = CS_HOUR + (30 * CS_MINUTE)
        result = ephem.cs2timestr(cs_value)
        assert "1:" in result
        assert ":30:" in result

    def test_cs2timestr_with_seconds(self):
        """Test conversion with seconds component."""
        # 1 hour 0 minutes 45 seconds = 360000 + 0 + 4500 = 364500 cs
        cs_value = CS_HOUR + (45 * CS_SECOND)
        result = ephem.cs2timestr(cs_value)
        assert "1:" in result
        assert ":45" in result

    def test_cs2timestr_complex_value(self):
        """Test conversion of a complex value with all components."""
        # 12 hours 34 minutes 56 seconds
        # = 12*360000 + 34*6000 + 56*100
        # = 4320000 + 204000 + 5600
        # = 4529600
        cs_value = 12 * CS_HOUR + 34 * CS_MINUTE + 56 * CS_SECOND
        result = ephem.cs2timestr(cs_value)
        assert "12:" in result
        assert ":34:" in result
        assert ":56" in result

    def test_cs2timestr_midnight(self):
        """Test conversion of midnight (0:00:00)."""
        result = ephem.cs2timestr(0)
        assert ":00:00" in result


class TestCs2timestrNegativeValues:
    """Tests for negative value handling (hours wrap mod 24)."""

    def test_cs2timestr_negative_1_hour(self):
        """Test conversion of -1 hour (wraps to 23:00:00 mod 24)."""
        result = ephem.cs2timestr(-CS_HOUR)
        assert "23:00:00" in result

    def test_cs2timestr_negative_12_hours(self):
        """Test conversion of -12 hours (wraps to 12:00:00 mod 24)."""
        result = ephem.cs2timestr(-CS_HOUR * 12)
        assert "12:" in result

    def test_cs2timestr_negative_with_minutes(self):
        """Test conversion of negative value with minutes (-1h30m wraps to 22:30)."""
        cs_value = -(CS_HOUR + 30 * CS_MINUTE)  # -1 hour 30 minutes
        result = ephem.cs2timestr(cs_value)
        # The WHOLE value wraps modulo 24h (documented contract), so
        # -1.5h == 22.5h -> 22:30:00. (Applying the sign to the hours field
        # alone would give 23:30 but breaks sub-hour negatives; the reference
        # API's output for negative cs is undefined garbage, so the total
        # wrap is this library's sane, documented behavior.)
        assert "22:30:00" in result


class TestCs2timestrLargeValues:
    """Tests for large value handling (hours wrap mod 24)."""

    def test_cs2timestr_24_hours(self):
        """Test conversion of 24 hours (wraps to 00:00:00 mod 24)."""
        result = ephem.cs2timestr(CS_HOUR * 24)
        assert "00:00:00" in result

    def test_cs2timestr_greater_than_24(self):
        """Test conversion of value greater than 24 hours (30h wraps to 06:00:00)."""
        result = ephem.cs2timestr(CS_HOUR * 30)  # 30 hours
        assert "06:" in result

    def test_cs2timestr_100_hours(self):
        """Test conversion of 100 hours (wraps to 04:00:00 mod 24)."""
        result = ephem.cs2timestr(CS_HOUR * 100)
        assert "04:" in result


class TestCs2timestrEdgeCases:
    """Edge case tests for cs2timestr."""

    def test_cs2timestr_1_centisecond(self):
        """Test conversion of 1 centisecond (rounds to 0 seconds)."""
        result = ephem.cs2timestr(1)
        assert ":00:00" in result

    def test_cs2timestr_50_centiseconds(self):
        """Test conversion of 50 centiseconds (rounds to 1 second)."""
        result = ephem.cs2timestr(50)
        # 50 cs rounds to 1 second
        assert ":00:01" in result

    def test_cs2timestr_99_centiseconds(self):
        """Test conversion of 99 centiseconds (rounds to 1 second)."""
        result = ephem.cs2timestr(99)
        # 99 + 50 = 149, 149 // 100 = 1 second
        assert ":00:01" in result

    def test_cs2timestr_100_centiseconds(self):
        """Test conversion of 100 centiseconds (exactly 1 second)."""
        result = ephem.cs2timestr(100)
        assert ":00:01" in result

    def test_cs2timestr_5950_centiseconds(self):
        """Test conversion near minute boundary (rounds up to 1 minute)."""
        # 5950 cs = 59.50 seconds, rounds to 60 seconds = 1 minute
        result = ephem.cs2timestr(5950)
        assert ":01:00" in result

    def test_cs2timestr_5949_centiseconds(self):
        """Test conversion near minute boundary (does not round up)."""
        # 5949 cs = 59.49 seconds, rounds to 59 seconds
        result = ephem.cs2timestr(5949)
        assert ":00:59" in result

    def test_cs2timestr_6000_centiseconds(self):
        """Test conversion of 6000 centiseconds (exactly 1 minute)."""
        result = ephem.cs2timestr(6000)
        assert ":01:00" in result


class TestCs2timestrFormat:
    """Tests for output format consistency."""

    def test_cs2timestr_format_has_two_colons(self):
        """Test that the format has two colons (HH:MM:SS)."""
        result = ephem.cs2timestr(12 * CS_HOUR + 34 * CS_MINUTE + 56 * CS_SECOND)
        assert result.count(":") == 2

    def test_cs2timestr_minutes_have_leading_zero(self):
        """Test that minutes always have leading zeros when < 10."""
        result = ephem.cs2timestr(CS_HOUR + 5 * CS_MINUTE)  # 1 hour 5 minutes
        assert ":05:" in result

    def test_cs2timestr_seconds_have_leading_zero(self):
        """Test that seconds always have leading zeros when < 10."""
        result = ephem.cs2timestr(CS_HOUR + 5 * CS_SECOND)  # 1 hour 0 min 5 sec
        assert ":00:05" in result


class TestCs2timestrRounding:
    """Tests for rounding behavior."""

    def test_rounding_at_second_boundary(self):
        """Test that centiseconds are rounded to nearest second."""
        # 49 centiseconds should round to 0 seconds
        result = ephem.cs2timestr(49)
        assert ":00:00" in result

        # 50 centiseconds should round to 1 second
        result = ephem.cs2timestr(50)
        assert ":00:01" in result

    def test_rounding_cascades_to_minute(self):
        """Test that rounding can cascade from seconds to minutes."""
        # 59.5 seconds = 5950 cs, should become 1:00
        result = ephem.cs2timestr(5950)
        assert ":01:00" in result

    def test_rounding_cascades_to_hour(self):
        """Test that rounding can cascade from minutes to hours."""
        # 59 minutes 59.5 seconds = 59*6000 + 5950 = 354000 + 5950 = 359950 cs
        # Should become 1:00:00
        result = ephem.cs2timestr(359950)
        assert "1:00:00" in result
