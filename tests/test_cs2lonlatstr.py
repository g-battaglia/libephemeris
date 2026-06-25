"""
Tests for the cs2lonlatstr centiseconds to longitude/latitude string conversion function.

Tests verify that cs2lonlatstr correctly converts angles in centiseconds to
the compact reference ephemeris format (e.g., "45N30", "122E15'30").
"""

import libephemeris as ephem


# Constants for centiseconds calculations
CS1 = 3600 * 100  # 360000 centiseconds per degree


class TestCs2lonlatstrBasic:
    """Basic functionality tests for cs2lonlatstr."""

    def test_cs2lonlatstr_exported(self):
        """Test that cs2lonlatstr is exported from the package."""
        assert hasattr(ephem, "cs2lonlatstr")
        assert callable(ephem.cs2lonlatstr)

    def test_cs2lonlatstr_zero(self):
        """Test conversion of zero."""
        result = ephem.cs2lonlatstr(0, "N", "S")
        assert result == "0N00"

    def test_cs2lonlatstr_returns_string(self):
        """Test that cs2lonlatstr returns a string."""
        result = ephem.cs2lonlatstr(CS1, "E", "W")
        assert isinstance(result, str)

    def test_cs2lonlatstr_one_degree(self):
        """Test conversion of exactly 1 degree."""
        result = ephem.cs2lonlatstr(CS1, "N", "S")
        assert result == "1N00"


class TestCs2lonlatstrDirectionalCharacters:
    """Tests for directional character handling."""

    def test_positive_uses_plus_char(self):
        """Test that positive values use plus_char."""
        result = ephem.cs2lonlatstr(CS1 * 45, "N", "S")
        assert "N" in result
        assert "S" not in result

    def test_negative_uses_minus_char(self):
        """Test that negative values use minus_char."""
        result = ephem.cs2lonlatstr(-CS1 * 45, "N", "S")
        assert "S" in result
        assert "N" not in result

    def test_zero_uses_plus_char(self):
        """Test that zero uses plus_char."""
        result = ephem.cs2lonlatstr(0, "N", "S")
        assert "N" in result
        assert "S" not in result

    def test_longitude_east_west(self):
        """Test longitude with E/W characters."""
        result_east = ephem.cs2lonlatstr(CS1 * 122, "E", "W")
        assert "E" in result_east
        assert "W" not in result_east

        result_west = ephem.cs2lonlatstr(-CS1 * 122, "E", "W")
        assert "W" in result_west
        assert "E" not in result_west

    def test_custom_direction_chars(self):
        """Test with custom directional characters."""
        result = ephem.cs2lonlatstr(CS1 * 30, "+", "-")
        assert "+" in result

        result_neg = ephem.cs2lonlatstr(-CS1 * 30, "+", "-")
        assert "-" in result_neg


class TestCs2lonlatstrValues:
    """Tests for various value conversions."""

    def test_45_degrees_30_minutes(self):
        """Test conversion of 45°30'."""
        cs_value = 45 * CS1 + 30 * 6000
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert result == "45N30"

    def test_122_degrees_15_minutes_30_seconds(self):
        """Test conversion of 122°15'30"."""
        cs_value = 122 * CS1 + 15 * 6000 + 30 * 100
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert result == "122E15'30"

    def test_negative_longitude(self):
        """Test negative longitude (West)."""
        cs_value = -(122 * CS1 + 15 * 6000 + 30 * 100)
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert result == "122W15'30"

    def test_90_degrees_latitude(self):
        """Test 90° latitude (North Pole)."""
        result = ephem.cs2lonlatstr(90 * CS1, "N", "S")
        assert result == "90N00"

    def test_negative_90_degrees_latitude(self):
        """Test -90° latitude (South Pole)."""
        result = ephem.cs2lonlatstr(-90 * CS1, "N", "S")
        assert result == "90S00"


class TestCs2lonlatstrRounding:
    """Tests for seconds rounding behavior."""

    def test_round_up_centiseconds(self):
        """Test that centiseconds >= 50 round up to next second."""
        cs_value = CS1 + 50
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert result == "1N00'01"

    def test_round_down_centiseconds(self):
        """Test that centiseconds < 50 round down."""
        cs_value = CS1 + 49
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert result == "1N00"

    def test_seconds_carry_to_minutes(self):
        """Test that 60 seconds rounds to next minute."""
        cs_value = CS1 + 59 * 6000 + 59 * 100 + 50
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert result == "2N00"

    def test_no_centiseconds_in_output(self):
        """Test that output does not contain centisecond decimals."""
        cs_value = CS1 + 45 * 100 + 78  # 1°0'45.78"
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        # Should show "46" (rounded from 45.78) not "45.78"
        assert ".78" not in result
        assert result == "1N00'46"


class TestCs2lonlatstrFormat:
    """Tests for output format consistency."""

    def test_format_direction_embedded(self):
        """Test that direction char is embedded between degree and minutes."""
        result = ephem.cs2lonlatstr(CS1 * 45, "N", "S")
        assert result == "45N00"

        result_neg = ephem.cs2lonlatstr(-CS1 * 45, "N", "S")
        assert result_neg == "45S00"

    def test_format_with_seconds(self):
        """Test format when seconds are present."""
        cs_value = 45 * CS1 + 30 * 6000 + 15 * 100
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert result == "45N30'15"


class TestCs2lonlatstrEdgeCases:
    """Edge case tests for cs2lonlatstr."""

    def test_small_positive_value(self):
        """Test very small positive value."""
        result = ephem.cs2lonlatstr(1, "N", "S")
        assert result == "0N00"

    def test_small_negative_value(self):
        """Test very small negative value."""
        result = ephem.cs2lonlatstr(-1, "N", "S")
        assert result == "0S00"

    def test_180_degrees_longitude(self):
        """Test 180° longitude (dateline)."""
        result = ephem.cs2lonlatstr(180 * CS1, "E", "W")
        assert result == "180E00"

    def test_large_value(self):
        """Test large value beyond normal coordinate range."""
        result = ephem.cs2lonlatstr(360 * CS1, "E", "W")
        assert result == "360E00"

    def test_30_minutes_only(self):
        """Test 0 degrees 30 minutes."""
        result = ephem.cs2lonlatstr(30 * 6000, "N", "S")
        assert result == "0N30"


class TestCs2lonlatstrTypicalUseCase:
    """Tests for typical geographic coordinate use cases."""

    def test_rome_latitude(self):
        """Test latitude of Rome (~41°54' N)."""
        cs_value = 41 * CS1 + 54 * 6000
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert result == "41N54"

    def test_rome_longitude(self):
        """Test longitude of Rome (~12°30' E)."""
        cs_value = 12 * CS1 + 30 * 6000
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert result == "12E30"

    def test_new_york_longitude(self):
        """Test longitude of New York (~74°00' W)."""
        cs_value = -(74 * CS1)
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert result == "74W00"

    def test_sydney_latitude(self):
        """Test latitude of Sydney (~33°52' S)."""
        cs_value = -(33 * CS1 + 52 * 6000)
        result = ephem.cs2lonlatstr(cs_value, "N", "S")
        assert result == "33S52"

    def test_sydney_longitude(self):
        """Test longitude of Sydney (~151°12' E)."""
        cs_value = 151 * CS1 + 12 * 6000
        result = ephem.cs2lonlatstr(cs_value, "E", "W")
        assert result == "151E12"
