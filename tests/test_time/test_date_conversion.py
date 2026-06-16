"""
Tests for date_conversion / date_conversion function.

After Phase 2, `date_conversion` is an alias for `date_conversion` and
returns ``(valid, jd, (year, month, day, hour))`` matching the reference ephemeris.

The old calendar-conversion helper is available as
``libephemeris.time_utils.date_conversion`` (renamed import in __init__).

Tests cover:
- reference-API-compatible return format
- Validation of dates
- Calendar parameter handling
- The old calendar-conversion helper (via _date_conversion_calendars)
"""

import pytest
import libephemeris as ephem
from libephemeris.time_utils import date_conversion as calendar_convert


class TestSWEDateConversion:
    """Test date_conversion / date_conversion reference-API-compatible API."""

    @pytest.mark.unit
    def test_returns_valid_jd_tuple(self):
        """date_conversion should return (valid, jd, (y, m, d, h))."""
        result = ephem.date_conversion(2000, 1, 1, 12.0, "g")
        assert isinstance(result, tuple)
        assert len(result) == 3
        valid, jd, date_tuple = result
        assert isinstance(valid, bool)
        assert isinstance(jd, float)
        assert isinstance(date_tuple, tuple)
        assert len(date_tuple) == 4

    @pytest.mark.unit
    def test_valid_gregorian_date(self):
        """A valid Gregorian date should return valid=True."""
        valid, jd, (y, m, d, h) = ephem.date_conversion(2000, 1, 1, 12.0, "g")
        assert valid is True
        assert y == 2000
        assert m == 1
        assert d == 1
        assert h == pytest.approx(12.0, abs=1e-10)
        assert jd == pytest.approx(2451545.0, abs=0.5)

    @pytest.mark.unit
    def test_valid_julian_date(self):
        """A valid Julian date should return valid=True."""
        valid, jd, (y, m, d, h) = ephem.date_conversion(2000, 1, 1, 12.0, "j")
        assert valid is True
        assert y == 2000
        assert m == 1
        assert d == 1

    @pytest.mark.unit
    def test_bytes_calendar(self):
        """Calendar parameter should accept bytes (reference ephemeris style)."""
        valid, jd, _ = ephem.date_conversion(2000, 1, 1, 12.0, b"g")
        assert valid is True
        assert isinstance(jd, float)

    @pytest.mark.unit
    def test_defaults(self):
        """hour defaults to 12.0, calendar defaults to b'g'."""
        valid, jd, (y, m, d, h) = ephem.date_conversion(2000, 1, 1)
        assert valid is True
        assert h == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_invalid_calendar_parameter(self):
        """Invalid calendar parameter should raise ValueError."""
        with pytest.raises(ValueError, match="calendar must be 'j' or 'g'"):
            ephem.date_conversion(2000, 1, 1, 12.0, "x")

    @pytest.mark.unit
    def test_case_insensitive_calendar(self):
        """Calendar parameter should be case-insensitive."""
        result_lower = ephem.date_conversion(2000, 1, 1, 12.0, "g")
        result_upper = ephem.date_conversion(2000, 1, 1, 12.0, "G")
        assert result_lower[1] == pytest.approx(result_upper[1])

    @pytest.mark.unit
    def test_swe_date_conversion_alias(self):
        """date_conversion should be the same function as date_conversion."""
        result_bare = ephem.date_conversion(2000, 6, 15, 8.5, "g")
        result_swe = ephem.date_conversion(2000, 6, 15, 8.5, "g")
        assert result_bare == result_swe


class TestCalendarConversionHelper:
    """Test the old calendar-conversion helper (time_utils.date_conversion)."""

    @pytest.mark.unit
    def test_gregorian_reform_date(self):
        """Oct 15, 1582 Gregorian = Oct 5, 1582 Julian."""
        year, month, day, hour = calendar_convert(1582, 10, 15, 12.0, "j")
        assert year == 1582
        assert month == 10
        assert day == 5
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_modern_date(self):
        """2000-01-01 Gregorian to Julian."""
        year, month, day, hour = calendar_convert(2000, 1, 1, 12.0, "j")
        # In 2000, Julian calendar is 13 days behind Gregorian
        assert year == 1999
        assert month == 12
        assert day == 19
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_julian_reform_date(self):
        """Oct 5, 1582 Julian = Oct 15, 1582 Gregorian."""
        year, month, day, hour = calendar_convert(1582, 10, 5, 12.0, "g")
        assert year == 1582
        assert month == 10
        assert day == 15
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_gregorian_to_gregorian(self):
        """Modern date stays same when converting Gregorian to Gregorian."""
        year, month, day, hour = calendar_convert(2020, 6, 15, 8.5, "g")
        assert year == 2020
        assert month == 6
        assert day == 15
        assert hour == pytest.approx(8.5, abs=1e-10)

    @pytest.mark.unit
    def test_julian_to_julian(self):
        """Pre-reform date stays same when converting Julian to Julian."""
        year, month, day, hour = calendar_convert(1400, 3, 21, 18.0, "j")
        assert year == 1400
        assert month == 3
        assert day == 21
        assert hour == pytest.approx(18.0, abs=1e-10)

    @pytest.mark.unit
    def test_roundtrip_pre_reform_date(self):
        """Pre-reform roundtrip via JD: verify Julian→Gregorian→Julian preserves JD.

        Note: date_conversion auto-detects pre-1582 dates as Julian, so the
        direct roundtrip calendar_convert(greg, 'j') treats the Gregorian
        output as Julian (since it's before 1582). We verify correctness via
        JD equivalence instead.
        """
        from libephemeris.time_utils import julday
        from libephemeris.constants import JUL_CAL, GREG_CAL

        original = (1400, 3, 15, 12.0)
        greg = calendar_convert(*original, "g")

        # Verify the JD is the same in both calendars
        jd_julian = julday(*original, JUL_CAL)
        jd_greg = julday(*greg, GREG_CAL)
        assert jd_julian == pytest.approx(jd_greg, abs=1e-10)

    @pytest.mark.unit
    def test_invalid_calendar_parameter(self):
        """Invalid calendar parameter should raise ValueError."""
        with pytest.raises(ValueError, match="calendar must be 'j' or 'g'"):
            calendar_convert(2000, 1, 1, 12.0, "x")

    @pytest.mark.unit
    def test_newton_birth(self):
        """Gregorian Dec 25, 1642 to Julian Dec 15, 1642."""
        year, month, day, hour = calendar_convert(1642, 12, 25, 12.0, "j")
        assert year == 1642
        assert month == 12
        assert day == 15

    @pytest.mark.unit
    def test_october_revolution(self):
        """Gregorian Nov 7, 1917 to Julian Oct 25, 1917."""
        year, month, day, hour = calendar_convert(1917, 11, 7, 0.0, "j")
        assert year == 1917
        assert month == 10
        assert day == 25
