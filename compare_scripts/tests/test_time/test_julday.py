"""
Comprehensive tests for Julian Day conversion (julday).

Tests cover:
- Standard epochs (J2000, J1900, Unix)
- Edge cases (leap years, century boundaries, Gregorian reform)
- Comparison with pyswisseph
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import GREG_CAL, JUL_CAL


class TestJuldayStandardEpochs:
    """Test Julian Day calculation for standard astronomical epochs."""

    @pytest.mark.unit
    def test_j2000_epoch(self):
        """J2000.0 = 2000-01-01 12:00 TT = JD 2451545.0"""
        jd = ephem.julday(2000, 1, 1, 12.0)
        assert abs(jd - 2451545.0) < 1e-10

    @pytest.mark.unit
    def test_j1900_epoch(self):
        """J1900.0 = 1899-12-31 12:00 TT = JD 2415020.0"""
        jd = ephem.julday(1899, 12, 31, 12.0)
        assert abs(jd - 2415020.0) < 1e-10

    @pytest.mark.unit
    def test_unix_epoch(self):
        """Unix epoch = 1970-01-01 00:00 UT = JD 2440587.5"""
        jd = ephem.julday(1970, 1, 1, 0.0)
        assert abs(jd - 2440587.5) < 1e-10

    @pytest.mark.unit
    def test_j2050_epoch(self):
        """J2050.0 for DE421 end range testing."""
        jd = ephem.julday(2050, 1, 1, 12.0)
        expected = 2469807.5  # Approximate
        assert abs(jd - expected) < 1.0  # Within 1 day


class TestJuldayTimeOfDay:
    """Test different times of day are correctly handled."""

    @pytest.mark.unit
    def test_midnight(self):
        """Midnight should give .5 fractional day."""
        jd = ephem.julday(2000, 1, 1, 0.0)
        # JD starts at noon, so midnight is 0.5 days before/after
        assert jd == pytest.approx(2451544.5, abs=1e-10)

    @pytest.mark.unit
    def test_noon(self):
        """Noon should give .0 fractional day."""
        jd = ephem.julday(2000, 1, 1, 12.0)
        # Check it ends in .0
        assert jd == pytest.approx(2451545.0, abs=1e-10)

    @pytest.mark.unit
    def test_6am(self):
        """6 AM = 0.25 days from midnight."""
        jd = ephem.julday(2000, 1, 1, 6.0)
        assert jd == pytest.approx(2451544.75, abs=1e-10)

    @pytest.mark.unit
    def test_6pm(self):
        """6 PM = 0.75 days from midnight."""
        jd = ephem.julday(2000, 1, 1, 18.0)
        assert jd == pytest.approx(2451545.25, abs=1e-10)

    @pytest.mark.unit
    def test_fractional_hours(self):
        """Test fractional hour conversion."""
        # 6:30:45 = 6 + 30/60 + 45/3600 = 6.5125 hours
        jd = ephem.julday(2000, 1, 1, 6.5125)
        expected_fraction = 6.5125 / 24
        jd_midnight = ephem.julday(2000, 1, 1, 0.0)
        assert jd - jd_midnight == pytest.approx(expected_fraction, abs=1e-6)


class TestJuldayLeapYears:
    """Test leap year handling in Julian Day calculation."""

    @pytest.mark.unit
    def test_leap_year_2000(self):
        """2000 is a leap year (divisible by 400)."""
        jd_feb28 = ephem.julday(2000, 2, 28, 12.0)
        jd_feb29 = ephem.julday(2000, 2, 29, 12.0)
        jd_mar1 = ephem.julday(2000, 3, 1, 12.0)
        assert jd_feb29 - jd_feb28 == pytest.approx(1.0, abs=1e-10)
        assert jd_mar1 - jd_feb29 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.unit
    def test_leap_year_2004(self):
        """2004 is a regular leap year."""
        jd_feb28 = ephem.julday(2004, 2, 28, 12.0)
        jd_feb29 = ephem.julday(2004, 2, 29, 12.0)
        assert jd_feb29 - jd_feb28 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.unit
    def test_non_leap_year_1900(self):
        """1900 is NOT a leap year (divisible by 100 but not 400)."""
        jd_feb28 = ephem.julday(1900, 2, 28, 12.0)
        jd_mar1 = ephem.julday(1900, 3, 1, 12.0)
        # Should be 1 day apart (no Feb 29)
        assert jd_mar1 - jd_feb28 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.unit
    def test_non_leap_year_2100(self):
        """2100 is NOT a leap year."""
        jd_feb28 = ephem.julday(2100, 2, 28, 12.0)
        jd_mar1 = ephem.julday(2100, 3, 1, 12.0)
        assert jd_mar1 - jd_feb28 == pytest.approx(1.0, abs=1e-10)


class TestJuldayMonthBoundaries:
    """Test month and year boundary handling."""

    @pytest.mark.unit
    def test_month_end_january(self):
        """January has 31 days."""
        jd_jan31 = ephem.julday(2000, 1, 31, 12.0)
        jd_feb1 = ephem.julday(2000, 2, 1, 12.0)
        assert jd_feb1 - jd_jan31 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.unit
    def test_month_end_april(self):
        """April has 30 days."""
        jd_apr30 = ephem.julday(2000, 4, 30, 12.0)
        jd_may1 = ephem.julday(2000, 5, 1, 12.0)
        assert jd_may1 - jd_apr30 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.unit
    def test_year_boundary(self):
        """Test year rollover."""
        jd_dec31 = ephem.julday(1999, 12, 31, 12.0)
        jd_jan1 = ephem.julday(2000, 1, 1, 12.0)
        assert jd_jan1 - jd_dec31 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.unit
    def test_y2k_boundary(self):
        """Test Y2K transition."""
        jd_before = ephem.julday(1999, 12, 31, 23.99)
        jd_after = ephem.julday(2000, 1, 1, 0.01)
        # Should be about 0.02/24 = 0.0008 days apart (+ 1 day for date change)
        expected_diff = 1.0 + (0.01 - 23.99) / 24
        # Actually: 0.01/24 + (24-23.99)/24 = difference
        assert jd_after > jd_before


class TestJuldayGregorianJulian:
    """Test Gregorian vs Julian calendar handling."""

    @pytest.mark.unit
    def test_gregorian_default(self):
        """Gregorian calendar should be default."""
        jd_greg = ephem.julday(2000, 1, 1, 12.0, GREG_CAL)
        jd_default = ephem.julday(2000, 1, 1, 12.0)
        assert jd_greg == jd_default

    @pytest.mark.unit
    def test_julian_calendar_differs(self):
        """Julian calendar should give different result for modern dates."""
        jd_greg = ephem.julday(2000, 1, 1, 12.0, GREG_CAL)
        jd_jul = ephem.julday(2000, 1, 1, 12.0, JUL_CAL)
        # Julian calendar is ~13 days behind Gregorian in 2000
        assert abs(jd_greg - jd_jul) > 10


class TestJuldayVsPyswisseph:
    """Compare results with pyswisseph for exact compatibility."""

    @pytest.mark.comparison
    def test_j2000_matches_swe(self):
        """J2000 should match pyswisseph exactly."""
        jd_lib = ephem.julday(2000, 1, 1, 12.0)
        jd_swe = swe.julday(2000, 1, 1, 12.0)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (1900, 1, 1, 0.0),
            (1950, 6, 15, 12.0),
            (1980, 3, 21, 6.5),
            (2000, 1, 1, 12.0),
            (2020, 12, 21, 18.0),
            (2050, 7, 4, 0.0),
        ],
    )
    def test_various_dates_match_swe(self, year, month, day, hour):
        """Various dates should match pyswisseph."""
        jd_lib = ephem.julday(year, month, day, hour)
        jd_swe = swe.julday(year, month, day, hour)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    def test_100_random_dates(self, random_dates_in_de421_range):
        """100 random dates should all match pyswisseph."""
        dates = random_dates_in_de421_range(100)
        for year, month, day, hour, _ in dates:
            jd_lib = ephem.julday(year, month, day, hour)
            jd_swe = swe.julday(year, month, day, hour)
            assert jd_lib == pytest.approx(jd_swe, abs=1e-10), (
                f"Mismatch for {year}-{month}-{day} {hour}h"
            )


class TestJuldayConsistency:
    """Test internal consistency of Julian Day calculations."""

    @pytest.mark.unit
    def test_consecutive_days_differ_by_one(self):
        """Consecutive days at same time should differ by exactly 1."""
        for day in range(1, 28):
            jd1 = ephem.julday(2000, 1, day, 12.0)
            jd2 = ephem.julday(2000, 1, day + 1, 12.0)
            assert jd2 - jd1 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.unit
    def test_year_has_365_or_366_days(self):
        """A year should have 365 or 366 days."""
        jd_start = ephem.julday(2000, 1, 1, 0.0)
        jd_end = ephem.julday(2001, 1, 1, 0.0)
        days = jd_end - jd_start
        assert days in [365.0, 366.0]  # 2000 is leap year
        assert days == 366.0  # 2000 is leap year

    @pytest.mark.unit
    def test_chronological_order(self):
        """Later dates should have larger JD values."""
        jd1 = ephem.julday(1990, 1, 1, 0.0)
        jd2 = ephem.julday(2000, 1, 1, 0.0)
        jd3 = ephem.julday(2010, 1, 1, 0.0)
        assert jd1 < jd2 < jd3
