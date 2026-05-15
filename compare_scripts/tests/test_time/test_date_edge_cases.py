"""
Tests for date edge cases in Julian Day conversion.

Tests cover:
- Julian/Gregorian calendar transition (October 4-15, 1582)
- Year 0 handling (which doesn't exist in historical calendar)
- Negative years (BCE dates)
- Very ancient dates (-3000)
- Very future dates (+3000)

All tests verify compatibility with pyswisseph.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import GREG_CAL, JUL_CAL


class TestGregorianJulianTransition:
    """Test the Gregorian calendar reform transition (October 1582).

    The Gregorian calendar was introduced on October 15, 1582, following
    October 4, 1582 in the Julian calendar. The 10 days from October 5-14
    were skipped in the Gregorian calendar.
    """

    @pytest.mark.comparison
    def test_last_julian_day_oct_4_1582(self):
        """October 4, 1582 was the last day of the Julian calendar."""
        jd_lib = ephem.julday(1582, 10, 4, 12.0, GREG_CAL)
        jd_swe = swe.julday(1582, 10, 4, 12.0, GREG_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    def test_first_gregorian_day_oct_15_1582(self):
        """October 15, 1582 was the first day of the Gregorian calendar."""
        jd_lib = ephem.julday(1582, 10, 15, 12.0, GREG_CAL)
        jd_swe = swe.julday(1582, 10, 15, 12.0, GREG_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "day",
        [5, 6, 7, 8, 9, 10, 11, 12, 13, 14],
    )
    def test_skipped_days_oct_5_14_1582_gregorian(self, day):
        """October 5-14, 1582 don't exist in Gregorian but julday() still computes them.

        Swiss Ephemeris doesn't validate dates - it just computes the JD mathematically.
        These dates are technically invalid in the Gregorian calendar but the function
        should still return consistent results matching pyswisseph.
        """
        jd_lib = ephem.julday(1582, 10, day, 12.0, GREG_CAL)
        jd_swe = swe.julday(1582, 10, day, 12.0, GREG_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "day",
        [5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15],
    )
    def test_julian_calendar_oct_1582(self, day):
        """October 5-15, 1582 exist normally in the Julian calendar."""
        jd_lib = ephem.julday(1582, 10, day, 12.0, JUL_CAL)
        jd_swe = swe.julday(1582, 10, day, 12.0, JUL_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    def test_transition_day_difference(self):
        """Oct 15 Gregorian should be Oct 5 Julian (10 day difference)."""
        jd_greg = ephem.julday(1582, 10, 15, 12.0, GREG_CAL)
        jd_jul = ephem.julday(1582, 10, 5, 12.0, JUL_CAL)
        # Same astronomical moment
        assert jd_greg == pytest.approx(jd_jul, abs=1e-10)

    @pytest.mark.comparison
    def test_revjul_around_transition_gregorian(self):
        """Test revjul around the transition using Gregorian calendar."""
        # JD 2299160.5 = October 15, 1582 00:00 Gregorian
        jd = 2299160.5
        result_lib = ephem.revjul(jd, GREG_CAL)
        result_swe = swe.revjul(jd, GREG_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)

    @pytest.mark.comparison
    def test_revjul_around_transition_julian(self):
        """Test revjul around the transition using Julian calendar."""
        # JD 2299160.5 = October 5, 1582 00:00 Julian
        jd = 2299160.5
        result_lib = ephem.revjul(jd, JUL_CAL)
        result_swe = swe.revjul(jd, JUL_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)


class TestYearZero:
    """Test year 0 handling.

    Year 0 doesn't exist in the historical/proleptic Julian calendar
    (1 BCE is followed by 1 CE), but astronomical year numbering
    uses year 0 = 1 BCE. Swiss Ephemeris uses astronomical year numbering.
    """

    @pytest.mark.comparison
    def test_year_zero_jan_1(self):
        """Year 0 (1 BCE) January 1 should match pyswisseph."""
        jd_lib = ephem.julday(0, 1, 1, 12.0, JUL_CAL)
        jd_swe = swe.julday(0, 1, 1, 12.0, JUL_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    def test_year_zero_mid_year(self):
        """Year 0 mid-year should match pyswisseph."""
        jd_lib = ephem.julday(0, 6, 15, 12.0, JUL_CAL)
        jd_swe = swe.julday(0, 6, 15, 12.0, JUL_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    def test_year_zero_dec_31(self):
        """Year 0 December 31 should match pyswisseph."""
        jd_lib = ephem.julday(0, 12, 31, 12.0, JUL_CAL)
        jd_swe = swe.julday(0, 12, 31, 12.0, JUL_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    def test_transition_year_minus_1_to_year_0(self):
        """Transition from year -1 to year 0 should be 1 year."""
        jd_minus1 = ephem.julday(-1, 1, 1, 12.0, JUL_CAL)
        jd_zero = ephem.julday(0, 1, 1, 12.0, JUL_CAL)
        # Should be approximately 365 or 366 days
        days_diff = jd_zero - jd_minus1
        assert 365 <= days_diff <= 366

    @pytest.mark.comparison
    def test_transition_year_0_to_year_1(self):
        """Transition from year 0 to year 1 should be 1 year."""
        jd_zero = ephem.julday(0, 1, 1, 12.0, JUL_CAL)
        jd_one = ephem.julday(1, 1, 1, 12.0, JUL_CAL)
        days_diff = jd_one - jd_zero
        assert 365 <= days_diff <= 366

    @pytest.mark.comparison
    def test_revjul_year_zero(self):
        """Revjul should return year 0 for appropriate JD."""
        # First, get JD for year 0
        jd = ephem.julday(0, 6, 15, 12.0, JUL_CAL)
        result_lib = ephem.revjul(jd, JUL_CAL)
        result_swe = swe.revjul(jd, JUL_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)


class TestNegativeYearsBCE:
    """Test negative years (BCE dates).

    In astronomical year numbering:
    - Year -1 = 2 BCE
    - Year -2 = 3 BCE
    - etc.
    """

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (-1, 1, 1, 12.0),  # 2 BCE
            (-1, 6, 15, 0.0),
            (-1, 12, 31, 23.5),
            (-10, 3, 21, 12.0),  # 11 BCE, vernal equinox date
            (-44, 3, 15, 12.0),  # 45 BCE, around Caesar's death
            (-100, 1, 1, 12.0),  # 101 BCE
            (-500, 7, 4, 12.0),  # 501 BCE
            (-753, 4, 21, 12.0),  # 754 BCE, traditional Rome founding
            (-1000, 1, 1, 12.0),  # 1001 BCE
        ],
    )
    def test_negative_years_match_pyswisseph_julian(self, year, month, day, hour):
        """Negative years (BCE) should match pyswisseph using Julian calendar."""
        jd_lib = ephem.julday(year, month, day, hour, JUL_CAL)
        jd_swe = swe.julday(year, month, day, hour, JUL_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year",
        [-1, -10, -100, -500, -1000],
    )
    def test_negative_years_revjul_roundtrip(self, year):
        """Negative years should roundtrip correctly."""
        jd = ephem.julday(year, 6, 15, 12.0, JUL_CAL)
        result_lib = ephem.revjul(jd, JUL_CAL)
        result_swe = swe.revjul(jd, JUL_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)

    @pytest.mark.comparison
    def test_chronological_order_bce(self):
        """Earlier BCE dates should have smaller JD values."""
        jd_500bce = ephem.julday(-499, 1, 1, 12.0, JUL_CAL)  # 500 BCE
        jd_100bce = ephem.julday(-99, 1, 1, 12.0, JUL_CAL)  # 100 BCE
        jd_1bce = ephem.julday(0, 1, 1, 12.0, JUL_CAL)  # 1 BCE = year 0
        assert jd_500bce < jd_100bce < jd_1bce


class TestVeryAncientDates:
    """Test very ancient dates (around -3000, i.e., ~3001 BCE)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (-3000, 1, 1, 12.0),  # 3001 BCE
            (-3000, 6, 21, 12.0),  # Summer solstice
            (-3000, 12, 21, 12.0),  # Winter solstice
            (-2500, 3, 21, 12.0),  # 2501 BCE, pyramid era
            (-4000, 1, 1, 12.0),  # 4001 BCE
            (-4712, 1, 1, 12.0),  # Near JD epoch start
            (-4713, 11, 24, 12.0),  # JD 0 is noon on this date (Julian)
        ],
    )
    def test_ancient_dates_match_pyswisseph(self, year, month, day, hour):
        """Very ancient dates should match pyswisseph."""
        jd_lib = ephem.julday(year, month, day, hour, JUL_CAL)
        jd_swe = swe.julday(year, month, day, hour, JUL_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    def test_jd_epoch_origin(self):
        """JD 0.0 should correspond to noon on Nov 24, 4713 BCE (Julian)."""
        # -4712 in astronomical year numbering = 4713 BCE
        jd_lib = ephem.julday(-4712, 1, 1, 12.0, JUL_CAL)
        jd_swe = swe.julday(-4712, 1, 1, 12.0, JUL_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)
        # JD should be small positive number near the epoch
        assert jd_lib < 100

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year",
        [-3000, -4000, -4712],
    )
    def test_ancient_dates_revjul_roundtrip(self, year):
        """Ancient dates should roundtrip correctly."""
        jd = ephem.julday(year, 6, 15, 12.0, JUL_CAL)
        result_lib = ephem.revjul(jd, JUL_CAL)
        result_swe = swe.revjul(jd, JUL_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)

    @pytest.mark.comparison
    def test_ancient_consecutive_days(self):
        """Consecutive days in ancient times should differ by 1."""
        jd1 = ephem.julday(-3000, 6, 15, 12.0, JUL_CAL)
        jd2 = ephem.julday(-3000, 6, 16, 12.0, JUL_CAL)
        assert jd2 - jd1 == pytest.approx(1.0, abs=1e-10)


class TestVeryFutureDates:
    """Test very future dates (around +3000)."""

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (3000, 1, 1, 12.0),
            (3000, 6, 21, 12.0),  # Summer solstice
            (3000, 12, 21, 12.0),  # Winter solstice
            (2500, 7, 4, 12.0),
            (4000, 1, 1, 12.0),
            (5000, 1, 1, 12.0),
            (10000, 1, 1, 12.0),
        ],
    )
    def test_future_dates_match_pyswisseph(self, year, month, day, hour):
        """Very future dates should match pyswisseph."""
        jd_lib = ephem.julday(year, month, day, hour, GREG_CAL)
        jd_swe = swe.julday(year, month, day, hour, GREG_CAL)
        assert jd_lib == pytest.approx(jd_swe, abs=1e-10)

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "year",
        [3000, 4000, 5000, 10000],
    )
    def test_future_dates_revjul_roundtrip(self, year):
        """Future dates should roundtrip correctly."""
        jd = ephem.julday(year, 6, 15, 12.0, GREG_CAL)
        result_lib = ephem.revjul(jd, GREG_CAL)
        result_swe = swe.revjul(jd, GREG_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)

    @pytest.mark.comparison
    def test_future_consecutive_days(self):
        """Consecutive days in the future should differ by 1."""
        jd1 = ephem.julday(3000, 6, 15, 12.0, GREG_CAL)
        jd2 = ephem.julday(3000, 6, 16, 12.0, GREG_CAL)
        assert jd2 - jd1 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.comparison
    def test_future_leap_year_3000(self):
        """Year 3000 is NOT a leap year (divisible by 100 but not 400)."""
        jd_feb28 = ephem.julday(3000, 2, 28, 12.0, GREG_CAL)
        jd_mar1 = ephem.julday(3000, 3, 1, 12.0, GREG_CAL)
        # Should be 1 day apart (no Feb 29)
        assert jd_mar1 - jd_feb28 == pytest.approx(1.0, abs=1e-10)

    @pytest.mark.comparison
    def test_future_leap_year_3200(self):
        """Year 3200 IS a leap year (divisible by 400)."""
        jd_feb28 = ephem.julday(3200, 2, 28, 12.0, GREG_CAL)
        jd_feb29 = ephem.julday(3200, 2, 29, 12.0, GREG_CAL)
        jd_mar1 = ephem.julday(3200, 3, 1, 12.0, GREG_CAL)
        assert jd_feb29 - jd_feb28 == pytest.approx(1.0, abs=1e-10)
        assert jd_mar1 - jd_feb29 == pytest.approx(1.0, abs=1e-10)


class TestExtremeRanges:
    """Test extreme date ranges to verify numerical stability."""

    @pytest.mark.comparison
    def test_full_range_chronological_order(self):
        """Dates across the full range should be in chronological order."""
        jd_ancient = ephem.julday(-4000, 1, 1, 12.0, JUL_CAL)
        jd_bce = ephem.julday(-500, 1, 1, 12.0, JUL_CAL)
        jd_year0 = ephem.julday(0, 1, 1, 12.0, JUL_CAL)
        jd_ce = ephem.julday(500, 1, 1, 12.0, JUL_CAL)
        jd_modern = ephem.julday(2000, 1, 1, 12.0, GREG_CAL)
        jd_future = ephem.julday(4000, 1, 1, 12.0, GREG_CAL)

        assert jd_ancient < jd_bce < jd_year0 < jd_ce < jd_modern < jd_future

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "jd",
        [
            0.0,  # JD epoch
            0.5,  # Half day after epoch
            1.0,  # One day after epoch
            1000000.0,  # Large round number
            2000000.0,
            3000000.0,
        ],
    )
    def test_round_jd_values_match_pyswisseph(self, jd):
        """Round JD values should match pyswisseph for revjul."""
        result_lib = ephem.revjul(jd, JUL_CAL)
        result_swe = swe.revjul(jd, JUL_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)

    @pytest.mark.comparison
    def test_jd_zero_revjul(self):
        """JD 0.0 should return the correct date (Julian calendar)."""
        result_lib = ephem.revjul(0.0, JUL_CAL)
        result_swe = swe.revjul(0.0, JUL_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)

    @pytest.mark.comparison
    def test_negative_jd_values(self):
        """Negative JD values (before the epoch) should match pyswisseph."""
        jd = -1000.0
        result_lib = ephem.revjul(jd, JUL_CAL)
        result_swe = swe.revjul(jd, JUL_CAL)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)
