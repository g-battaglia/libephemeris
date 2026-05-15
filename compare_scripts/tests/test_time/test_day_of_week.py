"""
Tests for day_of_week function.

Tests cover:
- Known calendar dates and their corresponding day of week
- The formula: floor(jd + 1.5) % 7 returns 0=Monday, 1=Tuesday, ..., 6=Sunday
- Edge cases around midnight and fractional days
- Comparison with pyswisseph results
"""

import pytest
import libephemeris as ephem


class TestDayOfWeekKnownDates:
    """Test day_of_week with known calendar dates."""

    @pytest.mark.unit
    def test_j2000_epoch_saturday(self):
        """J2000.0 epoch: Jan 1, 2000 was a Saturday (5)."""
        jd = ephem.julday(2000, 1, 1, 12.0)
        assert ephem.day_of_week(jd) == 5  # Saturday

    @pytest.mark.unit
    def test_monday(self):
        """Test a known Monday."""
        # Jan 3, 2000 was a Monday
        jd = ephem.julday(2000, 1, 3, 12.0)
        assert ephem.day_of_week(jd) == 0  # Monday

    @pytest.mark.unit
    def test_tuesday(self):
        """Test a known Tuesday."""
        # Jan 4, 2000 was a Tuesday
        jd = ephem.julday(2000, 1, 4, 12.0)
        assert ephem.day_of_week(jd) == 1  # Tuesday

    @pytest.mark.unit
    def test_wednesday(self):
        """Test a known Wednesday."""
        # Jan 5, 2000 was a Wednesday
        jd = ephem.julday(2000, 1, 5, 12.0)
        assert ephem.day_of_week(jd) == 2  # Wednesday

    @pytest.mark.unit
    def test_thursday(self):
        """Test a known Thursday."""
        # Jan 6, 2000 was a Thursday
        jd = ephem.julday(2000, 1, 6, 12.0)
        assert ephem.day_of_week(jd) == 3  # Thursday

    @pytest.mark.unit
    def test_friday(self):
        """Test a known Friday."""
        # Jan 7, 2000 was a Friday
        jd = ephem.julday(2000, 1, 7, 12.0)
        assert ephem.day_of_week(jd) == 4  # Friday

    @pytest.mark.unit
    def test_sunday(self):
        """Test a known Sunday."""
        # Jan 2, 2000 was a Sunday
        jd = ephem.julday(2000, 1, 2, 12.0)
        assert ephem.day_of_week(jd) == 6  # Sunday


class TestDayOfWeekFormula:
    """Test the mathematical formula floor(jd + 1.5) % 7."""

    @pytest.mark.unit
    def test_formula_directly(self):
        """Verify the formula matches expected output."""
        import math

        jd = ephem.julday(2000, 1, 1, 12.0)
        expected = int(math.floor(jd + 0.5)) % 7
        assert ephem.day_of_week(jd) == expected

    @pytest.mark.unit
    def test_consecutive_days_increment(self):
        """Verify consecutive days increment by 1 (mod 7)."""
        jd_start = ephem.julday(2000, 1, 1, 12.0)
        for i in range(7):
            jd = jd_start + i
            expected_dow = (ephem.day_of_week(jd_start) + i) % 7
            assert ephem.day_of_week(jd) == expected_dow

    @pytest.mark.unit
    def test_week_cycle(self):
        """Verify a full week cycle returns to the same day."""
        jd = ephem.julday(2000, 1, 1, 12.0)
        dow_start = ephem.day_of_week(jd)
        dow_after_week = ephem.day_of_week(jd + 7)
        assert dow_start == dow_after_week


class TestDayOfWeekEdgeCases:
    """Test edge cases for day_of_week."""

    @pytest.mark.unit
    def test_midnight_beginning(self):
        """Test at midnight (0 hours)."""
        jd = ephem.julday(2000, 1, 1, 0.0)
        # At midnight the day is still Jan 1, 2000 = Saturday
        assert ephem.day_of_week(jd) == 5

    @pytest.mark.unit
    def test_near_midnight_end(self):
        """Test just before midnight (23.999 hours)."""
        jd = ephem.julday(2000, 1, 1, 23.999)
        # Still Jan 1, 2000 = Saturday
        assert ephem.day_of_week(jd) == 5

    @pytest.mark.unit
    def test_noon(self):
        """Test at noon (12.0 hours)."""
        jd = ephem.julday(2000, 1, 1, 12.0)
        assert ephem.day_of_week(jd) == 5  # Saturday

    @pytest.mark.unit
    def test_ancient_date(self):
        """Test with an ancient date."""
        # Julius Caesar's death: March 15, 44 BC (Julian calendar)
        jd = ephem.julday(-43, 3, 15, 12.0, ephem.JUL_CAL)
        dow = ephem.day_of_week(jd)
        # Just verify it's a valid day (0-6)
        assert 0 <= dow <= 6

    @pytest.mark.unit
    def test_far_future_date(self):
        """Test with a far future date."""
        jd = ephem.julday(3000, 6, 15, 12.0)
        dow = ephem.day_of_week(jd)
        # Just verify it's a valid day (0-6)
        assert 0 <= dow <= 6


class TestDayOfWeekComparisonWithPyswisseph:
    """Compare day_of_week results with pyswisseph."""

    @pytest.mark.unit
    def test_matches_pyswisseph_j2000(self):
        """Verify result matches pyswisseph for J2000.0 epoch."""
        try:
            import swisseph as swe

            jd = ephem.julday(2000, 1, 1, 12.0)
            pyswisseph_result = swe.day_of_week(jd)
            assert ephem.day_of_week(jd) == pyswisseph_result
        except ImportError:
            pytest.skip("pyswisseph not installed")

    @pytest.mark.unit
    def test_matches_pyswisseph_various_dates(self):
        """Verify result matches pyswisseph for various dates."""
        try:
            import swisseph as swe

            test_dates = [
                (1985, 7, 23, 8.5),
                (2020, 12, 25, 0.0),
                (1969, 7, 20, 20.0),  # Moon landing
                (2024, 2, 29, 12.0),  # Leap year
            ]
            for year, month, day, hour in test_dates:
                jd = ephem.julday(year, month, day, hour)
                pyswisseph_result = swe.day_of_week(jd)
                assert ephem.day_of_week(jd) == pyswisseph_result, (
                    f"Mismatch for {year}-{month}-{day}"
                )
        except ImportError:
            pytest.skip("pyswisseph not installed")


class TestDayOfWeekHistoricalDates:
    """Test with historically significant dates."""

    @pytest.mark.unit
    def test_moon_landing(self):
        """Apollo 11 Moon landing: July 20, 1969 was a Sunday."""
        jd = ephem.julday(1969, 7, 20, 20.0)
        assert ephem.day_of_week(jd) == 6  # Sunday

    @pytest.mark.unit
    def test_gregorian_reform(self):
        """Gregorian reform: Oct 15, 1582 was a Friday."""
        jd = ephem.julday(1582, 10, 15, 12.0)
        assert ephem.day_of_week(jd) == 4  # Friday
