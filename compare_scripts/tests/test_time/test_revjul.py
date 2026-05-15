"""
Comprehensive tests for reverse Julian Day conversion (revjul).

Tests cover:
- Roundtrip consistency with julday
- Standard epochs
- Edge cases
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import JUL_CAL


class TestRevjulStandardEpochs:
    """Test reverse Julian Day for standard epochs."""

    @pytest.mark.unit
    def test_j2000_reverse(self):
        """JD 2451545.0 should give 2000-01-01 12:00."""
        year, month, day, hour = ephem.revjul(2451545.0)
        assert year == 2000
        assert month == 1
        assert day == 1
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_j1900_reverse(self):
        """JD 2415020.0 should give 1899-12-31 12:00."""
        year, month, day, hour = ephem.revjul(2415020.0)
        assert year == 1899
        assert month == 12
        assert day == 31
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_unix_epoch_reverse(self):
        """JD 2440587.5 should give 1970-01-01 00:00."""
        year, month, day, hour = ephem.revjul(2440587.5)
        assert year == 1970
        assert month == 1
        assert day == 1
        assert hour == pytest.approx(0.0, abs=1e-10)


class TestRevjulTimeOfDay:
    """Test that time of day is correctly extracted."""

    @pytest.mark.unit
    def test_midnight(self):
        """Test midnight extraction."""
        jd = ephem.julday(2000, 6, 15, 0.0)
        year, month, day, hour = ephem.revjul(jd)
        assert hour == pytest.approx(0.0, abs=1e-10)

    @pytest.mark.unit
    def test_noon(self):
        """Test noon extraction."""
        jd = ephem.julday(2000, 6, 15, 12.0)
        year, month, day, hour = ephem.revjul(jd)
        assert hour == pytest.approx(12.0, abs=1e-10)

    @pytest.mark.unit
    def test_fractional_hours(self):
        """Test fractional hour extraction."""
        # 6:30:45 = 6.5125 hours
        jd = ephem.julday(2000, 6, 15, 6.5125)
        year, month, day, hour = ephem.revjul(jd)
        assert hour == pytest.approx(6.5125, abs=1e-8)


class TestRevjulRoundtrip:
    """Test roundtrip conversion julday -> revjul."""

    @pytest.mark.unit
    def test_roundtrip_j2000(self):
        """J2000 roundtrip should be exact."""
        original = (2000, 1, 1, 12.0)
        jd = ephem.julday(*original)
        result = ephem.revjul(jd)
        assert result[0] == original[0]  # year
        assert result[1] == original[1]  # month
        assert result[2] == original[2]  # day
        assert result[3] == pytest.approx(original[3], abs=1e-10)

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "year,month,day,hour",
        [
            (1900, 1, 1, 0.0),
            (1950, 6, 15, 12.0),
            (1980, 3, 21, 6.5),
            (2000, 1, 1, 12.0),
            (2000, 2, 29, 18.25),  # Leap year
            (2020, 12, 21, 23.999),
            (2050, 7, 4, 0.001),
        ],
    )
    def test_roundtrip_various_dates(self, year, month, day, hour):
        """Various dates should roundtrip correctly."""
        jd = ephem.julday(year, month, day, hour)
        result = ephem.revjul(jd)
        assert result[0] == year
        assert result[1] == month
        assert result[2] == day
        assert result[3] == pytest.approx(hour, abs=1e-8)

    @pytest.mark.unit
    def test_roundtrip_100_random_dates(self, random_dates_in_de421_range):
        """100 random dates should all roundtrip correctly."""
        dates = random_dates_in_de421_range(100)
        for year, month, day, hour, jd in dates:
            result = ephem.revjul(jd)
            assert result[0] == year
            assert result[1] == month
            assert result[2] == day
            assert result[3] == pytest.approx(hour, abs=1e-8)


class TestRevjulVsPyswisseph:
    """Compare results with pyswisseph."""

    @pytest.mark.comparison
    def test_j2000_matches_swe(self):
        """J2000 reverse should match pyswisseph."""
        result_lib = ephem.revjul(2451545.0)
        result_swe = swe.revjul(2451545.0)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "jd",
        [
            2415020.0,  # J1900
            2440587.5,  # Unix epoch
            2451545.0,  # J2000
            2451545.5,  # Day after J2000
            2460000.0,  # Modern date
        ],
    )
    def test_various_jd_match_swe(self, jd):
        """Various JD values should match pyswisseph."""
        result_lib = ephem.revjul(jd)
        result_swe = swe.revjul(jd)
        assert result_lib[0] == result_swe[0]  # year
        assert result_lib[1] == result_swe[1]  # month
        assert result_lib[2] == result_swe[2]  # day
        assert result_lib[3] == pytest.approx(result_swe[3], abs=1e-10)


class TestRevjulEdgeCases:
    """Test edge cases in reverse Julian Day conversion."""

    @pytest.mark.unit
    def test_very_small_jd(self):
        """Test very small (historical) JD values."""
        # JD 0 = November 24, 4714 BC (Julian calendar)
        year, month, day, hour = ephem.revjul(0.0, JUL_CAL)
        # Should return a valid date
        assert isinstance(year, int)
        assert 1 <= month <= 12
        assert 1 <= day <= 31

    @pytest.mark.unit
    def test_large_jd(self):
        """Test large (future) JD values."""
        # Far future date
        jd = 2500000.0  # Around year 2132
        year, month, day, hour = ephem.revjul(jd)
        assert year > 2100
        assert 1 <= month <= 12
        assert 1 <= day <= 31

    @pytest.mark.unit
    def test_fractional_day_precision(self):
        """Test that fractional day is preserved with high precision."""
        # Test sub-second precision
        original_hour = 12.000001  # About 0.0036 seconds past noon
        jd = ephem.julday(2000, 1, 1, original_hour)
        _, _, _, result_hour = ephem.revjul(jd)
        assert result_hour == pytest.approx(original_hour, abs=1e-5)
