"""
Comprehensive tests for TAI (International Atomic Time) time scale functions.

Tests cover:
- TAI-UTC offset (leap seconds) retrieval
- UTC to TAI conversion
- TAI to UTC conversion
- TT to TAI conversion (fixed 32.184s offset)
- TAI to TT conversion
- Round-trip conversions
- Historical dates and leap second boundaries
"""

import pytest
import libephemeris as ephem


class TestTAIConstants:
    """Test TAI-related constants."""

    @pytest.mark.unit
    def test_tt_tai_offset_seconds(self):
        """TT-TAI offset should be exactly 32.184 seconds."""
        assert ephem.TT_TAI_OFFSET_SECONDS == 32.184

    @pytest.mark.unit
    def test_tt_tai_offset_days(self):
        """TT-TAI offset in days should be 32.184/86400."""
        expected = 32.184 / 86400.0
        assert ephem.TT_TAI_OFFSET_DAYS == pytest.approx(expected, rel=1e-10)


class TestGetTaiUtcForJd:
    """Test the get_tai_utc_for_jd function."""

    @pytest.mark.unit
    def test_modern_date_2020(self):
        """TAI-UTC should be 37 seconds for dates after 2017."""
        jd = ephem.julday(2020, 1, 1, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        assert tai_utc == 37.0

    @pytest.mark.unit
    def test_modern_date_2024(self):
        """TAI-UTC should be 37 seconds for 2024."""
        jd = ephem.julday(2024, 6, 15, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        assert tai_utc == 37.0

    @pytest.mark.unit
    def test_before_2017_leap_second(self):
        """TAI-UTC should be 36 seconds before 2017-01-01."""
        jd = ephem.julday(2016, 12, 31, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        assert tai_utc == 36.0

    @pytest.mark.unit
    def test_2015_leap_second(self):
        """TAI-UTC should be 35 seconds between 2015-07-01 and 2017-01-01."""
        jd = ephem.julday(2015, 8, 1, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        # After 2015-07-01 leap second, before 2017-01-01
        assert tai_utc == 36.0

    @pytest.mark.unit
    def test_year_2000(self):
        """TAI-UTC should be 32 seconds around J2000."""
        jd = ephem.julday(2000, 1, 1, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        assert tai_utc == 32.0

    @pytest.mark.unit
    def test_year_1990(self):
        """TAI-UTC should be 25 seconds in 1990."""
        jd = ephem.julday(1990, 6, 1, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        assert tai_utc == 25.0

    @pytest.mark.unit
    def test_returns_float(self):
        """Function should return a float."""
        jd = ephem.julday(2020, 1, 1, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        assert isinstance(tai_utc, float)


class TestUtcToTaiJd:
    """Test UTC to TAI Julian Day conversion."""

    @pytest.mark.unit
    def test_basic_conversion(self):
        """Basic UTC to TAI conversion should work."""
        jd_tai = ephem.utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0)
        # TAI should be ahead of UTC by 37 seconds in 2020
        jd_utc = ephem.julday(2020, 1, 1, 0.0)
        diff_seconds = (jd_tai - jd_utc) * 86400
        assert diff_seconds == pytest.approx(37.0, abs=0.1)

    @pytest.mark.unit
    def test_j2000_epoch(self):
        """Test TAI for J2000 epoch."""
        jd_tai = ephem.utc_to_tai_jd(2000, 1, 1, 12, 0, 0.0)
        # TAI should be ahead of UTC by 32 seconds in 2000
        jd_utc = ephem.julday(2000, 1, 1, 12.0)
        diff_seconds = (jd_tai - jd_utc) * 86400
        assert diff_seconds == pytest.approx(32.0, abs=0.1)

    @pytest.mark.unit
    def test_returns_float(self):
        """Function should return a float."""
        jd_tai = ephem.utc_to_tai_jd(2020, 1, 1, 12, 0, 0.0)
        assert isinstance(jd_tai, float)

    @pytest.mark.unit
    def test_tai_always_greater_than_utc(self):
        """TAI should always be greater than UTC (after 1972)."""
        dates = [
            (1980, 1, 1, 0, 0, 0.0),
            (2000, 6, 15, 12, 30, 45.0),
            (2020, 12, 31, 23, 59, 59.0),
        ]
        for date in dates:
            jd_tai = ephem.utc_to_tai_jd(*date)
            jd_utc = ephem.julday(
                date[0], date[1], date[2], date[3] + date[4] / 60.0 + date[5] / 3600.0
            )
            assert jd_tai > jd_utc, f"TAI should be > UTC for {date}"

    @pytest.mark.unit
    def test_with_fractional_seconds(self):
        """Test with fractional seconds."""
        jd_tai = ephem.utc_to_tai_jd(2020, 1, 1, 12, 0, 30.5)
        assert jd_tai > 0


class TestTaiJdToUtc:
    """Test TAI Julian Day to UTC conversion."""

    @pytest.mark.unit
    def test_basic_conversion(self):
        """Basic TAI to UTC conversion should work."""
        jd_tai = ephem.utc_to_tai_jd(2020, 1, 1, 12, 0, 0.0)
        year, month, day, hour, minute, second = ephem.tai_jd_to_utc(jd_tai)
        assert year == 2020
        assert month == 1
        assert day == 1
        assert hour == 12
        assert minute == 0
        assert second == pytest.approx(0.0, abs=0.001)

    @pytest.mark.unit
    def test_round_trip(self):
        """Round-trip UTC -> TAI -> UTC should preserve values."""
        orig_date = (2020, 6, 15, 14, 30, 45.5)
        jd_tai = ephem.utc_to_tai_jd(*orig_date)
        result = ephem.tai_jd_to_utc(jd_tai)

        assert result[0] == orig_date[0]  # year
        assert result[1] == orig_date[1]  # month
        assert result[2] == orig_date[2]  # day
        assert result[3] == orig_date[3]  # hour
        assert result[4] == orig_date[4]  # minute
        assert result[5] == pytest.approx(orig_date[5], abs=0.01)  # second

    @pytest.mark.unit
    def test_returns_tuple(self):
        """Function should return a 6-tuple."""
        jd_tai = ephem.utc_to_tai_jd(2020, 1, 1, 12, 0, 0.0)
        result = ephem.tai_jd_to_utc(jd_tai)
        assert isinstance(result, tuple)
        assert len(result) == 6


class TestTtToTaiJd:
    """Test TT to TAI Julian Day conversion."""

    @pytest.mark.unit
    def test_j2000_epoch(self):
        """Test TT to TAI conversion for J2000."""
        jd_tt = 2451545.0  # J2000.0 in TT
        jd_tai = ephem.tt_to_tai_jd(jd_tt)
        # TAI = TT - 32.184 seconds
        expected_diff = 32.184 / 86400.0
        assert jd_tt - jd_tai == pytest.approx(expected_diff, rel=1e-6)

    @pytest.mark.unit
    def test_tai_is_less_than_tt(self):
        """TAI should always be less than TT."""
        jd_tt = 2451545.0
        jd_tai = ephem.tt_to_tai_jd(jd_tt)
        assert jd_tai < jd_tt

    @pytest.mark.unit
    def test_offset_is_constant(self):
        """TT-TAI offset should be constant for all dates."""
        dates_tt = [2440000.0, 2451545.0, 2460000.0]
        expected_diff = 32.184 / 86400.0
        for jd_tt in dates_tt:
            jd_tai = ephem.tt_to_tai_jd(jd_tt)
            diff = jd_tt - jd_tai
            assert diff == pytest.approx(expected_diff, rel=1e-6)

    @pytest.mark.unit
    def test_returns_float(self):
        """Function should return a float."""
        jd_tai = ephem.tt_to_tai_jd(2451545.0)
        assert isinstance(jd_tai, float)


class TestTaiToTtJd:
    """Test TAI to TT Julian Day conversion."""

    @pytest.mark.unit
    def test_basic_conversion(self):
        """Basic TAI to TT conversion should work."""
        jd_tai = 2451545.0
        jd_tt = ephem.tai_to_tt_jd(jd_tai)
        # TT = TAI + 32.184 seconds
        expected_diff = 32.184 / 86400.0
        assert jd_tt - jd_tai == pytest.approx(expected_diff, rel=1e-6)

    @pytest.mark.unit
    def test_tt_is_greater_than_tai(self):
        """TT should always be greater than TAI."""
        jd_tai = 2451545.0
        jd_tt = ephem.tai_to_tt_jd(jd_tai)
        assert jd_tt > jd_tai

    @pytest.mark.unit
    def test_round_trip_tt_tai(self):
        """Round-trip TT -> TAI -> TT should preserve value."""
        jd_tt_orig = 2451545.0
        jd_tai = ephem.tt_to_tai_jd(jd_tt_orig)
        jd_tt_back = ephem.tai_to_tt_jd(jd_tai)
        assert jd_tt_back == pytest.approx(jd_tt_orig, rel=1e-15)

    @pytest.mark.unit
    def test_returns_float(self):
        """Function should return a float."""
        jd_tt = ephem.tai_to_tt_jd(2451545.0)
        assert isinstance(jd_tt, float)


class TestTaiToUtcJd:
    """Test TAI to UTC Julian Day conversion."""

    @pytest.mark.unit
    def test_basic_conversion(self):
        """Basic TAI to UTC JD conversion."""
        jd_utc = ephem.julday(2020, 1, 1, 12.0)
        jd_tai = ephem.utc_jd_to_tai(jd_utc)
        jd_utc_back = ephem.tai_to_utc_jd(jd_tai)
        assert jd_utc_back == pytest.approx(jd_utc, abs=1e-9)

    @pytest.mark.unit
    def test_utc_less_than_tai(self):
        """UTC JD should be less than TAI JD."""
        jd_tai = 2451545.0
        jd_utc = ephem.tai_to_utc_jd(jd_tai)
        assert jd_utc < jd_tai

    @pytest.mark.unit
    def test_round_trip(self):
        """Round-trip UTC -> TAI -> UTC should preserve value."""
        jd_utc_orig = ephem.julday(2020, 6, 15, 12.0)
        jd_tai = ephem.utc_jd_to_tai(jd_utc_orig)
        jd_utc_back = ephem.tai_to_utc_jd(jd_tai)
        assert jd_utc_back == pytest.approx(jd_utc_orig, abs=1e-9)


class TestUtcJdToTai:
    """Test UTC Julian Day to TAI conversion."""

    @pytest.mark.unit
    def test_basic_conversion(self):
        """Basic UTC JD to TAI conversion."""
        jd_utc = ephem.julday(2020, 1, 1, 12.0)
        jd_tai = ephem.utc_jd_to_tai(jd_utc)
        # Should be 37 seconds ahead
        diff_seconds = (jd_tai - jd_utc) * 86400
        assert diff_seconds == pytest.approx(37.0, abs=0.1)

    @pytest.mark.unit
    def test_tai_greater_than_utc(self):
        """TAI JD should be greater than UTC JD."""
        jd_utc = ephem.julday(2020, 1, 1, 12.0)
        jd_tai = ephem.utc_jd_to_tai(jd_utc)
        assert jd_tai > jd_utc

    @pytest.mark.unit
    def test_j2000(self):
        """Test for J2000 epoch."""
        jd_utc = ephem.julday(2000, 1, 1, 12.0)
        jd_tai = ephem.utc_jd_to_tai(jd_utc)
        diff_seconds = (jd_tai - jd_utc) * 86400
        assert diff_seconds == pytest.approx(32.0, abs=0.1)


class TestTAITimeScaleRelationships:
    """Test the relationships between TAI, TT, and UTC."""

    @pytest.mark.unit
    def test_tt_equals_tai_plus_offset(self):
        """TT = TAI + 32.184 seconds exactly."""
        jd_tai = 2451545.0
        jd_tt = ephem.tai_to_tt_jd(jd_tai)
        diff_seconds = (jd_tt - jd_tai) * 86400
        assert diff_seconds == pytest.approx(32.184, rel=1e-6)

    @pytest.mark.unit
    def test_tai_equals_utc_plus_leap_seconds(self):
        """TAI = UTC + leap_seconds."""
        jd_utc = ephem.julday(2020, 1, 1, 12.0)
        jd_tai = ephem.utc_jd_to_tai(jd_utc)
        leap_seconds = ephem.get_tai_utc_for_jd(jd_utc)
        expected_tai = jd_utc + leap_seconds / 86400.0
        assert jd_tai == pytest.approx(expected_tai, rel=1e-10)

    @pytest.mark.unit
    def test_tt_minus_utc_equals_delta_t(self):
        """TT - UT1 should approximately equal Delta T."""
        # Get JD values
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 1, 1, 12, 0, 0.0)
        jd_utc = ephem.julday(2020, 1, 1, 12.0)

        # Get Delta T
        delta_t_days = ephem.deltat(jd_ut1)
        delta_t_seconds = delta_t_days * 86400

        # TT - UT1 should be Delta T
        computed_delta_t = (jd_tt - jd_ut1) * 86400
        assert computed_delta_t == pytest.approx(delta_t_seconds, abs=0.1)

    @pytest.mark.unit
    def test_full_chain_conversion(self):
        """Test UTC -> TAI -> TT chain."""
        # Start with UTC
        jd_utc = ephem.julday(2020, 1, 1, 12.0)

        # Convert to TAI
        jd_tai = ephem.utc_jd_to_tai(jd_utc)

        # Convert to TT
        jd_tt = ephem.tai_to_tt_jd(jd_tai)

        # Verify the expected relationships
        # TAI - UTC = leap seconds (37 in 2020)
        tai_utc = (jd_tai - jd_utc) * 86400
        assert tai_utc == pytest.approx(37.0, abs=0.1)

        # TT - TAI = 32.184 seconds (always)
        tt_tai = (jd_tt - jd_tai) * 86400
        assert tt_tai == pytest.approx(32.184, rel=1e-6)

        # TT - UTC = TAI-UTC + 32.184
        tt_utc = (jd_tt - jd_utc) * 86400
        expected_tt_utc = 37.0 + 32.184
        assert tt_utc == pytest.approx(expected_tt_utc, abs=0.1)


class TestTAILeapSecondBoundaries:
    """Test behavior around leap second boundaries."""

    @pytest.mark.unit
    def test_before_and_after_2017_leap_second(self):
        """TAI-UTC should change across the 2017-01-01 leap second."""
        # Just before leap second
        jd_before = ephem.julday(2016, 12, 31, 23.0)
        tai_utc_before = ephem.get_tai_utc_for_jd(jd_before)

        # Just after leap second
        jd_after = ephem.julday(2017, 1, 1, 1.0)
        tai_utc_after = ephem.get_tai_utc_for_jd(jd_after)

        assert tai_utc_before == 36.0
        assert tai_utc_after == 37.0

    @pytest.mark.unit
    def test_before_and_after_2012_leap_second(self):
        """TAI-UTC should change across the 2012-07-01 leap second."""
        jd_before = ephem.julday(2012, 6, 30, 12.0)
        jd_after = ephem.julday(2012, 7, 1, 12.0)

        tai_utc_before = ephem.get_tai_utc_for_jd(jd_before)
        tai_utc_after = ephem.get_tai_utc_for_jd(jd_after)

        assert tai_utc_before == 34.0
        assert tai_utc_after == 35.0


class TestTAIHistoricalDates:
    """Test TAI functions with historical dates."""

    @pytest.mark.unit
    def test_year_1972_start_of_leap_seconds(self):
        """Test the beginning of the leap second era (1972)."""
        jd = ephem.julday(1972, 1, 1, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        assert tai_utc == 10.0  # TAI-UTC was 10 seconds on 1972-01-01

    @pytest.mark.unit
    def test_year_1980(self):
        """Test TAI-UTC for 1980."""
        jd = ephem.julday(1980, 6, 1, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)
        assert tai_utc == 19.0

    @pytest.mark.unit
    def test_utc_to_tai_1980(self):
        """Test UTC to TAI conversion for 1980."""
        jd_utc = ephem.julday(1980, 6, 1, 12.0)
        jd_tai = ephem.utc_jd_to_tai(jd_utc)
        diff_seconds = (jd_tai - jd_utc) * 86400
        assert diff_seconds == pytest.approx(19.0, abs=0.1)
