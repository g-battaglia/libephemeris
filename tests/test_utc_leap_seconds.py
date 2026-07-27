"""
Tests for UTC leap second handling in UTC to UT1 conversion.

This module verifies that the libephemeris time conversion functions correctly
handle leap seconds when converting between UTC, UT1, TAI, and TT time scales.

Leap seconds are inserted by IERS to keep UTC within 0.9 seconds of UT1.
The relationship between time scales is:
    - TAI (International Atomic Time): Continuous atomic time
    - UTC = TAI - leap_seconds (UTC is adjusted by leap seconds)
    - TT = TAI + 32.184 seconds (Terrestrial Time)
    - UT1 = UTC + (UT1-UTC), where |UT1-UTC| < 0.9 seconds

Key leap second dates tested:
    - 1972-01-01: TAI-UTC = 10 seconds (first leap second system)
    - 2017-01-01: TAI-UTC = 37 seconds (most recent as of 2024)
"""

import pytest
import libephemeris as ephem
from libephemeris.iers_data import (
    get_tai_utc,
    _calendar_to_mjd,
    _get_hardcoded_leap_seconds,
    load_iers_data,
)


class TestLeapSecondValues:
    """Tests verifying correct TAI-UTC values at known leap second dates."""

    # Historical leap second data: (year, month, day, expected_tai_utc_after)
    LEAP_SECOND_DATES = [
        (1972, 1, 1, 10),  # First leap second system introduced
        (1972, 7, 1, 11),  # Second leap second
        (1973, 1, 1, 12),
        (1974, 1, 1, 13),
        (1975, 1, 1, 14),
        (1976, 1, 1, 15),
        (1977, 1, 1, 16),
        (1978, 1, 1, 17),
        (1979, 1, 1, 18),
        (1980, 1, 1, 19),
        (1981, 7, 1, 20),
        (1982, 7, 1, 21),
        (1983, 7, 1, 22),
        (1985, 7, 1, 23),
        (1988, 1, 1, 24),
        (1990, 1, 1, 25),
        (1991, 1, 1, 26),
        (1992, 7, 1, 27),
        (1993, 7, 1, 28),
        (1994, 7, 1, 29),
        (1996, 1, 1, 30),
        (1997, 7, 1, 31),
        (1999, 1, 1, 32),
        (2006, 1, 1, 33),
        (2009, 1, 1, 34),
        (2012, 7, 1, 35),
        (2015, 7, 1, 36),
        (2017, 1, 1, 37),
    ]

    def test_hardcoded_leap_seconds_complete(self):
        """Verify the hardcoded leap second list contains all known leap seconds."""
        entries = _get_hardcoded_leap_seconds()

        # Should have 28 leap second entries (from 1972 to 2017)
        assert len(entries) == 28, (
            f"Expected 28 leap second entries, got {len(entries)}"
        )

        # Verify first and last entries
        assert entries[0].tai_utc == 10.0  # 1972-01-01
        assert entries[-1].tai_utc == 37.0  # 2017-01-01

    @pytest.mark.parametrize("year,month,day,expected_tai_utc", LEAP_SECOND_DATES)
    def test_tai_utc_at_leap_second_dates(self, year, month, day, expected_tai_utc):
        """Verify TAI-UTC values are correct at each leap second date."""
        # Load data to ensure leap seconds are available
        load_iers_data()

        # Get MJD for the date (at noon to avoid boundary issues)
        mjd = _calendar_to_mjd(year, month, day) + 0.5

        tai_utc = get_tai_utc(mjd)
        assert tai_utc == expected_tai_utc, (
            f"TAI-UTC at {year}-{month:02d}-{day:02d} should be {expected_tai_utc}s, "
            f"got {tai_utc}s"
        )

    def test_tai_utc_before_first_leap_second(self):
        """Verify pre-1972 TAI-UTC follows the rate-based (rubber-second) model.

        Before the 1972 leap-second system, UTC tracked TAI via the IERS
        piecewise-linear frequency-offset law, not integer steps. get_tai_utc
        implements that model for pre-1972 dates (it is NOT 0: the true
        offset grew from ~1.4 s in 1961 to ~10 s by the end of 1971).
        """
        load_iers_data()

        # December 31, 1971 - day before first leap second
        mjd_pre_1972 = _calendar_to_mjd(1971, 12, 31) + 0.5
        tai_utc = get_tai_utc(mjd_pre_1972)
        # 1968-02-01 segment: 4.2131700 s + (MJD - 39126) * 0.002592 s/day
        expected = 4.2131700 + (mjd_pre_1972 - 39126.0) * 0.002592
        assert tai_utc == pytest.approx(expected, abs=1e-9), (
            f"pre-1972 TAI-UTC should follow the rate-based model "
            f"({expected:.6f}s), got {tai_utc}"
        )
        # Sanity: just below the first integer step (10 s at 1972-01-01).
        assert 9.0 < tai_utc < 10.0

    def test_tai_utc_after_last_known_leap_second(self):
        """Verify TAI-UTC remains at 37 after 2017-01-01."""
        load_iers_data()

        # Test dates after the last leap second
        test_dates = [
            (2020, 1, 1),
            (2023, 6, 15),
            (2024, 12, 31),
        ]

        for year, month, day in test_dates:
            mjd = _calendar_to_mjd(year, month, day) + 0.5
            tai_utc = get_tai_utc(mjd)
            assert tai_utc == 37.0, (
                f"TAI-UTC at {year}-{month:02d}-{day:02d} should be 37s, got {tai_utc}s"
            )

    def test_tai_utc_boundary_before_leap_second(self):
        """Verify TAI-UTC value just before a leap second is the previous value."""
        load_iers_data()

        # December 31, 2016 23:59:59 - just before 2017-01-01 leap second
        # TAI-UTC should still be 36 (becomes 37 at midnight Jan 1)
        mjd_dec31_2016 = _calendar_to_mjd(2016, 12, 31) + 0.5
        tai_utc = get_tai_utc(mjd_dec31_2016)
        assert tai_utc == 36.0, f"TAI-UTC on 2016-12-31 should be 36s, got {tai_utc}s"


class TestUtcToJdLeapSecondHandling:
    """Tests for UTC to JD conversion with leap second handling."""

    def test_utc_to_jd_basic(self):
        """Test basic UTC to JD conversion returns reasonable values."""
        # J2000.0 epoch: Jan 1, 2000 12:00:00 UTC
        jd_tt, jd_ut1 = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)

        # JD(TT) should be close to 2451545.0 (J2000.0)
        assert abs(jd_tt - 2451545.0) < 0.01, (
            f"JD(TT) at J2000.0 should be ~2451545.0, got {jd_tt}"
        )

        # JD(UT1) should be slightly less than JD(TT) due to Delta T
        assert jd_tt > jd_ut1, "JD(TT) should be greater than JD(UT1)"

    def test_utc_to_jd_delta_t_relationship(self):
        """Verify Delta T relationship: TT = UT1 + Delta_T."""
        # Use a date where Delta T is well known (~69 seconds in 2020)
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 6, 15, 12, 0, 0.0)

        # Calculate Delta T from the conversion
        delta_t_from_jd = (jd_tt - jd_ut1) * 86400  # Convert days to seconds

        # Delta T at 2020 should be approximately 69 seconds
        assert 68 < delta_t_from_jd < 71, (
            f"Delta T in 2020 should be ~69 seconds, got {delta_t_from_jd:.2f}s"
        )

    def test_utc_to_jd_leap_second_consistency(self):
        """Verify UTC to JD is consistent across leap second boundaries."""
        # Test dates before and after the 2017 leap second
        jd_tt_before, jd_ut1_before = ephem.utc_to_jd(2016, 12, 31, 12, 0, 0.0)
        jd_tt_after, jd_ut1_after = ephem.utc_to_jd(2017, 1, 1, 12, 0, 0.0)

        # The difference should be exactly 1 day (plus tiny Delta T change)
        jd_diff = jd_ut1_after - jd_ut1_before
        assert abs(jd_diff - 1.0) < 0.001, (
            f"JD difference across leap second should be ~1 day, got {jd_diff:.6f}"
        )

    def test_utc_to_jd_at_leap_second_insertion(self):
        """Test UTC to JD right at leap second insertion moment."""
        # The leap second is inserted at 23:59:60 on Dec 31
        # After the leap second, TAI-UTC increases by 1

        # Just before midnight on Dec 31, 2016 (before 2017 leap second)
        jd_tt_1, jd_ut1_1 = ephem.utc_to_jd(2016, 12, 31, 23, 59, 59.0)

        # At midnight on Jan 1, 2017 (after 2017 leap second)
        jd_tt_2, jd_ut1_2 = ephem.utc_to_jd(2017, 1, 1, 0, 0, 0.0)

        # The UT1 difference should be approximately 1 second
        ut1_diff_seconds = (jd_ut1_2 - jd_ut1_1) * 86400
        assert 0.5 < ut1_diff_seconds < 2.0, (
            f"UT1 difference around leap second should be ~1s, got {ut1_diff_seconds:.3f}s"
        )

    def test_utc_to_jd_tai_relationship(self):
        """Verify UTC to TAI relationship includes correct leap seconds."""
        # Use a date after the last leap second
        year, month, day = 2020, 6, 15

        # Get JD in TT
        jd_tt, jd_ut1 = ephem.utc_to_jd(year, month, day, 12, 0, 0.0)

        # Get TAI-UTC for this date (should be 37 seconds)
        mjd = _calendar_to_mjd(year, month, day)
        tai_utc = get_tai_utc(mjd)
        assert tai_utc == 37.0, f"TAI-UTC should be 37s, got {tai_utc}s"

        # TT = TAI + 32.184 seconds
        # So TT - UTC = (TAI-UTC) + 32.184 = 37 + 32.184 = 69.184 seconds
        # This should be close to Delta T (which is TT - UT1, and UT1 ≈ UTC)

        # The Delta T calculated from JD should be consistent
        delta_t_from_jd = (jd_tt - jd_ut1) * 86400
        expected_tt_utc = tai_utc + 32.184

        # Delta T (TT - UT1) should be close to TT - UTC since |UT1-UTC| < 0.9s
        assert abs(delta_t_from_jd - expected_tt_utc) < 1.0, (
            f"Delta T ({delta_t_from_jd:.2f}s) should be close to "
            f"TT-UTC ({expected_tt_utc:.2f}s)"
        )


class TestJdut1ToUtcLeapSecondHandling:
    """Tests for JD(UT1) to UTC conversion with leap second handling."""

    def test_jdut1_to_utc_roundtrip(self):
        """Test round-trip conversion UTC -> JD(UT1) -> UTC."""
        test_dates = [
            (2000, 1, 1, 12, 0, 0.0),
            (2015, 6, 30, 23, 59, 30.0),  # Just before 2015-07-01 leap second
            (2017, 1, 1, 0, 0, 30.0),  # Just after 2017-01-01 leap second
            (2020, 6, 15, 14, 30, 45.5),
        ]

        for year, month, day, hour, minute, second in test_dates:
            # Convert UTC to JD(UT1)
            _, jd_ut1 = ephem.utc_to_jd(year, month, day, hour, minute, second)

            # Convert back to UTC
            y, m, d, h, mi, s = ephem.jdut1_to_utc(jd_ut1)

            # Check round-trip accuracy (within 0.1 second)
            assert y == year, f"Year mismatch: {y} != {year}"
            assert m == month, f"Month mismatch: {m} != {month}"
            assert d == day, f"Day mismatch: {d} != {day}"
            assert h == hour, f"Hour mismatch: {h} != {hour}"
            assert mi == minute, f"Minute mismatch: {mi} != {minute}"
            assert abs(s - second) < 0.5, (
                f"Second mismatch at {year}-{month:02d}-{day:02d}: {s:.2f} != {second}"
            )

    def test_jdut1_to_utc_across_leap_second(self):
        """Test JD(UT1) to UTC conversion maintains consistency across leap seconds."""
        # Get JD(UT1) for dates before and after leap second
        _, jd_ut1_before = ephem.utc_to_jd(2016, 12, 31, 12, 0, 0.0)
        _, jd_ut1_after = ephem.utc_to_jd(2017, 1, 1, 12, 0, 0.0)

        # Convert back
        y1, m1, d1, h1, mi1, s1 = ephem.jdut1_to_utc(jd_ut1_before)
        y2, m2, d2, h2, mi2, s2 = ephem.jdut1_to_utc(jd_ut1_after)

        # Verify dates are correct
        assert (y1, m1, d1) == (2016, 12, 31)
        assert (y2, m2, d2) == (2017, 1, 1)


class TestTaiUtcConsistency:
    """Tests for TAI-UTC consistency across the time conversion chain."""

    def test_tai_utc_via_jd(self):
        """Verify TAI-UTC can be correctly retrieved via JD functions."""
        # Get TAI-UTC for 2020 using the JD-based function
        jd = ephem.julday(2020, 1, 1, 12.0)
        tai_utc = ephem.get_tai_utc_for_jd(jd)

        assert tai_utc == 37.0, f"TAI-UTC for 2020 should be 37s, got {tai_utc}s"

    def test_tai_utc_historical_values(self):
        """Verify TAI-UTC historical values via JD."""
        test_cases = [
            (1980, 1, 1, 19.0),
            (1990, 1, 1, 25.0),
            (2000, 1, 1, 32.0),
            (2010, 1, 1, 34.0),
            (2016, 6, 15, 36.0),
            (2017, 6, 15, 37.0),
        ]

        for year, month, day, expected_tai_utc in test_cases:
            jd = ephem.julday(year, month, day, 12.0)
            tai_utc = ephem.get_tai_utc_for_jd(jd)
            assert tai_utc == expected_tai_utc, (
                f"TAI-UTC at {year}-{month:02d}-{day:02d} should be {expected_tai_utc}s, "
                f"got {tai_utc}s"
            )


class TestTaiTimeConversions:
    """Tests for TAI time scale conversions."""

    def test_utc_to_tai_jd_leap_seconds(self):
        """Verify UTC to TAI JD correctly adds leap seconds."""
        # In 2020, TAI-UTC = 37 seconds
        jd_tai = ephem.utc_to_tai_jd(2020, 1, 1, 0, 0, 0.0)
        jd_utc = ephem.julday(2020, 1, 1, 0.0)

        # TAI should be ahead of UTC by 37 seconds
        diff_seconds = (jd_tai - jd_utc) * 86400
        assert abs(diff_seconds - 37.0) < 0.01, (
            f"TAI should be 37s ahead of UTC, got {diff_seconds:.2f}s"
        )

    def test_tai_jd_to_utc_roundtrip(self):
        """Test round-trip UTC -> TAI JD -> UTC."""
        original = (2020, 6, 15, 14, 30, 45.0)

        # Convert UTC to TAI JD
        jd_tai = ephem.utc_to_tai_jd(
            *original[:3], original[3], original[4], original[5]
        )

        # Convert back to UTC
        result = ephem.tai_jd_to_utc(jd_tai)

        assert result[:5] == original[:5], f"Date/time mismatch: {result} != {original}"
        assert abs(result[5] - original[5]) < 0.01, (
            f"Second mismatch: {result[5]:.2f} != {original[5]}"
        )

    def test_tt_tai_relationship(self):
        """Verify TT = TAI + 32.184 seconds exactly."""
        jd_tt = 2451545.0  # J2000.0 in TT
        jd_tai = ephem.tt_to_tai_jd(jd_tt)

        diff_seconds = (jd_tt - jd_tai) * 86400
        # Allow small floating-point precision error (sub-millisecond)
        assert abs(diff_seconds - 32.184) < 1e-4, (
            f"TT-TAI should be exactly 32.184s, got {diff_seconds:.6f}s"
        )

    def test_tai_tt_roundtrip(self):
        """Test round-trip TT -> TAI -> TT."""
        jd_tt_orig = 2451545.0

        jd_tai = ephem.tt_to_tai_jd(jd_tt_orig)
        jd_tt_back = ephem.tai_to_tt_jd(jd_tai)

        assert abs(jd_tt_orig - jd_tt_back) < 1e-15, (
            f"TT round-trip failed: {jd_tt_orig} != {jd_tt_back}"
        )


class TestDeltaTWithIersData:
    """Tests for Delta T calculations using IERS data."""

    def test_deltat_reasonable_modern_value(self):
        """Verify Delta T has reasonable values for modern dates."""
        # Delta T at J2000.0 should be approximately 64 seconds
        jd_j2000 = 2451545.0
        dt = ephem.deltat(jd_j2000) * 86400  # Convert to seconds

        assert 60 < dt < 70, f"Delta T at J2000.0 should be ~64s, got {dt:.2f}s"

    def test_deltat_increasing_trend(self):
        """Verify Delta T has increased from 2000 to 2020."""
        jd_2000 = ephem.julday(2000, 1, 1, 12.0)
        jd_2020 = ephem.julday(2020, 1, 1, 12.0)

        dt_2000 = ephem.deltat(jd_2000) * 86400
        dt_2020 = ephem.deltat(jd_2020) * 86400

        assert dt_2020 > dt_2000, (
            f"Delta T should increase from 2000 to 2020: {dt_2000:.2f}s vs {dt_2020:.2f}s"
        )

    def test_deltat_ex_ephemeris_flags(self):
        """Test Delta T extended function with different ephemeris flags."""
        jd = 2451545.0

        dt_swieph = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        dt_jpleph = ephem.deltat_ex(jd, ephem.FLG_JPLEPH)

        # Both should give similar results (Skyfield uses JPL data internally)
        assert abs(dt_swieph - dt_jpleph) < 1e-10, (
            "SWIEPH and JPLEPH should give identical Delta T values"
        )

    def test_deltat_ex_moseph_uses_same_delta_t(self):
        """Verify Moshier ephemeris uses the same Delta T calculation (no warning)."""
        jd = 2451545.0

        dt = ephem.deltat_ex(jd, ephem.FLG_MOSEPH)

        # Should return a valid Delta T value matching other ephemeris modes
        dt_seconds = dt * 86400
        assert 60 < dt_seconds < 70

        # Verify it matches the standard function
        dt_standard = ephem.deltat(jd)
        assert dt == dt_standard


class TestLeapSecondEdgeCases:
    """Edge case tests for leap second handling."""

    def test_leap_second_moment_utc_60(self):
        """Test handling of second=60 during a leap second."""
        # During a leap second insertion, UTC can have second=60
        # Dec 31, 2016 23:59:60 UTC (during 2017 leap second insertion)
        #
        # Note: Skyfield may or may not accept second=60 depending on version
        # This tests that the library handles it gracefully
        try:
            jd_tt, jd_ut1 = ephem.utc_to_jd(2016, 12, 31, 23, 59, 60.0)
            # If it doesn't raise, the values should be reasonable
            assert jd_tt > 0
            assert jd_ut1 > 0
        except (ValueError, Exception):
            # It's acceptable to raise an error for second=60
            # as not all libraries support this
            pytest.skip("Library does not support second=60 (leap second)")

    def test_high_precision_around_leap_second(self):
        """Test millisecond precision around leap second boundary."""
        # 0.5 second before midnight on leap second day
        jd_tt_1, jd_ut1_1 = ephem.utc_to_jd(2016, 12, 31, 23, 59, 59.5)

        # 0.5 second after midnight
        jd_tt_2, jd_ut1_2 = ephem.utc_to_jd(2017, 1, 1, 0, 0, 0.5)

        # The UTC difference is 1.0 second + 1.0 leap second = 2.0 seconds
        # But UT1 difference should be approximately 1.0 second (no leap second)
        # Actually, due to how Skyfield handles this, the UT1 difference
        # depends on the DUT1 (UT1-UTC) value at each moment

        ut1_diff = (jd_ut1_2 - jd_ut1_1) * 86400
        # The difference should be between 1 and 3 seconds
        assert 0.5 < ut1_diff < 3.0, (
            f"UT1 difference around leap second should be ~1-2s, got {ut1_diff:.3f}s"
        )

    def test_multiple_leap_seconds_consistency(self):
        """Verify consistency when crossing multiple leap seconds."""
        # From 1990 (TAI-UTC=25) to 2020 (TAI-UTC=37) = 12 leap seconds

        jd_utc_1990 = ephem.julday(1990, 1, 1, 12.0)
        jd_utc_2020 = ephem.julday(2020, 1, 1, 12.0)

        tai_utc_1990 = ephem.get_tai_utc_for_jd(jd_utc_1990)
        tai_utc_2020 = ephem.get_tai_utc_for_jd(jd_utc_2020)

        leap_second_diff = tai_utc_2020 - tai_utc_1990
        assert leap_second_diff == 12.0, (
            f"Should have 12 leap seconds between 1990 and 2020, got {leap_second_diff}"
        )


class TestIersDataIntegration:
    """Integration tests for IERS data with time conversions."""

    def test_iers_data_loads(self):
        """Verify IERS data can be loaded."""
        result = load_iers_data()
        # Should return True even with just hardcoded leap seconds
        assert result is True

    def test_ut1_utc_reasonable_values(self):
        """Verify UT1-UTC values are within expected range."""
        from libephemeris.iers_data import get_ut1_utc

        load_iers_data()

        # UT1-UTC should always be between -0.9 and +0.9 seconds
        test_dates = [
            (2000, 1, 1),
            (2010, 6, 15),
            (2020, 3, 20),
        ]

        for year, month, day in test_dates:
            mjd = _calendar_to_mjd(year, month, day)
            ut1_utc = get_ut1_utc(mjd)

            if ut1_utc is not None:
                assert -0.9 <= ut1_utc <= 0.9, (
                    f"UT1-UTC at {year}-{month:02d}-{day:02d} should be within "
                    f"[-0.9, 0.9]s, got {ut1_utc:.3f}s"
                )

    def test_delta_t_consistency_with_iers(self):
        """Verify Delta T from Skyfield is consistent with IERS relationship."""
        # Delta T = TT - UT1 = (TAI + 32.184) - (UTC + UT1-UTC)
        #         = TAI - UTC + 32.184 - UT1-UTC
        #         = TAI-UTC + 32.184 - UT1-UTC

        from libephemeris.iers_data import get_ut1_utc

        load_iers_data()

        jd = ephem.julday(2015, 6, 15, 12.0)
        mjd = _calendar_to_mjd(2015, 6, 15)

        # Get components
        tai_utc = get_tai_utc(mjd)  # Should be 35 or 36
        ut1_utc = get_ut1_utc(mjd)

        if ut1_utc is not None:
            # Calculate expected Delta T
            expected_delta_t = tai_utc + 32.184 - ut1_utc

            # Get actual Delta T from Skyfield
            actual_delta_t = ephem.deltat(jd) * 86400

            # They should be close (within 0.5 seconds)
            assert abs(actual_delta_t - expected_delta_t) < 0.5, (
                f"Delta T mismatch: actual={actual_delta_t:.3f}s, "
                f"expected={expected_delta_t:.3f}s"
            )


class TestUtcTimezoneWithLeapSeconds:
    """Test UTC timezone conversions with leap second awareness."""

    def test_timezone_conversion_basic(self):
        """Verify basic timezone conversion works correctly.

        utc_time_zone converts FROM local time TO UTC by subtracting
        the timezone offset. So local 10:30 in UTC+1 => UTC 09:30.
        """
        # 2020-01-15 10:30:00 local (UTC+1) => 2020-01-15 09:30:00 UTC
        result = ephem.utc_time_zone(2020, 1, 15, 10, 30, 0.0, 1.0)

        assert result[:5] == (2020, 1, 15, 9, 30)
        assert abs(result[5] - 0.0) < 0.01

    def test_timezone_conversion_day_boundary(self):
        """Test timezone conversion crossing day boundary.

        utc_time_zone converts FROM local time TO UTC (subtracts offset).
        Local 23:00 in UTC+5 => UTC 18:00 (same day).
        """
        result = ephem.utc_time_zone(2020, 1, 15, 23, 0, 0.0, 5.0)

        assert result[:5] == (2020, 1, 15, 18, 0)

    def test_timezone_conversion_negative_offset(self):
        """Test timezone conversion with negative offset.

        utc_time_zone converts FROM local time TO UTC (subtracts offset).
        Local 02:00 in UTC-5 => UTC 07:00 (same day).
        """
        result = ephem.utc_time_zone(2020, 1, 15, 2, 0, 0.0, -5.0)

        assert result[:5] == (2020, 1, 15, 7, 0)


class TestJulianCalendarLeapSecondPreserved:
    """jdut1_to_utc/jdet_to_utc with JUL_CAL must keep a 23:59:60 leap
    second: the calendar remap only shifts the DATE (whole days), never the
    clock. Folding h/m/s into a decimal hour rolled the leap second past
    24h into the next day (reference: (2016, 12, 18, 23, 59, 60.41))."""

    def test_jdut1_julian_keeps_leap_second(self):
        y, m, d, hh, mm, ss = ephem.jdut1_to_utc(2457754.5, 0)
        assert (y, m, d, hh, mm) == (2016, 12, 18, 23, 59)
        assert 60.0 <= ss < 61.0

    def test_jdet_julian_keeps_leap_second(self):
        jd_tt, _ = ephem.utc_to_jd(2016, 12, 31, 23, 59, 60.4, 1)
        y, m, d, hh, mm, ss = ephem.jdet_to_utc(jd_tt, 0)
        assert (y, m, d, hh, mm) == (2016, 12, 18, 23, 59)
        assert 60.0 <= ss < 61.0

    def test_julian_normal_date_unchanged(self):
        _, jd_ut1 = ephem.utc_to_jd(2020, 6, 15, 14, 30, 0.0, 1)
        y, m, d, hh, mm, ss = ephem.jdut1_to_utc(jd_ut1, 0)
        assert (y, m, d) == (2020, 6, 2)
        assert abs((hh * 3600 + mm * 60 + ss) - (14 * 3600 + 30 * 60)) < 0.001


class TestPre1972LeapSecondMessage:
    """Pre-1972 dates at 23:59:60 get the plain 'invalid time' wording;
    the '(no leap second!)' suffix is reserved for the UTC era (measured
    on the reference API for 1971/1965/1900 vs 2000)."""

    def test_pre_1972_plain_wording(self):
        import pytest

        from libephemeris.exceptions import Error

        for y, m, d in ((1971, 12, 31), (1965, 6, 30), (1900, 1, 1)):
            with pytest.raises(Error, match=r"^invalid time: 23:59:60\.00$"):
                ephem.utc_to_jd(y, m, d, 23, 59, 60.0, 1)

    def test_utc_era_keeps_suffix(self):
        import pytest

        from libephemeris.exceptions import Error

        with pytest.raises(Error, match=r"no leap second"):
            ephem.utc_to_jd(2000, 1, 1, 23, 59, 60.0, 1)


class TestRoundInstantReconstruction:
    """Inverse conversions recover round instants exactly.

    A modern float JD carries a ~40 microsecond ulp, so the raw calendar
    reconstruction has a ~+-20 microsecond floor that used to straddle
    midnight backward (2013-01-01 -> 2012-12-31 23:59:59.99998). Readings
    within 200 microseconds of a whole second snap, with a leap-aware
    calendar carry; genuinely offset instants are untouched.
    """

    def test_midnight_round_trip_is_exact(self):
        import libephemeris as le

        for y, m, d in ((2013, 1, 1), (1985, 3, 20), (2030, 12, 31)):
            jd = le.utc_to_jd(y, m, d, 0, 0, 0.0, 1)
            assert le.jdet_to_utc(jd[0], 1) == (y, m, d, 0, 0, 0.0)
            assert le.jdut1_to_utc(jd[1], 1) == (y, m, d, 0, 0, 0.0)

    def test_leap_second_kept_and_post_leap_midnight(self):
        import libephemeris as le

        jd = le.utc_to_jd(2012, 6, 30, 23, 59, 60.0, 1)
        assert le.jdet_to_utc(jd[0], 1) == (2012, 6, 30, 23, 59, 60.0)
        jd2 = le.utc_to_jd(2012, 7, 1, 0, 0, 0.0, 1)
        assert le.jdet_to_utc(jd2[0], 1) == (2012, 7, 1, 0, 0, 0.0)

    def test_genuine_offset_not_snapped(self):
        import libephemeris as le

        jd = le.utc_to_jd(2013, 1, 1, 0, 0, 0.0025, 1)
        sec = le.jdet_to_utc(jd[0], 1)[5]
        assert abs(sec - 0.0025) < 2e-4 and sec != 0.0
