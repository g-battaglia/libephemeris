"""
Tests for TDB (Barycentric Dynamical Time) and TT (Terrestrial Time) handling.

This test module verifies that the library correctly handles the distinction
between TDB and TT where relevant.

Background:
-----------
- TT (Terrestrial Time): A theoretical ideal atomic clock on Earth's geoid.
  Used as input for the calc() function (pyswisseph-compatible API).
- TDB (Barycentric Dynamical Time): Time used for planetary ephemerides.
  It's the independent variable in the equations of motion.
- The difference TDB - TT is a periodic function with maximum amplitude of
  approximately 1.66 milliseconds (about 0.001657 seconds).

Design Decision:
----------------
This library treats TDB ≈ TT for the following reasons:
1. The maximum difference (~1.7 ms) is negligible for typical astrological
   calculations which require arcminute precision at best.
2. Skyfield uses TDB internally for ephemeris lookups (JPL DE4xx files).
3. For sub-arcsecond precision work, the difference matters but such
   precision is beyond typical use cases.

These tests document and verify this design decision.

References:
-----------
- IAU SOFA: Time Scale Transformations
- USNO Circular 179: The Astronomical Almanac
- JPL: Explanatory Supplement to the Astronomical Almanac
"""

import math
import pytest
import libephemeris as ephem
from libephemeris.state import get_timescale


class TestTDBTTBasics:
    """Test basic TDB/TT relationship and properties."""

    @pytest.mark.unit
    def test_skyfield_provides_both_tt_and_tdb(self):
        """Skyfield Time objects provide both TT and TDB Julian Day."""
        ts = get_timescale()
        t = ts.tt_jd(2451545.0)  # J2000.0 in TT

        # Skyfield should provide both time scales
        assert hasattr(t, "tt")
        assert hasattr(t, "tdb")
        assert isinstance(float(t.tt), float)
        assert isinstance(float(t.tdb), float)

    @pytest.mark.unit
    def test_tdb_tt_difference_small(self):
        """TDB - TT difference should be less than 2 milliseconds."""
        ts = get_timescale()

        # Test several dates throughout a year (TDB-TT varies periodically)
        test_dates = [
            (2000, 1, 1, 12.0),
            (2000, 4, 1, 12.0),
            (2000, 7, 1, 12.0),
            (2000, 10, 1, 12.0),
        ]

        for year, month, day, hour in test_dates:
            jd_tt = ephem.julday(year, month, day, hour)
            t = ts.tt_jd(jd_tt)

            # TDB and TT values in Julian Days
            tdb_jd = float(t.tdb)
            tt_jd = float(t.tt)

            # Difference in seconds
            diff_seconds = abs(tdb_jd - tt_jd) * 86400

            # Maximum TDB - TT is approximately 1.66 ms
            # Allow up to 2 ms to account for any variations
            assert diff_seconds < 0.002, (
                f"TDB - TT = {diff_seconds * 1000:.4f} ms at {year}-{month:02d}-{day:02d}, "
                f"exceeds expected maximum of ~1.7 ms"
            )

    @pytest.mark.unit
    def test_tdb_tt_periodic_variation(self):
        """TDB - TT varies periodically over the year."""
        ts = get_timescale()

        # Sample multiple points through a year
        jd_start = ephem.julday(2020, 1, 1, 12.0)
        differences = []

        for day_offset in range(0, 365, 30):  # Sample every 30 days
            jd = jd_start + day_offset
            t = ts.tt_jd(jd)
            diff_seconds = (float(t.tdb) - float(t.tt)) * 86400
            differences.append(diff_seconds)

        # The differences should vary (not constant)
        # Check that range is non-zero
        diff_range = max(differences) - min(differences)
        assert diff_range > 0.0001, (
            "TDB - TT should show periodic variation throughout the year"
        )

    @pytest.mark.unit
    def test_tdb_tt_difference_sub_millisecond_impact(self):
        """Verify the impact of TDB/TT difference on positions is sub-arcsecond."""
        # The TDB-TT difference of ~1.7 ms means:
        # - Moon moves about 0.5 arcseconds per second, so ~0.001 arcsec in 1.7 ms
        # - Fastest planet (Mercury) moves about 4 degrees/day = 0.17 arcsec/sec
        #   In 1.7 ms: ~0.0003 arcseconds
        # This is negligible for astrological purposes (arcminute precision)

        # Calculate position at same instant using TT and TT+1.7ms
        jd_tt = 2451545.0  # J2000.0
        dt_tdb_tt = 0.0017 / 86400  # 1.7 ms in days

        # Get Moon position (fastest moving body)
        pos_tt, _ = ephem.calc(jd_tt, ephem.MOON, ephem.FLG_SPEED)
        pos_offset, _ = ephem.calc(
            jd_tt + dt_tdb_tt, ephem.MOON, ephem.FLG_SPEED
        )

        # Position difference in arcseconds
        lon_diff_arcsec = abs(pos_tt[0] - pos_offset[0]) * 3600

        # Should be well under 1 arcsecond
        assert lon_diff_arcsec < 1.0, (
            f"Moon position difference for 1.7ms time offset = {lon_diff_arcsec:.4f} arcsec"
        )


class TestSweCalcTimeScales:
    """Test that calc and calc_ut use correct time scales."""

    @pytest.mark.unit
    def test_swe_calc_uses_tt_input(self):
        """calc() should accept TT as input."""
        jd_tt = 2451545.0  # J2000.0 in TT

        # This should work without error
        pos, flag = ephem.calc(jd_tt, ephem.SUN, ephem.FLG_SPEED)

        # Should return valid position
        assert len(pos) == 6
        assert 0 <= pos[0] < 360  # Longitude in range
        assert -90 <= pos[1] <= 90  # Latitude in range
        assert pos[2] > 0  # Distance positive

    @pytest.mark.unit
    def test_swe_calc_ut_uses_ut1_input(self):
        """calc_ut() should accept UT1 as input."""
        jd_ut = 2451545.0  # J2000.0 as UT1

        # This should work without error
        pos, flag = ephem.calc_ut(jd_ut, ephem.SUN, ephem.FLG_SPEED)

        # Should return valid position
        assert len(pos) == 6
        assert 0 <= pos[0] < 360  # Longitude in range
        assert -90 <= pos[1] <= 90  # Latitude in range
        assert pos[2] > 0  # Distance positive

    @pytest.mark.unit
    def test_swe_calc_vs_swe_calc_ut_consistency(self):
        """calc(TT) should match calc_ut(UT) when adjusted by Delta T."""
        jd_ut = ephem.julday(2000, 1, 1, 12.0)
        delta_t = ephem.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        # Calculate using both methods
        pos_tt, _ = ephem.calc(jd_tt, ephem.SUN, ephem.FLG_SPEED)
        pos_ut, _ = ephem.calc_ut(jd_ut, ephem.SUN, ephem.FLG_SPEED)

        # Positions should match closely (within a few arcseconds due to
        # the way Delta T is computed at slightly different times)
        lon_diff_arcsec = abs(pos_tt[0] - pos_ut[0]) * 3600
        assert lon_diff_arcsec < 1.0, (
            f"calc(TT) and calc_ut(UT) differ by {lon_diff_arcsec:.4f} arcsec"
        )

    @pytest.mark.unit
    def test_swe_calc_and_swe_calc_ut_differ_by_delta_t(self):
        """Using same JD with calc and calc_ut should show Delta T offset."""
        jd = 2451545.0  # Same numeric value for both
        delta_t = ephem.deltat(jd)
        delta_t_seconds = delta_t * 86400

        # Calculate Sun position with both functions using same numeric JD
        pos_tt, _ = ephem.calc(jd, ephem.SUN, ephem.FLG_SPEED)
        pos_ut, _ = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)

        # They should differ because one interprets JD as TT, other as UT1
        # Sun moves about 1 degree per day, Delta T is ~64 seconds
        # Expected difference: 1 deg/day * 64s / 86400s = ~0.00074 deg = 2.7 arcsec
        lon_diff_arcsec = abs(pos_tt[0] - pos_ut[0]) * 3600

        # Should show measurable difference from Delta T
        expected_diff_min = 1.0  # At least 1 arcsec difference
        expected_diff_max = 10.0  # But less than 10 arcsec

        assert expected_diff_min < lon_diff_arcsec < expected_diff_max, (
            f"Expected Sun position difference due to Delta T = ~2.7 arcsec, "
            f"got {lon_diff_arcsec:.2f} arcsec"
        )


class TestUtcToJdTimeScales:
    """Test time scale handling in utc_to_jd."""

    @pytest.mark.unit
    def test_utc_to_jd_returns_both_tt_and_ut1(self):
        """utc_to_jd() returns tuple (jd_tt, jd_ut1)."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2000, 1, 1, 12, 0, 0.0)

        assert isinstance(jd_tt, float)
        assert isinstance(jd_ut1, float)
        # TT is always ahead of UT1 for modern dates
        assert jd_tt > jd_ut1

    @pytest.mark.unit
    def test_utc_to_jd_tt_for_swe_calc(self):
        """jd_tt from utc_to_jd is suitable for calc()."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 6, 15, 12, 0, 0.0)

        # Should work correctly with calc
        pos, flag = ephem.calc(jd_tt, ephem.VENUS, ephem.FLG_SPEED)
        assert len(pos) == 6
        assert 0 <= pos[0] < 360

    @pytest.mark.unit
    def test_utc_to_jd_ut1_for_swe_calc_ut(self):
        """jd_ut1 from utc_to_jd is suitable for calc_ut()."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 6, 15, 12, 0, 0.0)

        # Should work correctly with calc_ut
        pos, flag = ephem.calc_ut(jd_ut1, ephem.VENUS, ephem.FLG_SPEED)
        assert len(pos) == 6
        assert 0 <= pos[0] < 360

    @pytest.mark.unit
    def test_utc_to_jd_consistency_with_swe_calc_methods(self):
        """Using jd_tt with calc and jd_ut1 with calc_ut should match."""
        jd_tt, jd_ut1 = ephem.utc_to_jd(2020, 6, 15, 12, 0, 0.0)

        pos_from_tt, _ = ephem.calc(jd_tt, ephem.VENUS, ephem.FLG_SPEED)
        pos_from_ut, _ = ephem.calc_ut(jd_ut1, ephem.VENUS, ephem.FLG_SPEED)

        # Both should give same position (within floating-point tolerance)
        lon_diff_arcsec = abs(pos_from_tt[0] - pos_from_ut[0]) * 3600
        assert lon_diff_arcsec < 0.1, (
            f"calc(jd_tt) and calc_ut(jd_ut1) should match, "
            f"diff = {lon_diff_arcsec:.6f} arcsec"
        )


class TestTDBApproximation:
    """Test and document the TDB ≈ TT approximation used in the library."""

    @pytest.mark.unit
    def test_tdb_tt_max_difference(self):
        """Document the maximum TDB - TT difference throughout a year."""
        ts = get_timescale()
        jd_start = ephem.julday(2020, 1, 1, 0.0)

        max_diff_ms = 0
        max_diff_date = None

        # Sample every day for a year
        for day_offset in range(366):
            jd = jd_start + day_offset
            t = ts.tt_jd(jd)
            diff_ms = abs(float(t.tdb) - float(t.tt)) * 86400 * 1000

            if diff_ms > max_diff_ms:
                max_diff_ms = diff_ms
                max_diff_date = ephem.revjul(jd)

        # Maximum should be around 1.6-1.7 ms
        assert 1.0 < max_diff_ms < 2.0, (
            f"Maximum TDB - TT = {max_diff_ms:.4f} ms, expected ~1.66 ms"
        )

    @pytest.mark.unit
    def test_tdb_approximation_acceptable_for_planets(self):
        """Verify TDB ≈ TT is acceptable for planetary positions."""
        ts = get_timescale()
        jd_tt = 2451545.0

        # Get position using TT (as library does)
        t = ts.tt_jd(jd_tt)
        tdb_jd = float(t.tdb)
        tt_jd = float(t.tt)

        # Calculate positions at both times
        for planet, max_speed_arcsec_per_second in [
            (ephem.MOON, 0.55),  # Moon: ~0.55 arcsec/s
            (ephem.MERCURY, 0.05),  # Mercury: ~4 deg/day max
            (ephem.SUN, 0.01),  # Sun: ~1 deg/day
        ]:
            pos_tt, _ = ephem.calc(tt_jd, planet, ephem.FLG_SPEED)
            pos_tdb, _ = ephem.calc(tdb_jd, planet, ephem.FLG_SPEED)

            lon_diff_arcsec = abs(pos_tt[0] - pos_tdb[0]) * 3600

            # Even for Moon, difference should be well under 0.01 arcsec
            assert lon_diff_arcsec < 0.1, (
                f"Planet {planet}: position difference = {lon_diff_arcsec:.6f} arcsec "
                f"for TDB vs TT input"
            )


class TestSkyfieldTimeIntegration:
    """Test Skyfield's time handling integration with the library."""

    @pytest.mark.unit
    def test_skyfield_ut1_to_tt_conversion(self):
        """Verify Skyfield converts UT1 to TT using Delta T."""
        ts = get_timescale()
        jd_ut = 2451545.0

        # Create time from UT1
        t = ts.ut1_jd(jd_ut)

        # TT should be UT1 + Delta T
        expected_delta_t = ephem.deltat(jd_ut)  # in days
        actual_delta_t = float(t.tt) - float(t.ut1)

        # Should match within 1 second
        diff_seconds = abs(expected_delta_t - actual_delta_t) * 86400
        assert diff_seconds < 1.0, (
            f"Delta T mismatch: expected {expected_delta_t * 86400:.2f}s, "
            f"got {actual_delta_t * 86400:.2f}s"
        )

    @pytest.mark.unit
    def test_skyfield_uses_tdb_for_ephemeris(self):
        """Document that Skyfield uses TDB for ephemeris calculations."""
        ts = get_timescale()
        t = ts.tt_jd(2451545.0)

        # Skyfield internally uses TDB for looking up positions in JPL ephemerides
        # The .tdb attribute gives the TDB Julian Day
        tdb_jd = float(t.tdb)

        # TDB should be close to but not exactly equal to TT
        tt_jd = float(t.tt)
        assert tdb_jd != tt_jd, "TDB and TT should differ slightly"
        assert abs(tdb_jd - tt_jd) * 86400 * 1000 < 2, "Difference should be < 2 ms"


class TestJdetToUtcTimeScales:
    """Test jdet_to_utc time scale handling."""

    @pytest.mark.unit
    def test_jdet_to_utc_roundtrip(self):
        """Converting UTC -> TT -> UTC should be reversible."""
        original = (2020, 6, 15, 14, 30, 45.5)
        year, month, day, hour, minute, second = original

        # UTC to JD
        jd_tt, jd_ut = ephem.utc_to_jd(year, month, day, hour, minute, second)

        # JD back to UTC
        result = ephem.jdet_to_utc(jd_tt)
        r_year, r_month, r_day, r_hour, r_minute, r_second = result

        # Should match original
        assert r_year == year
        assert r_month == month
        assert r_day == day
        assert r_hour == hour
        assert r_minute == minute
        assert r_second == pytest.approx(second, abs=0.01)

    @pytest.mark.unit
    def test_jdut1_to_utc_roundtrip(self):
        """Converting UTC -> UT1 -> UTC should be reversible."""
        original = (2020, 6, 15, 14, 30, 45.5)
        year, month, day, hour, minute, second = original

        # UTC to JD
        jd_tt, jd_ut = ephem.utc_to_jd(year, month, day, hour, minute, second)

        # JD back to UTC
        result = ephem.jdut1_to_utc(jd_ut)
        r_year, r_month, r_day, r_hour, r_minute, r_second = result

        # Should match original (within UTC-UT1 tolerance)
        assert r_year == year
        assert r_month == month
        assert r_day == day
        assert r_hour == hour
        assert r_minute == minute
        # UTC and UT1 differ by < 0.9 seconds by definition
        assert r_second == pytest.approx(second, abs=1.0)
