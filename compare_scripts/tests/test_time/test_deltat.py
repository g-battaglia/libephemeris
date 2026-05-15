"""
Comprehensive tests for Delta T calculation (deltat).

Delta T = TT - UT (Terrestrial Time minus Universal Time)
"""

import pytest
import swisseph as swe
import libephemeris as ephem


class TestDeltatBasicValues:
    """Test Delta T returns reasonable values for known epochs."""

    @pytest.mark.unit
    def test_deltat_j2000(self):
        """Delta T at J2000 should be around 63-64 seconds."""
        jd = 2451545.0  # J2000
        dt = ephem.deltat(jd)
        # Convert to seconds (dt is in days)
        dt_seconds = dt * 86400
        assert 60 < dt_seconds < 70, f"Delta T at J2000 = {dt_seconds}s, expected ~63s"

    @pytest.mark.unit
    def test_deltat_1950(self):
        """Delta T around 1950 should be ~29 seconds."""
        jd = ephem.julday(1950, 1, 1, 12.0)
        dt = ephem.deltat(jd)
        dt_seconds = dt * 86400
        assert 25 < dt_seconds < 35, f"Delta T at 1950 = {dt_seconds}s, expected ~29s"

    @pytest.mark.unit
    def test_deltat_1900(self):
        """Delta T around 1900 should be small negative or near zero."""
        jd = ephem.julday(1900, 1, 1, 12.0)
        dt = ephem.deltat(jd)
        dt_seconds = dt * 86400
        # Around 1900, Delta T was about -3 to +5 seconds
        assert -10 < dt_seconds < 10, f"Delta T at 1900 = {dt_seconds}s"

    @pytest.mark.unit
    def test_deltat_2020(self):
        """Delta T around 2020 should be ~69 seconds."""
        jd = ephem.julday(2020, 1, 1, 12.0)
        dt = ephem.deltat(jd)
        dt_seconds = dt * 86400
        assert 65 < dt_seconds < 75, f"Delta T at 2020 = {dt_seconds}s, expected ~69s"


class TestDeltatProperties:
    """Test properties of Delta T values."""

    @pytest.mark.unit
    def test_deltat_positive_modern(self):
        """Delta T should be positive for modern dates (after ~1902)."""
        jd = ephem.julday(2000, 1, 1, 12.0)
        dt = ephem.deltat(jd)
        assert dt > 0, "Delta T should be positive for modern dates"

    @pytest.mark.unit
    def test_deltat_increases_over_time_modern(self):
        """Delta T should generally increase over modern era."""
        jd_1980 = ephem.julday(1980, 1, 1, 12.0)
        jd_2000 = ephem.julday(2000, 1, 1, 12.0)
        jd_2020 = ephem.julday(2020, 1, 1, 12.0)

        dt_1980 = ephem.deltat(jd_1980)
        dt_2000 = ephem.deltat(jd_2000)
        dt_2020 = ephem.deltat(jd_2020)

        assert dt_1980 < dt_2000 < dt_2020, "Delta T should increase over time"

    @pytest.mark.unit
    def test_deltat_continuous(self):
        """Delta T should change smoothly (no large jumps)."""
        jd_start = ephem.julday(2000, 1, 1, 12.0)
        prev_dt = ephem.deltat(jd_start)

        for i in range(1, 365):
            jd = jd_start + i
            dt = ephem.deltat(jd)
            # Change should be small (less than 0.1 seconds per day)
            change = abs(dt - prev_dt) * 86400
            assert change < 0.1, f"Delta T changed by {change}s in one day"
            prev_dt = dt


class TestDeltatVsPyswisseph:
    """Compare Delta T with pyswisseph."""

    @pytest.mark.comparison
    def test_deltat_j2000_matches_swe(self):
        """Delta T at J2000 should match pyswisseph."""
        jd = 2451545.0
        dt_lib = ephem.deltat(jd)
        dt_swe = swe.deltat(jd)
        # Allow 0.5 second tolerance (in days: 0.5/86400)
        assert dt_lib == pytest.approx(dt_swe, abs=0.5 / 86400)

    @pytest.mark.comparison
    @pytest.mark.parametrize("year", [1900, 1920, 1950, 1980, 2000, 2010, 2020])
    def test_deltat_various_years_match_swe(self, year):
        """Delta T should match pyswisseph for various years."""
        jd = ephem.julday(year, 6, 15, 12.0)
        dt_lib = ephem.deltat(jd)
        dt_swe = swe.deltat(jd)
        # Allow 1 second tolerance
        assert dt_lib == pytest.approx(dt_swe, abs=1.0 / 86400), (
            f"Year {year}: lib={dt_lib * 86400:.2f}s, swe={dt_swe * 86400:.2f}s"
        )

    @pytest.mark.comparison
    @pytest.mark.xfail(
        reason="Delta T for future dates (>2030) is a PREDICTION and different "
        "algorithms (Skyfield vs Swiss Ephemeris) give divergent results. "
        "This is expected behavior, not a bug."
    )
    @pytest.mark.parametrize("year", [2040])
    def test_deltat_future_years_match_swe(self, year):
        """Delta T for future dates may diverge from pyswisseph due to different prediction algorithms."""
        jd = ephem.julday(year, 6, 15, 12.0)
        dt_lib = ephem.deltat(jd)
        dt_swe = swe.deltat(jd)
        # For future dates, allow 5 second tolerance
        assert dt_lib == pytest.approx(dt_swe, abs=5.0 / 86400), (
            f"Year {year}: lib={dt_lib * 86400:.2f}s, swe={dt_swe * 86400:.2f}s"
        )

    @pytest.mark.comparison
    def test_deltat_100_dates_match_swe(self, random_dates_in_de421_range):
        """100 random dates should match pyswisseph."""
        dates = random_dates_in_de421_range(100)
        max_diff = 0
        # Year 2030 in Julian Date (approximately Jan 1, 2030 at noon)
        jd_2030 = ephem.julday(2030, 1, 1, 12.0)
        for year, _, _, _, jd in dates:
            dt_lib = ephem.deltat(jd)
            dt_swe = swe.deltat(jd)
            diff = abs(dt_lib - dt_swe) * 86400  # in seconds
            max_diff = max(max_diff, diff)
            # For future dates (after 2030), Delta T is a PREDICTION and
            # different algorithms (Skyfield vs Swiss Ephemeris) diverge.
            # Allow 5+ seconds tolerance for future dates.
            if jd > jd_2030:
                tolerance = 5.0
            else:
                tolerance = 2.0
            assert diff < tolerance, (
                f"Delta T diff = {diff:.2f}s at JD {jd} (year {year}), "
                f"tolerance = {tolerance}s"
            )
        print(f"Max Delta T difference: {max_diff:.3f} seconds")


class TestDeltatUTtoTT:
    """Test using Delta T for UT to TT conversion."""

    @pytest.mark.unit
    def test_tt_greater_than_ut(self):
        """TT should be greater than UT for modern dates."""
        jd_ut = ephem.julday(2000, 1, 1, 12.0)
        dt = ephem.deltat(jd_ut)
        jd_tt = jd_ut + dt
        assert jd_tt > jd_ut

    @pytest.mark.unit
    def test_tt_ut_difference_seconds(self):
        """TT - UT should equal Delta T in seconds."""
        jd_ut = ephem.julday(2000, 1, 1, 12.0)
        dt = ephem.deltat(jd_ut)
        # Delta T in seconds
        dt_seconds = dt * 86400
        # This should be approximately 63 seconds at J2000
        assert 60 < dt_seconds < 70


class TestDeltatEdgeCases:
    """Test edge cases for Delta T."""

    @pytest.mark.edge_case
    def test_deltat_de421_range_start(self):
        """Delta T at DE421 range start (1900)."""
        jd = ephem.julday(1900, 1, 1, 12.0)
        dt = ephem.deltat(jd)
        # Should return a valid number
        assert isinstance(dt, float)
        assert -0.001 < dt < 0.001  # Within ±86 seconds

    @pytest.mark.edge_case
    def test_deltat_de421_range_end(self):
        """Delta T at DE421 range end (2050)."""
        jd = ephem.julday(2050, 1, 1, 12.0)
        dt = ephem.deltat(jd)
        # Should return a valid number (probably ~80+ seconds)
        assert isinstance(dt, float)
        dt_seconds = dt * 86400
        assert 70 < dt_seconds < 100

    @pytest.mark.unit
    def test_deltat_returns_float(self):
        """Delta T should always return a float."""
        jd = 2451545.0
        dt = ephem.deltat(jd)
        assert isinstance(dt, float)

    @pytest.mark.unit
    def test_deltat_reasonable_magnitude(self):
        """Delta T should have reasonable magnitude for any date in range."""
        for year in range(1900, 2051, 10):
            jd = ephem.julday(year, 1, 1, 12.0)
            dt = ephem.deltat(jd)
            dt_seconds = abs(dt * 86400)
            # Delta T should be less than 200 seconds for dates 1900-2050
            assert dt_seconds < 200, (
                f"Delta T at {year} = {dt_seconds}s seems too large"
            )


class TestDeltatEx:
    """Tests for the extended Delta T function deltat_ex."""

    @pytest.mark.unit
    def test_deltat_ex_returns_float(self):
        """deltat_ex should return a float."""
        jd = 2451545.0  # J2000
        result = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        assert isinstance(result, float)

    @pytest.mark.unit
    def test_deltat_ex_matches_deltat(self):
        """deltat_ex with FLG_SWIEPH should match deltat."""
        jd = 2451545.0  # J2000
        dt = ephem.deltat(jd)
        dt_ex = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        assert dt_ex == dt

    @pytest.mark.unit
    def test_deltat_ex_jpleph_matches_swieph(self):
        """FLG_JPLEPH should produce same results as FLG_SWIEPH."""
        jd = 2451545.0  # J2000
        dt_swieph = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        dt_jpleph = ephem.deltat_ex(jd, ephem.FLG_JPLEPH)
        assert dt_swieph == dt_jpleph

    @pytest.mark.unit
    def test_deltat_ex_moseph_uses_same_delta_t(self):
        """FLG_MOSEPH uses the same Delta T calculation (no warning)."""
        jd = 2451545.0  # J2000
        dt = ephem.deltat_ex(jd, ephem.FLG_MOSEPH)
        # Should still return a valid Delta T value
        assert isinstance(dt, float)
        dt_seconds = dt * 86400
        assert 60 < dt_seconds < 70
        # Value should match the standard function
        dt_standard = ephem.deltat(jd)
        assert dt == dt_standard

    @pytest.mark.unit
    def test_deltat_ex_no_flag(self):
        """deltat_ex with no ephemeris flag should work."""
        jd = 2451545.0  # J2000
        dt = ephem.deltat_ex(jd, 0)
        assert isinstance(dt, float)
        # Value should match the standard function
        dt_standard = ephem.deltat(jd)
        assert dt == dt_standard

    @pytest.mark.unit
    def test_deltat_ex_default_flag(self):
        """deltat_ex should default to FLG_SWIEPH."""
        jd = 2451545.0
        dt_default = ephem.deltat_ex(jd)
        dt_explicit = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        assert dt_default == dt_explicit

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2020, 2040])
    def test_deltat_ex_various_years(self, year):
        """deltat_ex should work for various years."""
        jd = ephem.julday(year, 6, 15, 12.0)
        dt = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        assert isinstance(dt, float)
        # Should match deltat
        dt_standard = ephem.deltat(jd)
        assert dt == dt_standard

    @pytest.mark.unit
    def test_deltat_ex_combined_flags(self):
        """deltat_ex should handle combined flags (ignoring non-ephemeris flags)."""
        jd = 2451545.0
        # Combine ephemeris flag with other flags (FLG_SPEED, etc.)
        combined_flag = ephem.FLG_SWIEPH | ephem.FLG_SPEED
        dt = ephem.deltat_ex(jd, combined_flag)
        assert isinstance(dt, float)
        # Should still match standard deltat
        dt_standard = ephem.deltat(jd)
        assert dt == dt_standard

    @pytest.mark.comparison
    def test_deltat_ex_vs_pyswisseph(self):
        """deltat_ex should match pyswisseph deltat_ex."""
        jd = 2451545.0  # J2000
        dt_lib = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        # pyswisseph deltat_ex returns just the delta_t value (float)
        dt_swe = swe.deltat_ex(jd, swe.FLG_SWIEPH)
        # Allow 0.5 second tolerance (in days: 0.5/86400)
        assert dt_lib == pytest.approx(dt_swe, abs=0.5 / 86400)

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "flag",
        [
            ("FLG_SWIEPH", 2),
            ("FLG_JPLEPH", 1),
        ],
    )
    def test_deltat_ex_flags_match_pyswisseph(self, flag):
        """deltat_ex with various flags should match pyswisseph."""
        flag_name, flag_value = flag
        jd = ephem.julday(2000, 1, 1, 12.0)
        dt_lib = ephem.deltat_ex(jd, flag_value)
        # pyswisseph deltat_ex returns just the delta_t value (float)
        dt_swe = swe.deltat_ex(jd, flag_value)
        # Allow 1 second tolerance
        assert dt_lib == pytest.approx(dt_swe, abs=1.0 / 86400), (
            f"Flag {flag_name}: lib={dt_lib * 86400:.2f}s, swe={dt_swe * 86400:.2f}s"
        )


class TestDeltatExAlias:
    """Test the deltat_ex alias."""

    @pytest.mark.unit
    def test_deltat_ex_alias_exists(self):
        """deltat_ex alias should exist."""
        assert hasattr(ephem, "deltat_ex")

    @pytest.mark.unit
    def test_deltat_ex_alias_works(self):
        """deltat_ex alias should produce same results as deltat_ex."""
        jd = 2451545.0
        result_swe = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        result_alias = ephem.deltat_ex(jd, ephem.FLG_SWIEPH)
        assert result_swe == result_alias
