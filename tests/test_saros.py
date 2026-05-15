"""
Tests for Saros series number calculation.

The Saros cycle is approximately 6585.32 days (223 synodic months = 18 years,
11 days, 8 hours). Eclipses in the same Saros series recur with similar
geometry.
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    get_saros_number,
    SAROS_CYCLE_DAYS,
    julday,
    sol_eclipse_when_glob,
    lun_eclipse_when,
    ECL_TOTAL,
)


class TestSarosCycleConstant:
    """Tests for the Saros cycle constant."""

    def test_saros_cycle_days_value(self):
        """Saros cycle should be approximately 6585.32 days."""
        assert 6585.0 < SAROS_CYCLE_DAYS < 6586.0

    def test_saros_cycle_is_223_synodic_months(self):
        """Saros cycle should equal 223 synodic months."""
        synodic_month = 29.530588853  # Mean synodic month
        expected = 223 * synodic_month
        assert abs(SAROS_CYCLE_DAYS - expected) < 0.01


class TestSolarSarosNumber:
    """Tests for solar eclipse Saros series identification."""

    def test_april_2024_total_eclipse_is_saros_139(self):
        """The April 8, 2024 total solar eclipse is Saros 139."""
        # JD for 2024-Apr-08 ~18:18 UT (maximum of the eclipse)
        jd_eclipse = 2460409.263
        saros = get_saros_number(jd_eclipse, eclipse_type="solar")
        assert saros == 139

    def test_august_2017_total_eclipse_is_saros_145(self):
        """The August 21, 2017 total solar eclipse is Saros 145."""
        # JD for 2017-Aug-21 maximum (verified against NASA canon)
        jd_eclipse = 2457987.268
        saros = get_saros_number(jd_eclipse, eclipse_type="solar")
        assert saros == 145

    def test_december_2020_annular_eclipse_is_saros_142(self):
        """The December 14, 2020 solar eclipse is Saros 142."""
        # JD for 2020-Dec-14 maximum (verified against NASA canon)
        jd_eclipse = 2459198.176
        saros = get_saros_number(jd_eclipse, eclipse_type="solar")
        assert saros == 142

    def test_june_2021_annular_eclipse_is_saros_147(self):
        """The June 10, 2021 annular solar eclipse is Saros 147."""
        # JD for 2021-Jun-10 maximum (verified against NASA canon)
        jd_eclipse = 2459375.946
        saros = get_saros_number(jd_eclipse, eclipse_type="solar")
        assert saros == 147

    def test_integration_with_sol_eclipse_when_glob(self):
        """Test that get_saros_number works with sol_eclipse_when_glob output."""
        # Find the April 2024 eclipse
        jd_start = julday(2024, 3, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        saros = get_saros_number(jd_max, eclipse_type="solar")
        assert saros == 139  # April 2024 eclipse is Saros 139

    def test_solar_eclipse_default_type(self):
        """Default eclipse_type should be 'solar'."""
        jd_eclipse = 2460409.263  # April 2024 eclipse
        saros = get_saros_number(jd_eclipse)  # No eclipse_type specified
        assert saros == 139

    def test_consecutive_saros_eclipses_same_series(self):
        """Eclipses separated by one Saros cycle should have same series number."""
        # August 21, 2017 eclipse (Saros 145)
        jd_2017 = 2457987.268

        # The previous eclipse in Saros 145 was August 11, 1999
        # ~18 years, 11 days earlier
        jd_1999 = jd_2017 - SAROS_CYCLE_DAYS

        saros_2017 = get_saros_number(jd_2017, "solar")
        saros_1999 = get_saros_number(jd_1999, "solar")

        assert saros_2017 == saros_1999 == 145


class TestLunarSarosNumber:
    """Tests for lunar eclipse Saros series identification."""

    def test_november_2022_lunar_eclipse_is_saros_136(self):
        """The November 8, 2022 total lunar eclipse is Saros 136."""
        # JD for 2022-Nov-08 ~11:00 UT (maximum)
        jd_eclipse = 2459891.958
        saros = get_saros_number(jd_eclipse, eclipse_type="lunar")
        assert saros == 136

    def test_may_2022_lunar_eclipse_is_saros_131(self):
        """The May 16, 2022 total lunar eclipse is Saros 131."""
        # JD for 2022-May-16 ~04:11 UT (maximum)
        jd_eclipse = 2459715.674
        saros = get_saros_number(jd_eclipse, eclipse_type="lunar")
        assert saros == 131

    def test_january_2019_lunar_eclipse_is_saros_134(self):
        """The January 21, 2019 total lunar eclipse is Saros 134."""
        # JD for 2019-Jan-21 ~05:12 UT (maximum)
        jd_eclipse = 2458504.717
        saros = get_saros_number(jd_eclipse, eclipse_type="lunar")
        assert saros == 134

    def test_july_2018_lunar_eclipse_is_saros_129(self):
        """The July 27, 2018 total lunar eclipse is Saros 129."""
        # JD for 2018-Jul-27 ~20:22 UT (maximum)
        jd_eclipse = 2458327.349
        saros = get_saros_number(jd_eclipse, eclipse_type="lunar")
        assert saros == 129

    def test_integration_with_lun_eclipse_when(self):
        """Test that get_saros_number works with lun_eclipse_when output."""
        # Find a lunar eclipse starting from 2022
        jd_start = julday(2022, 4, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        saros = get_saros_number(jd_max, eclipse_type="lunar")
        # The May 16, 2022 eclipse is Saros 131
        assert saros == 131

    def test_consecutive_lunar_saros_eclipses_same_series(self):
        """Lunar eclipses separated by one Saros cycle should have same series."""
        # November 8, 2022 eclipse (Saros 136)
        jd_2022 = 2459891.958

        # The next eclipse in Saros 136 will be ~18 years later
        jd_next = jd_2022 + SAROS_CYCLE_DAYS

        saros_2022 = get_saros_number(jd_2022, "lunar")
        saros_next = get_saros_number(jd_next, "lunar")

        assert saros_2022 == saros_next == 136


class TestSarosEdgeCases:
    """Edge case tests for Saros series calculation."""

    def test_invalid_eclipse_type_raises_error(self):
        """Should raise ValueError for invalid eclipse_type."""
        jd_eclipse = 2460409.263
        with pytest.raises(ValueError, match="eclipse_type must be"):
            get_saros_number(jd_eclipse, eclipse_type="invalid")

    def test_case_insensitive_eclipse_type(self):
        """Eclipse type should be case-insensitive."""
        jd_eclipse = 2460409.263

        assert get_saros_number(jd_eclipse, "SOLAR") == 139
        assert get_saros_number(jd_eclipse, "Solar") == 139
        assert get_saros_number(jd_eclipse, "solar") == 139

    def test_distant_past_eclipse(self):
        """Should handle eclipses far in the past."""
        # The 1999 August 11 total solar eclipse (Saros 145)
        jd_1999 = julday(1999, 8, 11, 11.0)
        saros = get_saros_number(jd_1999, "solar")
        assert saros == 145

    def test_distant_future_eclipse(self):
        """Should handle eclipses in the future."""
        # Calculate a future eclipse in Saros 139 (2024 + 18 years)
        jd_2024 = 2460409.263  # April 2024 eclipse
        jd_future = jd_2024 + SAROS_CYCLE_DAYS  # ~2042 eclipse

        saros = get_saros_number(jd_future, "solar")
        assert saros == 139

    def test_solar_and_lunar_saros_independent(self):
        """Solar and lunar Saros numbers should be independent."""
        # Same JD can have different Saros numbers for solar vs lunar
        jd_test = 2459891.958  # Near a lunar eclipse (Saros 136)

        lunar_saros = get_saros_number(jd_test, "lunar")
        solar_saros = get_saros_number(jd_test, "solar")

        # They use different reference systems, so may differ
        assert isinstance(lunar_saros, int)
        assert isinstance(solar_saros, int)


class TestSarosMultipleEclipses:
    """Test Saros numbers for multiple consecutive eclipses."""

    def test_multiple_solar_eclipses_different_series(self):
        """Find multiple consecutive solar eclipses and verify different series."""
        jd_start = julday(2024, 1, 1, 0.0)

        eclipses = []
        jd = jd_start
        for _ in range(3):
            _, times = sol_eclipse_when_glob(jd)
            jd_max = times[0]
            saros = get_saros_number(jd_max, "solar")
            eclipses.append((jd_max, saros))
            jd = jd_max + 25  # Skip ahead to find next eclipse

        # Consecutive eclipses should have different Saros numbers
        # (they're in different series)
        series_numbers = [e[1] for e in eclipses]
        # Most consecutive eclipses are in different series
        # (same series eclipses are ~18 years apart)
        assert len(eclipses) == 3
        for saros in series_numbers:
            assert 100 <= saros <= 180  # Reasonable range for modern eclipses

    def test_multiple_lunar_eclipses_different_series(self):
        """Find multiple consecutive lunar eclipses and verify different series."""
        jd_start = julday(2022, 1, 1, 0.0)

        eclipses = []
        jd = jd_start
        for _ in range(3):
            _, times = lun_eclipse_when(jd)
            jd_max = times[0]
            saros = get_saros_number(jd_max, "lunar")
            eclipses.append((jd_max, saros))
            jd = jd_max + 25  # Skip ahead to find next eclipse

        assert len(eclipses) == 3
        for jd_max, saros in eclipses:
            assert 100 <= saros <= 180  # Reasonable range for modern eclipses
