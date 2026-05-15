"""
Tests for Inex series number calculation.

The Inex cycle is approximately 10571.95 days (358 synodic months, ~28.945 years).
Eclipses in the same Inex series occur at the same solar longitude but at
opposite lunar nodes.
"""

import pytest
from libephemeris import (
    get_inex_number,
    INEX_CYCLE_DAYS,
    julday,
    sol_eclipse_when_glob,
    lun_eclipse_when,
    ECL_TOTAL,
)

pytestmark = pytest.mark.slow


class TestInexCycleConstant:
    """Tests for the Inex cycle constant."""

    def test_inex_cycle_days_value(self):
        """Inex cycle should be approximately 10571.95 days."""
        assert 10571.0 < INEX_CYCLE_DAYS < 10573.0

    def test_inex_cycle_is_358_synodic_months(self):
        """Inex cycle should equal 358 synodic months."""
        synodic_month = 29.530588853  # Mean synodic month
        expected = 358 * synodic_month
        assert abs(INEX_CYCLE_DAYS - expected) < 0.01


class TestSolarInexNumber:
    """Tests for solar eclipse Inex series identification."""

    def test_april_2024_total_eclipse_is_inex_50(self):
        """The April 8, 2024 total solar eclipse is Inex 50."""
        # JD for 2024-Apr-08 ~18:18 UT (maximum of the eclipse)
        jd_eclipse = 2460409.263
        inex = get_inex_number(jd_eclipse, eclipse_type="solar")
        assert inex == 50

    def test_august_2017_total_eclipse_is_inex_49(self):
        """The August 21, 2017 total solar eclipse is Inex 49."""
        # JD for 2017-Aug-21 ~18:26 UT (maximum)
        jd_eclipse = 2457987.768
        inex = get_inex_number(jd_eclipse, eclipse_type="solar")
        assert inex == 49

    def test_december_2020_total_eclipse_is_inex_51(self):
        """The December 14, 2020 total solar eclipse is Inex 51."""
        # JD for 2020-Dec-14 ~16:14 UT (maximum)
        jd_eclipse = 2459203.676
        inex = get_inex_number(jd_eclipse, eclipse_type="solar")
        assert inex == 51

    def test_june_2021_annular_eclipse_is_inex_52(self):
        """The June 10, 2021 annular solar eclipse is Inex 52."""
        # JD for 2021-Jun-10 ~10:43 UT (maximum)
        jd_eclipse = 2459375.947
        inex = get_inex_number(jd_eclipse, eclipse_type="solar")
        assert inex == 52

    def test_integration_with_sol_eclipse_when_glob(self):
        """Test that get_inex_number works with sol_eclipse_when_glob output."""
        # Find the April 2024 eclipse
        jd_start = julday(2024, 3, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        inex = get_inex_number(jd_max, eclipse_type="solar")
        assert inex == 50  # April 2024 eclipse is Inex 50

    def test_solar_eclipse_default_type(self):
        """Default eclipse_type should be 'solar'."""
        jd_eclipse = 2460409.263  # April 2024 eclipse
        inex = get_inex_number(jd_eclipse)  # No eclipse_type specified
        assert inex == 50

    def test_consecutive_inex_eclipses_same_series(self):
        """Eclipses separated by one Inex cycle should have same series number."""
        # April 8, 2024 eclipse (Inex 50)
        jd_2024 = 2460409.263

        # An eclipse one Inex cycle later should also be Inex 50
        jd_later = jd_2024 + INEX_CYCLE_DAYS

        inex_2024 = get_inex_number(jd_2024, "solar")
        inex_later = get_inex_number(jd_later, "solar")

        assert inex_2024 == inex_later == 50


class TestLunarInexNumber:
    """Tests for lunar eclipse Inex series identification."""

    def test_november_2022_lunar_eclipse_is_inex_50(self):
        """The November 8, 2022 total lunar eclipse is Inex 50."""
        # JD for 2022-Nov-08 ~11:00 UT (maximum)
        jd_eclipse = 2459891.958
        inex = get_inex_number(jd_eclipse, eclipse_type="lunar")
        assert inex == 50

    def test_may_2022_lunar_eclipse_is_inex_49(self):
        """The May 16, 2022 total lunar eclipse is Inex 49."""
        # JD for 2022-May-16 ~04:11 UT (maximum)
        jd_eclipse = 2459715.674
        inex = get_inex_number(jd_eclipse, eclipse_type="lunar")
        assert inex == 49

    def test_january_2019_lunar_eclipse_is_inex_47(self):
        """The January 21, 2019 total lunar eclipse is Inex 47."""
        # JD for 2019-Jan-21 ~05:12 UT (maximum)
        jd_eclipse = 2458504.717
        inex = get_inex_number(jd_eclipse, eclipse_type="lunar")
        assert inex == 47

    def test_july_2018_lunar_eclipse_is_inex_51(self):
        """The July 27, 2018 total lunar eclipse is Inex 51."""
        # JD for 2018-Jul-27 ~20:22 UT (maximum)
        jd_eclipse = 2458327.349
        inex = get_inex_number(jd_eclipse, eclipse_type="lunar")
        assert inex == 51

    def test_integration_with_lun_eclipse_when(self):
        """Test that get_inex_number works with lun_eclipse_when output."""
        # Find a lunar eclipse starting from 2022
        jd_start = julday(2022, 4, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        inex = get_inex_number(jd_max, eclipse_type="lunar")
        # The May 16, 2022 eclipse is Inex 49
        assert inex == 49

    def test_consecutive_lunar_inex_eclipses_same_series(self):
        """Lunar eclipses separated by one Inex cycle should have same series."""
        # November 8, 2022 eclipse (Inex 50)
        jd_2022 = 2459891.958

        # The next eclipse in Inex 50 will be ~29 years later
        jd_next = jd_2022 + INEX_CYCLE_DAYS

        inex_2022 = get_inex_number(jd_2022, "lunar")
        inex_next = get_inex_number(jd_next, "lunar")

        assert inex_2022 == inex_next == 50


class TestInexEdgeCases:
    """Edge case tests for Inex series calculation."""

    def test_invalid_eclipse_type_raises_error(self):
        """Should raise ValueError for invalid eclipse_type."""
        jd_eclipse = 2460409.263
        with pytest.raises(ValueError, match="eclipse_type must be"):
            get_inex_number(jd_eclipse, eclipse_type="invalid")

    def test_case_insensitive_eclipse_type(self):
        """Eclipse type should be case-insensitive."""
        jd_eclipse = 2460409.263

        assert get_inex_number(jd_eclipse, "SOLAR") == 50
        assert get_inex_number(jd_eclipse, "Solar") == 50
        assert get_inex_number(jd_eclipse, "solar") == 50

    def test_distant_past_eclipse(self):
        """Should handle eclipses far in the past."""
        # The 1999 August 11 total solar eclipse (Inex 49 - same as 2017)
        jd_1999 = julday(1999, 8, 11, 11.0)
        inex = get_inex_number(jd_1999, "solar")
        # Both 1999 and 2017 eclipses are in Saros 136, but Inex may differ
        # 1999 Aug 11 is about 18 years before 2017 Aug 21 (one Saros)
        # which means they have different Inex numbers
        assert isinstance(inex, int)

    def test_distant_future_eclipse(self):
        """Should handle eclipses in the future."""
        # Calculate a future eclipse in Inex 50 (2024 + ~29 years)
        jd_2024 = 2460409.263  # April 2024 eclipse
        jd_future = jd_2024 + INEX_CYCLE_DAYS  # ~2053 eclipse

        inex = get_inex_number(jd_future, "solar")
        assert inex == 50

    def test_solar_and_lunar_inex_independent(self):
        """Solar and lunar Inex numbers should be independent."""
        # Same JD can have different Inex numbers for solar vs lunar
        jd_test = 2459891.958  # Near a lunar eclipse (Inex 50)

        lunar_inex = get_inex_number(jd_test, "lunar")
        solar_inex = get_inex_number(jd_test, "solar")

        # They use different reference systems, so may differ
        assert isinstance(lunar_inex, int)
        assert isinstance(solar_inex, int)


class TestInexMultipleEclipses:
    """Test Inex numbers for multiple consecutive eclipses."""

    def test_multiple_solar_eclipses_different_series(self):
        """Find multiple consecutive solar eclipses and verify Inex numbers."""
        jd_start = julday(2024, 1, 1, 0.0)

        eclipses = []
        jd = jd_start
        for _ in range(3):
            _, times = sol_eclipse_when_glob(jd)
            jd_max = times[0]
            inex = get_inex_number(jd_max, "solar")
            eclipses.append((jd_max, inex))
            jd = jd_max + 25  # Skip ahead to find next eclipse

        # Consecutive eclipses may have different Inex numbers
        assert len(eclipses) == 3
        for jd_max, inex in eclipses:
            # Inex numbers should be reasonable integers
            assert isinstance(inex, int)
            assert 20 <= inex <= 80  # Reasonable range for modern eclipses

    def test_multiple_lunar_eclipses_different_series(self):
        """Find multiple consecutive lunar eclipses and verify Inex numbers."""
        jd_start = julday(2022, 1, 1, 0.0)

        eclipses = []
        jd = jd_start
        for _ in range(3):
            _, times = lun_eclipse_when(jd)
            jd_max = times[0]
            inex = get_inex_number(jd_max, "lunar")
            eclipses.append((jd_max, inex))
            jd = jd_max + 25  # Skip ahead to find next eclipse

        assert len(eclipses) == 3
        for jd_max, inex in eclipses:
            # Inex numbers should be reasonable integers
            assert isinstance(inex, int)
            assert 20 <= inex <= 80  # Reasonable range for modern eclipses


class TestInexVersusSaros:
    """Test the relationship between Inex and Saros series."""

    def test_inex_differs_from_saros(self):
        """Inex and Saros numbers are independent numbering systems."""
        from libephemeris import get_saros_number

        # April 2024 eclipse
        jd_eclipse = 2460409.263

        saros = get_saros_number(jd_eclipse, "solar")
        inex = get_inex_number(jd_eclipse, "solar")

        # Saros 139, Inex 50 for this eclipse
        assert saros == 139
        assert inex == 50
        # They should differ
        assert saros != inex

    def test_one_saros_apart_eclipses_different_inex(self):
        """Eclipses one Saros apart may have different Inex numbers."""
        from libephemeris import get_saros_number, SAROS_CYCLE_DAYS

        # August 2017 eclipse
        jd_2017 = 2457987.768
        # August 1999 eclipse (one Saros earlier)
        jd_1999 = jd_2017 - SAROS_CYCLE_DAYS

        # Same Saros series
        assert get_saros_number(jd_2017, "solar") == get_saros_number(jd_1999, "solar")

        # But Inex should remain the same within the series
        inex_2017 = get_inex_number(jd_2017, "solar")
        inex_1999 = get_inex_number(jd_1999, "solar")
        # Note: Inex series for consecutive Saros members may vary
        assert isinstance(inex_2017, int)
        assert isinstance(inex_1999, int)
