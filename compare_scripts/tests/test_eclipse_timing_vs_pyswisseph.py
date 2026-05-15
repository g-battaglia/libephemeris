"""
Tests validating eclipse timing against pyswisseph.

Compares eclipse maximum times from libephemeris against pyswisseph for
a list of known eclipses. Verifies that timing agrees within 10 seconds.

Known eclipses tested:
- Total solar eclipses: 2017-08-21, 2024-04-08
- Lunar eclipses: 2018-01-31, 2022-11-08

Reference: Swiss Ephemeris (pyswisseph) is the authoritative source for
high-precision eclipse calculations.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_ALLTYPES_SOLAR,
    ECL_ALLTYPES_LUNAR,
)


# Tolerance for eclipse timing comparison (seconds)
TIMING_TOLERANCE_SECONDS = 10.0


def time_diff_seconds(jd1: float, jd2: float) -> float:
    """Calculate absolute time difference between two Julian Days in seconds."""
    return abs(jd1 - jd2) * 86400.0


class TestSolarEclipseTimingVsPyswisseph:
    """Compare solar eclipse timing between libephemeris and pyswisseph."""

    def test_august_2017_total_solar_eclipse(self):
        """Test August 21, 2017 total solar eclipse ('Great American Eclipse').

        This was a major total solar eclipse crossing the continental United States.
        Reference: NASA Eclipse Bulletin
        """
        # Start search from August 1, 2017
        jd_start = swe.julday(2017, 8, 1, 0.0)

        # Get eclipse time from pyswisseph
        ret_swe = swe.sol_eclipse_when_glob(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max_swe = ret_swe[1][0]

        # Get eclipse time from libephemeris
        # Returns (ecl_type, times) tuple
        ecl_type_py, times_py = ephem.sol_eclipse_when_glob(
            jd_start, FLG_SWIEPH, ECL_TOTAL
        )
        jd_max_py = times_py[0]

        # Verify we found the same eclipse (August 21, 2017)
        # pyswisseph jd should correspond to Aug 21, 2017
        year_swe, month_swe, day_swe, hour_swe = swe.revjul(jd_max_swe)
        assert year_swe == 2017, f"Expected year 2017, got {year_swe}"
        assert month_swe == 8, f"Expected month 8, got {month_swe}"
        assert day_swe == 21, f"Expected day 21, got {day_swe}"

        # Calculate time difference
        diff_seconds = time_diff_seconds(jd_max_py, jd_max_swe)

        assert diff_seconds < TIMING_TOLERANCE_SECONDS, (
            f"August 2017 total solar eclipse timing difference: {diff_seconds:.2f} seconds. "
            f"Expected < {TIMING_TOLERANCE_SECONDS} seconds. "
            f"pyswisseph JD: {jd_max_swe:.8f}, libephemeris JD: {jd_max_py:.8f}"
        )

    def test_april_2024_total_solar_eclipse(self):
        """Test April 8, 2024 total solar eclipse (North American eclipse).

        This was a highly anticipated total solar eclipse crossing Mexico,
        United States, and Canada.
        Reference: NASA Eclipse Bulletin
        """
        # Start search from April 1, 2024
        jd_start = swe.julday(2024, 4, 1, 0.0)

        # Get eclipse time from pyswisseph
        ret_swe = swe.sol_eclipse_when_glob(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max_swe = ret_swe[1][0]

        # Get eclipse time from libephemeris
        # Returns (ecl_type, times) tuple
        ecl_type_py, times_py = ephem.sol_eclipse_when_glob(
            jd_start, FLG_SWIEPH, ECL_TOTAL
        )
        jd_max_py = times_py[0]

        # Verify we found the same eclipse (April 8, 2024)
        year_swe, month_swe, day_swe, hour_swe = swe.revjul(jd_max_swe)
        assert year_swe == 2024, f"Expected year 2024, got {year_swe}"
        assert month_swe == 4, f"Expected month 4, got {month_swe}"
        assert day_swe == 8, f"Expected day 8, got {day_swe}"

        # Calculate time difference
        diff_seconds = time_diff_seconds(jd_max_py, jd_max_swe)

        assert diff_seconds < TIMING_TOLERANCE_SECONDS, (
            f"April 2024 total solar eclipse timing difference: {diff_seconds:.2f} seconds. "
            f"Expected < {TIMING_TOLERANCE_SECONDS} seconds. "
            f"pyswisseph JD: {jd_max_swe:.8f}, libephemeris JD: {jd_max_py:.8f}"
        )

    def test_solar_eclipse_any_type_consistency(self):
        """Test that both libraries find the same eclipses when searching for any type."""
        # Search for the next 3 solar eclipses starting from 2017
        jd = swe.julday(2017, 1, 1, 0.0)

        for i in range(3):
            # Get eclipse from pyswisseph
            ret_swe = swe.sol_eclipse_when_glob(jd, FLG_SWIEPH, ECL_ALLTYPES_SOLAR)
            jd_max_swe = ret_swe[1][0]

            # Get eclipse from libephemeris (returns (ecl_type, times))
            _, times_py = ephem.sol_eclipse_when_glob(jd, FLG_SWIEPH, 0)
            jd_max_py = times_py[0]

            # Calculate time difference
            diff_seconds = time_diff_seconds(jd_max_py, jd_max_swe)

            assert diff_seconds < TIMING_TOLERANCE_SECONDS, (
                f"Solar eclipse #{i + 1} timing difference: {diff_seconds:.2f} seconds. "
                f"Expected < {TIMING_TOLERANCE_SECONDS} seconds."
            )

            # Move to next eclipse
            jd = jd_max_swe + 30  # Skip ahead to find next


class TestLunarEclipseTimingVsPyswisseph:
    """Compare lunar eclipse timing between libephemeris and pyswisseph."""

    def test_january_2018_lunar_eclipse(self):
        """Test January 31, 2018 total lunar eclipse ('Super Blue Blood Moon').

        This was a total lunar eclipse visible from Asia, Australia, and
        western North America. It was also a supermoon and the second
        full moon of January 2018 (blue moon).
        Reference: NASA Eclipse Bulletin
        """
        # Start search from January 1, 2018
        jd_start = swe.julday(2018, 1, 1, 0.0)

        # Get eclipse time from pyswisseph
        ret_swe = swe.lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max_swe = ret_swe[1][0]

        # Get eclipse time from libephemeris (returns (ecl_type, times))
        ecl_type_py, times_py = ephem.lun_eclipse_when(
            jd_start, FLG_SWIEPH, ECL_TOTAL
        )
        jd_max_py = times_py[0]

        # Verify we found the same eclipse (January 31, 2018)
        year_swe, month_swe, day_swe, hour_swe = swe.revjul(jd_max_swe)
        assert year_swe == 2018, f"Expected year 2018, got {year_swe}"
        assert month_swe == 1, f"Expected month 1, got {month_swe}"
        assert day_swe == 31, f"Expected day 31, got {day_swe}"

        # Calculate time difference
        diff_seconds = time_diff_seconds(jd_max_py, jd_max_swe)

        assert diff_seconds < TIMING_TOLERANCE_SECONDS, (
            f"January 2018 lunar eclipse timing difference: {diff_seconds:.2f} seconds. "
            f"Expected < {TIMING_TOLERANCE_SECONDS} seconds. "
            f"pyswisseph JD: {jd_max_swe:.8f}, libephemeris JD: {jd_max_py:.8f}"
        )

    def test_november_2022_lunar_eclipse(self):
        """Test November 8, 2022 total lunar eclipse.

        This total lunar eclipse was visible from North America, Pacific,
        Asia, and Australia.
        Reference: NASA Eclipse Bulletin
        """
        # Start search from November 1, 2022
        jd_start = swe.julday(2022, 11, 1, 0.0)

        # Get eclipse time from pyswisseph
        ret_swe = swe.lun_eclipse_when(jd_start, FLG_SWIEPH, ECL_TOTAL)
        jd_max_swe = ret_swe[1][0]

        # Get eclipse time from libephemeris (returns (ecl_type, times))
        ecl_type_py, times_py = ephem.lun_eclipse_when(
            jd_start, FLG_SWIEPH, ECL_TOTAL
        )
        jd_max_py = times_py[0]

        # Verify we found the same eclipse (November 8, 2022)
        year_swe, month_swe, day_swe, hour_swe = swe.revjul(jd_max_swe)
        assert year_swe == 2022, f"Expected year 2022, got {year_swe}"
        assert month_swe == 11, f"Expected month 11, got {month_swe}"
        assert day_swe == 8, f"Expected day 8, got {day_swe}"

        # Calculate time difference
        diff_seconds = time_diff_seconds(jd_max_py, jd_max_swe)

        assert diff_seconds < TIMING_TOLERANCE_SECONDS, (
            f"November 2022 lunar eclipse timing difference: {diff_seconds:.2f} seconds. "
            f"Expected < {TIMING_TOLERANCE_SECONDS} seconds. "
            f"pyswisseph JD: {jd_max_swe:.8f}, libephemeris JD: {jd_max_py:.8f}"
        )

    def test_lunar_eclipse_any_type_consistency(self):
        """Test that both libraries find the same lunar eclipses when searching for any type."""
        # Search for the next 3 lunar eclipses starting from 2018
        jd = swe.julday(2018, 1, 1, 0.0)

        for i in range(3):
            # Get eclipse from pyswisseph
            ret_swe = swe.lun_eclipse_when(jd, FLG_SWIEPH, ECL_ALLTYPES_LUNAR)
            jd_max_swe = ret_swe[1][0]

            # Get eclipse from libephemeris (returns (ecl_type, times))
            _, times_py = ephem.lun_eclipse_when(jd, FLG_SWIEPH, 0)
            jd_max_py = times_py[0]

            # Calculate time difference
            diff_seconds = time_diff_seconds(jd_max_py, jd_max_swe)

            assert diff_seconds < TIMING_TOLERANCE_SECONDS, (
                f"Lunar eclipse #{i + 1} timing difference: {diff_seconds:.2f} seconds. "
                f"Expected < {TIMING_TOLERANCE_SECONDS} seconds."
            )

            # Move to next eclipse
            jd = jd_max_swe + 30  # Skip ahead to find next


class TestEclipseTimingSummary:
    """Summary tests comparing multiple eclipses against pyswisseph."""

    @pytest.mark.parametrize(
        "year,month,day,eclipse_type,description",
        [
            (2017, 8, 1, ECL_TOTAL, "Aug 2017 Total Solar"),
            (2024, 4, 1, ECL_TOTAL, "Apr 2024 Total Solar"),
        ],
    )
    def test_solar_eclipses_parametrized(
        self, year, month, day, eclipse_type, description
    ):
        """Parametrized test for known solar eclipses."""
        jd_start = swe.julday(year, month, day, 0.0)

        # Get eclipse from pyswisseph
        ret_swe = swe.sol_eclipse_when_glob(jd_start, FLG_SWIEPH, eclipse_type)
        jd_max_swe = ret_swe[1][0]

        # Get eclipse from libephemeris (returns (ecl_type, times))
        _, times_py = ephem.sol_eclipse_when_glob(jd_start, FLG_SWIEPH, eclipse_type)
        jd_max_py = times_py[0]

        diff_seconds = time_diff_seconds(jd_max_py, jd_max_swe)

        assert diff_seconds < TIMING_TOLERANCE_SECONDS, (
            f"{description}: timing difference {diff_seconds:.2f}s exceeds tolerance"
        )

    @pytest.mark.parametrize(
        "year,month,day,eclipse_type,description",
        [
            (2018, 1, 1, ECL_TOTAL, "Jan 2018 Total Lunar"),
            (2022, 11, 1, ECL_TOTAL, "Nov 2022 Total Lunar"),
        ],
    )
    def test_lunar_eclipses_parametrized(
        self, year, month, day, eclipse_type, description
    ):
        """Parametrized test for known lunar eclipses."""
        jd_start = swe.julday(year, month, day, 0.0)

        # Get eclipse from pyswisseph
        ret_swe = swe.lun_eclipse_when(jd_start, FLG_SWIEPH, eclipse_type)
        jd_max_swe = ret_swe[1][0]

        # Get eclipse from libephemeris (returns (ecl_type, times))
        _, times_py = ephem.lun_eclipse_when(jd_start, FLG_SWIEPH, eclipse_type)
        jd_max_py = times_py[0]

        diff_seconds = time_diff_seconds(jd_max_py, jd_max_swe)

        assert diff_seconds < TIMING_TOLERANCE_SECONDS, (
            f"{description}: timing difference {diff_seconds:.2f}s exceeds tolerance"
        )
