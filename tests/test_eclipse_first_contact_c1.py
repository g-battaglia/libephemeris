"""
Tests for calc_eclipse_first_contact_c1 function in libephemeris.

Tests the calculation of first external contact (C1) for solar eclipses
using Besselian elements.

C1 is the moment when the Moon's disk first externally touches the Sun's disk,
marking the beginning of a solar eclipse. At this instant, the penumbral
shadow cone first touches Earth's surface.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
- Meeus "Astronomical Algorithms" Chapter 54
"""

import pytest
from libephemeris import (
    julday,
    calc_eclipse_first_contact_c1,
    sol_eclipse_when_glob,
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_ANNULAR,
)

pytestmark = pytest.mark.slow


class TestEclipseFirstContactC1BasicFunctionality:
    """Test basic functionality of calc_eclipse_first_contact_c1."""

    def test_function_exists(self):
        """Test that calc_eclipse_first_contact_c1 function exists and is callable."""
        from libephemeris.eclipse import calc_eclipse_first_contact_c1

        assert callable(calc_eclipse_first_contact_c1)

    def test_function_exported_from_main_module(self):
        """Test that function is exported from main libephemeris module."""
        from libephemeris import calc_eclipse_first_contact_c1

        assert callable(calc_eclipse_first_contact_c1)

    def test_returns_float(self):
        """Test that function returns a float value."""
        # Get a known eclipse maximum
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_eclipse_first_contact_c1(jd_max)

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts optional flags parameter."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        result = calc_eclipse_first_contact_c1(jd_max, flags=FLG_SWIEPH)

        assert isinstance(result, float)


class TestEclipseFirstContactC1TimingPrecision:
    """Test C1 timing precision against known reference data."""

    def test_april_2024_total_eclipse_c1(self):
        """Test April 8, 2024 total solar eclipse first contact timing.

        NASA Reference (from eclipse.gsfc.nasa.gov/SEpath/SEpath2001/SE2024Apr08Tpath.html):
        - Penumbral Eclipse Begins: 15:42:10.0 TD (Greatest eclipse: 18:17:18.3 UT)
        - Delta T ~69 seconds, so C1 ~ 15:41:01 UT

        Note: NASA times are in TD (Dynamical Time), we compute in UT.
        There is inherent uncertainty in the TD/UT conversion and reference data,
        so we allow 120 seconds tolerance.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_c1 = calc_eclipse_first_contact_c1(jd_max)

        # NASA reference: Penumbral begins ~15:42:10 TD = ~15:41:01 UT
        # 15 + 41/60 + 1/3600 = 15.6836 hours
        jd_reference = julday(2024, 4, 8, 15.6836)

        # Calculate time difference in seconds
        time_diff_days = abs(jd_c1 - jd_reference)
        time_diff_seconds = time_diff_days * 86400

        # Precision requirement: < 120 seconds (2 minute tolerance for reference data uncertainty)
        # This accounts for TD/UT conversion uncertainty and reference data precision
        assert time_diff_seconds < 120, (
            f"C1 timing error: {time_diff_seconds:.1f} seconds. "
            f"Expected < 120 seconds from reference. "
            f"Calculated JD: {jd_c1:.6f}, Reference JD: {jd_reference:.6f}"
        )

    def test_october_2023_annular_eclipse_c1(self):
        """Test October 14, 2023 annular solar eclipse first contact timing.

        This eclipse crossed North, Central, and South America.
        """
        jd_start = julday(2023, 9, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)
        jd_max = times[0]

        jd_c1 = calc_eclipse_first_contact_c1(jd_max)

        # C1 should be several hours before maximum
        # Maximum is around 18:00 UT, so C1 should be around 15:00-16:00 UT
        # (typical global eclipse duration is 4-6 hours total)
        time_before_max_hours = (jd_max - jd_c1) * 24

        assert 1.5 < time_before_max_hours < 3.5, (
            f"C1 should be 1.5-3.5 hours before maximum. "
            f"Got {time_before_max_hours:.2f} hours."
        )

    def test_december_2021_total_eclipse_c1(self):
        """Test December 4, 2021 total solar eclipse (Antarctica) first contact."""
        jd_start = julday(2021, 11, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_c1 = calc_eclipse_first_contact_c1(jd_max)

        # C1 should be before maximum
        assert jd_c1 < jd_max, "C1 should be before maximum eclipse"

        # And within reasonable range (< 3 hours before)
        time_before_max_hours = (jd_max - jd_c1) * 24
        assert 1.0 < time_before_max_hours < 3.0, (
            f"C1 should be 1-3 hours before maximum. "
            f"Got {time_before_max_hours:.2f} hours."
        )


class TestEclipseFirstContactC1ConsistencyWithSolEclipseWhenGlob:
    """Test that C1 calculation is consistent with sol_eclipse_when_glob results."""

    def test_c1_matches_sol_eclipse_when_glob_first_contact(self):
        """Test that our C1 matches the first contact from sol_eclipse_when_glob.

        sol_eclipse_when_glob returns contact times in its result array.
        Our dedicated C1 function should produce the same result.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]
        jd_first_glob = times[2]  # Eclipse begin from sol_eclipse_when_glob

        jd_c1 = calc_eclipse_first_contact_c1(jd_max)

        # Both should give the same first contact time (within numerical precision)
        time_diff_seconds = abs(jd_c1 - jd_first_glob) * 86400

        assert time_diff_seconds < 10, (
            f"C1 calculation differs from sol_eclipse_when_glob by {time_diff_seconds:.1f} seconds. "
            f"C1: {jd_c1:.6f}, sol_eclipse_when_glob: {jd_first_glob:.6f}"
        )

    def test_c1_consistency_across_multiple_eclipses(self):
        """Test C1 calculation consistency across multiple eclipses."""
        # Test several eclipses
        eclipse_starts = [
            (julday(2023, 9, 1, 0.0), ECL_ANNULAR),  # Oct 2023 annular
            (julday(2024, 1, 1, 0.0), ECL_TOTAL),  # Apr 2024 total
            (julday(2025, 1, 1, 0.0), ECL_TOTAL),  # Future eclipse
        ]

        for jd_start, ecl_type in eclipse_starts:
            _, times = sol_eclipse_when_glob(jd_start, ecltype=ecl_type)
            jd_max = times[0]
            jd_first_glob = times[2]

            jd_c1 = calc_eclipse_first_contact_c1(jd_max)

            time_diff_seconds = abs(jd_c1 - jd_first_glob) * 86400

            assert time_diff_seconds < 10, (
                f"C1 mismatch for eclipse at JD {jd_max:.2f}: "
                f"diff = {time_diff_seconds:.1f} seconds"
            )


class TestEclipseFirstContactC1PhysicalProperties:
    """Test physical properties and constraints of C1 calculations."""

    def test_c1_is_before_maximum(self):
        """Test that C1 is always before eclipse maximum."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_c1 = calc_eclipse_first_contact_c1(jd_max)

        assert jd_c1 < jd_max, (
            f"C1 should be before maximum. C1: {jd_c1:.6f}, Max: {jd_max:.6f}"
        )

    def test_c1_duration_is_physically_reasonable(self):
        """Test that time from C1 to maximum is physically reasonable.

        For a global solar eclipse, the time from first contact to maximum
        is typically 1.5-3 hours, depending on the eclipse geometry.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_c1 = calc_eclipse_first_contact_c1(jd_max)

        duration_hours = (jd_max - jd_c1) * 24

        # Physical constraints: penumbral contact to maximum typically 1.5-3 hours
        assert 1.0 < duration_hours < 4.0, (
            f"Duration from C1 to max ({duration_hours:.2f} hours) seems unreasonable. "
            f"Expected 1-4 hours for typical global eclipse geometry."
        )

    def test_c1_returns_nonzero_for_valid_eclipse(self):
        """Test that C1 returns a non-zero value for valid eclipses."""
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)
        jd_max = times[0]

        jd_c1 = calc_eclipse_first_contact_c1(jd_max)

        assert jd_c1 > 0, "C1 should be non-zero for a valid eclipse"
        assert jd_c1 > 2400000, "C1 should be a valid Julian Day (> 2400000)"


class TestEclipseFirstContactC1EdgeCases:
    """Test edge cases and boundary conditions."""

    def test_c1_with_partial_eclipse(self):
        """Test C1 calculation for a partial solar eclipse.

        Note: Pure partial-only solar eclipses are rare. This test uses any
        eclipse type to verify the function works with different eclipse geometries.
        """
        # Find any eclipse (partial eclipses are rare, so use ALLTYPES if needed)
        jd_start = julday(2025, 3, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start)  # Any type

        if times[0] > 0:  # If an eclipse was found
            jd_max = times[0]
            jd_c1 = calc_eclipse_first_contact_c1(jd_max)

            # C1 should work for any eclipse type
            assert jd_c1 > 0, "C1 should be calculable for all eclipse types"
            assert jd_c1 < jd_max, "C1 should be before maximum"

    def test_c1_symmetry_with_fourth_contact(self):
        """Test that C1 and C4 are roughly symmetric around maximum.

        For a global eclipse, the time from C1 to max should be roughly equal
        to the time from max to C4, unless the eclipse is asymmetric.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        _, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]
        jd_first = times[2]  # Eclipse begin (C1)
        jd_fourth = times[3]  # Eclipse end (C4)

        jd_c1 = calc_eclipse_first_contact_c1(jd_max)

        # Time from C1 to max
        dt_before = jd_max - jd_c1
        # Time from max to C4
        dt_after = jd_fourth - jd_max

        # These should be within 30% of each other for most eclipses
        ratio = dt_before / dt_after if dt_after > 0 else 0
        assert 0.6 < ratio < 1.7, (
            f"C1-to-max / max-to-C4 ratio ({ratio:.2f}) seems asymmetric. "
            f"Expected between 0.6 and 1.7 for typical eclipses."
        )
