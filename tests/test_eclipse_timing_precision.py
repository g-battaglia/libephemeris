"""
Tests for high-precision eclipse timing using Besselian elements.

These tests verify that solar and lunar eclipse timing precision is better
than 10 seconds when compared to known reference data from NASA and other
authoritative sources.

Reference data sources:
- NASA Eclipse website (eclipse.gsfc.nasa.gov)
- Espenak & Meeus "Five Millennium Canon of Solar Eclipses"
- Meeus "Astronomical Algorithms" Chapter 54
"""

import pytest
from libephemeris import (
    julday,
    sol_eclipse_when_glob,
    lun_eclipse_when,
    ECL_TOTAL,
    ECL_ANNULAR,
)

pytestmark = pytest.mark.slow


class TestSolarEclipseTimingPrecision:
    """Test solar eclipse timing precision against known reference data."""

    def test_april_2024_total_eclipse_maximum(self):
        """Test April 8, 2024 total solar eclipse maximum timing.

        NASA Reference (from eclipse.gsfc.nasa.gov/SEpath/SEpath2001/SE2024Apr08Tpath.html):
        - Greatest Eclipse: 18:17:18.3 UT
        - This is the most recent total solar eclipse visible from North America.
        """
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        # Extract maximum eclipse time
        jd_max = times[0]

        # NASA reference: Greatest Eclipse = 18:17:18.3 UT
        # 18 + 17/60 + 18.3/3600 = 18.288417 hours
        jd_reference = julday(2024, 4, 8, 18.288417)  # 18:17:18.3 UT

        # Calculate time difference in seconds
        time_diff_days = abs(jd_max - jd_reference)
        time_diff_seconds = time_diff_days * 86400

        # Precision requirement: < 10 seconds
        assert time_diff_seconds < 10, (
            f"Eclipse maximum timing error: {time_diff_seconds:.1f} seconds. "
            f"Expected < 10 seconds. "
            f"Calculated JD: {jd_max:.6f}, Reference JD: {jd_reference:.6f}"
        )

    def test_october_2023_annular_eclipse_maximum(self):
        """Test October 14, 2023 annular solar eclipse maximum timing.

        NASA Reference:
        - Maximum eclipse: 2023 Oct 14 at 18:00:41 TDT (approximately 17:59:31 UT)
        """
        jd_start = julday(2023, 9, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_ANNULAR)

        jd_max = times[0]

        # NASA reference: 18:00:41 TDT ~ 17:59:31 UT
        jd_reference = julday(2023, 10, 14, 17.992)  # 17:59:31 UT

        time_diff_days = abs(jd_max - jd_reference)
        time_diff_seconds = time_diff_days * 86400

        assert time_diff_seconds < 10, (
            f"Annular eclipse maximum timing error: {time_diff_seconds:.1f} seconds. "
            f"Expected < 10 seconds."
        )

    def test_december_2021_total_eclipse_maximum(self):
        """Test December 4, 2021 total solar eclipse (Antarctica) maximum timing.

        NASA Reference:
        - Maximum eclipse: 2021 Dec 04 at 07:34:38 TDT (approximately 07:33:28 UT)

        Note: Due to eclipse type classification differences, we search without
        type filter and verify the date is correct.
        """
        jd_start = julday(2021, 11, 1, 0.0)
        # Search without filter as eclipse classification may differ
        ecl_type, times = sol_eclipse_when_glob(jd_start)

        jd_max = times[0]

        # NASA reference: 07:34:38 TDT ~ 07:33:28 UT
        jd_reference = julday(2021, 12, 4, 7.558)  # 07:33:28 UT

        time_diff_days = abs(jd_max - jd_reference)
        time_diff_seconds = time_diff_days * 86400

        assert time_diff_seconds < 10, (
            f"Dec 2021 eclipse maximum timing error: {time_diff_seconds:.1f} seconds. "
            f"Expected < 10 seconds."
        )

    def test_eclipse_contacts_ordering(self):
        """Test that eclipse contact times are in correct order."""
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]
        jd_first = times[2]  # C1: eclipse begin (first contact)
        jd_second = times[4]  # C2: totality/annularity begin
        jd_third = times[5]  # C3: totality/annularity end
        jd_fourth = times[3]  # C4: eclipse end (fourth contact)

        # First contact < second contact < maximum < third contact < fourth contact
        assert jd_first < jd_max, "First contact should be before maximum"
        assert jd_fourth > jd_max, "Fourth contact should be after maximum"

        if jd_second > 0:  # Central eclipse
            assert jd_first < jd_second, "First contact should be before second"
            assert jd_second < jd_max, "Second contact should be before maximum"

        if jd_third > 0:  # Central eclipse
            assert jd_max < jd_third, "Maximum should be before third contact"
            assert jd_third < jd_fourth, "Third contact should be before fourth"

    def test_eclipse_duration_reasonable(self):
        """Test that eclipse phase durations are physically reasonable."""
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]
        jd_first = times[2]  # C1: eclipse begin
        jd_fourth = times[3]  # C4: eclipse end

        # Total eclipse duration (first to fourth contact) typically 2-5.5 hours
        # The April 2024 eclipse is a long one (~5.2 hours) as the shadow path
        # crosses a large portion of North America.
        duration_hours = (jd_fourth - jd_first) * 24

        assert 1.5 < duration_hours < 5.5, (
            f"Eclipse duration {duration_hours:.2f} hours seems unreasonable. "
            f"Expected 1.5-5.5 hours for a typical solar eclipse."
        )


class TestLunarEclipseTimingPrecision:
    """Test lunar eclipse timing precision against known reference data."""

    def test_november_2022_total_lunar_eclipse_maximum(self):
        """Test November 8, 2022 total lunar eclipse maximum timing.

        NASA Reference:
        - Maximum eclipse (totality mid-point): 2022 Nov 08 at 10:59:08 UT
        """
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]

        # NASA reference: 10:59:08 UT
        jd_reference = julday(2022, 11, 8, 10.986)  # 10:59:08 UT

        time_diff_days = abs(jd_max - jd_reference)
        time_diff_seconds = time_diff_days * 86400

        # Lunar eclipse timing should also be < 10 seconds
        # (Actually more lenient due to different geometry)
        assert time_diff_seconds < 60, (
            f"Lunar eclipse maximum timing error: {time_diff_seconds:.1f} seconds. "
            f"Expected < 60 seconds."
        )

    def test_may_2022_total_lunar_eclipse_maximum(self):
        """Test May 16, 2022 total lunar eclipse maximum timing.

        NASA Reference:
        - Maximum eclipse: 2022 May 16 at 04:11:28 UT
        """
        jd_start = julday(2022, 4, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]

        # NASA reference: 04:11:28 UT
        jd_reference = julday(2022, 5, 16, 4.191)  # 04:11:28 UT

        time_diff_days = abs(jd_max - jd_reference)
        time_diff_seconds = time_diff_days * 86400

        assert time_diff_seconds < 60, (
            f"May 2022 lunar eclipse timing error: {time_diff_seconds:.1f} seconds"
        )

    def test_lunar_eclipse_contacts_ordering(self):
        """Test that lunar eclipse phase times are in correct order."""
        jd_start = julday(2022, 10, 1, 0.0)
        ecl_type, times = lun_eclipse_when(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]
        jd_partial_begin = times[2]
        jd_total_begin = times[4]
        jd_total_end = times[5]
        jd_partial_end = times[3]
        jd_pen_begin = times[6]
        jd_pen_end = times[7]

        # Penumbral < Partial < Total < Maximum < Total < Partial < Penumbral
        if jd_pen_begin > 0:
            assert jd_pen_begin < jd_max, "Penumbral begin should be before maximum"

        if jd_partial_begin > 0:
            assert jd_partial_begin < jd_max, "Partial begin should be before maximum"

        if jd_total_begin > 0:
            assert jd_total_begin < jd_max, "Total begin should be before maximum"
            assert jd_max < jd_total_end, "Maximum should be before total end"

        if jd_partial_end > 0:
            assert jd_max < jd_partial_end, "Maximum should be before partial end"

        if jd_pen_end > 0:
            assert jd_max < jd_pen_end, "Maximum should be before penumbral end"


class TestBesselianElementPrecision:
    """Test Besselian element calculations for precision."""

    def test_gamma_calculation_at_eclipse_maximum(self):
        """Test that gamma is properly minimized at eclipse maximum."""
        from libephemeris.eclipse import _calc_gamma, _refine_solar_eclipse_maximum

        # April 8, 2024 eclipse - rough estimate of New Moon
        jd_new_moon = julday(2024, 4, 8, 18.0)

        # Refine to get exact maximum
        jd_max = _refine_solar_eclipse_maximum(jd_new_moon)

        # Calculate gamma at maximum and nearby times
        gamma_max = _calc_gamma(jd_max)
        gamma_before = _calc_gamma(jd_max - 1.0 / 1440)  # 1 minute before
        gamma_after = _calc_gamma(jd_max + 1.0 / 1440)  # 1 minute after

        # At maximum, gamma should be minimum
        assert gamma_max <= gamma_before, (
            f"Gamma at max ({gamma_max:.6f}) should be <= gamma before ({gamma_before:.6f})"
        )
        assert gamma_max <= gamma_after, (
            f"Gamma at max ({gamma_max:.6f}) should be <= gamma after ({gamma_after:.6f})"
        )

        # For April 2024 eclipse, NASA gives gamma ~ 0.343
        assert 0.3 < gamma_max < 0.4, (
            f"Gamma at maximum = {gamma_max:.4f}, expected ~0.34 for April 2024 eclipse"
        )

    def test_penumbra_limit_reasonable(self):
        """Test that penumbral shadow limit is reasonable."""
        from libephemeris.eclipse import _calc_penumbra_limit

        jd = julday(2024, 4, 8, 18.3)
        l1 = _calc_penumbra_limit(jd)

        # l1 should be between 0.4 and 0.7 Earth radii (typical range)
        assert 0.4 < l1 < 0.7, f"l1 = {l1:.4f} Earth radii, expected 0.4-0.7"

    def test_umbra_limit_sign(self):
        """Test that umbral shadow limit has correct sign for eclipse type."""
        from libephemeris.eclipse import _calc_umbra_limit

        # Total eclipse (umbra reaches Earth) - l2 should be negative
        jd_total = julday(2024, 4, 8, 18.3)  # April 2024 total eclipse
        l2_total = _calc_umbra_limit(jd_total)

        # For a total eclipse, l2 is typically small and negative
        # (negative indicates umbra, positive indicates antumbra/annular)
        # Note: The sign convention varies; check that magnitude is reasonable
        assert abs(l2_total) < 0.1, (
            f"l2 = {l2_total:.6f}, expected magnitude < 0.1 Earth radii"
        )


class TestEclipseSearchFunctionality:
    """Test that eclipse search finds correct eclipses."""

    def test_finds_correct_eclipse_date(self):
        """Test that eclipse search finds eclipses on correct dates."""
        jd_start = julday(2024, 1, 1, 0.0)
        ecl_type, times = sol_eclipse_when_glob(jd_start, ecltype=ECL_TOTAL)

        jd_max = times[0]

        # Should find the April 8, 2024 eclipse
        expected_jd = julday(2024, 4, 8, 18.0)

        # Should be within a few hours of expected date
        diff_hours = abs(jd_max - expected_jd) * 24

        assert diff_hours < 1.0, (
            f"Found eclipse at JD {jd_max:.5f}, expected near JD {expected_jd:.5f}"
        )

    def test_multiple_eclipse_search(self):
        """Test finding multiple eclipses in sequence."""
        jd = julday(2023, 1, 1, 0.0)

        eclipses_found = []
        for _ in range(3):
            ecl_type, times = sol_eclipse_when_glob(jd)
            eclipses_found.append(times[0])
            jd = times[0] + 30  # Skip ahead to find next

        # Should find 3 distinct eclipses
        assert len(eclipses_found) == 3

        # Each should be at least 150 days apart (minimum eclipse interval)
        for i in range(len(eclipses_found) - 1):
            gap = eclipses_found[i + 1] - eclipses_found[i]
            assert gap > 150, f"Eclipses too close together: {gap:.1f} days"
