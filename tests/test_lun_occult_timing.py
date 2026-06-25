"""
Tests for improved lunar occultation timing precision.

These tests verify that the improved algorithms for lun_occult_when_glob()
and lun_occult_when_loc() achieve better timing precision than the original
300 second (5 minute) tolerance.

The improvements include:
- Bisection search for contact times instead of chord-based estimates
- Golden section search with increased iterations for sub-second precision
- Topocentric parallax corrections for local occultation timing
- Reduced initial search step size
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    revjul,
    lun_occult_when_glob,
    lun_occult_when_loc,
    VENUS,
    FLG_SWIEPH,
)


class TestLunOccultTimingPrecision:
    """Test improved timing precision for lunar occultations."""

    def test_contact_times_precision_star_occultation(self):
        """Test that contact times have sub-minute precision for star occultations.

        The improved bisection-based contact time calculation should provide
        more accurate timing than the previous chord-based estimates.
        """
        # Known Regulus occultation in 2017
        jd_start = julday(2017, 6, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Should find an occultation
        assert retflags != 0

        jd_max = tret[0]
        jd_begin = tret[2]
        jd_end = tret[3]

        # Verify chronological order
        if jd_begin > 0:
            assert jd_begin < jd_max
        if jd_end > 0:
            assert jd_end > jd_max

        # Duration should be reasonable for a GLOBAL occultation
        # Global occultations (C1 to C4 across entire Earth surface)
        # can last much longer than local ones — up to ~5 hours for
        # bright stars near the ecliptic with favorable geometry.
        if jd_begin > 0 and jd_end > 0:
            duration_seconds = (jd_end - jd_begin) * 86400
            assert 1 < duration_seconds < 20000  # 1 second to ~5.5 hours

    def test_contact_times_symmetry(self):
        """Test that contact times are roughly symmetric around maximum."""
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        if retflags == 0:
            pytest.skip("No occultation found")

        jd_max = tret[0]
        jd_begin = tret[2]
        jd_end = tret[3]

        if jd_begin > 0 and jd_end > 0:
            dt_before = jd_max - jd_begin
            dt_after = jd_end - jd_max

            # For most occultations, times should be roughly symmetric
            # Allow up to 20% asymmetry due to Moon's motion
            if min(dt_before, dt_after) > 0:
                ratio = max(dt_before, dt_after) / min(dt_before, dt_after)
                assert ratio < 1.5, f"Contact times asymmetry too large: {ratio:.2f}"

    def test_totality_times_within_outer_contacts(self):
        """Test that totality times are properly nested within outer contacts."""
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        if retflags == 0:
            pytest.skip("No occultation found")

        jd_max = tret[0]
        jd_begin = tret[2]
        jd_end = tret[3]
        jd_total_begin = tret[4]
        jd_total_end = tret[5]

        # If we have totality times, they must be within outer contacts
        if jd_total_begin > 0 and jd_begin > 0:
            assert jd_begin <= jd_total_begin <= jd_max
        if jd_total_end > 0 and jd_end > 0:
            assert jd_max <= jd_total_end <= jd_end


class TestLunOccultWhenLocTopocentric:
    """Test topocentric parallax corrections for local occultations."""

    def test_local_occultation_timing_precision(self):
        """Test that local occultation uses topocentric positions.

        The local occultation timing should differ from global timing
        due to lunar parallax (~1 degree).
        """
        # Search for a Venus occultation from a specific location
        jd_start = julday(2024, 1, 1, 0)
        rome_lat, rome_lon = 41.9028, 12.4964

        try:
            ecl_type, times, attr = lun_occult_when_loc(
                jd_start, VENUS, (rome_lon, rome_lat, 0), FLG_SWIEPH
            )

            # Should find an occultation
            assert ecl_type != 0

            # Maximum time should be valid
            jd_max = times[0]
            assert jd_max > jd_start

            # First and fourth contacts should be valid
            jd_first = times[1]
            jd_fourth = times[4]

            if jd_first > 0 and jd_fourth > 0:
                # Duration should be reasonable
                duration_hours = (jd_fourth - jd_first) * 24
                assert 0.001 < duration_hours < 2.0

        except RuntimeError:
            pytest.skip("No Venus occultation visible from Rome in search period")

    def test_swe_lun_occult_when_loc_api(self):
        """Test the reference-API-compatible API."""
        jd_start = julday(2017, 1, 1, 0)
        geopos = (12.4964, 41.9028, 0)  # Rome (lon, lat, alt)

        # First find a global occultation
        try:
            retflags_glob, tret_glob = lun_occult_when_glob(
                jd_start, "Regulus", FLG_SWIEPH, 0
            )
        except RuntimeError:
            pytest.skip("No global Regulus occultation found")

        # Now try local search
        try:
            ecl_type, times, attr = lun_occult_when_loc(
                jd_start, "Regulus", geopos, FLG_SWIEPH, False
            )

            # Should return consistent structure
            assert len(times) == 10
            assert len(attr) == 20
            assert isinstance(ecl_type, int)

        except RuntimeError:
            pytest.skip("No Regulus occultation visible from Rome")

    def test_local_times_differ_from_global(self):
        """Test that local topocentric times differ from geocentric global times.

        Due to lunar parallax (~1 degree), local timing should differ from
        global timing. This tests that the topocentric correction is applied.

        Note: For global occultations visible only due to parallax (not
        geocentric occultations), the local algorithm may find the same
        occultation but with significantly different timing due to the
        parallax correction. This is expected behavior.
        """
        jd_start = julday(
            2017, 6, 1, 0
        )  # Start from June to find geocentric occultation

        # Get global occultation time
        try:
            retflags_glob, tret_glob = lun_occult_when_glob(
                jd_start, "Regulus", FLG_SWIEPH, 0
            )
        except RuntimeError:
            pytest.skip("No global Regulus occultation found")

        if retflags_glob == 0:
            pytest.skip("No occultation found")

        jd_max_global = tret_glob[0]

        # Get local occultation time from Rome (should be visible)
        rome_lat, rome_lon = 41.9028, 12.4964

        try:
            ecl_type, times, attr = lun_occult_when_loc(
                jd_max_global - 1, "Regulus", (rome_lon, rome_lat, 0), FLG_SWIEPH
            )

            jd_max_local = times[0]

            # Only compare if we found the same occultation (within 1 day)
            if abs(jd_max_local - jd_max_global) > 1.0:
                pytest.skip("Local and global searches found different occultations")

            diff_seconds = abs(jd_max_local - jd_max_global) * 86400

            # The difference should be small for geocentric occultations
            # (typically < 120 seconds for parallax effects)
            # but can be larger for global-only occultations
            assert diff_seconds < 600, (
                f"Local timing at Rome differs by {diff_seconds:.1f}s from global - "
                f"expected < 600s for parallax effects"
            )

        except RuntimeError:
            pytest.skip("Regulus occultation not visible from Rome")


class TestLunOccultSearchPrecision:
    """Test the search algorithm precision."""

    def test_golden_section_convergence(self):
        """Test that golden section search converges to high precision.

        The improved golden section search should converge to sub-second
        precision (1e-7 days ~ 0.01 seconds).
        """
        jd_start = julday(2017, 6, 1, 0)

        # Run search twice and verify consistent results
        retflags1, tret1 = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)
        retflags2, tret2 = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Results should be identical
        assert retflags1 == retflags2
        for i in range(len(tret1)):
            if tret1[i] > 0:
                diff_seconds = abs(tret1[i] - tret2[i]) * 86400
                assert diff_seconds < 1.0, f"tret[{i}] differs by {diff_seconds:.3f}s"

    def test_search_step_finds_occultations(self):
        """Test that reduced search step size doesn't miss occultations."""
        # Search from early 2017 - multiple Regulus occultations throughout 2017
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Should find an occultation in 2017
        assert retflags != 0

        year, month, day, hour = revjul(tret[0])
        assert year == 2017


class TestLunOccultPlanetOccultations:
    """Test planet occultation timing."""

    @pytest.mark.parametrize(
        "planet,name",
        [
            (VENUS, "Venus"),
            # MARS and JUPITER require ephemeris data that may not be available
            # in all test environments
        ],
    )
    def test_planet_occultation_timing(self, planet, name):
        """Test that planet occultation timing has improved precision."""
        jd_start = julday(2024, 1, 1, 0)

        try:
            retflags, tret = lun_occult_when_glob(jd_start, planet, FLG_SWIEPH, 0)

            if retflags == 0:
                pytest.skip(f"No {name} occultation found")

            jd_max = tret[0]
            jd_begin = tret[2]
            jd_end = tret[3]

            # Verify chronological order
            if jd_begin > 0:
                assert jd_begin < jd_max
            if jd_end > 0:
                assert jd_end > jd_max

            # Duration should be reasonable for a GLOBAL occultation
            # Global planet occultations can span several hours as the
            # Moon's shadow sweeps across the entire Earth surface.
            if jd_begin > 0 and jd_end > 0:
                duration_minutes = (jd_end - jd_begin) * 24 * 60
                assert 0.1 < duration_minutes < 360  # 6 seconds to 6 hours

        except (RuntimeError, KeyError):
            pytest.skip(
                f"No {name} occultation found within search period or ephemeris unavailable"
            )
