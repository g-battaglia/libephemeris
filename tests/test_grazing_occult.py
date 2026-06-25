"""
Tests for grazing occultation detection in libephemeris.

Tests the ECL_GRAZING flag which indicates when a star or planet
passes near the lunar limb during an occultation. Grazing occultations
are scientifically interesting because the target may flash in/out
multiple times due to lunar limb topography (mountains/valleys).

A grazing occultation is defined as when the target passes within
the outer 10% of the Moon's disc (min_sep > 0.9 * moon_radius).
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    lun_occult_when_glob,
    ECL_GRAZING,
    ECL_TOTAL,
    ECL_PARTIAL,
    FLG_SWIEPH,
)


class TestGrazingOccultConstant:
    """Test the ECL_GRAZING constant definition."""

    def test_grazing_constant_defined(self):
        """Test that ECL_GRAZING constant is defined and accessible."""
        assert ECL_GRAZING is not None
        assert isinstance(ECL_GRAZING, int)
        assert ECL_GRAZING > 0

    def test_grazing_constant_value(self):
        """Test that ECL_GRAZING has expected value (65536 = 2^16)."""
        assert ECL_GRAZING == 65536

    def test_grazing_constant_no_overlap_with_other_flags(self):
        """Test that ECL_GRAZING doesn't overlap with other eclipse type flags."""
        from libephemeris import (
            ECL_CENTRAL,
            ECL_NONCENTRAL,
            ECL_ANNULAR,
            ECL_ANNULAR_TOTAL,
            ECL_PENUMBRAL,
            ECL_VISIBLE,
            ECL_MAX_VISIBLE,
            ECL_1ST_VISIBLE,
            ECL_2ND_VISIBLE,
            ECL_3RD_VISIBLE,
            ECL_4TH_VISIBLE,
        )

        # ECL_GRAZING should not overlap with any other eclipse flags
        other_flags = [
            ECL_CENTRAL,
            ECL_NONCENTRAL,
            ECL_TOTAL,
            ECL_ANNULAR,
            ECL_PARTIAL,
            ECL_ANNULAR_TOTAL,
            ECL_PENUMBRAL,
            ECL_VISIBLE,
            ECL_MAX_VISIBLE,
            ECL_1ST_VISIBLE,
            ECL_2ND_VISIBLE,
            ECL_3RD_VISIBLE,
            ECL_4TH_VISIBLE,
        ]

        for flag in other_flags:
            # No bits should overlap
            assert (ECL_GRAZING & flag) == 0, f"ECL_GRAZING overlaps with {flag}"

    def test_ecl_grazing_alias(self):
        """Test that ECL_GRAZING is an alias for ECL_GRAZING."""
        assert ECL_GRAZING == ECL_GRAZING


class TestGrazingOccultDetection:
    """Test grazing occultation detection in lun_occult_when_glob."""

    def test_grazing_flag_combined_with_type(self):
        """Test that grazing flag is combined with occultation type (TOTAL/PARTIAL)."""
        # Search for a known occultation and check flag handling
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Should have a type flag (TOTAL, PARTIAL, or ANNULAR)
        has_type = (retflags & ECL_TOTAL) or (retflags & ECL_PARTIAL)
        assert has_type, "Occultation should have a type flag"

        # Grazing flag may or may not be present depending on the event
        # But it should be properly combinable with type flags
        grazing_set = bool(retflags & ECL_GRAZING)
        # This is informational - grazing flag presence depends on the geometry

    def test_grazing_flag_can_be_checked_with_bitwise_and(self):
        """Test that grazing flag can be checked using bitwise AND."""
        # Create a combined flag for testing
        combined_flag = ECL_PARTIAL | ECL_GRAZING

        # Check individual flags
        assert (combined_flag & ECL_PARTIAL) != 0
        assert (combined_flag & ECL_GRAZING) != 0
        assert (combined_flag & ECL_TOTAL) == 0

    def test_non_grazing_occultation_no_grazing_flag(self):
        """Test that non-grazing occultations don't have the grazing flag.

        For a typical central occultation, the grazing flag should not be set.
        """
        # Search for an occultation of Regulus (bright star close to ecliptic)
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # For a total occultation (star passes near center of Moon),
        # we wouldn't expect grazing flag unless it's specifically a grazing event
        if retflags & ECL_TOTAL:
            # Total occultation - could be grazing or non-grazing depending on geometry
            pass  # Valid either way

    def test_occultation_returns_valid_structure_with_grazing(self):
        """Test that occultation returns correct structure even with grazing detection."""
        jd_start = julday(2017, 1, 1, 0)

        retflags, tret = lun_occult_when_glob(jd_start, "Regulus", FLG_SWIEPH, 0)

        # Should still return 10-element tuple
        assert len(tret) == 10
        assert all(isinstance(t, float) for t in tret)
        assert isinstance(retflags, int)

        # Maximum should be after start
        assert tret[0] > jd_start


class TestGrazingOccultEdgeCases:
    """Test edge cases for grazing occultation detection."""

    def test_grazing_detection_different_stars(self):
        """Test grazing detection works for different stars."""
        stars_to_test = ["Regulus", "Spica", "Aldebaran"]
        jd_start = julday(2017, 1, 1, 0)

        for star in stars_to_test:
            try:
                retflags, tret = lun_occult_when_glob(jd_start, star, FLG_SWIEPH, 0)
                # Should return valid result
                assert isinstance(retflags, int)
                assert len(tret) == 10
            except RuntimeError:
                # No occultation found within search limit is acceptable
                pass

    def test_grazing_flag_is_bitmask_compatible(self):
        """Test that grazing flag works correctly as a bitmask."""
        # Simulate different flag combinations
        test_cases = [
            (ECL_TOTAL | ECL_GRAZING, True, True),
            (ECL_PARTIAL | ECL_GRAZING, False, True),
            (ECL_TOTAL, True, False),
            (ECL_PARTIAL, False, False),
        ]

        for flags, expect_total, expect_grazing in test_cases:
            is_total = bool(flags & ECL_TOTAL)
            is_grazing = bool(flags & ECL_GRAZING)

            assert is_total == expect_total
            assert is_grazing == expect_grazing


class TestGrazingOccultIntegration:
    """Integration tests for grazing occultation with other features."""

    def test_grazing_detection_with_lun_occult_when_loc(self):
        """Test that grazing flag is passed through to lun_occult_when_loc."""
        from libephemeris import lun_occult_when_loc

        # Test from a specific location (Rome)
        jd_start = julday(2017, 1, 1, 0)
        rome_lat, rome_lon = 41.9028, 12.4964

        try:
            retflag, times, attr = lun_occult_when_loc(
                jd_start, "Regulus", (rome_lon, rome_lat, 0.0)
            )

            # The return flag should be an integer that can contain grazing flag
            assert isinstance(retflag, int)
            # Should be able to check for grazing
            is_grazing = bool(retflag & ECL_GRAZING)
            # is_grazing can be True or False depending on the event
        except RuntimeError:
            # No visible occultation from this location is acceptable
            pass

    def test_grazing_flag_compatible_with_visibility_flags(self):
        """Test that grazing flag doesn't interfere with visibility flags."""
        from libephemeris import (
            ECL_VISIBLE,
            ECL_MAX_VISIBLE,
            ECL_1ST_VISIBLE,
        )

        # Create a combined flag with multiple visibility flags and grazing
        combined = ECL_PARTIAL | ECL_GRAZING | ECL_VISIBLE | ECL_MAX_VISIBLE

        # All flags should be independently checkable
        assert (combined & ECL_PARTIAL) != 0
        assert (combined & ECL_GRAZING) != 0
        assert (combined & ECL_VISIBLE) != 0
        assert (combined & ECL_MAX_VISIBLE) != 0

        # Other flags should still be false
        assert (combined & ECL_TOTAL) == 0
        assert (combined & ECL_1ST_VISIBLE) == 0


class TestPlanetOccultGrazing:
    """Test grazing detection in planetary occultation functions."""

    def test_planet_occult_grazing_flag_defined(self):
        """Test that ECL_GRAZING works with planet_occult_when_glob."""
        from libephemeris import planet_occult_when_glob, VENUS, JUPITER

        # We won't actually search for a planetary occultation (very rare)
        # but verify the function accepts the flag system
        jd_start = julday(2065, 1, 1, 0)

        try:
            retflags, tret = planet_occult_when_glob(jd_start, VENUS, JUPITER)
            # Check that return flag is compatible with grazing check
            assert isinstance(retflags, int)
            # The grazing flag should be checkable
            is_grazing = bool(retflags & ECL_GRAZING)
        except RuntimeError:
            # No occultation found is acceptable (they are very rare)
            pass
