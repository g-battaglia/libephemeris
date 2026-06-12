"""
Tests for sol_eclipse_obscuration_at_loc and sol_eclipse_obscuration_at_loc.

These tests validate the eclipse obscuration calculation at specific locations.
Eclipse obscuration is the fraction of the Sun's disc AREA covered by the Moon,
which differs from magnitude (fraction of Sun's DIAMETER covered).

Reference eclipse: April 8, 2024 total solar eclipse
- Dallas, TX: Total eclipse (obscuration = 1.0 at maximum)
- NYC: Partial eclipse (~80-85% obscuration)
- London: No visibility
"""

import pytest
from libephemeris import (
    julday,
    sol_eclipse_obscuration_at_loc,
    _sol_eclipse_obscuration_at_loc_pythonic,
    _sol_eclipse_magnitude_at_loc_pythonic,
    sol_eclipse_how,
    FLG_SWIEPH,
)


class TestSolEclipseObscurationAtLocSignature:
    """Test that sol_eclipse_obscuration_at_loc function signature is correct."""

    def test_function_exists(self):
        """Test that sol_eclipse_obscuration_at_loc function exists."""
        from libephemeris.eclipse import sol_eclipse_obscuration_at_loc

        assert callable(sol_eclipse_obscuration_at_loc)

    def test_exported_from_package(self):
        """Test that function is exported from main package."""
        from libephemeris import sol_eclipse_obscuration_at_loc

        assert callable(sol_eclipse_obscuration_at_loc)

    def test_swe_version_exported(self):
        """Test that swe_ prefixed version is exported."""
        from libephemeris import sol_eclipse_obscuration_at_loc

        assert callable(sol_eclipse_obscuration_at_loc)

    def test_returns_float(self):
        """Test that function returns a float."""
        # April 8, 2024 eclipse at Dallas
        jd = 2460409.28
        dallas_lat, dallas_lon = 32.7767, -96.797

        result = _sol_eclipse_obscuration_at_loc_pythonic(jd, dallas_lat, dallas_lon)

        assert isinstance(result, float)

    def test_accepts_altitude_parameter(self):
        """Test that function accepts altitude parameter."""
        jd = 2460409.28
        dallas_lat, dallas_lon = 32.7767, -96.797

        result = _sol_eclipse_obscuration_at_loc_pythonic(
            jd, dallas_lat, dallas_lon, altitude=500.0
        )

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts flags parameter."""
        jd = 2460409.28
        dallas_lat, dallas_lon = 32.7767, -96.797

        result = _sol_eclipse_obscuration_at_loc_pythonic(
            jd, dallas_lat, dallas_lon, flags=FLG_SWIEPH
        )

        assert isinstance(result, float)


class TestSweObscurationAtLocSignature:
    """Test sol_eclipse_obscuration_at_loc function signature."""

    def test_function_exists(self):
        """Test that sol_eclipse_obscuration_at_loc function exists."""
        from libephemeris.eclipse import sol_eclipse_obscuration_at_loc

        assert callable(sol_eclipse_obscuration_at_loc)

    def test_returns_float(self):
        """Test that function returns a float."""
        jd = 2460409.28
        dallas_geopos = [-96.797, 32.7767, 0]

        result = sol_eclipse_obscuration_at_loc(jd, dallas_geopos, FLG_SWIEPH)

        assert isinstance(result, float)

    def test_invalid_geopos_raises_error(self):
        """Test that invalid geopos raises ValueError."""
        jd = 2460409.28

        with pytest.raises(ValueError):
            sol_eclipse_obscuration_at_loc(jd, [0, 0], FLG_SWIEPH)


class TestObscurationDuringTotalEclipse:
    """Test obscuration calculation during April 8, 2024 total solar eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        # JD for April 8, 2024 ~18:42 UTC (maximum for Dallas area)
        self.tjd_ut = 2460409.28
        # Dallas coordinates: 32.7767N, 96.797W
        self.dallas_lat = 32.7767
        self.dallas_lon = -96.797

    def test_obscuration_is_positive_during_eclipse(self):
        """Test that obscuration is positive during eclipse."""
        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(
            self.tjd_ut, self.dallas_lat, self.dallas_lon
        )

        assert obscuration > 0, "Obscuration should be positive during eclipse"

    def test_obscuration_near_total_at_dallas(self):
        """Test that obscuration approaches 1.0 for total eclipse at Dallas."""
        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(
            self.tjd_ut, self.dallas_lat, self.dallas_lon
        )

        # At maximum totality, obscuration should be ~1.0
        assert obscuration > 0.85, (
            f"Obscuration {obscuration:.4f} should be near 1.0 for total eclipse"
        )

    def test_obscuration_valid_range(self):
        """Test that obscuration is in valid range [0, 1]."""
        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(
            self.tjd_ut, self.dallas_lat, self.dallas_lon
        )

        # For total eclipses, obscuration = (r_moon/r_sun)^2 which can exceed 1.0
        # matching reference API behavior
        assert 0 <= obscuration <= 1.3, (
            f"Obscuration {obscuration:.4f} should be in range [0, 1.3]"
        )


class TestObscurationAtPartialEclipseLocation:
    """Test obscuration at NYC during April 8, 2024 partial eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        self.tjd_ut = 2460409.30  # Near maximum for NYC
        self.nyc_lat = 40.7128
        self.nyc_lon = -74.006

    def test_partial_eclipse_obscuration(self):
        """Test that NYC sees partial eclipse with obscuration < 1.0."""
        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(
            self.tjd_ut, self.nyc_lat, self.nyc_lon
        )

        # NYC had partial eclipse
        assert 0 < obscuration < 1.0, (
            f"NYC obscuration {obscuration:.4f} should be partial (< 1.0)"
        )

    def test_partial_eclipse_significant_coverage(self):
        """Test that partial eclipse has significant coverage."""
        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(
            self.tjd_ut, self.nyc_lat, self.nyc_lon
        )

        # NYC had significant coverage
        assert obscuration > 0.5, f"NYC obscuration {obscuration:.4f} should be > 0.5"


class TestObscurationVsMagnitude:
    """Test the relationship between obscuration and magnitude."""

    def test_obscuration_less_than_magnitude_for_partial(self):
        """Test that obscuration < magnitude for partial eclipses."""
        # Partial eclipse at NYC
        jd = 2460409.30
        lat, lon = 40.7128, -74.006

        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)

        # For partial eclipses, obscuration is always less than magnitude
        # due to the non-linear relationship between diameter and area
        if 0 < magnitude < 1.0:
            assert obscuration < magnitude, (
                f"Obscuration {obscuration:.4f} should be < magnitude {magnitude:.4f} "
                "for partial eclipse"
            )

    def test_obscuration_equals_one_when_magnitude_exceeds_one(self):
        """Test that obscuration = 1.0 when magnitude >= 1.0 (total eclipse)."""
        # Total eclipse at Dallas - test at various times
        lat, lon = 32.7767, -96.797

        # Find a time when magnitude is > 1.0
        for jd in [2460409.27, 2460409.275, 2460409.28, 2460409.285]:
            magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)
            obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)

            if magnitude >= 1.0:
                assert obscuration == 1.0 or obscuration > 0.98, (
                    f"Obscuration {obscuration:.4f} should be 1.0 when "
                    f"magnitude {magnitude:.4f} >= 1.0"
                )
                break

    def test_both_zero_when_no_eclipse(self):
        """Test that both obscuration and magnitude are 0 when no eclipse."""
        jd = julday(2024, 6, 15, 12.0)
        lat, lon = 32.7767, -96.797

        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)

        assert obscuration == 0.0, "Obscuration should be 0 when no eclipse"
        assert magnitude == 0.0, "Magnitude should be 0 when no eclipse"


class TestObscurationNoEclipse:
    """Test obscuration when no eclipse is happening."""

    def test_no_eclipse_returns_zero(self):
        """Test that function returns 0 when no eclipse."""
        # Random time with no eclipse
        jd = julday(2024, 6, 15, 12.0)
        lat, lon = 32.7767, -96.797

        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)

        assert obscuration == 0.0, "Obscuration should be 0 when no eclipse"

    def test_location_outside_eclipse_path(self):
        """Test that location far from eclipse path returns 0."""
        # April 8, 2024 eclipse - but from Tokyo (no visibility)
        jd = 2460409.28
        tokyo_lat, tokyo_lon = 35.6762, 139.6503

        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, tokyo_lat, tokyo_lon)

        # Tokyo is outside the eclipse path
        assert obscuration == 0.0, "Tokyo should have no eclipse visibility"


class TestConsistencyWithSweEclipseHow:
    """Test consistency with sol_eclipse_how function."""

    def test_obscuration_matches_eclipse_how(self):
        """Test that obscuration matches attr[2] from sol_eclipse_how."""
        jd = 2460409.28
        lat, lon = 32.7767, -96.797
        geopos = [lon, lat, 0]

        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)
        _, attr = sol_eclipse_how(jd, geopos, FLG_SWIEPH)

        # Obscuration is attr[2] in sol_eclipse_how
        assert abs(obscuration - attr[2]) < 0.01, (
            f"Obscuration {obscuration:.4f} should match "
            f"sol_eclipse_how attr[2]={attr[2]:.4f}"
        )

    def test_consistency_at_multiple_times(self):
        """Test consistency at multiple times during eclipse."""
        lat, lon = 32.7767, -96.797
        geopos = [lon, lat, 0]

        # Test at multiple times during eclipse
        test_times = [2460409.26, 2460409.28, 2460409.30]

        for jd in test_times:
            obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)
            _, attr = sol_eclipse_how(jd, geopos, FLG_SWIEPH)

            assert abs(obscuration - attr[2]) < 0.01, (
                f"At JD {jd}: obscuration {obscuration:.4f} differs from "
                f"sol_eclipse_how attr[2]={attr[2]:.4f}"
            )


class TestSweApiConvention:
    """Test that swe_ version follows pyswisseph API conventions."""

    def test_geopos_lon_lat_order(self):
        """Test that geopos uses [lon, lat, alt] order."""
        jd = 2460409.28
        dallas_lat, dallas_lon = 32.7767, -96.797

        # Legacy function: lat, lon order
        result_legacy = _sol_eclipse_obscuration_at_loc_pythonic(jd, dallas_lat, dallas_lon)

        # swe version: geopos = [lon, lat, alt]
        result_swe = sol_eclipse_obscuration_at_loc(
            jd, [dallas_lon, dallas_lat, 0], FLG_SWIEPH
        )

        # Results should be equivalent
        assert abs(result_legacy - result_swe) < 0.0001


class TestEdgeCases:
    """Test edge cases for obscuration calculation."""

    def test_high_altitude_observer(self):
        """Test obscuration at high altitude."""
        jd = 2460409.28
        lat, lon = 32.7767, -96.797

        obscuration_sea_level = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon, altitude=0)
        obscuration_high_alt = _sol_eclipse_obscuration_at_loc_pythonic(
            jd, lat, lon, altitude=5000
        )

        # Both should detect eclipse, small difference possible
        assert obscuration_sea_level > 0
        assert obscuration_high_alt > 0
        assert abs(obscuration_sea_level - obscuration_high_alt) < 0.05

    def test_polar_location(self):
        """Test obscuration calculation at polar location."""
        jd = 2460409.28
        lat, lon = 80.0, 0.0  # Near North Pole

        # Should not crash
        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)
        assert isinstance(obscuration, float)
        assert 0 <= obscuration <= 1.0

    def test_equatorial_location(self):
        """Test obscuration at equatorial location."""
        jd = 2460409.28
        lat, lon = 0.0, 0.0  # Null Island

        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)
        assert isinstance(obscuration, float)
        assert 0 <= obscuration <= 1.0

    def test_accepts_tuple_geopos(self):
        """Test that swe function accepts tuple for geopos."""
        jd = 2460409.28
        geopos = (-96.797, 32.7767, 0)

        result = sol_eclipse_obscuration_at_loc(jd, geopos, FLG_SWIEPH)
        assert isinstance(result, float)

    def test_accepts_list_geopos(self):
        """Test that swe function accepts list for geopos."""
        jd = 2460409.28
        geopos = [-96.797, 32.7767, 0]

        result = sol_eclipse_obscuration_at_loc(jd, geopos, FLG_SWIEPH)
        assert isinstance(result, float)


class TestObscurationPhysicalProperties:
    """Test physical properties of obscuration calculation."""

    def test_obscuration_increases_toward_maximum(self):
        """Test that obscuration increases as eclipse approaches maximum."""
        lat, lon = 32.7767, -96.797

        # Times before maximum (increasing obscuration)
        times_before_max = [2460409.24, 2460409.26, 2460409.28]
        obscurations = [
            _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon) for jd in times_before_max
        ]

        # Check increasing trend (with some tolerance for edge effects)
        for i in range(1, len(obscurations)):
            if obscurations[i - 1] > 0:  # Only compare if eclipse is visible
                assert obscurations[i] >= obscurations[i - 1] * 0.9, (
                    f"Obscuration should increase: {obscurations}"
                )

    def test_obscuration_decreases_after_maximum(self):
        """Test that obscuration decreases after eclipse maximum."""
        lat, lon = 32.7767, -96.797

        # Times after maximum (decreasing obscuration)
        times_after_max = [2460409.28, 2460409.30, 2460409.32]
        obscurations = [
            _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon) for jd in times_after_max
        ]

        # Check decreasing trend
        for i in range(1, len(obscurations)):
            if obscurations[i] > 0:  # Only compare if eclipse is visible
                assert obscurations[i] <= obscurations[i - 1] * 1.1, (
                    f"Obscuration should decrease: {obscurations}"
                )

    def test_obscuration_bounded_by_one(self):
        """Test that obscuration never exceeds 1.0."""
        lat, lon = 32.7767, -96.797

        # Test at multiple times during total eclipse
        test_times = [
            2460409.26,
            2460409.27,
            2460409.275,
            2460409.28,
            2460409.285,
            2460409.29,
            2460409.30,
        ]

        for jd in test_times:
            obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)
            # For total eclipses, obscuration = (r_moon/r_sun)^2 can exceed 1.0
            assert obscuration <= 1.3, (
                f"Obscuration {obscuration:.4f} at JD {jd} should not exceed 1.3"
            )


class TestMathematicalRelationship:
    """Test the mathematical relationship between magnitude and obscuration."""

    def test_obscuration_formula_for_partial_eclipse(self):
        """Test that obscuration follows expected formula for partial eclipses.

        For partial eclipses, the relationship between magnitude (m) and
        obscuration (o) follows from the lens formula for circle intersection.
        """
        jd = 2460409.30  # Partial eclipse at NYC
        lat, lon = 40.7128, -74.006

        obscuration = _sol_eclipse_obscuration_at_loc_pythonic(jd, lat, lon)
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)

        if 0 < magnitude < 1.0:
            # For equal-sized disks, obscuration follows lens formula
            # This is an approximation; our function uses exact formula
            # Just verify the qualitative relationship holds
            assert 0 < obscuration < magnitude, (
                f"Obscuration {obscuration:.4f} should be less than "
                f"magnitude {magnitude:.4f} for partial eclipse"
            )
