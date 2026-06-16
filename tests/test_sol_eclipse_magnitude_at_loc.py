"""
Tests for sol_eclipse_magnitude_at_loc and sol_eclipse_magnitude_at_loc functions.

These tests validate the eclipse magnitude calculation at specific geographic locations.
Eclipse magnitude is defined as the fraction of the Sun's diameter covered by the Moon.

Reference eclipse: April 8, 2024 total solar eclipse
- Dallas, TX: Total eclipse (magnitude > 1.0 at maximum)
- NYC: Partial eclipse (~90% magnitude)
- London: No visibility
"""

import pytest
from libephemeris import (
    julday,
    sol_eclipse_magnitude_at_loc,
    _sol_eclipse_magnitude_at_loc_pythonic,
    sol_eclipse_how,
    FLG_SWIEPH,
)


class TestSolEclipseMagnitudeAtLocSignature:
    """Test that sol_eclipse_magnitude_at_loc function signature is correct."""

    def test_function_exists(self):
        """Test that sol_eclipse_magnitude_at_loc function exists."""
        from libephemeris.eclipse import sol_eclipse_magnitude_at_loc

        assert callable(sol_eclipse_magnitude_at_loc)

    def test_exported_from_package(self):
        """Test that function is exported from main package."""
        from libephemeris import sol_eclipse_magnitude_at_loc

        assert callable(sol_eclipse_magnitude_at_loc)

    def test_swe_version_exported(self):
        """Test that swe_ prefixed version is exported."""
        from libephemeris import sol_eclipse_magnitude_at_loc

        assert callable(sol_eclipse_magnitude_at_loc)

    def test_returns_float(self):
        """Test that function returns a float."""
        # April 8, 2024 eclipse at Dallas
        jd = 2460409.28
        dallas_lat, dallas_lon = 32.7767, -96.797

        result = _sol_eclipse_magnitude_at_loc_pythonic(jd, dallas_lat, dallas_lon)

        assert isinstance(result, float)

    def test_accepts_altitude_parameter(self):
        """Test that function accepts altitude parameter."""
        jd = 2460409.28
        dallas_lat, dallas_lon = 32.7767, -96.797

        result = _sol_eclipse_magnitude_at_loc_pythonic(
            jd, dallas_lat, dallas_lon, altitude=500.0
        )

        assert isinstance(result, float)

    def test_accepts_flags_parameter(self):
        """Test that function accepts flags parameter."""
        jd = 2460409.28
        dallas_lat, dallas_lon = 32.7767, -96.797

        result = _sol_eclipse_magnitude_at_loc_pythonic(
            jd, dallas_lat, dallas_lon, flags=FLG_SWIEPH
        )

        assert isinstance(result, float)


class TestSweSolEclipseMagnitudeAtLocSignature:
    """Test sol_eclipse_magnitude_at_loc function signature."""

    def test_function_exists(self):
        """Test that sol_eclipse_magnitude_at_loc function exists."""
        from libephemeris.eclipse import sol_eclipse_magnitude_at_loc

        assert callable(sol_eclipse_magnitude_at_loc)

    def test_returns_float(self):
        """Test that function returns a float."""
        jd = 2460409.28
        dallas_geopos = [-96.797, 32.7767, 0]

        result = sol_eclipse_magnitude_at_loc(jd, dallas_geopos, FLG_SWIEPH)

        assert isinstance(result, float)

    def test_invalid_geopos_raises_error(self):
        """Test that invalid geopos raises ValueError."""
        jd = 2460409.28

        with pytest.raises(ValueError):
            sol_eclipse_magnitude_at_loc(jd, [0, 0], FLG_SWIEPH)


class TestMagnitudeDuringTotalEclipse:
    """Test magnitude calculation during April 8, 2024 total solar eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        # JD for April 8, 2024 ~18:42 UTC (maximum for Dallas area)
        self.tjd_ut = 2460409.28
        # Dallas coordinates: 32.7767N, 96.797W
        self.dallas_lat = 32.7767
        self.dallas_lon = -96.797

    def test_magnitude_is_positive_during_eclipse(self):
        """Test that magnitude is positive during eclipse."""
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(
            self.tjd_ut, self.dallas_lat, self.dallas_lon
        )

        assert magnitude > 0, "Magnitude should be positive during eclipse"

    def test_magnitude_near_total_at_dallas(self):
        """Test that magnitude approaches or exceeds 1.0 for total eclipse at Dallas."""
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(
            self.tjd_ut, self.dallas_lat, self.dallas_lon
        )

        # At maximum totality, magnitude should be >= 1.0
        # Allow some tolerance since we may not be at exact maximum
        assert magnitude > 0.9, (
            f"Magnitude {magnitude:.4f} should be near 1.0 for total eclipse"
        )

    def test_magnitude_reasonable_range(self):
        """Test that magnitude is in a reasonable range."""
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(
            self.tjd_ut, self.dallas_lat, self.dallas_lon
        )

        # Magnitude should be between 0 and ~1.1 (can exceed 1.0 for total)
        assert 0 <= magnitude <= 1.5, f"Magnitude {magnitude:.4f} out of expected range"


class TestMagnitudeAtPartialEclipseLocation:
    """Test magnitude at NYC during April 8, 2024 partial eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        self.tjd_ut = 2460409.30  # Near maximum for NYC
        self.nyc_lat = 40.7128
        self.nyc_lon = -74.006

    def test_partial_eclipse_magnitude(self):
        """Test that NYC sees partial eclipse with magnitude < 1.0."""
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(
            self.tjd_ut, self.nyc_lat, self.nyc_lon
        )

        # NYC had partial eclipse with ~90% magnitude
        assert 0 < magnitude < 1.0, (
            f"NYC magnitude {magnitude:.4f} should be partial (< 1.0)"
        )

    def test_partial_eclipse_significant_coverage(self):
        """Test that partial eclipse has significant coverage."""
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(
            self.tjd_ut, self.nyc_lat, self.nyc_lon
        )

        # NYC had significant coverage
        assert magnitude > 0.7, f"NYC magnitude {magnitude:.4f} should be > 0.7"


class TestMagnitudeNoEclipse:
    """Test magnitude when no eclipse is happening."""

    def test_no_eclipse_returns_zero(self):
        """Test that function returns 0 when no eclipse."""
        # Random time with no eclipse
        jd = julday(2024, 6, 15, 12.0)
        lat, lon = 32.7767, -96.797

        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)

        assert magnitude == 0.0, "Magnitude should be 0 when no eclipse"

    def test_location_outside_eclipse_path(self):
        """Test that location far from eclipse path returns 0."""
        # April 8, 2024 eclipse - but from Tokyo (no visibility)
        jd = 2460409.28
        tokyo_lat, tokyo_lon = 35.6762, 139.6503

        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, tokyo_lat, tokyo_lon)

        # Tokyo is outside the eclipse path
        assert magnitude == 0.0, "Tokyo should have no eclipse visibility"


class TestConsistencyWithSweEclipseHow:
    """Test consistency with sol_eclipse_how function."""

    def test_magnitude_matches_eclipse_how(self):
        """Test that magnitude matches attr[0] from sol_eclipse_how."""
        jd = 2460409.28
        lat, lon = 32.7767, -96.797
        geopos = [lon, lat, 0]

        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)
        _, attr = sol_eclipse_how(jd, geopos, FLG_SWIEPH)

        # Magnitudes should match closely
        assert abs(magnitude - attr[0]) < 0.01, (
            f"Magnitude {magnitude:.4f} should match sol_eclipse_how {attr[0]:.4f}"
        )

    def test_consistency_at_multiple_times(self):
        """Test consistency at multiple times during eclipse."""
        lat, lon = 32.7767, -96.797
        geopos = [lon, lat, 0]

        # Test at multiple times during eclipse
        test_times = [2460409.26, 2460409.28, 2460409.30]

        for jd in test_times:
            magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)
            _, attr = sol_eclipse_how(jd, geopos, FLG_SWIEPH)

            assert abs(magnitude - attr[0]) < 0.01, (
                f"At JD {jd}: magnitude {magnitude:.4f} differs from "
                f"sol_eclipse_how {attr[0]:.4f}"
            )


class TestSweApiConvention:
    """Test that swe_ version follows the reference ephemeris API conventions."""

    def test_geopos_lon_lat_order(self):
        """Test that geopos uses [lon, lat, alt] order."""
        jd = 2460409.28
        dallas_lat, dallas_lon = 32.7767, -96.797

        # Legacy function: lat, lon order
        result_legacy = _sol_eclipse_magnitude_at_loc_pythonic(jd, dallas_lat, dallas_lon)

        # swe version: geopos = [lon, lat, alt]
        result_swe = sol_eclipse_magnitude_at_loc(
            jd, [dallas_lon, dallas_lat, 0], FLG_SWIEPH
        )

        # Results should be equivalent
        assert abs(result_legacy - result_swe) < 0.0001


class TestEdgeCases:
    """Test edge cases for magnitude calculation."""

    def test_high_altitude_observer(self):
        """Test magnitude at high altitude."""
        jd = 2460409.28
        lat, lon = 32.7767, -96.797

        magnitude_sea_level = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon, altitude=0)
        magnitude_high_alt = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon, altitude=5000)

        # Both should detect eclipse, small difference possible
        assert magnitude_sea_level > 0
        assert magnitude_high_alt > 0
        assert abs(magnitude_sea_level - magnitude_high_alt) < 0.05

    def test_polar_location(self):
        """Test magnitude calculation at polar location."""
        # During eclipse, polar regions may or may not see it
        jd = 2460409.28
        lat, lon = 80.0, 0.0  # Near North Pole

        # Should not crash
        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)
        assert isinstance(magnitude, float)
        assert 0 <= magnitude <= 1.5

    def test_equatorial_location(self):
        """Test magnitude at equatorial location."""
        jd = 2460409.28
        lat, lon = 0.0, 0.0  # Null Island

        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)
        assert isinstance(magnitude, float)
        assert 0 <= magnitude <= 1.5

    def test_accepts_tuple_geopos(self):
        """Test that swe function accepts tuple for geopos."""
        jd = 2460409.28
        geopos = (-96.797, 32.7767, 0)

        result = sol_eclipse_magnitude_at_loc(jd, geopos, FLG_SWIEPH)
        assert isinstance(result, float)

    def test_accepts_list_geopos(self):
        """Test that swe function accepts list for geopos."""
        jd = 2460409.28
        geopos = [-96.797, 32.7767, 0]

        result = sol_eclipse_magnitude_at_loc(jd, geopos, FLG_SWIEPH)
        assert isinstance(result, float)


class TestMagnitudeProgression:
    """Test that magnitude changes appropriately during eclipse."""

    def test_magnitude_increases_to_maximum(self):
        """Test that magnitude increases as eclipse approaches maximum."""
        lat, lon = 32.7767, -96.797

        # Times before maximum (increasing magnitude)
        times_before_max = [2460409.24, 2460409.26, 2460409.27]

        magnitudes = [
            _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon) for jd in times_before_max
        ]

        # Magnitudes should generally increase
        for i in range(1, len(magnitudes)):
            # Allow some noise but trend should be upward
            if magnitudes[i - 1] > 0:  # Only compare if eclipse is visible
                assert magnitudes[i] >= magnitudes[i - 1] * 0.9, (
                    f"Magnitude should increase: {magnitudes}"
                )

    def test_magnitude_decreases_after_maximum(self):
        """Test that magnitude decreases after eclipse maximum."""
        lat, lon = 32.7767, -96.797

        # Times after maximum (decreasing magnitude)
        times_after_max = [2460409.29, 2460409.30, 2460409.31]

        magnitudes = [
            _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon) for jd in times_after_max
        ]

        # Magnitudes should generally decrease
        for i in range(1, len(magnitudes)):
            if magnitudes[i] > 0:  # Only compare if eclipse still visible
                assert magnitudes[i] <= magnitudes[i - 1] * 1.1, (
                    f"Magnitude should decrease: {magnitudes}"
                )


class TestMagnitudeValueAccuracy:
    """Test the accuracy of magnitude values."""

    def test_magnitude_zero_outside_eclipse(self):
        """Test that magnitude is exactly 0 when disks don't overlap."""
        # Time well before any eclipse
        jd = julday(2024, 1, 15, 12.0)
        lat, lon = 32.7767, -96.797

        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)

        assert magnitude == 0.0, "Magnitude should be exactly 0 outside eclipse"

    def test_magnitude_definition(self):
        """Test that magnitude follows the standard definition.

        Eclipse magnitude = (overlap) / (sun_diameter)
        where overlap is the linear distance by which the Moon covers the Sun.
        """
        jd = 2460409.28
        lat, lon = 32.7767, -96.797

        magnitude = _sol_eclipse_magnitude_at_loc_pythonic(jd, lat, lon)

        # For a total eclipse, magnitude >= 1.0 is expected
        # For partial, 0 < magnitude < 1.0
        assert magnitude > 0, "Should detect eclipse at this time/location"
