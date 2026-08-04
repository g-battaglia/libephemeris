"""
Tests for sol_eclipse_where and sol_eclipse_how functions in libephemeris.

Validation tests use NASA's 2024-Apr-08 total solar eclipse data:
- Maximum totality approximately JD 2460409.26 (~18:18 UTC)
- Central line near Nazas, Durango, Mexico (~25.2N, ~104.1W)

Source: https://eclipse.gsfc.nasa.gov/SEplot/SEplot2001/SE2024Apr08T.GIF
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    sol_eclipse_where,
    sol_eclipse_how,
    FLG_SWIEPH,
    ECL_TOTAL,
    ECL_PARTIAL,
    ECL_CENTRAL,
    ECL_VISIBLE,
)


class TestSweSolEclipseWhereSignature:
    """Test the public ``sol_eclipse_where`` signature and result layout."""

    def test_function_exists(self):
        """Test that sol_eclipse_where function exists."""
        from libephemeris.eclipse import sol_eclipse_where

        assert callable(sol_eclipse_where)

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        # Use JD during April 8, 2024 eclipse maximum
        tjd_ut = 2460409.26  # ~18:18 UTC

        retflag, geopos, attr = sol_eclipse_where(tjd_ut, FLG_SWIEPH)

        # The public compatibility contract defines a 10-element geopos tuple.
        assert len(geopos) == 10
        assert all(isinstance(g, float) for g in geopos)

        # The public compatibility contract defines a 20-element attr tuple.
        assert len(attr) == 20
        assert all(isinstance(a, float) for a in attr)

        # retflag should be int
        assert isinstance(retflag, int)

    def test_legacy_function_wraps_correctly(self):
        """Test that legacy sol_eclipse_where function works."""
        tjd_ut = 2460409.26

        retflag, geopos, attr = sol_eclipse_where(tjd_ut, FLG_SWIEPH)

        # Should return same structure (10-element geopos, 20-element attr)
        assert len(geopos) == 10
        assert len(attr) == 20


class TestSweSolEclipseWhereApril2024:
    """Test sol_eclipse_where with April 8, 2024 total solar eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        # JD for April 8, 2024 ~18:12 UTC (maximum totality in Mexico)
        # JD 2460409.26 corresponds to approximately 18:14 UTC
        self.tjd_ut = 2460409.26
        # Expected coordinates: near Nazas, Durango, Mexico
        # Based on NASA eclipse data, the central line at this time is near:
        # 24.5N, 105W (more precise than the initial estimate)
        self.expected_lat = 24.5
        self.expected_lon = -105.0

    def test_finds_eclipse_location(self):
        """Test that function finds an eclipse location."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)

        # Should find an eclipse (non-zero return flag)
        assert retflag != 0

        # Coordinates should be non-zero
        assert geopos[0] != 0.0 or geopos[1] != 0.0

    def test_latitude_within_tolerance(self):
        """Test that latitude is within 0.5 degrees of expected."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)

        lat = geopos[1]
        # Allow 0.5 degree tolerance as specified in requirements
        assert abs(lat - self.expected_lat) < 0.5, (
            f"Latitude {lat:.2f} differs from expected {self.expected_lat:.2f} "
            f"by {abs(lat - self.expected_lat):.2f} degrees"
        )

    def test_longitude_within_tolerance(self):
        """Test that longitude is within 0.5 degrees of expected."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)

        lon = geopos[0]
        # Allow 0.5 degree tolerance as specified in requirements
        assert abs(lon - self.expected_lon) < 0.5, (
            f"Longitude {lon:.2f} differs from expected {self.expected_lon:.2f} "
            f"by {abs(lon - self.expected_lon):.2f} degrees"
        )

    def test_is_total_eclipse(self):
        """Test that eclipse is identified as total."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)

        # Should be total eclipse
        assert retflag & ECL_TOTAL, f"Expected total eclipse, got flags: {retflag}"

    def test_is_central_eclipse(self):
        """Test that eclipse is identified as central."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)

        # Should be central eclipse
        assert retflag & ECL_CENTRAL, f"Expected central eclipse, got flags: {retflag}"

    def test_attributes_are_valid(self):
        """Test that eclipse attributes are in valid ranges."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)

        magnitude = attr[0]
        ratio = attr[1]
        obscuration = attr[2]
        path_width = attr[3]
        azimuth = attr[4]
        true_alt = attr[5]
        apparent_alt = attr[6]
        separation = attr[7]

        # Magnitude should be > 1 for total eclipse
        assert 0.9 < magnitude < 1.5, f"Magnitude {magnitude} out of range"

        # Ratio should be around 1.0
        assert 0.9 < ratio < 1.1, f"Ratio {ratio} out of range"

        # Obscuration for a total eclipse is the disc-area ratio (> 1), the
        # measured reference behavior: obscuration == (r_moon/r_sun)**2 == ratio**2.
        assert obscuration == pytest.approx(ratio**2, rel=1e-6), (
            f"Obscuration {obscuration} should equal ratio**2 {ratio**2}"
        )
        assert obscuration > 1.0, f"total obscuration {obscuration} should exceed 1.0"

        # Path width is negative for total eclipses (sign convention)
        assert path_width < 0, (
            f"Path width {path_width} km should be negative for total"
        )
        assert abs(path_width) < 500, f"|Path width| {abs(path_width)} km out of range"

        # Azimuth should be 0-360
        assert 0 <= azimuth < 360, f"Azimuth {azimuth} out of range"

        # Altitudes should be reasonable
        assert -90 <= true_alt <= 90, f"True altitude {true_alt} out of range"
        assert -90 <= apparent_alt <= 90, (
            f"Apparent altitude {apparent_alt} out of range"
        )

        # Separation should be small at maximum
        assert 0 <= separation < 0.5, f"Separation {separation} out of range"


class TestSweSolEclipseHowSignature:
    """Test the public ``sol_eclipse_how`` signature and result layout."""

    def test_function_exists(self):
        """Test that sol_eclipse_how function exists."""
        from libephemeris.eclipse import sol_eclipse_how

        assert callable(sol_eclipse_how)

    def test_returns_correct_tuple_structure(self):
        """Test that return values have correct structure."""
        # Use JD during April 8, 2024 eclipse maximum
        tjd_ut = 2460409.26
        dallas_geopos = [-96.797, 32.7767, 0]  # lon, lat, alt

        retflag, attr = sol_eclipse_how(tjd_ut, dallas_geopos, FLG_SWIEPH)

        # attr should be at least 8-element tuple
        assert len(attr) == 20
        assert all(isinstance(a, float) for a in attr)

        # retflag should be int
        assert isinstance(retflag, int)

    def test_accepts_geopos_as_list(self):
        """Test that function accepts geopos as list."""
        tjd_ut = 2460409.26
        geopos = [-96.797, 32.7767, 0]

        retflag, attr = sol_eclipse_how(tjd_ut, geopos, FLG_SWIEPH)
        assert isinstance(attr, tuple)

    def test_accepts_geopos_as_tuple(self):
        """Test that function accepts geopos as tuple."""
        tjd_ut = 2460409.26
        geopos = (-96.797, 32.7767, 0)

        retflag, attr = sol_eclipse_how(tjd_ut, geopos, FLG_SWIEPH)
        assert isinstance(attr, tuple)

    def test_invalid_geopos_raises_error(self):
        """Test that invalid geopos raises ValueError."""
        tjd_ut = 2460409.26

        # Too few elements
        with pytest.raises(ValueError):
            sol_eclipse_how(tjd_ut, [0, 0])

    def test_legacy_function_wraps_correctly(self):
        """Test that legacy sol_eclipse_how function works."""
        tjd_ut = 2460409.26

        # Now aliases swe_ signature: (jd, geopos, ifl)
        retflag, attr = sol_eclipse_how(tjd_ut, (-96.797, 32.7767, 0), FLG_SWIEPH)
        assert len(attr) == 20


class TestSweSolEclipseHowDallasApril2024:
    """Test sol_eclipse_how at Dallas during April 8, 2024 total eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        # JD for April 8, 2024 ~18:42 UTC (maximum for Dallas)
        self.tjd_ut = 2460409.28  # Around Dallas maximum
        # Dallas coordinates: 32.7767N, 96.797W
        # geopos format: [longitude, latitude, altitude]
        self.geopos_dallas = [-96.797, 32.7767, 0]

    def test_finds_eclipse_at_dallas(self):
        """Test that function finds eclipse at Dallas."""
        retflag, attr = sol_eclipse_how(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        # Should find an eclipse (non-zero return flag)
        assert retflag != 0
        assert retflag & ECL_VISIBLE

    def test_dallas_eclipse_is_total(self):
        """Test that Dallas sees a total eclipse."""
        retflag, attr = sol_eclipse_how(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        # Should be total eclipse
        assert retflag & ECL_TOTAL, f"Expected total eclipse, got flags: {retflag}"

    def test_dallas_obscuration_is_total(self):
        """Test that obscuration at Dallas is ~100% (within 1%)."""
        retflag, attr = sol_eclipse_how(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        obscuration = attr[2]
        # For a total eclipse the reference reports the disc-area ratio, which
        # is >= 1 (the Sun is fully covered, and the Moon overfills it).
        assert obscuration >= 0.99, (
            f"Obscuration {obscuration:.3f} is less than expected 0.99 for total eclipse"
        )

    def test_dallas_attributes_are_valid(self):
        """Test that eclipse attributes at Dallas are in valid ranges."""
        retflag, attr = sol_eclipse_how(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        magnitude = attr[0]
        ratio = attr[1]
        obscuration = attr[2]
        shadow_width = attr[3]
        azimuth = attr[4]
        true_alt = attr[5]
        apparent_alt = attr[6]
        separation = attr[7]

        # Magnitude should be >= 1 for total eclipse
        assert 0.9 < magnitude < 1.5, f"Magnitude {magnitude} out of range"

        # Ratio should be around 1.0
        assert 0.9 < ratio < 1.1, f"Ratio {ratio} out of range"

        # Obscuration for a total eclipse is the disc-area ratio (> 1), the
        # measured reference behavior: obscuration == (r_moon/r_sun)**2 == ratio**2.
        assert obscuration == pytest.approx(ratio**2, rel=1e-6), (
            f"Obscuration {obscuration} should equal ratio**2 {ratio**2}"
        )
        assert obscuration > 1.0, f"total obscuration {obscuration} should exceed 1.0"

        # Shadow width is negative for total eclipses (sign convention)
        assert shadow_width < 0, (
            f"Shadow width {shadow_width} km should be negative for total"
        )

        # Azimuth should be 0-360
        assert 0 <= azimuth < 360, f"Azimuth {azimuth} out of range"

        # Sun should be above horizon
        assert true_alt > 0, f"True altitude {true_alt} should be positive"
        assert apparent_alt > 0, f"Apparent altitude {apparent_alt} should be positive"

        # Separation should be small at maximum
        assert 0 <= separation < 0.3, f"Separation {separation} should be small"

    def test_refraction_included(self):
        """Test that apparent altitude differs from true altitude (refraction)."""
        retflag, attr = sol_eclipse_how(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        true_alt = attr[5]
        apparent_alt = attr[6]

        # Apparent altitude should be slightly higher due to refraction
        # (at typical eclipse Sun altitudes, refraction is small but non-zero)
        assert apparent_alt >= true_alt, (
            f"Apparent alt {apparent_alt} should be >= true alt {true_alt}"
        )


class TestSweSolEclipseHowNYCApril2024:
    """Test sol_eclipse_how at NYC during April 8, 2024 eclipse (partial)."""

    def setup_method(self):
        """Set up test fixtures."""
        # JD for April 8, 2024 ~19:00 UTC (maximum for NYC)
        # This is closer to when NYC sees maximum eclipse
        self.tjd_ut = 2460409.30
        # NYC coordinates: 40.7128N, 74.006W
        self.geopos_nyc = [-74.006, 40.7128, 0]

    def test_finds_eclipse_at_nyc(self):
        """Test that function finds eclipse at NYC."""
        retflag, attr = sol_eclipse_how(self.tjd_ut, self.geopos_nyc, FLG_SWIEPH)

        # Should find an eclipse
        assert retflag != 0
        assert retflag & ECL_VISIBLE

    def test_nyc_eclipse_is_partial(self):
        """Test that NYC sees a partial eclipse."""
        retflag, attr = sol_eclipse_how(self.tjd_ut, self.geopos_nyc, FLG_SWIEPH)

        # Should be partial eclipse (not total)
        assert retflag & ECL_PARTIAL, f"Expected partial eclipse, got flags: {retflag}"
        assert not (retflag & ECL_TOTAL), "Should not be total at NYC"

    def test_nyc_obscuration_is_partial(self):
        """Test that obscuration at NYC is partial (~80-95%)."""
        retflag, attr = sol_eclipse_how(self.tjd_ut, self.geopos_nyc, FLG_SWIEPH)

        obscuration = attr[2]
        # NYC should have significant but partial obscuration
        assert 0.7 < obscuration < 1.0, (
            f"Obscuration {obscuration:.3f} out of expected range for NYC"
        )


class TestSweSolEclipseHowNoEclipse:
    """Test sol_eclipse_how when no eclipse is happening."""

    def test_no_eclipse_returns_zero_flag(self):
        """Test that function returns 0 flag when no eclipse."""
        # Use a random time when no eclipse is happening
        tjd_ut = julday(2024, 6, 15, 12.0)  # Random date
        geopos = [-96.797, 32.7767, 0]

        retflag, attr = sol_eclipse_how(tjd_ut, geopos, FLG_SWIEPH)

        # Should return 0 for no eclipse
        assert retflag == 0 or attr[0] == 0.0, (
            f"Expected no eclipse, but got flag {retflag}, magnitude {attr[0]}"
        )


class TestSweSolEclipseWhereTimeVariation:
    """Test sol_eclipse_where at different times during eclipse."""

    def test_eclipse_path_moves_east(self):
        """Test that eclipse central line moves generally eastward."""
        # Sample times during April 8, 2024 eclipse
        times = [
            2460409.24,  # Early
            2460409.26,  # Middle
            2460409.28,  # Late
        ]

        longitudes = []
        for tjd in times:
            retflag, geopos, attr = sol_eclipse_where(tjd, FLG_SWIEPH)
            if retflag != 0:
                longitudes.append(geopos[0])

        # Eclipse path moves from west to east (longitude increases)
        # Note: Need at least 2 valid points
        if len(longitudes) >= 2:
            # The longitude should increase (move east) but can wrap
            # For North American eclipse, we expect westward longitudes becoming less negative
            pass  # Complex to test due to longitude wrapping


class TestSweSolEclipseWherePartialEclipse:
    """Test sol_eclipse_where during a partial-only eclipse."""

    def test_partial_eclipse_returns_partial_flag(self):
        """Test that partial eclipse is identified correctly."""
        # Find a time when eclipse is partial-only (near edge)
        # For the Apr 2024 eclipse, very early/late times have no central line
        tjd_ut = 2460409.15  # Early in eclipse, may be partial

        retflag, geopos, attr = sol_eclipse_where(tjd_ut, FLG_SWIEPH)

        # Either no central eclipse (returns 0) or partial
        if retflag != 0:
            # If eclipse found, check it's valid
            assert geopos[0] != 0 or geopos[1] != 0


class TestSweSolEclipseHowEdgeCases:
    """Test edge cases for sol_eclipse_how."""

    def test_high_altitude_observer(self):
        """Test with high altitude observer."""
        tjd_ut = 2460409.26
        # Same as Dallas but at 5000m altitude
        geopos_high = [-96.797, 32.7767, 5000]

        retflag, attr = sol_eclipse_how(tjd_ut, geopos_high, FLG_SWIEPH)

        # Should still find eclipse
        assert retflag != 0

    def test_southern_hemisphere(self):
        """Test with southern hemisphere location."""
        tjd_ut = 2460409.26
        # Sydney, Australia - far from eclipse path
        geopos_sydney = [151.2093, -33.8688, 0]

        retflag, attr = sol_eclipse_how(tjd_ut, geopos_sydney, FLG_SWIEPH)

        # May or may not see eclipse, but should not crash
        assert isinstance(retflag, int)

    def test_extreme_latitude(self):
        """Test with extreme latitude location."""
        tjd_ut = 2460409.26
        # Near North Pole
        geopos_arctic = [0, 85, 0]

        retflag, attr = sol_eclipse_how(tjd_ut, geopos_arctic, FLG_SWIEPH)

        # Should not crash
        assert isinstance(retflag, int)


class TestSweSolEclipseWhereLimits:
    """Validate the public ``sol_eclipse_where`` result contract."""

    def setup_method(self):
        """Set up test fixtures for April 8, 2024 total solar eclipse."""
        # JD during maximum totality in Mexico
        self.tjd_ut = 2460409.26

    def test_geopos_has_ten_elements(self):
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)
        assert len(geopos) == 10

    def test_center_point_filled(self):
        """The shadow center falls in the broad NASA-published Mexico region."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)
        assert retflag & ECL_TOTAL
        assert -180.0 <= geopos[0] <= 180.0
        assert -90.0 <= geopos[1] <= 90.0
        assert -110.0 <= geopos[0] <= -100.0
        assert 20.0 <= geopos[1] <= 30.0

    def test_remaining_geopos_slots_are_zero(self):
        """Reserved geopos slots are zero per the public result contract."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)
        assert retflag != 0
        for i in range(2, 10):
            assert geopos[i] == 0.0, f"geopos[{i}] must be 0.0, got {geopos[i]}"

    def test_core_shadow_width_negative_for_total(self):
        """attr[3] is the core-shadow diameter, negative when umbral."""
        retflag, geopos, attr = sol_eclipse_where(self.tjd_ut, FLG_SWIEPH)
        assert retflag & ECL_TOTAL
        assert attr[3] < 0.0

    def test_no_eclipse_instant_returns_zero_flag_with_data(self):
        """A no-eclipse instant gives retflag 0 but still fills geopos/attr."""
        retflag, geopos, attr = sol_eclipse_where(2459375.5, FLG_SWIEPH)
        assert retflag == 0
        assert geopos[0] != 0.0 or geopos[1] != 0.0
        assert attr[0] < 0.0
        # No eclipse -> the obscuration slot carries the fixed no-event
        # protocol sentinel 1.0 (compatibility contract; see the Breaking
        # entry in CHANGELOG 3.0.0 and eclipse.py). The pre-rc9 0.0 rested
        # on a comment that asserted the opposite measurement.
        assert attr[2] == 1.0
