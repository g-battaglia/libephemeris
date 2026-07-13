"""
Tests for sol_eclipse_how_details and sol_eclipse_how_details functions.

These tests validate the comprehensive eclipse circumstances calculation
including contact times, position angles, obscuration percentages, and
Sun altitude/azimuth during eclipse phases.

Reference eclipse: April 8, 2024 total solar eclipse
- Dallas, TX: Total eclipse with ~4 minutes of totality
- NYC: Partial eclipse with ~90% obscuration
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris import (
    julday,
    sol_eclipse_how_details,
    _sol_eclipse_how_details_pythonic,
    sol_eclipse_how,
    sol_eclipse_when_loc,
    FLG_SWIEPH,
)


class TestSweSolEclipseHowDetailsSignature:
    """Test that sol_eclipse_how_details function signature and return structure."""

    def test_function_exists(self):
        """Test that sol_eclipse_how_details function exists."""
        from libephemeris.eclipse import sol_eclipse_how_details

        assert callable(sol_eclipse_how_details)

    def test_exported_from_package(self):
        """Test that function is exported from main package."""
        from libephemeris import sol_eclipse_how_details

        assert callable(sol_eclipse_how_details)

    def test_returns_dict(self):
        """Test that function returns a dictionary."""
        tjd_ut = 2460409.28  # April 8, 2024 eclipse
        dallas_geopos = [-96.797, 32.7767, 0]

        result = sol_eclipse_how_details(tjd_ut, dallas_geopos, FLG_SWIEPH)

        assert isinstance(result, dict)

    def test_dict_has_required_keys(self):
        """Test that returned dict has all required keys."""
        tjd_ut = 2460409.28
        dallas_geopos = [-96.797, 32.7767, 0]

        result = sol_eclipse_how_details(tjd_ut, dallas_geopos, FLG_SWIEPH)

        required_keys = [
            # Eclipse type
            "eclipse_type",
            "is_visible",
            "is_total",
            "is_annular",
            "is_partial",
            # Contact times
            "jd_c1",
            "jd_c2",
            "jd_max",
            "jd_c3",
            "jd_c4",
            # Magnitude and obscuration
            "magnitude",
            "max_magnitude",
            "obscuration",
            "max_obscuration",
            "obscuration_percent",
            "max_obscuration_percent",
            "ratio",
            "shadow_width_km",
            # Position angles
            "position_angle_c1",
            "position_angle_c2",
            "position_angle_c3",
            "position_angle_c4",
            # Sun positions
            "sun_alt_c1",
            "sun_az_c1",
            "sun_alt_c2",
            "sun_az_c2",
            "sun_alt_max",
            "sun_az_max",
            "sun_alt_c3",
            "sun_az_c3",
            "sun_alt_c4",
            "sun_az_c4",
            # Durations
            "duration_partial_minutes",
            "duration_total_minutes",
            # Angular sizes
            "sun_angular_radius",
            "moon_angular_radius",
            "separation",
        ]

        for key in required_keys:
            assert key in result, f"Missing key: {key}"

    def test_invalid_geopos_raises_error(self):
        """Test that invalid geopos raises ValueError."""
        tjd_ut = 2460409.28

        with pytest.raises(ValueError):
            sol_eclipse_how_details(tjd_ut, [0, 0], FLG_SWIEPH)


class TestSweSolEclipseHowDetailsDallasApril2024:
    """Test sol_eclipse_how_details at Dallas during April 8, 2024 total eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        # JD for April 8, 2024 ~18:42 UTC (maximum for Dallas area)
        self.tjd_ut = 2460409.28
        # Dallas coordinates: 32.7767N, 96.797W
        self.geopos_dallas = [-96.797, 32.7767, 0]

    def test_identifies_total_eclipse(self):
        """Test that Dallas sees a total eclipse."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        assert result["is_visible"] is True
        assert result["is_total"] is True
        assert result["is_annular"] is False

    def test_contact_times_are_valid(self):
        """Test that all four contact times are returned."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        # All contact times should be non-zero for total eclipse
        assert result["jd_c1"] > 0, "First contact time should be > 0"
        assert result["jd_c2"] > 0, "Second contact time should be > 0"
        assert result["jd_c3"] > 0, "Third contact time should be > 0"
        assert result["jd_c4"] > 0, "Fourth contact time should be > 0"

    def test_contact_times_in_correct_order(self):
        """Test that contact times are in chronological order."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        assert result["jd_c1"] < result["jd_c2"], "C1 should be before C2"
        assert result["jd_c2"] < result["jd_max"], "C2 should be before max"
        assert result["jd_max"] < result["jd_c3"], "Max should be before C3"
        assert result["jd_c3"] < result["jd_c4"], "C3 should be before C4"

    def test_maximum_obscuration_is_total(self):
        """Test that maximum obscuration is ~100% for total eclipse."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        # For total eclipse, max obscuration should be 100%
        assert result["max_obscuration"] >= 0.99, (
            f"Max obscuration {result['max_obscuration']:.3f} should be >= 0.99"
        )
        assert result["max_obscuration_percent"] >= 99.0

    def test_totality_duration_reasonable(self):
        """Test that totality duration is reasonable for Dallas."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        # Dallas had approximately 3-4 minutes of totality
        duration = result["duration_total_minutes"]
        assert 1.0 < duration < 6.0, (
            f"Totality duration {duration:.1f} min seems unreasonable"
        )

    def test_partial_phase_duration_reasonable(self):
        """Test that partial phase duration is reasonable."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        # Total eclipse duration (C1 to C4) typically 2-3 hours
        duration = result["duration_partial_minutes"]
        assert 60 < duration < 200, (
            f"Partial duration {duration:.1f} min seems unreasonable"
        )

    def test_sun_above_horizon_at_all_contacts(self):
        """Test that Sun is above horizon at all contacts."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        assert result["sun_alt_c1"] > 0, "Sun should be above horizon at C1"
        assert result["sun_alt_c2"] > 0, "Sun should be above horizon at C2"
        assert result["sun_alt_max"] > 0, "Sun should be above horizon at max"
        assert result["sun_alt_c3"] > 0, "Sun should be above horizon at C3"
        assert result["sun_alt_c4"] > 0, "Sun should be above horizon at C4"

    def test_position_angles_in_valid_range(self):
        """Test that position angles are in 0-360 degree range."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        for key in [
            "position_angle_c1",
            "position_angle_c2",
            "position_angle_c3",
            "position_angle_c4",
        ]:
            pa = result[key]
            if pa > 0:  # Only check non-zero values
                assert 0 <= pa <= 360, f"{key} = {pa} out of range"

    def test_azimuths_in_valid_range(self):
        """Test that azimuths are in 0-360 degree range."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        for key in ["sun_az_c1", "sun_az_c2", "sun_az_max", "sun_az_c3", "sun_az_c4"]:
            az = result[key]
            if az != 0:  # Only check non-zero values
                assert 0 <= az < 360, f"{key} = {az} out of range"

    def test_angular_sizes_reasonable(self):
        """Test that angular sizes are in reasonable range."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        # Sun angular radius ~0.26-0.27 degrees (~16 arcmin)
        assert 0.25 < result["sun_angular_radius"] < 0.28

        # Moon angular radius varies more, ~0.24-0.29 degrees
        assert 0.23 < result["moon_angular_radius"] < 0.30


class TestSweSolEclipseHowDetailsNYCApril2024:
    """Test sol_eclipse_how_details at NYC during April 8, 2024 partial eclipse."""

    def setup_method(self):
        """Set up test fixtures."""
        # JD for April 8, 2024 ~19:00 UTC (near maximum for NYC)
        self.tjd_ut = 2460409.30
        # NYC coordinates: 40.7128N, 74.006W
        self.geopos_nyc = [-74.006, 40.7128, 0]

    def test_identifies_partial_eclipse(self):
        """Test that NYC sees a partial eclipse."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_nyc, FLG_SWIEPH)

        assert result["is_visible"] is True
        assert result["is_partial"] is True
        assert result["is_total"] is False

    def test_c2_c3_are_zero_for_partial(self):
        """Test that C2 and C3 are zero for partial eclipse."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_nyc, FLG_SWIEPH)

        # For partial eclipse, C2 and C3 don't exist
        assert result["jd_c2"] == 0.0
        assert result["jd_c3"] == 0.0

    def test_c1_c4_exist_for_partial(self):
        """Test that C1 and C4 exist for partial eclipse."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_nyc, FLG_SWIEPH)

        assert result["jd_c1"] > 0
        assert result["jd_c4"] > 0
        assert result["jd_c1"] < result["jd_max"] < result["jd_c4"]

    def test_obscuration_is_partial(self):
        """Test that obscuration is partial (not 100%)."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_nyc, FLG_SWIEPH)

        # NYC should have significant but not total obscuration
        assert 0.7 < result["max_obscuration"] < 1.0
        assert result["max_obscuration_percent"] < 100.0

    def test_totality_duration_is_zero(self):
        """Test that totality duration is zero for partial eclipse."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_nyc, FLG_SWIEPH)

        assert result["duration_total_minutes"] == 0.0


class TestSweSolEclipseHowDetailsNoEclipse:
    """Test sol_eclipse_how_details when no eclipse is happening."""

    def test_no_eclipse_returns_zero_type(self):
        """Test that function returns zero eclipse type when no eclipse."""
        # Random time with no eclipse
        tjd_ut = julday(2024, 6, 15, 12.0)
        geopos = [-96.797, 32.7767, 0]

        result = sol_eclipse_how_details(tjd_ut, geopos, FLG_SWIEPH)

        # Should indicate no visible eclipse
        assert result["is_visible"] is False
        assert result["magnitude"] == 0.0 or result["eclipse_type"] == 0


class TestSolEclipseHowDetailsLegacy:
    """Test legacy sol_eclipse_how_details function."""

    def test_function_exists(self):
        """Test that legacy function exists."""
        from libephemeris.eclipse import sol_eclipse_how_details

        assert callable(sol_eclipse_how_details)

    def test_legacy_wraps_correctly(self):
        """Test that legacy function wraps swe version correctly."""
        tjd_ut = 2460409.28

        # Legacy: lat, lon order (pythonic helper)
        result_legacy = _sol_eclipse_how_details_pythonic(
            tjd_ut, 32.7767, -96.797, 0, FLG_SWIEPH
        )

        # reference-API style: geopos = [lon, lat, alt]
        result_swe = sol_eclipse_how_details(tjd_ut, [-96.797, 32.7767, 0], FLG_SWIEPH)

        # Results should be equivalent
        assert result_legacy["is_total"] == result_swe["is_total"]
        assert (
            abs(result_legacy["max_obscuration"] - result_swe["max_obscuration"]) < 0.01
        )


class TestSweSolEclipseHowDetailsConsistency:
    """Test consistency with sol_eclipse_how."""

    def test_magnitude_matches_base_function(self):
        """Test that magnitude matches sol_eclipse_how."""
        tjd_ut = 2460409.28
        geopos = [-96.797, 32.7767, 0]

        details = sol_eclipse_how_details(tjd_ut, geopos, FLG_SWIEPH)
        _, attr = sol_eclipse_how(tjd_ut, geopos, FLG_SWIEPH)

        # Magnitudes should be close (might differ slightly due to different JD)
        assert abs(details["magnitude"] - attr[0]) < 0.1

    def test_obscuration_matches_base_function(self):
        """Test that obscuration matches sol_eclipse_how."""
        tjd_ut = 2460409.28
        geopos = [-96.797, 32.7767, 0]

        details = sol_eclipse_how_details(tjd_ut, geopos, FLG_SWIEPH)
        _, attr = sol_eclipse_how(tjd_ut, geopos, FLG_SWIEPH)

        # Obscurations should be close
        assert abs(details["obscuration"] - attr[2]) < 0.1


class TestSweSolEclipseHowDetailsEdgeCases:
    """Test edge cases for sol_eclipse_how_details."""

    def test_high_altitude_observer(self):
        """Test with high altitude observer."""
        tjd_ut = 2460409.28
        geopos_high = [-96.797, 32.7767, 5000]  # 5000m altitude

        result = sol_eclipse_how_details(tjd_ut, geopos_high, FLG_SWIEPH)

        # Should still work and find eclipse
        assert result["is_visible"] is True

    def test_near_eclipse_boundary(self):
        """Test location near eclipse visibility boundary."""
        tjd_ut = 2460409.28
        # Location far from central line
        geopos_edge = [-70.0, 50.0, 0]  # Eastern Canada

        result = sol_eclipse_how_details(tjd_ut, geopos_edge, FLG_SWIEPH)

        # Should return valid data without crashing
        assert isinstance(result, dict)

    def test_accepts_tuple_geopos(self):
        """Test that function accepts geopos as tuple."""
        tjd_ut = 2460409.28
        geopos = (-96.797, 32.7767, 0)

        result = sol_eclipse_how_details(tjd_ut, geopos, FLG_SWIEPH)

        assert isinstance(result, dict)

    def test_accepts_list_geopos(self):
        """Test that function accepts geopos as list."""
        tjd_ut = 2460409.28
        geopos = [-96.797, 32.7767, 0]

        result = sol_eclipse_how_details(tjd_ut, geopos, FLG_SWIEPH)

        assert isinstance(result, dict)


class TestSweSolEclipseHowDetailsContactTimes:
    """Detailed tests for contact time calculations."""

    def setup_method(self):
        """Set up test fixtures."""
        self.tjd_ut = 2460409.28
        self.geopos_dallas = [-96.797, 32.7767, 0]

    def test_contact_times_close_to_when_loc(self):
        """Test that contact times are close to sol_eclipse_when_loc results."""
        # Get contact times from when_loc
        jd_start = julday(2024, 1, 1, 0)
        ecl_type, tret, attr = sol_eclipse_when_loc(
            jd_start, self.geopos_dallas, FLG_SWIEPH
        )

        # Get details at the time of the eclipse
        result = sol_eclipse_how_details(tret[0], self.geopos_dallas, FLG_SWIEPH)

        # Contact times should be reasonably close (within 5 minutes)
        tolerance_jd = 5.0 / (24 * 60)  # 5 minutes in JD

        if tret[1] > 0 and result["jd_c1"] > 0:
            assert abs(result["jd_c1"] - tret[1]) < tolerance_jd, (
                f"C1 differs by {abs(result['jd_c1'] - tret[1]) * 24 * 60:.1f} minutes"
            )

        if tret[4] > 0 and result["jd_c4"] > 0:
            assert abs(result["jd_c4"] - tret[4]) < tolerance_jd, (
                f"C4 differs by {abs(result['jd_c4'] - tret[4]) * 24 * 60:.1f} minutes"
            )

    def test_c2_c3_timing_symmetric(self):
        """Test that C2 and C3 are roughly symmetric around maximum."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        if result["jd_c2"] > 0 and result["jd_c3"] > 0:
            # Time from C2 to max and max to C3 should be similar
            time_c2_max = result["jd_max"] - result["jd_c2"]
            time_max_c3 = result["jd_c3"] - result["jd_max"]

            # Should be within ~30 seconds of each other
            assert abs(time_c2_max - time_max_c3) < 0.5 / (24 * 60), (
                "C2/C3 should be roughly symmetric around maximum"
            )


class TestSweSolEclipseHowDetailsPositionAngles:
    """Tests for position angle calculations."""

    def setup_method(self):
        """Set up test fixtures."""
        self.tjd_ut = 2460409.28
        self.geopos_dallas = [-96.797, 32.7767, 0]

    def test_position_angles_complementary(self):
        """Test that C1 and C4 position angles are roughly opposite."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        pa_c1 = result["position_angle_c1"]
        pa_c4 = result["position_angle_c4"]

        if pa_c1 > 0 and pa_c4 > 0:
            # C1 and C4 should be on opposite sides of the Sun
            # The difference should be roughly 180 degrees
            diff = abs(pa_c4 - pa_c1)
            if diff > 180:
                diff = 360 - diff

            # Allow some tolerance (within 90 degrees of being opposite)
            assert diff > 90, (
                f"C1 ({pa_c1:.1f}°) and C4 ({pa_c4:.1f}°) should be roughly opposite"
            )

    def test_c2_c3_position_angles(self):
        """Test that C2 and C3 position angles are calculated for total eclipse."""
        result = sol_eclipse_how_details(self.tjd_ut, self.geopos_dallas, FLG_SWIEPH)

        if result["is_total"]:
            # For total eclipse, C2 and C3 should have position angles
            assert result["position_angle_c2"] > 0 or result["jd_c2"] == 0
            assert result["position_angle_c3"] > 0 or result["jd_c3"] == 0
