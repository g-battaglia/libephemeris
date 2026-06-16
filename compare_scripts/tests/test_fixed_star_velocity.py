"""
Unit tests for fixed star velocity calculations.

Tests verify that:
1. FLG_SPEED flag enables velocity computation
2. Velocities match pyswisseph within 10%
3. Without FLG_SPEED, velocities remain zero (backward compatibility)
4. Speed is approximately the precession rate (~50.3 arcsec/year)
"""

import pytest
import swisseph as swe

import libephemeris as ephem
from libephemeris.constants import FLG_SPEED, FLG_SWIEPH
from libephemeris.fixed_stars import (
    calc_fixed_star_velocity,
    fixstar_ut,
    fixstar,
    fixstar2_ut,
    fixstar2,
)


@pytest.mark.unit
class TestFixedStarVelocity:
    """Tests for fixed star velocity calculations."""

    def test_velocity_calculation_regulus(self):
        """Test velocity calculation for Regulus."""
        from libephemeris.constants import REGULUS

        # J2000.0 epoch
        jd_tt = 2451545.0

        lon, lat, dist, speed_lon, speed_lat, speed_dist = calc_fixed_star_velocity(
            REGULUS, jd_tt
        )

        # Position should be valid
        assert 149 < lon < 151, f"Regulus lon: {lon:.4f}"
        assert -1 < lat < 2, f"Regulus lat: {lat:.4f}"
        assert dist > 1000, f"Regulus dist: {dist}"

        # Velocity should be around precession rate: ~50.3 arcsec/year
        # = 50.3 / 3600 / 365.25 = 0.0000378 deg/day
        # But actual value includes proper motion and nutation, so expect ~0.0001 deg/day range
        assert 0.00001 < abs(speed_lon) < 0.001, (
            f"Unexpected speed_lon: {speed_lon:.6f} deg/day"
        )

        # Latitude velocity should be small
        assert abs(speed_lat) < 0.0001, f"Unexpected speed_lat: {speed_lat:.6f} deg/day"

        # Distance "velocity" is a small numerical wobble of the catalogue
        # distance (parallax/aberration over the finite-difference step), not
        # exactly zero -- pyswisseph returns the same ~0.01-0.04 AU/day.
        assert abs(speed_dist) < 0.1, f"speed_dist unexpectedly large: {speed_dist}"

    def test_swe_fixstar_ut_without_speed_flag(self, standard_jd):
        """Test that without FLG_SPEED, velocities are zero."""
        pos, name, retflag = fixstar_ut("Regulus", standard_jd, 0)

        # Velocities should all be zero
        assert pos[3] == 0.0, f"speed_lon should be 0, got {pos[3]}"
        assert pos[4] == 0.0, f"speed_lat should be 0, got {pos[4]}"
        assert pos[5] == 0.0, f"speed_dist should be 0, got {pos[5]}"

        # Position should still be valid
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.4f}"
        assert "Regulus" in name, f"Unexpected name: {name}"

    def test_swe_fixstar_ut_with_speed_flag(self, standard_jd):
        """Test that with FLG_SPEED, velocities are computed."""
        pos, name, retflag = fixstar_ut("Regulus", standard_jd, FLG_SPEED)

        # Velocities should be non-zero
        assert pos[3] != 0.0, f"speed_lon should be non-zero, got {pos[3]}"

        # Position should still be valid
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.4f}"
        assert "Regulus" in name, f"Unexpected name: {name}"

    def test_swe_fixstar_with_speed_flag(self, standard_jd):
        """Test fixstar (TT version) with FLG_SPEED."""
        pos, name, retflag = fixstar("Regulus", standard_jd, FLG_SPEED)

        # Velocities should be non-zero
        assert pos[3] != 0.0, f"speed_lon should be non-zero, got {pos[3]}"

        # Position should still be valid
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.4f}"

    def test_swe_fixstar2_ut_with_speed_flag(self, standard_jd):
        """Test fixstar2_ut with FLG_SPEED."""
        pos, name, retflag = fixstar2_ut("Regulus", standard_jd, FLG_SPEED)

        # Name should be returned
        assert "Regulus" in name, f"Expected Regulus in name, got {name}"

        # Velocities should be non-zero
        assert pos[3] != 0.0, f"speed_lon should be non-zero, got {pos[3]}"

        # Position should still be valid
        assert 149 < pos[0] < 151, f"Regulus lon: {pos[0]:.4f}"

    def test_swe_fixstar2_with_speed_flag(self, standard_jd):
        """Test fixstar2 (TT version) with FLG_SPEED."""
        pos, name, retflag = fixstar2("Regulus", standard_jd, FLG_SPEED)

        # Name should be returned
        assert "Regulus" in name, f"Expected Regulus in name, got {name}"

        # Velocities should be non-zero
        assert pos[3] != 0.0, f"speed_lon should be non-zero, got {pos[3]}"

    @pytest.mark.parametrize(
        "star_name",
        [
            "Regulus",
            "Spica",
            "Sirius",
            "Aldebaran",
            "Antares",
        ],
    )
    def test_velocity_for_multiple_stars(self, standard_jd, star_name):
        """Test velocity calculation for multiple stars."""
        pos, name, retflag = fixstar_ut(star_name, standard_jd, FLG_SPEED)

        # Longitude velocity should be primarily due to precession and aberration
        # ~50.3 arcsec/year = 0.0000378 deg/day for precession
        # Plus proper motion, nutation, and aberration effects
        # Aberration can cause near-zero velocities at specific epochs
        # when its rate of change cancels other effects
        assert abs(pos[3]) < 0.001, f"{star_name} unexpected speed_lon: {pos[3]:.6f}"

        # Latitude velocity should be small
        assert abs(pos[4]) < 0.0001, f"{star_name} unexpected speed_lat: {pos[4]:.6f}"

        # Distance "velocity" is a small numerical wobble (parallax/aberration
        # over the finite-difference step), not exactly zero -- as in pyswisseph.
        assert abs(pos[5]) < 0.1, f"{star_name} speed_dist unexpectedly large: {pos[5]}"


@pytest.mark.unit
class TestFixedStarVelocityVsPyswisseph:
    """Tests comparing velocity with pyswisseph."""

    def test_regulus_velocity_vs_pyswisseph(self, standard_jd):
        """Compare Regulus velocity with pyswisseph - sign and magnitude.

        Note: Latitude velocity magnitude differs due to different calculation
        methods (finite difference vs. SE's analytical proper motion transformation),
        but signs now match after implementing SE-compatible velocity sign logic.
        """
        # Get libephemeris velocity
        pos_lib, name, retflag = fixstar_ut(
            "Regulus", standard_jd, FLG_SWIEPH | FLG_SPEED
        )

        # Check if pyswisseph has fixstar_ut
        if not hasattr(swe, "fixstar_ut"):
            pytest.skip("pyswisseph does not have fixstar_ut")

        # Get pyswisseph velocity
        try:
            pos_swe, name_swe, err_swe = swe.fixstar_ut(
                "Regulus", standard_jd, swe.FLG_SWIEPH | swe.FLG_SPEED
            )
        except Exception as e:
            pytest.skip(f"pyswisseph fixstar_ut failed: {e}")

        # Compare longitude velocity - should match within 10%
        speed_lon_lib = pos_lib[3]
        speed_lon_swe = pos_swe[3]

        # Calculate percentage difference
        if abs(speed_lon_swe) > 1e-10:
            pct_diff = abs((speed_lon_lib - speed_lon_swe) / speed_lon_swe) * 100
            assert pct_diff < 10, (
                f"speed_lon diff {pct_diff:.1f}%: lib={speed_lon_lib:.8f}, swe={speed_lon_swe:.8f}"
            )
        else:
            # If pyswisseph returns near-zero, check absolute difference
            assert abs(speed_lon_lib - speed_lon_swe) < 0.0001

        # Check latitude velocity sign matches (magnitude may differ due to
        # different calculation methods - SE uses analytical proper motion
        # transformation while we use finite differences)
        speed_lat_lib = pos_lib[4]
        speed_lat_swe = pos_swe[4]

        # Signs must match (unless both are near-zero)
        if abs(speed_lat_lib) > 1e-5 and abs(speed_lat_swe) > 1e-5:
            assert (speed_lat_lib >= 0) == (speed_lat_swe >= 0), (
                f"speed_lat sign mismatch: lib={speed_lat_lib:.8f}, swe={speed_lat_swe:.8f}"
            )

    def test_spica_velocity_vs_pyswisseph(self, standard_jd):
        """Compare Spica velocity with pyswisseph within 10%."""
        # Get libephemeris velocity
        pos_lib, name, retflag = fixstar_ut(
            "Spica", standard_jd, FLG_SWIEPH | FLG_SPEED
        )

        # Check if pyswisseph has fixstar_ut
        if not hasattr(swe, "fixstar_ut"):
            pytest.skip("pyswisseph does not have fixstar_ut")

        # Get pyswisseph velocity
        try:
            pos_swe, name_swe, err_swe = swe.fixstar_ut(
                "Spica", standard_jd, swe.FLG_SWIEPH | swe.FLG_SPEED
            )
        except Exception as e:
            pytest.skip(f"pyswisseph fixstar_ut failed: {e}")

        # Compare velocities - should match within 10%
        speed_lon_lib = pos_lib[3]
        speed_lon_swe = pos_swe[3]

        if abs(speed_lon_swe) > 1e-10:
            pct_diff = abs((speed_lon_lib - speed_lon_swe) / speed_lon_swe) * 100
            assert pct_diff < 10, (
                f"speed_lon diff {pct_diff:.1f}%: lib={speed_lon_lib:.8f}, swe={speed_lon_swe:.8f}"
            )

    def test_precession_rate_magnitude(self, standard_jd):
        """Verify velocity magnitude matches expected precession rate."""
        pos, name, retflag = fixstar_ut(
            "Regulus", standard_jd, FLG_SWIEPH | FLG_SPEED
        )

        speed_lon = pos[3]

        # Expected precession rate: ~50.3 arcsec/year
        # = 50.3 / 3600 degrees/year
        # = 50.3 / 3600 / 365.25 degrees/day
        # = 0.0000378 deg/day
        # With proper motion and nutation, expect range of 0.00003 to 0.0002 deg/day
        expected_min = 0.00003
        expected_max = 0.0002

        assert expected_min < abs(speed_lon) < expected_max, (
            f"speed_lon {speed_lon:.8f} outside expected range [{expected_min}, {expected_max}]"
        )


@pytest.mark.unit
class TestBackwardCompatibility:
    """Tests ensuring backward compatibility."""

    def test_no_velocity_when_flag_not_set(self, standard_jd):
        """Ensure velocities are zero when FLG_SPEED not set."""
        # Test all four fixstar functions
        pos1, _, _ = fixstar_ut("Regulus", standard_jd, 0)
        pos2, _, _ = fixstar("Regulus", standard_jd, 0)
        pos3, name3, _ = fixstar2_ut("Regulus", standard_jd, 0)
        pos4, name4, _ = fixstar2("Regulus", standard_jd, 0)

        for i, (pos, func_name) in enumerate(
            [
                (pos1, "fixstar_ut"),
                (pos2, "fixstar"),
                (pos3, "fixstar2_ut"),
                (pos4, "fixstar2"),
            ]
        ):
            assert pos[3] == 0.0, f"{func_name} speed_lon should be 0"
            assert pos[4] == 0.0, f"{func_name} speed_lat should be 0"
            assert pos[5] == 0.0, f"{func_name} speed_dist should be 0"

    def test_position_unchanged_with_speed_flag(self, standard_jd):
        """Verify position is the same with or without FLG_SPEED."""
        pos_without, _, _ = fixstar_ut("Regulus", standard_jd, 0)
        pos_with, _, _ = fixstar_ut("Regulus", standard_jd, FLG_SPEED)

        # Position components should be identical
        assert abs(pos_without[0] - pos_with[0]) < 1e-10, "Longitude should be same"
        assert abs(pos_without[1] - pos_with[1]) < 1e-10, "Latitude should be same"
        assert abs(pos_without[2] - pos_with[2]) < 1e-10, "Distance should be same"
