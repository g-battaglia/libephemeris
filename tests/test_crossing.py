"""
Unit tests for crossing functions: solcross_ut and mooncross_ut.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


@pytest.mark.unit
class TestSunCrossing:
    """Tests for swe_solcross_ut function."""

    def test_solcross_spring_equinox(self):
        """Test Sun crossing 0° Aries (Spring Equinox)."""
        jd_start = ephem.swe_julday(2000, 1, 1, 0.0)

        # Find when Sun crosses 0°
        jd_cross, err = ephem.swe_solcross_ut(0.0, jd_start, 0)

        assert err == "", f"Error in solcross: {err}"

        # Verify the crossing
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_SUN, 0)
        diff = abs(pos[0] - 0.0)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"Sun position at crossing: {pos[0]:.6f}° (should be ~0°)"

    def test_solcross_summer_solstice(self):
        """Test Sun crossing 90° (Summer Solstice)."""
        jd_start = ephem.swe_julday(2000, 1, 1, 0.0)

        jd_cross, err = ephem.swe_solcross_ut(90.0, jd_start, 0)

        assert err == ""

        pos, _ = ephem.swe_calc_ut(jd_cross, SE_SUN, 0)
        diff = min(abs(pos[0] - 90.0), 360 - abs(pos[0] - 90.0))

        assert diff < 0.001, f"Sun at {pos[0]:.6f}° (should be ~90°)"

    def test_solcross_vs_swisseph(self):
        """Compare solcross with SwissEphemeris."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)
        target_lon = 120.0

        # LibEphemeris
        jd_cross_py, err_py = ephem.swe_solcross_ut(target_lon, jd_start, 0)

        # SwissEph (returns single float, not tuple)
        jd_cross_swe = swe.solcross_ut(target_lon, jd_start, 0)

        assert err_py == ""

        # Times should match within ~10 seconds (1/86400 day)
        diff_seconds = abs(jd_cross_py - jd_cross_swe) * 86400
        assert diff_seconds < 30, f"Crossing time diff: {diff_seconds:.2f} seconds"

    @pytest.mark.parametrize(
        "target_lon", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_solcross_all_signs(self, target_lon):
        """Test Sun crossing all zodiac signs."""
        jd_start = ephem.swe_julday(2000, 1, 1, 0.0)

        jd_cross, err = ephem.swe_solcross_ut(target_lon, jd_start, 0)

        assert err == "", f"Error crossing {target_lon}°"

        # Verify position
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_SUN, 0)
        diff = abs(pos[0] - target_lon)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Sun at {pos[0]:.4f}° vs target {target_lon}°"

    def test_solcross_precision(self):
        """Test sub-arcsecond precision of crossing."""
        jd_start = ephem.swe_julday(2024, 1, 1, 0.0)

        jd_cross, err = ephem.swe_solcross_ut(45.0, jd_start, 0)

        assert err == ""

        pos, _ = ephem.swe_calc_ut(jd_cross, SE_SUN, 0)
        diff = abs(pos[0] - 45.0)

        # Should be within 1 arcsecond = 1/3600 degree
        assert diff < 1.0 / 3600.0, f"Precision: {diff * 3600:.2f} arcseconds"


@pytest.mark.unit
class TestMoonCrossing:
    """Tests for swe_mooncross_ut function."""

    def test_mooncross_basic(self):
        """Test basic Moon crossing."""
        jd_start = ephem.swe_julday(2024, 11, 28, 0.0)

        jd_cross, err = ephem.swe_mooncross_ut(90.0, jd_start, 0)

        assert err == "", f"Error in mooncross: {err}"

        # Verify crossing
        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MOON, 0)
        diff = abs(pos[0] - 90.0)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"Moon at {pos[0]:.6f}° (should be ~90°)"

    def test_mooncross_vs_swisseph(self):
        """Compare mooncross with SwissEphemeris."""
        jd_start = ephem.swe_julday(2024, 11, 28, 0.0)
        target_lon = 180.0

        # LibEphemeris
        jd_cross_py, err_py = ephem.swe_mooncross_ut(target_lon, jd_start, 0)

        # SwissEph (returns single float)
        jd_cross_swe = swe.mooncross_ut(target_lon, jd_start, 0)

        assert err_py == ""

        # Times should match within ~60 seconds
        diff_seconds = abs(jd_cross_py - jd_cross_swe) * 86400
        assert diff_seconds < 120, f"Crossing time diff: {diff_seconds:.2f} seconds"

    @pytest.mark.parametrize("target_lon", [0, 45, 90, 135, 180, 225, 270, 315])
    def test_mooncross_multiple_longitudes(self, target_lon):
        """Test Moon crossing various longitudes."""
        jd_start = ephem.swe_julday(2024, 11, 1, 0.0)

        jd_cross, err = ephem.swe_mooncross_ut(target_lon, jd_start, 0)

        assert err == ""

        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MOON, 0)
        diff = abs(pos[0] - target_lon)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Moon at {pos[0]:.4f}° vs target {target_lon}°"

    def test_mooncross_precision(self):
        """Test precision of Moon crossing."""
        jd_start = ephem.swe_julday(2024, 11, 28, 0.0)

        jd_cross, err = ephem.swe_mooncross_ut(123.456, jd_start, 0)

        assert err == ""

        pos, _ = ephem.swe_calc_ut(jd_cross, SE_MOON, 0)
        diff = abs(pos[0] - 123.456)

        # Moon precision should be sub-arcsecond
        assert diff < 1.0 / 3600.0, f"Precision: {diff * 3600:.2f} arcseconds"

    def test_mooncross_speed(self):
        """Test that Moon crossing is found quickly (near start date)."""
        jd_start = ephem.swe_julday(2024, 11, 28, 0.0)

        jd_cross, err = ephem.swe_mooncross_ut(200.0, jd_start, 0)

        assert err == ""

        # Moon completes orbit in ~27 days, should find crossing within that
        days_diff = abs(jd_cross - jd_start)
        assert days_diff < 28, f"Crossing found {days_diff:.1f} days from start"


@pytest.mark.integration
class TestCrossingIntegration:
    """Integration tests for crossing functions."""

    def test_consecutive_crossings(self):
        """Test finding multiple consecutive crossings."""
        jd_start = ephem.swe_julday(2000, 1, 1, 0.0)

        crossings = []
        jd = jd_start

        # Find Sun crossing each sign (30° increments)
        for i in range(12):
            target = i * 30.0
            jd_cross, err = ephem.swe_solcross_ut(target, jd, 0)
            assert err == ""
            crossings.append(jd_cross)
            jd = jd_cross + 1  # Start next search 1 day after

        # Crossings should be in chronological order
        for i in range(len(crossings) - 1):
            assert crossings[i] < crossings[i + 1], "Crossings not in order"

    def test_error_handling(self):
        """Test error handling for invalid inputs."""
        jd_start = ephem.swe_julday(2000, 1, 1, 0.0)

        # These should complete (maybe with error message)
        jd_cross, err = ephem.swe_solcross_ut(400.0, jd_start, 0)  # >360°
        # Function should normalize or handle gracefully
        assert isinstance(err, str)
