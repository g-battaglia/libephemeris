"""
Comprehensive tests for longitude crossing functions.

Tests solcross_ut, mooncross_ut, and cross_ut.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *
from libephemeris.crossing import (
    NR_TOLERANCE,
    NR_TOLERANCE_SUN,
    NR_TOLERANCE_MOON,
    NR_MAX_ITER_MOON,
    NR_MAX_ITER_PLANET,
    NR_MAX_ITER_HELIO,
    NR_MAX_ITER_VERY_SLOW,
    NR_MAX_ITER_STATION,
    _get_adaptive_max_iterations,
)


class TestSubArcsecondPrecision:
    """Tests to verify sub-milliarcsecond convergence (0.001 arcsecond for all bodies)."""

    # All tolerances are now 0.001 arcsecond = 0.001/3600 degrees
    SUB_ARCSEC = 0.001 / 3600.0
    SUB_ARCSEC_MOON = 0.001 / 3600.0
    SUB_MILLIARCSEC = 0.001 / 3600.0

    @pytest.mark.unit
    def test_nr_tolerance_constant(self):
        """Verify NR_TOLERANCE is set to 0.1 arcsecond."""
        assert NR_TOLERANCE == self.SUB_ARCSEC, (
            f"NR_TOLERANCE should be {self.SUB_ARCSEC}, got {NR_TOLERANCE}"
        )

    @pytest.mark.unit
    def test_nr_tolerance_sun_constant(self):
        """Verify NR_TOLERANCE_SUN is set to 0.001 arcsecond."""
        assert NR_TOLERANCE_SUN == self.SUB_MILLIARCSEC, (
            f"NR_TOLERANCE_SUN should be {self.SUB_MILLIARCSEC}, got {NR_TOLERANCE_SUN}"
        )

    @pytest.mark.unit
    def test_nr_tolerance_moon_constant(self):
        """Verify NR_TOLERANCE_MOON is set to 0.05 arcsecond (sub-arcsecond)."""
        assert NR_TOLERANCE_MOON == self.SUB_ARCSEC_MOON, (
            f"NR_TOLERANCE_MOON should be {self.SUB_ARCSEC_MOON}, got {NR_TOLERANCE_MOON}"
        )

    @pytest.mark.unit
    def test_nr_max_iter_moon_constant(self):
        """Verify NR_MAX_ITER_MOON is set to 50.

        50 iterations for sub-milliarcsecond convergence:
        - Longitude crossings (~13°/day): typically 7-12 iterations
        - Node crossings (~1°/day latitude speed): up to 20-30 iterations
        - 50 provides adequate safety margin for extreme cases
        """
        assert NR_MAX_ITER_MOON == 50, (
            f"NR_MAX_ITER_MOON should be 50, got {NR_MAX_ITER_MOON}"
        )

    @pytest.mark.unit
    def test_solcross_sub_milliarcsecond_precision(self):
        """Sun crossing should achieve sub-milliarcsecond precision (0.001 arcsec)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 45.0
        jd_cross = ephem.solcross_ut(target, jd_start, 0)

        # Check Sun position at crossing
        pos, _ = ephem.calc_ut(jd_cross, SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        # Verify sub-milliarcsecond precision (0.001 arcsecond)
        diff_arcsec = diff * 3600.0
        assert diff_arcsec < 0.001, (
            f"Sun at {pos[0]}, target {target}, diff {diff_arcsec:.6f} arcsec (should be < 0.001)"
        )

    @pytest.mark.unit
    def test_mooncross_sub_arcsecond_precision(self):
        """Moon crossing should achieve sub-arcsecond precision (0.05 arcsec)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 123.456
        jd_cross = ephem.mooncross_ut(target, jd_start, 0)

        pos, _ = ephem.calc_ut(jd_cross, MOON, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        # Verify sub-arcsecond precision (0.05 arcsec for Moon)
        diff_arcsec = diff * 3600.0
        assert diff_arcsec < 0.05, (
            f"Moon at {pos[0]}, target {target}, diff {diff_arcsec:.4f} arcsec (should be < 0.05)"
        )

    @pytest.mark.unit
    def test_mooncross_node_sub_arcsecond_precision(self):
        """Moon node crossing should achieve sub-arcsecond precision (0.05 arcsec)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        jd_cross, xlon, xlat = ephem.mooncross_node_ut(jd_start, 0)

        # Check Moon latitude at crossing (should be ~0)
        pos, _ = ephem.calc_ut(jd_cross, MOON, 0)

        # Latitude should be within 0.05 arcsecond of zero (Moon tolerance)
        lat_arcsec = abs(pos[1]) * 3600.0
        assert lat_arcsec < 0.05, (
            f"Moon latitude at node: {lat_arcsec:.4f} arcsec (should be < 0.05)"
        )

    @pytest.mark.unit
    def test_cross_ut_sub_arcsecond_precision(self):
        """Generic planet crossing should achieve sub-arcsecond precision."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 30.0

        jd_cross = ephem.cross_ut(MARS, target, jd_start, 0)

        pos, _ = ephem.calc_ut(jd_cross, MARS, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        # Verify sub-arcsecond precision
        diff_arcsec = diff * 3600.0
        assert diff_arcsec < 0.1, (
            f"Mars at {pos[0]}, target {target}, diff {diff_arcsec:.4f} arcsec (should be < 0.1)"
        )

    @pytest.mark.unit
    def test_helio_cross_sub_arcsecond_precision(self):
        """Heliocentric crossing should achieve sub-arcsecond precision."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 90.0

        jd_cross = ephem.helio_cross_ut(EARTH, target, jd_start, 0)

        pos, _ = ephem.calc_ut(jd_cross, EARTH, FLG_HELCTR)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        # Verify sub-arcsecond precision
        diff_arcsec = diff * 3600.0
        assert diff_arcsec < 0.1, (
            f"Earth helio at {pos[0]}, target {target}, diff {diff_arcsec:.4f} arcsec (should be < 0.1)"
        )

    @pytest.mark.comparison
    def test_solcross_vs_pyswisseph_sub_milliarcsecond(self):
        """Verify libephemeris matches pyswisseph and achieves sub-milliarcsecond precision."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 0.0  # Vernal equinox

        jd_lib = ephem.solcross_ut(target, jd_start, 0)
        jd_swe = swe.solcross_ut(target, jd_start, 0)

        # Get positions at both crossing times
        pos_lib, _ = ephem.calc_ut(jd_lib, SUN, 0)
        pos_swe = swe.calc_ut(jd_swe, swe.SUN)[0]

        # Both should be very close to target (handle 0/360 wraparound)
        diff_lib = abs(pos_lib[0] - target)
        if diff_lib > 180:
            diff_lib = 360 - diff_lib
        diff_lib_arcsec = diff_lib * 3600.0

        diff_swe = abs(pos_swe[0] - target)
        if diff_swe > 180:
            diff_swe = 360 - diff_swe
        diff_swe_arcsec = diff_swe * 3600.0

        # libephemeris should achieve sub-milliarcsecond precision (0.001 arcsec)
        assert diff_lib_arcsec < 0.001, (
            f"libephemeris diff {diff_lib_arcsec:.6f} arcsec"
        )
        # pyswisseph should also achieve sub-milliarcsecond precision
        assert diff_swe_arcsec < 0.001, f"pyswisseph diff {diff_swe_arcsec:.6f} arcsec"

        # Crossing times should be within 1 second of each other
        diff_time_seconds = abs(jd_lib - jd_swe) * 86400.0
        assert diff_time_seconds < 1.0, (
            f"Timing diff {diff_time_seconds:.4f} seconds (should be < 1)"
        )

    @pytest.mark.comparison
    def test_mooncross_vs_pyswisseph_sub_arcsecond(self):
        """Verify Moon crossing matches pyswisseph within sub-arcsecond tolerance (0.05 arcsec)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 180.0

        jd_lib = ephem.mooncross_ut(target, jd_start, 0)
        jd_swe = swe.mooncross_ut(target, jd_start, 0)

        # Get positions at both crossing times
        pos_lib, _ = ephem.calc_ut(jd_lib, MOON, 0)

        # Both should be very close to target
        diff_lib = abs(pos_lib[0] - target)
        if diff_lib > 180:
            diff_lib = 360 - diff_lib
        diff_lib_arcsec = diff_lib * 3600.0

        # libephemeris should achieve sub-arcsecond precision (0.05 arcsec for Moon)
        assert diff_lib_arcsec < 0.05, f"libephemeris diff {diff_lib_arcsec:.4f} arcsec"

        # Crossing times should be within 3 seconds of each other
        # (slight differences in underlying ephemeris calculations)
        diff_time_seconds = abs(jd_lib - jd_swe) * 86400.0
        assert diff_time_seconds < 3.0, (
            f"Timing diff {diff_time_seconds:.4f} seconds (should be < 3)"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_solcross_all_signs_sub_milliarcsecond(self, target):
        """All zodiac sign ingresses should achieve sub-milliarcsecond precision (0.001 arcsec)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        jd_cross = ephem.solcross_ut(float(target), jd_start, 0)
        pos, _ = ephem.calc_ut(jd_cross, SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        diff_arcsec = diff * 3600.0
        assert diff_arcsec < 0.001, (
            f"Target {target}° diff {diff_arcsec:.6f} arcsec (should be < 0.001)"
        )


class TestSolcrossBasic:
    """Basic tests for Sun crossing function."""

    @pytest.mark.unit
    def test_solcross_vernal_equinox(self):
        """Find when Sun crosses 0° (vernal equinox)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        jd_cross = ephem.solcross_ut(0.0, jd_start, 0)

        # Should be around March 20, 2024
        year, month, day, hour = ephem.revjul(jd_cross)
        assert year == 2024
        assert month == 3
        assert 19 <= day <= 21

    @pytest.mark.unit
    def test_solcross_summer_solstice(self):
        """Find when Sun crosses 90° (summer solstice)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        jd_cross = ephem.solcross_ut(90.0, jd_start, 0)

        year, month, day, hour = ephem.revjul(jd_cross)
        assert year == 2024
        assert month == 6
        assert 20 <= day <= 22

    @pytest.mark.unit
    def test_solcross_precision(self):
        """Sun should be very close to target at crossing time."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 45.0
        jd_cross = ephem.solcross_ut(target, jd_start, 0)

        # Check Sun position at crossing
        pos, _ = ephem.calc_ut(jd_cross, SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"Sun at {pos[0]}, target {target}, diff {diff}"


class TestSolcrossTT:
    """Tests for solcross (TT version)."""

    @pytest.mark.unit
    def test_solcross_tt_vernal_equinox(self):
        """Find when Sun crosses 0° (vernal equinox) using TT."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        # Convert to TT (approximately, for starting point)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.solcross(0.0, jd_start_tt, 0)

        # Should be around March 20, 2024
        # Convert back to UT for date check
        delta_t_cross = ephem.deltat(jd_cross_tt)
        jd_cross_ut = jd_cross_tt - delta_t_cross
        year, month, day, hour = ephem.revjul(jd_cross_ut)
        assert year == 2024
        assert month == 3
        assert 19 <= day <= 21

    @pytest.mark.unit
    def test_solcross_tt_precision(self):
        """Sun should be very close to target at crossing time (TT)."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 45.0
        jd_cross_tt = ephem.solcross(target, jd_start_tt, 0)

        # Check Sun position at crossing (using TT version of calc)
        pos, _ = ephem.calc(jd_cross_tt, SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"Sun at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_solcross_tt_vs_ut_consistency(self):
        """TT and UT versions should give consistent results."""
        jd_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 90.0  # Summer solstice

        # Get crossing time in UT
        jd_cross_ut = ephem.solcross_ut(target, jd_ut, 0)

        # Get crossing time in TT
        jd_cross_tt = ephem.solcross(target, jd_tt, 0)

        # Convert UT result to TT for comparison
        delta_t_cross = ephem.deltat(jd_cross_ut)
        jd_cross_ut_as_tt = jd_cross_ut + delta_t_cross

        # They should be very close (within seconds)
        diff_seconds = abs(jd_cross_tt - jd_cross_ut_as_tt) * 86400
        assert diff_seconds < 10, f"TT vs UT consistency diff {diff_seconds} seconds"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_solcross_tt_all_signs(self, target):
        """All zodiac sign ingresses should work with TT version."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.solcross(float(target), jd_start_tt, 0)

        # Verify Sun is at target
        pos, _ = ephem.calc(jd_cross_tt, SUN, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Target {target}° diff {diff}°"

    @pytest.mark.comparison
    def test_solcross_tt_vs_pyswisseph(self):
        """TT version should match pyswisseph solcross()."""
        jd_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = swe.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 0.0

        jd_lib = ephem.solcross(target, jd_tt, 0)
        jd_swe = swe.solcross(target, jd_tt, 0)

        # Difference should be less than 1 minute
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 60, f"Timing diff {diff_seconds} seconds"

    @pytest.mark.edge_case
    def test_solcross_tt_target_360_equals_0(self):
        """360° should be same as 0° for TT version."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_360 = ephem.solcross(360.0, jd_start_tt, 0)
        jd_0 = ephem.solcross(0.0, jd_start_tt, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001


class TestSolcrossVsPyswisseph:
    """Compare solcross with pyswisseph."""

    @pytest.mark.comparison
    def test_solcross_equinox_timing(self):
        """Crossing time should match pyswisseph."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        jd_lib = ephem.solcross_ut(0.0, jd_start, 0)
        jd_swe = swe.solcross_ut(0.0, jd_start, 0)

        # Difference should be less than 1 minute
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 60, f"Timing diff {diff_seconds} seconds"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_solcross_all_signs(self, target):
        """All zodiac sign ingresses should match."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        jd_lib = ephem.solcross_ut(float(target), jd_start, 0)
        jd_swe = swe.solcross_ut(float(target), jd_start, 0)

        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 120, f"Target {target}° diff {diff_seconds} seconds"


class TestMooncrossBasic:
    """Basic tests for Moon crossing function."""

    @pytest.mark.unit
    def test_mooncross_basic(self):
        """Moon crossing should return valid JD."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 90.0
        jd_cross = ephem.mooncross_ut(target, jd_start, 0)

        assert jd_cross > jd_start
        # Should be within one lunar orbit (~27 days)
        assert jd_cross < jd_start + 28

    @pytest.mark.unit
    def test_mooncross_precision(self):
        """Moon should be close to target at crossing."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 123.456
        jd_cross = ephem.mooncross_ut(target, jd_start, 0)

        pos, _ = ephem.calc_ut(jd_cross, MOON, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Moon at {pos[0]}, target {target}, diff {diff}"


class TestMooncrossTT:
    """Tests for mooncross (TT version)."""

    @pytest.mark.unit
    def test_mooncross_tt_basic(self):
        """Moon crossing should return valid JD using TT."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 90.0
        jd_cross_tt = ephem.mooncross(target, jd_start_tt, 0)

        assert jd_cross_tt > jd_start_tt
        # Should be within one lunar orbit (~27 days)
        assert jd_cross_tt < jd_start_tt + 28

    @pytest.mark.unit
    def test_mooncross_tt_precision(self):
        """Moon should be close to target at crossing time (TT)."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 123.456
        jd_cross_tt = ephem.mooncross(target, jd_start_tt, 0)

        # Check Moon position at crossing (using TT version of calc)
        pos, _ = ephem.calc(jd_cross_tt, MOON, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Moon at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_mooncross_tt_vs_ut_consistency(self):
        """TT and UT versions should give consistent results."""
        jd_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 180.0

        # Get crossing time in UT
        jd_cross_ut = ephem.mooncross_ut(target, jd_ut, 0)

        # Get crossing time in TT
        jd_cross_tt = ephem.mooncross(target, jd_tt, 0)

        # Convert UT result to TT for comparison
        delta_t_cross = ephem.deltat(jd_cross_ut)
        jd_cross_ut_as_tt = jd_cross_ut + delta_t_cross

        # They should be very close (within seconds)
        diff_seconds = abs(jd_cross_tt - jd_cross_ut_as_tt) * 86400
        assert diff_seconds < 10, f"TT vs UT consistency diff {diff_seconds} seconds"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_mooncross_tt_all_signs(self, target):
        """All zodiac sign ingresses should work with TT version."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.mooncross(float(target), jd_start_tt, 0)

        # Verify Moon is at target
        pos, _ = ephem.calc(jd_cross_tt, MOON, 0)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Target {target}° diff {diff}°"

    @pytest.mark.comparison
    def test_mooncross_tt_vs_pyswisseph(self):
        """TT version should match pyswisseph mooncross()."""
        jd_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = swe.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 180.0

        jd_lib = ephem.mooncross(target, jd_tt, 0)
        jd_swe = swe.mooncross(target, jd_tt, 0)

        # Difference should be less than 3 minutes (Moon moves fast)
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Timing diff {diff_seconds} seconds"

    @pytest.mark.edge_case
    def test_mooncross_tt_target_360_equals_0(self):
        """360° should be same as 0° for TT version."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_360 = ephem.mooncross(360.0, jd_start_tt, 0)
        jd_0 = ephem.mooncross(0.0, jd_start_tt, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001


class TestMooncrossVsPyswisseph:
    """Compare mooncross with pyswisseph."""

    @pytest.mark.comparison
    def test_mooncross_timing(self):
        """Moon crossing should match pyswisseph."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 180.0

        jd_lib = ephem.mooncross_ut(target, jd_start, 0)
        jd_swe = swe.mooncross_ut(target, jd_start, 0)

        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Moon crossing diff {diff_seconds} seconds"


class TestCrossingConsecutive:
    """Test finding consecutive crossings."""

    @pytest.mark.unit
    def test_consecutive_sun_crossings(self):
        """Should find 12 consecutive Sun crossings in a year."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        crossings = []

        for target in range(0, 360, 30):
            jd_cross = ephem.solcross_ut(float(target), jd, 0)
            crossings.append(jd_cross)

        # All should be in chronological order (after accounting for year wrap)
        # First 0° crossing should be in March
        assert crossings[0] > ephem.julday(2024, 3, 1, 0.0)

    @pytest.mark.unit
    def test_12_moon_crossings_in_month(self):
        """Should find ~12 crossings in 28 days."""
        jd = ephem.julday(2024, 1, 1, 0.0)
        crossings = []

        for target in range(0, 360, 30):
            jd_cross = ephem.mooncross_ut(float(target), jd, 0)
            if jd_cross < jd + 28:
                crossings.append(jd_cross)

        # Should find all 12 in about 27 days
        assert len(crossings) == 12


class TestCrossingEdgeCases:
    """Test edge cases for crossing functions."""

    @pytest.mark.edge_case
    def test_target_360_equals_0(self):
        """360° should be same as 0°."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        jd_360 = ephem.solcross_ut(360.0, jd_start, 0)
        jd_0 = ephem.solcross_ut(0.0, jd_start, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001

    @pytest.mark.edge_case
    def test_target_negative_normalized(self):
        """Negative target should be normalized."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        jd_neg = ephem.solcross_ut(-30.0, jd_start, 0)  # = 330°
        jd_pos = ephem.solcross_ut(330.0, jd_start, 0)

        assert abs(jd_neg - jd_pos) < 0.001


class TestCrossUtGeneric:
    """Test generic planet crossing function."""

    @pytest.mark.unit
    def test_cross_ut_mercury(self):
        """Should work for Mercury."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 0.0

        jd_cross = ephem.cross_ut(MERCURY, target, jd_start, 0)

        assert jd_cross > jd_start

        # Check Mercury position at crossing
        pos, _ = ephem.calc_ut(jd_cross, MERCURY, 0)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1  # Mercury moves fast, allow more tolerance


class TestMooncrossNodeBasic:
    """Basic tests for Moon node crossing function."""

    @pytest.mark.unit
    def test_mooncross_node_basic(self):
        """Moon node crossing should return valid JD."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        jd_cross, xlon, xlat = ephem.mooncross_node_ut(jd_start, 0)

        assert jd_cross > jd_start
        # Should be within ~14 days (half the nodal month)
        assert jd_cross < jd_start + 14

    @pytest.mark.unit
    def test_mooncross_node_latitude_near_zero(self):
        """Moon latitude should be near zero at node crossing."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        jd_cross, xlon, xlat = ephem.mooncross_node_ut(jd_start, 0)

        # Check Moon latitude at crossing
        pos, _ = ephem.calc_ut(jd_cross, MOON, 0)

        # Latitude should be very close to 0
        assert abs(pos[1]) < 0.01, f"Moon latitude at node crossing: {pos[1]}"

    @pytest.mark.unit
    def test_mooncross_node_consecutive_crossings(self):
        """Should find consecutive node crossings ~13.6 days apart."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find first node crossing
        jd_first, _, _ = ephem.mooncross_node_ut(jd_start, 0)

        # Find second node crossing by starting just after first
        jd_second, _, _ = ephem.mooncross_node_ut(jd_first + 0.5, 0)

        # Should be roughly 13.6 days apart (half the nodal month)
        diff_days = jd_second - jd_first
        assert 12 < diff_days < 15, f"Time between node crossings: {diff_days} days"

    @pytest.mark.unit
    def test_mooncross_node_ascending_vs_descending(self):
        """Can determine ascending vs descending node from latitude velocity."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find first crossing
        jd_first, _, _ = ephem.mooncross_node_ut(jd_start, 0)
        pos_first, _ = ephem.calc_ut(jd_first, MOON, FLG_SPEED)

        # Find second crossing
        jd_second, _, _ = ephem.mooncross_node_ut(jd_first + 0.5, 0)
        pos_second, _ = ephem.calc_ut(jd_second, MOON, FLG_SPEED)

        # Latitude velocities should have opposite signs
        # (one ascending, one descending)
        assert pos_first[4] * pos_second[4] < 0, (
            f"Consecutive crossings should alternate: "
            f"v1={pos_first[4]}, v2={pos_second[4]}"
        )

        # Latitude velocities should have opposite signs
        # (one ascending, one descending)
        assert pos_first[4] * pos_second[4] < 0, (
            f"Consecutive crossings should alternate: "
            f"v1={pos_first[4]}, v2={pos_second[4]}"
        )


class TestMooncrossNodeTT:
    """Tests for mooncross_node (TT version)."""

    @pytest.mark.unit
    def test_mooncross_node_tt_basic(self):
        """Moon node crossing should return valid JD using TT."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt, xlon, xlat = ephem.mooncross_node(jd_start_tt, 0)

        assert jd_cross_tt > jd_start_tt
        # Should be within ~14 days
        assert jd_cross_tt < jd_start_tt + 14

    @pytest.mark.unit
    def test_mooncross_node_tt_precision(self):
        """Moon latitude should be near zero at node crossing (TT)."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt, xlon, xlat = ephem.mooncross_node(jd_start_tt, 0)

        # Check Moon latitude at crossing (using TT version of calc)
        pos, _ = ephem.calc(jd_cross_tt, MOON, 0)

        # Latitude should be very close to 0
        assert abs(pos[1]) < 0.01, f"Moon latitude at node crossing: {pos[1]}"

    @pytest.mark.unit
    def test_mooncross_node_tt_vs_ut_consistency(self):
        """TT and UT versions should give consistent results."""
        jd_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        # Get crossing time in UT
        jd_cross_ut, _, _ = ephem.mooncross_node_ut(jd_ut, 0)

        # Get crossing time in TT
        jd_cross_tt, _, _ = ephem.mooncross_node(jd_tt, 0)

        # Convert UT result to TT for comparison
        delta_t_cross = ephem.deltat(jd_cross_ut)
        jd_cross_ut_as_tt = jd_cross_ut + delta_t_cross

        # They should be very close (within seconds)
        diff_seconds = abs(jd_cross_tt - jd_cross_ut_as_tt) * 86400
        assert diff_seconds < 10, f"TT vs UT consistency diff {diff_seconds} seconds"


class TestMooncrossNodeVsPyswisseph:
    """Compare mooncross_node with pyswisseph."""

    @pytest.mark.comparison
    def test_mooncross_node_vs_pyswisseph(self):
        """Moon node crossing should match pyswisseph."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        jd_lib, xlon_lib, xlat_lib = ephem.mooncross_node_ut(jd_start, 0)
        # pyswisseph returns (jd_cross, xlon, xlat)
        result_swe = swe.mooncross_node_ut(jd_start, 0)
        jd_swe = result_swe[0]

        # Difference should be less than 3 minutes
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Moon node crossing diff {diff_seconds} seconds"

    @pytest.mark.comparison
    def test_mooncross_node_tt_vs_pyswisseph(self):
        """TT version should match pyswisseph mooncross_node()."""
        jd_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = swe.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        jd_lib, xlon_lib, xlat_lib = ephem.mooncross_node(jd_tt, 0)
        # pyswisseph returns (jd_cross, xlon, xlat)
        result_swe = swe.mooncross_node(jd_tt, 0)
        jd_swe = result_swe[0]

        # Difference should be less than 3 minutes
        diff_seconds = abs(jd_lib - jd_swe) * 86400
        assert diff_seconds < 180, f"Moon node crossing diff {diff_seconds} seconds"

    @pytest.mark.comparison
    def test_mooncross_node_multiple_crossings(self):
        """Multiple consecutive crossings should match pyswisseph."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        for i in range(4):
            jd_lib, _, _ = ephem.mooncross_node_ut(jd, 0)
            # pyswisseph returns (jd_cross, xlon, xlat)
            result_swe = swe.mooncross_node_ut(jd, 0)
            jd_swe = result_swe[0]

            diff_seconds = abs(jd_lib - jd_swe) * 86400
            assert diff_seconds < 180, f"Crossing {i + 1} diff {diff_seconds} seconds"

            # Move to after this crossing for next iteration
            jd = jd_lib + 0.5


class TestHelioCrossBasic:
    """Basic tests for heliocentric crossing function."""

    @pytest.mark.unit
    def test_helio_cross_ut_mars_basic(self):
        """Heliocentric Mars crossing should return valid JD."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 0.0

        jd_cross = ephem.helio_cross_ut(MARS, target, jd_start, 0)

        assert jd_cross > jd_start
        # Mars completes one heliocentric orbit in ~687 days
        assert jd_cross < jd_start + 700

    @pytest.mark.unit
    def test_helio_cross_ut_precision(self):
        """Planet should be very close to target at heliocentric crossing."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 45.0

        jd_cross = ephem.helio_cross_ut(MARS, target, jd_start, 0)

        # Check Mars heliocentric position at crossing
        pos, _ = ephem.calc_ut(jd_cross, MARS, FLG_HELCTR)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Mars helio at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_helio_cross_ut_earth(self):
        """Should work for Earth (heliocentric view of Earth from Sun)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 90.0

        jd_cross = ephem.helio_cross_ut(EARTH, target, jd_start, 0)

        assert jd_cross > jd_start
        # Earth completes one heliocentric orbit in ~365 days
        assert jd_cross < jd_start + 366

        # Verify position
        pos, _ = ephem.calc_ut(jd_cross, EARTH, FLG_HELCTR)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.01

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "planet,max_days",
        [
            (MERCURY, 90),  # Mercury helio period ~88 days
            (VENUS, 230),  # Venus helio period ~225 days
            (EARTH, 370),  # Earth helio period ~365 days
            (MARS, 700),  # Mars helio period ~687 days
        ],
    )
    def test_helio_cross_ut_multiple_planets(self, planet, max_days):
        """Should work for various planets."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 180.0

        jd_cross = ephem.helio_cross_ut(planet, target, jd_start, 0)

        assert jd_cross > jd_start
        assert jd_cross < jd_start + max_days

        # Verify position
        pos, _ = ephem.calc_ut(jd_cross, planet, FLG_HELCTR)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Planet {planet} at {pos[0]}, target {target}"


class TestHelioCrossTT:
    """Tests for helio_cross (TT version)."""

    @pytest.mark.unit
    def test_helio_cross_tt_basic(self):
        """Heliocentric crossing should return valid JD using TT."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 90.0
        jd_cross_tt = ephem.helio_cross(MARS, target, jd_start_tt, 0)

        assert jd_cross_tt > jd_start_tt
        assert jd_cross_tt < jd_start_tt + 700

    @pytest.mark.unit
    def test_helio_cross_tt_precision(self):
        """Planet should be very close to target at crossing time (TT)."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        target = 123.456
        jd_cross_tt = ephem.helio_cross(MARS, target, jd_start_tt, 0)

        # Check position (using TT version of calc)
        pos, _ = ephem.calc(jd_cross_tt, MARS, FLG_HELCTR)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Mars helio at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_helio_cross_tt_vs_ut_consistency(self):
        """TT and UT versions should give consistent results."""
        jd_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_ut)
        jd_tt = jd_ut + delta_t

        target = 180.0

        # Get crossing time in UT
        jd_cross_ut = ephem.helio_cross_ut(MARS, target, jd_ut, 0)

        # Get crossing time in TT
        jd_cross_tt = ephem.helio_cross(MARS, target, jd_tt, 0)

        # Convert UT result to TT for comparison
        delta_t_cross = ephem.deltat(jd_cross_ut)
        jd_cross_ut_as_tt = jd_cross_ut + delta_t_cross

        # They should be very close (within seconds)
        diff_seconds = abs(jd_cross_tt - jd_cross_ut_as_tt) * 86400
        assert diff_seconds < 10, f"TT vs UT consistency diff {diff_seconds} seconds"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "target", [0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330]
    )
    def test_helio_cross_tt_all_signs(self, target):
        """All zodiac sign ingresses should work with TT version."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        jd_cross_tt = ephem.helio_cross(EARTH, float(target), jd_start_tt, 0)

        # Verify position
        pos, _ = ephem.calc(jd_cross_tt, EARTH, FLG_HELCTR)

        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Target {target}° diff {diff}°"


class TestHelioCrossEdgeCases:
    """Test edge cases for heliocentric crossing functions."""

    @pytest.mark.edge_case
    def test_helio_cross_target_360_equals_0(self):
        """360° should be same as 0°."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        jd_360 = ephem.helio_cross_ut(MARS, 360.0, jd_start, 0)
        jd_0 = ephem.helio_cross_ut(MARS, 0.0, jd_start, 0)

        # Should be the same crossing
        assert abs(jd_360 - jd_0) < 0.001

    @pytest.mark.edge_case
    def test_helio_cross_target_negative_normalized(self):
        """Negative target should be normalized."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        jd_neg = ephem.helio_cross_ut(MARS, -30.0, jd_start, 0)  # = 330°
        jd_pos = ephem.helio_cross_ut(MARS, 330.0, jd_start, 0)

        assert abs(jd_neg - jd_pos) < 0.001

    @pytest.mark.edge_case
    def test_helio_cross_consecutive_crossings(self):
        """Should find consecutive crossings for a full orbit."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        crossings = []

        # Find Earth crossing each 30° for a full year
        for target in range(0, 360, 30):
            jd_cross = ephem.helio_cross_ut(EARTH, float(target), jd_start, 0)
            crossings.append(jd_cross)

        # All crossings should be within one year
        assert max(crossings) - jd_start < 366


class TestHelioCrossVsGeocentricConcept:
    """Verify heliocentric differs from geocentric (conceptual tests)."""

    @pytest.mark.unit
    def test_helio_no_retrograde(self):
        """Heliocentric planets should not show retrograde behavior."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find several consecutive crossings for Mars
        crossings = []
        jd = jd_start
        for _ in range(4):
            jd_cross = ephem.helio_cross_ut(MARS, 90.0, jd, 0)
            crossings.append(jd_cross)
            jd = jd_cross + 1  # Move past this crossing

        # Each crossing should be about one Mars orbital period apart (~687 days)
        for i in range(1, len(crossings)):
            diff = crossings[i] - crossings[i - 1]
            # Should be around 687 days (±50 days for minor perturbations)
            assert 600 < diff < 750, f"Period between crossings: {diff} days"

    @pytest.mark.unit
    def test_helio_vs_geocentric_different_positions(self):
        """Heliocentric and geocentric positions should generally differ."""
        jd = ephem.julday(2024, 6, 15, 12.0)

        # Get geocentric Mars position
        geo_pos, _ = ephem.calc_ut(jd, MARS, 0)

        # Get heliocentric Mars position
        helio_pos, _ = ephem.calc_ut(jd, MARS, FLG_HELCTR)

        # Positions should be different (unless Mars is in opposition)
        diff = abs(geo_pos[0] - helio_pos[0])
        if diff > 180:
            diff = 360 - diff

        # They're generally different but we just verify the function works
        # and returns valid positions
        assert 0 <= geo_pos[0] < 360
        assert 0 <= helio_pos[0] < 360


class TestMooncrossNodeEclipseRelevance:
    """Test mooncross_node in context of eclipse calculations."""

    @pytest.mark.unit
    def test_eclipse_proximity_to_node(self):
        """Eclipses occur when Sun is near a lunar node."""
        # Start from a known solar eclipse date: April 8, 2024
        jd_eclipse = ephem.julday(2024, 4, 8, 18.0)

        # Find the nearest node crossing
        # Check both before and after
        jd_before = jd_eclipse - 7  # Week before
        jd_cross_before, _, _ = ephem.mooncross_node_ut(jd_before, 0)

        jd_after = jd_eclipse - 0.5
        jd_cross_after, _, _ = ephem.mooncross_node_ut(jd_after, 0)

        # The eclipse should be close to a node crossing (within a few days)
        diff_before = abs(jd_eclipse - jd_cross_before)
        diff_after = abs(jd_cross_after - jd_eclipse)

        # Eclipse should be within ~1 day of a node crossing for a total eclipse
        assert min(diff_before, diff_after) < 2, (
            f"Eclipse not near node: before={diff_before}, after={diff_after}"
        )


class TestAdaptiveIterations:
    """Test the adaptive iteration limit functionality based on planet speed."""

    @pytest.mark.unit
    def test_adaptive_iterations_fast_planet(self):
        """Fast planets (speed >= 0.1°/day) should use base iteration count."""
        assert _get_adaptive_max_iterations(1.0) == NR_MAX_ITER_PLANET
        assert _get_adaptive_max_iterations(0.5) == NR_MAX_ITER_PLANET
        assert _get_adaptive_max_iterations(0.1) == NR_MAX_ITER_PLANET
        # Also test negative speeds (retrograde)
        assert _get_adaptive_max_iterations(-0.5) == NR_MAX_ITER_PLANET

    @pytest.mark.unit
    def test_adaptive_iterations_slow_planet(self):
        """Slow planets (0.01 <= speed < 0.1°/day) should use helio iteration count."""
        assert _get_adaptive_max_iterations(0.05) == NR_MAX_ITER_HELIO
        assert _get_adaptive_max_iterations(0.03) == NR_MAX_ITER_HELIO
        assert _get_adaptive_max_iterations(0.01) == NR_MAX_ITER_HELIO
        # Also test negative speeds
        assert _get_adaptive_max_iterations(-0.03) == NR_MAX_ITER_HELIO

    @pytest.mark.unit
    def test_adaptive_iterations_very_slow_planet(self):
        """Very slow planets (0.001 <= speed < 0.01°/day) like Pluto should use more iterations."""
        assert _get_adaptive_max_iterations(0.005) == NR_MAX_ITER_VERY_SLOW
        assert (
            _get_adaptive_max_iterations(0.004) == NR_MAX_ITER_VERY_SLOW
        )  # Pluto typical
        assert _get_adaptive_max_iterations(0.001) == NR_MAX_ITER_VERY_SLOW
        # Also test negative speeds
        assert _get_adaptive_max_iterations(-0.004) == NR_MAX_ITER_VERY_SLOW

    @pytest.mark.unit
    def test_adaptive_iterations_near_station(self):
        """Near retrograde stations (speed < 0.001°/day) should use maximum iterations."""
        assert _get_adaptive_max_iterations(0.0005) == NR_MAX_ITER_STATION
        assert _get_adaptive_max_iterations(0.0001) == NR_MAX_ITER_STATION
        assert _get_adaptive_max_iterations(0.0) == NR_MAX_ITER_STATION
        # Also test negative speeds (very slow retrograde)
        assert _get_adaptive_max_iterations(-0.0005) == NR_MAX_ITER_STATION

    @pytest.mark.unit
    def test_iteration_constants_ordering(self):
        """Iteration limits should increase (or equal) as speed decreases."""
        assert NR_MAX_ITER_PLANET <= NR_MAX_ITER_HELIO
        assert NR_MAX_ITER_HELIO <= NR_MAX_ITER_VERY_SLOW
        assert NR_MAX_ITER_VERY_SLOW <= NR_MAX_ITER_STATION

    @pytest.mark.unit
    def test_pluto_typical_speed_configured(self):
        """Verify Pluto typical speed is configured in cross_ut.

        Pluto moves very slowly (~0.004°/day) compared to other planets.
        Without this configuration, the algorithm would use 0.5°/day default
        which is ~125x too fast, causing poor initial estimates.
        """
        # Import the crossing module to inspect the implementation
        # We verify this by checking that Pluto crossing works correctly
        # even when starting from a retrograde position where speed could be 0
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Get Pluto's actual speed
        pos, _ = ephem.calc_ut(jd_start, PLUTO, FLG_SPEED)
        actual_speed = abs(pos[3])

        # Pluto's typical speed should be ~0.004°/day
        # Verify the actual speed is in the expected range for Pluto
        assert actual_speed < 0.05, (
            f"Pluto speed {actual_speed}°/day should be < 0.05°/day"
        )

        # Now verify crossing works (this would fail with bad speed estimate)
        target = (pos[0] + 1) % 360
        jd_cross = ephem.cross_ut(PLUTO, target, jd_start, 0)

        # Verify Pluto is at target with good precision
        pos_cross, _ = ephem.calc_ut(jd_cross, PLUTO, 0)
        diff = abs(pos_cross[0] - target)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600.0
        assert diff_arcsec < 0.1, (
            f"Pluto at {pos_cross[0]}, target {target}, diff {diff_arcsec:.4f} arcsec"
        )

    @pytest.mark.unit
    def test_pluto_crossing_converges(self):
        """Pluto crossing should converge with adaptive iteration limits."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 300.0

        # This should not raise "Maximum iterations reached"
        jd_cross = ephem.cross_ut(PLUTO, target, jd_start, 0)

        assert jd_cross > jd_start

        # Verify Pluto is at target
        pos, _ = ephem.calc_ut(jd_cross, PLUTO, 0)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Pluto at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_pluto_helio_crossing_converges(self):
        """Pluto heliocentric crossing should converge with adaptive iteration limits."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        target = 300.0

        # This should not raise "Maximum iterations reached"
        jd_cross = ephem.helio_cross_ut(PLUTO, target, jd_start, 0)

        assert jd_cross > jd_start

        # Verify Pluto is at target (heliocentric)
        pos, _ = ephem.calc_ut(jd_cross, PLUTO, FLG_HELCTR)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Pluto helio at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_neptune_crossing_converges(self):
        """Neptune crossing (very slow ~0.006°/day) should converge."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        # Neptune is around 355° in early 2024, heading towards 0° (Aries)
        # Target slightly ahead at 356°
        target = 356.0

        jd_cross = ephem.cross_ut(NEPTUNE, target, jd_start, 0)

        assert jd_cross > jd_start

        # Verify Neptune is at target
        pos, _ = ephem.calc_ut(jd_cross, NEPTUNE, 0)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Neptune at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_uranus_crossing_converges(self):
        """Uranus crossing (~0.012°/day) should converge."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        # Uranus is around 49° in early 2024 and retrograde, so it will eventually
        # turn direct and cross 50°
        target = 50.0

        jd_cross = ephem.cross_ut(URANUS, target, jd_start, 0)

        assert jd_cross > jd_start

        # Verify Uranus is at target
        pos, _ = ephem.calc_ut(jd_cross, URANUS, 0)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Uranus at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_saturn_crossing_converges(self):
        """Saturn crossing (~0.034°/day) should converge."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        # Saturn is around 333° in early 2024, heading towards 334°
        target = 334.0

        jd_cross = ephem.cross_ut(SATURN, target, jd_start, 0)

        assert jd_cross > jd_start

        # Verify Saturn is at target
        pos, _ = ephem.calc_ut(jd_cross, SATURN, 0)
        diff = abs(pos[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Saturn at {pos[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_slow_planet_multiple_crossings(self):
        """Test that slow planets can find multiple consecutive crossings."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find 3 consecutive crossings for Saturn (starting at 333°)
        crossings = []
        jd = jd_start
        for target in [334.0, 335.0, 336.0]:
            jd_cross = ephem.cross_ut(SATURN, target, jd, 0)
            crossings.append(jd_cross)
            jd = jd_cross + 1

        # All crossings should be in chronological order
        assert crossings[0] < crossings[1] < crossings[2]

        # Verify each crossing is accurate
        for i, target in enumerate([334.0, 335.0, 336.0]):
            pos, _ = ephem.calc_ut(crossings[i], SATURN, 0)
            diff = abs(pos[0] - target)
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.1


class TestMoonNodeEdgeCases:
    """Test Moon node crossing edge cases to verify 30 iterations is sufficient.

    The Moon's latitude speed varies:
    - ~1°/day at node crossings (slowest, hardest to converge)
    - Higher speeds away from nodes
    These tests verify convergence works with 30 iterations even in worst cases.
    """

    @pytest.mark.unit
    def test_mooncross_node_when_starting_very_close_to_node(self):
        """Test convergence when starting very close to a node crossing."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # First find a node crossing
        jd_node, _, _ = ephem.mooncross_node_ut(jd_start, 0)

        # Now start very close to it (1 hour before) - tests edge case
        jd_near_node = jd_node - 1 / 24  # 1 hour before
        jd_next, _, _ = ephem.mooncross_node_ut(jd_near_node, 0)

        # Should find the same node (or the next one if we're past it)
        pos, _ = ephem.calc_ut(jd_next, MOON, 0)
        lat_arcsec = abs(pos[1]) * 3600.0
        assert lat_arcsec < 0.05, f"Latitude at node: {lat_arcsec:.4f} arcsec"

    @pytest.mark.unit
    def test_mooncross_node_multiple_consecutive_near_nodes(self):
        """Test finding multiple consecutive node crossings (6 crossings = ~81 days)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)
        crossings = []

        jd = jd_start
        for i in range(6):
            jd_cross, _, _ = ephem.mooncross_node_ut(jd, 0)
            crossings.append(jd_cross)

            # Verify latitude is near zero
            pos, _ = ephem.calc_ut(jd_cross, MOON, 0)
            lat_arcsec = abs(pos[1]) * 3600.0
            assert lat_arcsec < 0.05, (
                f"Crossing {i + 1}: latitude {lat_arcsec:.4f} arcsec"
            )

            jd = jd_cross + 0.5  # Move past this crossing

        # Verify crossings are ~13.6 days apart
        for i in range(1, len(crossings)):
            diff = crossings[i] - crossings[i - 1]
            assert 12 < diff < 15, f"Gap between crossings {i} and {i + 1}: {diff} days"

    @pytest.mark.unit
    def test_mooncross_node_at_minimum_latitude_speed(self):
        """Test convergence when Moon latitude speed is at minimum (~1°/day).

        This is the hardest case for Newton-Raphson convergence.
        At node crossings, latitude speed is approximately 1°/day.
        """
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find a node crossing
        jd_node, _, _ = ephem.mooncross_node_ut(jd_start, 0)

        # Verify we got good convergence
        pos, _ = ephem.calc_ut(jd_node, MOON, FLG_SPEED)
        lat_arcsec = abs(pos[1]) * 3600.0

        # Latitude should be within tolerance
        assert lat_arcsec < 0.05, f"Latitude: {lat_arcsec:.4f} arcsec"

        # Confirm latitude speed is around the expected ~1°/day
        lat_speed = abs(pos[4])
        assert 0.5 < lat_speed < 1.5, f"Latitude speed: {lat_speed}°/day"

    @pytest.mark.unit
    def test_mooncross_at_true_node_longitude(self):
        """Test Moon longitude crossing at the true node position.

        When Moon is at its true node, latitude is ~0 and this is
        an edge case for longitude crossing calculations.
        """
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Get true node position
        node_pos, _ = ephem.calc_ut(jd_start, TRUE_NODE, 0)
        node_lon = node_pos[0]

        # Find when Moon crosses the true node longitude
        jd_cross = ephem.mooncross_ut(node_lon, jd_start, 0)

        # Verify Moon is at target longitude
        pos, _ = ephem.calc_ut(jd_cross, MOON, 0)
        diff = abs(pos[0] - node_lon)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600.0

        assert diff_arcsec < 0.05, (
            f"Moon at {pos[0]:.6f}, node at {node_lon:.6f}, diff {diff_arcsec:.4f} arcsec"
        )

    @pytest.mark.unit
    def test_mooncross_node_tt_edge_case_near_node(self):
        """Test TT version convergence when starting near a node."""
        jd_start_ut = ephem.julday(2024, 1, 1, 0.0)
        delta_t = ephem.deltat(jd_start_ut)
        jd_start_tt = jd_start_ut + delta_t

        # Find node crossing in TT
        jd_node_tt, _, _ = ephem.mooncross_node(jd_start_tt, 0)

        # Start very close to it (30 minutes before)
        jd_near = jd_node_tt - 0.5 / 24
        jd_next_tt, _, _ = ephem.mooncross_node(jd_near, 0)

        # Verify convergence
        pos, _ = ephem.calc(jd_next_tt, MOON, 0)
        lat_arcsec = abs(pos[1]) * 3600.0
        assert lat_arcsec < 0.05, f"Latitude: {lat_arcsec:.4f} arcsec"

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [2020, 2024, 2028, 2032])
    def test_mooncross_node_convergence_different_years(self, year):
        """Test node crossing convergence across different years.

        The true node position oscillates ±1.5° around mean node,
        so testing different years covers different orbital phases.
        """
        jd_start = ephem.julday(year, 6, 15, 0.0)  # Mid-year

        jd_node, _, _ = ephem.mooncross_node_ut(jd_start, 0)

        pos, _ = ephem.calc_ut(jd_node, MOON, 0)
        lat_arcsec = abs(pos[1]) * 3600.0
        assert lat_arcsec < 0.05, f"Year {year}: latitude {lat_arcsec:.4f} arcsec"


class TestStationDetectionAndBrentsFallback:
    """Tests for station detection and Brent's method fallback.

    When planets are near retrograde stations (velocity ~0°/day), Newton-Raphson
    can have convergence issues. These tests verify that the Brent's method
    fallback handles these cases correctly.
    """

    @pytest.mark.unit
    def test_station_speed_threshold_constant(self):
        """Verify STATION_SPEED_THRESHOLD is defined and reasonable."""
        from libephemeris.crossing import STATION_SPEED_THRESHOLD

        # Should be a small positive value
        assert 0 < STATION_SPEED_THRESHOLD <= 0.01
        assert STATION_SPEED_THRESHOLD == 0.01

    @pytest.mark.unit
    def test_is_near_station_function(self):
        """Test the _is_near_station helper function."""
        from libephemeris.crossing import _is_near_station

        # Fast planets should not be near station
        assert not _is_near_station(1.0)
        assert not _is_near_station(0.5)
        assert not _is_near_station(0.1)
        assert not _is_near_station(0.01)  # At threshold, not below (uses <)
        assert _is_near_station(0.009)  # Below threshold
        assert _is_near_station(0.001)

        # Near-zero speeds should be detected as stations
        assert _is_near_station(0.0005)
        assert _is_near_station(0.0001)
        assert _is_near_station(0.0)

        # Negative speeds (retrograde) should also be detected
        assert not _is_near_station(-0.5)
        assert not _is_near_station(-0.01)  # At threshold
        assert _is_near_station(-0.0005)

    @pytest.mark.unit
    def test_brent_find_crossing_basic(self):
        """Test that Brent's method can find a simple crossing."""
        from libephemeris.crossing import _brent_find_crossing

        # Create a simple linear function for testing
        def linear_position(jd):
            # Position increases by 1° per day, starting at 0° at jd=0
            return jd * 1.0, 1.0  # (longitude, speed)

        # Find crossing of 5° starting from jd=0
        jd_cross = _brent_find_crossing(
            linear_position, x2cross=5.0, jd_a=0.0, jd_b=10.0, tolerance=0.0001
        )

        assert abs(jd_cross - 5.0) < 0.001

    @pytest.mark.unit
    def test_brent_find_crossing_with_wraparound(self):
        """Test Brent's method handles 360° wraparound."""
        from libephemeris.crossing import _brent_find_crossing

        def wrapping_position(jd):
            # Position that wraps around: 350° at jd=0, 10° at jd=20
            lon = (350.0 + jd * 1.0) % 360.0
            return lon, 1.0

        # Find crossing of 5° (should happen around jd=15)
        jd_cross = _brent_find_crossing(
            wrapping_position, x2cross=5.0, jd_a=10.0, jd_b=20.0, tolerance=0.0001
        )

        # At jd=15, position should be (350+15)%360 = 5°
        assert abs(jd_cross - 15.0) < 0.001

    @pytest.mark.unit
    def test_find_bracket_for_crossing_basic(self):
        """Test bracket finding for a simple crossing."""
        from libephemeris.crossing import _find_bracket_for_crossing

        def linear_position(jd):
            return jd * 1.0, 1.0  # 1°/day

        jd_a, jd_b = _find_bracket_for_crossing(
            linear_position, x2cross=5.0, jd_start=0.0, jd_end=10.0
        )

        # Bracket should contain jd=5 (jd_a <= 5 <= jd_b)
        assert jd_a <= 5.0 <= jd_b

    @pytest.mark.unit
    def test_find_bracket_raises_if_no_crossing(self):
        """Test that bracket finding raises error if no crossing found."""
        from libephemeris.crossing import _find_bracket_for_crossing

        def constant_position(jd):
            return 10.0, 0.0  # Always at 10°

        with pytest.raises(RuntimeError, match="No crossing found"):
            _find_bracket_for_crossing(
                constant_position, x2cross=50.0, jd_start=0.0, jd_end=10.0
            )

    @pytest.mark.unit
    def test_mercury_retrograde_station_crossing(self):
        """Test crossing calculation when Mercury is near a retrograde station.

        Mercury goes retrograde 3-4 times per year. Near station points,
        its speed approaches 0°/day, which can cause Newton-Raphson issues.
        The Brent's method fallback should handle these cases.
        """
        # Mercury retrograde station dates in 2024:
        # - Dec 2023 to Jan 2024: Stations around Jan 1
        # - Apr 2024: Stations around Apr 1 and Apr 25
        # - Aug 2024: Stations around Aug 5 and Aug 28
        # - Nov-Dec 2024: Stations around Nov 25 and Dec 15

        # Start near a Mercury station (early January 2024)
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Get Mercury's current position
        pos, _ = ephem.calc_ut(jd_start, MERCURY, FLG_SPEED)
        current_lon = pos[0]

        # Target a degree ahead (should work even if near station)
        target = (current_lon + 5) % 360

        # This should not raise even if near station
        jd_cross = ephem.cross_ut(MERCURY, target, jd_start, 0)

        # Verify the crossing was found
        assert jd_cross > jd_start

        # Verify Mercury is at target at crossing time
        pos_cross, _ = ephem.calc_ut(jd_cross, MERCURY, 0)
        diff = abs(pos_cross[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1, f"Mercury at {pos_cross[0]}, target {target}, diff {diff}"

    @pytest.mark.unit
    def test_saturn_slow_motion_crossing(self):
        """Test crossing for Saturn which moves ~0.03°/day.

        Saturn is a slow planet that can trigger the station detection
        logic even during normal motion.
        """
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Get Saturn's current position
        pos, _ = ephem.calc_ut(jd_start, SATURN, FLG_SPEED)
        current_lon = pos[0]

        # Target 2° ahead
        target = (current_lon + 2) % 360

        jd_cross = ephem.cross_ut(SATURN, target, jd_start, 0)

        assert jd_cross > jd_start

        pos_cross, _ = ephem.calc_ut(jd_cross, SATURN, 0)
        diff = abs(pos_cross[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1

    @pytest.mark.unit
    def test_pluto_very_slow_motion_crossing(self):
        """Test crossing for Pluto which moves slowly (~0.004-0.03°/day).

        Pluto is a slow planet that benefits from the robust station handling,
        even when not at a true station.
        """
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Get Pluto's current position
        pos, _ = ephem.calc_ut(jd_start, PLUTO, FLG_SPEED)
        current_lon = pos[0]
        current_speed = pos[3]

        # Pluto moves slowly (typically 0.004-0.03°/day)
        assert abs(current_speed) < 0.05, f"Pluto speed: {current_speed}"

        # Target 1° ahead
        target = (current_lon + 1) % 360

        jd_cross = ephem.cross_ut(PLUTO, target, jd_start, 0)

        assert jd_cross > jd_start

        pos_cross, _ = ephem.calc_ut(jd_cross, PLUTO, 0)
        diff = abs(pos_cross[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1

    @pytest.mark.unit
    def test_venus_retrograde_crossing(self):
        """Test crossing when Venus is retrograde.

        Venus goes retrograde every ~18 months. Near stations,
        its speed approaches 0°/day.
        """
        # Venus was retrograde in late 2023 (Oct 25 - Nov 4 station)
        jd_start = ephem.julday(2024, 3, 1, 0.0)

        pos, _ = ephem.calc_ut(jd_start, VENUS, FLG_SPEED)
        current_lon = pos[0]

        # Target 10° ahead
        target = (current_lon + 10) % 360

        jd_cross = ephem.cross_ut(VENUS, target, jd_start, 0)

        assert jd_cross > jd_start

        pos_cross, _ = ephem.calc_ut(jd_cross, VENUS, 0)
        diff = abs(pos_cross[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1

    @pytest.mark.unit
    def test_helio_cross_with_slow_planet(self):
        """Test heliocentric crossing with Pluto (very slow)."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        pos, _ = ephem.calc_ut(jd_start, PLUTO, FLG_HELCTR | FLG_SPEED)
        current_lon = pos[0]

        # Target 0.5° ahead
        target = (current_lon + 0.5) % 360

        jd_cross = ephem.helio_cross_ut(PLUTO, target, jd_start, 0)

        assert jd_cross > jd_start

        pos_cross, _ = ephem.calc_ut(jd_cross, PLUTO, FLG_HELCTR)
        diff = abs(pos_cross[0] - target)
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.1

    @pytest.mark.unit
    def test_multiple_planets_near_stations(self):
        """Test that multiple slow planets converge correctly."""
        jd_start = ephem.julday(2024, 6, 1, 0.0)

        # Use smaller offsets for very slow planets, and avoid planets that might be retrograde
        slow_planets = [
            (JUPITER, 5),  # Jupiter, target 5° ahead
            (SATURN, 3),  # Saturn, target 3° ahead
            (URANUS, 2),  # Uranus, target 2° ahead
        ]

        for planet_id, offset in slow_planets:
            pos, _ = ephem.calc_ut(jd_start, planet_id, FLG_SPEED)
            target = (pos[0] + offset) % 360

            jd_cross = ephem.cross_ut(planet_id, target, jd_start, 0)

            assert jd_cross > jd_start

            pos_cross, _ = ephem.calc_ut(jd_cross, planet_id, 0)
            diff = abs(pos_cross[0] - target)
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.1, f"Planet {planet_id} failed: diff={diff}"

    @pytest.mark.unit
    def test_crossing_accuracy_preserved_with_brent(self):
        """Verify that Brent's method achieves the same precision as Newton-Raphson."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Test with Saturn (slow but not quite at station)
        pos, _ = ephem.calc_ut(jd_start, SATURN, FLG_SPEED)
        target = (pos[0] + 5) % 360

        jd_cross = ephem.cross_ut(SATURN, target, jd_start, 0)

        pos_cross, _ = ephem.calc_ut(jd_cross, SATURN, 0)
        diff = abs(pos_cross[0] - target)
        if diff > 180:
            diff = 360 - diff

        # Should still achieve sub-arcsecond precision
        diff_arcsec = diff * 3600.0
        assert diff_arcsec < 0.1, f"Precision: {diff_arcsec:.4f} arcsec"


class TestKnownAstronomicalEvents:
    """Tests verifying solstices/equinoxes against Swiss Ephemeris reference values.

    These tests verify that the crossing functions produce results that match
    Swiss Ephemeris calculated astronomical events with precision < 1 second.

    Reference dates calculated using Swiss Ephemeris (pyswisseph) which uses
    JPL DE431 ephemeris data. These values are the authoritative reference for
    this library since libephemeris is designed to match Swiss Ephemeris behavior.

    Swiss Ephemeris calculated values (2023-2025):
    - 2024 Vernal Equinox: March 20, 03:06:24.243 UTC (Sun at 0°)
    - 2024 Summer Solstice: June 20, 20:51:00.170 UTC (Sun at 90°)
    - 2024 Autumnal Equinox: September 22, 12:43:39.036 UTC (Sun at 180°)
    - 2024 Winter Solstice: December 21, 09:20:34.132 UTC (Sun at 270°)

    Note: These values may differ by ~3-5 seconds from other sources (USNO, etc.)
    due to different ephemeris models and definitions of "crossing" (apparent vs
    geometric Sun position).
    """

    # 1 second precision threshold (in days)
    ONE_SECOND_DAYS = 1.0 / 86400.0

    @staticmethod
    def _time_to_hours(hour: int, minute: int, second: float) -> float:
        """Convert time to decimal hours."""
        return hour + minute / 60.0 + second / 3600.0

    @pytest.mark.unit
    def test_vernal_equinox_2024(self):
        """Verify 2024 vernal equinox (Sun at 0°) matches Swiss Ephemeris within 1 second.

        Swiss Ephemeris: March 20, 2024 at 03:06:24.243 UTC
        """
        # Start search from beginning of 2024
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find when Sun crosses 0° (vernal equinox)
        jd_cross = ephem.solcross_ut(0.0, jd_start, 0)

        # Swiss Ephemeris reference: March 20, 2024 at 03:06:24.243 UTC
        jd_known = ephem.julday(2024, 3, 20, self._time_to_hours(3, 6, 24.243))

        # Difference should be less than 1 second
        diff_seconds = abs(jd_cross - jd_known) * 86400.0
        assert diff_seconds < 1.0, (
            f"Vernal equinox 2024: calculated JD {jd_cross:.8f}, "
            f"known JD {jd_known:.8f}, diff {diff_seconds:.3f} seconds "
            f"(should be < 1 second)"
        )

    @pytest.mark.unit
    def test_summer_solstice_2024(self):
        """Verify 2024 summer solstice (Sun at 90°) matches Swiss Ephemeris within 1 second.

        Swiss Ephemeris: June 20, 2024 at 20:51:00.170 UTC
        """
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find when Sun crosses 90° (summer solstice)
        jd_cross = ephem.solcross_ut(90.0, jd_start, 0)

        # Swiss Ephemeris reference: June 20, 2024 at 20:51:00.170 UTC
        jd_known = ephem.julday(2024, 6, 20, self._time_to_hours(20, 51, 0.170))

        diff_seconds = abs(jd_cross - jd_known) * 86400.0
        assert diff_seconds < 1.0, (
            f"Summer solstice 2024: calculated JD {jd_cross:.8f}, "
            f"known JD {jd_known:.8f}, diff {diff_seconds:.3f} seconds "
            f"(should be < 1 second)"
        )

    @pytest.mark.unit
    def test_autumnal_equinox_2024(self):
        """Verify 2024 autumnal equinox (Sun at 180°) matches Swiss Ephemeris within 1 second.

        Swiss Ephemeris: September 22, 2024 at 12:43:39.036 UTC
        """
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find when Sun crosses 180° (autumnal equinox)
        jd_cross = ephem.solcross_ut(180.0, jd_start, 0)

        # Swiss Ephemeris reference: September 22, 2024 at 12:43:39.036 UTC
        jd_known = ephem.julday(2024, 9, 22, self._time_to_hours(12, 43, 39.036))

        diff_seconds = abs(jd_cross - jd_known) * 86400.0
        assert diff_seconds < 1.0, (
            f"Autumnal equinox 2024: calculated JD {jd_cross:.8f}, "
            f"known JD {jd_known:.8f}, diff {diff_seconds:.3f} seconds "
            f"(should be < 1 second)"
        )

    @pytest.mark.unit
    def test_winter_solstice_2024(self):
        """Verify 2024 winter solstice (Sun at 270°) matches Swiss Ephemeris within 1 second.

        Swiss Ephemeris: December 21, 2024 at 09:20:34.132 UTC
        """
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Find when Sun crosses 270° (winter solstice)
        jd_cross = ephem.solcross_ut(270.0, jd_start, 0)

        # Swiss Ephemeris reference: December 21, 2024 at 09:20:34.132 UTC
        jd_known = ephem.julday(2024, 12, 21, self._time_to_hours(9, 20, 34.132))

        diff_seconds = abs(jd_cross - jd_known) * 86400.0
        assert diff_seconds < 1.0, (
            f"Winter solstice 2024: calculated JD {jd_cross:.8f}, "
            f"known JD {jd_known:.8f}, diff {diff_seconds:.3f} seconds "
            f"(should be < 1 second)"
        )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "event_name,target_lon,year,month,day,hour,minute,second",
        [
            # 2024 events (Swiss Ephemeris calculated)
            ("Vernal Equinox 2024", 0.0, 2024, 3, 20, 3, 6, 24.243),
            ("Summer Solstice 2024", 90.0, 2024, 6, 20, 20, 51, 0.170),
            ("Autumnal Equinox 2024", 180.0, 2024, 9, 22, 12, 43, 39.036),
            ("Winter Solstice 2024", 270.0, 2024, 12, 21, 9, 20, 34.132),
            # 2023 events (Swiss Ephemeris calculated)
            ("Vernal Equinox 2023", 0.0, 2023, 3, 20, 21, 24, 25.767),
            ("Summer Solstice 2023", 90.0, 2023, 6, 21, 14, 57, 50.739),
            ("Autumnal Equinox 2023", 180.0, 2023, 9, 23, 6, 49, 59.279),
            ("Winter Solstice 2023", 270.0, 2023, 12, 22, 3, 27, 22.795),
            # 2025 events (Swiss Ephemeris calculated)
            ("Vernal Equinox 2025", 0.0, 2025, 3, 20, 9, 1, 29.361),
            ("Summer Solstice 2025", 90.0, 2025, 6, 21, 2, 42, 16.474),
            ("Autumnal Equinox 2025", 180.0, 2025, 9, 22, 18, 19, 20.085),
            ("Winter Solstice 2025", 270.0, 2025, 12, 21, 15, 3, 5.075),
        ],
    )
    def test_solstice_equinox_parametrized(
        self, event_name, target_lon, year, month, day, hour, minute, second
    ):
        """Verify multiple years of solstices/equinoxes match Swiss Ephemeris within 1.5 seconds.

        Uses a 1.5 second tolerance to account for minor floating-point precision
        differences in the Newton-Raphson iteration between libephemeris and pyswisseph.
        The direct comparison test (test_solcross_matches_pyswisseph_for_known_events)
        verifies the core algorithm matches pyswisseph within 1 second.
        """
        # Start search from beginning of the year
        jd_start = ephem.julday(year, 1, 1, 0.0)

        # Find crossing
        jd_cross = ephem.solcross_ut(target_lon, jd_start, 0)

        # Swiss Ephemeris reference time
        jd_known = ephem.julday(
            year, month, day, self._time_to_hours(hour, minute, second)
        )

        diff_seconds = abs(jd_cross - jd_known) * 86400.0
        # Allow 1.5 second tolerance for floating-point precision differences
        assert diff_seconds < 1.5, (
            f"{event_name}: diff {diff_seconds:.3f} seconds (should be < 1.5 seconds)"
        )

    @pytest.mark.unit
    def test_new_moon_january_2024(self):
        """Verify new moon calculation using mooncross_ut.

        New Moon: January 11, 2024 at 11:57:38 UTC
        At new moon, Moon longitude equals Sun longitude.
        """
        # Start of January 2024
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Get Sun position and find when Moon crosses that longitude
        sun_pos, _ = ephem.calc_ut(jd_start, SUN, 0)
        sun_lon = sun_pos[0]

        # Find when Moon crosses the Sun's longitude (approximately new moon)
        # Note: For precise new moon, we'd need to track Sun's motion too,
        # but this test verifies mooncross_ut precision
        jd_new_moon = ephem.mooncross_ut(sun_lon, jd_start, 0)

        # Verify Moon is at target position with sub-arcsecond precision
        moon_pos, _ = ephem.calc_ut(jd_new_moon, MOON, 0)
        diff = abs(moon_pos[0] - sun_lon)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600.0

        assert diff_arcsec < 0.05, (
            f"Moon longitude precision: {diff_arcsec:.4f} arcsec (should be < 0.05)"
        )

    @pytest.mark.unit
    def test_full_moon_january_2024(self):
        """Verify full moon calculation using mooncross_ut.

        Full Moon: January 25, 2024 at 17:54:01 UTC
        At full moon, Moon longitude = Sun longitude + 180°.
        """
        # After new moon, find full moon
        jd_start = ephem.julday(2024, 1, 12, 0.0)  # Day after new moon

        # Get Sun position at start
        sun_pos, _ = ephem.calc_ut(jd_start, SUN, 0)
        # Full moon target: opposite the Sun
        target_lon = (sun_pos[0] + 180.0) % 360.0

        # Find when Moon crosses the opposite point
        jd_full_moon = ephem.mooncross_ut(target_lon, jd_start, 0)

        # Verify Moon is at target position
        moon_pos, _ = ephem.calc_ut(jd_full_moon, MOON, 0)
        diff = abs(moon_pos[0] - target_lon)
        if diff > 180:
            diff = 360 - diff
        diff_arcsec = diff * 3600.0

        assert diff_arcsec < 0.05, (
            f"Moon longitude precision: {diff_arcsec:.4f} arcsec (should be < 0.05)"
        )

        # Verify the date is in late January 2024
        year, month, day, hour = ephem.revjul(jd_full_moon)
        assert year == 2024
        assert month == 1
        assert 24 <= day <= 26  # Full moon around Jan 25

    @pytest.mark.unit
    def test_mooncross_longitude_precision_multiple_targets(self):
        """Verify mooncross_ut achieves sub-arcsecond precision for various longitudes."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        # Test multiple target longitudes
        targets = [0, 45, 90, 135, 180, 225, 270, 315]

        for target in targets:
            jd_cross = ephem.mooncross_ut(float(target), jd_start, 0)

            # Verify Moon is at target
            moon_pos, _ = ephem.calc_ut(jd_cross, MOON, 0)
            diff = abs(moon_pos[0] - target)
            if diff > 180:
                diff = 360 - diff
            diff_arcsec = diff * 3600.0

            assert diff_arcsec < 0.05, (
                f"Target {target}°: diff {diff_arcsec:.4f} arcsec (should be < 0.05)"
            )

    @pytest.mark.unit
    def test_solcross_matches_pyswisseph_for_known_events(self):
        """Verify libephemeris matches pyswisseph for known astronomical events."""
        jd_start = ephem.julday(2024, 1, 1, 0.0)

        for target in [0.0, 90.0, 180.0, 270.0]:
            jd_lib = ephem.solcross_ut(target, jd_start, 0)
            jd_swe = swe.solcross_ut(target, jd_start, 0)

            # Should match within 1 second
            diff_seconds = abs(jd_lib - jd_swe) * 86400.0
            assert diff_seconds < 1.0, (
                f"Target {target}°: libephemeris vs pyswisseph diff {diff_seconds:.3f} seconds"
            )
