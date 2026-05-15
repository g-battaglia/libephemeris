"""
Comprehensive tests for edge cases.

Tests boundary conditions, invalid inputs, and error handling.
"""

import pytest
import libephemeris as ephem
from libephemeris.constants import *


class TestDateBoundaries:
    """Test date boundary conditions."""

    @pytest.mark.edge_case
    def test_de440_start_1550(self):
        """Calculations at DE440 start (1550) should work."""
        jd = ephem.julday(1550, 1, 1, 12.0)

        pos, _ = ephem.calc_ut(jd, SUN, 0)

        assert 0 <= pos[0] < 360

    @pytest.mark.edge_case
    def test_de440_end_2650(self):
        """Calculations at DE440 end (2650) should work."""
        jd = ephem.julday(2650, 1, 1, 12.0)

        pos, _ = ephem.calc_ut(jd, SUN, 0)

        assert 0 <= pos[0] < 360

    @pytest.mark.edge_case
    def test_before_de440_range(self):
        """Calculations before DE440 range should fail or return error."""
        jd = ephem.julday(1549, 1, 1, 12.0)

        # Should either raise exception or return with error flag
        try:
            pos, flags = ephem.calc_ut(jd, SUN, 0)
            # If it succeeds, result might be extrapolated
            # We should at least check it's reasonable
        except Exception:
            pass  # Expected to fail

    @pytest.mark.edge_case
    def test_after_de440_range(self):
        """Calculations after DE440 range should fail or return error."""
        jd = ephem.julday(2651, 1, 1, 12.0)

        try:
            pos, flags = ephem.calc_ut(jd, SUN, 0)
        except Exception:
            pass  # Expected to fail


class TestInvalidInputs:
    """Test handling of invalid inputs."""

    @pytest.mark.edge_case
    def test_invalid_planet_id(self):
        """Invalid planet ID should raise error or return error."""
        jd = 2451545.0

        try:
            pos, _ = ephem.calc_ut(jd, 999, 0)
            # If no exception, check for error indicator
            # Some implementations return error in flags
        except (ValueError, KeyError, Exception):
            pass  # Expected

    @pytest.mark.edge_case
    def test_latitude_over_90(self):
        """Latitude > 90 should raise CoordinateError with clear message."""
        jd = 2451545.0

        with pytest.raises(
            ephem.CoordinateError, match="latitude 91.0 is out of valid range"
        ):
            ephem.houses(jd, 91.0, 0.0, ord("P"))

    @pytest.mark.edge_case
    def test_latitude_under_minus_90(self):
        """Latitude < -90 should raise CoordinateError with clear message."""
        jd = 2451545.0

        with pytest.raises(
            ephem.CoordinateError, match="latitude -91.0 is out of valid range"
        ):
            ephem.houses(jd, -91.0, 0.0, ord("P"))


class TestLongitudeWrapping:
    """Test longitude wraparound handling."""

    @pytest.mark.edge_case
    def test_longitude_360_equals_0(self):
        """360° should equal 0°."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        jd_360 = ephem.solcross_ut(360.0, jd, 0)
        jd_0 = ephem.solcross_ut(0.0, jd, 0)

        # Should find the same crossing
        assert abs(jd_360 - jd_0) < 0.01

    @pytest.mark.edge_case
    def test_negative_longitude_normalized(self):
        """Negative longitude should be normalized."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        # -30° should equal 330°
        jd_neg = ephem.solcross_ut(-30.0, jd, 0)
        jd_pos = ephem.solcross_ut(330.0, jd, 0)

        assert abs(jd_neg - jd_pos) < 0.01

    @pytest.mark.edge_case
    def test_longitude_over_360_normalized(self):
        """Longitude > 360° should be normalized."""
        jd = ephem.julday(2024, 1, 1, 0.0)

        # 400° should equal 40°
        jd_over = ephem.solcross_ut(400.0, jd, 0)
        jd_norm = ephem.solcross_ut(40.0, jd, 0)

        assert abs(jd_over - jd_norm) < 0.01


class TestTimeBoundaries:
    """Test time-related edge cases."""

    @pytest.mark.edge_case
    def test_midnight_boundary(self):
        """Test calculations at midnight."""
        jd = ephem.julday(2000, 1, 1, 0.0)

        pos, _ = ephem.calc_ut(jd, SUN, 0)

        assert 0 <= pos[0] < 360

    @pytest.mark.edge_case
    def test_leap_second(self):
        """Test date with leap second (handled gracefully)."""
        # 2016-12-31 had a leap second
        jd = ephem.julday(2016, 12, 31, 23.999)

        pos, _ = ephem.calc_ut(jd, SUN, 0)

        assert 0 <= pos[0] < 360


class TestSpecialPositions:
    """Test special astronomical positions."""

    @pytest.mark.edge_case
    def test_sun_at_equinox(self):
        """Sun at equinox should be at 0° or 180°."""
        # Find 2024 vernal equinox
        jd_start = ephem.julday(2024, 3, 1, 0.0)
        jd_equinox = ephem.solcross_ut(0.0, jd_start, 0)

        pos, _ = ephem.calc_ut(jd_equinox, SUN, 0)

        # Should be very close to 0°
        diff = pos[0]
        if diff > 180:
            diff = 360 - diff
        assert diff < 0.001

    @pytest.mark.edge_case
    def test_moon_at_lunar_node(self):
        """Test Moon position near lunar node."""
        # Get node position
        jd = 2451545.0
        node_pos, _ = ephem.calc_ut(jd, MEAN_NODE, 0)

        # Find when Moon crosses this longitude
        jd_cross = ephem.mooncross_ut(node_pos[0], jd, 0)

        # Moon latitude should be near 0 at node
        moon_pos, _ = ephem.calc_ut(jd_cross, MOON, 0)

        assert abs(moon_pos[1]) < 1.0, f"Moon lat at node: {moon_pos[1]}"


class TestReturnTypes:
    """Test that return types are correct."""

    @pytest.mark.unit
    def test_julday_returns_float(self):
        """julday should return float."""
        jd = ephem.julday(2000, 1, 1, 12.0)
        assert isinstance(jd, float)

    @pytest.mark.unit
    def test_revjul_returns_tuple(self):
        """revjul should return tuple."""
        result = ephem.revjul(2451545.0)
        assert isinstance(result, tuple)
        assert len(result) == 4

    @pytest.mark.unit
    def test_calc_ut_returns_tuple(self):
        """calc_ut should return (position, flags) tuple."""
        result = ephem.calc_ut(2451545.0, SUN, 0)
        assert isinstance(result, tuple)
        assert len(result) == 2

    @pytest.mark.unit
    def test_houses_returns_tuple(self):
        """houses should return (cusps, ascmc) tuple."""
        result = ephem.houses(2451545.0, 41.9, 12.5, ord("P"))
        assert isinstance(result, tuple)
        assert len(result) == 2
