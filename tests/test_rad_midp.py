"""
Tests for the rad_midp angle midpoint function.

Tests verify that rad_midp correctly calculates the midpoint between two angles
handling wraparound at 2*pi radians.
"""

import math
import pytest
import libephemeris as ephem


TWO_PI = 2.0 * math.pi


class TestRadMidpBasic:
    """Basic functionality tests for rad_midp."""

    def test_rad_midp_exported(self):
        """Test that rad_midp is exported from the package."""
        assert hasattr(ephem, "rad_midp")
        assert callable(ephem.rad_midp)

    def test_rad_midp_simple(self):
        """Test midpoint of simple angles."""
        result = ephem.rad_midp(0, math.pi / 2)
        assert result == pytest.approx(math.pi / 4)

    def test_rad_midp_same_angle(self):
        """Test midpoint when both angles are the same."""
        result = ephem.rad_midp(math.pi / 4, math.pi / 4)
        assert result == pytest.approx(math.pi / 4)

    def test_rad_midp_returns_float(self):
        """Test that rad_midp returns a float."""
        result = ephem.rad_midp(0, math.pi / 2)
        assert isinstance(result, float)


class TestRadMidpWraparound:
    """Tests for wraparound handling at 2*pi radians."""

    def test_rad_midp_wraparound_near_2pi(self):
        """Test midpoint between angles near 2*pi and 0."""
        # Equivalent to 350 and 10 degrees
        a = math.radians(350)
        b = math.radians(10)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_rad_midp_wraparound_reverse(self):
        """Test midpoint between 10 and 350 degrees (in radians)."""
        a = math.radians(10)
        b = math.radians(350)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_rad_midp_wraparound_355_5(self):
        """Test midpoint between 355 and 5 degrees (in radians)."""
        a = math.radians(355)
        b = math.radians(5)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_rad_midp_wraparound_340_20(self):
        """Test midpoint between 340 and 20 degrees (in radians)."""
        a = math.radians(340)
        b = math.radians(20)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_rad_midp_wraparound_270_90(self):
        """Test midpoint between 270 and 90 degrees (in radians)."""
        # 3*pi/2 to pi/2 is 180 degrees either way
        a = 3 * math.pi / 2
        b = math.pi / 2
        result = ephem.rad_midp(a, b)
        # Both paths are equal length, algorithm chooses one consistently
        assert result == pytest.approx(0.0) or result == pytest.approx(math.pi)


class TestRadMidpOppositeAngles:
    """Tests for opposite angles (pi radians apart)."""

    def test_rad_midp_0_pi(self):
        """Test midpoint between 0 and pi."""
        result = ephem.rad_midp(0, math.pi)
        assert result == pytest.approx(math.pi / 2)

    def test_rad_midp_pi_0(self):
        """Test midpoint between pi and 0.

        When both arcs are equally long (π), the reference ephemeris chooses the
        positive (clockwise) arc, yielding 3π/2 not π/2.
        """
        result = ephem.rad_midp(math.pi, 0)
        assert result == pytest.approx(3 * math.pi / 2)

    def test_rad_midp_pi2_3pi2(self):
        """Test midpoint between pi/2 and 3*pi/2."""
        result = ephem.rad_midp(math.pi / 2, 3 * math.pi / 2)
        # Could be 0 or pi, both are equidistant
        assert result == pytest.approx(0.0) or result == pytest.approx(math.pi)

    def test_rad_midp_170_190_degrees(self):
        """Test midpoint between 170 and 190 degrees (in radians)."""
        a = math.radians(170)
        b = math.radians(190)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(math.pi)


class TestRadMidpNegativeAngles:
    """Tests for negative angle handling."""

    def test_rad_midp_negative_small(self):
        """Test midpoint between negative and positive small angles."""
        a = math.radians(-10)
        b = math.radians(10)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_rad_midp_negative_90_0(self):
        """Test midpoint between -90 and 0 degrees (in radians)."""
        a = math.radians(-90)
        b = 0.0
        result = ephem.rad_midp(a, b)
        # -90 degrees = 270 degrees, midpoint with 0 is 315 degrees
        assert result == pytest.approx(math.radians(315))

    def test_rad_midp_both_negative(self):
        """Test midpoint between two negative angles."""
        a = math.radians(-45)
        b = math.radians(-90)
        # -45 = 315 degrees, -90 = 270 degrees
        # Midpoint should be 292.5 degrees
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(math.radians(292.5))


class TestRadMidpLargeAngles:
    """Tests for large angle handling."""

    def test_rad_midp_greater_than_2pi(self):
        """Test midpoint with angle > 2*pi."""
        # 370 degrees = 10 degrees
        a = math.radians(370)
        b = math.radians(10)
        result = ephem.rad_midp(a, b)
        # Both normalize to ~10 degrees
        assert result == pytest.approx(math.radians(10))

    def test_rad_midp_multiple_rotations(self):
        """Test midpoint with multiple rotations."""
        # 4*pi = 0
        a = 4 * math.pi
        b = 0.0
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_rad_midp_large_angles(self):
        """Test midpoint with large angles."""
        # 710 = 350 degrees, 730 = 10 degrees, midpoint should be 0
        a = math.radians(710)
        b = math.radians(730)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=1e-10)


class TestRadMidpSpecificValues:
    """Tests for specific known values."""

    def test_rad_midp_0_0(self):
        """Test midpoint of 0 and 0."""
        result = ephem.rad_midp(0, 0)
        assert result == pytest.approx(0.0)

    def test_rad_midp_2pi_0(self):
        """Test midpoint of 2*pi and 0."""
        result = ephem.rad_midp(TWO_PI, 0)
        assert result == pytest.approx(0.0, abs=1e-10)

    def test_rad_midp_45_135_degrees(self):
        """Test midpoint of 45 and 135 degrees (in radians)."""
        a = math.radians(45)
        b = math.radians(135)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(math.radians(90))

    def test_rad_midp_315_45_degrees(self):
        """Test midpoint of 315 and 45 degrees (in radians)."""
        a = math.radians(315)
        b = math.radians(45)
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=1e-10)


class TestRadMidpSymmetry:
    """Tests for symmetry properties."""

    @pytest.mark.parametrize(
        "a,b",
        [
            (0, math.pi / 2),
            (math.pi / 2, math.pi),
            (math.pi, 3 * math.pi / 2),
            (3 * math.pi / 2, TWO_PI),
            (math.radians(10), math.radians(350)),
            (math.radians(30), math.radians(60)),
            (math.radians(100), math.radians(200)),
        ],
    )
    def test_rad_midp_order_invariance(self, a, b):
        """Test that rad_midp(a, b) == rad_midp(b, a)."""
        result_ab = ephem.rad_midp(a, b)
        result_ba = ephem.rad_midp(b, a)
        assert result_ab == pytest.approx(result_ba, abs=1e-10)


class TestRadMidpEdgeCases:
    """Edge case tests for rad_midp."""

    def test_rad_midp_very_small_difference(self):
        """Test with very small angle difference."""
        a = 0.0001
        b = 0.0002
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.00015, abs=1e-10)

    def test_rad_midp_near_2pi_boundary(self):
        """Test near the 2*pi boundary."""
        a = TWO_PI - 0.001
        b = 0.001
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(0.0, abs=0.01)

    def test_rad_midp_fractional_angles(self):
        """Test with fractional angles."""
        a = 0.5
        b = 1.5
        result = ephem.rad_midp(a, b)
        assert result == pytest.approx(1.0)


class TestRadMidpConsistencyWithDegMidp:
    """Tests that rad_midp is consistent with deg_midp."""

    @pytest.mark.parametrize(
        "deg_a,deg_b",
        [
            (0, 90),
            (350, 10),
            (180, 0),
            (45, 135),
            (270, 90),
            (315, 45),
        ],
    )
    def test_rad_midp_matches_deg_midp(self, deg_a, deg_b):
        """Test that rad_midp gives same result as deg_midp when converted."""
        deg_result = ephem.deg_midp(deg_a, deg_b)
        rad_result = ephem.rad_midp(math.radians(deg_a), math.radians(deg_b))
        # Convert rad_result to degrees for comparison
        rad_result_deg = math.degrees(rad_result)
        # Normalize to [0, 360)
        rad_result_deg = rad_result_deg % 360.0
        assert rad_result_deg == pytest.approx(deg_result, abs=1e-10)
