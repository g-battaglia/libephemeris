"""
Tests for the deg_midp angle midpoint function.

Tests verify that deg_midp correctly calculates the midpoint between two angles
handling wraparound at 360 degrees.
"""

import pytest
import libephemeris as ephem


class TestDegMidpBasic:
    """Basic functionality tests for deg_midp."""

    def test_deg_midp_exported(self):
        """Test that deg_midp is exported from the package."""
        assert hasattr(ephem, "deg_midp")
        assert callable(ephem.deg_midp)

    def test_deg_midp_simple(self):
        """Test midpoint of simple angles."""
        result = ephem.deg_midp(0, 90)
        assert result == pytest.approx(45.0)

    def test_deg_midp_same_angle(self):
        """Test midpoint when both angles are the same."""
        result = ephem.deg_midp(45, 45)
        assert result == pytest.approx(45.0)

    def test_deg_midp_returns_float(self):
        """Test that deg_midp returns a float."""
        result = ephem.deg_midp(0, 90)
        assert isinstance(result, float)


class TestDegMidpWraparound:
    """Tests for wraparound handling at 360 degrees."""

    def test_deg_midp_wraparound_350_10(self):
        """Test midpoint between 350 and 10 is 0 (not 180)."""
        result = ephem.deg_midp(350, 10)
        assert result == pytest.approx(0.0)

    def test_deg_midp_wraparound_10_350(self):
        """Test midpoint between 10 and 350 is 0 (order shouldn't matter for result)."""
        result = ephem.deg_midp(10, 350)
        assert result == pytest.approx(0.0)

    def test_deg_midp_wraparound_355_5(self):
        """Test midpoint between 355 and 5."""
        result = ephem.deg_midp(355, 5)
        assert result == pytest.approx(0.0)

    def test_deg_midp_wraparound_340_20(self):
        """Test midpoint between 340 and 20."""
        result = ephem.deg_midp(340, 20)
        assert result == pytest.approx(0.0)

    def test_deg_midp_wraparound_270_90(self):
        """Test midpoint between 270 and 90 via 0/360 is 0 or 180 (equidistant)."""
        # 270 to 90 is 180 degrees either way, so either 0 or 180 is valid
        # Our algorithm should return 0 (going via shorter path)
        result = ephem.deg_midp(270, 90)
        # Both paths are equal length, algorithm chooses one consistently
        assert result == pytest.approx(0.0) or result == pytest.approx(180.0)


class TestDegMidpOppositeAngles:
    """Tests for opposite angles (180 degrees apart)."""

    def test_deg_midp_0_180(self):
        """Test midpoint between 0 and 180."""
        result = ephem.deg_midp(0, 180)
        assert result == pytest.approx(90.0)

    def test_deg_midp_180_0(self):
        """Test midpoint between 180 and 0.

        When both arcs are equally long (180°), the reference ephemeris chooses the
        positive (clockwise) arc, yielding 270° not 90°.
        """
        result = ephem.deg_midp(180, 0)
        assert result == pytest.approx(270.0)

    def test_deg_midp_90_270(self):
        """Test midpoint between 90 and 270."""
        result = ephem.deg_midp(90, 270)
        # Could be 0 or 180, both are equidistant
        assert result == pytest.approx(0.0) or result == pytest.approx(180.0)

    def test_deg_midp_170_190(self):
        """Test midpoint between 170 and 190."""
        result = ephem.deg_midp(170, 190)
        assert result == pytest.approx(180.0)


class TestDegMidpNegativeAngles:
    """Tests for negative angle handling."""

    def test_deg_midp_negative_10_10(self):
        """Test midpoint between -10 and 10."""
        result = ephem.deg_midp(-10, 10)
        assert result == pytest.approx(0.0)

    def test_deg_midp_negative_90_0(self):
        """Test midpoint between -90 and 0."""
        result = ephem.deg_midp(-90, 0)
        assert result == pytest.approx(315.0)

    def test_deg_midp_both_negative(self):
        """Test midpoint between two negative angles."""
        result = ephem.deg_midp(-45, -90)
        # -45 = 315, -90 = 270
        # Midpoint should be 292.5
        assert result == pytest.approx(292.5)


class TestDegMidpLargeAngles:
    """Tests for large angle handling."""

    def test_deg_midp_370_10(self):
        """Test midpoint with angle > 360."""
        result = ephem.deg_midp(370, 10)
        # 370 = 10, so midpoint of 10 and 10 is 10
        assert result == pytest.approx(10.0)

    def test_deg_midp_720_0(self):
        """Test midpoint with multiple rotations."""
        result = ephem.deg_midp(720, 0)
        # 720 = 0, so midpoint of 0 and 0 is 0
        assert result == pytest.approx(0.0)

    def test_deg_midp_large_angles(self):
        """Test midpoint with large angles."""
        result = ephem.deg_midp(710, 730)
        # 710 = 350, 730 = 10, midpoint should be 0
        assert result == pytest.approx(0.0)


class TestDegMidpSpecificValues:
    """Tests for specific known values."""

    def test_deg_midp_0_0(self):
        """Test midpoint of 0 and 0."""
        result = ephem.deg_midp(0, 0)
        assert result == pytest.approx(0.0)

    def test_deg_midp_360_0(self):
        """Test midpoint of 360 and 0."""
        result = ephem.deg_midp(360, 0)
        assert result == pytest.approx(0.0)

    def test_deg_midp_45_135(self):
        """Test midpoint of 45 and 135."""
        result = ephem.deg_midp(45, 135)
        assert result == pytest.approx(90.0)

    def test_deg_midp_315_45(self):
        """Test midpoint of 315 and 45."""
        result = ephem.deg_midp(315, 45)
        assert result == pytest.approx(0.0)


class TestDegMidpSymmetry:
    """Tests for symmetry properties."""

    @pytest.mark.parametrize(
        "a,b",
        [
            (0, 90),
            (90, 180),
            (180, 270),
            (270, 360),
            (10, 350),
            (30, 60),
            (100, 200),
        ],
    )
    def test_deg_midp_order_invariance(self, a, b):
        """Test that deg_midp(a, b) == deg_midp(b, a)."""
        result_ab = ephem.deg_midp(a, b)
        result_ba = ephem.deg_midp(b, a)
        assert result_ab == pytest.approx(result_ba, abs=1e-10)


class TestDegMidpEdgeCases:
    """Edge case tests for deg_midp."""

    def test_deg_midp_very_small_difference(self):
        """Test with very small angle difference."""
        result = ephem.deg_midp(0.0001, 0.0002)
        assert result == pytest.approx(0.00015, abs=1e-10)

    def test_deg_midp_near_360_boundary(self):
        """Test near the 360 boundary."""
        result = ephem.deg_midp(359.9, 0.1)
        assert result == pytest.approx(0.0, abs=0.1)

    def test_deg_midp_fractional_angles(self):
        """Test with fractional angles."""
        result = ephem.deg_midp(30.5, 60.5)
        assert result == pytest.approx(45.5)
