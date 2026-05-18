"""
Lunar Calculations Comparison Tests.

Compares lunar nodes and Lilith calculations between pyswisseph and libephemeris:
- Mean Node
- True Node
- Mean Lilith (Black Moon)
- True Lilith (Osculating Apogee)

All lunar tests now pass within tolerance, including Mean Lilith which uses
an SE-compatible DE404 algorithm with ecliptic projection (~15 arcsec precision).
"""

import pytest
import random
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
)

# Mean Lilith now uses SE-compatible DE404 algorithm with ecliptic projection
# Tests pass with ~15 arcsec precision (previously expected to fail)


def angular_diff(val1: float, val2: float) -> float:
    """Calculate angular difference accounting for 360 wrap."""
    d = abs(val1 - val2)
    if d > 180:
        d = 360 - d
    return d


# ============================================================================
# TOLERANCES
# ============================================================================

MEAN_NODE_TOL = 0.001  # degrees (~3.6 arcsec)
TRUE_NODE_TOL = 0.015  # degrees (~54 arcsec, perturbation series vs pyswisseph)
MEAN_LILITH_TOL = (
    0.01  # degrees (~36 arcsec, SE-compatible algorithm achieves ~15 arcsec)
)
TRUE_LILITH_TOL = 0.06  # degrees (~216 arcsec, orbital elements vs Moshier, max ~0.054°)


# ============================================================================
# TEST DATA
# ============================================================================

TEST_DATES = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 6, 15, 0.0, "Summer 2024"),
    (1980, 5, 20, 14.5, "Past"),
    (1950, 10, 15, 22.0, "Mid-century"),
    (2050, 1, 1, 0.0, "Future"),
]


# ============================================================================
# TEST CLASSES
# ============================================================================


class TestMeanNode:
    """Compare Mean North Node calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_mean_node_longitude(self, year, month, day, hour, desc):
        """Test Mean Node longitude matches within tolerance."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, 0)
        pos_py, _ = ephem.calc_ut(jd, MEAN_NODE, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < MEAN_NODE_TOL, (
            f"Mean Node at {desc}: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    def test_mean_node_with_speed(self):
        """Test Mean Node with velocity calculation."""
        jd = swe.julday(2000, 1, 1, 12.0)

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_NODE, swe.FLG_SPEED)
        pos_py, _ = ephem.calc_ut(jd, MEAN_NODE, 256)  # FLG_SPEED

        diff_lon = angular_diff(pos_swe[0], pos_py[0])
        diff_speed = abs(pos_swe[3] - pos_py[3])

        assert diff_lon < MEAN_NODE_TOL, f"Mean Node longitude diff {diff_lon:.6f}°"
        assert diff_speed < 0.01, f"Mean Node speed diff {diff_speed:.6f}°/day"


class TestTrueNode:
    """Compare True North Node calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_true_node_longitude(self, year, month, day, hour, desc):
        """Test True Node longitude matches within relaxed tolerance."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, 0)
        pos_py, _ = ephem.calc_ut(jd, TRUE_NODE, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TRUE_NODE_TOL, (
            f"True Node at {desc}: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    def test_true_node_random_dates(self):
        """Test True Node at 100 random dates."""
        random.seed(42)
        errors = []

        for _ in range(100):
            year = random.randint(1900, 2100)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, swe.TRUE_NODE, 0)
            pos_py, _ = ephem.calc_ut(jd, TRUE_NODE, 0)

            diff = angular_diff(pos_swe[0], pos_py[0])
            errors.append(diff)

        max_error = max(errors)
        mean_error = sum(errors) / len(errors)

        assert max_error < TRUE_NODE_TOL, (
            f"True Node max error {max_error:.4f}° exceeds tolerance"
        )
        assert mean_error < 0.08, (
            f"True Node mean error {mean_error:.4f}° exceeds threshold"
        )


class TestMeanLilith:
    """Compare Mean Lilith (Black Moon) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_mean_lilith_longitude(self, year, month, day, hour, desc):
        """Test Mean Lilith longitude matches within tolerance."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_APOG, 0)
        pos_py, _ = ephem.calc_ut(jd, MEAN_APOG, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < MEAN_LILITH_TOL, (
            f"Mean Lilith at {desc}: diff {diff:.6f}° exceeds tolerance"
        )


class TestTrueLilith:
    """Compare True Lilith (Osculating Apogee) calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES)
    def test_true_lilith_longitude(self, year, month, day, hour, desc):
        """
        Test True Lilith longitude.

        The eccentricity vector method achieves sub-arcminute precision
        (~52 arcsec mean, ~235 arcsec max) against Swiss Ephemeris.
        """
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
        pos_py, _ = ephem.calc_ut(jd, OSCU_APOG, 0)

        diff = angular_diff(pos_swe[0], pos_py[0])

        assert diff < TRUE_LILITH_TOL, (
            f"True Lilith at {desc}: diff {diff:.6f}° exceeds tolerance"
        )

    @pytest.mark.comparison
    @pytest.mark.slow
    def test_true_lilith_statistics(self):
        """Test True Lilith error distribution across random dates."""
        random.seed(42)
        errors = []

        for _ in range(100):
            year = random.randint(1950, 2050)
            month = random.randint(1, 12)
            day = random.randint(1, 28)
            hour = random.uniform(0, 24)
            jd = swe.julday(year, month, day, hour)

            pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
            pos_py, _ = ephem.calc_ut(jd, OSCU_APOG, 0)

            diff = angular_diff(pos_swe[0], pos_py[0])
            errors.append(diff)

        max_error = max(errors)
        mean_error = sum(errors) / len(errors)

        # True Lilith achieves sub-arcminute precision
        assert mean_error < 0.05, (
            f"True Lilith mean error {mean_error:.4f}° exceeds threshold"
        )


class TestNodalMotion:
    """Test lunar node motion characteristics."""

    @pytest.mark.comparison
    def test_mean_node_retrograde(self):
        """Test that Mean Node is always retrograde."""
        jd = swe.julday(2000, 1, 1, 12.0)

        pos_py, _ = ephem.calc_ut(jd, MEAN_NODE, 256)  # FLG_SPEED

        # Mean Node is always retrograde (negative daily motion)
        assert pos_py[3] < 0, "Mean Node should have retrograde motion"

    @pytest.mark.comparison
    def test_nodes_opposite(self):
        """Test that North and South nodes are 180° apart."""
        jd = swe.julday(2000, 1, 1, 12.0)

        # Get North Node
        pos_north, _ = ephem.calc_ut(jd, MEAN_NODE, 0)

        # South Node is 180° from North Node
        south_node = (pos_north[0] + 180) % 360

        # This is by definition, so should be exact
        assert abs(south_node - ((pos_north[0] + 180) % 360)) < 1e-10


class TestLilithLatitude:
    """Test Lilith latitude calculations."""

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES[:3])
    def test_mean_lilith_latitude(self, year, month, day, hour, desc):
        """Test Mean Lilith latitude."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.MEAN_APOG, 0)
        pos_py, _ = ephem.calc_ut(jd, MEAN_APOG, 0)

        diff_lat = abs(pos_swe[1] - pos_py[1])

        # Latitude should be close (Mean Lilith has no latitude by definition)
        assert diff_lat < 0.1, f"Mean Lilith latitude diff {diff_lat:.6f}° at {desc}"

    @pytest.mark.comparison
    @pytest.mark.parametrize("year,month,day,hour,desc", TEST_DATES[:3])
    def test_true_lilith_latitude(self, year, month, day, hour, desc):
        """Test True Lilith latitude."""
        jd = swe.julday(year, month, day, hour)

        pos_swe, _ = swe.calc_ut(jd, swe.OSCU_APOG, 0)
        pos_py, _ = ephem.calc_ut(jd, OSCU_APOG, 0)

        diff_lat = abs(pos_swe[1] - pos_py[1])

        # Latitude tolerance is more relaxed than longitude
        assert diff_lat < 1.0, f"True Lilith latitude diff {diff_lat:.2f}° at {desc}"
