"""
Unit tests for lunar module: Mean/True Nodes and Lilith.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


@pytest.mark.unit
class TestLunarNodes:
    """Tests for Mean and True Lunar Nodes."""

    def test_mean_node_j2000(self, standard_jd):
        """Test Mean Node at J2000."""
        mean_node, _ = ephem.swe_calc_ut(standard_jd, SE_MEAN_NODE, 0)
        assert 124 < mean_node[0] < 126  # Expected ~125° at J2000
        assert mean_node[1] == 0.0  # Latitude always 0

    def test_mean_node_vs_swisseph(self, standard_jd):
        """Compare Mean Node with SwissEph."""
        node_py, _ = ephem.swe_calc_ut(standard_jd, SE_MEAN_NODE, 0)
        node_swe, _ = swe.calc_ut(standard_jd, swe.MEAN_NODE, 0)

        diff = abs(node_py[0] - node_swe[0])
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.01, f"Mean Node diff: {diff}°"

    def test_true_node_j2000(self, standard_jd):
        """Test True Node at J2000."""
        true_node, _ = ephem.swe_calc_ut(standard_jd, SE_TRUE_NODE, 0)
        assert 0 <= true_node[0] < 360

    def test_true_node_vs_swisseph(self, standard_jd):
        """Compare True Node with SwissEph."""
        node_py, _ = ephem.swe_calc_ut(standard_jd, SE_TRUE_NODE, 0)
        node_swe, _ = swe.calc_ut(standard_jd, swe.TRUE_NODE, 0)

        diff = abs(node_py[0] - node_swe[0])
        # Proper wrap-around handling
        if diff > 180:
            diff = 360 - diff

        # True node calculation may have small differences
        assert diff < 2.0, f"True Node diff: {diff}°"

    def test_mean_vs_true_node(self, standard_jd):
        """Mean and True nodes should be different."""
        mean_node, _ = ephem.swe_calc_ut(standard_jd, SE_MEAN_NODE, 0)
        true_node, _ = ephem.swe_calc_ut(standard_jd, SE_TRUE_NODE, 0)

        diff = abs(mean_node[0] - true_node[0])
        assert diff > 0.1, "Mean and True should differ"

    def test_south_node(self, standard_jd):
        """Test South Node is 180° from North Node."""
        # Note: SE_MEAN_SOUTH_NODE would be -SE_MEAN_NODE
        # For now we test the logic is sound
        mean_north, _ = ephem.swe_calc_ut(standard_jd, SE_MEAN_NODE, 0)
        expected_south = (mean_north[0] + 180.0) % 360.0

        # South node calculation
        south_id = -SE_MEAN_NODE
        mean_south, _ = ephem.swe_calc_ut(standard_jd, south_id, 0)

        diff = abs(mean_south[0] - expected_south)
        if diff > 180:
            diff = 360 - diff

        assert diff < 0.001, f"South node should be 180° from North, diff={diff}°"


@pytest.mark.unit
class TestLilith:
    """Tests for Mean and True Lilith (Lunar Apogee)."""

    def test_mean_lilith_j2000(self, standard_jd):
        """Test Mean Lilith at J2000."""
        lilith, _ = ephem.swe_calc_ut(standard_jd, SE_MEAN_APOG, 0)
        assert 0 <= lilith[0] < 360
        assert lilith[1] == 0.0  # Mean Lilith has 0 latitude

    def test_mean_lilith_vs_swisseph(self, standard_jd):
        """Compare Mean Lilith with SwissEph."""
        lilith_py, _ = ephem.swe_calc_ut(standard_jd, SE_MEAN_APOG, 0)
        lilith_swe, _ = swe.calc_ut(standard_jd, swe.MEAN_APOG, 0)

        diff = abs(lilith_py[0] - lilith_swe[0])
        # Proper wrap-around
        if diff > 180:
            diff = 360 - diff

        assert diff < 5.0, f"Mean Lilith diff: {diff}°"

    def test_true_lilith_j2000(self, standard_jd):
        """Test True Lilith at J2000."""
        lilith, _ = ephem.swe_calc_ut(standard_jd, SE_OSCU_APOG, 0)
        assert 0 <= lilith[0] < 360

    def test_true_lilith_vs_swisseph(self, standard_jd):
        """Compare True Lilith with SwissEph."""
        lilith_py, _ = ephem.swe_calc_ut(standard_jd, SE_OSCU_APOG, 0)
        lilith_swe, _ = swe.calc_ut(standard_jd, swe.OSCU_APOG, 0)

        diff = abs(lilith_py[0] - lilith_swe[0])
        # Proper wrap-around
        if diff > 180:
            diff = 360 - diff

        # Osculating apogee may have larger differences
        assert diff < 10.0, f"True Lilith diff: {diff}°"

    def test_mean_vs_true_lilith(self, standard_jd):
        """Mean and True Lilith should be different."""
        mean_lilith, _ = ephem.swe_calc_ut(standard_jd, SE_MEAN_APOG, 0)
        true_lilith, _ = ephem.swe_calc_ut(standard_jd, SE_OSCU_APOG, 0)

        diff = abs(mean_lilith[0] - true_lilith[0])
        assert diff > 0.1, "Mean and True Lilith should differ"

    @pytest.mark.parametrize(
        "year,month,day",
        [
            (2000, 1, 1),
            (1980, 6, 15),
            (2024, 11, 28),
        ],
    )
    def test_lilith_multiple_dates(self, year, month, day):
        """Test Lilith calculation at various dates."""
        jd = ephem.swe_julday(year, month, day, 12.0)

        mean_lilith, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, 0)
        true_lilith, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, 0)

        # Both should return valid longitudes
        assert 0 <= mean_lilith[0] < 360
        assert 0 <= true_lilith[0] < 360


@pytest.mark.integration
class TestLunarIntegration:
    """Integration tests for lunar calculations across date ranges."""

    def test_nodes_over_time(self, test_dates):
        """Test nodes across multiple dates."""
        for year, month, day, hour, name in test_dates:
            jd = ephem.swe_julday(year, month, day, hour)

            mean_node, _ = ephem.swe_calc_ut(jd, SE_MEAN_NODE, 0)
            true_node, _ = ephem.swe_calc_ut(jd, SE_TRUE_NODE, 0)

            # Should always be valid positions
            assert 0 <= mean_node[0] < 360, f"Invalid Mean Node for {name}"
            assert 0 <= true_node[0] < 360, f"Invalid True Node for {name}"

    def test_lilith_over_time(self, test_dates):
        """Test Lilith across multiple dates."""
        for year, month, day, hour, name in test_dates:
            jd = ephem.swe_julday(year, month, day, hour)

            mean_lilith, _ = ephem.swe_calc_ut(jd, SE_MEAN_APOG, 0)
            true_lilith, _ = ephem.swe_calc_ut(jd, SE_OSCU_APOG, 0)

            # Should always be valid positions
            assert 0 <= mean_lilith[0] < 360, f"Invalid Mean Lilith for {name}"
            assert 0 <= true_lilith[0] < 360, f"Invalid True Lilith for {name}"
