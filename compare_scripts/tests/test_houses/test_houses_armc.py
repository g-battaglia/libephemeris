"""
Tests for houses_armc function.

Tests the house calculation from ARMC (Right Ascension of Medium Coeli)
instead of from Julian Day.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import *


class TestHousesArmcBasic:
    """Basic tests for houses_armc function."""

    @pytest.mark.unit
    def test_houses_armc_returns_correct_structure(self):
        """houses_armc() should return (cusps, ascmc) tuples."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # Should have 12 cusps
        assert len(cusps) == 12
        # Should have 8 ASCMC elements
        assert len(ascmc) == 8

    @pytest.mark.unit
    def test_houses_armc_cusps_are_valid_longitudes(self):
        """All cusps should be valid longitudes (0-360)."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        for i, cusp in enumerate(cusps):
            assert 0 <= cusp < 360, f"Cusp {i + 1} = {cusp} out of range"

    @pytest.mark.unit
    def test_houses_armc_asc_mc_valid(self):
        """ASC and MC should be valid longitudes."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        assert 0 <= ascmc[0] < 360, f"ASC = {ascmc[0]} out of range"
        assert 0 <= ascmc[1] < 360, f"MC = {ascmc[1]} out of range"

    @pytest.mark.unit
    def test_houses_armc_returns_correct_armc(self):
        """The returned ARMC in ascmc[2] should match input."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # ascmc[2] is ARMC
        assert abs(ascmc[2] - armc) < 0.001, f"ARMC mismatch: {ascmc[2]} vs {armc}"


class TestHousesArmcVsPyswisseph:
    """Compare houses_armc with pyswisseph."""

    @pytest.mark.comparison
    def test_houses_armc_matches_pyswisseph_placidus(self):
        """houses_armc should match pyswisseph for Placidus."""
        armc = 292.957072438026
        lat = 41.9
        eps = 23.43767671605485

        cusps_lib, ascmc_lib = ephem.houses_armc(armc, lat, eps, ord("P"))
        cusps_swe, ascmc_swe = swe.houses_armc(armc, lat, eps, b"P")

        # Compare ASC
        asc_diff = abs(ascmc_lib[0] - ascmc_swe[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        assert asc_diff < 0.01, f"Placidus ASC diff {asc_diff}"

        # Compare MC
        mc_diff = abs(ascmc_lib[1] - ascmc_swe[1])
        if mc_diff > 180:
            mc_diff = 360 - mc_diff
        assert mc_diff < 0.01, f"Placidus MC diff {mc_diff}"

        # Compare all cusps
        for i in range(12):
            diff = abs(cusps_lib[i] - cusps_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.1, f"Placidus cusp {i + 1} diff {diff}"

    @pytest.mark.comparison
    @pytest.mark.parametrize(
        "hsys,name,tolerance",
        [
            (ord("P"), "Placidus", 0.1),
            (ord("K"), "Koch", 0.1),
            (ord("O"), "Porphyry", 0.1),
            (ord("R"), "Regiomontanus", 0.1),
            (ord("C"), "Campanus", 0.1),
            (ord("E"), "Equal", 0.1),
            (ord("W"), "Whole Sign", 0.1),
            (ord("B"), "Alcabitius", 0.1),
            (ord("M"), "Morinus", 0.1),
        ],
    )
    def test_houses_armc_various_systems(self, hsys, name, tolerance):
        """Test houses_armc for various house systems."""
        armc = 150.0
        lat = 45.0
        eps = 23.44

        cusps_lib, ascmc_lib = ephem.houses_armc(armc, lat, eps, hsys)
        cusps_swe, ascmc_swe = swe.houses_armc(armc, lat, eps, bytes([hsys]))

        # Compare ASC
        asc_diff = abs(ascmc_lib[0] - ascmc_swe[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        assert asc_diff < tolerance, f"{name} ASC diff {asc_diff}"

        # Compare MC
        mc_diff = abs(ascmc_lib[1] - ascmc_swe[1])
        if mc_diff > 180:
            mc_diff = 360 - mc_diff
        assert mc_diff < tolerance, f"{name} MC diff {mc_diff}"


class TestHousesArmcConsistencyWithHouses:
    """Test that houses_armc is consistent with houses()."""

    @pytest.mark.unit
    def test_houses_armc_matches_houses_with_correct_obliquity(self):
        """houses_armc should produce same results as houses() with same inputs."""
        jd = 2451545.0
        lat = 41.9
        lon = 12.5

        # First, calculate using houses() to get ARMC and obliquity
        cusps_houses, ascmc_houses = ephem.houses(jd, lat, lon, ord("P"))

        armc = ascmc_houses[2]  # Get ARMC from houses()

        # Get the obliquity from pyswisseph for consistency
        eps = swe.calc_ut(jd, swe.ECL_NUT)[0][0]

        # Now calculate using houses_armc with the same ARMC and obliquity
        cusps_armc, ascmc_armc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # Compare ASC (might have small differences due to obliquity calc)
        asc_diff = abs(ascmc_houses[0] - ascmc_armc[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        # Allow for some difference due to obliquity calculation differences
        assert asc_diff < 0.01, f"ASC diff {asc_diff}"

        # Compare MC
        mc_diff = abs(ascmc_houses[1] - ascmc_armc[1])
        if mc_diff > 180:
            mc_diff = 360 - mc_diff
        assert mc_diff < 0.01, f"MC diff {mc_diff}"


class TestHousesArmcEdgeCases:
    """Test edge cases for houses_armc."""

    @pytest.mark.edge_case
    def test_houses_armc_at_equator(self):
        """houses_armc should work at equator (lat=0)."""
        armc = 180.0
        lat = 0.0
        eps = 23.44

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        assert 0 <= ascmc[0] < 360
        assert 0 <= ascmc[1] < 360

    @pytest.mark.edge_case
    @pytest.mark.parametrize("lat", [60.0, 65.0])
    def test_houses_armc_high_latitude(self, lat):
        """houses_armc should handle high latitudes below polar circle."""
        armc = 100.0
        eps = 23.44

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # Should return valid values (even if approximation)
        assert 0 <= ascmc[0] < 360
        assert 0 <= ascmc[1] < 360

    @pytest.mark.edge_case
    def test_houses_armc_polar_circle_raises_error(self):
        """houses_armc should raise Error for Placidus at polar latitudes."""
        armc = 100.0
        eps = 23.44

        # 70° is within the polar circle - Placidus should raise an error
        with pytest.raises(ephem.Error) as exc_info:
            ephem.houses_armc(armc, 70.0, eps, ord("P"))

        assert "polar circle" in str(exc_info.value).lower()

        # But Porphyry should work
        cusps, ascmc = ephem.houses_armc(armc, 70.0, eps, ord("O"))
        assert 0 <= ascmc[0] < 360

    @pytest.mark.edge_case
    @pytest.mark.parametrize("armc", [0.0, 90.0, 180.0, 270.0, 359.99])
    def test_houses_armc_various_armc_values(self, armc):
        """houses_armc should handle various ARMC values."""
        lat = 45.0
        eps = 23.44

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        assert 0 <= ascmc[0] < 360
        assert 0 <= ascmc[1] < 360
        assert len(cusps) == 12

    @pytest.mark.edge_case
    def test_houses_armc_southern_hemisphere(self):
        """houses_armc should work in southern hemisphere."""
        armc = 200.0
        lat = -33.9  # Sydney
        eps = 23.44

        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        assert 0 <= ascmc[0] < 360
        assert 0 <= ascmc[1] < 360

    @pytest.mark.edge_case
    def test_houses_armc_armc_normalization(self):
        """ARMC values outside 0-360 should be normalized."""
        lat = 45.0
        eps = 23.44

        # Test with ARMC > 360
        cusps1, ascmc1 = ephem.houses_armc(400.0, lat, eps, ord("P"))
        cusps2, ascmc2 = ephem.houses_armc(40.0, lat, eps, ord("P"))

        assert abs(ascmc1[0] - ascmc2[0]) < 0.001
        assert abs(ascmc1[1] - ascmc2[1]) < 0.001

        # Test with negative ARMC
        cusps3, ascmc3 = ephem.houses_armc(-40.0, lat, eps, ord("P"))
        cusps4, ascmc4 = ephem.houses_armc(320.0, lat, eps, ord("P"))

        assert abs(ascmc3[0] - ascmc4[0]) < 0.001
        assert abs(ascmc3[1] - ascmc4[1]) < 0.001


class TestHousesArmcAlias:
    """Test the houses_armc alias (without swe_ prefix)."""

    @pytest.mark.unit
    def test_houses_armc_alias_exists(self):
        """houses_armc alias should exist."""
        assert hasattr(ephem, "houses_armc")

    @pytest.mark.unit
    def test_houses_armc_alias_works(self):
        """houses_armc alias should work correctly."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # Both should return the same result
        cusps1, ascmc1 = ephem.houses_armc(armc, lat, eps, ord("P"))
        cusps2, ascmc2 = ephem.houses_armc(armc, lat, eps, ord("P"))

        assert cusps1 == cusps2
        assert ascmc1 == ascmc2


class TestHousesArmcEx2Basic:
    """Basic tests for houses_armc_ex2 function."""

    @pytest.mark.unit
    def test_houses_armc_ex2_returns_correct_structure(self):
        """houses_armc_ex2() should return (cusps, ascmc, cusps_speed, ascmc_speed) tuples."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P")
        )

        # Should have 12 cusps
        assert len(cusps) == 12
        # Should have 8 ASCMC elements
        assert len(ascmc) == 8
        # Should have 12 cusp velocities
        assert len(cusps_speed) == 12
        # Should have 8 ASCMC velocities
        assert len(ascmc_speed) == 8

    @pytest.mark.unit
    def test_houses_armc_ex2_cusps_match_houses_armc(self):
        """Cusps and ascmc from houses_armc_ex2 should match houses_armc."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        cusps_ex2, ascmc_ex2, _, _ = ephem.houses_armc_ex2(armc, lat, eps, ord("P"))
        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))

        # Positions should be identical
        assert cusps_ex2 == cusps
        assert ascmc_ex2 == ascmc

    @pytest.mark.unit
    def test_houses_armc_ex2_velocities_are_reasonable(self):
        """Velocities should be in reasonable range for house cusps when FLG_SPEED is set."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        _, _, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P"), FLG_SPEED
        )

        # House cusps and angles rotate with Earth: ~360°/day
        # Velocities should be approximately 360°/day (±50% due to projection effects)
        for i, speed in enumerate(cusps_speed):
            assert 100 < abs(speed) < 600, f"Cusp {i + 1} speed {speed} out of range"

        # ASC velocity
        assert 100 < abs(ascmc_speed[0]) < 600, (
            f"ASC speed {ascmc_speed[0]} out of range"
        )
        # MC velocity
        assert 100 < abs(ascmc_speed[1]) < 600, (
            f"MC speed {ascmc_speed[1]} out of range"
        )

    @pytest.mark.unit
    def test_houses_armc_ex2_velocities_are_positive(self):
        """Velocities should generally be positive (prograde motion) when FLG_SPEED is set."""
        armc = 150.0
        lat = 45.0
        eps = 23.44

        _, _, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P"), FLG_SPEED
        )

        # Most cusps should have positive velocity (increasing longitude)
        positive_cusps = sum(1 for s in cusps_speed if s > 0)
        assert positive_cusps >= 8, (
            f"Only {positive_cusps} cusps have positive velocity"
        )

        # ARMC velocity (ascmc_speed[2]) should be the sidereal rotation
        # rate (~360.9856°/day), not the solar rate (360°/day).
        assert abs(ascmc_speed[2] - 360.9856) < 0.01, (
            f"ARMC speed {ascmc_speed[2]} != ~360.986"
        )


class TestHousesArmcEx2VariousSystems:
    """Test houses_armc_ex2 for various house systems."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys,name",
        [
            (ord("P"), "Placidus"),
            (ord("K"), "Koch"),
            (ord("O"), "Porphyry"),
            (ord("R"), "Regiomontanus"),
            (ord("C"), "Campanus"),
            (ord("E"), "Equal"),
            (ord("W"), "Whole Sign"),
            (ord("B"), "Alcabitius"),
            (ord("M"), "Morinus"),
        ],
    )
    def test_houses_armc_ex2_various_systems(self, hsys, name):
        """Test houses_armc_ex2 works for various house systems."""
        armc = 150.0
        lat = 45.0
        eps = 23.44

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, hsys
        )

        # All should return valid structure
        assert len(cusps) == 12, f"{name} cusps count wrong"
        assert len(ascmc) == 8, f"{name} ascmc count wrong"
        assert len(cusps_speed) == 12, f"{name} cusps_speed count wrong"
        assert len(ascmc_speed) == 8, f"{name} ascmc_speed count wrong"

        # All positions should be valid longitudes
        for i, cusp in enumerate(cusps):
            assert 0 <= cusp < 360, f"{name} cusp {i + 1} = {cusp} out of range"


class TestHousesArmcEx2EdgeCases:
    """Test edge cases for houses_armc_ex2."""

    @pytest.mark.edge_case
    def test_houses_armc_ex2_at_equator(self):
        """houses_armc_ex2 should work at equator (lat=0)."""
        armc = 180.0
        lat = 0.0
        eps = 23.44

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P")
        )

        assert len(cusps) == 12
        assert len(cusps_speed) == 12

    @pytest.mark.edge_case
    @pytest.mark.parametrize("lat", [60.0, 65.0])
    def test_houses_armc_ex2_high_latitude(self, lat):
        """houses_armc_ex2 should handle high latitudes."""
        armc = 100.0
        eps = 23.44

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P")
        )

        assert len(cusps) == 12
        assert len(cusps_speed) == 12

    @pytest.mark.edge_case
    def test_houses_armc_ex2_southern_hemisphere(self):
        """houses_armc_ex2 should work in southern hemisphere."""
        armc = 200.0
        lat = -33.9  # Sydney
        eps = 23.44

        cusps, ascmc, cusps_speed, ascmc_speed = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P")
        )

        assert len(cusps) == 12
        assert len(cusps_speed) == 12


class TestHousesArmcEx2Alias:
    """Test the houses_armc_ex2 alias (without swe_ prefix)."""

    @pytest.mark.unit
    def test_houses_armc_ex2_alias_exists(self):
        """houses_armc_ex2 alias should exist."""
        assert hasattr(ephem, "houses_armc_ex2")

    @pytest.mark.unit
    def test_houses_armc_ex2_alias_works(self):
        """houses_armc_ex2 alias should work correctly."""
        armc = 292.957
        lat = 41.9
        eps = 23.4393

        # Both should return the same result
        result1 = ephem.houses_armc_ex2(armc, lat, eps, ord("P"))
        result2 = ephem.houses_armc_ex2(armc, lat, eps, ord("P"))

        assert result1 == result2


class TestHousesArmcEx2VsPyswisseph:
    """Compare houses_armc_ex2 with pyswisseph houses_armc_ex2."""

    @pytest.mark.comparison
    def test_houses_armc_ex2_cusps_match_pyswisseph(self):
        """houses_armc_ex2 cusps should match pyswisseph."""
        armc = 292.957072438026
        lat = 41.9
        eps = 23.43767671605485

        cusps_lib, ascmc_lib, cusps_speed_lib, ascmc_speed_lib = (
            ephem.houses_armc_ex2(armc, lat, eps, ord("P"))
        )
        cusps_swe, ascmc_swe, cusps_speed_swe, ascmc_speed_swe = swe.houses_armc_ex2(
            armc, lat, eps, b"P"
        )

        # Compare cusps
        for i in range(12):
            diff = abs(cusps_lib[i] - cusps_swe[i])
            if diff > 180:
                diff = 360 - diff
            assert diff < 0.1, f"Cusp {i + 1} diff {diff}"

        # Compare ASC
        asc_diff = abs(ascmc_lib[0] - ascmc_swe[0])
        if asc_diff > 180:
            asc_diff = 360 - asc_diff
        assert asc_diff < 0.01, f"ASC diff {asc_diff}"

        # Compare MC
        mc_diff = abs(ascmc_lib[1] - ascmc_swe[1])
        if mc_diff > 180:
            mc_diff = 360 - mc_diff
        assert mc_diff < 0.01, f"MC diff {mc_diff}"

    @pytest.mark.comparison
    def test_houses_armc_ex2_velocities_match_pyswisseph(self):
        """houses_armc_ex2 velocities should be close to pyswisseph when FLG_SPEED is set."""
        armc = 150.0
        lat = 45.0
        eps = 23.44

        _, _, cusps_speed_lib, ascmc_speed_lib = ephem.houses_armc_ex2(
            armc, lat, eps, ord("P"), FLG_SPEED
        )
        _, _, cusps_speed_swe, ascmc_speed_swe = swe.houses_armc_ex2(
            armc, lat, eps, b"P"
        )

        # Compare cusp velocities (allow 6% tolerance due to different finite diff methods)
        for i in range(12):
            if abs(cusps_speed_swe[i]) > 1:  # Avoid division by very small values
                rel_diff = abs(cusps_speed_lib[i] - cusps_speed_swe[i]) / abs(
                    cusps_speed_swe[i]
                )
                assert rel_diff < 0.06, (
                    f"Cusp {i + 1} speed rel diff {rel_diff * 100:.1f}%"
                )

        # Compare ASC velocity
        if abs(ascmc_speed_swe[0]) > 1:
            rel_diff = abs(ascmc_speed_lib[0] - ascmc_speed_swe[0]) / abs(
                ascmc_speed_swe[0]
            )
            assert rel_diff < 0.06, f"ASC speed rel diff {rel_diff * 100:.1f}%"

        # Compare MC velocity
        if abs(ascmc_speed_swe[1]) > 1:
            rel_diff = abs(ascmc_speed_lib[1] - ascmc_speed_swe[1]) / abs(
                ascmc_speed_swe[1]
            )
            assert rel_diff < 0.06, f"MC speed rel diff {rel_diff * 100:.1f}%"
