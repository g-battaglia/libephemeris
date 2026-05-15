"""
Tests for planet center vs barycenter positions.

This module verifies that libephemeris uses planet centers (NAIF ID x99)
rather than system barycenters (NAIF ID x) for gas giant positions.

The difference between planet center and system barycenter can be significant:
- Jupiter: up to ~7000 km (due to Galilean moons)
- Saturn: up to ~3000 km (due to Titan and other moons)
- Uranus: up to ~500 km
- Neptune: up to ~300 km (due to Triton)
- Pluto: up to ~2000 km (due to Charon)

These differences cause position errors of several arcseconds when viewed from Earth.
"""

import pytest
import swisseph as swe
import libephemeris as ephem
from libephemeris.constants import (
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    MARS,
    FLG_SWIEPH,
    FLG_SPEED,
)


@pytest.mark.unit
class TestPlanetCenters:
    """Tests verifying planet centers are used instead of barycenters."""

    @pytest.fixture
    def gas_giants(self):
        """Gas giant planets that have significant barycenter offset."""
        return [
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ]

    def test_jupiter_position_matches_swisseph(self, standard_jd):
        """Test that Jupiter position matches SwissEph (which uses planet center)."""
        res_swe, _ = swe.calc_ut(standard_jd, JUPITER, FLG_SWIEPH)
        res_lib, _ = ephem.calc_ut(standard_jd, JUPITER, FLG_SWIEPH)

        lon_diff = abs(res_swe[0] - res_lib[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        # Position should match within 0.001 degrees (3.6 arcsec)
        # Previously with barycenter, this could differ by ~0.01 degrees
        assert lon_diff < 0.001, (
            f"Jupiter longitude diff: {lon_diff}° ({lon_diff * 3600:.1f} arcsec)"
        )

    def test_saturn_position_matches_swisseph(self, standard_jd):
        """Test that Saturn position matches SwissEph (which uses planet center)."""
        res_swe, _ = swe.calc_ut(standard_jd, SATURN, FLG_SWIEPH)
        res_lib, _ = ephem.calc_ut(standard_jd, SATURN, FLG_SWIEPH)

        lon_diff = abs(res_swe[0] - res_lib[0])
        if lon_diff > 180:
            lon_diff = 360 - lon_diff

        assert lon_diff < 0.001, (
            f"Saturn longitude diff: {lon_diff}° ({lon_diff * 3600:.1f} arcsec)"
        )

    def test_all_gas_giants_positions(self, standard_jd, gas_giants):
        """Test all gas giant positions against SwissEph."""
        for planet_id, planet_name in gas_giants:
            res_swe, _ = swe.calc_ut(standard_jd, planet_id, FLG_SWIEPH)
            res_lib, _ = ephem.calc_ut(standard_jd, planet_id, FLG_SWIEPH)

            lon_diff = abs(res_swe[0] - res_lib[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            assert lon_diff < 0.001, (
                f"{planet_name} longitude diff: {lon_diff}° ({lon_diff * 3600:.1f} arcsec)"
            )

    def test_gas_giants_with_speed(self, standard_jd, gas_giants):
        """Test gas giant positions and velocities with FLG_SPEED."""
        for planet_id, planet_name in gas_giants:
            res_swe, _ = swe.calc_ut(standard_jd, planet_id, FLG_SWIEPH | FLG_SPEED)
            res_lib, _ = ephem.calc_ut(
                standard_jd, planet_id, FLG_SWIEPH | FLG_SPEED
            )

            lon_diff = abs(res_swe[0] - res_lib[0])
            if lon_diff > 180:
                lon_diff = 360 - lon_diff

            speed_diff = abs(res_swe[3] - res_lib[3])

            assert lon_diff < 0.001, f"{planet_name} longitude diff: {lon_diff}°"
            assert speed_diff < 0.01, f"{planet_name} speed diff: {speed_diff}°/day"

    def test_gas_giants_over_date_range(self, gas_giants):
        """Test gas giant positions across multiple dates."""
        test_dates = [
            2451545.0,  # J2000
            2458849.5,  # 2020-01-01
            2460310.5,  # 2024-01-01
            2440587.5,  # 1970-01-01
        ]

        for jd in test_dates:
            for planet_id, planet_name in gas_giants:
                res_swe, _ = swe.calc_ut(jd, planet_id, FLG_SWIEPH)
                res_lib, _ = ephem.calc_ut(jd, planet_id, FLG_SWIEPH)

                lon_diff = abs(res_swe[0] - res_lib[0])
                if lon_diff > 180:
                    lon_diff = 360 - lon_diff

                assert lon_diff < 0.001, (
                    f'{planet_name} @ JD {jd}: diff {lon_diff}° ({lon_diff * 3600:.1f}")'
                )


@pytest.mark.unit
class TestPlanetMapContents:
    """Tests verifying the _PLANET_MAP dictionary contents."""

    def test_planet_map_uses_planet_centers(self):
        """Verify _PLANET_MAP doesn't use 'barycenter' for any planets."""
        from libephemeris.planets import _PLANET_MAP

        for planet_id, target_name in _PLANET_MAP.items():
            assert "barycenter" not in target_name.lower(), (
                f"Planet ID {planet_id} uses barycenter: '{target_name}'. "
                "Should use planet center instead."
            )

    def test_gas_giant_names_correct(self):
        """Verify gas giant target names are correct (without 'barycenter')."""
        from libephemeris.planets import _PLANET_MAP

        expected = {
            MARS: "mars",
            JUPITER: "jupiter",
            SATURN: "saturn",
            URANUS: "uranus",
            NEPTUNE: "neptune",
            PLUTO: "pluto",
        }

        for planet_id, expected_name in expected.items():
            actual_name = _PLANET_MAP.get(planet_id)
            assert actual_name == expected_name, (
                f"Planet ID {planet_id}: expected '{expected_name}', got '{actual_name}'"
            )
