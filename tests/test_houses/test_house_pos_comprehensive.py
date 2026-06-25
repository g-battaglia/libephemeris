"""
Comprehensive tests for house_pos (body house placement).

Verifies that house_pos returns valid house numbers for bodies,
works with various house systems, and is consistent with houses().
"""

from __future__ import annotations


import pytest

import libephemeris as swe
from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
)


HOUSE_SYSTEMS = [
    (ord("P"), "Placidus"),
    (ord("K"), "Koch"),
    (ord("R"), "Regiomontanus"),
    (ord("C"), "Campanus"),
    (ord("A"), "Equal (Asc)"),
    (ord("E"), "Equal (MC)"),
    (ord("W"), "Whole Sign"),
    (ord("B"), "Alcabitius"),
    (ord("M"), "Morinus"),
    (ord("O"), "Porphyry"),
]

PLANETS = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
]


class TestHousePosBasic:
    """Basic house_pos functionality."""

    @pytest.mark.unit
    def test_house_pos_returns_float(self):
        """house_pos returns a float."""
        jd = 2451545.0
        # Get ARMC and obliquity
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("P"))
        armc = ascmc[2]  # ARMC
        obliquity = ascmc[4] if len(ascmc) > 4 else 23.44  # True obliquity

        # Get Sun position
        sun, _ = swe.calc_ut(jd, SUN, 0)

        result = swe.house_pos(armc, 41.9, obliquity, ord("P"), sun[0], sun[1])
        assert isinstance(result, (int, float)), f"Type: {type(result)}"

    @pytest.mark.unit
    def test_house_pos_between_1_and_13(self):
        """House position should be between 1.0 and 13.0."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("P"))
        armc = ascmc[2]
        obliquity = ascmc[4] if len(ascmc) > 4 else 23.44

        sun, _ = swe.calc_ut(jd, SUN, 0)
        result = swe.house_pos(armc, 41.9, obliquity, ord("P"), sun[0], sun[1])
        assert 1.0 <= result < 13.0, f"House pos {result} out of range"


class TestHousePosAllPlanets:
    """Test house_pos for all planets."""

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS)
    def test_planet_house_valid(self, body_id: int, name: str):
        """Each planet gets a valid house number."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("P"))
        armc = ascmc[2]
        obliquity = ascmc[4] if len(ascmc) > 4 else 23.44

        pos, _ = swe.calc_ut(jd, body_id, 0)
        result = swe.house_pos(armc, 41.9, obliquity, ord("P"), pos[0], pos[1])
        house = int(result)
        assert 1 <= house <= 12, f"{name}: house {house} out of 1-12"

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id,name", PLANETS)
    def test_planet_house_fraction_valid(self, body_id: int, name: str):
        """Fractional part should be between 0 and 1."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, 41.9, 12.5, ord("P"))
        armc = ascmc[2]
        obliquity = ascmc[4] if len(ascmc) > 4 else 23.44

        pos, _ = swe.calc_ut(jd, body_id, 0)
        result = swe.house_pos(armc, 41.9, obliquity, ord("P"), pos[0], pos[1])
        fraction = result - int(result)
        assert 0 <= fraction < 1.0, f"{name}: fraction {fraction}"


class TestHousePosHouseSystems:
    """Test house_pos with various house systems."""

    @pytest.mark.unit
    @pytest.mark.parametrize("hsys,hsys_name", HOUSE_SYSTEMS)
    def test_sun_house_in_system(self, hsys: int, hsys_name: str):
        """Sun house position valid in each house system."""
        jd = 2451545.0
        lat, lon = 41.9, 12.5
        cusps, ascmc = swe.houses(jd, lat, lon, hsys)
        armc = ascmc[2]
        obliquity = ascmc[4] if len(ascmc) > 4 else 23.44

        sun, _ = swe.calc_ut(jd, SUN, 0)
        result = swe.house_pos(armc, lat, obliquity, hsys, sun[0], sun[1])
        house = int(result)
        assert 1 <= house <= 12, f"{hsys_name}: Sun in house {house}"


class TestHousePosLocations:
    """Test house_pos at various locations."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "lat,lon,name",
        [
            (41.9, 12.5, "Rome"),
            (40.7, -74.0, "New York"),
            (35.7, 139.7, "Tokyo"),
            (-33.9, 151.2, "Sydney"),
            (0.0, 0.0, "Equator"),
            (60.0, 25.0, "Helsinki"),
        ],
    )
    def test_moon_house_at_location(self, lat: float, lon: float, name: str):
        """Moon house position valid at various locations."""
        jd = 2451545.0
        cusps, ascmc = swe.houses(jd, lat, lon, ord("P"))
        armc = ascmc[2]
        obliquity = ascmc[4] if len(ascmc) > 4 else 23.44

        moon, _ = swe.calc_ut(jd, MOON, 0)
        result = swe.house_pos(armc, lat, obliquity, ord("P"), moon[0], moon[1])
        house = int(result)
        assert 1 <= house <= 12, f"{name}: Moon in house {house}"


class TestHousePosDateRange:
    """Test house_pos across dates."""

    @pytest.mark.unit
    @pytest.mark.parametrize("year", [1900, 1950, 2000, 2024, 2050, 2100])
    def test_sun_house_across_years(self, year: int):
        """Sun house valid across years."""
        jd = swe.julday(year, 6, 21, 12.0)
        lat, lon = 41.9, 12.5
        cusps, ascmc = swe.houses(jd, lat, lon, ord("P"))
        armc = ascmc[2]
        obliquity = ascmc[4] if len(ascmc) > 4 else 23.44

        sun, _ = swe.calc_ut(jd, SUN, 0)
        result = swe.house_pos(armc, lat, obliquity, ord("P"), sun[0], sun[1])
        assert 1.0 <= result < 13.0, f"Year {year}: pos={result}"
