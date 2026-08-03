"""
Tests for house_name and get_planet_name functions.

Verifies correct name lookup for all house systems and body IDs.
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
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    OSCU_APOG,
    EARTH,
    CHIRON,
    CERES,
    PALLAS,
    JUNO,
    VESTA,
)


class TestHouseName:
    """Test house_name function."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys_char,expected",
        [
            ("P", "Placidus"),
            ("K", "Koch"),
            ("R", "Regiomontanus"),
            ("C", "Campanus"),
            ("B", "Alcabitius"),
            ("O", "Porphyry"),
            ("M", "Morinus"),
            ("W", None),  # Whole Sign or Equal
        ],
    )
    def test_known_systems_return_names(self, hsys_char: str, expected):
        """Known house system chars return correct names."""
        result = swe.house_name(ord(hsys_char))
        assert isinstance(result, str)
        assert len(result) > 0
        if expected:
            assert expected in result, (
                f"'{hsys_char}': got '{result}', expected '{expected}'"
            )

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "hsys_char",
        list("PKROCABMWETUGHFVXSLQNYDI"),
    )
    def test_all_systems_return_nonempty(self, hsys_char: str):
        """All supported system chars return non-empty strings."""
        result = swe.house_name(ord(hsys_char))
        assert isinstance(result, str)
        assert len(result) > 0
        assert result != "Unknown", f"System '{hsys_char}' returned 'Unknown'"

    @pytest.mark.unit
    def test_savard_a_system(self):
        """Savard-A house system ('J') should return its name."""
        result = swe.house_name(ord("J"))
        assert isinstance(result, str)
        assert len(result) > 0
        assert result != "Unknown"

    @pytest.mark.unit
    def test_unknown_system_returns_empty_string(self):
        """An unknown system selector returns the empty string.

        Measured contract (matched exactly): every unrecognized selector —
        letters outside the table and non-letters alike — maps to '', not to
        a placeholder such as "Unknown".
        """
        for selector in (ord("Z"), ord("z"), ord("5")):
            result = swe.house_name(selector)
            assert isinstance(result, str)
            assert result == ""

    @pytest.mark.unit
    def test_alias_works(self):
        """house_name alias should work the same."""
        r1 = swe.house_name(ord("P"))
        r2 = swe.house_name(ord("P"))
        assert r1 == r2


class TestGetPlanetName:
    """Test get_planet_name function."""

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,expected_substr",
        [
            (SUN, "Sun"),
            (MOON, "Moon"),
            (MERCURY, "Mercury"),
            (VENUS, "Venus"),
            (MARS, "Mars"),
            (JUPITER, "Jupiter"),
            (SATURN, "Saturn"),
            (URANUS, "Uranus"),
            (NEPTUNE, "Neptune"),
            (PLUTO, "Pluto"),
        ],
    )
    def test_major_planets(self, body_id: int, expected_substr: str):
        """Major planets return correct names."""
        result = swe.get_planet_name(body_id)
        assert isinstance(result, str)
        assert expected_substr in result, f"Body {body_id}: got '{result}'"

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,expected_substr",
        [
            (MEAN_NODE, "Node"),
            (TRUE_NODE, "Node"),
            (MEAN_APOG, "Apog"),
            (OSCU_APOG, "Apog"),
        ],
    )
    def test_lunar_nodes_apsides(self, body_id: int, expected_substr: str):
        """Lunar nodes and apsides return appropriate names."""
        result = swe.get_planet_name(body_id)
        assert isinstance(result, str)
        assert expected_substr.lower() in result.lower(), (
            f"Body {body_id}: got '{result}', expected substring '{expected_substr}'"
        )

    @pytest.mark.unit
    def test_earth_name(self):
        """Earth body ID returns 'Earth'."""
        result = swe.get_planet_name(EARTH)
        assert "Earth" in result

    @pytest.mark.unit
    @pytest.mark.parametrize(
        "body_id,expected_substr",
        [
            (CHIRON, "Chiron"),
            (CERES, "Ceres"),
            (PALLAS, "Pallas"),
            (JUNO, "Juno"),
            (VESTA, "Vesta"),
        ],
    )
    def test_asteroids(self, body_id: int, expected_substr: str):
        """Asteroid names are correct."""
        result = swe.get_planet_name(body_id)
        assert isinstance(result, str)
        assert expected_substr in result, f"Body {body_id}: got '{result}'"

    @pytest.mark.unit
    def test_unknown_body_returns_string(self):
        """Unknown body IDs return the compatibility API's empty string."""
        result = swe.get_planet_name(8888)
        assert isinstance(result, str)
        assert result == ""

    @pytest.mark.unit
    def test_unknown_asteroid_returns_empty(self):
        """Unknown AST_OFFSET ids return '' (reference parity)."""
        # swe.get_planet_name(99999) == '' in the reference ephemeris (asteroid 89999,
        # no name without its ephemeris file)
        assert swe.get_planet_name(99999) == ""

    @pytest.mark.unit
    def test_get_planet_name_alias(self):
        """get_planet_name alias should work the same."""
        r1 = swe.get_planet_name(MARS)
        r2 = swe.get_planet_name(MARS)
        assert r1 == r2

    @pytest.mark.unit
    @pytest.mark.parametrize("body_id", range(40, 49))
    def test_uranian_bodies_named(self, body_id: int):
        """Uranian body IDs 40-48 should return names."""
        result = swe.get_planet_name(body_id)
        assert isinstance(result, str)
        assert len(result) > 0
        # Should not be just "Unknown (40)" etc
        assert "Unknown" not in result or body_id == 48, (
            f"Body {body_id}: got '{result}'"
        )
