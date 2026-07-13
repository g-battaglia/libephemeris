"""
Tests for get_planet_name function.
"""

from libephemeris import (
    get_planet_name,
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
    VULCAN,
    WHITE_MOON,
    PROSERPINA,
    WALDEMATH,
)


class TestGetPlanetName:
    """Tests for the get_planet_name function."""

    def test_sun(self):
        """Test Sun name lookup."""
        assert get_planet_name(SUN) == "Sun"
        assert get_planet_name(0) == "Sun"

    def test_moon(self):
        """Test Moon name lookup."""
        assert get_planet_name(MOON) == "Moon"
        assert get_planet_name(1) == "Moon"

    def test_planets(self):
        """Test all main planet name lookups."""
        assert get_planet_name(MERCURY) == "Mercury"
        assert get_planet_name(VENUS) == "Venus"
        assert get_planet_name(MARS) == "Mars"
        assert get_planet_name(JUPITER) == "Jupiter"
        assert get_planet_name(SATURN) == "Saturn"
        assert get_planet_name(URANUS) == "Uranus"
        assert get_planet_name(NEPTUNE) == "Neptune"
        assert get_planet_name(PLUTO) == "Pluto"

    def test_earth(self):
        """Test Earth name lookup."""
        assert get_planet_name(EARTH) == "Earth"

    def test_lunar_nodes(self):
        """Test lunar node name lookups (the reference ephemeris uses lowercase)."""
        assert get_planet_name(MEAN_NODE) == "mean Node"
        assert get_planet_name(TRUE_NODE) == "true Node"

    def test_lunar_apogee(self):
        """Test lunar apogee name lookups (the reference ephemeris uses lowercase/abbreviated)."""
        assert get_planet_name(MEAN_APOG) == "mean Apogee"
        assert get_planet_name(OSCU_APOG) == "osc. Apogee"

    def test_unknown_planet_id(self):
        """Unknown IDs return the compatibility API's empty string."""
        result = get_planet_name(9999)
        assert result == ""

    def test_negative_planet_id(self):
        """Test negative planet ID returns descriptive string."""
        result = get_planet_name(-1)
        assert result == ""

    def test_constants_without_builtin_elements_have_no_name(self):
        """IDs 55-58 are constants only and have no built-in names."""
        for body_id in (VULCAN, WHITE_MOON, PROSERPINA, WALDEMATH):
            assert get_planet_name(body_id) == ""

    def test_return_type(self):
        """Test that get_planet_name always returns a string."""
        assert isinstance(get_planet_name(SUN), str)
        assert isinstance(get_planet_name(MOON), str)
        assert isinstance(get_planet_name(9999), str)
