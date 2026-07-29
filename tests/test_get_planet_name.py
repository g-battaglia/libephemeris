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
    VULKANUS,
    NEPTUNE_LEVERRIER,
    NEPTUNE_ADAMS,
    PLUTO_LOWELL,
    PLUTO_PICKERING,
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

    def test_symbolic_points_have_names(self):
        """IDs 55-58 are computed symbolic points, so they carry names.

        The runtime returns positions for Vulcan, the White Moon (Selena),
        Proserpina and Waldemath, so the metadata contract exposes their
        published names rather than an empty string.
        """
        assert get_planet_name(VULCAN) == "Vulcan"
        assert get_planet_name(WHITE_MOON) == "Selena/White Moon"
        assert get_planet_name(PROSERPINA) == "Proserpina"
        assert get_planet_name(WALDEMATH) == "Waldemath"

    def test_vulcanus_hamburg_spelling(self):
        """The seventh Uranian planet displays as 'Vulcanus' (Hamburg School)."""
        assert get_planet_name(VULKANUS) == "Vulcanus"

    def test_historical_predictions_annotated(self):
        """Historical predictions carry the parenthetical target planet."""
        assert get_planet_name(NEPTUNE_LEVERRIER) == "Leverrier (Neptune)"
        assert get_planet_name(NEPTUNE_ADAMS) == "Adams (Neptune)"
        assert get_planet_name(PLUTO_LOWELL) == "Lowell (Pluto)"
        assert get_planet_name(PLUTO_PICKERING) == "Pickering (Pluto)"

    def test_return_type(self):
        """Test that get_planet_name always returns a string."""
        assert isinstance(get_planet_name(SUN), str)
        assert isinstance(get_planet_name(MOON), str)
        assert isinstance(get_planet_name(9999), str)
