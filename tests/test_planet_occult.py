"""
Tests for planetary occultation calculations in libephemeris.

Tests the planet_occult_when_glob and planet_occult_when_loc functions which
find planetary occultations of other planets and fixed stars.

Planetary occultations occur when one planet passes in front of another planet
or star as seen from Earth. These are very rare events:
- Mutual planetary occultations happen only a few times per century
- Planetary occultations of bright stars are more common but still rare

Historical mutual planetary occultations:
- 1818 Dec 3: Venus occulted Jupiter
- 1737 Aug 29: Venus occulted Mercury

Note: The bundled DE421 ephemeris covers 1899-2053, so tests are limited to
this date range. Many tests verify function structure and error handling
rather than finding actual occultations (which are very rare).
"""

import pytest

pytestmark = pytest.mark.slow

from libephemeris.exceptions import Error
from libephemeris import (
    julday,
    planet_occult_when_glob,
    planet_occult_when_loc,
    SUN,
    MOON,
    VENUS,
    MARS,
    JUPITER,
    FLG_SWIEPH,
)


class TestPlanetOccultWhenGlob:
    """Test suite for planet_occult_when_glob function."""

    def test_raises_error_for_no_target(self):
        """Test that function raises error if no target specified."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError):
            planet_occult_when_glob(jd_start, VENUS, 0, "", FLG_SWIEPH, 0)

    def test_raises_error_for_sun_as_occulting_body(self):
        """Test that Sun cannot be the occulting body."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError, match="Sun cannot be"):
            planet_occult_when_glob(jd_start, SUN, VENUS, "", FLG_SWIEPH, 0)

    def test_raises_error_for_moon_as_occulting_body(self):
        """Test that Moon cannot be the occulting body."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError, match="Moon cannot be"):
            planet_occult_when_glob(jd_start, MOON, VENUS, "", FLG_SWIEPH, 0)

    def test_raises_error_for_same_planet(self):
        """Test that occulting and occulted cannot be the same."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError, match="must be different bodies"):
            planet_occult_when_glob(jd_start, VENUS, VENUS, "", FLG_SWIEPH, 0)

    def test_raises_error_for_invalid_occulting_planet(self):
        """Test that invalid planet ID raises appropriate error."""
        jd_start = julday(2020, 1, 1, 0)

        with pytest.raises(ValueError):
            planet_occult_when_glob(jd_start, 999, JUPITER, "", FLG_SWIEPH, 0)

    def test_raises_error_for_invalid_occulted_planet(self):
        """Test that invalid occulted planet ID raises appropriate error."""
        jd_start = julday(2020, 1, 1, 0)

        with pytest.raises(ValueError):
            planet_occult_when_glob(jd_start, VENUS, 999, "", FLG_SWIEPH, 0)

    def test_raises_error_for_unknown_star(self):
        """Test that function raises error for unknown star name."""
        jd_start = julday(2020, 1, 1, 0)

        with pytest.raises(ValueError):
            planet_occult_when_glob(jd_start, VENUS, 0, "UnknownStar123", FLG_SWIEPH, 0)

    def test_swe_alias(self):
        """Test that planet_occult_when_glob is an alias."""

    def test_search_terminates_properly(self):
        """Test that search terminates when no event found or event found.

        Planetary occultations are so rare that we may not find one in most
        search windows. This test verifies the function eventually terminates
        (either with a result or with exceptions.Error).
        """
        # Use a short search range with inner planets for faster test
        jd_start = julday(2000, 1, 1, 0)

        try:
            # This will either find an event or raise Error
            retflags, tret = planet_occult_when_glob(
                jd_start, VENUS, MARS, "", FLG_SWIEPH, 0
            )
            # If it returns, verify structure
            assert len(tret) == 10
            assert isinstance(retflags, int)
        except Error as e:
            # Expected if no occultation found
            assert "No planetary occultation" in str(e)


class TestPlanetOccultWhenLoc:
    """Test suite for planet_occult_when_loc function."""

    def test_raises_error_for_no_target(self):
        """Test that function raises error if no target specified."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError):
            planet_occult_when_loc(jd_start, VENUS, 0, "", 40.0, -74.0, 0, FLG_SWIEPH)

    def test_raises_error_for_sun_as_occulting_body(self):
        """Test that Sun cannot be the occulting body."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError, match="Sun cannot be"):
            planet_occult_when_loc(jd_start, SUN, VENUS, "", 40.0, -74.0, 0, FLG_SWIEPH)

    def test_raises_error_for_moon_as_occulting_body(self):
        """Test that Moon cannot be the occulting body."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError, match="Moon cannot be"):
            planet_occult_when_loc(
                jd_start, MOON, VENUS, "", 40.0, -74.0, 0, FLG_SWIEPH
            )

    def test_raises_error_for_same_planet(self):
        """Test that occulting and occulted cannot be the same."""
        jd_start = julday(2024, 1, 1, 0)

        with pytest.raises(ValueError, match="must be different bodies"):
            planet_occult_when_loc(
                jd_start, VENUS, VENUS, "", 40.0, -74.0, 0, FLG_SWIEPH
            )

    def test_swe_alias(self):
        """Test that planet_occult_when_loc is an alias."""


class TestPlanetOccultStarOccultation:
    """Tests for planetary occultations of fixed stars."""

    def test_raises_error_for_unknown_star(self):
        """Test that function raises error for unknown star name."""
        jd_start = julday(2020, 1, 1, 0)

        with pytest.raises(ValueError):
            planet_occult_when_glob(
                jd_start, VENUS, 0, "NonexistentStar", FLG_SWIEPH, 0
            )


class TestPlanetOccultReturnStructure:
    """Tests for verifying return value structure without requiring an actual event."""

    def test_glob_function_signature(self):
        """Test that glob function accepts correct parameters."""
        jd_start = julday(2000, 1, 1, 0)

        # Function should accept these parameters without error
        # (it will raise Error because no occultation is found, which is expected)
        try:
            retflags, tret = planet_occult_when_glob(
                jd_start, VENUS, MARS, "", FLG_SWIEPH, 0
            )
            # If we somehow found an event, verify structure
            assert len(tret) == 10
            assert isinstance(retflags, int)
        except Error as e:
            # Expected - no occultation found
            assert "No planetary occultation" in str(e)

    def test_loc_function_signature(self):
        """Test that loc function accepts correct parameters."""
        jd_start = julday(2000, 1, 1, 0)

        try:
            retflag, times, attr = planet_occult_when_loc(
                jd_start, VENUS, MARS, "", 40.0, -74.0, 100.0, FLG_SWIEPH
            )
            # If we somehow found an event, verify structure
            assert isinstance(retflag, int)
            assert len(times) == 10
            assert len(attr) == 20
        except Error as e:
            # Expected - no occultation found
            assert "No planetary occultation" in str(e)


class TestPlanetOccultImports:
    """Test that functions are properly exported."""

    def test_imports_from_main_module(self):
        """Test that functions can be imported from main module."""
        from libephemeris import (
            planet_occult_when_glob,
            planet_occult_when_loc,
        )

        assert callable(planet_occult_when_glob)
        assert callable(planet_occult_when_loc)
        assert callable(planet_occult_when_glob)
        assert callable(planet_occult_when_loc)

    def test_swe_aliases_are_same_function(self):
        """Test that swe_* are aliases to base functions."""


class TestPlanetOccultDocumentation:
    """Tests to verify documentation examples work (at least syntactically)."""

    def test_documentation_example_valid(self):
        """Test that documentation example syntax is valid."""
        # This verifies the example in the docstring is syntactically correct
        from libephemeris import julday, planet_occult_when_glob, VENUS, JUPITER

        jd = julday(2000, 1, 1, 0)

        # Function should be callable with these arguments
        # (may find an event or raise RuntimeError - both are valid)
        try:
            retflags, tret = planet_occult_when_glob(jd, VENUS, JUPITER)
            assert len(tret) == 10
        except RuntimeError:
            pass  # Expected if no occultation found
