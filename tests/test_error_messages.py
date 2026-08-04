"""
Test that error messages match the reference ephemeris format for client code compatibility.

This module tests that libephemeris error messages use the same format as the reference ephemeris,
allowing existing client code that does pattern matching on error messages to work
correctly with libephemeris.

Key formats verified:
- "illegal planet number {n}." for invalid planet IDs
- "could not find star name {name}" for stars not found
- "within polar circle, switched to Porphyry" for houses in polar regions
"""

import pytest
import re
from libephemeris.exceptions import Error


class TestErrorMessageFormat:
    """Test that error messages match the reference ephemeris format."""

    def test_illegal_planet_number_format_in_message(self):
        """
        Error messages for invalid planets should contain 'illegal planet number'.

        The reference ephemeris uses format: "illegal planet number {n}."
        """

        # These should fail silently (return zeros) rather than raise
        # because invalid planets in _calc_body just return zeros
        # But functions that do validation should use "illegal planet number"
        pass  # Validation happens in specific functions, tested below

    def test_star_not_found_format(self):
        """
        Star-not-found messages use the measured per-family reference formats.

        v1 (fixstar/fixstar_ut/fixstar_mag): "star {name} not found", echoing
        the raw search string. v2 (fixstar2*): "could not find star name
        {key}", echoing the lowercased space-stripped key.
        """
        from libephemeris.fixed_stars import fixstar_ut, fixstar2_ut

        # Test fixstar_ut (v1 message shape, raw echo)
        with pytest.raises(Error, match="star NonExistentStar123 not found"):
            fixstar_ut("NonExistentStar123", 2451545.0, 0)

        # Test fixstar2_ut (v2 message shape, normalized echo)
        with pytest.raises(Error, match="could not find star name nonexistentstar123"):
            fixstar2_ut("NonExistentStar123", 2451545.0, 0)

    def test_polar_circle_error_format_houses(self):
        """
        Error messages for polar circle houses should contain 'polar circle'.

        libephemeris uses a more descriptive format than the reference ephemeris:
        "(within Northern polar circle)" with suggestions for alternatives.
        """
        import libephemeris as ephem
        from libephemeris.exceptions import PolarCircleError

        # Test houses with high latitude (polar circle condition)
        # For Placidus system at high latitude
        with pytest.raises(PolarCircleError) as excinfo:
            ephem.houses(2451545.0, 85.0, 0.0, ord("P"))  # Placidus at 85° latitude

        error_msg = str(excinfo.value)
        assert "polar circle" in error_msg
        assert "Porphyry" in error_msg  # Suggested as alternative

    def test_polar_circle_error_format_houses_armc(self):
        """
        Error messages for polar circle in houses_armc should indicate polar region.
        """
        import libephemeris as ephem
        from libephemeris.exceptions import PolarCircleError

        # Test houses_armc with high latitude
        with pytest.raises(PolarCircleError) as excinfo:
            ephem.houses_armc(0.0, 85.0, 23.44, ord("P"))  # Placidus at 85° latitude

        error_msg = str(excinfo.value)
        assert "polar circle" in error_msg
        assert "Porphyry" in error_msg  # Suggested as alternative


class TestIllegalPlanetMessages:
    """Test that 'illegal planet number' format is used for invalid planets."""

    def test_minor_body_not_found_message(self):
        """Minor body not found should use 'illegal planet number' format."""
        from libephemeris.minor_bodies import calc_minor_body_heliocentric

        with pytest.raises(ValueError) as excinfo:
            calc_minor_body_heliocentric(999999, 2451545.0)

        assert "illegal planet number" in str(excinfo.value)

    def test_star_resolve_error_formats(self):
        """Each resolver family uses its measured reference message shape."""
        from libephemeris.fixed_stars import _resolve_star_id, _resolve_star2

        # _resolve_star_id (v1): "star {raw} not found"
        star_id, error, _ = _resolve_star_id("NonExistentStar")
        assert star_id == -1
        assert error == "star NonExistentStar not found"

        # _resolve_star2 (v2): "could not find star name {normalized}"
        entry, error = _resolve_star2("NonExistentStar")
        assert entry is None
        assert error == "could not find star name nonexistentstar"


class TestPatternMatchingCompatibility:
    """
    Test that client code can pattern match on error messages.

    These tests simulate what existing reference ephemeris client code might do.
    """

    def test_can_detect_illegal_planet_pattern(self):
        """Client code should be able to detect 'illegal planet number' pattern."""
        from libephemeris.minor_bodies import calc_minor_body_heliocentric

        try:
            calc_minor_body_heliocentric(999999, 2451545.0)
        except ValueError as e:
            error_msg = str(e)
            # Pattern matching that client code might use
            illegal_pattern = re.search(r"illegal planet number (\d+)", error_msg)
            assert illegal_pattern is not None, f"Pattern not found in: {error_msg}"
            assert illegal_pattern.group(1) == "999999"

    def test_can_detect_star_not_found_pattern(self):
        """Client code can pattern-match the per-family not-found messages."""
        from libephemeris.fixed_stars import fixstar_ut, fixstar2_ut

        # v1: "star {raw} not found"
        with pytest.raises(Error) as excinfo:
            fixstar_ut("FakeStar", 2451545.0, 0)
        star_pattern = re.search(r"star (\w+) not found", str(excinfo.value))
        assert star_pattern is not None, f"Pattern not found in: {excinfo.value}"
        assert star_pattern.group(1) == "FakeStar"

        # v2: "could not find star name {normalized}"
        with pytest.raises(Error) as excinfo:
            fixstar2_ut("FakeStar", 2451545.0, 0)
        star_pattern = re.search(r"could not find star name (\w+)", str(excinfo.value))
        assert star_pattern is not None, f"Pattern not found in: {excinfo.value}"
        assert star_pattern.group(1) == "fakestar"

    def test_can_detect_polar_circle_pattern(self):
        """Client code should be able to detect 'polar circle' pattern."""
        import libephemeris as ephem
        from libephemeris.exceptions import PolarCircleError

        try:
            ephem.houses(2451545.0, 85.0, 0.0, ord("P"))
        except PolarCircleError as e:
            error_msg = str(e)
            # Pattern matching that client code might use
            assert "polar circle" in error_msg
            assert "Porphyry" in error_msg
