"""
Test that error messages keep a stable, pattern-matchable shape.

Client code sometimes inspects the text of an exception to decide what went
wrong. These tests pin the shape of the messages that such code is most likely
to look at, so a rewording is a deliberate choice rather than an accident:

- "body id {n} ..." for a body id that has no model or no elements
- "no fixed star matches the search string '{key}'" for stars not found
- "polar circle" and a suggested fallback system for houses in polar regions
"""

import pytest
import re
from libephemeris.exceptions import Error


class TestErrorMessageFormat:
    """Test the shape of the messages client code may pattern-match."""

    def test_star_not_found_format(self):
        """
        Star-not-found messages echo the search key of each family.

        v1 (fixstar/fixstar_ut/fixstar_mag) echoes the raw search string;
        v2 (fixstar2*) echoes the lowercased, space-stripped key.
        """
        from libephemeris.fixed_stars import fixstar_ut, fixstar2_ut

        # Test fixstar_ut (v1 message shape, raw echo)
        with pytest.raises(
            Error, match="no fixed star matches the search string 'NonExistentStar123'"
        ):
            fixstar_ut("NonExistentStar123", 2451545.0, 0)

        # Test fixstar2_ut (v2 message shape, normalized echo)
        with pytest.raises(
            Error, match="no fixed star matches the search string 'nonexistentstar123'"
        ):
            fixstar2_ut("NonExistentStar123", 2451545.0, 0)

    def test_polar_circle_error_format_houses(self):
        """
        Error messages for polar circle houses should contain 'polar circle'
        and name Porphyry as the suggested alternative.
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


class TestUnknownBodyMessages:
    """Test that unknown-body messages name the offending id."""

    def test_minor_body_not_found_message(self):
        """A body id outside the element table is named in the ValueError."""
        from libephemeris.minor_bodies import calc_minor_body_heliocentric

        with pytest.raises(ValueError) as excinfo:
            calc_minor_body_heliocentric(999999, 2451545.0)

        assert "body id 999999 has no entry in the minor-body element table" == str(
            excinfo.value
        )

    def test_star_resolve_error_formats(self):
        """Each resolver family echoes its own search key."""
        from libephemeris.fixed_stars import _resolve_star_id, _resolve_star2

        # _resolve_star_id (v1): raw echo
        star_id, error, _ = _resolve_star_id("NonExistentStar")
        assert star_id == -1
        assert error == "no fixed star matches the search string 'NonExistentStar'"

        # _resolve_star2 (v2): normalized echo
        entry, error = _resolve_star2("NonExistentStar")
        assert entry is None
        assert error == "no fixed star matches the search string 'nonexistentstar'"


class TestPatternMatchingCompatibility:
    """
    Test that client code can pattern match on error messages.

    These tests simulate what client code that inspects messages might do.
    """

    def test_can_detect_unknown_body_pattern(self):
        """Client code can extract the offending body id from the message."""
        from libephemeris.minor_bodies import calc_minor_body_heliocentric

        with pytest.raises(ValueError) as excinfo:
            calc_minor_body_heliocentric(999999, 2451545.0)
        error_msg = str(excinfo.value)
        # Pattern matching that client code might use
        body_pattern = re.search(r"body id (\d+) has no entry", error_msg)
        assert body_pattern is not None, f"Pattern not found in: {error_msg}"
        assert body_pattern.group(1) == "999999"

    def test_can_detect_star_not_found_pattern(self):
        """Client code can extract the search key from the not-found messages."""
        from libephemeris.fixed_stars import fixstar_ut, fixstar2_ut

        # v1: raw echo
        with pytest.raises(Error) as excinfo:
            fixstar_ut("FakeStar", 2451545.0, 0)
        star_pattern = re.search(
            r"no fixed star matches the search string '(\w+)'", str(excinfo.value)
        )
        assert star_pattern is not None, f"Pattern not found in: {excinfo.value}"
        assert star_pattern.group(1) == "FakeStar"

        # v2: normalized echo
        with pytest.raises(Error) as excinfo:
            fixstar2_ut("FakeStar", 2451545.0, 0)
        star_pattern = re.search(
            r"no fixed star matches the search string '(\w+)'", str(excinfo.value)
        )
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
