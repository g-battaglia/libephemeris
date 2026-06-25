"""Tests for UnknownBodyError exception handling.

These tests verify that clear error messages are provided when unknown
body IDs are requested in calculation functions.
"""

import pytest
import libephemeris as ephem
from libephemeris.exceptions import UnknownBodyError


class TestUnknownBodyError:
    """Test UnknownBodyError exception class."""

    def test_exception_attributes(self):
        """Test that exception has correct attributes."""
        err = UnknownBodyError(message="Test message", body_id=99999)
        assert err.body_id == 99999
        assert err.message == "Test message"
        assert str(err) == "Test message"

    def test_exception_repr(self):
        """Test exception repr format."""
        err = UnknownBodyError(message="Test message", body_id=99999)
        repr_str = repr(err)
        assert "UnknownBodyError" in repr_str
        assert "99999" in repr_str

    def test_inherits_from_error(self):
        """Test that UnknownBodyError inherits from Error."""
        err = UnknownBodyError(message="Test", body_id=1)
        assert isinstance(err, ephem.Error)


class TestCalcUtUnknownBody:
    """Test calc_ut with unknown body IDs."""

    def test_unknown_body_raises_error(self):
        """Test that unknown body ID raises UnknownBodyError."""
        jd = 2451545.0  # J2000.0
        unknown_id = 8888

        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc_ut(jd, unknown_id, 0)

        assert exc_info.value.body_id == unknown_id
        assert "8888" in str(exc_info.value)
        assert "Unknown body ID" in str(exc_info.value)

    def test_negative_unknown_id(self):
        """Test that negative unknown body ID raises error."""
        jd = 2451545.0
        # -100 is not a valid ID (south nodes use small negatives like -10, -11)
        unknown_id = -100

        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc_ut(jd, unknown_id, 0)

        assert exc_info.value.body_id == unknown_id

    def test_gap_between_ranges(self):
        """Test that IDs in gaps between valid ranges raise error."""
        jd = 2451545.0
        # IDs 23-39 are in a gap (between NPLANETS=23 and FICT_OFFSET=40)
        gap_id = 30

        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc_ut(jd, gap_id, 0)

        assert exc_info.value.body_id == gap_id

    def test_error_message_contains_guidance(self):
        """Test that error message provides useful guidance."""
        jd = 2451545.0
        unknown_id = 8888

        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc_ut(jd, unknown_id, 0)

        msg = str(exc_info.value)
        # Should mention supported body types
        assert "standard planets" in msg.lower() or "planets" in msg.lower()
        assert "constants" in msg.lower()


class TestCalcUnknownBody:
    """Test calc with unknown body IDs."""

    def test_unknown_body_raises_error(self):
        """Test that unknown body ID raises UnknownBodyError in calc.

        88888 would be asteroid 78888 (any AST_OFFSET + N is attempted
        via SPK/SBDB since WS7), so use an id below the asteroid range
        that no dispatch claims.
        """
        jd = 2451545.0
        unknown_id = 8888

        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc(jd, unknown_id, 0)

        assert exc_info.value.body_id == unknown_id


class TestCalcPctrUnknownBody:
    """Test calc_pctr with unknown body IDs."""

    def test_unknown_target_raises_error(self):
        """Test that unknown target body ID raises UnknownBodyError."""
        jd = 2451545.0
        unknown_id = 77777

        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc_pctr(jd, unknown_id, ephem.EARTH, 0)

        assert exc_info.value.body_id == unknown_id
        assert "target" in str(exc_info.value).lower()

    def test_unknown_observer_raises_error(self):
        """Test that unknown observer body ID raises UnknownBodyError."""
        jd = 2451545.0
        unknown_id = 66666

        with pytest.raises(UnknownBodyError) as exc_info:
            ephem.calc_pctr(jd, ephem.MOON, unknown_id, 0)

        assert exc_info.value.body_id == unknown_id
        assert "observer" in str(exc_info.value).lower()


class TestValidBodiesStillWork:
    """Ensure valid body calculations still work correctly."""

    def test_sun_calculation(self):
        """Test that Sun calculation works."""
        jd = 2451545.0
        pos, flags = ephem.calc_ut(jd, ephem.SUN, ephem.FLG_SPEED)
        # Sun should be near 280 degrees at J2000.0
        assert 270 < pos[0] < 290

    def test_moon_calculation(self):
        """Test that Moon calculation works."""
        jd = 2451545.0
        pos, flags = ephem.calc_ut(jd, ephem.MOON, ephem.FLG_SPEED)
        # Moon should have valid longitude
        assert 0 <= pos[0] < 360

    def test_mars_calculation(self):
        """Test that Mars calculation works."""
        jd = 2451545.0
        pos, flags = ephem.calc_ut(jd, ephem.MARS, ephem.FLG_SPEED)
        assert 0 <= pos[0] < 360

    def test_lunar_node_calculation(self):
        """Test that lunar node calculation works."""
        jd = 2451545.0
        pos, flags = ephem.calc_ut(jd, ephem.MEAN_NODE, 0)
        assert 0 <= pos[0] < 360

    def test_pctr_valid_bodies(self):
        """Test that planet-centric calculation works for valid bodies."""
        jd = 2451545.0
        pos, flags = ephem.calc_pctr(jd, ephem.MOON, ephem.MARS, 0)
        assert 0 <= pos[0] < 360


class TestExceptionCanBeCaught:
    """Test that exception can be caught properly."""

    def test_catch_as_unknown_body_error(self):
        """Test catching as UnknownBodyError."""
        jd = 2451545.0
        caught = False
        try:
            ephem.calc_ut(jd, 8888, 0)
        except UnknownBodyError:
            caught = True
        assert caught

    def test_catch_as_generic_error(self):
        """Test catching as generic Error."""
        jd = 2451545.0
        caught = False
        try:
            ephem.calc_ut(jd, 8888, 0)
        except ephem.Error:
            caught = True
        assert caught

    def test_access_body_id_from_exception(self):
        """Test accessing body_id attribute from caught exception."""
        jd = 2451545.0
        try:
            ephem.calc_ut(jd, 12345, 0)
        except UnknownBodyError as e:
            assert e.body_id == 12345
