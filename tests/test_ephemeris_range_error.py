"""
Tests for EphemerisRangeError improved error messages.

This module tests that when calculations fail due to dates outside the
ephemeris range, clear error messages are provided with:
- The requested Julian Day number
- The supported date range in both JD and calendar format
- The body being calculated
- The ephemeris file in use
"""

import pytest
import libephemeris as eph
from libephemeris.exceptions import EphemerisRangeError


class TestEphemerisRangeErrorClass:
    """Test the EphemerisRangeError exception class itself."""

    def test_exception_has_all_attributes(self):
        """EphemerisRangeError should have all documented attributes."""
        err = EphemerisRangeError(
            message="Test error message",
            requested_jd=3000000.0,
            start_jd=2287184.5,
            end_jd=2688976.5,
            start_date="1549-12-31",
            end_date="2650-01-25",
            body_id=0,
            body_name="Sun",
            ephemeris_file="de440.bsp",
        )

        assert err.requested_jd == 3000000.0
        assert err.start_jd == 2287184.5
        assert err.end_jd == 2688976.5
        assert err.start_date == "1549-12-31"
        assert err.end_date == "2650-01-25"
        assert err.body_id == 0
        assert err.body_name == "Sun"
        assert err.ephemeris_file == "de440.bsp"
        assert str(err) == "Test error message"

    def test_exception_inherits_from_error(self):
        """EphemerisRangeError should inherit from Error."""
        err = EphemerisRangeError("Test")
        assert isinstance(err, eph.Error)
        assert isinstance(err, Exception)

    def test_exception_repr(self):
        """EphemerisRangeError should have informative repr."""
        err = EphemerisRangeError(
            message="Test",
            requested_jd=3000000.0,
            body_id=0,
        )
        repr_str = repr(err)
        assert "EphemerisRangeError" in repr_str
        assert "3000000.0" in repr_str


class TestCalcUtRangeError:
    """Test that calc_ut raises EphemerisRangeError for out-of-range dates."""

    def test_far_future_date_raises_error(self):
        """Calculation for far future date should raise EphemerisRangeError."""
        # JD 3000000 is approximately year 3501 - way beyond DE440 range
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_ut(far_future_jd, eph.SUN, 0)

        err = exc_info.value
        # Check that error contains useful information
        assert err.requested_jd == far_future_jd
        assert err.body_id == eph.SUN
        assert err.body_name == "Sun"
        assert err.start_jd is not None
        assert err.end_jd is not None
        assert err.ephemeris_file is not None

    def test_far_past_date_raises_error(self):
        """Calculation for far past date should raise EphemerisRangeError."""
        # JD 1000000 is approximately year -1975 - before DE440 range
        far_past_jd = 1000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_ut(far_past_jd, eph.MOON, 0)

        err = exc_info.value
        assert err.requested_jd == far_past_jd
        assert err.body_id == eph.MOON
        assert err.body_name == "Moon"

    def test_error_message_contains_jd_info(self):
        """Error message should contain requested JD and supported range."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_ut(far_future_jd, eph.MARS, 0)

        message = str(exc_info.value)
        # Check message contains key information
        assert "3000000" in message  # Requested JD
        assert "Mars" in message  # Body name
        assert "outside ephemeris range" in message.lower()
        assert "JD" in message  # Should include JD in supported range

    def test_error_message_contains_date_strings(self):
        """Error message should contain date strings for readability."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_ut(far_future_jd, eph.SUN, 0)

        message = str(exc_info.value)
        # Should contain calendar dates
        assert "1549-12-31" in message or "1549" in message
        assert "2650" in message

    def test_error_message_contains_ephemeris_file(self):
        """Error message should indicate which ephemeris file is in use."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_ut(far_future_jd, eph.SUN, 0)

        message = str(exc_info.value)
        # Should contain ephemeris file name
        assert ".bsp" in message or "de" in message.lower()


class TestCalcRangeError:
    """Test that calc raises EphemerisRangeError for out-of-range dates."""

    def test_far_future_date_raises_error(self):
        """calc should also raise EphemerisRangeError."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc(far_future_jd, eph.JUPITER, 0)

        err = exc_info.value
        assert err.requested_jd == far_future_jd
        assert err.body_id == eph.JUPITER
        assert err.body_name == "Jupiter"


class TestCalcPctrRangeError:
    """Test that calc_pctr raises EphemerisRangeError for out-of-range dates."""

    def test_far_future_date_raises_error(self):
        """calc_pctr should also raise EphemerisRangeError."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_pctr(far_future_jd, eph.MOON, eph.MARS, 0)

        err = exc_info.value
        assert err.requested_jd == far_future_jd
        assert err.body_id == eph.MOON


class TestExceptionCanBeCaught:
    """Test that the exception can be caught in various ways."""

    def test_catch_as_ephemeris_range_error(self):
        """Can catch as EphemerisRangeError."""
        caught = False
        try:
            eph.calc_ut(3000000.0, eph.SUN, 0)
        except EphemerisRangeError:
            caught = True

        assert caught

    def test_catch_as_error(self):
        """Can catch as Error (base class)."""
        caught = False
        try:
            eph.calc_ut(3000000.0, eph.SUN, 0)
        except eph.Error:
            caught = True

        assert caught

    def test_catch_as_exception(self):
        """Can catch as Exception."""
        caught = False
        try:
            eph.calc_ut(3000000.0, eph.SUN, 0)
        except Exception:
            caught = True

        assert caught


class TestEdgeCases:
    """Test edge cases for date range errors."""

    def test_date_just_inside_range_works(self):
        """Date just inside ephemeris range should work."""
        # J2000.0 is well within DE440 range
        jd = 2451545.0  # 2000-01-01

        pos, _ = eph.calc_ut(jd, eph.SUN, 0)
        assert 0 <= pos[0] < 360  # Valid longitude

    def test_various_planets_get_correct_names(self):
        """Different planets should get correct names in error messages."""
        far_future_jd = 3000000.0

        test_cases = [
            (eph.SUN, "Sun"),
            (eph.MOON, "Moon"),
            (eph.MERCURY, "Mercury"),
            (eph.VENUS, "Venus"),
            (eph.MARS, "Mars"),
            (eph.JUPITER, "Jupiter"),
            (eph.SATURN, "Saturn"),
        ]

        for body_id, expected_name in test_cases:
            with pytest.raises(EphemerisRangeError) as exc_info:
                eph.calc_ut(far_future_jd, body_id, 0)

            assert exc_info.value.body_name == expected_name
            assert expected_name in str(exc_info.value)


class TestExportedFromModule:
    """Test that EphemerisRangeError is properly exported."""

    def test_exported_from_libephemeris(self):
        """EphemerisRangeError should be accessible from main module."""
        assert hasattr(eph, "EphemerisRangeError")
        assert eph.EphemerisRangeError is EphemerisRangeError

    def test_in_all_list(self):
        """EphemerisRangeError should be in __all__."""
        assert "EphemerisRangeError" in eph.__all__
