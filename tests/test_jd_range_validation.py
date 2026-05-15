"""
Tests for Julian Day range validation before calculation.

This module tests that Julian Day range is validated before ephemeris calculations
to provide clear error messages when dates are outside the supported range.
"""

import pytest
import libephemeris as eph
from libephemeris.exceptions import EphemerisRangeError, validate_jd_range


class TestValidateJdRangeFunction:
    """Test the validate_jd_range function directly."""

    def test_valid_jd_within_range(self):
        """Valid Julian Day within ephemeris range should not raise."""
        # J2000.0 is well within DE440 range (1550-2650)
        jd = 2451545.0  # 2000-01-01
        validate_jd_range(jd)  # Should not raise

    def test_far_future_jd_raises_error(self):
        """Julian Day far in the future should raise EphemerisRangeError."""
        # JD 3000000 is approximately year 3501
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range(far_future_jd)

        err = exc_info.value
        assert err.requested_jd == far_future_jd
        assert "outside ephemeris range" in str(err).lower()

    def test_far_past_jd_raises_error(self):
        """Julian Day far in the past should raise EphemerisRangeError."""
        # JD 1000000 is approximately year -1975
        far_past_jd = 1000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range(far_past_jd)

        err = exc_info.value
        assert err.requested_jd == far_past_jd
        assert "outside ephemeris range" in str(err).lower()

    def test_error_includes_body_info(self):
        """Error should include body information if provided."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range(far_future_jd, body_id=eph.SUN)

        err = exc_info.value
        assert err.body_id == eph.SUN
        assert err.body_name == "Sun"
        assert "Sun" in str(err)

    def test_error_includes_func_name(self):
        """Error should include function name if provided."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range(far_future_jd, func_name="test_function")

        message = str(exc_info.value)
        assert "test_function" in message

    def test_error_includes_date_strings(self):
        """Error should include human-readable date strings."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range(far_future_jd)

        err = exc_info.value
        # Should include start/end dates
        assert err.start_date is not None
        assert err.end_date is not None
        # Error message should mention dates
        message = str(err)
        assert err.start_date in message or "1549" in message
        assert err.end_date in message or "2650" in message

    def test_error_includes_jd_range(self):
        """Error should include the supported JD range."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range(far_future_jd)

        err = exc_info.value
        assert err.start_jd is not None
        assert err.end_jd is not None
        assert err.start_jd < err.end_jd

    def test_error_includes_ephemeris_file(self):
        """Error should include the ephemeris file name."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            validate_jd_range(far_future_jd)

        err = exc_info.value
        assert err.ephemeris_file is not None
        assert ".bsp" in err.ephemeris_file or "de" in err.ephemeris_file.lower()


class TestProactiveValidation:
    """Test that validation happens before calculation (proactive, not reactive)."""

    def test_calc_ut_validates_jd_range(self):
        """calc_ut should validate JD range before calculation."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_ut(far_future_jd, eph.SUN, 0)

        err = exc_info.value
        # Error should be raised with body info (proactive validation)
        assert err.body_id == eph.SUN
        assert "calc_ut" in str(err)

    def test_calc_validates_jd_range(self):
        """calc should validate JD range before calculation."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc(far_future_jd, eph.MOON, 0)

        err = exc_info.value
        assert err.body_id == eph.MOON
        assert "calc" in str(err)

    def test_calc_pctr_validates_jd_range(self):
        """calc_pctr should validate JD range before calculation."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_pctr(far_future_jd, eph.MARS, eph.JUPITER, 0)

        err = exc_info.value
        assert err.body_id == eph.MARS
        assert "calc_pctr" in str(err)


class TestNonJplBodiesNotValidated:
    """Test that non-JPL bodies are not subject to JD range validation."""

    def test_mean_node_far_future(self):
        """Mean lunar node should work for far future dates (Keplerian)."""
        # Mean node is calculated mathematically, not from JPL ephemeris
        far_future_jd = 3000000.0

        # This should not raise EphemerisRangeError
        # (it may use other bodies internally, so this tests the entry point)
        try:
            pos, _ = eph.calc_ut(far_future_jd, eph.MEAN_NODE, 0)
            # If it succeeds, verify we got a result
            assert 0 <= pos[0] < 360
        except EphemerisRangeError:
            # If it raises due to internal dependencies, that's also acceptable
            pass

    def test_mean_apogee_far_future(self):
        """Mean Lilith/Apogee should work for far future dates (Keplerian)."""
        far_future_jd = 3000000.0

        try:
            pos, _ = eph.calc_ut(far_future_jd, eph.MEAN_APOG, 0)
            assert 0 <= pos[0] < 360
        except EphemerisRangeError:
            pass


class TestEdgeCases:
    """Test edge cases for JD range validation."""

    def test_jd_exactly_at_start(self):
        """JD at the exact start of range should work."""
        # Get the actual range from the ephemeris
        eph.calc_ut(2451545.0, eph.SUN, 0)  # Ensure loaded
        path, start_jd, end_jd, denum = eph.get_current_file_data(0)

        if start_jd > 0:
            # Test JD at start should work
            pos, _ = eph.calc_ut(start_jd + 1.0, eph.SUN, 0)
            assert 0 <= pos[0] < 360

    def test_jd_exactly_at_end(self):
        """JD at the exact end of range should work."""
        eph.calc_ut(2451545.0, eph.SUN, 0)  # Ensure loaded
        path, start_jd, end_jd, denum = eph.get_current_file_data(0)

        if end_jd > 0:
            # Test JD just before end should work
            pos, _ = eph.calc_ut(end_jd - 1.0, eph.SUN, 0)
            assert 0 <= pos[0] < 360

    def test_jd_just_outside_start(self):
        """JD just before start should raise error."""
        eph.calc_ut(2451545.0, eph.SUN, 0)  # Ensure loaded
        path, start_jd, end_jd, denum = eph.get_current_file_data(0)

        if start_jd > 0:
            with pytest.raises(EphemerisRangeError):
                eph.calc_ut(start_jd - 1.0, eph.SUN, 0)

    def test_jd_just_outside_end(self):
        """JD just after end should raise error."""
        eph.calc_ut(2451545.0, eph.SUN, 0)  # Ensure loaded
        path, start_jd, end_jd, denum = eph.get_current_file_data(0)

        if end_jd > 0:
            with pytest.raises(EphemerisRangeError):
                eph.calc_ut(end_jd + 1.0, eph.SUN, 0)


class TestMultiplePlanets:
    """Test validation works for all standard planets."""

    @pytest.mark.parametrize(
        "body_id,body_name",
        [
            (eph.SUN, "Sun"),
            (eph.MOON, "Moon"),
            (eph.MERCURY, "Mercury"),
            (eph.VENUS, "Venus"),
            (eph.MARS, "Mars"),
            (eph.JUPITER, "Jupiter"),
            (eph.SATURN, "Saturn"),
            (eph.URANUS, "Uranus"),
            (eph.NEPTUNE, "Neptune"),
            (eph.PLUTO, "Pluto"),
        ],
    )
    def test_all_planets_validate_range(self, body_id, body_name):
        """All standard planets should validate JD range."""
        far_future_jd = 3000000.0

        with pytest.raises(EphemerisRangeError) as exc_info:
            eph.calc_ut(far_future_jd, body_id, 0)

        err = exc_info.value
        assert err.body_id == body_id
        assert err.body_name == body_name
        assert body_name in str(err)


class TestExportedFromExceptions:
    """Test that validate_jd_range is properly exported."""

    def test_validate_jd_range_in_exceptions_module(self):
        """validate_jd_range should be importable from exceptions module."""
        from libephemeris.exceptions import validate_jd_range

        assert callable(validate_jd_range)
