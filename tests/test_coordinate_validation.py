"""Tests for geographic coordinate validation.

This module tests the validation of latitude and longitude coordinates
throughout the libephemeris library.
"""

import pytest
import libephemeris as ephem
from libephemeris.exceptions import (
    CoordinateError,
    validate_latitude,
    validate_longitude,
    validate_coordinates,
)


class TestCoordinateError:
    """Tests for CoordinateError exception class."""

    def test_coordinate_error_attributes(self):
        """CoordinateError should store all relevant attributes."""
        error = CoordinateError(
            message="test message",
            coordinate_name="latitude",
            value=91.0,
            min_value=-90.0,
            max_value=90.0,
        )
        assert str(error) == "test message"
        assert error.coordinate_name == "latitude"
        assert error.value == 91.0
        assert error.min_value == -90.0
        assert error.max_value == 90.0

    def test_coordinate_error_repr(self):
        """CoordinateError should have useful repr."""
        error = CoordinateError(
            message="test",
            coordinate_name="longitude",
            value=200.0,
            min_value=-180.0,
            max_value=180.0,
        )
        repr_str = repr(error)
        assert "CoordinateError" in repr_str
        assert "longitude" in repr_str
        assert "200.0" in repr_str

    def test_coordinate_error_is_subclass_of_error(self):
        """CoordinateError should be a subclass of libephemeris.Error."""
        assert issubclass(CoordinateError, ephem.Error)


class TestValidateLatitude:
    """Tests for validate_latitude function."""

    def test_valid_latitudes(self):
        """Valid latitudes should not raise exceptions."""
        valid_values = [0.0, 45.0, -45.0, 90.0, -90.0, 89.999, -89.999]
        for lat in valid_values:
            validate_latitude(lat)  # Should not raise

    def test_invalid_latitude_too_high(self):
        """Latitude above 90 should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_latitude(91.0)
        error = exc_info.value
        assert error.coordinate_name == "latitude"
        assert error.value == 91.0
        assert error.min_value == -90.0
        assert error.max_value == 90.0
        assert "91" in str(error)
        assert "-90" in str(error) and "90" in str(error)

    def test_invalid_latitude_too_low(self):
        """Latitude below -90 should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_latitude(-91.0)
        error = exc_info.value
        assert error.coordinate_name == "latitude"
        assert error.value == -91.0

    def test_latitude_with_func_name(self):
        """Error message should include function name when provided."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_latitude(100.0, "my_function")
        assert "my_function" in str(exc_info.value)


class TestValidateLongitude:
    """Tests for validate_longitude function."""

    def test_valid_longitudes(self):
        """Valid longitudes should not raise exceptions."""
        valid_values = [0.0, 90.0, -90.0, 180.0, -180.0, 179.999, -179.999]
        for lon in valid_values:
            validate_longitude(lon)  # Should not raise

    def test_valid_east_longitudes(self):
        """East-positive longitudes in (180, 360] should not raise."""
        for lon in [181.0, 200.0, 270.0, 359.999, 360.0]:
            validate_longitude(lon)  # Should not raise

    def test_invalid_longitude_too_high(self):
        """Longitude above 360 should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_longitude(361.0)
        error = exc_info.value
        assert error.coordinate_name == "longitude"
        assert error.value == 361.0
        assert error.min_value == -180.0
        assert error.max_value == 360.0

    def test_invalid_longitude_too_low(self):
        """Longitude below -180 should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_longitude(-181.0)
        error = exc_info.value
        assert error.coordinate_name == "longitude"
        assert error.value == -181.0

    def test_longitude_with_func_name(self):
        """Error message should include function name when provided."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_longitude(400.0, "test_func")
        assert "test_func" in str(exc_info.value)


class TestValidateCoordinates:
    """Tests for validate_coordinates function."""

    def test_valid_coordinates(self):
        """Valid coordinate pairs should not raise exceptions."""
        valid_pairs = [
            (0.0, 0.0),
            (41.9, 12.5),  # Rome
            (-33.9, 18.4),  # Cape Town
            (90.0, 180.0),  # Edge case
            (-90.0, -180.0),  # Edge case
        ]
        for lat, lon in valid_pairs:
            validate_coordinates(lat, lon)  # Should not raise

    def test_invalid_latitude_in_coordinates(self):
        """Invalid latitude should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_coordinates(100.0, 0.0)
        assert exc_info.value.coordinate_name == "latitude"

    def test_invalid_longitude_in_coordinates(self):
        """Invalid longitude should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_coordinates(0.0, 400.0)
        assert exc_info.value.coordinate_name == "longitude"

    def test_latitude_checked_before_longitude(self):
        """Latitude should be validated before longitude."""
        # Both are invalid, but latitude error should be raised first
        with pytest.raises(CoordinateError) as exc_info:
            validate_coordinates(100.0, 400.0)
        assert exc_info.value.coordinate_name == "latitude"


class TestSetTopoValidation:
    """Tests for coordinate validation in set_topo()."""

    def test_set_topo_valid_coordinates(self):
        """set_topo with valid coordinates should work."""
        # Rome
        ephem.set_topo(12.5, 41.9, 0)  # Should not raise

    def test_set_topo_invalid_latitude(self):
        """set_topo with invalid latitude should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            ephem.set_topo(0.0, 91.0, 0)
        error = exc_info.value
        assert error.coordinate_name == "latitude"
        assert "set_topo" in str(error)

    def test_set_topo_invalid_longitude(self):
        """set_topo with invalid longitude should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            ephem.set_topo(400.0, 45.0, 0)
        error = exc_info.value
        assert error.coordinate_name == "longitude"
        assert "set_topo" in str(error)

    def test_set_topo_boundary_values(self):
        """set_topo should accept boundary values."""
        ephem.set_topo(180.0, 90.0, 0)  # Should not raise
        ephem.set_topo(-180.0, -90.0, 0)  # Should not raise

    def test_set_topo_east_longitude(self):
        """set_topo should accept east-positive longitude in (180, 360]."""
        ephem.set_topo(200.0, 45.0, 0)  # Should not raise (= -160 deg)


class TestSweHousesValidation:
    """Tests for coordinate validation in houses()."""

    def test_swe_houses_valid_coordinates(self):
        """houses with valid coordinates should work."""
        jd = 2451545.0  # J2000
        cusps, ascmc = ephem.houses(jd, 41.9, 12.5, ord("P"))
        assert len(cusps) == 12
        assert len(ascmc) == 8

    def test_swe_houses_invalid_latitude(self):
        """houses with invalid latitude should raise CoordinateError."""
        jd = 2451545.0
        with pytest.raises(CoordinateError) as exc_info:
            ephem.houses(jd, 91.0, 0.0, ord("P"))
        error = exc_info.value
        assert error.coordinate_name == "latitude"
        assert "houses" in str(error)

    def test_swe_houses_invalid_longitude(self):
        """houses with invalid longitude should raise CoordinateError."""
        jd = 2451545.0
        with pytest.raises(CoordinateError) as exc_info:
            ephem.houses(jd, 45.0, 400.0, ord("P"))
        error = exc_info.value
        assert error.coordinate_name == "longitude"
        assert "houses" in str(error)

    def test_swe_houses_east_longitude(self):
        """houses accepts east-positive longitude (0-360) like the reference.

        200 deg E is the same meridian as -160 deg, so the cusps must match
        those computed with the signed longitude.
        """
        jd = 2451545.0
        cusps_east, ascmc_east = ephem.houses(jd, 45.0, 200.0, ord("P"))
        cusps_signed, ascmc_signed = ephem.houses(jd, 45.0, -160.0, ord("P"))
        assert len(cusps_east) == 12
        for ce, cs in zip(cusps_east, cusps_signed, strict=True):
            assert abs(ce - cs) < 1e-6, (ce, cs)
        for ae, a_s in zip(ascmc_east, ascmc_signed, strict=True):
            assert abs(ae - a_s) < 1e-6, (ae, a_s)


class TestSweHousesArmcValidation:
    """Tests for coordinate validation in houses_armc()."""

    def test_swe_houses_armc_valid_latitude(self):
        """houses_armc with valid latitude should work."""
        armc = 45.0
        lat = 41.9
        eps = 23.44
        cusps, ascmc = ephem.houses_armc(armc, lat, eps, ord("P"))
        assert len(cusps) == 12

    def test_swe_houses_armc_invalid_latitude(self):
        """houses_armc with invalid latitude should raise CoordinateError."""
        with pytest.raises(CoordinateError) as exc_info:
            ephem.houses_armc(45.0, 91.0, 23.44, ord("P"))
        error = exc_info.value
        assert error.coordinate_name == "latitude"
        assert "houses_armc" in str(error)


class TestContextSetTopoValidation:
    """Tests for coordinate validation in EphemerisContext.set_topo()."""

    def test_context_set_topo_valid_coordinates(self):
        """Context set_topo with valid coordinates should work."""
        ctx = ephem.EphemerisContext()
        ctx.set_topo(12.5, 41.9, 0)  # Should not raise
        assert ctx.topo is not None

    def test_context_set_topo_invalid_latitude(self):
        """Context set_topo with invalid latitude should raise CoordinateError."""
        ctx = ephem.EphemerisContext()
        with pytest.raises(CoordinateError) as exc_info:
            ctx.set_topo(0.0, 95.0, 0)
        assert exc_info.value.coordinate_name == "latitude"

    def test_context_set_topo_invalid_longitude(self):
        """Context set_topo with invalid longitude should raise CoordinateError."""
        ctx = ephem.EphemerisContext()
        with pytest.raises(CoordinateError) as exc_info:
            ctx.set_topo(400.0, 45.0, 0)
        assert exc_info.value.coordinate_name == "longitude"


class TestErrorMessageClarity:
    """Tests for error message clarity and usefulness."""

    def test_latitude_error_message_includes_value(self):
        """Error message should include the invalid value."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_latitude(123.456)
        assert "123.456" in str(exc_info.value)

    def test_latitude_error_message_includes_range(self):
        """Error message should include the valid range."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_latitude(100.0)
        msg = str(exc_info.value)
        assert "-90" in msg
        assert "90" in msg

    def test_longitude_error_message_includes_range(self):
        """Error message should include the valid range."""
        with pytest.raises(CoordinateError) as exc_info:
            validate_longitude(400.0)
        msg = str(exc_info.value)
        assert "-180" in msg
        assert "360" in msg
