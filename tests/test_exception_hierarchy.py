"""
Tests for the improved exception hierarchy.

This module tests:
1. Exception class inheritance relationships
2. Exception attributes and initialization
3. Category-based exception catching
4. Backward compatibility with existing code
5. Export from libephemeris module
"""

import pytest
import libephemeris as ephem
from libephemeris.exceptions import (
    # Base
    Error,
    # Input validation category
    InputValidationError,
    CoordinateError,
    InvalidBodyError,
    # Data not found category
    DataNotFoundError,
    UnknownBodyError,
    StarNotFoundError,
    SPKNotFoundError,
    # Calculation category
    CalculationError,
    PolarCircleError,
    EphemerisRangeError,
    ConvergenceError,
    # Configuration category
    ConfigurationError,
)


class TestExceptionHierarchyStructure:
    """Test the inheritance structure of the exception hierarchy."""

    @pytest.mark.unit
    def test_error_is_base_of_all(self):
        """All custom exceptions should inherit from Error."""
        assert issubclass(InputValidationError, Error)
        assert issubclass(DataNotFoundError, Error)
        assert issubclass(CalculationError, Error)
        assert issubclass(ConfigurationError, Error)

    @pytest.mark.unit
    def test_error_is_exception(self):
        """Error should inherit from Exception (reference-API compatibility)."""
        assert issubclass(Error, Exception)

    @pytest.mark.unit
    def test_input_validation_category(self):
        """Input validation exceptions inherit from InputValidationError."""
        assert issubclass(CoordinateError, InputValidationError)
        assert issubclass(InvalidBodyError, InputValidationError)
        # And also from Error
        assert issubclass(CoordinateError, Error)
        assert issubclass(InvalidBodyError, Error)

    @pytest.mark.unit
    def test_data_not_found_category(self):
        """Data not found exceptions inherit from DataNotFoundError."""
        assert issubclass(UnknownBodyError, DataNotFoundError)
        assert issubclass(StarNotFoundError, DataNotFoundError)
        assert issubclass(SPKNotFoundError, DataNotFoundError)
        # And also from Error
        assert issubclass(UnknownBodyError, Error)
        assert issubclass(StarNotFoundError, Error)
        assert issubclass(SPKNotFoundError, Error)

    @pytest.mark.unit
    def test_calculation_category(self):
        """Calculation exceptions inherit from CalculationError."""
        assert issubclass(PolarCircleError, CalculationError)
        assert issubclass(EphemerisRangeError, CalculationError)
        assert issubclass(ConvergenceError, CalculationError)
        # And also from Error
        assert issubclass(PolarCircleError, Error)
        assert issubclass(EphemerisRangeError, Error)
        assert issubclass(ConvergenceError, Error)

    @pytest.mark.unit
    def test_configuration_category(self):
        """ConfigurationError inherits from Error."""
        assert issubclass(ConfigurationError, Error)


class TestCategoryBasedCatching:
    """Test catching exceptions by category."""

    @pytest.mark.unit
    def test_catch_input_validation_category(self):
        """InputValidationError catches all input validation errors."""
        # CoordinateError
        with pytest.raises(InputValidationError):
            raise CoordinateError("test")

        # InvalidBodyError
        with pytest.raises(InputValidationError):
            raise InvalidBodyError("test")

    @pytest.mark.unit
    def test_catch_data_not_found_category(self):
        """DataNotFoundError catches all data not found errors."""
        # UnknownBodyError
        with pytest.raises(DataNotFoundError):
            raise UnknownBodyError("test")

        # StarNotFoundError
        with pytest.raises(DataNotFoundError):
            raise StarNotFoundError("test")

        # SPKNotFoundError
        with pytest.raises(DataNotFoundError):
            raise SPKNotFoundError("test")

    @pytest.mark.unit
    def test_catch_calculation_category(self):
        """CalculationError catches all calculation errors."""
        # PolarCircleError
        with pytest.raises(CalculationError):
            raise PolarCircleError("test")

        # EphemerisRangeError
        with pytest.raises(CalculationError):
            raise EphemerisRangeError("test")

        # ConvergenceError
        with pytest.raises(CalculationError):
            raise ConvergenceError("test")

    @pytest.mark.unit
    def test_catch_all_with_error(self):
        """Error catches all library exceptions."""
        exceptions_to_test = [
            InputValidationError("test"),
            CoordinateError("test"),
            InvalidBodyError("test"),
            DataNotFoundError("test"),
            UnknownBodyError("test"),
            StarNotFoundError("test"),
            SPKNotFoundError("test"),
            CalculationError("test"),
            PolarCircleError("test"),
            EphemerisRangeError("test"),
            ConvergenceError("test"),
            ConfigurationError("test"),
        ]
        for exc in exceptions_to_test:
            with pytest.raises(Error):
                raise exc


class TestNewExceptionClasses:
    """Test the new exception classes added in this improvement."""

    @pytest.mark.unit
    def test_invalid_body_error_attributes(self):
        """InvalidBodyError stores all expected attributes."""
        err = InvalidBodyError(
            message="Sun cannot be used for heliacal rising",
            body_id=0,
            body_name="Sun",
            operation="heliacal_rising",
        )
        assert err.message == "Sun cannot be used for heliacal rising"
        assert err.body_id == 0
        assert err.body_name == "Sun"
        assert err.operation == "heliacal_rising"
        assert str(err) == "Sun cannot be used for heliacal rising"

    @pytest.mark.unit
    def test_invalid_body_error_repr(self):
        """InvalidBodyError has informative repr."""
        err = InvalidBodyError("test", body_id=1, body_name="Moon")
        r = repr(err)
        assert "InvalidBodyError" in r
        assert "body_id=1" in r
        assert "Moon" in r

    @pytest.mark.unit
    def test_star_not_found_error_attributes(self):
        """StarNotFoundError stores all expected attributes."""
        err = StarNotFoundError(
            message="Star 'NonexistentStar' not found in catalog",
            star_id="NonexistentStar",
            search_type="name",
        )
        assert err.star_id == "NonexistentStar"
        assert err.search_type == "name"
        assert str(err) == "Star 'NonexistentStar' not found in catalog"

    @pytest.mark.unit
    def test_star_not_found_error_repr(self):
        """StarNotFoundError has informative repr."""
        err = StarNotFoundError("test", star_id="Sirius", search_type="hip")
        r = repr(err)
        assert "StarNotFoundError" in r
        assert "Sirius" in r
        assert "hip" in r

    @pytest.mark.unit
    def test_convergence_error_attributes(self):
        """ConvergenceError stores all expected attributes."""
        err = ConvergenceError(
            message="Brent's method failed to converge",
            algorithm="brent",
            iterations=100,
            max_iterations=100,
            tolerance=1e-10,
            last_value=0.0001,
        )
        assert err.algorithm == "brent"
        assert err.iterations == 100
        assert err.max_iterations == 100
        assert err.tolerance == 1e-10
        assert err.last_value == 0.0001
        assert str(err) == "Brent's method failed to converge"

    @pytest.mark.unit
    def test_convergence_error_repr(self):
        """ConvergenceError has informative repr."""
        err = ConvergenceError("test", algorithm="newton", iterations=50)
        r = repr(err)
        assert "ConvergenceError" in r
        assert "newton" in r
        assert "iterations=50" in r

    @pytest.mark.unit
    def test_configuration_error_attributes(self):
        """ConfigurationError stores all expected attributes."""
        err = ConfigurationError(
            message="Observer location not set",
            missing_config="observer_location",
            suggestion="Call set_topo(lon, lat, alt) first",
        )
        assert err.missing_config == "observer_location"
        assert err.suggestion == "Call set_topo(lon, lat, alt) first"
        assert str(err) == "Observer location not set"

    @pytest.mark.unit
    def test_configuration_error_repr(self):
        """ConfigurationError has informative repr."""
        err = ConfigurationError(
            "test", missing_config="ephe_path", suggestion="Set path"
        )
        r = repr(err)
        assert "ConfigurationError" in r
        assert "ephe_path" in r


class TestCategoryExceptions:
    """Test the category base exception classes."""

    @pytest.mark.unit
    def test_input_validation_error_basic(self):
        """InputValidationError can be created and raised."""
        err = InputValidationError("Invalid input provided")
        assert str(err) == "Invalid input provided"
        with pytest.raises(InputValidationError):
            raise err

    @pytest.mark.unit
    def test_data_not_found_error_basic(self):
        """DataNotFoundError can be created and raised."""
        err = DataNotFoundError("Data not found")
        assert str(err) == "Data not found"
        with pytest.raises(DataNotFoundError):
            raise err

    @pytest.mark.unit
    def test_calculation_error_basic(self):
        """CalculationError can be created and raised."""
        err = CalculationError("Calculation failed")
        assert str(err) == "Calculation failed"
        with pytest.raises(CalculationError):
            raise err


class TestBackwardCompatibility:
    """Test backward compatibility with existing exception usage."""

    @pytest.mark.unit
    def test_existing_exceptions_still_work(self):
        """Existing exceptions work exactly as before."""
        # CoordinateError with attributes
        coord_err = CoordinateError(
            "Invalid latitude",
            coordinate_name="latitude",
            value=91.0,
            min_value=-90.0,
            max_value=90.0,
        )
        assert coord_err.coordinate_name == "latitude"
        assert coord_err.value == 91.0

        # PolarCircleError with attributes
        polar_err = PolarCircleError(
            "Polar circle error",
            latitude=70.0,
            threshold=66.5,
            obliquity=23.4,
            house_system="P",
        )
        assert polar_err.latitude == 70.0
        assert polar_err.house_system == "P"

        # SPKNotFoundError with factory method
        spk_err = SPKNotFoundError.from_filepath(
            "/path/to/file.bsp", body_name="Chiron", body_id="2060"
        )
        assert spk_err.filepath == "/path/to/file.bsp"
        assert "download_spk" in str(spk_err)

        # UnknownBodyError with body_id
        unknown_err = UnknownBodyError("Unknown body", body_id=99999)
        assert unknown_err.body_id == 99999

        # EphemerisRangeError with attributes
        range_err = EphemerisRangeError(
            "Date out of range",
            requested_jd=5000000.0,
            start_jd=2287184.5,
            end_jd=2688976.5,
        )
        assert range_err.requested_jd == 5000000.0

    @pytest.mark.unit
    def test_catch_by_original_type_still_works(self):
        """Catching by the original exception type still works."""
        with pytest.raises(CoordinateError):
            raise CoordinateError("test")

        with pytest.raises(PolarCircleError):
            raise PolarCircleError("test")

        with pytest.raises(UnknownBodyError):
            raise UnknownBodyError("test")

        with pytest.raises(SPKNotFoundError):
            raise SPKNotFoundError("test")

        with pytest.raises(EphemerisRangeError):
            raise EphemerisRangeError("test")

    @pytest.mark.unit
    def test_catch_by_error_base_still_works(self):
        """Catching by Error base class still works for all."""
        exceptions = [
            CoordinateError("test"),
            PolarCircleError("test"),
            UnknownBodyError("test"),
            SPKNotFoundError("test"),
            EphemerisRangeError("test"),
        ]
        for exc in exceptions:
            with pytest.raises(Error):
                raise exc


class TestModuleExports:
    """Test that exceptions are properly exported from libephemeris."""

    @pytest.mark.unit
    def test_all_exceptions_in_module(self):
        """All exception classes are accessible via libephemeris."""
        # Base
        assert hasattr(ephem, "Error")
        # Category classes
        assert hasattr(ephem, "InputValidationError")
        assert hasattr(ephem, "DataNotFoundError")
        assert hasattr(ephem, "CalculationError")
        assert hasattr(ephem, "ConfigurationError")
        # Specific exceptions
        assert hasattr(ephem, "CoordinateError")
        assert hasattr(ephem, "InvalidBodyError")
        assert hasattr(ephem, "UnknownBodyError")
        assert hasattr(ephem, "StarNotFoundError")
        assert hasattr(ephem, "SPKNotFoundError")
        assert hasattr(ephem, "PolarCircleError")
        assert hasattr(ephem, "EphemerisRangeError")
        assert hasattr(ephem, "ConvergenceError")

    @pytest.mark.unit
    def test_all_exceptions_in_all(self):
        """All exception classes are in __all__."""
        expected = [
            "Error",
            "InputValidationError",
            "DataNotFoundError",
            "CalculationError",
            "ConfigurationError",
            "CoordinateError",
            "InvalidBodyError",
            "UnknownBodyError",
            "StarNotFoundError",
            "SPKNotFoundError",
            "PolarCircleError",
            "EphemerisRangeError",
            "ConvergenceError",
        ]
        for name in expected:
            assert name in ephem.__all__, f"{name} not in __all__"

    @pytest.mark.unit
    def test_exception_identity(self):
        """Exceptions from module match those from exceptions module."""
        assert ephem.Error is Error
        assert ephem.InputValidationError is InputValidationError
        assert ephem.DataNotFoundError is DataNotFoundError
        assert ephem.CalculationError is CalculationError
        assert ephem.ConfigurationError is ConfigurationError
        assert ephem.CoordinateError is CoordinateError
        assert ephem.InvalidBodyError is InvalidBodyError
        assert ephem.UnknownBodyError is UnknownBodyError
        assert ephem.StarNotFoundError is StarNotFoundError
        assert ephem.SPKNotFoundError is SPKNotFoundError
        assert ephem.PolarCircleError is PolarCircleError
        assert ephem.EphemerisRangeError is EphemerisRangeError
        assert ephem.ConvergenceError is ConvergenceError


class TestExceptionDocstrings:
    """Test that all exceptions have proper docstrings."""

    @pytest.mark.unit
    def test_all_exceptions_have_docstrings(self):
        """All exception classes have non-empty docstrings."""
        exceptions = [
            Error,
            InputValidationError,
            DataNotFoundError,
            CalculationError,
            ConfigurationError,
            CoordinateError,
            InvalidBodyError,
            UnknownBodyError,
            StarNotFoundError,
            SPKNotFoundError,
            PolarCircleError,
            EphemerisRangeError,
            ConvergenceError,
        ]
        for exc in exceptions:
            assert exc.__doc__ is not None, f"{exc.__name__} has no docstring"
            assert len(exc.__doc__) > 50, f"{exc.__name__} docstring too short"


class TestExceptionUsagePatterns:
    """Test typical usage patterns for the exception hierarchy."""

    @pytest.mark.unit
    def test_try_except_with_category_fallback(self):
        """Demonstrate category-based error handling."""

        def mock_calculation():
            raise PolarCircleError("Cannot calculate at polar latitude")

        result = None
        try:
            mock_calculation()
        except CalculationError as e:
            # Can catch any calculation error
            result = "calculation_failed"

        assert result == "calculation_failed"

    @pytest.mark.unit
    def test_multiple_category_handling(self):
        """Handle different error categories differently."""

        def handle_errors(exc):
            try:
                raise exc
            except InputValidationError:
                return "validation_error"
            except DataNotFoundError:
                return "not_found_error"
            except CalculationError:
                return "calculation_error"
            except ConfigurationError:
                return "config_error"
            except Error:
                return "unknown_error"

        assert handle_errors(CoordinateError("bad coords")) == "validation_error"
        assert handle_errors(InvalidBodyError("bad body")) == "validation_error"
        assert handle_errors(UnknownBodyError("unknown")) == "not_found_error"
        assert handle_errors(StarNotFoundError("no star")) == "not_found_error"
        assert handle_errors(SPKNotFoundError("no spk")) == "not_found_error"
        assert handle_errors(PolarCircleError("polar")) == "calculation_error"
        assert handle_errors(EphemerisRangeError("range")) == "calculation_error"
        assert handle_errors(ConvergenceError("no converge")) == "calculation_error"
        assert handle_errors(ConfigurationError("no config")) == "config_error"
