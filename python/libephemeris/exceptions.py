"""Public exception classes created by the native extension."""

from ._core import (
    CalculationError,
    ConfigurationError,
    ConvergenceError,
    CoordinateError,
    DataNotFoundError,
    EphemerisRangeError,
    Error,
    IllegalBodyError,
    InputValidationError,
    InvalidBodyError,
    LEBCorruptionError,
    FictitiousRuntimeDispatch,
    NetworkSealedError,
    PolarCircleError,
    SPKEvaluationError,
    SPKNotFoundError,
    SPKRequiredError,
    StarNotFoundError,
    UnknownBodyError,
    validate_coordinates,
    validate_latitude,
    validate_longitude,
)

__all__ = [
    "Error", "InputValidationError", "CoordinateError", "InvalidBodyError",
    "DataNotFoundError", "UnknownBodyError", "IllegalBodyError", "StarNotFoundError",
    "SPKNotFoundError", "SPKRequiredError", "CalculationError", "PolarCircleError",
    "EphemerisRangeError", "ConvergenceError", "SPKEvaluationError",
    "ConfigurationError", "NetworkSealedError",
    "LEBCorruptionError", "FictitiousRuntimeDispatch",
    "validate_coordinates", "validate_latitude", "validate_longitude", "validate_jd_range",
]

from ._core import validate_jd_range
