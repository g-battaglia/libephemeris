class Error(Exception): ...
def validate_latitude(lat: float, func_name: str = '') -> None: ...
def validate_longitude(lon: float, func_name: str = '') -> None: ...
def validate_coordinates(lat: float, lon: float, func_name: str = '') -> None: ...
class InputValidationError(Error): ...
class CoordinateError(InputValidationError):
    message: str
    coordinate_name: str | None
    value: float | None
    min_value: float | None
    max_value: float | None
    def __init__(self, message: str, coordinate_name: str | None = None, value: float | None = None, min_value: float | None = None, max_value: float | None = None) -> None: ...
class InvalidBodyError(InputValidationError):
    message: str
    body_id: int | None
    body_name: str | None
    operation: str | None
    def __init__(self, message: str, body_id: int | None = None, body_name: str | None = None, operation: str | None = None) -> None: ...
class DataNotFoundError(Error): ...
class UnknownBodyError(DataNotFoundError):
    message: str
    body_id: int | None
    def __init__(self, message: str, body_id: int | None = None) -> None: ...
class IllegalBodyError(UnknownBodyError, ValueError): ...
class StarNotFoundError(DataNotFoundError):
    message: str
    star_id: str | None
    search_type: str | None
    def __init__(self, message: str, star_id: str | None = None, search_type: str | None = None) -> None: ...
class SPKNotFoundError(DataNotFoundError):
    message: str
    filepath: str | None
    body_name: str | None
    body_id: str | None
    def __init__(self, message: str, filepath: str | None = None, body_name: str | None = None, body_id: str | None = None) -> None: ...
    @classmethod
    def from_filepath(cls, filepath: str, body_name: str | None = None, body_id: str | None = None) -> SPKNotFoundError: ...
class SPKRequiredError(DataNotFoundError):
    message: str
    body_id: int | None
    body_name: str | None
    def __init__(self, message: str, body_id: int | None = None, body_name: str | None = None) -> None: ...
    @classmethod
    def for_body(cls, body_id: int, body_name: str, horizons_id: str) -> SPKRequiredError: ...
class CalculationError(Error): ...
class PolarCircleError(CalculationError):
    message: str
    latitude: float | None
    threshold: float | None
    obliquity: float | None
    house_system: str | None
    def __init__(self, message: str, latitude: float | None = None, threshold: float | None = None, obliquity: float | None = None, house_system: str | None = None) -> None: ...
class EphemerisRangeError(CalculationError):
    message: str
    requested_jd: float | None
    start_jd: float | None
    end_jd: float | None
    start_date: str | None
    end_date: str | None
    body_id: int | None
    body_name: str | None
    ephemeris_file: str | None
    def __init__(self, message: str, requested_jd: float | None = None, start_jd: float | None = None, end_jd: float | None = None, start_date: str | None = None, end_date: str | None = None, body_id: int | None = None, body_name: str | None = None, ephemeris_file: str | None = None) -> None: ...
class ConvergenceError(CalculationError):
    message: str
    algorithm: str | None
    iterations: int | None
    max_iterations: int | None
    tolerance: float | None
    last_value: float | None
    def __init__(self, message: str, algorithm: str | None = None, iterations: int | None = None, max_iterations: int | None = None, tolerance: float | None = None, last_value: float | None = None) -> None: ...
class SPKEvaluationError(CalculationError):
    message: str
    body_id: int | None
    naif_id: int | None
    spk_file: str | None
    requested_jd: float | None
    def __init__(self, message: str, body_id: int | None = None, naif_id: int | None = None, spk_file: str | None = None, requested_jd: float | None = None) -> None: ...
class ConfigurationError(Error):
    message: str
    missing_config: str | None
    suggestion: str | None
    def __init__(self, message: str, missing_config: str | None = None, suggestion: str | None = None) -> None: ...
class NetworkSealedError(ConfigurationError, RuntimeError):
    purpose: str | None
    def __init__(self, message: str, purpose: str | None = None) -> None: ...

class LEBCorruptionError(ValueError): ...
class FictitiousRuntimeDispatch(KeyError): ...

def validate_jd_range(jd: float, body_id: int | None=None, func_name: str="") -> None: ...
