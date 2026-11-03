from .constants import INTP_APOG, INTP_PERG
from .exceptions import CalculationError, EphemerisRangeError

StateVector = tuple[float, float, float, float, float, float]
CartesianVector = tuple[float, float, float]
SUPPORTED_BODIES: tuple[int, int]
def state_to_polar(state: StateVector) -> StateVector: ...
