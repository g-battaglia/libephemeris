"""Native scalar lunar-state conversion; use Session for interpolated apsides."""

from .constants import INTP_APOG, INTP_PERG
from .exceptions import CalculationError, EphemerisRangeError
from ._core import _nf_state_to_polar as state_to_polar
from ._core import _nf_SUPPORTED_LUNAR_APSIDES as SUPPORTED_BODIES
