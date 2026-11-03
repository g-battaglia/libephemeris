from typing import Tuple, List
from typing import Any

from .._core import (
    _module_moon_theories_tass17_TASS_EPOCH_JD as TASS_EPOCH_JD,
    _module_moon_theories_tass17_AU_KM as AU_KM,
    _module_moon_theories_tass17_SATELLITE_NAMES as SATELLITE_NAMES,
)
def saturn_moon_position(jd: float, body: int=5) -> Tuple[float, float, float]: ...
def saturn_moon_position_velocity(jd: float, body: int=5) -> Tuple[float, float, float, float, float, float]: ...
