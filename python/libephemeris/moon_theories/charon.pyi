from typing import Tuple, List
from typing import Any

from .._core import (
    _module_moon_theories_charon_CHARON_A_KM as CHARON_A_KM,
    _module_moon_theories_charon_CHARON_PERIOD_DAYS as CHARON_PERIOD_DAYS,
    _module_moon_theories_charon_CHARON_N_RAD_DAY as CHARON_N_RAD_DAY,
    _module_moon_theories_charon_CHARON_E as CHARON_E,
    _module_moon_theories_charon_EPOCH_JD as EPOCH_JD,
    _module_moon_theories_charon_M0_RAD as M0_RAD,
    _module_moon_theories_charon_OMEGA_RAD as OMEGA_RAD,
    _module_moon_theories_charon_I_RAD as I_RAD,
    _module_moon_theories_charon_NODE_RAD as NODE_RAD,
    _module_moon_theories_charon_NODE_RATE_RAD_DAY as NODE_RATE_RAD_DAY,
    _module_moon_theories_charon_AU_KM as AU_KM,
)
def charon_position(jd: float) -> Tuple[float, float, float]: ...
