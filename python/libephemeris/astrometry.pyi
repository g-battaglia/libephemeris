from typing import Any, Dict, Optional, Tuple
import numpy as np

C_LIGHT_AU_DAY: Any
J2000: Any
JD_PER_CENTURY: Any
TWO_PI: Any
DEG_TO_RAD: Any
RAD_TO_DEG: Any
ARCSEC_TO_RAD: Any
OBLIQUITY_J2000_ARCSEC: Any
OBLIQUITY_COEFFS: Any
def nutation_angles(jd_tt: float) -> Tuple[float, float]: ...
ZETA_A_COEFFS: Any
Z_A_COEFFS: Any
THETA_A_COEFFS: Any
def precess_from_j2000(lon: float, lat: float, jd_tt: float) -> Tuple[float, float]: ...
def apply_aberration_to_position(position: Tuple[float, float, float], earth_velocity: Tuple[float, float, float]) -> Tuple[float, float, float]: ...
