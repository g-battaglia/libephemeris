from typing import Any, Dict, Optional, Tuple
import numpy as np

from .constants import ASCENDANT, DESCENDANT, IC, ANTIVERTEX, MC, VERTEX
from .houses import houses_with_fallback
def calc_angles(jd_ut: float, lat: float, lon: float) -> Dict[str, float]: ...
def get_angle_value(angle_id: int, jd_ut: float, lat: float, lon: float) -> float: ...
