from typing import Any, Dict, Optional, Tuple
import numpy as np

from .constants import EARTH, JUPITER, MARS, MERCURY, NEPTUNE, SATURN, URANUS, VENUS
def mean_orbital_elements(ipl: int, jd_tt: float) -> tuple[float, float, float, float, float] | None: ...
