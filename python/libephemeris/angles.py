"""Native numeric utilities for angles."""

from .constants import ASCENDANT, DESCENDANT, IC, ANTIVERTEX, MC, VERTEX
from .houses import houses_with_fallback
from ._core import _nf_angles_calc_angles as calc_angles
from ._core import _nf_get_angle_value as get_angle_value
