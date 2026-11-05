"""Native numeric utilities for cache."""

from .precession_vondrak import vondrak_mean_obliquity_deg
from ._core import _nf_get_cached_nutation as get_cached_nutation
from ._core import _nf_get_nutation_degrees as get_nutation_degrees
from ._core import _nf_get_cached_obliquity as get_cached_obliquity
from ._core import _nf_get_mean_obliquity as get_mean_obliquity
from ._core import _nf_get_true_obliquity as get_true_obliquity
from ._core import _nf_get_cached_time_ut1 as get_cached_time_ut1
from ._core import _nf_get_cached_time_tt as get_cached_time_tt
from ._core import _nf_clear_caches as clear_caches
from ._core import _nf_get_cache_info as get_cache_info
