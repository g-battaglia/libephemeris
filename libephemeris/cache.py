# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
Caching utilities for libephemeris hot path optimization.

This module provides LRU caches for expensive calculations that are called
repeatedly with the same inputs. Key optimizations include:

1. Nutation caching - IAU 2006/2000A nutation calculations (via pyerfa)
2. Obliquity caching - True obliquity of the ecliptic (IAU 2006)
3. Timescale object caching - Skyfield timescale objects

These caches significantly improve performance when calculating multiple
planetary positions for the same Julian Day (common in chart calculations).

Precision:
    - Nutation: IAU 2006/2000A model via pyerfa (~0.01-0.05 mas precision)
    - Obliquity: IAU 2006 model via pyerfa (consistent across all code paths)

Thread Safety:
    The caches use functools.lru_cache which is thread-safe for reads but
    not for simultaneous writes with the same key. In practice, this is
    acceptable since we're caching pure functions with deterministic outputs.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Tuple

import erfa


# Cache size limits
# - NUTATION_CACHE_SIZE: Caches nutation angles (dpsi, deps) for Julian Days
#   Most chart calculations use 1-2 JDs, but time ranges may use more
# - OBLIQUITY_CACHE_SIZE: Caches obliquity values, same usage pattern
# - TIME_CACHE_SIZE: Caches Skyfield Time objects per JD to avoid recreating
#   nutation/precession matrices for repeated calls at the same JD
_NUTATION_CACHE_SIZE = 16
_OBLIQUITY_CACHE_SIZE = 16
_TIME_CACHE_SIZE = 8
_OBSERVER_CACHE_MAX = 32

# J2000.0 epoch in Julian Days
_J2000_JD = 2451545.0


@lru_cache(maxsize=_NUTATION_CACHE_SIZE)
def get_cached_nutation(jd_tt: float) -> Tuple[float, float]:
    """
    Get cached nutation angles (dpsi, deps) for a Julian Day.

    Uses the IAU 2006/2000A nutation model via pyerfa (erfa.nut06a), which
    is the most precise nutation model currently adopted by the IAU.
    This includes the IAU 2000A lunisolar+planetary nutation (1365 terms)
    with IAU 2006 J2 secular variation correction.

    Precision: ~0.01-0.05 milliarcsecond

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (dpsi_radians, deps_radians):
            - dpsi: Nutation in longitude (radians)
            - deps: Nutation in obliquity (radians)

    Performance:
        - Uncached: ~0.02ms per call
        - Cached hit: ~0.0001ms per call (200x faster)
    """
    # erfa.nut06a expects (jd1, jd2) where jd1 + jd2 = JD in TT
    # Using 2-part JD for maximum numerical precision (SOFA convention)
    dpsi, deps = erfa.nut06a(_J2000_JD, jd_tt - _J2000_JD)
    return float(dpsi), float(deps)


def get_nutation_degrees(jd_tt: float) -> Tuple[float, float]:
    """
    Get nutation in degrees (convenience wrapper).

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (dpsi_degrees, deps_degrees)
    """
    dpsi_rad, deps_rad = get_cached_nutation(jd_tt)
    return math.degrees(dpsi_rad), math.degrees(deps_rad)


@lru_cache(maxsize=_OBLIQUITY_CACHE_SIZE)
def get_cached_obliquity(jd_tt: float) -> Tuple[float, float]:
    """
    Get cached mean and true obliquity for a Julian Day.

    Uses IAU 2006 mean obliquity (via erfa.obl06) and IAU 2006/2000A
    nutation (via erfa.nut06a) for true obliquity, providing consistency
    with all other code paths in libephemeris.

    Mean obliquity formula (IAU 2006, Capitaine et al. 2003):
        eps0 = 84381.406″ - 46.836769″T - 0.0001831″T² + 0.00200340″T³
               - 0.000000576″T⁴ - 0.0000000434″T⁵

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (mean_obliquity, true_obliquity) in degrees
    """
    # Mean obliquity via IAU 2006 (erfa.obl06)
    # Returns mean obliquity in radians
    eps0_rad = erfa.obl06(_J2000_JD, jd_tt - _J2000_JD)
    eps0 = math.degrees(eps0_rad)

    # Nutation in obliquity from IAU 2006/2000A
    _, deps_rad = get_cached_nutation(jd_tt)
    deps = math.degrees(deps_rad)

    # True obliquity = mean + nutation in obliquity
    eps = eps0 + deps

    return eps0, eps


def get_mean_obliquity(jd_tt: float) -> float:
    """Get mean obliquity in degrees (IAU 2006)."""
    eps0, _ = get_cached_obliquity(jd_tt)
    return eps0


def get_true_obliquity(jd_tt: float) -> float:
    """Get true obliquity in degrees (IAU 2006 + nutation)."""
    _, eps = get_cached_obliquity(jd_tt)
    return eps


@lru_cache(maxsize=_TIME_CACHE_SIZE)
def get_cached_time_ut1(jd: float):
    """
    Get a cached Skyfield Time object for a UT1 Julian Day.

    Avoids recreating Time objects (and their lazy nutation/precession
    matrices) when the same JD is used repeatedly, e.g. when computing
    positions for multiple planets in a single chart.

    Args:
        jd: Julian Day in UT1

    Returns:
        Skyfield Time object
    """
    from .state import get_timescale

    return get_timescale().ut1_jd(jd)


@lru_cache(maxsize=_TIME_CACHE_SIZE)
def get_cached_time_tt(jd: float):
    """
    Get a cached Skyfield Time object for a TT Julian Day.

    Args:
        jd: Julian Day in Terrestrial Time (TT)

    Returns:
        Skyfield Time object
    """
    from .state import get_timescale

    return get_timescale().tt_jd(jd)


# Observer position cache: stores observer.at(t) results keyed on
# (id(observer), jd_tt). This avoids recomputing Earth's position
# for every planet in a multi-planet chart at the same JD.
#
# To guard against Python reusing memory addresses after deallocation
# (which would cause cache hits returning positions computed for a
# *different* observer), we store the observer object alongside the
# result and validate identity on lookup.
_observer_at_cache: dict = {}


def get_cached_observer_at(observer, t):
    """
    Get a cached result of observer.at(t).

    For geocentric calculations, Earth's position at a given JD is
    identical for all planets. This cache prevents ~10x redundant
    evaluations per chart.

    Args:
        observer: Skyfield VectorSum (e.g. earth, sun)
        t: Skyfield Time object

    Returns:
        Skyfield position object (Geocentric, Barycentric, etc.)
    """
    key = (id(observer), float(t.tt))
    entry = _observer_at_cache.get(key)
    # Validate that the cached observer is the *same* object, not just
    # a new object that was allocated at the same memory address.
    if entry is not None:
        cached_observer, result = entry
        if cached_observer is observer:
            return result
    if len(_observer_at_cache) > _OBSERVER_CACHE_MAX:
        _observer_at_cache.clear()
    result = observer.at(t)
    _observer_at_cache[key] = (observer, result)
    return result


def clear_observer_cache():
    """Clear the observer-at-time cache.

    Called when a Skyfield Time reify descriptor corruption is detected,
    forcing fresh position computations on the next call.
    """
    _observer_at_cache.clear()


def clear_caches() -> None:
    """
    Clear all computation caches.

    Call this when:
    - Changing ephemeris files
    - After calling close()
    - When memory needs to be freed
    """
    get_cached_nutation.cache_clear()
    get_cached_obliquity.cache_clear()
    get_cached_time_ut1.cache_clear()
    get_cached_time_tt.cache_clear()
    _observer_at_cache.clear()

    # Clear fast_calc frame-data caches (imported lazily to avoid circular deps)
    try:
        from .fast_calc import (
            _get_skyfield_frame_data,
            _get_precession_matrix,
            _leb_frame_cache,
        )

        _get_skyfield_frame_data.cache_clear()
        _get_precession_matrix.cache_clear()
        # LEB frame data is keyed by jd_tt but depends on the reader's
        # nutation table; clear on close to avoid serving stale values
        # if the active ephemeris is swapped.
        _leb_frame_cache.clear()
    except ImportError:
        pass

    # Clear the refraction atmosphere-profile cache. The cached profile is
    # keyed by (altitude, pressure, temperature, lapse_rate); clearing on
    # close ensures a clean baseline for the next session.
    try:
        from . import refraction as _refraction

        _refraction._refr_cache_entry = None
    except ImportError:
        pass


def get_cache_info() -> dict:
    """
    Get cache statistics for debugging and optimization.

    Returns:
        Dictionary with cache statistics for each cached function
    """
    return {
        "nutation": get_cached_nutation.cache_info()._asdict(),
        "obliquity": get_cached_obliquity.cache_info()._asdict(),
        "time_ut1": get_cached_time_ut1.cache_info()._asdict(),
        "time_tt": get_cached_time_tt.cache_info()._asdict(),
        "observer_at": {
            "size": len(_observer_at_cache),
            "maxsize": _OBSERVER_CACHE_MAX,
        },
    }
