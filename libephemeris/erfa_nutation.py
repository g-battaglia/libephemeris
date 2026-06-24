# SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-LibEphemeris-Commercial
# Copyright (c) 2025-2026 Giacomo Battaglia
"""
ERFA Nutation Integration Module.

This module provides pyerfa-based nutation functions as an optional enhancement
to the Skyfield-based IAU 2000A/B implementations.

Key functions:
    - erfa.nut00a(): IAU 2000A nutation (equivalent to Skyfield's iau2000a_radians)
    - erfa.nut06a(): IAU 2000A nutation with IAU 2006 precession adjustments
    - erfa.pnm06a(): Combined bias-precession-nutation matrix (IAU 2006/2000A)

Precision comparison:
    - IAU 2000B (Skyfield default): ~1 milliarcsecond precision (77 terms)
    - IAU 2000A (Skyfield high-precision): ~0.1 milliarcsecond precision (1365 terms)
    - IAU 2006/2000A (ERFA nut06a): ~0.01-0.05 milliarcsecond precision
      (IAU 2000A + J2 secular variation + obliquity frame correction)

References:
    - IERS Conventions 2010, Chapter 5
    - Capitaine, N. & Wallace, P.T., 2006, Astron.Astrophys. 450, 855
    - Mathews, P.M., Herring, T.A., Buffet, B.A. 2002, J.Geophys.Res.
"""

from __future__ import annotations

import math
from functools import lru_cache
from typing import Tuple, Optional

# Try to import pyerfa
try:
    import erfa

    _HAS_ERFA = True
except ImportError:
    _HAS_ERFA = False


def has_erfa() -> bool:
    """Check if pyerfa is available."""
    return _HAS_ERFA


def get_erfa_nutation_nut00a(jd_tt: float) -> Optional[Tuple[float, float]]:
    """
    Get nutation angles using ERFA's nut00a (IAU 2000A model).

    This is the same underlying model as Skyfield's iau2000a_radians but
    computed directly by the ERFA library.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (dpsi_radians, deps_radians) if pyerfa available, None otherwise

    Notes:
        - Uses IAU 2000A nutation model (MHB2000)
        - 1365 terms for luni-solar and planetary nutation
        - Precision: ~100 microarcseconds (0.1 milliarcseconds)
        - Free core nutation (FCN) is omitted
    """
    if not _HAS_ERFA:
        return None

    # ERFA uses 2-part Julian date for precision
    # Using J2000 method: date1 = 2451545.0, date2 = offset
    dpsi, deps = erfa.nut00a(2451545.0, jd_tt - 2451545.0)
    return float(dpsi), float(deps)


def get_erfa_nutation_nut06a(jd_tt: float) -> Optional[Tuple[float, float]]:
    """
    Get nutation angles using ERFA's nut06a (IAU 2006/2000A model).

    This is the most accurate nutation model available, combining:
    - IAU 2000A nutation (1365 terms)
    - Corrections for IAU 2006 precession frame
    - J2 secular variation correction

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Tuple of (dpsi_radians, deps_radians) if pyerfa available, None otherwise

    Notes:
        - Precision: few tens of microarcseconds (~0.01-0.05 mas)
        - The most rigorous nutation model currently available
        - Recommended for highest-precision applications
    """
    if not _HAS_ERFA:
        return None

    # Using J2000 method for optimal internal precision
    dpsi, deps = erfa.nut06a(2451545.0, jd_tt - 2451545.0)
    return float(dpsi), float(deps)


def get_erfa_pnm06a_matrix(jd_tt: float) -> Optional[Tuple]:
    """
    Get the bias-precession-nutation matrix using ERFA's pnm06a.

    This returns a 3x3 rotation matrix that combines:
    - Frame bias (GCRS to mean J2000)
    - IAU 2006 precession
    - IAU 2000A nutation (with IAU 2006 adjustments)

    The matrix transforms coordinates from GCRS to true equator and equinox
    of date in a single operation.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        3x3 numpy array if pyerfa available, None otherwise.
        The matrix operates as: V(date) = rbpn @ V(GCRS)

    Notes:
        - This is more rigorous than applying bias, precession, and nutation
          as separate rotations, avoiding cross-term errors.
        - After a century, the separate rotation approach introduces errors
          exceeding 1 mas; this matrix avoids that.
    """
    if not _HAS_ERFA:
        return None

    rbpn = erfa.pnm06a(2451545.0, jd_tt - 2451545.0)
    return rbpn


@lru_cache(maxsize=128)
def get_erfa_nutation_cached(
    jd_tt: float, model: str = "nut06a"
) -> Tuple[float, float]:
    """
    Get cached ERFA nutation angles.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)
        model: Either "nut00a" (IAU 2000A) or "nut06a" (IAU 2006/2000A)

    Returns:
        Tuple of (dpsi_radians, deps_radians)

    Raises:
        ImportError: If pyerfa is not available
        ValueError: If model is not recognized

    Notes:
        Falls back to Skyfield's iau2000a_radians if pyerfa unavailable.
    """
    if not _HAS_ERFA:
        # Fallback to Skyfield
        from skyfield.nutationlib import iau2000a_radians

        from .state import get_timescale

        ts = get_timescale()
        t = ts.tt_jd(jd_tt)
        dpsi_rad, deps_rad = iau2000a_radians(t)
        # Skyfield returns numpy.float64; coerce to native floats to match
        # the erfa paths and the library-wide "always return float" invariant.
        return float(dpsi_rad), float(deps_rad)

    if model == "nut00a":
        result = get_erfa_nutation_nut00a(jd_tt)
    elif model == "nut06a":
        result = get_erfa_nutation_nut06a(jd_tt)
    else:
        raise ValueError(f"Unknown nutation model: {model}. Use 'nut00a' or 'nut06a'.")

    if result is None:
        raise ImportError("pyerfa is required for ERFA nutation models")

    return result


def get_erfa_obliquity_iau2006(jd_tt: float) -> Optional[float]:
    """
    Get mean obliquity using IAU 2006 precession model (erfa.obl06).

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Mean obliquity in radians if pyerfa available, None otherwise

    Notes:
        This is slightly different from the Laskar 1986 formula used
        as fallback. The difference is typically < 0.001 arcseconds
        for dates within a few centuries of J2000.
    """
    if not _HAS_ERFA:
        return None

    eps0 = erfa.obl06(2451545.0, jd_tt - 2451545.0)
    return float(eps0)


def compare_nutation_models(jd_tt: float) -> dict:
    """
    Compare nutation values from different models for a given date.

    This is primarily for validation and research purposes.

    Args:
        jd_tt: Julian Day in Terrestrial Time (TT)

    Returns:
        Dictionary with nutation values from each model and their differences
    """
    from skyfield.nutationlib import iau2000a_radians, iau2000b_radians

    from .state import get_timescale

    ts = get_timescale()
    t = ts.tt_jd(jd_tt)

    # Skyfield models
    dpsi_2000b, deps_2000b = iau2000b_radians(t)
    dpsi_2000a, deps_2000a = iau2000a_radians(t)

    result = {
        "jd_tt": jd_tt,
        "skyfield_iau2000b": {
            "dpsi_arcsec": math.degrees(dpsi_2000b) * 3600,
            "deps_arcsec": math.degrees(deps_2000b) * 3600,
        },
        "skyfield_iau2000a": {
            "dpsi_arcsec": math.degrees(dpsi_2000a) * 3600,
            "deps_arcsec": math.degrees(deps_2000a) * 3600,
        },
    }

    # ERFA models (if available)
    if _HAS_ERFA:
        erfa_00a = get_erfa_nutation_nut00a(jd_tt)
        erfa_06a = get_erfa_nutation_nut06a(jd_tt)

        if erfa_00a:
            result["erfa_nut00a"] = {
                "dpsi_arcsec": math.degrees(erfa_00a[0]) * 3600,
                "deps_arcsec": math.degrees(erfa_00a[1]) * 3600,
            }

        if erfa_06a:
            result["erfa_nut06a"] = {
                "dpsi_arcsec": math.degrees(erfa_06a[0]) * 3600,
                "deps_arcsec": math.degrees(erfa_06a[1]) * 3600,
            }

        # Compute differences (in milliarcseconds)
        if erfa_00a and erfa_06a:
            result["differences_mas"] = {
                "skyfield_2000a_vs_erfa_00a_dpsi": (
                    math.degrees(dpsi_2000a - erfa_00a[0]) * 3600 * 1000
                ),
                "skyfield_2000a_vs_erfa_00a_deps": (
                    math.degrees(deps_2000a - erfa_00a[1]) * 3600 * 1000
                ),
                "erfa_00a_vs_erfa_06a_dpsi": (
                    math.degrees(erfa_00a[0] - erfa_06a[0]) * 3600 * 1000
                ),
                "erfa_00a_vs_erfa_06a_deps": (
                    math.degrees(erfa_00a[1] - erfa_06a[1]) * 3600 * 1000
                ),
                "skyfield_2000b_vs_2000a_dpsi": (
                    math.degrees(dpsi_2000b - dpsi_2000a) * 3600 * 1000
                ),
                "skyfield_2000b_vs_2000a_deps": (
                    math.degrees(deps_2000b - deps_2000a) * 3600 * 1000
                ),
            }

    return result


def clear_erfa_nutation_cache() -> None:
    """Clear the ERFA nutation cache."""
    get_erfa_nutation_cached.cache_clear()


def get_erfa_nutation_cache_info() -> dict:
    """Get ERFA nutation cache statistics."""
    return get_erfa_nutation_cached.cache_info()._asdict()
