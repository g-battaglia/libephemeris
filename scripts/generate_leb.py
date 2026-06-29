#!/usr/bin/env python3
"""
LEB file generator.

Produces a .leb binary ephemeris file by sampling Skyfield/libephemeris
analytical functions and fitting Chebyshev polynomials to the results.

Usage:
    python scripts/generate_leb.py --output ephemeris.leb
    python scripts/generate_leb.py --tier medium
    python scripts/generate_leb.py --tier extended --start -5000 --end 5000
    python scripts/generate_leb.py --start 1550 --end 2650 --output de440_full.leb
    python scripts/generate_leb.py --output ephemeris.leb --verify --verify-samples 500
    python scripts/generate_leb.py --output ephemeris.leb --workers 8
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import struct
import subprocess
import sys
import time
from typing import Callable, List, Optional, Tuple

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval


# =============================================================================
# PROGRESS BAR
# =============================================================================


class ProgressBar:
    """Lightweight terminal progress bar (no dependencies).

    Usage::

        bar = ProgressBar(total=100, label="Sun", width=30)
        for i in range(100):
            do_work()
            bar.update(i + 1)
        bar.finish()
    """

    def __init__(
        self,
        total: int,
        label: str = "",
        width: int = 30,
        unit: str = "seg",
        enabled: bool = True,
    ):
        self.total = max(total, 1)
        self.label = label
        self.width = width
        self.unit = unit
        self.enabled = enabled
        self.t0 = time.monotonic()
        self._last_draw = 0.0
        # Determine terminal width once
        try:
            self._term_cols = shutil.get_terminal_size().columns
        except Exception:
            self._term_cols = 80

    def update(self, current: int) -> None:
        if not self.enabled:
            return
        now = time.monotonic()
        # Throttle redraws to every 100ms (unless it's the last update)
        if current < self.total and now - self._last_draw < 0.1:
            return
        self._last_draw = now
        self._draw(current, now)

    def _draw(self, current: int, now: float) -> None:
        elapsed = now - self.t0
        frac = current / self.total
        filled = int(self.width * frac)
        bar = "\u2588" * filled + "\u2591" * (self.width - filled)
        pct = f"{frac * 100:5.1f}%"

        # ETA
        if current > 0 and frac < 1.0:
            eta = elapsed / frac * (1.0 - frac)
            time_str = f"ETA {_fmt_time(eta)}"
        else:
            time_str = f"{_fmt_time(elapsed)}"

        counter = f"{current}/{self.total} {self.unit}"
        line = f"\r  {self.label:20s} {bar} {pct} {counter:>14s} {time_str:>10s}"
        # Truncate to terminal width
        if len(line) > self._term_cols:
            line = line[: self._term_cols]
        sys.stdout.write(line)
        sys.stdout.flush()

    def finish(self, suffix: str = "") -> None:
        if not self.enabled:
            return
        self._draw(self.total, time.monotonic())
        if suffix:
            sys.stdout.write(f"  {suffix}")
        sys.stdout.write("\n")
        sys.stdout.flush()


def _fmt_time(seconds: float) -> str:
    """Format seconds as mm:ss or hh:mm:ss."""
    s = int(seconds)
    if s < 3600:
        return f"{s // 60}:{s % 60:02d}"
    h = s // 3600
    m = (s % 3600) // 60
    return f"{h}:{m:02d}:{s % 60:02d}"


# Add parent directory to path so we can import libephemeris
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

from libephemeris.leb_format import (
    BODY_ENTRY_SIZE,
    BODY_PARAMS,
    COORD_ECLIPTIC,
    COORD_HELIO_ECL,
    COORD_ICRS_BARY,
    COORD_ICRS_BARY_SYSTEM,
    DELTA_T_ENTRY_FMT,
    DELTA_T_ENTRY_SIZE,
    DELTA_T_HEADER_FMT,
    DELTA_T_HEADER_SIZE,
    HEADER_SIZE,
    MAGIC,
    NUTATION_HEADER_SIZE,
    SECTION_BODY_INDEX,
    SECTION_CHEBYSHEV,
    SECTION_DELTA_T,
    SECTION_DIR_SIZE,
    SECTION_NUTATION,
    SECTION_STARS,
    STAR_ENTRY_SIZE,
    VERSION,
    BodyEntry,
    FileHeader,
    NutationHeader,
    SectionEntry,
    StarEntry,
    segment_byte_size,
    write_body_entry,
    write_header,
    write_nutation_header,
    write_section_dir,
    write_star_entry,
)
from libephemeris.exotic_bodies import (
    EXOTIC_IDS,
    naif_map as _exotic_naif,
    name_map as _exotic_names,
)

# =============================================================================
# CONSTANTS
# =============================================================================

# Julian Day for J2000.0
J2000 = 2451545.0

# Default date range (200 years)
DEFAULT_START_YEAR = 1900
DEFAULT_END_YEAR = 2100

# Tier configurations: (ephemeris_file, default_start_year, default_end_year, output_name)
TIER_CONFIGS = {
    "base": ("de440s.bsp", 1850, 2150, "ephemeris_base.leb"),
    "medium": ("de440.bsp", 1550, 2650, "ephemeris_medium.leb"),
    "extended": ("de441.bsp", -5000, 5000, "ephemeris_extended.leb"),
}

# Default output directory for tier-based generation
DEFAULT_LEB_DIR = os.path.join(os.path.dirname(__file__), "..", "data", "leb")

# Nutation Chebyshev parameters
NUTATION_INTERVAL = 16.0  # days
NUTATION_DEGREE = 16
NUTATION_COMPONENTS = 2  # dpsi, deps

# Delta-T sampling interval
DELTA_T_INTERVAL = 30.0  # days

# Number of sections in the file
NUM_SECTIONS = 5  # body_index, chebyshev, nutation, delta_t, stars

# Body names for logging
BODY_NAMES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
    10: "Mean Node",
    11: "True Node",
    12: "Mean Apogee",
    13: "Oscu Apogee",
    14: "Earth",
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
    21: "Interp Apogee",
    22: "Interp Perigee",
    40: "Cupido",
    41: "Hades",
    42: "Zeus",
    43: "Kronos",
    44: "Apollon",
    45: "Admetos",
    46: "Vulkanus",
    47: "Poseidon",
    48: "Transpluto",
}
# Exotic minor bodies (centaurs/TNOs/NEAs) — labels from the registry.
BODY_NAMES.update(_exotic_names())

# Reverse lookup: case-insensitive name → body ID.
# Supports both canonical names ("Mean Node") and short aliases ("moon").
_NAME_TO_BODY_ID: dict[str, int] = {}
for _bid, _bname in BODY_NAMES.items():
    _NAME_TO_BODY_ID[_bname.lower()] = _bid
    # Also register without spaces (e.g. "meannode" → 10)
    _no_space = _bname.lower().replace(" ", "")
    if _no_space != _bname.lower():
        _NAME_TO_BODY_ID[_no_space] = _bid


def _resolve_body_token(token: str) -> int:
    """Resolve a body token (ID or name) to a numeric body ID.

    Accepts:
      - Numeric IDs: "1", "14"
      - Names (case-insensitive): "moon", "Moon", "MOON"
      - Names without spaces: "meannode", "oscuapogee", "interpapogee"

    Raises:
        ValueError: If the token cannot be resolved.
    """
    token = token.strip()
    # Try numeric first
    try:
        return int(token)
    except ValueError:
        pass
    # Try name lookup
    key = token.lower()
    if key in _NAME_TO_BODY_ID:
        return _NAME_TO_BODY_ID[key]
    raise ValueError(
        f"Unknown body '{token}'. Use a numeric ID or one of: "
        + ", ".join(sorted(_NAME_TO_BODY_ID.keys()))
    )


# Planet name map for Skyfield (body_id -> skyfield name)
_PLANET_MAP = {
    0: "sun",
    1: "moon",
    2: "mercury",
    3: "venus",
    4: "mars",
    5: "jupiter",
    6: "saturn",
    7: "uranus",
    8: "neptune",
    9: "pluto",
    14: "earth",
}

# Outer planet barycenter names for system-barycenter generation
_SYSTEM_BARY_MAP = {
    5: "jupiter barycenter",
    6: "saturn barycenter",
    7: "uranus barycenter",
    8: "neptune barycenter",
    9: "pluto barycenter",
}

# Asteroid NAIF IDs for Skyfield SPK lookup
_ASTEROID_NAIF = {
    15: 2060,  # Chiron
    17: 2000001,  # Ceres
    18: 2000002,  # Pallas
    19: 2000003,  # Juno
    20: 2000004,  # Vesta
}
# Exotic minor bodies (centaurs/TNOs/NEAs) — registry is the single source of
# truth. Membership in _ASTEROID_NAIF is the gate that routes a body through
# the SPK-Chebyshev generation path (download, per-body coverage clamp, verify).
_ASTEROID_NAIF.update(_exotic_naif())

# Body groups for independent generation + merge workflow.
# Each group can be generated as a standalone .leb file and later merged.
BODY_GROUPS: dict[str, List[int]] = {
    "planets": sorted(_PLANET_MAP.keys()),  # 11 planets (ICRS_BARY, vectorized)
    # Classic asteroids = SPK bodies minus the exotic registry (stays correct if
    # a 6th classic is ever added to _ASTEROID_NAIF without touching EXOTIC_IDS).
    "asteroids": sorted(set(_ASTEROID_NAIF) - set(EXOTIC_IDS)),
    "exotics": list(EXOTIC_IDS),  # centaurs/TNOs/NEAs (ICRS_BARY, spktype21)
    "analytical": sorted(
        bid
        for bid in BODY_PARAMS
        if bid not in _PLANET_MAP and bid not in _ASTEROID_NAIF
    ),  # ecliptic/helio analytical bodies
}


# =============================================================================
# CHEBYSHEV FITTING
# =============================================================================


def chebyshev_nodes(n: int) -> np.ndarray:
    """Compute n Chebyshev nodes (Type I) on [-1, 1]."""
    return np.cos(np.pi * (np.arange(n) + 0.5) / n)


def fit_segment(
    func: Callable[[float], np.ndarray],
    jd_start: float,
    jd_end: float,
    degree: int,
    components: int,
) -> np.ndarray:
    """Fit a multi-component function over [jd_start, jd_end] with Chebyshev polynomials.

    Args:
        func: Function taking JD and returning array of shape (components,).
        jd_start: Start of segment.
        jd_end: End of segment.
        degree: Polynomial degree.
        components: Number of output components.

    Returns:
        Coefficients array of shape (components, degree+1), stored as
        [c0_comp0, c1_comp0, ..., cN_comp0, c0_comp1, ...].
    """
    nodes = chebyshev_nodes(degree + 1)
    # Map nodes from [-1, 1] to [jd_start, jd_end]
    jd_nodes = 0.5 * (jd_end - jd_start) * nodes + 0.5 * (jd_start + jd_end)

    # Evaluate function at nodes
    values = np.array([func(jd) for jd in jd_nodes])  # shape: (degree+1, components)

    # Fit each component independently
    coeffs = np.zeros((components, degree + 1))
    for c in range(components):
        coeffs[c] = chebfit(nodes, values[:, c], degree)

    return coeffs


def verify_segment(
    func: Callable[[float], np.ndarray],
    coeffs: np.ndarray,
    seg_start: float,
    seg_end: float,
    components: int,
    n_test: int = 10,
    verify_end: float | None = None,
) -> float:
    """Verify a Chebyshev fit by evaluating at intermediate points.

    Args:
        seg_start: Segment start JD (defines the polynomial domain).
        seg_end: Segment end JD (defines the polynomial domain).
        verify_end: If given, only sample verification points up to this JD
            (for the last segment which may extend beyond the data range).
            Tau is always computed relative to [seg_start, seg_end] since
            that is the domain the polynomial was fitted on.

    Returns the maximum error across all components and test points.
    """
    if verify_end is None:
        verify_end = seg_end
    mid = 0.5 * (seg_start + seg_end)
    half = 0.5 * (seg_end - seg_start)
    max_error = 0.0
    for i in range(n_test):
        # Uniform test points (not Chebyshev nodes)
        frac = (i + 0.5) / n_test
        jd = seg_start + frac * (verify_end - seg_start)
        tau = (jd - mid) / half

        reference = func(jd)
        for c in range(components):
            fitted = chebval(tau, coeffs[c])
            error = abs(fitted - reference[c])
            if error > max_error:
                max_error = error

    return max_error


# =============================================================================
# VECTORIZED CHEBYSHEV FITTING
# =============================================================================

N_VERIFY = 10  # Number of verification points per segment


def _compute_all_segment_jds(
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    n_verify: int = N_VERIFY,
) -> Tuple[np.ndarray, int, int]:
    """Precompute all JDs (fit nodes + verification points) for all segments.

    Returns:
        (all_jds, n_segments, pts_per_segment)
        First (degree+1) JDs per segment are Chebyshev nodes.
        Next n_verify JDs per segment are uniform verification points.
    """
    n_segments = int(math.ceil((jd_end - jd_start) / interval_days))
    nodes_01 = chebyshev_nodes(degree + 1)  # in [-1, 1]
    pts_per_seg = degree + 1 + n_verify

    all_jds = np.empty(n_segments * pts_per_seg)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        seg_end = seg_start + interval_days
        mid = 0.5 * (seg_start + seg_end)
        half = 0.5 * (seg_end - seg_start)
        offset = i * pts_per_seg

        # Chebyshev nodes mapped to [seg_start, seg_end]
        all_jds[offset : offset + degree + 1] = half * nodes_01 + mid

        # Uniform verification points (not on Chebyshev nodes)
        verify_end = min(seg_end, jd_end)
        for v in range(n_verify):
            frac = (v + 0.5) / n_verify
            all_jds[offset + degree + 1 + v] = seg_start + frac * (
                verify_end - seg_start
            )

    return all_jds, n_segments, pts_per_seg


def _fit_and_verify_from_values(
    all_values: np.ndarray,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    components: int,
    n_segments: int,
    pts_per_seg: int,
    n_verify: int = N_VERIFY,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Fit Chebyshev segments and verify from pre-evaluated values.

    Args:
        all_values: shape (n_total, components) with pre-evaluated function values.
            Layout: [seg0_node0, ..., seg0_nodeN, seg0_verify0, ..., seg0_verifyM,
                     seg1_node0, ...]

    Returns:
        (list_of_coefficient_arrays, max_error)
    """
    nodes_01 = chebyshev_nodes(degree + 1)
    all_coeffs: List[np.ndarray] = []
    max_error = 0.0
    bar = ProgressBar(n_segments, label=label or "fit+verify", enabled=verbose)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        seg_end = seg_start + interval_days
        offset = i * pts_per_seg

        # Extract fitting values (Chebyshev nodes)
        fit_values = all_values[offset : offset + degree + 1]

        # Fit each component
        coeffs = np.zeros((components, degree + 1))
        for c in range(components):
            coeffs[c] = chebfit(nodes_01, fit_values[:, c], degree)

        # Verify using pre-computed verification points
        verify_end = min(seg_end, jd_end)
        mid = 0.5 * (seg_start + seg_end)
        half = 0.5 * (seg_end - seg_start)

        for v in range(n_verify):
            frac = (v + 0.5) / n_verify
            jd_v = seg_start + frac * (verify_end - seg_start)
            tau = (jd_v - mid) / half

            ref = all_values[offset + degree + 1 + v]
            for c in range(components):
                fitted = float(chebval(tau, coeffs[c]))
                error = abs(fitted - ref[c])
                if error > max_error:
                    max_error = error

        all_coeffs.append(coeffs)
        bar.update(i + 1)

    bar.finish()
    return all_coeffs, max_error


def _fit_and_verify_from_values_unwrap(
    all_values: np.ndarray,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    components: int,
    n_segments: int,
    pts_per_seg: int,
    n_verify: int = N_VERIFY,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Fit Chebyshev segments with longitude unwrapping from pre-evaluated values.

    Same as _fit_and_verify_from_values but unwraps component 0 (longitude)
    before fitting and re-wraps during verification. Used for ecliptic and
    heliocentric bodies.

    Args:
        all_values: shape (n_total, components) with pre-evaluated function values.
            Layout: [seg0_node0, ..., seg0_nodeN, seg0_verify0, ..., seg0_verifyM,
                     seg1_node0, ...]

    Returns:
        (list_of_coefficient_arrays, max_error)
    """
    nodes_01 = chebyshev_nodes(degree + 1)
    all_coeffs: List[np.ndarray] = []
    max_error = 0.0
    bar = ProgressBar(n_segments, label=label or "fit+verify", enabled=verbose)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        seg_end = seg_start + interval_days
        offset = i * pts_per_seg

        # Extract fitting values (Chebyshev nodes)
        fit_values = all_values[offset : offset + degree + 1].copy()

        # Unwrap longitude (component 0) to remove 360-degree jumps
        fit_values[:, 0] = np.unwrap(np.radians(fit_values[:, 0]))
        fit_values[:, 0] = np.degrees(fit_values[:, 0])

        # Fit each component
        coeffs = np.zeros((components, degree + 1))
        for c in range(components):
            coeffs[c] = chebfit(nodes_01, fit_values[:, c], degree)

        # Verify using pre-computed verification points (with re-wrapping)
        verify_end = min(seg_end, jd_end)
        mid = 0.5 * (seg_start + seg_end)
        half = 0.5 * (seg_end - seg_start)

        for v in range(n_verify):
            frac = (v + 0.5) / n_verify
            jd_v = seg_start + frac * (verify_end - seg_start)
            tau = (jd_v - mid) / half

            ref = all_values[offset + degree + 1 + v]
            for c in range(components):
                fitted = float(chebval(tau, coeffs[c]))
                if c == 0:
                    # Re-wrap longitude for comparison
                    fitted = fitted % 360.0
                    ref_val = ref[c] % 360.0
                    error = abs(fitted - ref_val)
                    if error > 180.0:
                        error = 360.0 - error
                else:
                    error = abs(fitted - ref[c])
                if error > max_error:
                    max_error = error

        all_coeffs.append(coeffs)
        bar.update(i + 1)

    bar.finish()
    return all_coeffs, max_error


# =============================================================================
# VECTORIZED ECLIPTIC BODY BATCH EVALUATION
# =============================================================================


def _calc_mean_lilith_batch(all_jds: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Vectorized version of calc_mean_lilith / _calc_mean_apse_analytical.

    Args:
        all_jds: Array of Julian Days (TT).
        T: Array of Julian centuries from J2000 (pre-computed).

    Returns:
        Array of mean apogee longitudes in degrees [0, 360).
    """
    T2 = T * T
    fracT = T % 1.0

    z_F_T2 = -1.312045233711e01
    z_F_T3 = -1.138215912580e-03
    z_F_T4 = -9.646018347184e-06
    z_MP_T2 = 3.146734198839e01
    z_MP_T3 = 4.768357585780e-02
    z_MP_T4 = -3.421689790404e-04
    z_LP_T2 = -5.663161722088e00
    z_LP_T3 = 5.722859298199e-03
    z_LP_T4 = -8.466472828815e-05

    NF = 1739232000.0 * fracT + 295263.0983 * T - 0.2079419901760 * T + 335779.55755
    NF = NF % 1296000.0
    NF += ((z_F_T4 * T + z_F_T3) * T + z_F_T2) * T2

    MP = 1717200000.0 * fracT + 715923.4728 * T - 0.2035946368532 * T + 485868.28096
    MP = MP % 1296000.0
    MP += ((z_MP_T4 * T + z_MP_T3) * T + z_MP_T2) * T2

    LP = 1731456000.0 * fracT + 1108372.83264 * T - 0.6784914260953 * T + 785939.95571
    LP = LP % 1296000.0
    LP += ((z_LP_T4 * T + z_LP_T3) * T + z_LP_T2) * T2

    STR = np.pi / (180.0 * 3600.0)

    apogee_rad = (LP - MP) * STR + np.pi
    apogee_rad = apogee_rad % (2.0 * np.pi)

    node_rad = (LP - NF) * STR
    node_rad = node_rad % (2.0 * np.pi)

    MOON_MEAN_INCL = 5.1453964
    lon_from_node = apogee_rad - node_rad

    x = np.cos(lon_from_node)
    y = np.sin(lon_from_node)

    incl_rad = -np.radians(MOON_MEAN_INCL)
    cos_incl = np.cos(incl_rad)
    # z=0 so y_new = y * cos_incl, x stays the same
    y_new = y * cos_incl

    lon_from_node_proj = np.arctan2(y_new, x)

    apogee_projected = (lon_from_node_proj + node_rad) % (2.0 * np.pi)
    return np.degrees(apogee_projected)


def _calc_lunar_fundamental_arguments_batch(
    all_jds: np.ndarray, T: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Vectorized version of _calc_lunar_fundamental_arguments.

    Args:
        all_jds: Array of Julian Days (TT).
        T: Array of Julian centuries from J2000.

    Returns:
        (D, M, M_prime, F) arrays in radians.
    """
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    T5 = T4 * T

    D = (
        297.8501921
        + 445267.1114034 * T
        - 0.0018819 * T2
        + T3 / 545868.0
        - T4 / 113065000.0
        + T5 / 18999000000.0
    )

    M = (
        357.5291092
        + 35999.0502909 * T
        - 0.0001536 * T2
        + T3 / 24490000.0
        - T4 / 992300000.0
        + T5 / 189900000000.0
    )

    M_prime = (
        134.9633964
        + 477198.8675055 * T
        + 0.0087414 * T2
        + T3 / 69699.0
        - T4 / 14712000.0
        + T5 / 2520410000.0
    )

    F = (
        93.2720950
        + 483202.0175233 * T
        - 0.0036539 * T2
        - T3 / 3526000.0
        + T4 / 863310000.0
        - T5 / 142650000000.0
    )

    D = np.radians(D % 360.0)
    M = np.radians(M % 360.0)
    M_prime = np.radians(M_prime % 360.0)
    F = np.radians(F % 360.0)

    return D, M, M_prime, F


def _calc_elp2000_apogee_perturbations_batch(
    all_jds: np.ndarray, T: np.ndarray
) -> np.ndarray:
    """Vectorized version of _calc_elp2000_apogee_perturbations.

    Returns array of perturbation corrections in degrees.
    """
    D, M, M_prime, F = _calc_lunar_fundamental_arguments_batch(all_jds, T)

    E = 1.0 - 0.002516 * T - 0.0000074 * T**2
    E2 = E * E

    p = np.zeros_like(T)

    # Primary evection harmonics
    p += 0.1892 * np.sin(D - M_prime)
    p += 4.6921 * np.sin(2.0 * D - 2.0 * M_prime)
    p += -0.0127 * np.sin(3.0 * D - 3.0 * M_prime)
    p += 0.7854 * np.sin(4.0 * D - 4.0 * M_prime)
    p += 0.0089 * np.sin(5.0 * D - 5.0 * M_prime)
    p += 0.1634 * np.sin(6.0 * D - 6.0 * M_prime)
    p += -0.0056 * np.sin(7.0 * D - 7.0 * M_prime)
    p += 0.0412 * np.sin(8.0 * D - 8.0 * M_prime)
    p += 0.0023 * np.sin(9.0 * D - 9.0 * M_prime)
    p += 0.0108 * np.sin(10.0 * D - 10.0 * M_prime)

    # Solar anomaly coupling
    p += 0.3847 * E * np.sin(M)
    p += 0.0198 * E2 * np.sin(2.0 * M)
    p += 0.5123 * E * np.sin(2.0 * D - 2.0 * M_prime - M)
    p += 0.1287 * E * np.sin(2.0 * D - 2.0 * M_prime + M)
    p += -0.0523 * E * np.sin(D - M_prime - M)
    p += 0.0412 * E * np.sin(D - M_prime + M)
    p += 0.0876 * E * np.sin(4.0 * D - 4.0 * M_prime - M)
    p += 0.0234 * E * np.sin(4.0 * D - 4.0 * M_prime + M)
    p += 0.0187 * E * np.sin(6.0 * D - 6.0 * M_prime - M)

    # Double solar anomaly coupling
    p += 0.0312 * E2 * np.sin(2.0 * D - 2.0 * M_prime - 2.0 * M)
    p += 0.0156 * E2 * np.sin(2.0 * D - 2.0 * M_prime + 2.0 * M)

    # Lunar anomaly harmonics
    p += -0.0234 * np.sin(M_prime)
    p += 0.0087 * np.sin(2.0 * M_prime)
    p += -0.0034 * np.sin(3.0 * M_prime)

    # Latitude coupling
    p += 0.2634 * np.sin(2.0 * F - 2.0 * M_prime)
    p += 0.0423 * np.sin(2.0 * F - 2.0 * D)
    p += -0.0289 * np.sin(2.0 * F)
    p += 0.0156 * np.sin(2.0 * F - 4.0 * M_prime + 2.0 * D)
    p += -0.0098 * np.sin(2.0 * F + 2.0 * M_prime - 2.0 * D)

    # Cross-coupling
    p += 0.0723 * np.sin(2.0 * D - M_prime)
    p += -0.0567 * np.sin(2.0 * D - 3.0 * M_prime)
    p += 0.0234 * np.sin(4.0 * D - 3.0 * M_prime)
    p += -0.0178 * np.sin(4.0 * D - 5.0 * M_prime)
    p += 0.0112 * np.sin(2.0 * D)
    p += -0.0089 * np.sin(4.0 * D)

    # Solar-lunar-latitude coupling
    p += 0.0178 * E * np.sin(2.0 * F - 2.0 * M_prime + M)
    p += -0.0134 * E * np.sin(2.0 * F - 2.0 * M_prime - M)
    p += 0.0089 * E * np.sin(2.0 * F - 2.0 * D + M)
    p += -0.0067 * E * np.sin(2.0 * F - 2.0 * D - M)

    # Secular corrections
    p += 0.0012 * T * np.sin(2.0 * D - 2.0 * M_prime)
    p += -0.0008 * T * np.sin(M)

    return p


def _calc_elp2000_perigee_perturbations_batch(
    all_jds: np.ndarray, T: np.ndarray
) -> np.ndarray:
    """Vectorized version of _calc_elp2000_perigee_perturbations.

    Returns array of perturbation corrections in degrees.
    """
    D, M, M_prime, F = _calc_lunar_fundamental_arguments_batch(all_jds, T)

    E = 1.0 - 0.002516 * T - 0.0000074 * T**2
    E2 = E * E

    p = np.zeros_like(T)

    # Polynomial mean perigee corrections
    p += -0.1749
    p += -0.1411 * T
    p += -0.0140 * T * T
    p += +0.0168 * T * T * T

    # Primary evection harmonics sin(kD - kM')
    p += +0.3002 * np.sin(D - M_prime)
    p += -22.2062 * np.sin(2.0 * D - 2.0 * M_prime)
    p += -0.1594 * np.sin(3.0 * D - 3.0 * M_prime)
    p += +6.4536 * np.sin(4.0 * D - 4.0 * M_prime)
    p += +0.0938 * np.sin(5.0 * D - 5.0 * M_prime)
    p += -2.2814 * np.sin(6.0 * D - 6.0 * M_prime)
    p += -0.0375 * np.sin(7.0 * D - 7.0 * M_prime)
    p += +0.4792 * np.sin(8.0 * D - 8.0 * M_prime)
    p += +0.0075 * np.sin(9.0 * D - 9.0 * M_prime)
    p += -0.0598 * np.sin(10.0 * D - 10.0 * M_prime)
    p += +0.0114 * np.sin(12.0 * D - 12.0 * M_prime)
    p += -0.0031 * np.sin(14.0 * D - 14.0 * M_prime)
    p += +0.0011 * np.sin(16.0 * D - 16.0 * M_prime)

    # Evection phase corrections cos(kD - kM')
    p += -0.0750 * np.cos(2.0 * D - 2.0 * M_prime)
    p += -0.0013 * np.cos(3.0 * D - 3.0 * M_prime)
    p += +0.0393 * np.cos(4.0 * D - 4.0 * M_prime)
    p += -0.0061 * np.cos(8.0 * D - 8.0 * M_prime)
    p += +0.0039 * np.cos(10.0 * D - 10.0 * M_prime)
    p += -0.0023 * np.cos(6.0 * D - 6.0 * M_prime)
    p += -0.0011 * np.cos(9.0 * D - 9.0 * M_prime)

    # Solar anomaly coupling
    p += +0.4684 * E * np.sin(M)
    p += -0.9747 * E * np.sin(2.0 * D - 2.0 * M_prime - M)
    p += +0.0935 * E * np.sin(2.0 * D - 2.0 * M_prime + M)
    p += -0.0266 * E * np.sin(D - M_prime - M)
    p += -0.0580 * E * np.sin(D - M_prime + M)
    p += +0.5348 * E * np.sin(4.0 * D - 4.0 * M_prime - M)
    p += -0.0829 * E * np.sin(4.0 * D - 4.0 * M_prime + M)
    p += -0.2059 * E * np.sin(6.0 * D - 6.0 * M_prime - M)
    p += +0.0586 * E * np.sin(6.0 * D - 6.0 * M_prime + M)

    # Solar double coupling (E² terms)
    p += +0.0016 * E2 * np.sin(2.0 * M)
    p += -0.0390 * E2 * np.sin(2.0 * D - 2.0 * M_prime - 2.0 * M)
    p += +0.0707 * E2 * np.sin(2.0 * D - 2.0 * M_prime + 2.0 * M)
    p += +0.0284 * E2 * np.sin(4.0 * D - 4.0 * M_prime - 2.0 * M)

    # Lunar anomaly harmonics
    p += +0.0106 * np.sin(M_prime)
    p += +0.0013 * np.sin(2.0 * M_prime)

    # Latitude coupling
    p += +0.1695 * np.sin(2.0 * F - 2.0 * M_prime)
    p += -0.0539 * np.sin(2.0 * F - 2.0 * D)
    p += -0.0258 * np.sin(2.0 * F - 4.0 * M_prime + 2.0 * D)

    # Cross-coupling
    p += -0.0354 * np.sin(2.0 * D - M_prime)
    p += +0.0039 * np.sin(2.0 * D - 3.0 * M_prime)
    p += +0.1551 * np.sin(4.0 * D - 3.0 * M_prime)
    p += +0.0067 * np.sin(4.0 * D - 5.0 * M_prime)
    p += -0.0024 * np.sin(2.0 * D)
    p += -0.4541 * np.sin(6.0 * D - 5.0 * M_prime)
    p += -0.0010 * np.sin(3.0 * D - 2.0 * M_prime)

    # Solar-latitude cross-coupling
    p += -0.0017 * E * np.sin(2.0 * F - 2.0 * M_prime + M)
    p += -0.0098 * E * np.sin(2.0 * F - 2.0 * M_prime - M)
    p += +0.0095 * E * np.sin(2.0 * F - 2.0 * D - M)

    # Higher-order evection-solar coupling
    p += +0.0376 * E * np.sin(8.0 * D - 8.0 * M_prime - M)
    p += -0.0209 * E * np.sin(8.0 * D - 8.0 * M_prime + M)
    p += -0.0066 * E * np.sin(10.0 * D - 10.0 * M_prime - M)

    # Secular and long-period corrections
    p += +0.0013 * T * np.sin(M)
    p += -0.0014 * T * np.sin(D - M_prime)
    p += -0.0042 * T * np.cos(2.0 * D - 2.0 * M_prime)

    # Cosine phase corrections
    p += +0.0168 * np.cos(M)
    p += +0.0217 * np.cos(2.0 * F - 2.0 * M_prime)

    # Sun-Moon anomaly coupling
    p += -0.0021 * E * np.sin(M - M_prime)
    p += -0.0012 * E * np.sin(2.0 * D + M - M_prime)

    return p


def _eval_ecliptic_bodies_batch(
    all_jds: np.ndarray,
    body_ids: List[int],
    verbose: bool = False,
) -> dict:
    """Evaluate all ecliptic bodies at all JDs in a single vectorized pass.

    Uses ONE Skyfield call for all bodies that need Moon ephemeris data,
    then computes each body's values vectorially with numpy.

    Args:
        all_jds: Array of Julian Days (TT) for evaluation.
        body_ids: List of ecliptic body IDs to evaluate (subset of [10-13, 21, 22]).
        verbose: Print progress.

    Returns:
        dict mapping body_id -> (N, 3) array of [lon, lat, dist].
    """
    N = len(all_jds)
    T = (all_jds - 2451545.0) / 36525.0
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T

    results: dict[int, np.ndarray] = {}

    # Bodies needing Skyfield: 11 (true node), 13 (oscu apogee),
    # 21 (interp apogee lat/dist), 22 (interp perigee lat/dist)
    skyfield_bodies = {11, 13, 21, 22}
    need_skyfield = bool(skyfield_bodies.intersection(body_ids))

    if need_skyfield:
        from skyfield.framelib import ecliptic_frame

        planets, ts = _init_skyfield()
        earth = planets["earth"]
        moon = planets["moon"]

        if verbose:
            print("    Vectorized Skyfield call for Moon ecliptic state...")
        t = ts.tt_jd(all_jds)
        moon_pos = (moon - earth).at(t)
        r = moon_pos.frame_xyz(ecliptic_frame).au  # (3, N)
        _, v_obj = moon_pos.frame_xyz_and_velocity(ecliptic_frame)
        v = v_obj.au_per_d  # (3, N)

        # Angular momentum h = r × v  (component-wise on (3, N) arrays)
        h_x = r[1] * v[2] - r[2] * v[1]
        h_y = r[2] * v[0] - r[0] * v[2]
        h_z = r[0] * v[1] - r[1] * v[0]
        h_mag = np.sqrt(h_x**2 + h_y**2 + h_z**2)

        # |r| for eccentricity vector
        r_mag = np.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2)

        # Gravitational parameter μ for Earth-Moon system in AU³/day²
        gm_earth = 398600.435436  # km³/s²
        earth_moon_mass_ratio = 81.3005691
        gm_moon = gm_earth / earth_moon_mass_ratio
        gm_earth_moon = gm_earth + gm_moon
        mu = gm_earth_moon / (149597870.7**3) * (86400**2)

        # v × h (for eccentricity vector)
        vxh_x = v[1] * h_z - v[2] * h_y
        vxh_y = v[2] * h_x - v[0] * h_z
        vxh_z = v[0] * h_y - v[1] * h_x

        # Eccentricity vector e = (v×h)/μ - r/|r|
        e_x = vxh_x / mu - r[0] / r_mag
        e_y = vxh_y / mu - r[1] / r_mag
        e_z = vxh_z / mu - r[2] / r_mag
        e_mag = np.sqrt(e_x**2 + e_y**2 + e_z**2)

        # Semi-latus rectum
        p = h_mag**2 / mu

    # --- Body 10: Mean Node (pure polynomial) ---
    if 10 in body_ids:
        Omega = (
            125.0445479
            - 1934.1362891 * T
            + 0.0020754 * T2
            + T3 / 467441.0
            - T4 / 60616000.0
        ) % 360.0
        results[10] = np.column_stack([Omega, np.zeros(N), np.zeros(N)])

    # --- Body 11: True Node (from angular momentum) ---
    if 11 in body_ids:
        node_lon = np.degrees(np.arctan2(h_x, -h_y)) % 360.0
        node_dist = h_mag * 1000.0
        results[11] = np.column_stack([node_lon, np.zeros(N), node_dist])

    # --- Body 12: Mean Lilith with latitude ---
    if 12 in body_ids:
        mean_lilith = _calc_mean_lilith_batch(all_jds, T)
        # Mean node (same formula as body 10)
        mean_node = (
            125.0445479
            - 1934.1362891 * T
            + 0.0020754 * T2
            + T3 / 467441.0
            - T4 / 60616000.0
        ) % 360.0
        omega = (mean_lilith - mean_node) % 360.0
        mean_lilith_lat = 5.145 * np.sin(np.radians(omega))
        results[12] = np.column_stack([mean_lilith, mean_lilith_lat, np.zeros(N)])
    else:
        # Still need mean_lilith for bodies 21 and 22
        mean_lilith = None

    # --- Body 13: Oscu Apogee / True Lilith (from eccentricity vector) ---
    if 13 in body_ids:
        apogee_lon = np.degrees(np.arctan2(-e_y, -e_x)) % 360.0
        apogee_lat = np.degrees(np.arcsin(-e_z / e_mag))
        apogee_dist = p / (1.0 - e_mag)
        results[13] = np.column_stack([apogee_lon, apogee_lat, apogee_dist])

    # --- Body 21: Interpolated Apogee ---
    if 21 in body_ids:
        if mean_lilith is None:
            mean_lilith = _calc_mean_lilith_batch(all_jds, T)
        apogee_pert = _calc_elp2000_apogee_perturbations_batch(all_jds, T)
        interp_apogee_lon = (mean_lilith + apogee_pert) % 360.0
        # Lat/dist from osculating apogee (True Lilith)
        ia_lat = np.degrees(np.arcsin(-e_z / e_mag))
        ia_dist = p / (1.0 - e_mag)
        results[21] = np.column_stack([interp_apogee_lon, ia_lat, ia_dist])

    # --- Body 22: Interpolated Perigee ---
    if 22 in body_ids:
        if mean_lilith is None:
            mean_lilith = _calc_mean_lilith_batch(all_jds, T)
        mean_perigee = (mean_lilith + 180.0) % 360.0
        perigee_pert = _calc_elp2000_perigee_perturbations_batch(all_jds, T)
        interp_perigee_lon = (mean_perigee + perigee_pert) % 360.0
        # Lat/dist from osculating perigee (+e direction)
        ip_lat = np.degrees(np.arcsin(e_z / e_mag))
        ip_dist = p / (1.0 + e_mag)
        results[22] = np.column_stack([interp_perigee_lon, ip_lat, ip_dist])

    return results


def generate_ecliptic_bodies_vectorized(
    body_ids: List[int],
    jd_start: float,
    jd_end: float,
    verbose: bool = False,
) -> dict:
    """Generate Chebyshev coefficients for all ecliptic bodies in one vectorized pass.

    Uses a single Skyfield call for all bodies, then numpy for all post-processing.
    This is ~100x faster than the scalar per-body path for bodies that need Skyfield.

    For large date ranges (e.g. extended tier, 10,000 years) the evaluation is
    split into chunks to keep peak memory bounded (~1 GB per chunk instead of
    ~8 GB for 11M JDs at once).

    Args:
        body_ids: List of ecliptic body IDs to generate (subset of [10-13, 21, 22]).
        jd_start: Start Julian Day.
        jd_end: End Julian Day.
        verbose: Print progress.

    Returns:
        dict mapping body_id -> (body_id, coefficients_list, max_error)
    """
    # All ecliptic bodies use the same parameters (interval=8, degree=13)
    interval_days = 8.0
    degree = 13
    components = 3

    # Maximum JDs per evaluation chunk.  Each JD consumes ~500-800 bytes of
    # peak memory in the vectorized Skyfield call, so 500K JDs ≈ 300 MB peak.
    # With ~5 GB already used by planet/asteroid coefficients, this keeps total
    # RSS under ~5.4 GB, well within macOS limits.
    max_chunk_jds = 500_000

    if verbose:
        print(f"    Vectorized ecliptic generation: {len(body_ids)} bodies")

    # 1. Precompute all JDs
    all_jds, n_segments, pts_per_seg = _compute_all_segment_jds(
        jd_start, jd_end, interval_days, degree
    )

    if verbose:
        print(
            f"    Total JDs: {len(all_jds):,} ({n_segments} segments × {pts_per_seg} pts)"
        )

    # 2. Evaluate all bodies — chunked if needed to avoid OOM
    if len(all_jds) <= max_chunk_jds:
        body_values = _eval_ecliptic_bodies_batch(all_jds, body_ids, verbose=verbose)
    else:
        # Split into chunks aligned to segment boundaries (pts_per_seg points each)
        chunk_segs = max_chunk_jds // pts_per_seg
        chunk_jds = chunk_segs * pts_per_seg
        n_chunks = math.ceil(len(all_jds) / chunk_jds)

        accumulated: dict[int, list] = {bid: [] for bid in body_ids}
        bar = ProgressBar(
            total=len(all_jds),
            label="Ecliptic eval",
            unit="JDs",
            enabled=verbose,
        )
        for ci in range(n_chunks):
            start_idx = ci * chunk_jds
            end_idx = min(start_idx + chunk_jds, len(all_jds))
            jds_chunk = all_jds[start_idx:end_idx]
            chunk_results = _eval_ecliptic_bodies_batch(
                jds_chunk, body_ids, verbose=False
            )
            for bid in body_ids:
                if bid in chunk_results:
                    accumulated[bid].append(chunk_results[bid])
            # Free chunk intermediates
            del chunk_results, jds_chunk
            bar.update(end_idx)
        bar.finish()

        body_values = {
            bid: np.concatenate(parts) for bid, parts in accumulated.items() if parts
        }
        del accumulated

    # 3. Fit and verify each body
    results = {}
    for bid in body_ids:
        if bid not in body_values:
            continue
        label = BODY_NAMES.get(bid, f"Body {bid}")
        coeffs, error = _fit_and_verify_from_values_unwrap(
            body_values[bid],
            jd_start,
            jd_end,
            interval_days,
            degree,
            components,
            n_segments,
            pts_per_seg,
            label=label,
            verbose=verbose,
        )
        results[bid] = (bid, coeffs, error)

    return results


# =============================================================================
# BODY GENERATORS
# =============================================================================


def _init_skyfield():
    """Initialize Skyfield resources (called once per process)."""
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()
    return planets, ts


def _get_spk_jd_range(planets) -> Tuple[float, float]:
    """Get the valid JD range of the loaded SPK ephemeris.

    Returns:
        (jd_min, jd_max) covering all segments in the ephemeris.
    """
    T0 = 2451545.0  # J2000 in JD
    spk = planets.spk
    starts = [T0 + seg.start_second / 86400.0 for seg in spk.segments]
    ends = [T0 + seg.end_second / 86400.0 for seg in spk.segments]
    return min(starts), max(ends)


def _spk_covers_range(
    spk_file: str,
    body_id: int,
    jd_start: float,
    jd_end: float,
) -> bool:
    """Check if an SPK type 21 file covers the required JD range for a body.

    Opens the SPK, finds the target NAIF ID, and checks whether the
    segments span [jd_start, jd_end].

    Returns:
        True if coverage is sufficient, False otherwise.
    """
    try:
        from spktype21 import SPKType21

        kernel = SPKType21.open(spk_file)
    except Exception:
        return False

    try:
        naif_id = _ASTEROID_NAIF.get(body_id)
        if naif_id is None:
            return False

        T0 = 2451545.0
        # Find segments matching this target (try both NAIF conventions)
        target_segs = [
            seg
            for seg in kernel.segments
            if seg.target == naif_id
            or seg.target == naif_id + 18000000  # Horizons 20000000+N convention
        ]
        if not target_segs:
            # Try any non-center segment
            target_segs = [seg for seg in kernel.segments if seg.target not in (0, 10)]

        if not target_segs:
            return False

        seg_start = min(T0 + seg.start_second / 86400.0 for seg in target_segs)
        seg_end = max(T0 + seg.end_second / 86400.0 for seg in target_segs)

        # Allow 2-day tolerance: Horizons SPK boundaries can be off by ~1 day
        # from the requested range, and generation clamps to SPK range anyway.
        return seg_start <= jd_start + 2.0 and seg_end >= jd_end - 2.0
    finally:
        kernel.close()


def _find_covering_cached_spk(
    body_id: int,
    jd_start: float,
    jd_end: float,
) -> Optional[str]:
    """Return an already-cached SPK path that covers [jd_start, jd_end] for
    this body, or None.

    The auto-download registers the first kernel it finds (often a narrower
    range from a prior tier). When the registered kernel is too narrow, a
    wider one may still be sitting in the cache from an earlier run — prefer
    it over a slow, timeout-prone forced re-download of the whole tier range.
    """
    import glob

    from libephemeris.constants import SPK_BODY_NAME_MAP
    from libephemeris.spk_auto import DEFAULT_AUTO_SPK_DIR

    entry = SPK_BODY_NAME_MAP.get(body_id)
    if entry is None:
        return None
    # Numbered asteroids cache as "{number}_{start}_{end}.bsp"; the horizons_id
    # is the number string (drop the trailing ";" used for SBDB lookups).
    prefix = str(entry[0]).rstrip(";")
    cache_dir = os.environ.get("LIBEPHEMERIS_SPK_DIR", DEFAULT_AUTO_SPK_DIR)
    for path in sorted(glob.glob(os.path.join(cache_dir, f"{prefix}_*.bsp"))):
        if _spk_covers_range(path, body_id, jd_start, jd_end):
            return path
    return None


def _get_asteroid_spk_range(
    spk_file: str,
    body_id: int,
) -> Optional[Tuple[float, float]]:
    """Get the JD coverage range of an asteroid SPK type 21 file.

    Opens the SPK, finds the target NAIF ID, and returns the full coverage
    range as (jd_start, jd_end).

    Args:
        spk_file: Path to the SPK type 21 file.
        body_id: Internal body ID (SE_* constant) for NAIF lookup.

    Returns:
        (jd_start, jd_end) or None if the file cannot be read.
    """
    try:
        from spktype21 import SPKType21

        kernel = SPKType21.open(spk_file)
    except Exception:
        return None

    try:
        naif_id = _ASTEROID_NAIF.get(body_id)
        if naif_id is None:
            return None

        T0 = 2451545.0
        # Find segments matching this target (try both NAIF conventions)
        target_segs = [
            seg
            for seg in kernel.segments
            if seg.target == naif_id
            or seg.target == naif_id + 18000000  # Horizons 20000000+N convention
        ]
        if not target_segs:
            # Try any non-center segment
            target_segs = [seg for seg in kernel.segments if seg.target not in (0, 10)]

        if not target_segs:
            return None

        seg_start = min(T0 + seg.start_second / 86400.0 for seg in target_segs)
        seg_end = max(T0 + seg.end_second / 86400.0 for seg in target_segs)

        return (seg_start, seg_end)
    finally:
        kernel.close()


def _eval_target_vectorized(
    target,
    all_jds: np.ndarray,
    ts,
    spk_jd_min: float,
    spk_jd_max: float,
) -> np.ndarray:
    """Evaluate a Skyfield target at all JDs, with linear extrapolation for
    JDs outside the SPK ephemeris range.

    Last Chebyshev segments may extend their fit nodes a few days beyond
    jd_end (to keep full-width segment alignment with the reader).
    Instead of clamping (which corrupts the polynomial fit), out-of-range
    nodes are extrapolated linearly using position + velocity at the boundary.
    This preserves fit quality for in-range dates.

    Args:
        target: Skyfield VectorFunction (planet, barycenter, etc.)
        all_jds: Array of Julian Days (TT)
        ts: Skyfield Timescale
        spk_jd_min: Start of SPK valid range
        spk_jd_max: End of SPK valid range

    Returns:
        Array of shape (N, 3) with positions in AU.
    """
    margin = 1.0  # days of safety margin from SPK boundary
    lo = spk_jd_min + margin
    hi = spk_jd_max - margin

    in_range = (all_jds >= lo) & (all_jds <= hi)

    if np.all(in_range):
        # Fast path: everything in range
        t_arr = ts.tt_jd(all_jds)
        return np.asarray(target.at(t_arr).position.au).T  # (N, 3)

    # Split: evaluate in-range vectorized, extrapolate out-of-range
    all_values = np.empty((len(all_jds), 3))

    # In-range: vectorized evaluation
    in_idx = np.where(in_range)[0]
    if len(in_idx) > 0:
        t_in = ts.tt_jd(all_jds[in_idx])
        pos_in = np.asarray(target.at(t_in).position.au)  # (3, N_in)
        all_values[in_idx] = pos_in.T

    # Out-of-range: linear extrapolation from boundary
    out_idx = np.where(~in_range)[0]
    if len(out_idx) > 0:
        # Determine if overshoot is at start or end
        over_end = all_jds[out_idx] > hi
        over_start = ~over_end

        # Extrapolate from end boundary
        end_idx = out_idx[over_end]
        if len(end_idx) > 0:
            t_bnd = ts.tt_jd(hi)
            bnd_pos = target.at(t_bnd)
            p = np.asarray(bnd_pos.position.au).ravel()  # (3,)
            v = np.asarray(bnd_pos.velocity.au_per_d).ravel()  # (3,)
            dt = all_jds[end_idx] - hi  # days past boundary
            all_values[end_idx, 0] = p[0] + v[0] * dt
            all_values[end_idx, 1] = p[1] + v[1] * dt
            all_values[end_idx, 2] = p[2] + v[2] * dt

        # Extrapolate from start boundary
        start_idx = out_idx[over_start]
        if len(start_idx) > 0:
            t_bnd = ts.tt_jd(lo)
            bnd_pos = target.at(t_bnd)
            p = np.asarray(bnd_pos.position.au).ravel()
            v = np.asarray(bnd_pos.velocity.au_per_d).ravel()
            dt = all_jds[start_idx] - lo
            all_values[start_idx, 0] = p[0] + v[0] * dt
            all_values[start_idx, 1] = p[1] + v[1] * dt
            all_values[start_idx, 2] = p[2] + v[2] * dt

    return all_values


def _eval_body_icrs_vectorized(
    target_name: str,
    all_jds: np.ndarray,
    planets,
    ts,
) -> np.ndarray:
    """Get vectorized ICRS barycentric positions for a planet.

    Handles inner planets (direct Skyfield target) and outer planets
    (barycenter + SPK center offset or COB correction) transparently.

    For outer planets, uses a **hybrid per-JD approach** that matches the
    runtime behavior of _SpkCenterTarget: SPK center offsets are used for
    JDs within the planet_centers.bsp coverage, and analytical COB corrections
    are used for JDs outside that range. This ensures the stored Chebyshev
    data matches what calc() produces at runtime.

    JDs that extend beyond the SPK ephemeris range (from last-segment
    overshoot) are linearly extrapolated using position + velocity at
    the boundary. This preserves Chebyshev fit quality for in-range dates.

    Args:
        target_name: Planet name from _PLANET_MAP (e.g., 'jupiter', 'sun')
        all_jds: Array of Julian Days (TT)
        planets: Skyfield SpiceKernel ephemeris
        ts: Skyfield Timescale

    Returns:
        Array of shape (N, 3) with ICRS barycentric positions in AU.
    """
    from libephemeris.planets import _PLANET_FALLBACK, _PLANET_CENTER_NAIF_IDS
    from libephemeris.state import get_planet_center_segment

    spk_min, spk_max = _get_spk_jd_range(planets)

    # Try direct target first (works for inner planets: sun, moon, mercury, venus, earth)
    try:
        target = planets[target_name]
        return _eval_target_vectorized(target, all_jds, ts, spk_min, spk_max)
    except KeyError:
        pass

    # Outer planet: needs barycenter + center offset
    bary_name = _PLANET_FALLBACK.get(target_name)
    if bary_name is None:
        raise ValueError(f"No target or fallback for {target_name}")

    barycenter = planets[bary_name]

    # Get barycenter positions for ALL JDs (vectorized, always available)
    bary_vals = _eval_target_vectorized(barycenter, all_jds, ts, spk_min, spk_max)

    # Hybrid SPK/COB approach: use SPK center where in range, analytical COB
    # where out of range. This matches _SpkCenterTarget's runtime behavior.
    if target_name in _PLANET_CENTER_NAIF_IDS:
        naif_id = _PLANET_CENTER_NAIF_IDS[target_name]
        center_segment = get_planet_center_segment(naif_id)
        if center_segment is not None:
            # Get the SPK center segment's JD range
            spk_seg = center_segment.spk_segment
            center_jd_min = spk_seg.start_jd
            center_jd_max = spk_seg.end_jd
            margin = 1.0  # days safety margin
            clo = center_jd_min + margin
            chi = center_jd_max - margin

            in_spk = (all_jds >= clo) & (all_jds <= chi)
            in_spk_idx = np.where(in_spk)[0]
            out_spk_idx = np.where(~in_spk)[0]

            # Apply SPK center offset for in-range JDs (vectorized)
            if len(in_spk_idx) > 0:
                t_in = ts.tt_jd(all_jds[in_spk_idx])
                offset_pos = np.asarray(center_segment.at(t_in).position.au)  # (3, N)
                bary_vals[in_spk_idx] += offset_pos.T

            # Apply analytical COB for out-of-range JDs (scalar)
            if len(out_spk_idx) > 0:
                from libephemeris.moon_theories import get_cob_offset

                for i in out_spk_idx:
                    t_single = ts.tt_jd(float(all_jds[i]))
                    offset = get_cob_offset(bary_name, t_single)
                    bary_vals[i, 0] += offset[0]
                    bary_vals[i, 1] += offset[1]
                    bary_vals[i, 2] += offset[2]

            return bary_vals

    # No SPK center available at all: pure analytical COB fallback
    from libephemeris.moon_theories import get_cob_offset

    for i in range(len(all_jds)):
        t_single = ts.tt_jd(float(all_jds[i]))
        offset = get_cob_offset(bary_name, t_single)
        bary_vals[i, 0] += offset[0]
        bary_vals[i, 1] += offset[1]
        bary_vals[i, 2] += offset[2]

    return bary_vals


def generate_body_icrs(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for an ICRS barycentric body.

    Uses vectorized Skyfield evaluation: all JDs (fit nodes + verification)
    are computed in a single batch call for ~150x speedup over scalar loops.

    Returns:
        (list_of_coefficient_arrays, max_error_au)
    """
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()

    target_name = _PLANET_MAP.get(body_id)
    if target_name is None:
        raise ValueError(f"No planet map entry for body_id={body_id}")

    # Precompute all JDs for all segments (fit + verify)
    all_jds, n_segments, pts_per_seg = _compute_all_segment_jds(
        jd_start, jd_end, interval_days, degree
    )

    # Single vectorized evaluation (handles COB correction for outer planets)
    all_values = _eval_body_icrs_vectorized(target_name, all_jds, planets, ts)

    return _fit_and_verify_from_values(
        all_values,
        jd_start,
        jd_end,
        interval_days,
        degree,
        3,
        n_segments,
        pts_per_seg,
        label=label,
        verbose=verbose,
    )


def generate_body_icrs_system_bary(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for a system barycenter (ICRS).

    Stores the pure system barycenter position (no COB correction),
    which is ultra-smooth and fits Chebyshev polynomials with negligible error.
    COB correction is applied at runtime by fast_calc._apply_cob_correction().

    This eliminates the high-frequency moon oscillations (e.g., Io's 1.77-day
    period for Jupiter, Charon's 6.4-day period for Pluto) that cause
    Chebyshev fitting errors of 0.01-0.15" when storing planet center positions.

    Returns:
        (list_of_coefficient_arrays, max_error_au)
    """
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()

    bary_name = _SYSTEM_BARY_MAP.get(body_id)
    if bary_name is None:
        raise ValueError(f"No system barycenter map entry for body_id={body_id}")

    spk_min, spk_max = _get_spk_jd_range(planets)
    barycenter = planets[bary_name]

    # Precompute all JDs for all segments (fit + verify)
    all_jds, n_segments, pts_per_seg = _compute_all_segment_jds(
        jd_start, jd_end, interval_days, degree
    )

    # Pure system barycenter evaluation (no COB — ultra-smooth)
    all_values = _eval_target_vectorized(barycenter, all_jds, ts, spk_min, spk_max)

    return _fit_and_verify_from_values(
        all_values,
        jd_start,
        jd_end,
        interval_days,
        degree,
        3,
        n_segments,
        pts_per_seg,
        label=label,
        verbose=verbose,
    )


def generate_body_icrs_asteroid(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for an asteroid (ICRS barycentric).

    Uses spktype21 to read the SPK type 21 file directly (~36x faster than
    scalar libephemeris calls). The SPK provides heliocentric positions
    (center=10/Sun), so the Sun's barycentric position is added via a single
    vectorized Skyfield call.

    Raises RuntimeError if the SPK file is not available or does not cover
    the requested date range. Never falls back to Keplerian — that would
    produce errors of degrees over decades.

    Returns:
        (list_of_coefficient_arrays, max_error_au)
    """
    from libephemeris.state import get_planets, get_timescale, _SPK_BODY_MAP

    planets = get_planets()
    ts = get_timescale()

    # Check if we have an SPK file for this asteroid
    spk_info = _SPK_BODY_MAP.get(body_id)
    ast_name = BODY_NAMES.get(body_id, f"Body {body_id}")

    if spk_info is None:
        raise RuntimeError(
            f"No SPK kernel registered for {ast_name} (body {body_id}). "
            f"Cannot generate LEB data without SPK — Keplerian fallback "
            f"produces errors of degrees over decades. Use "
            f"auto_download_asteroid_spk() first or exclude this body."
        )

    spk_file, naif_id = spk_info
    try:
        from spktype21 import SPKType21

        kernel = SPKType21.open(spk_file)
    except Exception as exc:
        raise RuntimeError(
            f"Cannot open SPK file for {ast_name}: {spk_file}: {exc}"
        ) from exc

    try:
        # Find the center ID and target from the kernel segments
        center_id = kernel.segments[0].center  # typically 10 (Sun)
        target_id = kernel.segments[0].target

        # Check if the asteroid SPK covers the full requested range
        # (with some margin for last-segment overshoot)
        T0_SPK = 2451545.0
        ast_jd_min = min(
            T0_SPK + seg.start_second / 86400.0
            for seg in kernel.segments
            if seg.target == target_id
        )
        ast_jd_max = max(
            T0_SPK + seg.end_second / 86400.0
            for seg in kernel.segments
            if seg.target == target_id
        )

        if ast_jd_min > jd_start + 2.0 or ast_jd_max < jd_end - 2.0:
            kernel.close()
            raise RuntimeError(
                f"SPK for {ast_name} covers JD {ast_jd_min:.1f}–{ast_jd_max:.1f} "
                f"but requested range is JD {jd_start:.1f}–{jd_end:.1f}. "
                f"Cannot generate LEB data — SPK coverage is insufficient. "
                f"Use a narrower date range or exclude this body."
            )
    except RuntimeError:
        raise
    except Exception as exc:
        kernel.close()
        raise RuntimeError(f"Error reading SPK segments for {ast_name}: {exc}") from exc

    try:
        center_id = kernel.segments[0].center
        target_id = kernel.segments[0].target

        # Recompute asteroid SPK range for clamping
        T0_SPK2 = 2451545.0
        ast_jd_lo = min(
            T0_SPK2 + seg.start_second / 86400.0
            for seg in kernel.segments
            if seg.target == target_id
        )
        ast_jd_hi = max(
            T0_SPK2 + seg.end_second / 86400.0
            for seg in kernel.segments
            if seg.target == target_id
        )

        # Precompute all JDs
        all_jds, n_segments, pts_per_seg = _compute_all_segment_jds(
            jd_start, jd_end, interval_days, degree
        )

        # Get Sun barycentric positions with extrapolation for overshoot
        spk_min, spk_max = _get_spk_jd_range(planets)
        sun = planets["sun"]
        sun_bary = _eval_target_vectorized(sun, all_jds, ts, spk_min, spk_max)  # (N, 3)

        # Clamp to asteroid SPK range for spktype21 calls
        ast_clamped = np.clip(all_jds, ast_jd_lo + 0.01, ast_jd_hi - 0.01)

        # Compute asteroid heliocentric positions via spktype21 (scalar loop)
        AU_KM = 149597870.7
        all_values = np.empty((len(all_jds), 3))
        spk_bar = ProgressBar(
            len(all_jds),
            label=label + " (spk21)" if label else "spk21 eval",
            unit="pt",
            enabled=verbose,
        )
        for i in range(len(all_jds)):
            pos_km, _ = kernel.compute_type21(
                center_id, target_id, float(ast_clamped[i])
            )
            # helio km -> helio AU -> + Sun bary = SSB bary (ICRS)
            all_values[i, 0] = pos_km[0] / AU_KM + sun_bary[i, 0]
            all_values[i, 1] = pos_km[1] / AU_KM + sun_bary[i, 1]
            all_values[i, 2] = pos_km[2] / AU_KM + sun_bary[i, 2]
            spk_bar.update(i + 1)
        spk_bar.finish()

        return _fit_and_verify_from_values(
            all_values,
            jd_start,
            jd_end,
            interval_days,
            degree,
            3,
            n_segments,
            pts_per_seg,
            label=label,
            verbose=verbose,
        )
    except Exception as exc:
        raise RuntimeError(
            f"SPK evaluation failed for {ast_name} (body {body_id}): {exc}"
        ) from exc
    finally:
        kernel.close()


def generate_body_ecliptic(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for an ecliptic-direct body.

    Handles longitude unwrapping before fitting.
    """
    from libephemeris.lunar import (
        calc_mean_lunar_node,
        calc_true_lunar_node,
        calc_mean_lilith_with_latitude,
        calc_true_lilith,
        calc_interpolated_apogee,
        calc_interpolated_perigee,
    )

    # Map body_id to evaluation function
    eval_funcs = {
        10: lambda jd: np.array([calc_mean_lunar_node(jd), 0.0, 0.0]),
        11: lambda jd: np.array(calc_true_lunar_node(jd)),
        12: lambda jd: np.array([*calc_mean_lilith_with_latitude(jd), 0.0]),
        13: lambda jd: np.array(calc_true_lilith(jd)),
        21: lambda jd: np.array(calc_interpolated_apogee(jd)),
        22: lambda jd: np.array(calc_interpolated_perigee(jd)),
    }

    if body_id not in eval_funcs:
        raise ValueError(f"No ecliptic eval function for body_id={body_id}")

    raw_func = eval_funcs[body_id]

    # Generate with longitude unwrapping
    return _generate_segments_unwrap(
        raw_func,
        jd_start,
        jd_end,
        interval_days,
        degree,
        3,
        label=label,
        verbose=verbose,
    )


def generate_body_helio(
    body_id: int,
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for a heliocentric ecliptic body.

    Handles Uranian hypotheticals and Transpluto.
    """
    from libephemeris.hypothetical import calc_uranian_planet, calc_transpluto

    if body_id == 48:  # ISIS / Transpluto

        def eval_func(jd: float) -> np.ndarray:
            result = calc_transpluto(jd)
            return np.array([result[0], result[1], result[2]])
    elif 40 <= body_id <= 47:  # Uranian planets

        def eval_func(jd: float) -> np.ndarray:
            result = calc_uranian_planet(body_id, jd)
            return np.array([result[0], result[1], result[2]])
    else:
        raise ValueError(f"No heliocentric eval function for body_id={body_id}")

    # Generate with longitude unwrapping
    return _generate_segments_unwrap(
        eval_func,
        jd_start,
        jd_end,
        interval_days,
        degree,
        3,
        label=label,
        verbose=verbose,
    )


def _generate_segments(
    func: Callable[[float], np.ndarray],
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    components: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev segments for a function (no unwrapping).

    Returns:
        (list_of_coefficient_arrays, max_error)
    """
    n_segments = int(math.ceil((jd_end - jd_start) / interval_days))
    all_coeffs = []
    max_error = 0.0
    bar = ProgressBar(n_segments, label=label or "segments", enabled=verbose)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        # IMPORTANT: Always use full-width segments. The reader maps tau
        # using interval_days, so truncating the last segment would cause
        # a mismatch between the fitted polynomial domain and the reader's
        # tau computation.  The underlying functions are valid beyond jd_end.
        seg_end = seg_start + interval_days

        coeffs = fit_segment(func, seg_start, seg_end, degree, components)
        # Verify only within the actual requested range, but tau is always
        # relative to the full segment [seg_start, seg_end].
        v_end = min(seg_end, jd_end)
        error = verify_segment(
            func, coeffs, seg_start, seg_end, components, verify_end=v_end
        )
        if error > max_error:
            max_error = error

        all_coeffs.append(coeffs)
        bar.update(i + 1)

    bar.finish()
    return all_coeffs, max_error


def _generate_segments_unwrap(
    func: Callable[[float], np.ndarray],
    jd_start: float,
    jd_end: float,
    interval_days: float,
    degree: int,
    components: int,
    label: str = "",
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev segments with longitude unwrapping (component 0).

    The first component (longitude) is unwrapped before fitting to handle
    the 0/360 discontinuity. After fitting, verification re-wraps with % 360.
    """
    n_segments = int(math.ceil((jd_end - jd_start) / interval_days))
    all_coeffs = []
    max_error = 0.0
    bar = ProgressBar(n_segments, label=label or "segments", enabled=verbose)

    for i in range(n_segments):
        seg_start = jd_start + i * interval_days
        # Always use full-width segments (see _generate_segments comment)
        seg_end = seg_start + interval_days

        nodes = chebyshev_nodes(degree + 1)
        jd_nodes = 0.5 * (seg_end - seg_start) * nodes + 0.5 * (seg_start + seg_end)

        # Evaluate function at nodes
        values = np.array([func(jd) for jd in jd_nodes])

        # Unwrap longitude (component 0) to remove 360-degree jumps
        values[:, 0] = np.unwrap(np.radians(values[:, 0]))
        values[:, 0] = np.degrees(values[:, 0])

        # Fit each component
        coeffs = np.zeros((components, degree + 1))
        for c in range(components):
            coeffs[c] = chebfit(nodes, values[:, c], degree)

        # Verify with re-wrapping (only within actual requested range)
        v_end = min(seg_end, jd_end)
        error = _verify_segment_unwrapped(
            func, coeffs, seg_start, seg_end, components, verify_end=v_end
        )
        if error > max_error:
            max_error = error

        all_coeffs.append(coeffs)
        bar.update(i + 1)

    bar.finish()
    return all_coeffs, max_error


def _verify_segment_unwrapped(
    func: Callable[[float], np.ndarray],
    coeffs: np.ndarray,
    seg_start: float,
    seg_end: float,
    components: int,
    n_test: int = 10,
    verify_end: float | None = None,
) -> float:
    """Verify a Chebyshev fit for ecliptic bodies (with longitude re-wrapping).

    Args:
        seg_start: Segment start JD (defines the polynomial domain).
        seg_end: Segment end JD (defines the polynomial domain).
        verify_end: If given, only sample verification points up to this JD.
            Tau is always computed relative to [seg_start, seg_end].
    """
    if verify_end is None:
        verify_end = seg_end
    mid = 0.5 * (seg_start + seg_end)
    half = 0.5 * (seg_end - seg_start)
    max_error = 0.0
    for i in range(n_test):
        frac = (i + 0.5) / n_test
        jd = seg_start + frac * (verify_end - seg_start)
        tau = (jd - mid) / half

        reference = func(jd)
        for c in range(components):
            fitted = chebval(tau, coeffs[c])
            if c == 0:
                # Re-wrap longitude
                fitted = fitted % 360.0
                ref = reference[c] % 360.0
                # Handle wrap-around comparison
                error = abs(fitted - ref)
                if error > 180.0:
                    error = 360.0 - error
            else:
                error = abs(fitted - reference[c])
            if error > max_error:
                max_error = error

    return max_error


# =============================================================================
# NUTATION GENERATOR
# =============================================================================


def generate_nutation(
    jd_start: float,
    jd_end: float,
    verbose: bool = False,
) -> Tuple[List[np.ndarray], float]:
    """Generate Chebyshev coefficients for nutation (dpsi, deps in radians).

    Uses vectorized erfa.nut06a (IAU 2006/2000A) for ~50x speedup over
    scalar get_cached_nutation() calls.
    """
    import erfa

    # Precompute all JDs for all segments
    all_jds, n_segments, pts_per_seg = _compute_all_segment_jds(
        jd_start, jd_end, NUTATION_INTERVAL, NUTATION_DEGREE
    )
    n_points = len(all_jds)

    # Convert JD to TT (J2000 epoch split for erfa)
    # erfa uses (jd1, jd2) split: jd1 = 2451545.0 (J2000), jd2 = jd - 2451545.0
    jd1 = np.full_like(all_jds, 2451545.0)
    jd2 = all_jds - 2451545.0

    # Chunked erfa computation with progress bar for large arrays
    CHUNK_SIZE = 500_000
    if n_points > CHUNK_SIZE:
        dpsi_list = []
        deps_list = []
        bar = ProgressBar(
            total=n_points,
            label="Nutation (erfa)",
            unit="pts",
            enabled=verbose,
        )
        for i in range(0, n_points, CHUNK_SIZE):
            end_i = min(i + CHUNK_SIZE, n_points)
            dpsi_chunk, deps_chunk = erfa.nut06a(jd1[i:end_i], jd2[i:end_i])
            dpsi_list.append(dpsi_chunk)
            deps_list.append(deps_chunk)
            bar.update(end_i)
        bar.finish()
        dpsi = np.concatenate(dpsi_list)
        deps = np.concatenate(deps_list)
    else:
        if verbose:
            print(f"Computing nutation for {n_points:,} points...")
        dpsi, deps = erfa.nut06a(jd1, jd2)

    # Stack into (N, 2) array
    all_values = np.column_stack([dpsi, deps])

    return _fit_and_verify_from_values(
        all_values,
        jd_start,
        jd_end,
        NUTATION_INTERVAL,
        NUTATION_DEGREE,
        NUTATION_COMPONENTS,
        n_segments,
        pts_per_seg,
        label="Nutation",
        verbose=verbose,
    )


# =============================================================================
# DELTA-T GENERATOR
# =============================================================================


def generate_delta_t(
    jd_start: float,
    jd_end: float,
) -> List[Tuple[float, float]]:
    """Generate Delta-T sparse table.

    Samples deltat() every DELTA_T_INTERVAL days.

    Returns:
        List of (jd, delta_t_days) tuples.
    """
    from libephemeris.time_utils import deltat

    table = []
    jd = jd_start
    while jd <= jd_end:
        dt = deltat(jd)
        table.append((jd, dt))
        jd += DELTA_T_INTERVAL

    # Ensure we include the end point
    if table[-1][0] < jd_end:
        table.append((jd_end, deltat(jd_end)))

    return table


# =============================================================================
# STAR CATALOG GENERATOR
# =============================================================================


def generate_star_catalog() -> List[StarEntry]:
    """Extract star catalog entries from libephemeris.

    Returns:
        List of StarEntry records.
    """
    from libephemeris.fixed_stars import STAR_CATALOG

    entries = []
    for star in STAR_CATALOG:
        entries.append(
            StarEntry(
                star_id=star.id,
                ra_j2000=star.data.ra_j2000,
                dec_j2000=star.data.dec_j2000,
                pm_ra=star.data.pm_ra / 3600.0,  # arcsec/yr -> deg/yr
                pm_dec=star.data.pm_dec / 3600.0,  # arcsec/yr -> deg/yr
                parallax=0.0,  # Not in StarData
                rv=0.0,  # Not in StarData
                magnitude=star.magnitude,
            )
        )

    return entries


# =============================================================================
# FILE ASSEMBLY
# =============================================================================


def _year_to_jd(year: int) -> float:
    """Convert a year to Julian Day (January 1.0)."""
    from libephemeris.time_utils import julday

    return julday(year, 1, 1, 0.0)


def generate_single_body(
    body_id: int,
    jd_start: float,
    jd_end: float,
    verbose: bool = False,
) -> Tuple[int, List[np.ndarray], float]:
    """Generate Chebyshev data for a single body.

    Returns:
        (body_id, coefficients_list, max_error)
    """
    params = BODY_PARAMS[body_id]
    interval_days, degree, coord_type, components = params
    label = BODY_NAMES.get(body_id, f"Body {body_id}")

    if coord_type == COORD_ICRS_BARY:
        if body_id in _PLANET_MAP:
            coeffs, error = generate_body_icrs(
                body_id,
                jd_start,
                jd_end,
                interval_days,
                degree,
                label=label,
                verbose=verbose,
            )
        else:
            coeffs, error = generate_body_icrs_asteroid(
                body_id,
                jd_start,
                jd_end,
                interval_days,
                degree,
                label=label,
                verbose=verbose,
            )
    elif coord_type == COORD_ICRS_BARY_SYSTEM:
        # System barycenter — store pure barycenter, COB applied at runtime
        coeffs, error = generate_body_icrs_system_bary(
            body_id,
            jd_start,
            jd_end,
            interval_days,
            degree,
            label=label,
            verbose=verbose,
        )
    elif coord_type == COORD_ECLIPTIC:
        coeffs, error = generate_body_ecliptic(
            body_id,
            jd_start,
            jd_end,
            interval_days,
            degree,
            label=label,
            verbose=verbose,
        )
    elif coord_type == COORD_HELIO_ECL:
        coeffs, error = generate_body_helio(
            body_id,
            jd_start,
            jd_end,
            interval_days,
            degree,
            label=label,
            verbose=verbose,
        )
    else:
        raise ValueError(f"Unknown coord_type {coord_type} for body {body_id}")

    return body_id, coeffs, error


def assemble_leb(
    output: str,
    jd_start: float,
    jd_end: float,
    bodies: Optional[List[int]] = None,
    workers: int = 1,
    verbose: bool = True,
    skip_aux: bool = False,
) -> None:
    """Assemble a complete .leb file.

    Args:
        output: Output file path.
        jd_start: Start Julian Day.
        jd_end: End Julian Day.
        bodies: List of body IDs to include (None = all from BODY_PARAMS).
        workers: Number of parallel workers for body generation.
        verbose: Print progress messages.
        skip_aux: If True, skip nutation, Delta-T and star catalog generation.
            The resulting file will only contain body coefficients.  Useful
            when generating single-body partial files that will be merged
            later (merge takes aux data from the first file that has it).
    """
    if bodies is None:
        bodies = sorted(BODY_PARAMS.keys())

    now_jd = J2000 + (time.time() / 86400.0 - 10957.5)  # Approximate current JD

    if verbose:
        print(f"Generating LEB file: {output}")
        print(f"  Date range: JD {jd_start:.1f} to {jd_end:.1f}")
        print(f"  Bodies: {len(bodies)}")
        print(f"  Workers: {workers}")
        print()

    # -------------------------------------------------------------------------
    # 0. Ensure SPK kernels for asteroids; discover per-body date ranges
    # -------------------------------------------------------------------------
    # Asteroids may have SPK coverage narrower than the tier range.
    # Instead of excluding them, we include them with their actual SPK range.
    # The LEB format already supports per-body jd_start/jd_end, and the
    # reader raises ValueError for out-of-range JDs, which calc_ut()
    # catches to fall through to Skyfield.
    asteroid_bodies = [b for b in bodies if b in _ASTEROID_NAIF]
    excluded_asteroids: List[int] = []
    # Per-body date ranges: body_id -> (jd_start, jd_end)
    # Non-asteroids use the global range; asteroids use their SPK range.
    body_jd_ranges: dict[int, Tuple[float, float]] = {}

    # Minimum useful SPK coverage (years). Asteroids with less than this
    # are excluded — too narrow to be worth including.
    MIN_ASTEROID_COVERAGE_DAYS = 365.25 * 20  # 20 years

    if asteroid_bodies:
        import libephemeris as ephem
        from libephemeris.minor_bodies import (
            auto_download_asteroid_spk,
            is_spk_available_for_body,
        )
        from libephemeris.state import _SPK_BODY_MAP

        # Enable auto-download for SPK acquisition
        os.environ["LIBEPHEMERIS_AUTO_SPK"] = "1"
        ephem.set_auto_spk_download(True)

        if verbose:
            print("  Preparing asteroid SPK kernels...")
            print(f"    Tier range: JD {jd_start:.1f} to {jd_end:.1f}")

        for bid in asteroid_bodies:
            name = BODY_NAMES.get(bid, f"Body {bid}")

            # Check if an already-cached SPK covers the full range.
            # auto_download_asteroid_spk() short-circuits when a body is
            # already registered, even if its coverage is too narrow (e.g.
            # a previously downloaded ±10yr SPK).  We detect that here and
            # force a re-download with the full tier range.
            need_force = False
            if bid in _SPK_BODY_MAP:
                cached_file, _ = _SPK_BODY_MAP[bid]
                if not _spk_covers_range(cached_file, bid, jd_start, jd_end):
                    # The registered kernel is too narrow — prefer a wider one
                    # already in the cache over a forced (slow, timeout-prone)
                    # re-download of the whole tier range.
                    covering = _find_covering_cached_spk(bid, jd_start, jd_end)
                    if covering is not None:
                        from libephemeris import spk as _spk

                        _spk.register_spk_body(bid, covering, _ASTEROID_NAIF[bid])
                        if verbose:
                            print(
                                f"    {name}: using cached wider SPK "
                                f"{os.path.basename(covering)}"
                            )
                    else:
                        need_force = True
                        if verbose:
                            print(f"    {name}: cached SPK too narrow, re-downloading...")

            try:
                auto_download_asteroid_spk(
                    bid,
                    jd_start=jd_start,
                    jd_end=jd_end,
                    force=need_force,
                )
            except Exception as exc:
                if verbose:
                    print(f"    {name}: SPK download failed: {exc}")

            if not is_spk_available_for_body(bid) or bid not in _SPK_BODY_MAP:
                excluded_asteroids.append(bid)
                if verbose:
                    print(f"    {name}: no SPK available (EXCLUDED)")
                continue

            spk_file, _ = _SPK_BODY_MAP[bid]

            # Full coverage?
            if _spk_covers_range(spk_file, bid, jd_start, jd_end):
                body_jd_ranges[bid] = (jd_start, jd_end)
                if verbose:
                    print(f"    {name}: SPK covers full tier range (spk21)")
                continue

            # Partial coverage — discover the actual SPK range
            spk_range = _get_asteroid_spk_range(spk_file, bid)
            if spk_range is None:
                excluded_asteroids.append(bid)
                if verbose:
                    print(f"    {name}: cannot read SPK range (EXCLUDED)")
                continue

            spk_jd_start, spk_jd_end = spk_range
            # Intersect SPK range with tier range
            eff_start = max(spk_jd_start, jd_start)
            eff_end = min(spk_jd_end, jd_end)
            eff_days = eff_end - eff_start

            if eff_days < MIN_ASTEROID_COVERAGE_DAYS:
                excluded_asteroids.append(bid)
                if verbose:
                    eff_years = eff_days / 365.25
                    print(
                        f"    {name}: SPK range JD {spk_jd_start:.1f}–{spk_jd_end:.1f} "
                        f"overlaps only {eff_years:.0f} years with tier (EXCLUDED)"
                    )
                continue

            # Use per-body range (with 1-day inward margin for SPK boundary safety)
            body_jd_ranges[bid] = (eff_start + 1.0, eff_end - 1.0)
            if verbose:
                eff_start_yr = 2000.0 + (eff_start - 2451545.0) / 365.25
                eff_end_yr = 2000.0 + (eff_end - 2451545.0) / 365.25
                print(
                    f"    {name}: per-body range JD {eff_start:.1f}–{eff_end:.1f} "
                    f"(~{eff_start_yr:.0f}–{eff_end_yr:.0f} CE)"
                )

        # Remove unavailable asteroids from the body list
        if excluded_asteroids:
            bodies = [b for b in bodies if b not in excluded_asteroids]
            if verbose:
                excluded_names = [BODY_NAMES.get(b, str(b)) for b in excluded_asteroids]
                print(
                    f"\n  WARNING: {len(excluded_asteroids)} asteroid(s) excluded "
                    f"(no SPK available): {', '.join(excluded_names)}"
                )
                print(
                    "  Excluded asteroids will NOT be in the LEB file. "
                    "They will use live Skyfield/SPK at runtime."
                )
            # Fail fast on an UNEXPECTED exotic drop: every registry body is
            # supposed to have an obtainable SPK over its coverage window, so a
            # silent exclusion would ship an incomplete file that looks fine.
            dropped_exotics = sorted(set(excluded_asteroids) & set(EXOTIC_IDS))
            if dropped_exotics:
                names = [BODY_NAMES.get(b, str(b)) for b in dropped_exotics]
                raise RuntimeError(
                    "Exotic minor bodies were dropped (no SPK obtainable): "
                    + ", ".join(names)
                    + ". Pre-download their SPK (scripts/download_max_range_spk.py) "
                    "or remove them from libephemeris.exotic_bodies."
                )

        if verbose:
            print()

    # -------------------------------------------------------------------------
    # 1. Generate body coefficients
    # -------------------------------------------------------------------------
    body_data: dict[int, List[np.ndarray]] = {}
    body_errors: dict[int, float] = {}

    # Categorize bodies by generation strategy:
    # - Planets (in _PLANET_MAP): vectorized Skyfield — ICRS_BARY
    # - Asteroids (in _ASTEROID_NAIF): SPK-based — ICRS_BARY
    # - Ecliptic/Helio: scalar analytical funcs (slow, parallelize across workers)
    planet_bodies = [b for b in bodies if b in _PLANET_MAP]
    asteroid_bodies_gen = [b for b in bodies if b in _ASTEROID_NAIF]
    analytical_bodies = [
        b for b in bodies if b not in _PLANET_MAP and b not in _ASTEROID_NAIF
    ]

    t0 = time.time()

    # 1a. Generate planets (vectorized Skyfield)
    if planet_bodies and verbose:
        print("  --- Planets (vectorized Skyfield) ---")
    for bid in planet_bodies:
        bid, coeffs, error = generate_single_body(
            bid, jd_start, jd_end, verbose=verbose
        )
        body_data[bid] = coeffs
        body_errors[bid] = error

    # 1b. Generate asteroids (SPK-based)
    # Each asteroid uses its own date range (from body_jd_ranges) which may
    # be narrower than the tier range when SPK coverage is partial.
    if asteroid_bodies_gen and verbose:
        print("  --- Asteroids ---")
    for bid in asteroid_bodies_gen:
        ast_start, ast_end = body_jd_ranges.get(bid, (jd_start, jd_end))
        bid, coeffs, error = generate_single_body(
            bid, ast_start, ast_end, verbose=verbose
        )
        body_data[bid] = coeffs
        body_errors[bid] = error

    # 1c. Generate analytical bodies
    # Split into ecliptic (vectorized) and heliocentric (scalar, already fast).
    # The vectorized ecliptic path assumes uniform params (8d/13), so bodies
    # with non-standard params are routed to the scalar path instead.
    _STD_ECL_PARAMS = (8.0, 13)  # (interval_days, degree) for vectorized path
    ecliptic_body_ids_vec = [
        b
        for b in analytical_bodies
        if BODY_PARAMS[b][2] == COORD_ECLIPTIC
        and (BODY_PARAMS[b][0], BODY_PARAMS[b][1]) == _STD_ECL_PARAMS
    ]
    ecliptic_body_ids_scalar = [
        b
        for b in analytical_bodies
        if BODY_PARAMS[b][2] == COORD_ECLIPTIC
        and (BODY_PARAMS[b][0], BODY_PARAMS[b][1]) != _STD_ECL_PARAMS
    ]
    helio_body_ids = [
        b for b in analytical_bodies if BODY_PARAMS[b][2] == COORD_HELIO_ECL
    ]

    # Ecliptic bodies with standard params: single vectorized Skyfield call + numpy
    # This replaces ~328K scalar Skyfield calls per body with ONE array call.
    if ecliptic_body_ids_vec:
        if verbose:
            print(
                f"  --- Ecliptic bodies ({len(ecliptic_body_ids_vec)} bodies, vectorized) ---"
            )
        vec_results = generate_ecliptic_bodies_vectorized(
            ecliptic_body_ids_vec, jd_start, jd_end, verbose=verbose
        )
        for bid, (_, coeffs, error) in vec_results.items():
            body_data[bid] = coeffs
            body_errors[bid] = error

    # Ecliptic bodies with non-standard params: scalar path (per-body params)
    if ecliptic_body_ids_scalar:
        if verbose:
            print(
                f"  --- Ecliptic bodies ({len(ecliptic_body_ids_scalar)} bodies, scalar) ---"
            )
        for bid in ecliptic_body_ids_scalar:
            bid, coeffs, error = generate_single_body(
                bid, jd_start, jd_end, verbose=verbose
            )
            body_data[bid] = coeffs
            body_errors[bid] = error

    # Heliocentric bodies (Uranians 40-48, Transpluto): scalar but pure math, fast
    if helio_body_ids:
        if verbose:
            print(
                f"  --- Heliocentric bodies ({len(helio_body_ids)} bodies, scalar) ---"
            )
        for bid in helio_body_ids:
            bid, coeffs, error = generate_single_body(
                bid, jd_start, jd_end, verbose=verbose
            )
            body_data[bid] = coeffs
            body_errors[bid] = error

    t_bodies = time.time() - t0
    if verbose:
        print(f"\n  Body generation: {t_bodies:.1f}s")

    # -------------------------------------------------------------------------
    # 2. Generate nutation
    # -------------------------------------------------------------------------
    nutation_coeffs: List[np.ndarray] = []
    nutation_error = 0.0
    if skip_aux:
        if verbose:
            print("  Skipping nutation (--skip-aux)")
    else:
        if verbose:
            print("  Generating nutation...")
        t0 = time.time()
        nutation_coeffs, nutation_error = generate_nutation(
            jd_start, jd_end, verbose=verbose
        )
        t_nut = time.time() - t0
        if verbose:
            print(
                f"{len(nutation_coeffs)} segments, max error={nutation_error:.2e} rad ({t_nut:.1f}s)"
            )

    # -------------------------------------------------------------------------
    # 3. Generate Delta-T
    # -------------------------------------------------------------------------
    delta_t_table: List[Tuple[float, float]] = []
    if skip_aux:
        if verbose:
            print("  Skipping Delta-T (--skip-aux)")
    else:
        if verbose:
            print("  Generating Delta-T...", end=" ", flush=True)
        t0 = time.time()
        delta_t_table = generate_delta_t(jd_start, jd_end)
        t_dt = time.time() - t0
        if verbose:
            print(f"{len(delta_t_table)} samples ({t_dt:.1f}s)")

    # -------------------------------------------------------------------------
    # 4. Generate star catalog
    # -------------------------------------------------------------------------
    star_entries: List[StarEntry] = []
    if skip_aux:
        if verbose:
            print("  Skipping star catalog (--skip-aux)")
    else:
        if verbose:
            print("  Generating star catalog...", end=" ", flush=True)
        star_entries = generate_star_catalog()
        if verbose:
            print(f"{len(star_entries)} stars")

    # -------------------------------------------------------------------------
    # 5. Calculate layout and write file
    # -------------------------------------------------------------------------
    if verbose:
        print("\n  Assembling file...", end=" ", flush=True)

    # Calculate section sizes
    body_count = len(bodies)
    body_index_size = body_count * BODY_ENTRY_SIZE

    # Chebyshev data size (sum of all body segments)
    chebyshev_size = 0
    for bid in bodies:
        params = BODY_PARAMS[bid]
        interval_days, degree, coord_type, components = params
        seg_size = segment_byte_size(degree, components)
        chebyshev_size += len(body_data[bid]) * seg_size

    # Nutation data size
    nut_data_size = NUTATION_HEADER_SIZE + len(nutation_coeffs) * segment_byte_size(
        NUTATION_DEGREE, NUTATION_COMPONENTS
    )

    # Delta-T size
    delta_t_size = DELTA_T_HEADER_SIZE + len(delta_t_table) * DELTA_T_ENTRY_SIZE

    # Star catalog size
    star_size = len(star_entries) * STAR_ENTRY_SIZE

    # Section directory
    section_dir_total = NUM_SECTIONS * SECTION_DIR_SIZE

    # Calculate offsets
    body_index_offset = HEADER_SIZE + section_dir_total
    chebyshev_offset = body_index_offset + body_index_size
    nutation_offset = chebyshev_offset + chebyshev_size
    delta_t_offset = nutation_offset + nut_data_size
    star_offset = delta_t_offset + delta_t_size

    total_size = star_offset + star_size

    # Allocate buffer
    buf = bytearray(total_size)

    # Write header
    header = FileHeader(
        magic=MAGIC,
        version=VERSION,
        section_count=NUM_SECTIONS,
        body_count=body_count,
        jd_start=jd_start,
        jd_end=jd_end,
        generation_epoch=now_jd,
        flags=0,
    )
    write_header(buf, header)

    # Write section directory
    sections = [
        SectionEntry(SECTION_BODY_INDEX, body_index_offset, body_index_size),
        SectionEntry(SECTION_CHEBYSHEV, chebyshev_offset, chebyshev_size),
        SectionEntry(SECTION_NUTATION, nutation_offset, nut_data_size),
        SectionEntry(SECTION_DELTA_T, delta_t_offset, delta_t_size),
        SectionEntry(SECTION_STARS, star_offset, star_size),
    ]
    for i, sec in enumerate(sections):
        write_section_dir(buf, HEADER_SIZE + i * SECTION_DIR_SIZE, sec)

    # Write body index and coefficient data
    # Asteroids may have per-body date ranges narrower than the global range.
    coeff_write_offset = chebyshev_offset
    for idx, bid in enumerate(sorted(bodies)):
        params = BODY_PARAMS[bid]
        interval_days, degree, coord_type, components = params
        seg_size = segment_byte_size(degree, components)
        n_segments = len(body_data[bid])

        # Use per-body range if available (asteroids with partial SPK coverage),
        # otherwise use global tier range.
        bid_jd_start, bid_jd_end = body_jd_ranges.get(bid, (jd_start, jd_end))

        entry = BodyEntry(
            body_id=bid,
            coord_type=coord_type,
            segment_count=n_segments,
            jd_start=bid_jd_start,
            jd_end=bid_jd_end,
            interval_days=interval_days,
            degree=degree,
            components=components,
            data_offset=coeff_write_offset,
        )
        write_body_entry(buf, body_index_offset + idx * BODY_ENTRY_SIZE, entry)

        # Write coefficient data for each segment
        for seg_coeffs in body_data[bid]:
            # seg_coeffs shape: (components, degree+1)
            # Flatten to [c0_comp0, c1_comp0, ..., cN_comp0, c0_comp1, ...]
            flat = seg_coeffs.flatten()
            struct.pack_into(f"<{len(flat)}d", buf, coeff_write_offset, *flat)
            coeff_write_offset += seg_size

    # Write nutation section
    nut_n_segments = len(nutation_coeffs)
    nut_header = NutationHeader(
        jd_start=jd_start,
        jd_end=jd_end,
        interval_days=NUTATION_INTERVAL,
        degree=NUTATION_DEGREE,
        components=NUTATION_COMPONENTS,
        segment_count=nut_n_segments,
        reserved=0,
    )
    write_nutation_header(buf, nutation_offset, nut_header)

    nut_data_offset = nutation_offset + NUTATION_HEADER_SIZE
    nut_seg_size = segment_byte_size(NUTATION_DEGREE, NUTATION_COMPONENTS)
    for seg_coeffs in nutation_coeffs:
        flat = seg_coeffs.flatten()
        struct.pack_into(f"<{len(flat)}d", buf, nut_data_offset, *flat)
        nut_data_offset += nut_seg_size

    # Write Delta-T section
    struct.pack_into(
        DELTA_T_HEADER_FMT,
        buf,
        delta_t_offset,
        len(delta_t_table),
        0,  # reserved
    )
    dt_data_offset = delta_t_offset + DELTA_T_HEADER_SIZE
    for jd, dt in delta_t_table:
        struct.pack_into(DELTA_T_ENTRY_FMT, buf, dt_data_offset, jd, dt)
        dt_data_offset += DELTA_T_ENTRY_SIZE

    # Write star catalog
    for i, star in enumerate(star_entries):
        write_star_entry(buf, star_offset + i * STAR_ENTRY_SIZE, star)

    # Write to file
    with open(output, "wb") as f:
        f.write(buf)

    if verbose:
        print("done!")
        print(f"\n  File: {output}")
        print(f"  Size: {total_size:,} bytes ({total_size / (1024 * 1024):.1f} MB)")
        print(f"  Bodies: {body_count}")
        if not skip_aux:
            print(f"  Nutation segments: {nut_n_segments}")
            print(f"  Delta-T samples: {len(delta_t_table)}")
            print(f"  Stars: {len(star_entries)}")
        print()

        # Print error summary
        # Approximate minimum geocentric distances (AU) for angular error
        # conversion. Using min distance gives the worst-case angular error.
        _MIN_GEO_DIST: dict[int, float] = {
            0: 0.98,  # Sun
            1: 0.0024,  # Moon
            2: 0.55,  # Mercury
            3: 0.26,  # Venus
            4: 0.37,  # Mars
            5: 3.9,  # Jupiter
            6: 8.0,  # Saturn
            7: 17.3,  # Uranus
            8: 28.8,  # Neptune
            9: 28.7,  # Pluto
            14: 0.0,  # Earth (geocentric = 0)
            15: 7.5,  # Chiron
            17: 1.6,  # Ceres
            18: 1.2,  # Pallas
            19: 1.0,  # Juno
            20: 1.1,  # Vesta
        }
        print("  Max fitting errors:")
        for bid in sorted(bodies):
            name = BODY_NAMES.get(bid, f"Body {bid}")
            error = body_errors[bid]
            params = BODY_PARAMS[bid]
            coord_type = params[2]
            # Show per-body range annotation if different from global
            range_note = ""
            if bid in body_jd_ranges:
                br_start, br_end = body_jd_ranges[bid]
                if abs(br_start - jd_start) > 1.0 or abs(br_end - jd_end) > 1.0:
                    yr_s = 2000.0 + (br_start - 2451545.0) / 365.25
                    yr_e = 2000.0 + (br_end - 2451545.0) / 365.25
                    range_note = f" [~{yr_s:.0f}-{yr_e:.0f}]"
            if coord_type in (COORD_ICRS_BARY, COORD_ICRS_BARY_SYSTEM):
                # Convert AU error to arcseconds using min geocentric distance
                geo_dist = _MIN_GEO_DIST.get(bid, 1.0)
                if geo_dist > 0.01:
                    arcsec = (error / geo_dist) * 206265.0
                else:
                    arcsec = error * 206265.0  # fallback for Earth
                print(f'    {name:20s}: {error:.2e} AU ({arcsec:.4f}"){range_note}')
            else:
                # ECLIPTIC, HELIO_ECL: already in degrees
                arcsec = error * 3600.0
                print(f'    {name:20s}: {error:.2e} deg ({arcsec:.4f}"){range_note}')
        if not skip_aux:
            nut_arcsec = math.degrees(nutation_error) * 3600.0
            print(f'    {"Nutation":20s}: {nutation_error:.2e} rad ({nut_arcsec:.4f}")')


# =============================================================================
# MERGE PARTIAL LEB FILES
# =============================================================================


def merge_leb_files(
    inputs: List[str],
    output: str,
    verbose: bool = True,
) -> None:
    """Merge multiple partial .leb files into a single complete file.

    Each input file must cover the same JD range but contain different bodies.
    Nutation, Delta-T, and star catalog are taken from the first input file
    that contains them.

    This allows generating body groups independently (e.g. planets, asteroids,
    analytical) and combining them afterward, which avoids the fork-deadlock
    issues of multiprocessing on macOS and gives finer control over
    regeneration.

    Args:
        inputs: List of paths to partial .leb files.
        output: Output path for the merged file.
        verbose: Print progress.

    Raises:
        ValueError: If inputs have mismatched JD ranges or overlapping bodies.
    """
    from libephemeris.leb_format import (
        read_body_entry,
        read_header,
        read_nutation_header,
        read_section_dir,
    )

    if not inputs:
        raise ValueError("No input files provided")

    if verbose:
        print(f"Merging {len(inputs)} LEB files -> {output}")

    # -------------------------------------------------------------------------
    # 1. Read all input files
    # -------------------------------------------------------------------------
    all_bodies: dict[int, tuple[str, int]] = {}  # body_id -> (source_file, idx_in_file)
    ref_jd_start: Optional[float] = None
    ref_jd_end: Optional[float] = None

    # Parsed data from each input
    input_data: List[dict] = []

    for path in inputs:
        with open(path, "rb") as f:
            data = f.read()

        hdr = read_header(data, 0)
        if hdr.magic != MAGIC:
            raise ValueError(f"Invalid LEB magic in {path}")
        if hdr.version != VERSION:
            raise ValueError(f"Unsupported LEB version {hdr.version} in {path}")

        # Validate JD range consistency
        if ref_jd_start is None:
            ref_jd_start = hdr.jd_start
            ref_jd_end = hdr.jd_end
        else:
            assert ref_jd_start is not None and ref_jd_end is not None
            if (
                abs(hdr.jd_start - ref_jd_start) > 0.5
                or abs(hdr.jd_end - ref_jd_end) > 0.5
            ):
                raise ValueError(
                    f"JD range mismatch: {path} has "
                    f"[{hdr.jd_start:.1f}, {hdr.jd_end:.1f}] but expected "
                    f"[{ref_jd_start:.1f}, {ref_jd_end:.1f}]"
                )

        # Parse sections
        sections: dict[int, SectionEntry] = {}
        for i in range(hdr.section_count):
            offset = HEADER_SIZE + i * SECTION_DIR_SIZE
            sec = read_section_dir(data, offset)
            sections[sec.section_id] = sec

        # Parse bodies
        bodies: dict[int, BodyEntry] = {}
        if SECTION_BODY_INDEX in sections:
            sec = sections[SECTION_BODY_INDEX]
            for i in range(hdr.body_count):
                off = sec.offset + i * BODY_ENTRY_SIZE
                entry = read_body_entry(data, off)
                bodies[entry.body_id] = entry

                # Check for duplicates
                if entry.body_id in all_bodies:
                    src, _ = all_bodies[entry.body_id]
                    raise ValueError(
                        f"Body {entry.body_id} ({BODY_NAMES.get(entry.body_id, '?')}) "
                        f"found in both {src} and {path}"
                    )
                all_bodies[entry.body_id] = (path, i)

        info = {
            "path": path,
            "data": data,
            "header": hdr,
            "sections": sections,
            "bodies": bodies,
        }
        input_data.append(info)

        if verbose:
            body_names = [
                BODY_NAMES.get(bid, str(bid)) for bid in sorted(bodies.keys())
            ]
            print(f"  {path}: {len(bodies)} bodies ({', '.join(body_names)})")

    assert ref_jd_start is not None and ref_jd_end is not None

    # -------------------------------------------------------------------------
    # 2. Collect body coefficient data (raw bytes)
    # -------------------------------------------------------------------------
    merged_bodies = sorted(all_bodies.keys())
    body_entries: List[BodyEntry] = []
    body_coeff_blobs: List[bytes] = []  # raw coefficient bytes per body

    for bid in merged_bodies:
        # Find the source file
        for info in input_data:
            if bid in info["bodies"]:
                entry = info["bodies"][bid]
                data = info["data"]
                seg_size = segment_byte_size(entry.degree, entry.components)
                total_bytes = entry.segment_count * seg_size
                blob = data[entry.data_offset : entry.data_offset + total_bytes]
                body_entries.append(entry)
                body_coeff_blobs.append(blob)
                break

    # -------------------------------------------------------------------------
    # 3. Collect nutation from first file that has REAL nutation data.
    # A --skip-aux partial writes an empty nutation section (header with
    # segment_count == 0); taking that one would poison the merged file.
    # -------------------------------------------------------------------------
    nutation_blob: Optional[bytes] = None
    for info in input_data:
        if SECTION_NUTATION in info["sections"]:
            sec = info["sections"][SECTION_NUTATION]
            nut_header = read_nutation_header(info["data"], sec.offset)
            if nut_header.segment_count <= 0:
                continue
            nutation_blob = info["data"][sec.offset : sec.offset + sec.size]
            break

    # -------------------------------------------------------------------------
    # 4. Collect delta-T from first file that has it
    # -------------------------------------------------------------------------
    delta_t_blob: Optional[bytes] = None
    for info in input_data:
        if SECTION_DELTA_T in info["sections"]:
            sec = info["sections"][SECTION_DELTA_T]
            delta_t_blob = info["data"][sec.offset : sec.offset + sec.size]
            break

    # -------------------------------------------------------------------------
    # 5. Collect star catalog from first file that has it
    # -------------------------------------------------------------------------
    star_blob: Optional[bytes] = None
    for info in input_data:
        if SECTION_STARS in info["sections"]:
            sec = info["sections"][SECTION_STARS]
            star_blob = info["data"][sec.offset : sec.offset + sec.size]
            break

    # -------------------------------------------------------------------------
    # 6. Calculate layout and write merged file
    # -------------------------------------------------------------------------
    body_count = len(merged_bodies)
    body_index_size = body_count * BODY_ENTRY_SIZE
    chebyshev_size = sum(len(b) for b in body_coeff_blobs)
    nut_size = len(nutation_blob) if nutation_blob else 0
    dt_size = len(delta_t_blob) if delta_t_blob else 0
    star_size = len(star_blob) if star_blob else 0

    section_dir_total = NUM_SECTIONS * SECTION_DIR_SIZE
    body_index_offset = HEADER_SIZE + section_dir_total
    chebyshev_offset = body_index_offset + body_index_size
    nutation_offset = chebyshev_offset + chebyshev_size
    delta_t_offset = nutation_offset + nut_size
    star_offset = delta_t_offset + dt_size
    total_size = star_offset + star_size

    buf = bytearray(total_size)

    # Header
    now_jd = J2000 + (time.time() / 86400.0 - 10957.5)
    header = FileHeader(
        magic=MAGIC,
        version=VERSION,
        section_count=NUM_SECTIONS,
        body_count=body_count,
        jd_start=ref_jd_start,
        jd_end=ref_jd_end,
        generation_epoch=now_jd,
        flags=0,
    )
    write_header(buf, header)

    # Section directory
    sections_list = [
        SectionEntry(SECTION_BODY_INDEX, body_index_offset, body_index_size),
        SectionEntry(SECTION_CHEBYSHEV, chebyshev_offset, chebyshev_size),
        SectionEntry(SECTION_NUTATION, nutation_offset, nut_size),
        SectionEntry(SECTION_DELTA_T, delta_t_offset, dt_size),
        SectionEntry(SECTION_STARS, star_offset, star_size),
    ]
    for i, sec in enumerate(sections_list):
        write_section_dir(buf, HEADER_SIZE + i * SECTION_DIR_SIZE, sec)

    # Body index + coefficient data
    coeff_write_offset = chebyshev_offset
    for idx, (entry, blob) in enumerate(zip(body_entries, body_coeff_blobs)):
        new_entry = BodyEntry(
            body_id=entry.body_id,
            coord_type=entry.coord_type,
            segment_count=entry.segment_count,
            jd_start=entry.jd_start,
            jd_end=entry.jd_end,
            interval_days=entry.interval_days,
            degree=entry.degree,
            components=entry.components,
            data_offset=coeff_write_offset,
        )
        write_body_entry(buf, body_index_offset + idx * BODY_ENTRY_SIZE, new_entry)
        buf[coeff_write_offset : coeff_write_offset + len(blob)] = blob
        coeff_write_offset += len(blob)

    # Nutation, Delta-T, stars (copy raw blobs)
    if nutation_blob:
        buf[nutation_offset : nutation_offset + len(nutation_blob)] = nutation_blob
    if delta_t_blob:
        buf[delta_t_offset : delta_t_offset + len(delta_t_blob)] = delta_t_blob
    if star_blob:
        buf[star_offset : star_offset + len(star_blob)] = star_blob

    with open(output, "wb") as f:
        f.write(buf)

    if verbose:
        print(f"\n  Merged file: {output}")
        print(f"  Size: {total_size:,} bytes ({total_size / (1024 * 1024):.1f} MB)")
        print(f"  Bodies: {body_count}")
        print(f"  JD range: {ref_jd_start:.1f} to {ref_jd_end:.1f}")
        body_list = [BODY_NAMES.get(b, str(b)) for b in merged_bodies]
        print(f"  Body list: {', '.join(body_list)}")
        print()


# =============================================================================
# VERIFICATION
# =============================================================================


def verify_leb(
    leb_path: str,
    n_samples: int = 500,
    verbose: bool = True,
) -> bool:
    """Post-generation validation of a .leb file.

    Compares every body in the LEB file against its reference source:
    - ICRS planets: direct Skyfield comparison (sub-arcsecond)
    - ICRS asteroids: spktype21 comparison (sub-arcsecond)
    - Ecliptic bodies: analytical function comparison (sub-arcsecond)
    - Heliocentric bodies: analytical function comparison (sub-arcsecond)

    Args:
        leb_path: Path to the .leb file.
        n_samples: Number of random JDs to test per body.
        verbose: Print progress.

    Returns:
        True if all bodies pass, False otherwise.
    """
    from libephemeris.leb_reader import LEBReader

    # Enable auto-download for asteroid SPK acquisition during verification
    import libephemeris as _ephem

    _ephem.set_auto_spk_download(True)

    reader = LEBReader(leb_path)
    jd_start, jd_end = reader.jd_range

    # Ensure asteroid SPKs cover each asteroid's actual date range in the LEB
    # (which may be narrower than the global file range for per-body coverage).
    asteroid_ids_in_file = [bid for bid in reader._bodies if bid in _ASTEROID_NAIF]
    if asteroid_ids_in_file:
        from libephemeris.minor_bodies import auto_download_asteroid_spk

        if verbose:
            print("  Preparing asteroid SPKs for verification...")
        for bid in asteroid_ids_in_file:
            body = reader._bodies[bid]
            try:
                auto_download_asteroid_spk(
                    bid, jd_start=body.jd_start, jd_end=body.jd_end, force=True
                )
            except Exception:
                pass  # Will show as FAIL if data is unavailable

    all_pass = True

    if verbose:
        print(f"Verifying {leb_path}")
        print(f"  Range: JD {jd_start:.1f} to {jd_end:.1f}")
        print(f"  Samples per body: {n_samples}")
        print()

    rng = np.random.default_rng(42)
    test_jds = rng.uniform(jd_start + 1, jd_end - 1, n_samples)

    # Build ecliptic/helio eval functions (same as used during generation)
    ecliptic_eval_funcs = _build_ecliptic_eval_funcs()
    helio_eval_funcs = _build_helio_eval_funcs()

    for body_id in sorted(reader._bodies.keys()):
        body = reader._bodies[body_id]
        name = BODY_NAMES.get(body_id, f"Body {body_id}")
        max_error = 0.0
        worst_dist = 1.0  # distance (AU) at sample with worst error
        error_unit = "AU"  # tracking what unit max_error is in

        # Use body-specific JD range if available
        body_jd_start = body.jd_start if hasattr(body, "jd_start") else jd_start
        body_jd_end = body.jd_end if hasattr(body, "jd_end") else jd_end
        body_test_jds = test_jds[
            (test_jds >= body_jd_start + 1) & (test_jds <= body_jd_end - 1)
        ]
        if len(body_test_jds) == 0:
            body_test_jds = rng.uniform(
                body_jd_start + 1, body_jd_end - 1, min(n_samples, 100)
            )

        for jd in body_test_jds:
            pos, vel = reader.eval_body(body_id, jd)

            if body.coord_type in (COORD_ICRS_BARY, COORD_ICRS_BARY_SYSTEM):
                if body_id in _PLANET_MAP:
                    # Planet: direct ICRS comparison with Skyfield
                    # For COORD_ICRS_BARY_SYSTEM, stored data is system
                    # barycenter — verification compares raw fit accuracy
                    # (COB correction is applied at runtime, not stored).
                    sample_err, dist = _verify_icrs_planet(body_id, jd, pos)
                    if sample_err > max_error:
                        max_error = sample_err
                        worst_dist = dist
                    error_unit = "AU"
                elif body_id in _ASTEROID_NAIF:
                    # Asteroid: ICRS comparison via spktype21
                    sample_err, dist = _verify_icrs_asteroid(body_id, jd, pos)
                    if sample_err > max_error:
                        max_error = sample_err
                        worst_dist = dist
                    error_unit = "AU"

            elif body.coord_type == COORD_ECLIPTIC:
                # Ecliptic body: compare with analytical function
                eval_func = ecliptic_eval_funcs.get(body_id)
                if eval_func is not None:
                    sample_err = _verify_ecliptic_body(eval_func, jd, pos)
                    if sample_err > max_error:
                        max_error = sample_err
                    error_unit = "deg"
                else:
                    # No eval function — just check range
                    if not (0.0 <= pos[0] < 360.0 or pos[0] == 0.0):
                        max_error = 999.0
                    error_unit = "deg"

            elif body.coord_type == COORD_HELIO_ECL:
                # Heliocentric body: compare with analytical function
                eval_func = helio_eval_funcs.get(body_id)
                if eval_func is not None:
                    sample_err = _verify_ecliptic_body(eval_func, jd, pos)
                    if sample_err > max_error:
                        max_error = sample_err
                    error_unit = "deg"
                else:
                    if not (0.0 <= pos[0] < 360.0 or pos[0] == 0.0):
                        max_error = 999.0
                    error_unit = "deg"

        # Convert to arcseconds and determine pass/fail
        if error_unit == "AU":
            # Convert AU error to angular error
            if worst_dist > 0.001:
                arcsec = (max_error / worst_dist) * 206265.0
            else:
                arcsec = max_error * 206265.0
            passed = arcsec < 1.0
        else:
            # Error is in degrees
            arcsec = max_error * 3600.0
            passed = arcsec < 1.0

        status = "PASS" if passed else "FAIL"
        if not passed:
            all_pass = False

        if verbose:
            # Show per-body range annotation if narrower than global
            range_note = ""
            if abs(body.jd_start - jd_start) > 1.0 or abs(body.jd_end - jd_end) > 1.0:
                yr_s = 2000.0 + (body.jd_start - 2451545.0) / 365.25
                yr_e = 2000.0 + (body.jd_end - 2451545.0) / 365.25
                range_note = f" [~{yr_s:.0f}-{yr_e:.0f}]"
            print(
                f'  {name:20s}: max error = {max_error:.2e} ({arcsec:.4f}") [{status}]{range_note}'
            )

    reader.close()

    if verbose:
        print()
        if all_pass:
            print("  ALL BODIES PASSED")
        else:
            print("  SOME BODIES FAILED")

    return all_pass


def _build_ecliptic_eval_funcs() -> dict:
    """Build evaluation functions for ecliptic bodies (used by verify_leb)."""
    from libephemeris.lunar import (
        calc_mean_lunar_node,
        calc_true_lunar_node,
        calc_mean_lilith_with_latitude,
        calc_true_lilith,
        calc_interpolated_apogee,
        calc_interpolated_perigee,
    )

    return {
        10: lambda jd: np.array([calc_mean_lunar_node(jd), 0.0, 0.0]),
        11: lambda jd: np.array(calc_true_lunar_node(jd)),
        12: lambda jd: np.array([*calc_mean_lilith_with_latitude(jd), 0.0]),
        13: lambda jd: np.array(calc_true_lilith(jd)),
        21: lambda jd: np.array(calc_interpolated_apogee(jd)),
        22: lambda jd: np.array(calc_interpolated_perigee(jd)),
    }


def _build_helio_eval_funcs() -> dict:
    """Build evaluation functions for heliocentric bodies (used by verify_leb)."""
    from libephemeris.hypothetical import calc_uranian_planet, calc_transpluto

    funcs: dict = {}
    # Uranian planets (body IDs 40-47)
    for bid in range(40, 48):
        _bid = bid  # capture for closure
        funcs[_bid] = lambda jd, b=_bid: np.array(calc_uranian_planet(b, jd)[:3])
    # Transpluto (body ID 48)
    funcs[48] = lambda jd: np.array(calc_transpluto(jd)[:3])
    return funcs


def _verify_icrs_planet(body_id: int, jd: float, leb_pos: tuple) -> Tuple[float, float]:
    """Verify an ICRS planet against Skyfield. Returns (max_err_au, dist_au)."""
    from libephemeris.state import get_planets, get_timescale
    from libephemeris.planets import get_planet_target

    planets = get_planets()
    ts = get_timescale()

    target_name = _PLANET_MAP[body_id]
    target = get_planet_target(planets, target_name)
    t = ts.tt_jd(jd)
    ref_pos = target.at(t).position.au
    sample_err = 0.0
    for c in range(3):
        err = abs(leb_pos[c] - float(ref_pos[c]))
        if err > sample_err:
            sample_err = err
    dist = math.sqrt(
        float(ref_pos[0]) ** 2 + float(ref_pos[1]) ** 2 + float(ref_pos[2]) ** 2
    )
    return sample_err, dist


def _verify_icrs_asteroid(
    body_id: int, jd: float, leb_pos: tuple
) -> Tuple[float, float]:
    """Verify an ICRS asteroid against spktype21. Returns (max_err_au, dist_au)."""
    from libephemeris.state import get_planets, get_timescale, _SPK_BODY_MAP

    planets = get_planets()
    ts = get_timescale()

    spk_info = _SPK_BODY_MAP.get(body_id)
    if spk_info is None:
        # No SPK available — report large error
        return 1.0, 1.0

    spk_file, naif_id = spk_info
    try:
        from spktype21 import SPKType21

        kernel = SPKType21.open(spk_file)
    except Exception:
        return 1.0, 1.0

    try:
        center_id = kernel.segments[0].center  # typically 10 (Sun)
        target_id = kernel.segments[0].target

        AU_KM = 149597870.7
        pos_km, _ = kernel.compute_type21(center_id, target_id, jd)

        # Get Sun barycentric position via Skyfield
        sun = planets["sun"]
        t = ts.tt_jd(jd)
        sun_bary = sun.at(t).position.au

        # helio km -> helio AU -> + Sun bary = SSB bary (ICRS)
        ref_x = pos_km[0] / AU_KM + float(sun_bary[0])
        ref_y = pos_km[1] / AU_KM + float(sun_bary[1])
        ref_z = pos_km[2] / AU_KM + float(sun_bary[2])

        err_x = abs(leb_pos[0] - ref_x)
        err_y = abs(leb_pos[1] - ref_y)
        err_z = abs(leb_pos[2] - ref_z)
        sample_err = max(err_x, err_y, err_z)
        dist = math.sqrt(ref_x**2 + ref_y**2 + ref_z**2)
        return sample_err, dist
    except Exception:
        return 1.0, 1.0
    finally:
        kernel.close()


def _verify_ecliptic_body(
    eval_func: Callable[[float], np.ndarray], jd: float, leb_pos: tuple
) -> float:
    """Verify an ecliptic/helio body against its analytical function.

    Returns max error in degrees.
    """
    ref = eval_func(jd)

    # Longitude comparison (handle 0/360 wrap)
    dlon = abs(float(leb_pos[0]) - float(ref[0]))
    if dlon > 180.0:
        dlon = 360.0 - dlon

    # Latitude comparison (no wrapping)
    dlat = abs(float(leb_pos[1]) - float(ref[1]))

    # Distance comparison (if available, convert to angular equivalent)
    # For ecliptic bodies, distance is typically small or zero
    return max(dlon, dlat)


# =============================================================================
# CLI
# =============================================================================


def _resolve_tier(args) -> Tuple[float, float, str]:
    """Resolve tier/start/end/output from CLI args.

    Returns:
        (jd_start, jd_end, output_path)
    """
    if args.tier:
        ephem_file, tier_start, tier_end, tier_output = TIER_CONFIGS[args.tier]

        from libephemeris import set_jpl_file

        set_jpl_file(ephem_file)

        start_year = args.start if args.start is not None else tier_start
        end_year = args.end if args.end is not None else tier_end

        if args.output is None:
            os.makedirs(DEFAULT_LEB_DIR, exist_ok=True)
            output = os.path.join(DEFAULT_LEB_DIR, tier_output)
        else:
            output = args.output
    else:
        if args.output is None:
            raise SystemExit("--output is required when --tier is not specified")
        output = args.output
        start_year = args.start if args.start is not None else DEFAULT_START_YEAR
        end_year = args.end if args.end is not None else DEFAULT_END_YEAR

    jd_start = _year_to_jd(start_year)
    jd_end = _year_to_jd(end_year)
    return jd_start, jd_end, output


def _group_output_path(base_output: str, group: str) -> str:
    """Derive the partial-file path for a body group.

    Example: ``data/leb/ephemeris_base.leb`` + ``planets``
             -> ``data/leb/ephemeris_base_planets.leb``
    """
    root, ext = os.path.splitext(base_output)
    return f"{root}_{group}{ext}"


def main():
    parser = argparse.ArgumentParser(
        description="Generate LEB binary ephemeris file",
        epilog=(
            "Body groups for --group:\n"
            "  planets    : Sun, Moon, Mercury-Pluto, Earth (vectorized Skyfield)\n"
            "  asteroids  : Chiron, Ceres, Pallas, Juno, Vesta (spktype21)\n"
            "  analytical : Lunar nodes, Lilith variants, Uranians, Transpluto\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output .leb file path (default: auto from --tier)",
    )
    parser.add_argument(
        "--tier",
        "-t",
        choices=["base", "medium", "extended"],
        default=None,
        help="Precision tier: base (de440s, 1850-2150), medium (de440, 1550-2650), "
        "extended (de441, -5000 to 5000). Sets ephemeris file, date range, "
        "and output path automatically.",
    )
    parser.add_argument(
        "--start",
        type=int,
        default=None,
        help=f"Start year (default: from --tier, or {DEFAULT_START_YEAR})",
    )
    parser.add_argument(
        "--end",
        type=int,
        default=None,
        help=f"End year (default: from --tier, or {DEFAULT_END_YEAR})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=os.cpu_count() or 1,
        help=f"Number of parallel workers (default: {os.cpu_count() or 1}, "
        "auto-detected CPU count)",
    )
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Run post-generation verification",
    )
    parser.add_argument(
        "--verify-samples",
        type=int,
        default=500,
        help="Number of verification samples per body (default: 500)",
    )
    parser.add_argument(
        "--bodies",
        type=str,
        default=None,
        help="Comma-separated list of body IDs or names (e.g. '1,2,3' or "
        "'moon,mercury,venus' or 'Moon,14'). Case-insensitive.",
    )
    parser.add_argument(
        "--group",
        choices=["planets", "asteroids", "exotics", "analytical"],
        default=None,
        help="Generate only a specific body group (partial file). "
        "Use --merge to combine partial files afterward.",
    )
    parser.add_argument(
        "--merge",
        nargs="+",
        metavar="FILE",
        default=None,
        help="Merge multiple partial .leb files into one. "
        "Requires --output (or --tier for auto path). "
        "Example: --merge planets.leb asteroids.leb analytical.leb",
    )
    parser.add_argument(
        "--single",
        action="store_true",
        help="Generate each body in its own subprocess (lowest memory usage). "
        "Each body runs sequentially, then all partial files are merged. "
        "Use this on memory-constrained machines.",
    )
    parser.add_argument(
        "--skip-aux",
        action="store_true",
        help="Skip nutation, Delta-T and star catalog generation. "
        "Used internally by --single mode for all bodies except the first.",
    )
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="Only verify an existing .leb file (no generation). "
        "Use with --tier to auto-resolve the file path, or --output to specify one.",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Suppress progress output",
    )

    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Mode 0: Verify-only (no generation)
    # ------------------------------------------------------------------
    if args.verify_only:
        if args.tier:
            ephem_file, _, _, tier_output = TIER_CONFIGS[args.tier]
            from libephemeris import set_jpl_file

            set_jpl_file(ephem_file)
            leb_path = (
                args.output
                if args.output
                else os.path.join(DEFAULT_LEB_DIR, tier_output)
            )
        elif args.output:
            leb_path = args.output
        else:
            parser.error("--verify-only requires --tier or --output")
            return  # unreachable

        if not os.path.exists(leb_path):
            print(f"Error: file not found: {leb_path}", file=sys.stderr)
            sys.exit(1)

        ok = verify_leb(
            leb_path,
            n_samples=args.verify_samples,
            verbose=not args.quiet,
        )
        sys.exit(0 if ok else 1)

    # ------------------------------------------------------------------
    # Mode 1: Merge existing partial files
    # ------------------------------------------------------------------
    if args.merge:
        if args.tier:
            _, _, output = _resolve_tier(args)
        elif args.output:
            output = args.output
        else:
            parser.error("--output (or --tier) is required with --merge")
            return  # unreachable

        t0 = time.time()
        merge_leb_files(args.merge, output, verbose=not args.quiet)
        elapsed = time.time() - t0
        if not args.quiet:
            print(f"  Merge time: {elapsed:.1f}s")

        if args.verify:
            print()
            ok = verify_leb(
                output,
                n_samples=args.verify_samples,
                verbose=not args.quiet,
            )
            if not ok:
                sys.exit(1)
        return

    # ------------------------------------------------------------------
    # Mode 2: Generate (full or group)
    # ------------------------------------------------------------------
    jd_start, jd_end, output = _resolve_tier(args)

    if not args.quiet and args.tier:
        ephem_file = TIER_CONFIGS[args.tier][0]
        print(f"  Tier: {args.tier} (ephemeris: {ephem_file})")

    # Resolve body list
    bodies = None
    if args.bodies:
        try:
            bodies = [_resolve_body_token(b) for b in args.bodies.split(",")]
        except ValueError as exc:
            parser.error(str(exc))
    elif args.group:
        if args.group not in BODY_GROUPS:
            parser.error(f"Unknown group: {args.group}")
        bodies = BODY_GROUPS[args.group]
        # Auto-suffix the output path when using --group without explicit --output
        if args.output is None:
            output = _group_output_path(output, args.group)
        if not args.quiet:
            print(f"  Group: {args.group} ({len(bodies)} bodies)")

    # ------------------------------------------------------------------
    # Mode 2a: Subprocess orchestration (full generation, no --group/--bodies)
    #
    # Two sub-modes:
    #   --single : one subprocess per body (lowest memory, ~31 subprocesses)
    #   default  : one subprocess per group (3 subprocesses)
    #
    # The partial files are then merged in-process and deleted.
    # ------------------------------------------------------------------
    if bodies is None and args.group is None:
        partial_files: list[str] = []

        # Build the base command shared by all subprocesses
        base_cmd = [sys.executable, os.path.abspath(__file__)]
        if args.tier:
            base_cmd += ["--tier", args.tier]
        if args.start is not None:
            base_cmd += ["--start", str(args.start)]
        if args.end is not None:
            base_cmd += ["--end", str(args.end)]
        base_cmd += ["--workers", str(args.workers)]
        if args.quiet:
            base_cmd.append("--quiet")

        t0 = time.time()

        if args.single:
            # Single-body mode: one subprocess per body, except
            # Skyfield-dependent ecliptic bodies (11, 13, 21, 22) are grouped
            # into a single subprocess to share the expensive Moon computation.
            all_bodies = []
            for group_bodies in BODY_GROUPS.values():
                all_bodies.extend(group_bodies)
            # Deduplicate preserving order
            seen: set[int] = set()
            unique_bodies: list[int] = []
            for bid in all_bodies:
                if bid not in seen:
                    seen.add(bid)
                    unique_bodies.append(bid)

            # Group Skyfield-dependent ecliptic bodies into one unit.
            # These all need (moon-earth).at(t) which is the dominant cost;
            # running them together computes it once instead of 4 times.
            _SKYFIELD_ECLIPTIC = {11, 13, 21, 22}
            units: list[tuple[list[int], str]] = []
            skyfield_group: list[int] = []
            for bid in unique_bodies:
                if bid in _SKYFIELD_ECLIPTIC:
                    skyfield_group.append(bid)
                else:
                    name = BODY_NAMES.get(bid, f"Body {bid}")
                    units.append(([bid], name))
            if skyfield_group:
                names = ", ".join(BODY_NAMES.get(b, str(b)) for b in skyfield_group)
                units.append((skyfield_group, names))

            if not args.quiet:
                print(
                    f"  Single-body mode: {len(units)} units "
                    f"({len(unique_bodies)} bodies, "
                    f"ecliptic group: {len(skyfield_group)} bodies)"
                )

            for i, (bids, name) in enumerate(units, 1):
                bid_str = ",".join(str(b) for b in bids)
                suffix = (
                    f"body{bids[0]}"
                    if len(bids) == 1
                    else f"bodies{'_'.join(str(b) for b in bids)}"
                )
                partial = _group_output_path(output, suffix)
                partial_files.append(partial)

                if not args.quiet:
                    print(f"\n[{i}/{len(units)}] Generating {name}...")

                cmd = base_cmd + [
                    "--bodies",
                    bid_str,
                    "--output",
                    partial,
                ]
                # Only the first unit generates nutation/delta-T/stars;
                # all others skip aux data (merge takes it from the first).
                if i > 1:
                    cmd.append("--skip-aux")
                result = subprocess.run(cmd)
                if result.returncode != 0:
                    print(
                        f"Error: {name} generation failed "
                        f"(exit code {result.returncode})",
                        file=sys.stderr,
                    )
                    sys.exit(result.returncode)
        else:
            # Group mode (default): one subprocess per group. Derive from
            # BODY_GROUPS so every defined group (incl. exotics) is generated
            # and merged — matching the --single branch above.
            groups = list(BODY_GROUPS.keys())

            for i, group in enumerate(groups, 1):
                partial = _group_output_path(output, group)
                partial_files.append(partial)

                if not args.quiet:
                    print(f"\n[{i}/{len(groups)}] Generating {group} group...")

                # Pass --output explicitly so the subprocess writes to the
                # exact partial path we expect (avoids auto-suffix mismatch).
                cmd = base_cmd + ["--group", group, "--output", partial]
                result = subprocess.run(cmd)
                if result.returncode != 0:
                    print(
                        f"Error: {group} group generation failed "
                        f"(exit code {result.returncode})",
                        file=sys.stderr,
                    )
                    sys.exit(result.returncode)

        # Merge partial files into final output
        if not args.quiet:
            print(f"\nMerging {len(partial_files)} partial files...")
        merge_leb_files(partial_files, output, verbose=not args.quiet)

        # Cleanup partial files
        for pf in partial_files:
            if os.path.exists(pf):
                os.remove(pf)
                if not args.quiet:
                    print(f"  Removed {os.path.basename(pf)}")

        elapsed = time.time() - t0
        if not args.quiet:
            print(f"\n  Total time: {elapsed:.1f}s")

    # ------------------------------------------------------------------
    # Mode 2b: Direct generation (--group or --bodies specified)
    # ------------------------------------------------------------------
    else:
        t0 = time.time()
        assemble_leb(
            output=output,
            jd_start=jd_start,
            jd_end=jd_end,
            bodies=bodies,
            workers=args.workers,
            verbose=not args.quiet,
            skip_aux=args.skip_aux,
        )
        elapsed = time.time() - t0

        if not args.quiet:
            print(f"  Total time: {elapsed:.1f}s")

    if args.verify:
        print()
        ok = verify_leb(
            output,
            n_samples=args.verify_samples,
            verbose=not args.quiet,
        )
        if not ok:
            sys.exit(1)


if __name__ == "__main__":
    main()
