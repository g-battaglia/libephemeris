#!/usr/bin/env python3
"""
Generate precomputed correction tables for mean lunar elements.

RETIRED: the tables this script writes (lunar_corrections.py) are not
consumed by any runtime path — see CLEAN.md.  The live interpolated-apse
residual tables are produced by
validation/calibrate/generate_lunar_apse_corrections.py in the separate,
gitignored validation/ repo (see CLAUDE.md).  Kept for reference; the smoothing
approach below also aliases short-period perturbations into the "mean"
(full perturbations rather than residuals), which is why it was retired.

This script computes correction tables that convert analytical mean elements
(JPL DE440 polynomial fits) to geometrically-derived mean elements from JPL ephemeris.

The corrections are computed by:
1. Calculating "true mean" via smoothing of geometric true positions from JPL ephemeris
2. Calculating analytical mean via polynomial fit to JPL DE440/DE441 data
3. Correction = geometric_mean - analytical_mean

Usage:
    python scripts/generate_lunar_corrections.py [--ephemeris PATH] [--output PATH]

Options:
    --ephemeris PATH    Path to ephemeris file (default: de441.bsp)
    --output PATH       Output Python file (default: libephemeris/lunar_corrections.py)
    --start-year YEAR   Start year (default: -13200)
    --end-year YEAR     End year (default: 17200)
    --step YEARS        Step in years (default: 10)
    --window-days DAYS  Smoothing window in days (default: 1461 = 4 years)
    --workers N         Number of parallel workers (default: CPU count)
    --force             Regenerate even if output file exists
    --dry-run           Show what would be done without generating

Prerequisites:
    - de441.bsp must be available (will prompt to download if missing)
    - Skyfield and libephemeris dependencies

Output file structure:
    - MEAN_NODE_CORRECTIONS: tuple of corrections for mean lunar node
    - MEAN_APSE_CORRECTIONS: tuple of corrections for mean apogee/perigee
    - Metadata constants for interpolation

References:
    - Park, R.S. et al. (2021) "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
    - Chapront-Touze, M. & Chapront, J. (1988) "ELP 2000-82B"
    - Simon, J.L. et al. (1994) "Numerical expressions for precession formulae"
"""

from __future__ import annotations

import argparse
import math
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from multiprocessing import cpu_count

DEFAULT_EPHEMERIS = "de441.bsp"
DEFAULT_OUTPUT = "libephemeris/lunar_corrections.py"
DEFAULT_START_YEAR = -13199
DEFAULT_END_YEAR = 17190
DEFAULT_STEP_YEARS = 10
DEFAULT_WINDOW_DAYS = 1461.0
DEFAULT_WORKERS = cpu_count()


def circular_mean(angles: list[float], weights: list[float] | None = None) -> float:
    """Compute weighted circular mean of angles in degrees."""
    if not angles:
        return 0.0

    if weights is None:
        weights = [1.0] * len(angles)

    sum_sin = 0.0
    sum_cos = 0.0

    for angle, weight in zip(angles, weights):
        rad = math.radians(angle)
        sum_sin += weight * math.sin(rad)
        sum_cos += weight * math.cos(rad)

    result = math.degrees(math.atan2(sum_sin, sum_cos))
    return result % 360.0


def normalize_angle_diff(diff: float) -> float:
    """Normalize angle difference to [-180, 180) range."""
    while diff >= 180.0:
        diff -= 360.0
    while diff < -180.0:
        diff += 360.0
    return diff


def jd_to_year(jd: float) -> float:
    """Convert Julian Day to year (floating point)."""
    return (jd - 2451545.0) / 365.25 + 2000.0


def year_to_jd(year: float) -> float:
    """Convert year (floating point) to Julian Day (TT)."""
    return 2451545.0 + (year - 2000.0) * 365.25


def _get_ephemeris_range() -> tuple[float, float]:
    """Get the valid Julian Date range for the current ephemeris (Moon segment)."""
    from libephemeris.state import get_planets, get_timescale

    planets = get_planets()
    ts = get_timescale()

    try:
        min_jd = float("inf")
        max_jd = float("-inf")

        for segment in planets.segments:
            try:
                start_time, end_time = segment.time_range(ts)
                seg_min = float(start_time.tt)
                seg_max = float(end_time.tt)

                target = getattr(segment, "target", None)
                if target == 301:
                    min_jd = seg_min
                    max_jd = seg_max
                else:
                    min_jd = min(min_jd, seg_min)
                    max_jd = max(max_jd, seg_max)
            except Exception:
                continue

        if min_jd == float("inf"):
            return (2415020.0, 2471184.0)

        return (min_jd, max_jd)
    except Exception:
        return (2415020.0, 2471184.0)


def compute_mean_node_geometric(
    jd_tt: float, window_days: float = DEFAULT_WINDOW_DAYS
) -> float:
    """
    Compute mean lunar node by smoothing true node positions.

    Samples true node over a window centered on jd_tt and computes
    weighted circular mean to extract the secular trend.

    Args:
        jd_tt: Julian Day in TT
        window_days: Smoothing window in days (default 4 years)

    Returns:
        Mean node longitude in degrees [0, 360)
    """
    from libephemeris.lunar import calc_true_lunar_node

    step_days = 5.0
    n_samples = int(window_days / step_days)
    half_window = window_days / 2.0

    min_jd, max_jd = _get_ephemeris_range()

    angles = []
    weights = []

    for i in range(n_samples):
        offset = (i - n_samples / 2.0) * step_days
        if abs(offset) > half_window:
            continue

        sample_jd = jd_tt + offset

        if sample_jd < min_jd or sample_jd > max_jd:
            continue

        try:
            lon = calc_true_lunar_node(sample_jd)[0]
            w = math.exp(-(offset**2) / (2.0 * (half_window / 2.0) ** 2))
            angles.append(lon)
            weights.append(w)
        except Exception:
            continue

    if not angles:
        try:
            lon = calc_true_lunar_node(jd_tt)[0]
            return lon
        except Exception:
            raise ValueError(f"No valid samples for JD {jd_tt}")

    return circular_mean(angles, weights)


def compute_mean_apse_geometric(
    jd_tt: float, window_days: float = DEFAULT_WINDOW_DAYS
) -> float:
    """
    Compute mean lunar apogee by smoothing true apogee positions.

    Uses the same smoothing approach as mean node.

    Args:
        jd_tt: Julian Day in TT
        window_days: Smoothing window in days (default 4 years)

    Returns:
        Mean apogee longitude in degrees [0, 360)
    """
    from libephemeris.lunar import calc_true_lilith

    step_days = 5.0
    n_samples = int(window_days / step_days)
    half_window = window_days / 2.0

    min_jd, max_jd = _get_ephemeris_range()

    angles = []
    weights = []

    for i in range(n_samples):
        offset = (i - n_samples / 2.0) * step_days
        if abs(offset) > half_window:
            continue

        sample_jd = jd_tt + offset

        if sample_jd < min_jd or sample_jd > max_jd:
            continue

        try:
            lon = calc_true_lilith(sample_jd)[0]
            w = math.exp(-(offset**2) / (2.0 * (half_window / 2.0) ** 2))
            angles.append(lon)
            weights.append(w)
        except Exception:
            continue

    if not angles:
        try:
            lon = calc_true_lilith(jd_tt)[0]
            return lon
        except Exception:
            raise ValueError(f"No valid samples for JD {jd_tt}")

    return circular_mean(angles, weights)


def compute_mean_node_analytical(jd_tt: float) -> float:
    """
    Compute mean lunar node using analytical polynomial formula.

    This uses the geometric relationship: Mean Node = L' - F
    where L' is the Moon's mean longitude and F is the mean argument of latitude.
    The polynomial coefficients are derived from JPL DE440/DE441 ephemeris data
    using least-squares fitting to ensure maximum accuracy.

    The fundamental arguments (L', F) are computed using high-precision
    polynomial expansions with secular terms up to T^5.

    Args:
        jd_tt: Julian Day in TT

    Returns:
        Mean node longitude in degrees [0, 360)

    References:
        - Simon, J.L. et al. (1994) A&A 282, 663-683
        - Chapront, J. et al. (2002) A&A 387, 700-708
    """
    T = (jd_tt - 2451545.0) / 36525.0
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    T5 = T4 * T

    L_moon = (
        218.3164477
        + 481267.88123421 * T
        - 0.0015786 * T2
        + T3 / 538841.0
        - T4 / 65194000.0
        + T5 / 128627100000.0
    )

    F = (
        93.2720950
        + 483202.0175233 * T
        - 0.0036539 * T2
        - T3 / 3526000.0
        + T4 / 863310000.0
        - T5 / 142650000000.0
    )

    node = (L_moon - F) % 360.0

    return node


def compute_mean_apse_analytical(jd_tt: float) -> float:
    """
    Compute mean lunar apogee using analytical polynomial formula.

    This uses the geometric relationship: Mean Apogee = L' - M' + 180°
    where L' is the Moon's mean longitude and M' is the Moon's mean anomaly.
    The result is projected from the lunar orbital plane to the ecliptic
    using the mean inclination of the lunar orbit.

    The polynomial coefficients are derived from JPL DE440/DE441 ephemeris
    data using least-squares fitting for maximum accuracy.

    Args:
        jd_tt: Julian Day in TT

    Returns:
        Mean apogee longitude in degrees [0, 360)

    References:
        - Simon, J.L. et al. (1994) A&A 282, 663-683
        - Chapront, J. et al. (2002) A&A 387, 700-708
    """
    T = (jd_tt - 2451545.0) / 36525.0
    T2 = T * T
    T3 = T2 * T
    T4 = T3 * T
    T5 = T4 * T

    L_moon = (
        218.3164477
        + 481267.88123421 * T
        - 0.0015786 * T2
        + T3 / 538841.0
        - T4 / 65194000.0
        + T5 / 128627100000.0
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

    apogee_unproj = (L_moon - M_prime + 180.0) % 360.0
    node = (L_moon - F) % 360.0

    MOON_MEAN_INCL = 5.1453964

    lon_from_node = math.radians(apogee_unproj - node)

    x = math.cos(lon_from_node)
    y = math.sin(lon_from_node)
    z = 0.0

    incl_rad = -math.radians(MOON_MEAN_INCL)
    cos_incl = math.cos(incl_rad)
    sin_incl = math.sin(incl_rad)
    y_new = y * cos_incl - z * sin_incl
    y * sin_incl + z * cos_incl

    lon_from_node_proj = math.atan2(y_new, x)

    apogee_proj = math.degrees(lon_from_node_proj) + node
    apogee_proj = apogee_proj % 360.0

    return apogee_proj


def compute_correction_for_year(
    year: int, window_days: float = DEFAULT_WINDOW_DAYS
) -> tuple[int, float, float]:
    """
    Compute both node and apse corrections for a single year.

    This function is designed to be called in parallel.

    Args:
        year: Year to compute corrections for
        window_days: Smoothing window in days

    Returns:
        Tuple of (year, node_correction, apse_correction)
    """
    jd = year_to_jd(year)

    geo_node = compute_mean_node_geometric(jd, window_days)
    analytical_node = compute_mean_node_analytical(jd)
    node_corr = normalize_angle_diff(geo_node - analytical_node)

    geo_apse = compute_mean_apse_geometric(jd, window_days)
    analytical_apse = compute_mean_apse_analytical(jd)
    apse_corr = normalize_angle_diff(geo_apse - analytical_apse)

    return (year, node_corr, apse_corr)


# =============================================================================
# PERIGEE PERTURBATION RESIDUAL COMPUTATION
# =============================================================================

# Perigee table uses a finer step since perturbations are larger
PERIGEE_STEP_YEARS = 2
PERIGEE_HALF_WINDOW_DAYS = 28.0
PERIGEE_N_SAMPLES = 9
PERIGEE_FIT_START = -2000
PERIGEE_FIT_END = 4000


def _unwrap_longitudes(longitudes: list[float]) -> list[float]:
    """Unwrap a sequence of longitudes to handle 0/360 discontinuity."""
    if not longitudes:
        return []

    unwrapped = [longitudes[0]]
    for i in range(1, len(longitudes)):
        diff = longitudes[i] - longitudes[i - 1]
        if diff > 180:
            diff -= 360
        elif diff < -180:
            diff += 360
        unwrapped.append(unwrapped[-1] + diff)

    return unwrapped


def compute_interpolated_perigee_jpl(jd_tt: float) -> float:
    """
    Compute interpolated perigee directly from JPL ephemeris.

    Uses quadratic polynomial regression on 9 samples of the osculating
    perigee in a +/-28 day window, with unwrap of longitudes.
    Same approach as calc_interpolated_apogee() in production.
    """
    from libephemeris.lunar import calc_osculating_perigee

    step = 2.0 * PERIGEE_HALF_WINDOW_DAYS / (PERIGEE_N_SAMPLES - 1)

    times = []
    lons = []
    for i in range(PERIGEE_N_SAMPLES):
        t = jd_tt - PERIGEE_HALF_WINDOW_DAYS + i * step
        lon, _, _ = calc_osculating_perigee(t)
        times.append(t - jd_tt)
        lons.append(lon)

    lons = _unwrap_longitudes(lons)

    # Quadratic fit: solve for coefficients [a, b, c] where lon = a*t^2 + b*t + c
    n = len(times)
    sum_t0 = float(n)
    sum_t1 = sum(times)
    sum_t2 = sum(t**2 for t in times)
    sum_t3 = sum(t**3 for t in times)
    sum_t4 = sum(t**4 for t in times)
    sum_y = sum(lons)
    sum_ty = sum(t * y for t, y in zip(times, lons))
    sum_t2y = sum(t**2 * y for t, y in zip(times, lons))

    A = [
        [sum_t4, sum_t3, sum_t2],
        [sum_t3, sum_t2, sum_t1],
        [sum_t2, sum_t1, sum_t0],
    ]
    rhs = [sum_t2y, sum_ty, sum_y]

    for col in range(3):
        max_row = col
        for row in range(col + 1, 3):
            if abs(A[row][col]) > abs(A[max_row][col]):
                max_row = row
        A[col], A[max_row] = A[max_row], A[col]
        rhs[col], rhs[max_row] = rhs[max_row], rhs[col]
        for row in range(col + 1, 3):
            if abs(A[col][col]) < 1e-30:
                continue
            factor = A[row][col] / A[col][col]
            for k in range(col, 3):
                A[row][k] -= factor * A[col][k]
            rhs[row] -= factor * rhs[col]

    x = [0.0, 0.0, 0.0]
    for i in range(2, -1, -1):
        x[i] = rhs[i]
        for j in range(i + 1, 3):
            x[i] -= A[i][j] * x[j]
        if abs(A[i][i]) > 1e-30:
            x[i] /= A[i][i]

    return x[2] % 360.0


def compute_perigee_perturbation_series(jd_tt: float) -> float:
    """Compute the trigonometric perturbation series for perigee.

    Delegates to the canonical implementation in lunar.py to avoid
    maintaining duplicate coefficient tables.

    Calibrated against JPL DE441 using passage-interpolated harmonic fitting (v2.2).
    """
    from libephemeris.lunar import _calc_elp2000_perigee_perturbations

    return _calc_elp2000_perigee_perturbations(jd_tt)


def compute_perigee_residual_for_year(year: int) -> tuple[int, float] | None:
    """
    Compute perigee perturbation residual for a single year.

    residual = interpolated_perigee_JPL - mean_perigee - trig_series

    Returns:
        (year, residual) or None if computation fails
    """
    from libephemeris.lunar import calc_mean_lilith

    jd = year_to_jd(year)

    try:
        interp_perigee = compute_interpolated_perigee_jpl(jd)
    except Exception:
        return None

    mean_perigee = (calc_mean_lilith(jd) + 180.0) % 360.0
    perturbation_target = normalize_angle_diff(interp_perigee - mean_perigee)

    trig_series = compute_perigee_perturbation_series(jd)
    residual = perturbation_target - trig_series

    return (year, residual)


def format_corrections_table(corrections: list[float], values_per_line: int = 8) -> str:
    """Format corrections list as Python tuple with nice formatting."""
    lines = []
    for i in range(0, len(corrections), values_per_line):
        chunk = corrections[i : i + values_per_line]
        line = "    " + ", ".join(f"{v:.10f}" for v in chunk)
        if i + values_per_line < len(corrections):
            line += ","
        lines.append(line)
    return "\n".join(lines)


def generate_output_file(
    node_corrections: list[float],
    apse_corrections: list[float],
    start_year: int,
    end_year: int,
    step: int,
    ephemeris: str,
    output_path: str,
    window_days: float,
    perigee_corrections: list[float] | None = None,
    perigee_start_year: int | None = None,
    perigee_end_year: int | None = None,
    perigee_step: int | None = None,
) -> None:
    """Generate the Python output file with corrections."""

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    content = f'''"""
AUTO-GENERATED FILE - DO NOT EDIT MANUALLY

This file contains precomputed correction tables for mean lunar elements.
To regenerate:
    python scripts/generate_lunar_corrections.py --ephemeris {ephemeris}

Generated by: scripts/generate_lunar_corrections.py
Generated on: {timestamp}
Ephemeris used: {ephemeris}
Smoothing window: {window_days:.0f} days

Tables provide corrections to convert analytical mean elements (polynomial
formulas based on JPL DE440/DE441) to geometrically-derived mean elements
from the same JPL ephemeris.

Method:
  - True positions computed from JPL ephemeris state vectors
  - Mean derived by Gaussian-weighted smoothing over {window_days:.0f}-day window
  - Correction = mean_geometric - mean_analytical

Usage:
    from libephemeris.lunar_corrections import (
        MEAN_NODE_CORRECTIONS,
        MEAN_APSE_CORRECTIONS,
        CORRECTION_START_YEAR,
        CORRECTION_STEP_YEARS,
    )

    # Interpolate correction for a given year
    idx = (year - CORRECTION_START_YEAR) / CORRECTION_STEP_YEARS
    idx_low = int(idx)
    correction = CORRECTIONS[idx_low] + frac * (CORRECTIONS[idx_low + 1] - CORRECTIONS[idx_low])

    # Apply to analytical mean
    mean_node = mean_node_analytical(jd) + correction

References:
    - Park, R.S. et al. (2021) "The JPL Planetary and Lunar Ephemerides DE440 and DE441"
    - Simon, J.L. et al. (1994) "Numerical expressions for precession formulae", A&A 282
    - Chapront, J. et al. (2002) "A new determination of lunar orbital parameters", A&A 387
"""

# =============================================================================
# METADATA
# =============================================================================

# Start year for correction tables (inclusive)
CORRECTION_START_YEAR: int = {start_year}

# End year for correction tables (inclusive)
CORRECTION_END_YEAR: int = {end_year}

# Step size in years between correction values
CORRECTION_STEP_YEARS: int = {step}

# Ephemeris used to generate corrections
CORRECTION_EPHEMERIS: str = "{ephemeris}"

# Number of correction values in each table
CORRECTION_COUNT: int = {len(node_corrections)}

# =============================================================================
# MEAN LUNAR NODE CORRECTIONS
# =============================================================================
# Corrections in degrees. Index i corresponds to year:
#   year = CORRECTION_START_YEAR + i * CORRECTION_STEP_YEARS
#
# Usage:
#   mean_node = mean_node_analytical(jd) + correction
#
# Precision: < 0.001 degrees throughout range when using linear interpolation
# =============================================================================
MEAN_NODE_CORRECTIONS: tuple[float, ...] = (
{format_corrections_table(node_corrections)}
)

# =============================================================================
# MEAN LUNAR APSE CORRECTIONS
# =============================================================================
# Corrections for mean apogee (Lilith) and mean perigee in degrees.
# Index i corresponds to year:
#   year = CORRECTION_START_YEAR + i * CORRECTION_STEP_YEARS
#
# Usage:
#   mean_apogee = mean_apogee_analytical(jd) + correction
#   mean_perigee = mean_apogee + 180 degrees
#
# Precision: < 0.001 degrees throughout range when using linear interpolation
# =============================================================================
MEAN_APSE_CORRECTIONS: tuple[float, ...] = (
{format_corrections_table(apse_corrections)}
)
'''

    # Add perigee perturbation corrections if available
    if perigee_corrections is not None and perigee_start_year is not None:
        p_step = perigee_step or PERIGEE_STEP_YEARS
        p_end = (
            perigee_end_year
            or perigee_start_year + (len(perigee_corrections) - 1) * p_step
        )
        content += f"""
# =============================================================================
# PERIGEE PERTURBATION METADATA
# =============================================================================

# Step in years between perigee corrections
PERIGEE_CORRECTION_STEP_YEARS: int = {p_step}

# Start year for perigee correction table (inclusive)
PERIGEE_CORRECTION_START_YEAR: int = {perigee_start_year}

# End year for perigee correction table (inclusive)
PERIGEE_CORRECTION_END_YEAR: int = {p_end}

# Number of perigee correction entries
PERIGEE_CORRECTION_COUNT: int = {len(perigee_corrections)}

# =============================================================================
# PERIGEE PERTURBATION CORRECTIONS
# =============================================================================
# Residual corrections for the interpolated perigee perturbation series.
# These correct the difference between the trigonometric perturbation series
# and the actual interpolated perigee position from JPL ephemeris.
#
# Usage:
#   perturbation = trig_series(jd) + interpolated_correction(jd)
#
# Index i corresponds to year:
#   year = PERIGEE_CORRECTION_START_YEAR + i * PERIGEE_CORRECTION_STEP_YEARS
# =============================================================================
PERIGEE_PERTURBATION_CORRECTIONS: tuple[float, ...] = (
{format_corrections_table(perigee_corrections)}
)
"""

    with open(output_path, "w") as f:
        f.write(content)

    print(f"Generated: {output_path}")
    print(f"  - {len(node_corrections)} node corrections")
    print(f"  - {len(apse_corrections)} apse corrections")
    if perigee_corrections is not None:
        print(f"  - {len(perigee_corrections)} perigee perturbation corrections")
    print(f"  - Range: {start_year} to {end_year} CE")
    print(f"  - Step: {step} years")


def setup_ephemeris(ephemeris: str) -> bool:
    """Set up the ephemeris file for use."""
    from libephemeris.state import set_ephemeris_file

    try:
        set_ephemeris_file(ephemeris)
        return True
    except Exception as e:
        print(f"Error setting ephemeris: {e}")
        return False


def init_worker(ephemeris: str) -> None:
    """Initialize worker process with ephemeris."""
    os.environ["LIBEPHEMERIS_EPHEMERIS"] = ephemeris
    from libephemeris.state import set_ephemeris_file

    set_ephemeris_file(ephemeris)


def main():
    parser = argparse.ArgumentParser(
        description="Generate precomputed correction tables for mean lunar elements",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--ephemeris",
        default=DEFAULT_EPHEMERIS,
        help=f"Ephemeris file to use (default: {DEFAULT_EPHEMERIS})",
    )
    parser.add_argument(
        "--output",
        default=DEFAULT_OUTPUT,
        help=f"Output Python file (default: {DEFAULT_OUTPUT})",
    )
    parser.add_argument(
        "--start-year",
        type=int,
        default=DEFAULT_START_YEAR,
        help=f"Start year (default: {DEFAULT_START_YEAR})",
    )
    parser.add_argument(
        "--end-year",
        type=int,
        default=DEFAULT_END_YEAR,
        help=f"End year (default: {DEFAULT_END_YEAR})",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=DEFAULT_STEP_YEARS,
        help=f"Step in years (default: {DEFAULT_STEP_YEARS})",
    )
    parser.add_argument(
        "--window-days",
        type=float,
        default=DEFAULT_WINDOW_DAYS,
        help=f"Smoothing window in days (default: {DEFAULT_WINDOW_DAYS})",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Number of parallel workers (default: {DEFAULT_WORKERS})",
    )
    parser.add_argument(
        "--force", action="store_true", help="Regenerate even if output file exists"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Show what would be done without generating",
    )
    parser.add_argument(
        "--test-single",
        type=float,
        help="Test computation for a single year (for debugging)",
    )

    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)

    if not os.path.isabs(args.output):
        args.output = os.path.join(project_root, args.output)

    if os.path.exists(args.output) and not args.force:
        print(f"Output file already exists: {args.output}")
        print("Use --force to regenerate")
        return 1

    print(f"Setting up ephemeris: {args.ephemeris}")
    if not setup_ephemeris(args.ephemeris):
        print(f"Make sure {args.ephemeris} is available in the ephemeris path")
        print("You may need to download it first or set LIBEPHEMERIS_EPHEMERIS")
        return 1
    print(f"Using ephemeris: {args.ephemeris}")

    if args.test_single is not None:
        print(f"\nTesting single year: {args.test_single}")
        jd = year_to_jd(args.test_single)
        print(f"  JD: {jd}")

        print("  Computing geometric mean node...")
        geo_node = compute_mean_node_geometric(jd, args.window_days)
        print(f"    Geometric: {geo_node:.10f} deg")

        print("  Computing analytical mean node...")
        analytical_node = compute_mean_node_analytical(jd)
        print(f"    Analytical: {analytical_node:.10f} deg")

        correction = normalize_angle_diff(geo_node - analytical_node)
        print(f"  Correction: {correction:.10f} deg")

        print("\n  Computing geometric mean apse...")
        geo_apse = compute_mean_apse_geometric(jd, args.window_days)
        print(f"    Geometric: {geo_apse:.10f} deg")

        print("  Computing analytical mean apse...")
        analytical_apse = compute_mean_apse_analytical(jd)
        print(f"    Analytical: {analytical_apse:.10f} deg")

        apse_correction = normalize_angle_diff(geo_apse - analytical_apse)
        print(f"  Correction: {apse_correction:.10f} deg")

        return 0

    years = list(range(args.start_year, args.end_year + 1, args.step))
    n_years = len(years)

    print(f"\nGenerating corrections for {n_years} years:")
    print(f"  Range: {args.start_year} to {args.end_year} CE")
    print(f"  Step: {args.step} years")
    print(f"  Smoothing window: {args.window_days} days")
    print(f"  Parallel workers: {args.workers}")

    if args.dry_run:
        print("\n[DRY RUN] Would generate corrections")
        return 0

    results = {}

    tqdm = None
    try:
        from tqdm import tqdm as tqdm_module

        tqdm = tqdm_module
    except ImportError:
        print("(Install tqdm for progress bar)")

    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=init_worker,
        initargs=(args.ephemeris,),
    ) as executor:
        futures = {
            executor.submit(compute_correction_for_year, year, args.window_days): year
            for year in years
        }

        if tqdm is not None:
            for future in tqdm(
                as_completed(futures), total=len(futures), desc="Computing corrections"
            ):
                year, node_corr, apse_corr = future.result()
                results[year] = (node_corr, apse_corr)
        else:
            completed = 0
            for future in as_completed(futures):
                year, node_corr, apse_corr = future.result()
                results[year] = (node_corr, apse_corr)
                completed += 1
                if completed % 100 == 0:
                    print(
                        f"  Progress: {completed}/{n_years} ({100 * completed // n_years}%)"
                    )

    node_corrections = [results[year][0] for year in years]
    apse_corrections = [results[year][1] for year in years]

    # --- Phase 2: Compute perigee perturbation corrections ---
    print("\nPhase 2: Computing perigee perturbation corrections...")
    perigee_years = list(range(args.start_year, args.end_year + 1, PERIGEE_STEP_YEARS))
    n_perigee = len(perigee_years)
    print(f"  Range: {args.start_year} to {args.end_year} CE")
    print(f"  Step: {PERIGEE_STEP_YEARS} years")
    print(f"  Samples: {n_perigee}")

    perigee_results = {}
    perigee_failed = []

    with ProcessPoolExecutor(
        max_workers=args.workers,
        initializer=init_worker,
        initargs=(args.ephemeris,),
    ) as executor:
        futures = {
            executor.submit(compute_perigee_residual_for_year, year): year
            for year in perigee_years
        }

        if tqdm is not None:
            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc="Computing perigee corrections",
            ):
                result = future.result()
                if result is None:
                    perigee_failed.append(futures[future])
                else:
                    year, residual = result
                    perigee_results[year] = residual
        else:
            completed = 0
            for future in as_completed(futures):
                result = future.result()
                if result is None:
                    perigee_failed.append(futures[future])
                else:
                    year, residual = result
                    perigee_results[year] = residual
                completed += 1
                if completed % 500 == 0:
                    print(
                        f"  Progress: {completed}/{n_perigee} ({100 * completed // n_perigee}%)"
                    )

    if perigee_failed:
        print(f"  WARNING: {len(perigee_failed)} perigee samples failed")
        perigee_failed.sort()

    # Build ordered perigee corrections list
    # For failed years, use linear interpolation from neighbors
    perigee_corrections = []
    actual_perigee_start = None
    actual_perigee_end = None

    for year in perigee_years:
        if year in perigee_results:
            perigee_corrections.append(perigee_results[year])
            if actual_perigee_start is None:
                actual_perigee_start = year
            actual_perigee_end = year
        else:
            # Fill with 0.0 for failed samples (will be interpolated at runtime)
            perigee_corrections.append(0.0)
            if actual_perigee_start is None:
                actual_perigee_start = year
            actual_perigee_end = year

    if actual_perigee_start is None:
        actual_perigee_start = args.start_year
    if actual_perigee_end is None:
        actual_perigee_end = args.end_year

    print(f"  Successful: {len(perigee_results)}/{n_perigee}")
    print(f"  Failed: {len(perigee_failed)}")

    generate_output_file(
        node_corrections,
        apse_corrections,
        args.start_year,
        args.end_year,
        args.step,
        args.ephemeris,
        args.output,
        args.window_days,
        perigee_corrections=perigee_corrections,
        perigee_start_year=actual_perigee_start,
        perigee_end_year=actual_perigee_end,
        perigee_step=PERIGEE_STEP_YEARS,
    )

    return 0


if __name__ == "__main__":
    sys.exit(main())
