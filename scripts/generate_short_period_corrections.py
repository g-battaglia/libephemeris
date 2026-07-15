#!/usr/bin/env python3
"""Generate empirical short-period perturbation correction coefficients.

Provenance:
    NON-PRODUCTION, REJECTED EXPERIMENT. Inputs are public NASA/JPL SPK states
    and the project's cited Keplerian fallback, but the residual fit is not a
    primary-source orbital theory. Its output is not shipped, imported, or
    accepted as runtime provenance. The script is retained to make the rejected
    methodology reproducible and to prevent it being rediscovered silently;
    ``docs/development/roadmap.md`` records why numerical integration won.

Status:
    Output from this script is neither packaged nor accepted as a runtime
    model. Reintroducing any generated coefficient requires a new, independently
    sourced scientific design and review of the provenance registry.

Fits Fourier series to the DETRENDED residuals (SPK - Keplerian) for each
major asteroid, capturing the dominant short-period perturbation terms from
Jupiter, Saturn, and (for centaurs) Uranus.

Key insight: residuals over centuries are dominated by secular drift from
imperfect secular perturbation theory. Polynomial detrending removes this
secular component, leaving the periodic oscillations that Fourier terms
can capture effectively.

The corrections are applied as delta_lon, delta_lat, delta_dist on top of
the existing Keplerian propagation with secular perturbations and multi-epoch
elements.

Usage:
    python scripts/generate_short_period_corrections.py

Output:
    Prints Python code for SHORT_PERIOD_CORRECTIONS dict to paste into
    minor_bodies.py, and writes it to data/short_period_corrections.py.
"""

from __future__ import annotations

import math
import os
import sys
import time
from pathlib import Path

import numpy as np

# Add project root to path
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from libephemeris.constants import (
    CERES,
    CHIRON,
    JUNO,
    PALLAS,
    PHOLUS,
    VESTA,
)
from libephemeris.minor_bodies import (
    JUPITER_N,
    MINOR_BODY_ELEMENTS,
    SATURN_N,
    URANUS_N,
    _get_closest_epoch_elements,
    calc_minor_body_position,
)

# =============================================================================
# CONSTANTS
# =============================================================================

# Reference epoch: J2000.0
J2000_JD = 2451545.0

# Mean anomalies at J2000.0 (degrees)
# Source: Simon et al. (1994), consistent with JPL DE440
JUPITER_M0_J2000 = 20.020
SATURN_M0_J2000 = 317.021
URANUS_M0_J2000 = 142.238

# SPK file directory
SPK_DIR = Path(os.path.expanduser("~/.libephemeris/spk"))

# Bodies to fit — (body_id, name, spk_prefix)
BODIES = [
    (CERES, "Ceres", "ceres"),
    (PALLAS, "Pallas", "pallas"),
    (JUNO, "Juno", "juno"),
    (VESTA, "Vesta", "vesta"),
    (CHIRON, "Chiron", "2060"),
    (PHOLUS, "Pholus", "5145"),
]

# Fitting parameters
SAMPLE_INTERVAL_DAYS = 5.0  # days between samples
# Use the wide-range SPK coverage (~1600-2500) but leave margins
FIT_START_JD = 2325020.5  # ~1653 CE
FIT_END_JD = 2614500.5  # ~2445 CE

# Detrending polynomial degree — captures secular drift without
# absorbing the periodic signal. Degree 5 handles up to quintic
# secular error growth over 800 years.
DETREND_DEGREE = 5

# Amplitude thresholds for pruning — keep terms above these
# Use lower thresholds since detrending isolates the periodic signal
AMPLITUDE_THRESHOLD_LON_ARCSEC = 2.0
AMPLITUDE_THRESHOLD_LAT_ARCSEC = 1.0
AMPLITUDE_THRESHOLD_DIST_AU = 0.00002

# Maximum number of terms per coordinate per body (to keep dict compact)
MAX_TERMS_PER_COORD = 40

# Obliquity of J2000 ecliptic (degrees)
OBLIQUITY_J2000 = 23.4392911


# =============================================================================
# SPK REFERENCE POSITIONS
# =============================================================================


def _find_spk_file(spk_prefix: str) -> Path | None:
    """Find the wide-range SPK file for a body."""
    pattern = f"{spk_prefix}_160001_250001.bsp"
    path = SPK_DIR / pattern
    if path.exists():
        return path
    # Try other wide-range patterns
    for f in sorted(SPK_DIR.glob(f"{spk_prefix}_*.bsp")):
        return f
    return None


def _compute_spk_positions_batch(
    spk_file: Path, jd_array: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute heliocentric ecliptic positions from SPK for all JDs.

    Returns (lon_deg, lat_deg, dist_au, valid_mask) arrays.
    """
    from spktype21 import SPKType21

    n = len(jd_array)
    lon = np.zeros(n)
    lat = np.zeros(n)
    dist = np.zeros(n)
    valid = np.ones(n, dtype=bool)

    eps = math.radians(OBLIQUITY_J2000)
    cos_eps = math.cos(eps)
    sin_eps = math.sin(eps)
    AU_KM = 149597870.7

    kernel = SPKType21.open(str(spk_file))
    try:
        center_id = kernel.segments[0].center
        target_id = kernel.segments[0].target

        for i in range(n):
            try:
                pos_km, _ = kernel.compute_type21(center_id, target_id, jd_array[i])

                # ICRS -> ecliptic J2000
                x_eq = pos_km[0] / AU_KM
                y_eq = pos_km[1] / AU_KM
                z_eq = pos_km[2] / AU_KM

                x_ecl = x_eq
                y_ecl = y_eq * cos_eps + z_eq * sin_eps
                z_ecl = -y_eq * sin_eps + z_eq * cos_eps

                r = math.sqrt(x_ecl**2 + y_ecl**2 + z_ecl**2)
                lon[i] = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
                lat[i] = math.degrees(math.asin(max(-1, min(1, z_ecl / r))))
                dist[i] = r
            except Exception:
                valid[i] = False
    finally:
        kernel.close()

    return lon, lat, dist, valid


def _compute_keplerian_positions_batch(
    body_id: int, jd_array: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute heliocentric ecliptic positions from Keplerian propagation."""
    n = len(jd_array)
    lon = np.zeros(n)
    lat = np.zeros(n)
    dist = np.zeros(n)

    for i in range(n):
        elements = _get_closest_epoch_elements(body_id, jd_array[i])
        x, y, z = calc_minor_body_position(elements, jd_array[i], body_id=body_id)
        r = math.sqrt(x**2 + y**2 + z**2)
        lon[i] = math.degrees(math.atan2(y, x)) % 360.0
        lat[i] = math.degrees(math.asin(max(-1, min(1, z / r)))) if r > 0 else 0.0
        dist[i] = r

    return lon, lat, dist


# =============================================================================
# ANGULAR VARIABLES
# =============================================================================


def _compute_mean_anomalies(
    body_id: int, jd_array: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Compute mean anomalies for asteroid, Jupiter, Saturn, Uranus.

    Returns (M_ast, M_J, M_S, M_U) arrays in degrees.
    """
    n = len(jd_array)
    M_ast = np.zeros(n)

    for i in range(n):
        elements = _get_closest_epoch_elements(body_id, jd_array[i])
        dt = jd_array[i] - elements.epoch
        M_ast[i] = (elements.M0 + elements.n * dt) % 360.0

    dt_j2000 = jd_array - J2000_JD
    M_J = (JUPITER_M0_J2000 + JUPITER_N * dt_j2000) % 360.0
    M_S = (SATURN_M0_J2000 + SATURN_N * dt_j2000) % 360.0
    M_U = (URANUS_M0_J2000 + URANUS_N * dt_j2000) % 360.0

    return M_ast, M_J, M_S, M_U


# =============================================================================
# FOURIER BASIS
# =============================================================================


def _build_basis_labels(include_uranus: bool = False) -> list[tuple[int, int, str]]:
    """Build the list of (j, k, planet) basis function labels.

    Each label identifies a term: angle = j*M_ast + k*M_planet.
    For planet='A' (asteroid only), k is ignored (always 0 in the angle).
    """
    labels = []

    # Pure asteroid harmonics: j*M_ast, j = 1..6
    for j in range(1, 7):
        labels.append((j, 0, "A"))

    # Jupiter coupling: j*M_ast + k*M_J
    for j in range(-5, 6):
        for k in range(-4, 5):
            if k == 0:
                continue
            labels.append((j, k, "J"))

    # Saturn coupling: j*M_ast + k*M_S
    for j in range(-4, 5):
        for k in range(-3, 4):
            if k == 0:
                continue
            labels.append((j, k, "S"))

    # Uranus coupling (for centaurs)
    if include_uranus:
        for j in range(-3, 4):
            for k in range(-2, 3):
                if k == 0:
                    continue
                labels.append((j, k, "U"))

    return labels


def _build_design_matrix(
    labels: list[tuple[int, int, str]],
    M_ast: np.ndarray,
    M_J: np.ndarray,
    M_S: np.ndarray,
    M_U: np.ndarray,
) -> np.ndarray:
    """Build the design matrix for least-squares fitting.

    Each basis function contributes 2 columns: cos(angle) and sin(angle).
    Returns matrix of shape (n_samples, 2 * n_terms).
    """
    n = len(M_ast)
    n_terms = len(labels)
    X = np.zeros((n, 2 * n_terms))

    M_ast_rad = np.radians(M_ast)
    M_J_rad = np.radians(M_J)
    M_S_rad = np.radians(M_S)
    M_U_rad = np.radians(M_U)

    for idx, (j, k, planet) in enumerate(labels):
        if planet == "A":
            angle = j * M_ast_rad
        elif planet == "J":
            angle = j * M_ast_rad + k * M_J_rad
        elif planet == "S":
            angle = j * M_ast_rad + k * M_S_rad
        elif planet == "U":
            angle = j * M_ast_rad + k * M_U_rad
        else:
            raise ValueError(f"Unknown planet: {planet}")

        X[:, 2 * idx] = np.cos(angle)
        X[:, 2 * idx + 1] = np.sin(angle)

    return X


# =============================================================================
# DETRENDING
# =============================================================================


def _detrend(jd_array: np.ndarray, residuals: np.ndarray, degree: int) -> np.ndarray:
    """Remove polynomial trend from residuals to isolate periodic component.

    Uses normalized time variable to avoid numerical issues with high-degree
    polynomials at large JD values.

    Returns detrended residuals (same shape as input).
    """
    # Normalize time to [-1, 1] for numerical stability
    t_center = (jd_array[0] + jd_array[-1]) / 2.0
    t_scale = (jd_array[-1] - jd_array[0]) / 2.0
    t_norm = (jd_array - t_center) / t_scale

    # Fit polynomial
    poly_coeffs = np.polyfit(t_norm, residuals, degree)
    trend = np.polyval(poly_coeffs, t_norm)

    return residuals - trend


# =============================================================================
# FITTING AND PRUNING
# =============================================================================


def _fit_and_prune(
    X: np.ndarray,
    residuals: np.ndarray,
    labels: list[tuple[int, int, str]],
    amplitude_threshold: float,
    max_terms: int,
) -> list[tuple[int, int, str, float, float]]:
    """Fit Fourier coefficients to detrended residuals and prune small terms.

    Args:
        X: Design matrix (n_samples, 2*n_terms)
        residuals: Detrended target values (n_samples,)
        labels: Basis function labels
        amplitude_threshold: Minimum amplitude to keep a term
        max_terms: Maximum number of terms to keep

    Returns:
        List of (j, k, planet, a_coeff, b_coeff) tuples, sorted by amplitude.
        The model is: sum(a * cos(angle) + b * sin(angle)).
    """
    # Solve least squares
    coeffs, _, _, _ = np.linalg.lstsq(X, residuals, rcond=None)

    # Extract (a, b) pairs and compute amplitudes
    significant = []
    for idx, (j, k, planet) in enumerate(labels):
        a = coeffs[2 * idx]
        b = coeffs[2 * idx + 1]
        amplitude = math.sqrt(a**2 + b**2)
        if amplitude >= amplitude_threshold:
            significant.append((j, k, planet, float(a), float(b), amplitude))

    # Sort by amplitude (descending) and limit
    significant.sort(key=lambda x: x[5], reverse=True)
    significant = significant[:max_terms]

    # Return without the amplitude column
    return [(j, k, p, a, b) for j, k, p, a, b, _ in significant]


def _evaluate_correction(
    terms: list[tuple[int, int, str, float, float]],
    M_ast: np.ndarray,
    M_J: np.ndarray,
    M_S: np.ndarray,
    M_U: np.ndarray,
) -> np.ndarray:
    """Evaluate the Fourier correction from fitted terms."""
    result = np.zeros(len(M_ast))

    M_ast_rad = np.radians(M_ast)
    M_J_rad = np.radians(M_J)
    M_S_rad = np.radians(M_S)
    M_U_rad = np.radians(M_U)

    for j, k, planet, a, b in terms:
        if planet == "A":
            angle = j * M_ast_rad
        elif planet == "J":
            angle = j * M_ast_rad + k * M_J_rad
        elif planet == "S":
            angle = j * M_ast_rad + k * M_S_rad
        elif planet == "U":
            angle = j * M_ast_rad + k * M_U_rad
        else:
            continue

        result += a * np.cos(angle) + b * np.sin(angle)

    return result


# =============================================================================
# VALIDATION
# =============================================================================


def _validate_on_windows(
    body_id: int,
    spk_file: Path,
    lon_terms: list,
    lat_terms: list,
    dist_terms: list,
    epoch_jd: float = 2461000.5,
) -> None:
    """Validate corrections at specific time offsets (matching benchmark test).

    Tests at ±1mo, ±6mo, ±1yr, ±5yr, ±10yr from the element epoch.
    """
    offsets_days = [
        (30.0, "1 month"),
        (182.6, "6 months"),
        (365.25, "1 year"),
        (5 * 365.25, "5 years"),
        (10 * 365.25, "10 years"),
        (25 * 365.25, "25 years"),
    ]

    print("\n  Validation at specific offsets from epoch:")
    print(f"    {'Offset':>12s}  {'Before':>10s}  {'After':>10s}  {'Improvement':>12s}")

    for dt_days, label in offsets_days:
        errors_before = []
        errors_after = []

        for sign in [1, -1]:
            jd = epoch_jd + sign * dt_days
            try:
                # SPK reference
                from spktype21 import SPKType21

                kernel = SPKType21.open(str(spk_file))
                try:
                    center_id = kernel.segments[0].center
                    target_id = kernel.segments[0].target
                    pos_km, _ = kernel.compute_type21(center_id, target_id, jd)
                    AU_KM = 149597870.7
                    eps = math.radians(OBLIQUITY_J2000)
                    x_eq = pos_km[0] / AU_KM
                    y_eq = pos_km[1] / AU_KM
                    z_eq = pos_km[2] / AU_KM
                    x_ecl = x_eq
                    y_ecl = y_eq * math.cos(eps) + z_eq * math.sin(eps)
                    z_ecl = -y_eq * math.sin(eps) + z_eq * math.cos(eps)
                    r_spk = math.sqrt(x_ecl**2 + y_ecl**2 + z_ecl**2)
                    lon_spk = math.degrees(math.atan2(y_ecl, x_ecl)) % 360.0
                    lat_spk = math.degrees(math.asin(z_ecl / r_spk))
                finally:
                    kernel.close()

                # Keplerian
                elements = _get_closest_epoch_elements(body_id, jd)
                x, y, z = calc_minor_body_position(elements, jd, body_id=body_id)
                r_kep = math.sqrt(x**2 + y**2 + z**2)
                lon_kep = math.degrees(math.atan2(y, x)) % 360.0
                lat_kep = math.degrees(math.asin(z / r_kep))

                # Angular separation before correction
                err_before = _angular_sep(lon_spk, lat_spk, lon_kep, lat_kep)
                errors_before.append(err_before)

                # Apply correction
                jd_arr = np.array([jd])
                M_a, M_j, M_s, M_u = _compute_mean_anomalies(body_id, jd_arr)
                dlon = _evaluate_correction(lon_terms, M_a, M_j, M_s, M_u)[0]
                dlat = _evaluate_correction(lat_terms, M_a, M_j, M_s, M_u)[0]

                lon_corr = lon_kep + dlon / 3600.0
                lat_corr = lat_kep + dlat / 3600.0

                err_after = _angular_sep(lon_spk, lat_spk, lon_corr, lat_corr)
                errors_after.append(err_after)

            except Exception as exc:
                print(f"    {label:>12s}  Error: {exc}")

        if errors_before and errors_after:
            max_before = max(errors_before)
            max_after = max(errors_after)
            improvement = max_before / max_after if max_after > 0 else float("inf")
            print(
                f"    {label:>12s}  {_format_arcsec(max_before):>10s}  "
                f"{_format_arcsec(max_after):>10s}  {improvement:>10.1f}x"
            )


def _angular_sep(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
    """Angular separation in arcseconds."""
    lon1_r, lat1_r = math.radians(lon1), math.radians(lat1)
    lon2_r, lat2_r = math.radians(lon2), math.radians(lat2)
    dlon = lon2_r - lon1_r
    cos1, cos2 = math.cos(lat1_r), math.cos(lat2_r)
    sin1, sin2 = math.sin(lat1_r), math.sin(lat2_r)
    num = math.sqrt(
        (cos2 * math.sin(dlon)) ** 2 + (cos1 * sin2 - sin1 * cos2 * math.cos(dlon)) ** 2
    )
    den = sin1 * sin2 + cos1 * cos2 * math.cos(dlon)
    return math.degrees(math.atan2(num, den)) * 3600.0


def _format_arcsec(val: float) -> str:
    """Format arcsecond value with appropriate units."""
    if val < 1.0:
        return f'{val:.3f}"'
    elif val < 60.0:
        return f'{val:.1f}"'
    elif val < 3600.0:
        return f"{val / 60.0:.1f}'"
    else:
        return f"{val / 3600.0:.2f}d"


# =============================================================================
# MAIN PROCESSING
# =============================================================================


def process_body(body_id: int, name: str, spk_prefix: str) -> dict | None:
    """Process one body: sample, detrend, fit, prune, validate."""
    spk_file = _find_spk_file(spk_prefix)
    if spk_file is None:
        print(f"  WARNING: No SPK file found for {name} ({spk_prefix}), skipping")
        return None

    print(f"\n{'=' * 70}")
    print(f"Processing {name} (body_id={body_id}, SPK={spk_file.name})")
    print(f"{'=' * 70}")

    elements = MINOR_BODY_ELEMENTS[body_id]
    include_uranus = elements.a > 10.0
    if include_uranus:
        print(f"  a = {elements.a:.2f} AU — including Uranus terms")

    # Generate sample JDs
    jd_array = np.arange(FIT_START_JD, FIT_END_JD, SAMPLE_INTERVAL_DAYS)
    n_samples = len(jd_array)
    print(f"  Samples: {n_samples} at {SAMPLE_INTERVAL_DAYS}-day intervals")

    # Compute SPK reference
    t0 = time.time()
    print("  Computing SPK reference positions...", end="", flush=True)
    spk_lon, spk_lat, spk_dist, valid = _compute_spk_positions_batch(spk_file, jd_array)
    n_valid = np.sum(valid)
    print(f" done ({time.time() - t0:.1f}s, {n_valid}/{n_samples} valid)")

    if n_valid < 100:
        print(f"  ERROR: Too few valid SPK points ({n_valid}), skipping")
        return None

    # Compute Keplerian positions
    t0 = time.time()
    print("  Computing Keplerian positions...", end="", flush=True)
    kep_lon, kep_lat, kep_dist = _compute_keplerian_positions_batch(body_id, jd_array)
    print(f" done ({time.time() - t0:.1f}s)")

    # Apply valid mask
    jd_valid = jd_array[valid]
    spk_lon = spk_lon[valid]
    spk_lat = spk_lat[valid]
    spk_dist = spk_dist[valid]
    kep_lon = kep_lon[valid]
    kep_lat = kep_lat[valid]
    kep_dist = kep_dist[valid]

    # Compute residuals
    dlon = (spk_lon - kep_lon + 180.0) % 360.0 - 180.0
    dlat = spk_lat - kep_lat
    ddist = spk_dist - kep_dist

    dlon_arcsec = dlon * 3600.0
    dlat_arcsec = dlat * 3600.0

    print("  Raw residuals:")
    print(
        f'    lon:  RMS={np.std(dlon_arcsec):.1f}", max={np.max(np.abs(dlon_arcsec)):.1f}"'
    )
    print(
        f'    lat:  RMS={np.std(dlat_arcsec):.1f}", max={np.max(np.abs(dlat_arcsec)):.1f}"'
    )
    print(f"    dist: RMS={np.std(ddist):.6f} AU, max={np.max(np.abs(ddist)):.6f} AU")

    # Detrend to remove secular drift
    print(f"  Detrending with polynomial degree {DETREND_DEGREE}...")
    dlon_detrended = _detrend(jd_valid, dlon_arcsec, DETREND_DEGREE)
    dlat_detrended = _detrend(jd_valid, dlat_arcsec, DETREND_DEGREE)
    ddist_detrended = _detrend(jd_valid, ddist, DETREND_DEGREE)

    print("  Detrended residuals:")
    print(
        f'    lon:  RMS={np.std(dlon_detrended):.1f}", max={np.max(np.abs(dlon_detrended)):.1f}"'
    )
    print(
        f'    lat:  RMS={np.std(dlat_detrended):.1f}", max={np.max(np.abs(dlat_detrended)):.1f}"'
    )
    print(
        f"    dist: RMS={np.std(ddist_detrended):.6f} AU, max={np.max(np.abs(ddist_detrended)):.6f} AU"
    )

    # Compute mean anomalies
    t0 = time.time()
    print("  Computing mean anomalies...", end="", flush=True)
    M_ast, M_J, M_S, M_U = _compute_mean_anomalies(body_id, jd_valid)
    print(f" done ({time.time() - t0:.1f}s)")

    # Build basis
    labels = _build_basis_labels(include_uranus=include_uranus)
    print(f"  Basis: {len(labels)} terms")

    # Build design matrix
    t0 = time.time()
    print("  Building design matrix...", end="", flush=True)
    X = _build_design_matrix(labels, M_ast, M_J, M_S, M_U)
    print(f" done ({time.time() - t0:.1f}s)")

    # Fit detrended residuals
    t0 = time.time()
    print("  Fitting detrended longitude...", end="", flush=True)
    lon_terms = _fit_and_prune(
        X,
        dlon_detrended,
        labels,
        AMPLITUDE_THRESHOLD_LON_ARCSEC,
        MAX_TERMS_PER_COORD,
    )
    print(f" {len(lon_terms)} terms")

    print("  Fitting detrended latitude...", end="", flush=True)
    lat_terms = _fit_and_prune(
        X,
        dlat_detrended,
        labels,
        AMPLITUDE_THRESHOLD_LAT_ARCSEC,
        MAX_TERMS_PER_COORD,
    )
    print(f" {len(lat_terms)} terms")

    print("  Fitting detrended distance...", end="", flush=True)
    dist_terms = _fit_and_prune(
        X,
        ddist_detrended,
        labels,
        AMPLITUDE_THRESHOLD_DIST_AU,
        MAX_TERMS_PER_COORD,
    )
    print(f" {len(dist_terms)} terms")
    print(f"  Fitting complete ({time.time() - t0:.1f}s)")

    # Evaluate correction quality on detrended data
    corr_lon = _evaluate_correction(lon_terms, M_ast, M_J, M_S, M_U)
    corr_lat = _evaluate_correction(lat_terms, M_ast, M_J, M_S, M_U)
    corr_dist = _evaluate_correction(dist_terms, M_ast, M_J, M_S, M_U)

    res_lon = dlon_detrended - corr_lon
    res_lat = dlat_detrended - corr_lat
    res_dist = ddist_detrended - corr_dist

    print("  After Fourier correction (detrended):")
    print(f'    lon:  RMS={np.std(res_lon):.1f}", max={np.max(np.abs(res_lon)):.1f}"')
    print(f'    lat:  RMS={np.std(res_lat):.1f}", max={np.max(np.abs(res_lat)):.1f}"')
    print(
        f"    dist: RMS={np.std(res_dist):.6f} AU, max={np.max(np.abs(res_dist)):.6f} AU"
    )

    rr_lon = (
        (1.0 - np.std(res_lon) / np.std(dlon_detrended)) * 100
        if np.std(dlon_detrended) > 0
        else 0
    )
    rr_lat = (
        (1.0 - np.std(res_lat) / np.std(dlat_detrended)) * 100
        if np.std(dlat_detrended) > 0
        else 0
    )
    print(f"  Periodic RMS reduction: lon={rr_lon:.1f}%, lat={rr_lat:.1f}%")

    # Print top terms
    print("\n  Top longitude terms:")
    for j, k, p, a, b in lon_terms[:10]:
        amp = math.sqrt(a**2 + b**2)
        phase = math.degrees(math.atan2(b, a))
        if p == "A":
            desc = f"{j}*M_ast"
        else:
            desc = f"{j}*M_ast + {k}*M_{p}"
        print(f'    {desc:25s}  amp={amp:8.2f}"  phase={phase:7.1f}')

    # Validate at specific time offsets
    _validate_on_windows(
        body_id,
        spk_file,
        lon_terms,
        lat_terms,
        dist_terms,
    )

    return {
        "lon_terms": lon_terms,
        "lat_terms": lat_terms,
        "dist_terms": dist_terms,
    }


def _format_terms(terms: list[tuple[int, int, str, float, float]], indent: str) -> str:
    """Format a list of terms as Python list of tuples."""
    if not terms:
        return f"{indent}[]"

    lines = [f"{indent}["]
    for j, k, p, a, b in terms:
        p_code = {"A": 0, "J": 1, "S": 2, "U": 3}[p]
        lines.append(f"{indent}    ({j}, {k}, {p_code}, {a:.6f}, {b:.6f}),")
    lines.append(f"{indent}]")
    return "\n".join(lines)


def generate_python_code(results: dict[int, dict]) -> str:
    """Generate Python code for the SHORT_PERIOD_CORRECTIONS dict."""
    lines = []
    lines.append(
        "# ============================================================================="
    )
    lines.append("# SHORT-PERIOD PERTURBATION CORRECTIONS (empirical Fourier fit)")
    lines.append(
        "# ============================================================================="
    )
    lines.append("# Generated by scripts/generate_short_period_corrections.py")
    lines.append(
        "# Each body has corrections for longitude (arcsec), latitude (arcsec),"
    )
    lines.append(
        "# and distance (AU). Each term is (j, k, planet, a_coeff, b_coeff) where:"
    )
    lines.append("#   angle = j * M_ast + k * M_planet (radians)")
    lines.append("#   correction += a_coeff * cos(angle) + b_coeff * sin(angle)")
    lines.append("#   planet: 0=asteroid_only, 1=Jupiter, 2=Saturn, 3=Uranus")
    lines.append("#")
    lines.append("# Mean anomaly reference values at J2000.0 (JD 2451545.0):")
    lines.append(
        f"# JUPITER_M0 = {JUPITER_M0_J2000} deg, SATURN_M0 = {SATURN_M0_J2000} deg,"
    )
    lines.append(f"# URANUS_M0 = {URANUS_M0_J2000} deg")
    lines.append(
        "# Mean motions: JUPITER_N, SATURN_N, URANUS_N from module-level constants"
    )
    lines.append("")
    lines.append("# Planet mean anomalies at J2000.0 (degrees)")
    lines.append(f"_JUPITER_M0_J2000 = {JUPITER_M0_J2000}")
    lines.append(f"_SATURN_M0_J2000 = {SATURN_M0_J2000}")
    lines.append(f"_URANUS_M0_J2000 = {URANUS_M0_J2000}")
    lines.append(f"_J2000_JD = {J2000_JD}")
    lines.append("")
    lines.append("")
    lines.append("# Short-period correction coefficients per body")
    lines.append("# Structure: body_id -> (lon_terms, lat_terms, dist_terms)")
    lines.append("# Each term: (j, k, planet_code, a_coeff, b_coeff)")
    lines.append("SHORT_PERIOD_CORRECTIONS: dict[int, tuple[list, list, list]] = {")

    for body_id, data in sorted(results.items()):
        body_name = MINOR_BODY_ELEMENTS[body_id].name
        n_lon = len(data["lon_terms"])
        n_lat = len(data["lat_terms"])
        n_dist = len(data["dist_terms"])
        lines.append(
            f"    # {body_name}: {n_lon} lon + {n_lat} lat + {n_dist} dist terms"
        )
        lines.append(f"    {body_id}: (")

        lines.append("        # Longitude corrections (arcseconds)")
        lines.append(_format_terms(data["lon_terms"], "        "))
        lines.append(",")
        lines.append("        # Latitude corrections (arcseconds)")
        lines.append(_format_terms(data["lat_terms"], "        "))
        lines.append(",")
        lines.append("        # Distance corrections (AU)")
        lines.append(_format_terms(data["dist_terms"], "        "))
        lines.append(",")
        lines.append("    ),")

    lines.append("}")
    lines.append("")

    return "\n".join(lines)


def main():
    """Fit and emit the experimental short-period correction tables.

    The command applies this module's declared sampling, detrending, spectral
    selection, and amplitude thresholds to public SPK-derived residuals.  All
    fit settings are printed with the result; no compatibility-target output
    is sampled or used as a fitting authority.
    """
    print("=" * 70)
    print("SHORT-PERIOD PERTURBATION CORRECTION GENERATOR (v2 — detrended)")
    print("=" * 70)
    print(f"Fit range: JD {FIT_START_JD:.1f} to {FIT_END_JD:.1f}")
    print(f"Sample interval: {SAMPLE_INTERVAL_DAYS} days")
    print(f"Detrend degree: {DETREND_DEGREE}")
    print(f"Max terms per coordinate: {MAX_TERMS_PER_COORD}")
    print(
        f'Thresholds: lon>{AMPLITUDE_THRESHOLD_LON_ARCSEC}", '
        f'lat>{AMPLITUDE_THRESHOLD_LAT_ARCSEC}", '
        f"dist>{AMPLITUDE_THRESHOLD_DIST_AU} AU"
    )
    print()

    results = {}
    total_t0 = time.time()

    for body_id, name, spk_prefix in BODIES:
        body_result = process_body(body_id, name, spk_prefix)
        if body_result is not None:
            results[body_id] = body_result

    total_time = time.time() - total_t0
    print(f"\n{'=' * 70}")
    print(f"All bodies processed in {total_time:.1f}s")
    print(f"{'=' * 70}")

    if results:
        code = generate_python_code(results)
        print("\n\n# ===== BEGIN GENERATED CODE =====\n")
        print(code)
        print("# ===== END GENERATED CODE =====\n")

        # Write to file
        out_path = Path(__file__).parent.parent / "data" / "short_period_corrections.py"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w") as f:
            f.write(
                "# Auto-generated by scripts/generate_short_period_corrections.py\n"
            )
            f.write("# Do not edit manually.\n\n")
            f.write(code)
        print(f"Code written to {out_path}")


if __name__ == "__main__":
    main()
