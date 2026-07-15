#!/usr/bin/env python3
"""Automated Chebyshev parameter sweep for LEB optimization.

Tests a matrix of (interval_days, degree) combinations for all LEB bodies,
measuring the maximum Chebyshev fitting error in arcseconds. Used to find
optimal parameters that minimize file size while maintaining <0.001" precision.

Usage:
    # Sweep all bodies with proposed optimized params (500 random segments each):
    python scripts/sweep_leb_params.py

    # Sweep a single body:
    python scripts/sweep_leb_params.py --body sun

    # Scan ALL segments (definitive, slower):
    python scripts/sweep_leb_params.py --scan-all

    # Custom scan count:
    python scripts/sweep_leb_params.py --scan 1000

Provenance:
    Project-authored search over explicit LEB representation parameters. Truth
    samples come from the registered JPL/IAU generator channel and errors are
    evaluated off-node. A winning interval/degree pair is an engineering choice
    accepted only after full-range validation; residuals never become physical
    coefficients or compatibility fits.
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval

sys.path.insert(0, ".")
from libephemeris import set_calc_mode, calc
from libephemeris.constants import FLG_SPEED
from libephemeris.leb_format import (
    BODY_PARAMS,
    COORD_ICRS_BARY,
    COORD_ICRS_BARY_SYSTEM,
)

# Bodies stored as ICRS barycentric in LEB — sweep cannot test these accurately
# because the sweep fits ecliptic lon/lat but LEB fits ICRS x,y,z (much smoother).
# These bodies are validated by regeneration + leb_accuracy_sweep instead.
_ICRS_BODIES = {
    bid
    for bid, p in BODY_PARAMS.items()
    if p[2] in (COORD_ICRS_BARY, COORD_ICRS_BARY_SYSTEM)
}


def _calc_flags(body_id: int) -> int:
    """Return calc flags appropriate for this body type."""
    del body_id
    return FLG_SPEED


# Base tier JD range (1850-2150)
BASE_JD_START = 2396758.5
BASE_JD_END = 2506331.5

# All body names
BODY_NAMES: dict[int, str] = {
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
}

# Reverse lookup
_NAME_TO_ID: dict[str, int] = {}
for _bid, _bname in BODY_NAMES.items():
    _NAME_TO_ID[_bname.lower()] = _bid
    _NAME_TO_ID[_bname.lower().replace(" ", "")] = _bid

# Current params from BODY_PARAMS and proposed optimized params.
# Format: body_id -> list of (interval, degree) to test.
# First entry is always the current params, rest are proposed alternatives.
#
# NOTE: Only ECLIPTIC and HELIO_ECL bodies are tested here. ICRS bodies
# (planets, asteroids) are validated by regeneration + leb_accuracy_sweep,
# because this script fits ecliptic lon/lat but LEB fits ICRS x,y,z.
SWEEP_MATRIX: dict[int, list[tuple[float, int]]] = {
    # ── ECLIPTIC bodies (lon/lat/dist fitting matches LEB) ──
    # Mean Node (nearly linear, 4th-degree polynomial in T)
    10: [(128, 7), (128, 9), (256, 7)],
    # True Node
    11: [(32, 11), (32, 13), (64, 11)],
    # Mean Apogee
    12: [(32, 9), (64, 9), (64, 11)],
    # Oscu Apogee (fast oscillation ~2.6 deg/day)
    13: [(6, 15), (8, 15), (8, 17)],
    # Interp Apogee
    21: [(8, 13), (16, 13), (16, 15)],
    # Interp Perigee
    22: [(8, 13), (16, 13), (16, 15)],
}

# Precision target: 0.0005" (half of 0.001" to leave 2x safety margin)
PRECISION_TARGET = 0.0005


def _chebyshev_nodes(n: int) -> np.ndarray:
    """Return the standard ``n`` first-kind Gauss-Chebyshev nodes."""
    return np.cos(np.pi * (np.arange(n) + 0.5) / n)


def _test_segment(
    body_id: int,
    seg_start: float,
    seg_end: float,
    degree: int,
    n_test: int = 50,
) -> tuple[float, float]:
    """Test Chebyshev fitting for a single segment.

    Returns (max_lon_err_arcsec, max_lat_err_arcsec).
    """
    nodes = _chebyshev_nodes(degree + 1)
    jd_nodes = 0.5 * (seg_end - seg_start) * nodes + 0.5 * (seg_start + seg_end)

    # Evaluate at Chebyshev nodes
    values = np.zeros((degree + 1, 3))
    for i, jd in enumerate(jd_nodes):
        result, _ = calc(float(jd), body_id, _calc_flags(body_id))
        values[i] = [result[0], result[1], result[2]]

    # Unwrap longitude to remove 360 deg jumps before fitting
    values[:, 0] = np.degrees(np.unwrap(np.radians(values[:, 0])))

    # Fit each component
    coeffs = np.zeros((3, degree + 1))
    for c in range(3):
        coeffs[c] = chebfit(nodes, values[:, c], degree)

    # Test at dense intermediate points
    mid = 0.5 * (seg_start + seg_end)
    half = 0.5 * (seg_end - seg_start)

    max_lon_err = 0.0
    max_lat_err = 0.0

    for i in range(n_test):
        frac = (i + 0.5) / n_test
        jd = seg_start + frac * (seg_end - seg_start)
        tau = (jd - mid) / half

        ref, _ = calc(float(jd), body_id, _calc_flags(body_id))

        # Longitude error
        fitted_lon = float(chebval(tau, coeffs[0])) % 360.0
        ref_lon = ref[0] % 360.0
        lon_err = abs(fitted_lon - ref_lon)
        if lon_err > 180.0:
            lon_err = 360.0 - lon_err
        lon_err *= 3600.0  # arcsec

        # Latitude error
        lat_err = abs(float(chebval(tau, coeffs[1])) - ref[1]) * 3600.0

        max_lon_err = max(max_lon_err, lon_err)
        max_lat_err = max(max_lat_err, lat_err)

    return max_lon_err, max_lat_err


def _build_segments(
    jd_start: float,
    jd_end: float,
    interval: float,
    scan: int,
    scan_all: bool,
) -> list[tuple[float, float]]:
    """Build segment list for testing."""
    n_total = int((jd_end - jd_start) / interval)
    if n_total == 0:
        n_total = 1

    if scan_all:
        return [
            (jd_start + i * interval, jd_start + (i + 1) * interval)
            for i in range(n_total)
        ]

    if scan > 0 and scan < n_total:
        rng = np.random.default_rng(42)
        indices = sorted(rng.choice(n_total, min(scan, n_total), replace=False))
        return [
            (jd_start + int(i) * interval, jd_start + (int(i) + 1) * interval)
            for i in indices
        ]

    # Default: all segments if n_total <= scan, else scan random
    return [
        (jd_start + i * interval, jd_start + (i + 1) * interval) for i in range(n_total)
    ]


def _estimate_size_mb(
    interval: float, degree: int, jd_start: float, jd_end: float
) -> float:
    """Estimate file size contribution in MB for one body."""
    n_segs = int((jd_end - jd_start) / interval) + 1
    seg_bytes = 3 * (degree + 1) * 8
    return n_segs * seg_bytes / (1024 * 1024)


def sweep_body(
    body_id: int,
    combos: list[tuple[float, int]],
    jd_start: float,
    jd_end: float,
    scan: int,
    scan_all: bool,
    n_test: int,
) -> list[dict]:
    """Sweep all param combos for a single body. Returns list of result dicts."""
    results = []
    name = BODY_NAMES.get(body_id, f"Body {body_id}")

    for interval, degree in combos:
        segments = _build_segments(jd_start, jd_end, interval, scan, scan_all)
        n_segs = len(segments)

        t0 = time.monotonic()
        max_lon = 0.0
        max_lat = 0.0

        for idx, (seg_s, seg_e) in enumerate(segments):
            lon_err, lat_err = _test_segment(body_id, seg_s, seg_e, degree, n_test)
            max_lon = max(max_lon, lon_err)
            max_lat = max(max_lat, lat_err)

            # Progress every 100 segments
            if (idx + 1) % 100 == 0:
                elapsed = time.monotonic() - t0
                pct = (idx + 1) / n_segs * 100
                print(
                    f"\r    {name:16s} {interval:g}d/{degree:2d}  "
                    f"[{idx + 1}/{n_segs} segs, {pct:.0f}%] "
                    f'max_err={max(max_lon, max_lat):.6f}"  '
                    f"{elapsed:.0f}s",
                    end="",
                    flush=True,
                )

        elapsed = time.monotonic() - t0
        worst = max(max_lon, max_lat)
        size_mb = _estimate_size_mb(interval, degree, jd_start, jd_end)

        # Clear progress line
        print(f"\r{' ' * 100}\r", end="", flush=True)

        results.append(
            {
                "body_id": body_id,
                "name": name,
                "interval": interval,
                "degree": degree,
                "max_lon": max_lon,
                "max_lat": max_lat,
                "worst": worst,
                "size_mb": size_mb,
                "pass": worst < PRECISION_TARGET,
                "n_segs_tested": n_segs,
                "elapsed": elapsed,
            }
        )

    return results


def main() -> None:
    """Screen declared interval/degree candidates on reproducible segments.

    Candidate matrices, target, scan size, and random seed are LibEphemeris
    engineering choices recorded by this script.  A finite random scan is a
    filter rather than a continuous accuracy proof; ``--scan-all`` enumerates
    segments but still uses finite evaluation points within each segment.
    """
    parser = argparse.ArgumentParser(description="LEB parameter sweep")
    parser.add_argument("--body", help="Single body name or ID (default: all)")
    parser.add_argument(
        "--scan",
        type=int,
        default=500,
        help="Random segments to scan per combo (default: 500)",
    )
    parser.add_argument(
        "--scan-all", action="store_true", help="Scan ALL segments (definitive, slow)"
    )
    parser.add_argument(
        "--n-test", type=int, default=50, help="Test points per segment (default: 50)"
    )
    parser.add_argument("--jd-start", type=float, default=BASE_JD_START)
    parser.add_argument("--jd-end", type=float, default=BASE_JD_END)
    args = parser.parse_args()

    set_calc_mode("skyfield")

    # Determine which bodies to sweep
    if args.body:
        try:
            body_id = int(args.body)
        except ValueError:
            body_id = _NAME_TO_ID.get(args.body.lower())
            if body_id is None:
                print(f"Unknown body: {args.body}")
                print(f"Available: {', '.join(sorted(_NAME_TO_ID.keys()))}")
                sys.exit(1)
        bodies_to_sweep = {
            body_id: SWEEP_MATRIX.get(
                body_id, [(BODY_PARAMS[body_id][0], BODY_PARAMS[body_id][1])]
            )
        }
    else:
        bodies_to_sweep = SWEEP_MATRIX

    scan_mode = "ALL" if args.scan_all else f"{args.scan} random"
    print(f"{'=' * 85}")
    print("  LEB Parameter Sweep — Base Tier")
    print(f"  JD range: {args.jd_start:.1f} – {args.jd_end:.1f}")
    print(f"  Segments: {scan_mode} | Test points/seg: {args.n_test}")
    print(f'  Precision target: <{PRECISION_TARGET}" (2x margin below 0.001")')
    print(f"{'=' * 85}")

    all_results = []
    for body_id in sorted(bodies_to_sweep.keys()):
        combos = bodies_to_sweep[body_id]
        name = BODY_NAMES.get(body_id, f"Body {body_id}")
        print(f"\n  --- {name} (id={body_id}) ---")
        results = sweep_body(
            body_id,
            combos,
            args.jd_start,
            args.jd_end,
            args.scan,
            args.scan_all,
            args.n_test,
        )
        all_results.extend(results)

        # Print results for this body
        for r in results:
            "PASS" if r["pass"] else "FAIL"
            icon = "+" if r["pass"] else "X"
            is_current = (
                r["interval"] == BODY_PARAMS[body_id][0]
                and r["degree"] == BODY_PARAMS[body_id][1]
            )
            tag = " [current]" if is_current else ""
            print(
                f"    [{icon}] {r['interval']:>5g}d / deg {r['degree']:2d}  "
                f'lon={r["max_lon"]:.6f}"  lat={r["max_lat"]:.6f}"  '
                f'max={r["worst"]:.6f}"  size={r["size_mb"]:.2f} MB  '
                f"{r['elapsed']:.1f}s{tag}"
            )

    # Summary table
    print(f"\n{'=' * 85}")
    print("  SUMMARY — Best passing combo per body")
    print(f"{'=' * 85}")
    print(
        f"  {'Body':16s}  {'Current':>10s}  {'Proposed':>10s}  "
        f"{'Curr MB':>8s}  {'Prop MB':>8s}  {'Saving':>8s}  {'Max Err':>10s}"
    )
    print(
        f"  {'─' * 16}  {'─' * 10}  {'─' * 10}  "
        f"{'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 10}"
    )

    total_current = 0.0
    total_proposed = 0.0

    for body_id in sorted(bodies_to_sweep.keys()):
        body_results = [r for r in all_results if r["body_id"] == body_id]
        if not body_results:
            continue

        name = BODY_NAMES.get(body_id, f"Body {body_id}")
        cur_params = BODY_PARAMS[body_id]
        cur_str = f"{cur_params[0]:g}d/{cur_params[1]}"
        cur_mb = _estimate_size_mb(
            cur_params[0], cur_params[1], args.jd_start, args.jd_end
        )
        total_current += cur_mb

        # Find best passing combo (smallest size among passing)
        passing = [r for r in body_results if r["pass"]]
        if passing:
            best = min(passing, key=lambda r: r["size_mb"])
            prop_str = f"{best['interval']:g}d/{best['degree']}"
            prop_mb = best["size_mb"]
            saving = cur_mb - prop_mb
            max_err = best["worst"]
            total_proposed += prop_mb
            print(
                f"  {name:16s}  {cur_str:>10s}  {prop_str:>10s}  "
                f'{cur_mb:7.2f}M  {prop_mb:7.2f}M  {saving:+7.2f}M  {max_err:.6f}"'
            )
        else:
            # No passing combo — keep current
            total_proposed += cur_mb
            print(
                f"  {name:16s}  {cur_str:>10s}  {'(keep)':>10s}  "
                f"{cur_mb:7.2f}M  {cur_mb:7.2f}M  {'+0.00':>8s}M  {'N/A':>10s}"
            )

    saving_total = total_current - total_proposed
    pct = saving_total / total_current * 100 if total_current > 0 else 0
    print(
        f"  {'─' * 16}  {'─' * 10}  {'─' * 10}  "
        f"{'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 10}"
    )
    print(
        f"  {'TOTAL':16s}  {'':>10s}  {'':>10s}  "
        f"{total_current:7.2f}M  {total_proposed:7.2f}M  {saving_total:+7.2f}M  ({pct:.0f}%)"
    )
    print()


if __name__ == "__main__":
    main()
