#!/usr/bin/env python3
"""Fast single-segment Chebyshev parameter test.

Tests Chebyshev fitting accuracy for a given body/interval/degree combination
without generating a full .leb file. Evaluates only a few segments (around the
known worst-case JD and/or random samples), so it runs in seconds instead of
the minutes required by a full generation.

Usage:
    # Test a single combo at the known worst-case JD (< 5 seconds):
    python scripts/test_chebyshev_params.py --body jupiter --interval 1 --degree 19 --jd 2458845.3

    # Compare multiple combos side by side (< 15 seconds):
    python scripts/test_chebyshev_params.py --body jupiter --combos "1/17,1/19,1/21,2/17,2/19" --jd 2458845.3

    # Scan 200 random segments for a more robust estimate (30-60 seconds):
    python scripts/test_chebyshev_params.py --body jupiter --interval 1 --degree 19 --scan 200

    # Scan random + worst-case combined:
    python scripts/test_chebyshev_params.py --body jupiter --interval 1 --degree 19 --jd 2458845.3 --scan 200

    # Find worst-case across ALL segments (slow but definitive):
    python scripts/test_chebyshev_params.py --body jupiter --interval 1 --degree 19 --scan-all

Provenance:
    Project-authored focused parameter test using the registered JPL/IAU source
    channel and public Chebyshev approximation. User-supplied intervals, degrees,
    and probe epochs affect only the experiment. Reported errors must be followed
    by full-range verification and are never emitted as runtime correction data.
"""

from __future__ import annotations

import argparse
import sys
import time

import numpy as np
from numpy.polynomial.chebyshev import chebfit, chebval

sys.path.insert(0, ".")
from libephemeris import set_calc_mode, calc

# ── Body name mapping ───────────────────────────────────────────────────────

_BODY_NAMES: dict[str, int] = {
    "sun": 0,
    "moon": 1,
    "mercury": 2,
    "venus": 3,
    "mars": 4,
    "jupiter": 5,
    "saturn": 6,
    "uranus": 7,
    "neptune": 8,
    "pluto": 9,
    "earth": 14,
}

_BODY_LABELS: dict[int, str] = {v: k.capitalize() for k, v in _BODY_NAMES.items()}

# Base tier range (default)
BASE_JD_START = 2396758.5
BASE_JD_END = 2506331.5


# ── Chebyshev helpers ───────────────────────────────────────────────────────


def _chebyshev_nodes(n: int) -> np.ndarray:
    """Chebyshev Type-I nodes on [-1, 1]."""
    return np.cos(np.pi * (np.arange(n) + 0.5) / n)


def _test_segment(
    body_id: int,
    seg_start: float,
    seg_end: float,
    degree: int,
    n_test: int = 50,
) -> tuple[float, float, float, float]:
    """Test Chebyshev fitting for a single segment.

    Returns (max_lon_err_arcsec, max_lat_err_arcsec, worst_lon_jd, worst_lat_jd).
    """
    nodes = _chebyshev_nodes(degree + 1)
    jd_nodes = 0.5 * (seg_end - seg_start) * nodes + 0.5 * (seg_start + seg_end)

    # Evaluate at Chebyshev nodes
    values = np.zeros((degree + 1, 3))
    for i, jd in enumerate(jd_nodes):
        result, _ = calc(float(jd), body_id, 256)
        values[i] = [result[0], result[1], result[2]]

    # Unwrap longitude to remove 360° jumps before fitting
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
    worst_lon_jd = seg_start
    worst_lat_jd = seg_start

    for i in range(n_test):
        frac = (i + 0.5) / n_test
        jd = seg_start + frac * (seg_end - seg_start)
        tau = (jd - mid) / half

        ref, _ = calc(float(jd), body_id, 256)

        # Longitude: re-wrap for comparison
        fitted_lon = float(chebval(tau, coeffs[0])) % 360.0
        ref_lon = ref[0] % 360.0
        lon_err = abs(fitted_lon - ref_lon)
        if lon_err > 180.0:
            lon_err = 360.0 - lon_err
        lon_err *= 3600.0  # arcsec

        # Latitude: direct comparison
        lat_err = abs(float(chebval(tau, coeffs[1])) - ref[1]) * 3600.0

        if lon_err > max_lon_err:
            max_lon_err = lon_err
            worst_lon_jd = jd
        if lat_err > max_lat_err:
            max_lat_err = lat_err
            worst_lat_jd = jd

    return max_lon_err, max_lat_err, worst_lon_jd, worst_lat_jd


# ── Segment selection ───────────────────────────────────────────────────────


def _build_segment_list(
    jd_start: float,
    jd_end: float,
    interval: float,
    worst_jd: float | None,
    scan: int,
    scan_all: bool,
) -> list[tuple[str, float, float]]:
    """Build a list of (label, seg_start, seg_end) to test."""
    n_total = int((jd_end - jd_start) / interval)
    segments: list[tuple[str, float, float]] = []
    seen: set[int] = set()

    def _add(label: str, idx: int) -> None:
        if 0 <= idx < n_total and idx not in seen:
            seen.add(idx)
            s = jd_start + idx * interval
            segments.append((label, s, s + interval))

    # Segment(s) around the known worst-case JD
    if worst_jd is not None:
        idx = int((worst_jd - jd_start) / interval)
        _add("worst-case", idx)
        _add("neighbor(-1)", idx - 1)
        _add("neighbor(+1)", idx + 1)

    if scan_all:
        # Every segment
        for i in range(n_total):
            _add(f"seg-{i}", i)
    elif scan > 0:
        # Random sample
        rng = np.random.default_rng(42)
        indices = rng.choice(n_total, min(scan, n_total), replace=False)
        for i in sorted(indices):
            _add(f"scan-{i}", int(i))

    # If nothing was requested, default to 200 evenly-spaced segments
    if not segments:
        step = max(1, n_total // 200)
        for i in range(0, n_total, step):
            _add(f"grid-{i}", i)

    return segments


# ── Run one param combo ─────────────────────────────────────────────────────


def _run_combo(
    body_id: int,
    interval: float,
    degree: int,
    segments: list[tuple[str, float, float]],
    n_test: int,
) -> tuple[float, float, str, str]:
    """Run a single interval/degree combo. Returns (max_lon, max_lat, lon_info, lat_info)."""
    max_lon = 0.0
    max_lat = 0.0
    lon_info = ""
    lat_info = ""

    for label, seg_s, seg_e in segments:
        lon_err, lat_err, lon_jd, lat_jd = _test_segment(
            body_id, seg_s, seg_e, degree, n_test
        )
        if lon_err > max_lon:
            max_lon = lon_err
            lon_info = f"{label} JD {lon_jd:.1f}"
        if lat_err > max_lat:
            max_lat = lat_err
            lat_info = f"{label} JD {lat_jd:.1f}"

    return max_lon, max_lat, lon_info, lat_info


# ── Main ────────────────────────────────────────────────────────────────────


def main() -> None:
    """Run focused Chebyshev interval/degree experiments for one body.

    User-supplied combinations and probe dates affect this diagnostic only.
    Random scans use the script's reproducible sampling policy, while an
    all-segment scan still tests finitely many points per segment; neither mode
    silently promotes measured residuals into runtime coefficients.
    """
    parser = argparse.ArgumentParser(
        description="Fast Chebyshev parameter tester (single-segment)"
    )
    parser.add_argument("--body", required=True, help="Body name or numeric ID")
    parser.add_argument("--interval", type=float, help="Interval in days")
    parser.add_argument("--degree", type=int, help="Chebyshev degree")
    parser.add_argument(
        "--combos",
        help='Comma-separated interval/degree combos, e.g. "1/17,1/19,2/19"',
    )
    parser.add_argument("--jd", type=float, help="Known worst-case JD to test around")
    parser.add_argument(
        "--scan", type=int, default=0, help="Number of random segments to scan"
    )
    parser.add_argument(
        "--scan-all", action="store_true", help="Scan ALL segments (definitive)"
    )
    parser.add_argument(
        "--n-test", type=int, default=50, help="Test points per segment (default: 50)"
    )
    parser.add_argument("--jd-start", type=float, default=BASE_JD_START)
    parser.add_argument("--jd-end", type=float, default=BASE_JD_END)
    args = parser.parse_args()

    # Resolve body
    try:
        body_id = int(args.body)
    except ValueError:
        body_id = _BODY_NAMES.get(args.body.lower())  # type: ignore[assignment]
        if body_id is None:
            print(f"Unknown body: {args.body}")
            sys.exit(1)

    body_label = _BODY_LABELS.get(body_id, f"Body {body_id}")

    # Build combo list
    combos: list[tuple[float, int]] = []
    if args.combos:
        for token in args.combos.split(","):
            token = token.strip()
            iv, dg = token.split("/")
            combos.append((float(iv), int(dg)))
    elif args.interval is not None and args.degree is not None:
        combos.append((args.interval, args.degree))
    else:
        print("Error: provide --interval/--degree or --combos")
        sys.exit(1)

    set_calc_mode("skyfield")

    print(f"{'=' * 65}")
    print(f"  {body_label} — Chebyshev parameter test")
    print(f"  Range: JD {args.jd_start:.1f} – {args.jd_end:.1f}")
    if args.jd:
        print(f"  Worst-case JD: {args.jd}")
    print(f"  Test points per segment: {args.n_test}")
    print(f"{'=' * 65}")

    # Header
    hdr_lon = 'lon"'
    hdr_lat = 'lat"'
    hdr_max = 'max"'
    print(
        f"\n  {'Params':>8s}  {hdr_lon:>8s}  {hdr_lat:>8s}  {hdr_max:>8s}  {'Status':>6s}  Time"
    )
    print(f"  {'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 8}  {'─' * 6}  {'─' * 6}")

    for interval, degree in combos:
        segments = _build_segment_list(
            args.jd_start, args.jd_end, interval, args.jd, args.scan, args.scan_all
        )

        t0 = time.monotonic()
        max_lon, max_lat, lon_info, lat_info = _run_combo(
            body_id, interval, degree, segments, args.n_test
        )
        elapsed = time.monotonic() - t0

        worst = max(max_lon, max_lat)
        status = "PASS" if worst < 0.5 else "FAIL"
        marker = "✅" if status == "PASS" else "❌"

        combo_str = f"{interval:g}d/{degree}"
        print(
            f"  {combo_str:>8s}  {max_lon:8.4f}  {max_lat:8.4f}  {worst:8.4f}  {marker} {status}  {elapsed:.1f}s"
        )

    print()


if __name__ == "__main__":
    main()
