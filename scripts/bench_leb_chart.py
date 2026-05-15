#!/usr/bin/env python3
"""Benchmark: LEB full-chart calculation speed.

Simulates realistic astrology usage: compute all bodies at the same JD,
repeated for many different dates.  This is exactly the access pattern
that the eval_cache and time-cache optimizations target.

Usage:
    # On each branch (main vs opt/leb-time-cache):
    python scripts/bench_leb_chart.py

    # Quick smoke test (fewer iterations):
    python scripts/bench_leb_chart.py --charts 50

    # Custom body set:
    python scripts/bench_leb_chart.py --charts 200 --bodies core
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time

# Ensure project root on path
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from libephemeris.constants import (
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    CHIRON,
    EARTH,
    FLG_SPEED,
    FLG_SWIEPH,
)
from libephemeris.time_utils import julday

# Body sets -------------------------------------------------------------------

BODIES_CORE = [
    SUN,
    MOON,
    MERCURY,
    VENUS,
    MARS,
    JUPITER,
    SATURN,
    URANUS,
    NEPTUNE,
    PLUTO,
    MEAN_NODE,
    TRUE_NODE,
    MEAN_APOG,
    CHIRON,
]

BODIES_MINIMAL = [SUN, MOON, MARS, JUPITER, SATURN]

BODY_SETS = {
    "core": BODIES_CORE,
    "minimal": BODIES_MINIMAL,
}


def generate_test_leb(bodies: list[int], years: int = 5) -> str:
    """Generate a temporary .leb file and return its path."""
    from scripts.generate_leb import assemble_leb

    jd_start = julday(2023, 1, 1, 0.0)
    jd_end = julday(2023 + years, 1, 1, 0.0)

    # Include EARTH — needed internally as observer
    all_bodies = sorted(set(bodies) | {EARTH})

    fd, path = tempfile.mkstemp(suffix=".leb")
    os.close(fd)

    print(f"Generating LEB file ({len(all_bodies)} bodies, {years} years)...")
    t0 = time.perf_counter()
    assemble_leb(
        output=path,
        jd_start=jd_start,
        jd_end=jd_end,
        bodies=all_bodies,
        workers=1,
        verbose=False,
    )
    print(f"  Generated in {time.perf_counter() - t0:.1f}s -> {path}")
    return path


def run_benchmark(
    leb_path: str,
    bodies: list[int],
    n_charts: int,
    iflag: int,
) -> dict:
    """Run the benchmark and return timing stats."""
    from libephemeris.leb_reader import LEBReader
    from libephemeris.fast_calc import fast_calc_ut

    reader = LEBReader(leb_path)

    # Generate evenly-spaced JDs within the file's range
    cov_start, cov_end = reader.jd_range
    jd_start = cov_start + 1.0  # small margin
    jd_end = cov_end - 1.0
    step = (jd_end - jd_start) / n_charts
    jds = [jd_start + i * step for i in range(n_charts)]

    n_calls = n_charts * len(bodies)

    # Warm-up (1 chart) to trigger any lazy loading
    for body_id in bodies:
        try:
            fast_calc_ut(reader, jds[0], body_id, iflag)
        except (KeyError, ValueError):
            pass

    # Clear any caches so the benchmark starts clean
    try:
        from libephemeris.cache import clear_caches

        clear_caches()
    except ImportError:
        pass
    reader._eval_cache.clear() if hasattr(reader, "_eval_cache") else None

    # --- Timed section -------------------------------------------------------
    errors = 0
    t_start = time.perf_counter()

    for jd in jds:
        for body_id in bodies:
            try:
                fast_calc_ut(reader, jd, body_id, iflag)
            except (KeyError, ValueError):
                errors += 1

    t_elapsed = time.perf_counter() - t_start
    # -------------------------------------------------------------------------

    reader.close()

    return {
        "charts": n_charts,
        "bodies_per_chart": len(bodies),
        "total_calls": n_calls,
        "errors": errors,
        "elapsed_s": t_elapsed,
        "ms_per_chart": (t_elapsed / n_charts) * 1000,
        "us_per_call": (t_elapsed / n_calls) * 1_000_000,
    }


def main():
    parser = argparse.ArgumentParser(description="Benchmark LEB chart calculations")
    parser.add_argument(
        "--charts",
        type=int,
        default=200,
        help="Number of charts (different JDs) to compute (default: 200)",
    )
    parser.add_argument(
        "--bodies",
        choices=list(BODY_SETS.keys()),
        default="core",
        help="Body set to use (default: core = 14 bodies)",
    )
    parser.add_argument(
        "--leb",
        type=str,
        default=None,
        help="Path to existing .leb file (skip generation)",
    )
    parser.add_argument(
        "--rounds",
        type=int,
        default=3,
        help="Number of benchmark rounds (default: 3, reports best)",
    )
    args = parser.parse_args()

    bodies = BODY_SETS[args.bodies]
    iflag = FLG_SPEED | FLG_SWIEPH

    # Get or generate LEB file
    if args.leb:
        leb_path = args.leb
        cleanup = False
    else:
        leb_path = generate_test_leb(bodies)
        cleanup = True

    print(
        f"\nBenchmark: {args.charts} charts x {len(bodies)} bodies "
        f"= {args.charts * len(bodies)} calls"
    )
    print(f"LEB file: {leb_path}")
    print(f"Rounds: {args.rounds} (reporting best)\n")

    results = []
    for i in range(args.rounds):
        r = run_benchmark(leb_path, bodies, args.charts, iflag)
        results.append(r)
        print(
            f"  Round {i + 1}: {r['elapsed_s']:.3f}s "
            f"({r['ms_per_chart']:.2f} ms/chart, "
            f"{r['us_per_call']:.1f} us/call)"
        )

    best = min(results, key=lambda r: r["elapsed_s"])

    print(f"\n{'=' * 60}")
    print(f"BEST of {args.rounds} rounds:")
    print(f"  Total time:    {best['elapsed_s']:.3f}s")
    print(f"  Per chart:     {best['ms_per_chart']:.2f} ms")
    print(f"  Per call:      {best['us_per_call']:.1f} us")
    print(f"  Errors/skips:  {best['errors']}")
    print(f"{'=' * 60}")

    if cleanup:
        os.unlink(leb_path)
        print(f"\nCleaned up temp file: {leb_path}")


if __name__ == "__main__":
    main()
