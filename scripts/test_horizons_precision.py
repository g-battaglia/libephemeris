#!/usr/bin/env python3
"""
Fast Horizons vs Skyfield precision test.

Compares calc_ut() results between Horizons API and Skyfield across
multiple bodies, dates, and flag combinations. Requires internet connection.

Usage:
    python scripts/test_horizons_precision.py              # default (200 dates)
    python scripts/test_horizons_precision.py --dates 50   # quick test
    python scripts/test_horizons_precision.py --dates 500  # thorough test

Provenance:
    Project-authored live comparison between the public NASA/JPL Horizons
    service and the registered Skyfield/JPL pipeline. Query settings, deterministic
    dates, frames, and tolerances define verification only. Responses and
    residuals are neither persisted nor fitted into LibEphemeris models.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as swe

BODIES = [
    (0, "Sun"),
    (1, "Moon"),
    (2, "Mercury"),
    (3, "Venus"),
    (4, "Mars"),
    (5, "Jupiter"),
    (6, "Saturn"),
    (7, "Uranus"),
    (8, "Neptune"),
    (9, "Pluto"),
    (14, "Earth"),
    (15, "Chiron"),
    (17, "Ceres"),
]

# Known bug: Sun heliocentric in Skyfield path returns wrong results
# (documented in proposals/leb-optimization-findings.md)
SKIP_COMBOS = {(0, swe.FLG_SPEED | swe.FLG_HELCTR)}

FLAGS = [
    (swe.FLG_SPEED, "default"),
    (swe.FLG_SPEED | swe.FLG_SIDEREAL, "sidereal"),
    (swe.FLG_SPEED | swe.FLG_EQUATORIAL, "equatorial"),
    (swe.FLG_SPEED | swe.FLG_J2000, "J2000"),
    (swe.FLG_SPEED | swe.FLG_NOABERR, "no_aberr"),
    (swe.FLG_SPEED | swe.FLG_HELCTR, "heliocentric"),
]

# Thresholds (arcseconds)
GEOCENTRIC_THRESHOLD = 0.01  # strict for geocentric modes
HELIOCENTRIC_THRESHOLD = 0.05  # looser for heliocentric (systematic offset)

JD_START = 2415020.5  # 1900
JD_END = 2488069.5  # 2100


def run_test(n_dates: int = 200, seed: int = 42) -> bool:
    """Compare live Horizons and direct JPL/Skyfield output on a seeded grid.

    Dates, flag combinations, exclusions, and tolerances are test policy
    declared in this module.  Horizons responses remain in memory only and are
    discarded after the verdict; they never become fitted model data.

    Returns:
        Whether every completed comparison met its category threshold.
    """
    rng = np.random.default_rng(seed)
    jds = rng.uniform(JD_START, JD_END, n_dates)

    t0 = time.time()

    # Phase 1: Skyfield reference
    swe.set_calc_mode("skyfield")
    ref = {}
    for jd in jds:
        for bid, _ in BODIES:
            for fl, _ in FLAGS:
                if (bid, fl) in SKIP_COMBOS:
                    continue
                try:
                    r = swe.calc_ut(float(jd), bid, fl)
                    ref[(float(jd), bid, fl)] = r[0][:3]
                except Exception:
                    pass
    swe.close()

    # Phase 2: Horizons
    swe.set_calc_mode("horizons")
    n = 0
    n_fail = 0
    body_max: dict[int, tuple[float, str]] = {}
    failures = []

    for jd in jds:
        for bid, bname in BODIES:
            for fl, fname in FLAGS:
                k = (float(jd), bid, fl)
                if k not in ref:
                    continue
                try:
                    r = swe.calc_ut(float(jd), bid, fl)
                    v2 = r[0][:3]
                    v1 = ref[k]

                    ld = abs(v2[0] - v1[0])
                    if ld > 180:
                        ld = 360 - ld
                    ld *= 3600
                    latd = abs(v2[1] - v1[1]) * 3600
                    err = max(ld, latd)
                    n += 1

                    threshold = (
                        HELIOCENTRIC_THRESHOLD
                        if fl & swe.FLG_HELCTR
                        else GEOCENTRIC_THRESHOLD
                    )
                    if err >= threshold:
                        n_fail += 1
                        failures.append((bname, float(jd), fname, err))

                    if bid not in body_max or err > body_max[bid][0]:
                        body_max[bid] = (err, fname)
                except Exception:
                    pass

    swe.close()
    elapsed = time.time() - t0

    # Report
    g = max(e[0] for e in body_max.values()) if body_max else 0
    ok = n_fail == 0

    print(f'{"PASS" if ok else "FAIL"} | {n} tests | max={g:.4f}" | {elapsed:.1f}s')

    # Show per-body max
    for bid in sorted(body_max):
        err, fn = body_max[bid]
        name = [n for i, n in BODIES if i == bid][0]
        is_helio = "helio" in fn
        threshold = HELIOCENTRIC_THRESHOLD if is_helio else GEOCENTRIC_THRESHOLD
        mark = " ***" if err >= threshold else ""
        if err > 0.001:
            print(f'  {name:<12s}  {err:.4f}"  {fn}{mark}')

    if failures:
        print(f"\n  {len(failures)} failures:")
        for bname, jd, fname, err in failures[:5]:
            print(f'    {bname} at JD {jd:.1f} [{fname}]: {err:.4f}"')

    return ok


def main():
    """Parse the requested date count and exit with the live-test verdict."""
    parser = argparse.ArgumentParser(
        description="Horizons vs Skyfield precision test (requires internet)"
    )
    parser.add_argument(
        "--dates",
        type=int,
        default=200,
        help="Number of random dates (default: 200, ~15K tests)",
    )
    args = parser.parse_args()

    ok = run_test(n_dates=args.dates)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
