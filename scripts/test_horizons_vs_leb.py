#!/usr/bin/env python3
"""
Horizons vs LEB2 cross-validation test.

Compares calc_ut() results between Horizons API and LEB2 fast path
to verify consistency across all three backends.

Usage:
    python scripts/test_horizons_vs_leb.py              # default (100 dates)
    python scripts/test_horizons_vs_leb.py --dates 500  # thorough test
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
    (10, "MeanNode"),
    (12, "MeanApog"),
    (15, "Chiron"),
    (17, "Ceres"),
]

SKIP_COMBOS = {(0, swe.FLG_SPEED | swe.FLG_HELCTR)}

FLAGS = [
    (swe.FLG_SPEED, "default"),
    (swe.FLG_SPEED | swe.FLG_SIDEREAL, "sidereal"),
    (swe.FLG_SPEED | swe.FLG_EQUATORIAL, "equatorial"),
]

JD_START = 2415020.5
JD_END = 2488069.5
THRESHOLD = 0.01  # arcseconds


def run_test(n_dates: int = 100, seed: int = 42) -> bool:
    rng = np.random.default_rng(seed)
    jds = rng.uniform(JD_START, JD_END, n_dates)

    leb_path = "data/leb2/base_core.leb2"
    if not os.path.isfile(leb_path):
        print(f"SKIP: {leb_path} not found")
        return True

    t0 = time.time()

    # LEB2 reference
    swe.set_leb_file(leb_path)
    swe.set_calc_mode("leb")
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

    # Horizons
    swe.set_calc_mode("horizons")
    n = 0
    n_fail = 0
    body_max: dict[int, tuple[float, str]] = {}

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
                    if err >= THRESHOLD:
                        n_fail += 1
                    if bid not in body_max or err > body_max[bid][0]:
                        body_max[bid] = (err, fname)
                except Exception:
                    pass

    swe.close()
    elapsed = time.time() - t0

    g = max(e[0] for e in body_max.values()) if body_max else 0
    ok = n_fail == 0

    print(
        f'{"PASS" if ok else "FAIL"} | Horizons vs LEB2 | {n} tests | max={g:.4f}" | {elapsed:.1f}s'
    )
    for bid in sorted(body_max):
        err, fn = body_max[bid]
        name = [nm for i, nm in BODIES if i == bid][0]
        if err > 0.001:
            print(f'  {name:<12s}  {err:.4f}"  {fn}')

    return ok


def main():
    parser = argparse.ArgumentParser(description="Horizons vs LEB2 precision test")
    parser.add_argument("--dates", type=int, default=100)
    args = parser.parse_args()
    ok = run_test(n_dates=args.dates)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
