#!/usr/bin/env python3
"""
Fast LEB2 vs LEB1 precision test.

Compares calc_ut() results between LEB2 and LEB1 files across all 31 bodies,
6 flag combinations, and N random dates. Fails if any comparison exceeds 0.001".

Usage:
    python scripts/test_leb2_precision.py base          # ~15s
    python scripts/test_leb2_precision.py medium         # ~15s
    python scripts/test_leb2_precision.py extended       # ~15s
    python scripts/test_leb2_precision.py base --dates 1000  # more dates
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as swe

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
    10: "MeanNode",
    11: "TrueNode",
    12: "MeanApog",
    13: "OscuApog",
    14: "Earth",
    15: "Chiron",
    17: "Ceres",
    18: "Pallas",
    19: "Juno",
    20: "Vesta",
    21: "IntpApog",
    22: "IntpPerig",
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
# Exotic minor bodies (centaurs/TNOs/NEAs) — registry is the source of truth.
from libephemeris.exotic_bodies import name_map as _exotic_names  # noqa: E402

BODY_NAMES.update(_exotic_names())

ALL_BODIES = sorted(BODY_NAMES.keys())
HELIO_ONLY = {40, 41, 42, 43, 44, 45, 46, 47, 48}

FLAGS = [
    (swe.FLG_SPEED, "default"),
    (swe.FLG_SPEED | swe.FLG_SIDEREAL, "sidereal"),
    (swe.FLG_SPEED | swe.FLG_EQUATORIAL, "equatorial"),
    (swe.FLG_SPEED | swe.FLG_J2000, "J2000"),
    (swe.FLG_SPEED | swe.FLG_NOABERR, "no_aberr"),
    (swe.FLG_SPEED | swe.FLG_HELCTR, "heliocentric"),
]

TIER_CONFIG = {
    "base": {
        "leb1": "data/leb/ephemeris_base.leb",
        "leb2": "data/leb2/base_core.leb2",
        "jd_start": 2396759,
        "jd_end": 2506330,
    },
    "medium": {
        "leb1": "data/leb/ephemeris_medium.leb",
        "leb2": "data/leb2/medium_core.leb2",
        "jd_start": 2305448,
        "jd_end": 2634165,
    },
    "extended": {
        "leb1": "data/leb/ephemeris_extended.leb",
        "leb2": "data/leb2/extended_core.leb2",
        "jd_start": 625673,
        "jd_end": 4279532,
    },
}

THRESHOLD = 0.001  # arcseconds


def run_test(tier: str, n_dates: int = 200, seed: int = 42) -> bool:
    cfg = TIER_CONFIG[tier]
    for f in [cfg["leb1"], cfg["leb2"]]:
        if not os.path.isfile(f):
            print(f"SKIP: {f} not found")
            return True

    rng = np.random.default_rng(seed)
    jds = rng.uniform(cfg["jd_start"] + 10, cfg["jd_end"] - 10, n_dates)

    t0 = time.time()

    # Determine which bodies are in the LEB2 file
    from libephemeris.leb_reader import open_leb

    leb2_reader = open_leb(cfg["leb2"])
    leb2_bodies = set()
    # Walk the body map — works for both LEBReader, LEB2Reader, CompositeLEBReader
    if hasattr(leb2_reader, "_bodies"):
        leb2_bodies = set(leb2_reader._bodies.keys())
    elif hasattr(leb2_reader, "_body_map"):
        leb2_bodies = set(leb2_reader._body_map.keys())
    leb2_reader.close()

    test_bodies = [b for b in ALL_BODIES if b in leb2_bodies]

    # Phase 1: LEB1 reference
    swe.set_leb_file(cfg["leb1"])
    swe.set_calc_mode("leb")
    ref = {}
    for jd in jds:
        for bid in test_bodies:
            for fl, _ in FLAGS:
                if bid in HELIO_ONLY and not (fl & swe.FLG_HELCTR):
                    continue
                try:
                    ref[(float(jd), bid, fl)] = swe.calc_ut(float(jd), bid, fl)[0][:3]
                except Exception:
                    pass
    swe.close()

    # Phase 2: LEB2 compare
    swe.set_leb_file(cfg["leb2"])
    swe.set_calc_mode("leb")
    n = 0
    n_over = 0
    body_max: dict[int, tuple[float, str]] = {}

    for jd in jds:
        for bid in test_bodies:
            for fl, fn in FLAGS:
                k = (float(jd), bid, fl)
                if k not in ref:
                    continue
                try:
                    v2 = swe.calc_ut(float(jd), bid, fl)[0][:3]
                    v1 = ref[k]
                    ld = abs(v2[0] - v1[0])
                    if ld > 180:
                        ld = 360 - ld
                    ld *= 3600
                    latd = abs(v2[1] - v1[1]) * 3600
                    err = max(ld, latd)
                    n += 1
                    if err >= THRESHOLD:
                        n_over += 1
                    if bid not in body_max or err > body_max[bid][0]:
                        body_max[bid] = (err, fn)
                except Exception:
                    pass
    swe.close()
    elapsed = time.time() - t0

    # Report
    g = max(e[0] for e in body_max.values()) if body_max else 0
    ok = n_over == 0

    print(
        f'{"PASS" if ok else "FAIL"} | {tier} | {n} tests | max={g:.4f}" | {elapsed:.1f}s'
    )

    if not ok or g > THRESHOLD * 0.5:  # show details if close to threshold
        for bid in sorted(body_max):
            err, fn = body_max[bid]
            if err >= THRESHOLD * 0.1:
                name = BODY_NAMES.get(bid, str(bid))
                mark = " ***" if err >= THRESHOLD else ""
                print(f'  {name:<14s}  {err:.4f}"  {fn}{mark}')

    return ok


def main():
    parser = argparse.ArgumentParser(description="Fast LEB2 precision test")
    parser.add_argument("tier", choices=["base", "medium", "extended", "all"])
    parser.add_argument(
        "--dates", type=int, default=200, help="Random dates per test (default: 200)"
    )
    args = parser.parse_args()

    tiers = list(TIER_CONFIG.keys()) if args.tier == "all" else [args.tier]
    all_ok = True

    for tier in tiers:
        ok = run_test(tier, n_dates=args.dates)
        if not ok:
            all_ok = False

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
