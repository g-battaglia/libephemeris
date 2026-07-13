#!/usr/bin/env python3
"""
Fast LEB2 vs LEB1 precision test.

Compares calc_ut() results between a tier's core LEB2 companion and LEB1 across
every body in that core file, 6 flag combinations, and N random dates. Fails
closed on missing inputs, calculation errors, absent coverage, or an angular
comparison at or above 0.001".

Usage:
    python scripts/test_leb2_precision.py base          # ~15s
    python scripts/test_leb2_precision.py medium         # ~15s
    python scripts/test_leb2_precision.py extended       # ~15s
    python scripts/test_leb2_precision.py base --dates 1000  # more dates
"""

from __future__ import annotations

import argparse
import math
import os
import sys
import time
import warnings

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as swe
from libephemeris.leb_reader import open_leb  # noqa: E402

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
    },
    "medium": {
        "leb1": "data/leb/ephemeris_medium.leb",
        "leb2": "data/leb2/medium_core.leb2",
    },
    "extended": {
        "leb1": "data/leb/ephemeris_extended.leb",
        "leb2": "data/leb2/extended_core.leb2",
    },
}

THRESHOLD = 0.001  # arcseconds
CORE_BODY_IDS = frozenset({*range(13), 14})


def _angular_error_arcsec(value: float, reference: float) -> float:
    """Return the shortest finite angular separation in arcseconds."""
    return abs((value - reference + 180.0) % 360.0 - 180.0) * 3600.0


def run_test(tier: str, n_dates: int = 200, seed: int = 42) -> bool:
    cfg = TIER_CONFIG[tier]
    for f in [cfg["leb1"], cfg["leb2"]]:
        if not os.path.isfile(f):
            print(f"FAIL: required file not found: {f}")
            return False

    if n_dates < 1:
        print("FAIL: --dates must be at least 1")
        return False

    t0 = time.time()

    # Read exactly the named files. Runtime companion discovery is not part of
    # this focused precision gate, and the sampled interval comes from their
    # actual headers rather than duplicated tier constants.
    try:
        leb2_reader = open_leb(cfg["leb2"])
        try:
            leb2_bodies = set(leb2_reader._bodies)
            leb2_range = leb2_reader.jd_range
        finally:
            leb2_reader.close()
        leb1_reader = open_leb(cfg["leb1"])
        try:
            leb1_bodies = set(leb1_reader._bodies)
            leb1_range = leb1_reader.jd_range
        finally:
            leb1_reader.close()
    except Exception as exc:
        print(f"FAIL: cannot inspect required ephemeris headers: {exc}")
        return False

    if leb2_bodies != CORE_BODY_IDS:
        missing = sorted(CORE_BODY_IDS - leb2_bodies)
        unexpected = sorted(leb2_bodies - CORE_BODY_IDS)
        print(
            "FAIL: invalid core companion body inventory: "
            f"missing={missing}, unexpected={unexpected}"
        )
        return False
    missing_reference = sorted(CORE_BODY_IDS - leb1_bodies)
    if missing_reference:
        print(f"FAIL: LEB1 reference is missing core bodies: {missing_reference}")
        return False

    jd_start = max(leb1_range[0], leb2_range[0]) + 10.0
    jd_end = min(leb1_range[1], leb2_range[1]) - 10.0
    if jd_start >= jd_end:
        print(
            "FAIL: LEB1/LEB2 coverage has no usable overlap: "
            f"LEB1={leb1_range}, LEB2={leb2_range}"
        )
        return False
    rng = np.random.default_rng(seed)
    jds = rng.uniform(jd_start, jd_end, n_dates)
    test_bodies = sorted(CORE_BODY_IDS)

    # Phase 1: LEB1 reference
    swe.set_leb_file(cfg["leb1"])
    swe.set_calc_mode("leb")
    ref: dict[tuple[float, int, int], tuple[float, ...]] = {}
    failures: list[str] = []
    reference_coverage = {bid: 0 for bid in test_bodies}
    for jd in jds:
        for bid in test_bodies:
            for fl, _ in FLAGS:
                if bid in HELIO_ONLY and not (fl & swe.FLG_HELCTR):
                    continue
                try:
                    value = tuple(swe.calc_ut(float(jd), bid, fl)[0][:3])
                    if len(value) != 3:
                        raise ValueError(
                            f"LEB1 result has {len(value)} position components; expected 3"
                        )
                    if not all(math.isfinite(component) for component in value):
                        raise ValueError("non-finite LEB1 result")
                    ref[(float(jd), bid, fl)] = value
                    reference_coverage[bid] += 1
                except Exception as exc:
                    failures.append(
                        f"LEB1 body={bid} jd={float(jd):.6f} flags={fl}: {exc}"
                    )
    swe.close()

    # Phase 2: LEB2 compare
    swe.set_leb_file(cfg["leb2"])
    swe.set_calc_mode("leb")
    n = 0
    n_over = 0
    body_max: dict[int, tuple[float, str]] = {}
    comparison_coverage = {bid: 0 for bid in test_bodies}

    for jd in jds:
        for bid in test_bodies:
            for fl, fn in FLAGS:
                k = (float(jd), bid, fl)
                if k not in ref:
                    continue
                try:
                    v2 = tuple(swe.calc_ut(float(jd), bid, fl)[0][:3])
                    v1 = ref[k]
                    if len(v2) != 3 or len(v1) != 3:
                        raise ValueError(
                            "LEB1/LEB2 result does not contain exactly three "
                            "position components"
                        )
                    if not all(math.isfinite(component) for component in (*v1, *v2)):
                        raise ValueError("non-finite LEB1/LEB2 result")
                    ld = _angular_error_arcsec(v2[0], v1[0])
                    latd = abs(v2[1] - v1[1]) * 3600
                    err = max(ld, latd)
                    if not math.isfinite(err):
                        raise ValueError("non-finite angular error")
                    n += 1
                    comparison_coverage[bid] += 1
                    if err >= THRESHOLD:
                        n_over += 1
                    if bid not in body_max or err > body_max[bid][0]:
                        body_max[bid] = (err, fn)
                except Exception as exc:
                    failures.append(
                        f"LEB2 body={bid} jd={float(jd):.6f} flags={fl}: {exc}"
                    )
    swe.close()
    elapsed = time.time() - t0

    # Report
    g = max(e[0] for e in body_max.values()) if body_max else 0
    uncovered = [
        bid
        for bid in test_bodies
        if reference_coverage[bid] == 0 or comparison_coverage[bid] == 0
    ]
    ok = n > 0 and n_over == 0 and not failures and not uncovered

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

    if uncovered:
        print(f"  Missing per-body coverage: {uncovered}")
    if failures:
        print(f"  Calculation errors: {len(failures)}")
        for failure in failures[:10]:
            print(f"    {failure}")
        if len(failures) > 10:
            print(f"    ... {len(failures) - 10} more")
    if n == 0:
        print("  No comparisons completed")

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

    # This gate compares LEB2 compression against the same LEB1 model. Its
    # extended-tier samples intentionally cross the analytical model's advisory
    # range, so those model-validity warnings are unrelated to compression loss.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", swe.MeeusPolynomialWarning)
        for tier in tiers:
            ok = run_test(tier, n_dates=args.dates)
            if not ok:
                all_ok = False

    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
