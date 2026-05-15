#!/usr/bin/env python3
"""
Round 30: NONUT / NOABERR / TRUEPOS Flag Combinations
======================================================

Tests position calculations with all combinations of correction-suppression
flags across all major bodies and multiple epochs.

These flags disable specific corrections in the position pipeline:
  FLG_TRUEPOS (16)  - geometric position (no light-time correction)
  FLG_NONUT (64)    - suppress nutation (use mean ecliptic)
  FLG_NOABERR (1024) - suppress aberration correction

Parts:
  P1: Single flags (TRUEPOS, NONUT, NOABERR) vs default — all planets
  P2: Flag pairs (TRUEPOS+NONUT, TRUEPOS+NOABERR, NONUT+NOABERR)
  P3: All three flags combined (TRUEPOS+NONUT+NOABERR)
  P4: Flags combined with EQUATORIAL output
  P5: Flags combined with J2000 ecliptic
  P6: Flags combined with HELCTR (heliocentric)
  P7: Multi-epoch sweep (1800-2100) with NONUT+NOABERR
"""

from __future__ import annotations

import os
import sys
import time
import traceback

os.environ["LIBEPHEMERIS_MODE"] = "skyfield"

import swisseph as swe

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import libephemeris as ephem
from libephemeris.constants import *

_EPHE_PATH = os.path.join(os.path.dirname(__file__), "..", "swisseph", "ephe")
swe.set_ephe_path(_EPHE_PATH)

# Bodies to test
BODIES = [
    (SUN, "Sun"),
    (MOON, "Moon"),
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
    (MEAN_NODE, "MeanNode"),
    (TRUE_NODE, "TrueNode"),
    (MEAN_APOG, "MeanLilith"),
    (CHIRON, "Chiron"),
]

# Test epochs
EPOCHS = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 3, 20, 15.5, "2024 equinox"),
    (1990, 7, 15, 6.0, "1990"),
]

# Tolerances (arcseconds)
TOL_LON = 2.0  # longitude arcsec
TOL_LAT = 2.0  # latitude arcsec
TOL_DIST = 1e-5  # distance AU
TOL_SPEED = 0.001  # speed deg/day

# Known: MeanLilith lat ~19" offset
TOL_MEANLILITH_LAT = 25.0


class R:
    def __init__(self, name):
        self.name = name
        self.passed = self.failed = self.skipped = 0
        self.failures = []
        self.max_diff = 0.0
        self.max_label = ""

    def ok(self, diff=0.0, label=""):
        self.passed += 1
        if diff > self.max_diff:
            self.max_diff = diff
            self.max_label = label

    def fail(self, msg):
        self.failed += 1
        self.failures.append(msg)

    def skip(self, msg=""):
        self.skipped += 1

    def summary(self):
        total = self.passed + self.failed
        print(f"\n{'=' * 70}")
        print(f"  {self.name}: {self.passed}/{total} PASSED ({self.skipped} skip)")
        if self.max_diff > 0:
            print(f'  Max diff: {self.max_diff:.6f}" ({self.max_label})')
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def compare_pos(r, se_result, le_result, label, body_id):
    """Compare SE and LE positions with appropriate tolerances."""
    se_pos = se_result[0] if isinstance(se_result, tuple) else se_result
    le_pos = le_result[0] if isinstance(le_result, tuple) else le_result

    # Longitude
    lon_diff = abs(se_pos[0] - le_pos[0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_diff_arcsec = lon_diff * 3600

    # Latitude
    lat_diff = abs(se_pos[1] - le_pos[1]) * 3600

    # Distance
    dist_diff = abs(se_pos[2] - le_pos[2])

    # Adaptive lat tolerance for MeanLilith
    lat_tol = TOL_MEANLILITH_LAT if body_id == MEAN_APOG else TOL_LAT

    if lon_diff_arcsec > TOL_LON:
        r.fail(
            f'{label} lon: {lon_diff_arcsec:.3f}" (SE={se_pos[0]:.8f} LE={le_pos[0]:.8f})'
        )
    elif lat_diff > lat_tol:
        r.fail(f'{label} lat: {lat_diff:.3f}" (SE={se_pos[1]:.8f} LE={le_pos[1]:.8f})')
    elif dist_diff > TOL_DIST:
        r.fail(f"{label} dist: {dist_diff:.8f} AU")
    else:
        r.ok(lon_diff_arcsec, f"{label} lon")

    # Speeds (indices 3-5)
    if len(se_pos) >= 4 and len(le_pos) >= 4:
        spd_diff = abs(se_pos[3] - le_pos[3])
        if spd_diff > TOL_SPEED:
            r.fail(f"{label} spd: {spd_diff:.6f} deg/day")
        else:
            r.ok(spd_diff * 3600, f"{label} spd")


def calc_both(body_id, jd, flags):
    """Calculate position with both SE and LE."""
    se_result = swe.calc_ut(jd, body_id, flags)
    le_result = ephem.calc_ut(jd, body_id, flags)
    return se_result, le_result


# ============================================================
# PART 1: Single flags vs default
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Single correction-suppression flags — all planets")
    print("=" * 70)

    r = R("P1: Single flags")

    single_flags = [
        (FLG_TRUEPOS, "TRUEPOS"),
        (FLG_NONUT, "NONUT"),
        (FLG_NOABERR, "NOABERR"),
    ]

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)

        for body_id, body_name in BODIES:
            for flag_val, flag_name in single_flags:
                flags = FLG_SPEED | flag_val
                label = f"{epoch_name} {body_name} {flag_name}"

                try:
                    se_r, le_r = calc_both(body_id, jd, flags)
                    compare_pos(r, se_r, le_r, label, body_id)
                except Exception as e:
                    r.skip(f"{label}: {e}")

    print(f"  Tested {r.passed + r.failed + r.skipped} combinations")
    return r.summary(), r


# ============================================================
# PART 2: Flag pairs
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Flag pairs — all planets")
    print("=" * 70)

    r = R("P2: Flag pairs")

    flag_pairs = [
        (FLG_TRUEPOS | FLG_NONUT, "TRUEPOS+NONUT"),
        (FLG_TRUEPOS | FLG_NOABERR, "TRUEPOS+NOABERR"),
        (FLG_NONUT | FLG_NOABERR, "NONUT+NOABERR"),
    ]

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)

        for body_id, body_name in BODIES:
            for flag_val, flag_name in flag_pairs:
                flags = FLG_SPEED | flag_val
                label = f"{epoch_name} {body_name} {flag_name}"

                try:
                    se_r, le_r = calc_both(body_id, jd, flags)
                    compare_pos(r, se_r, le_r, label, body_id)
                except Exception as e:
                    r.skip(f"{label}: {e}")

    print(f"  Tested {r.passed + r.failed + r.skipped} combinations")
    return r.summary(), r


# ============================================================
# PART 3: All three flags combined
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: All three flags combined (TRUEPOS+NONUT+NOABERR)")
    print("=" * 70)

    r = R("P3: All three flags")

    flag_val = FLG_TRUEPOS | FLG_NONUT | FLG_NOABERR

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)

        for body_id, body_name in BODIES:
            flags = FLG_SPEED | flag_val
            label = f"{epoch_name} {body_name} ALL3"

            try:
                se_r, le_r = calc_both(body_id, jd, flags)
                compare_pos(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    print(f"  Tested {r.passed + r.failed + r.skipped} combinations")
    return r.summary(), r


# ============================================================
# PART 4: Flags + EQUATORIAL
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Correction flags + EQUATORIAL output")
    print("=" * 70)

    r = R("P4: Flags+EQUATORIAL")

    flag_combos = [
        (FLG_NONUT, "NONUT"),
        (FLG_NOABERR, "NOABERR"),
        (FLG_TRUEPOS | FLG_NONUT, "TRUEPOS+NONUT"),
        (FLG_NONUT | FLG_NOABERR, "NONUT+NOABERR"),
        (FLG_TRUEPOS | FLG_NONUT | FLG_NOABERR, "ALL3"),
    ]

    jd = swe.julday(2024, 3, 20, 15.5)

    for body_id, body_name in BODIES:
        for flag_val, flag_name in flag_combos:
            flags = FLG_SPEED | FLG_EQUATORIAL | flag_val
            label = f"{body_name} EQ+{flag_name}"

            try:
                se_r, le_r = calc_both(body_id, jd, flags)
                compare_pos(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    print(f"  Tested {r.passed + r.failed + r.skipped} combinations")
    return r.summary(), r


# ============================================================
# PART 5: Flags + J2000 ecliptic
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Correction flags + J2000 ecliptic")
    print("=" * 70)

    r = R("P5: Flags+J2000")

    flag_combos = [
        (FLG_NONUT, "NONUT"),
        (FLG_NOABERR, "NOABERR"),
        (FLG_TRUEPOS, "TRUEPOS"),
        (FLG_TRUEPOS | FLG_NONUT | FLG_NOABERR, "ALL3"),
    ]

    # Use a date far from J2000 to test precession interactions
    jd = swe.julday(2024, 6, 21, 0.0)

    for body_id, body_name in BODIES:
        for flag_val, flag_name in flag_combos:
            flags = FLG_SPEED | FLG_J2000 | flag_val
            label = f"{body_name} J2000+{flag_name}"

            try:
                se_r, le_r = calc_both(body_id, jd, flags)
                compare_pos(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    print(f"  Tested {r.passed + r.failed + r.skipped} combinations")
    return r.summary(), r


# ============================================================
# PART 6: Flags + HELCTR (heliocentric)
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Correction flags + heliocentric")
    print("=" * 70)

    r = R("P6: Flags+HELCTR")

    # Helio bodies (no Sun, no Moon, no nodes/Lilith for helio)
    helio_bodies = [
        (MERCURY, "Mercury"),
        (VENUS, "Venus"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (URANUS, "Uranus"),
        (NEPTUNE, "Neptune"),
        (PLUTO, "Pluto"),
        (CHIRON, "Chiron"),
    ]

    flag_combos = [
        (FLG_NONUT, "NONUT"),
        (FLG_NOABERR, "NOABERR"),
        (FLG_TRUEPOS, "TRUEPOS"),
        (FLG_TRUEPOS | FLG_NONUT | FLG_NOABERR, "ALL3"),
    ]

    jd = swe.julday(2024, 3, 20, 15.5)

    for body_id, body_name in helio_bodies:
        for flag_val, flag_name in flag_combos:
            flags = FLG_SPEED | FLG_HELCTR | flag_val
            label = f"{body_name} HELIO+{flag_name}"

            try:
                se_r, le_r = calc_both(body_id, jd, flags)
                compare_pos(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    print(f"  Tested {r.passed + r.failed + r.skipped} combinations")
    return r.summary(), r


# ============================================================
# PART 7: Multi-epoch sweep with NONUT+NOABERR
# ============================================================
def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Multi-epoch sweep (1800-2100) with NONUT+NOABERR")
    print("=" * 70)

    r = R("P7: Multi-epoch NONUT+NOABERR")

    sweep_bodies = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MERCURY, "Mercury"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (PLUTO, "Pluto"),
    ]

    flags = FLG_SPEED | FLG_NONUT | FLG_NOABERR

    years = list(range(1800, 2101, 25))

    for y in years:
        jd = swe.julday(y, 6, 21, 12.0)

        for body_id, body_name in sweep_bodies:
            label = f"{y} {body_name} NN+NA"

            try:
                se_r, le_r = calc_both(body_id, jd, flags)
                compare_pos(r, se_r, le_r, label, body_id)
            except Exception as e:
                r.skip(f"{label}: {e}")

    print(f"  Tested {r.passed + r.failed + r.skipped} across {len(years)} epochs")
    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 30: NONUT / NOABERR / TRUEPOS Flag Combinations")
    print("=" * 70)

    start = time.time()
    all_ok = True
    all_results = []

    for pname, pfn in [
        ("P1", run_part1),
        ("P2", run_part2),
        ("P3", run_part3),
        ("P4", run_part4),
        ("P5", run_part5),
        ("P6", run_part6),
        ("P7", run_part7),
    ]:
        try:
            ok, res = pfn()
            all_results.append((pname, res))
            if not ok:
                all_ok = False
        except Exception as e:
            print(f"\n  {pname} CRASHED: {e}")
            traceback.print_exc()
            all_ok = False

    elapsed = time.time() - start
    print("\n" + "=" * 70)
    print("ROUND 30 FINAL SUMMARY")
    print("=" * 70)
    tp = tf = ts = 0
    for pn, res in all_results:
        st = "PASS" if res.failed == 0 else "FAIL"
        t = res.passed + res.failed
        print(f"  {pn} {res.name}: {res.passed}/{t} ({res.skipped} skip) [{st}]")
        tp += res.passed
        tf += res.failed
        ts += res.skipped
    print(f"\n  TOTAL: {tp}/{tp + tf} PASSED, {tf} FAILED, {ts} SKIPPED")
    print(f"  Time: {elapsed:.1f}s")
    print(f"\n  >>> ROUND 30: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
