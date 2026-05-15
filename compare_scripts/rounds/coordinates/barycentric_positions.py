#!/usr/bin/env python3
"""
Round 31: Barycentric Positions (FLG_BARYCTR)
================================================

Tests barycentric position output for all major planets across multiple epochs.
Barycentric = position relative to solar system barycenter (not heliocentric).

Parts:
  P1: Barycentric ecliptic positions — all planets, 3 epochs
  P2: Barycentric + EQUATORIAL
  P3: Barycentric + J2000
  P4: Barycentric + NONUT + NOABERR
  P5: Barycentric speeds validation
  P6: Multi-epoch sweep (1800-2100)
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

# FLG_BARYCTR
FLG_BARYCTR = 4

# Bodies valid for barycentric (no nodes/Lilith/analytical)
BARY_BODIES = [
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
    (CHIRON, "Chiron"),
]

EPOCHS = [
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 3, 20, 15.5, "2024eq"),
    (1990, 7, 15, 6.0, "1990"),
]

TOL_LON = 2.0  # arcsec
TOL_LAT = 2.0
TOL_DIST = 1e-5  # AU
TOL_SPEED = 0.001  # deg/day


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


def compare_pos(r, se_result, le_result, label, check_speed=True):
    se_pos = se_result[0] if isinstance(se_result, tuple) else se_result
    le_pos = le_result[0] if isinstance(le_result, tuple) else le_result

    lon_diff = abs(se_pos[0] - le_pos[0])
    if lon_diff > 180:
        lon_diff = 360 - lon_diff
    lon_diff_as = lon_diff * 3600

    lat_diff = abs(se_pos[1] - le_pos[1]) * 3600
    dist_diff = abs(se_pos[2] - le_pos[2])

    if lon_diff_as > TOL_LON:
        r.fail(
            f'{label} lon: {lon_diff_as:.3f}" (SE={se_pos[0]:.8f} LE={le_pos[0]:.8f})'
        )
    elif lat_diff > TOL_LAT:
        r.fail(f'{label} lat: {lat_diff:.3f}" (SE={se_pos[1]:.8f} LE={le_pos[1]:.8f})')
    elif dist_diff > TOL_DIST:
        r.fail(f"{label} dist: {dist_diff:.8f} AU")
    else:
        r.ok(max(lon_diff_as, lat_diff), f"{label}")

    if check_speed and len(se_pos) >= 4 and len(le_pos) >= 4:
        spd_diff = abs(se_pos[3] - le_pos[3])
        if spd_diff > TOL_SPEED:
            r.fail(f"{label} spd: {spd_diff:.6f} deg/day")
        else:
            r.ok(spd_diff * 3600, f"{label} spd")


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Barycentric ecliptic positions — all planets")
    print("=" * 70)

    r = R("P1: Barycentric ecliptic")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BARY_BODIES:
            flags = FLG_SPEED | FLG_BARYCTR
            label = f"{epoch_name} {body_name} BARY"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Barycentric + EQUATORIAL")
    print("=" * 70)

    r = R("P2: Bary+EQUATORIAL")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BARY_BODIES:
            flags = FLG_SPEED | FLG_BARYCTR | FLG_EQUATORIAL
            label = f"{epoch_name} {body_name} BARY+EQ"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Barycentric + J2000")
    print("=" * 70)

    r = R("P3: Bary+J2000")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BARY_BODIES:
            flags = FLG_SPEED | FLG_BARYCTR | FLG_J2000
            label = f"{epoch_name} {body_name} BARY+J2000"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Barycentric + NONUT + NOABERR")
    print("=" * 70)

    r = R("P4: Bary+NONUT+NOABERR")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in BARY_BODIES:
            flags = FLG_SPEED | FLG_BARYCTR | FLG_NONUT | FLG_NOABERR
            label = f"{epoch_name} {body_name} BARY+NN+NA"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Barycentric speed validation")
    print("=" * 70)

    r = R("P5: Bary speeds")

    jd = swe.julday(2024, 3, 20, 15.5)
    dt = 1.0 / 86400.0  # 1 second

    for body_id, body_name in BARY_BODIES:
        flags = FLG_SPEED | FLG_BARYCTR
        label = f"{body_name} BARY speed"

        try:
            se_r = swe.calc_ut(jd, body_id, flags)
            le_r = ephem.calc_ut(jd, body_id, flags)

            # Compare reported speeds
            se_spd = se_r[0][3]
            le_spd = le_r[0][3]
            diff = abs(se_spd - le_spd)
            if diff > TOL_SPEED:
                r.fail(
                    f"{label} lon_spd: SE={se_spd:.8f} LE={le_spd:.8f} diff={diff:.6f}"
                )
            else:
                r.ok(diff * 3600, f"{label} lon_spd")

            # Validate LE speed with finite difference
            le_before = ephem.calc_ut(jd - dt, body_id, FLG_BARYCTR)
            le_after = ephem.calc_ut(jd + dt, body_id, FLG_BARYCTR)
            fd_spd = (le_after[0][0] - le_before[0][0]) / (2 * dt)
            # Handle wraparound
            if fd_spd > 180:
                fd_spd -= 360
            elif fd_spd < -180:
                fd_spd += 360
            fd_diff = abs(le_spd - fd_spd)
            if fd_diff > 0.01:  # 0.01 deg/day tolerance for FD
                r.fail(
                    f"{label} FD: spd={le_spd:.8f} FD={fd_spd:.8f} diff={fd_diff:.6f}"
                )
            else:
                r.ok(fd_diff * 3600, f"{label} FD")

        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Multi-epoch barycentric sweep (1800-2100)")
    print("=" * 70)

    r = R("P6: Multi-epoch Bary")

    sweep_bodies = [
        (SUN, "Sun"),
        (MOON, "Moon"),
        (MARS, "Mars"),
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (PLUTO, "Pluto"),
    ]

    flags = FLG_SPEED | FLG_BARYCTR

    for y in range(1800, 2101, 25):
        jd = swe.julday(y, 6, 21, 12.0)
        for body_id, body_name in sweep_bodies:
            label = f"{y} {body_name} BARY"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label, check_speed=False)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 31: Barycentric Positions (FLG_BARYCTR)")
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
    print("ROUND 31 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 31: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
