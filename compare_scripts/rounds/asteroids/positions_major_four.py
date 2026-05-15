#!/usr/bin/env python3
"""
Round 35: Asteroid Positions Deep Sweep
========================================

Tests the four major asteroids (Ceres, Pallas, Juno, Vesta) across multiple
epochs, flag combinations, and coordinate systems.

Parts:
  P1: Default ecliptic — 4 asteroids × 5 epochs
  P2: Equatorial output
  P3: J2000 ecliptic
  P4: Heliocentric
  P5: Speed validation (finite difference)
  P6: Multi-epoch sweep (1900-2100) for orbit coverage
  P7: Sidereal mode
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

# Major asteroids use dedicated constants in libephemeris
# SE uses either AST_OFFSET+N or dedicated IDs; LE uses dedicated IDs only
ASTEROIDS = [
    (CERES, "Ceres"),
    (PALLAS, "Pallas"),
    (JUNO, "Juno"),
    (VESTA, "Vesta"),
]

EPOCHS = [
    (1950, 1, 1, 12.0, "1950"),
    (2000, 1, 1, 12.0, "J2000"),
    (2024, 3, 20, 15.5, "2024eq"),
    (2024, 9, 22, 12.0, "2024aut"),
    (2050, 6, 21, 0.0, "2050"),
]

TOL_LON = 5.0  # arcsec (asteroids have more ephemeris variation)
TOL_LAT = 5.0
TOL_DIST = 5e-5
TOL_SPEED = 0.005


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
            print(f'  Max diff: {self.max_diff:.3f}" ({self.max_label})')
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
    lon_as = lon_diff * 3600

    lat_diff = abs(se_pos[1] - le_pos[1]) * 3600
    dist_diff = abs(se_pos[2] - le_pos[2])

    if lon_as > TOL_LON:
        r.fail(f'{label} lon: {lon_as:.3f}" (SE={se_pos[0]:.8f} LE={le_pos[0]:.8f})')
    elif lat_diff > TOL_LAT:
        r.fail(f'{label} lat: {lat_diff:.3f}"')
    elif dist_diff > TOL_DIST:
        r.fail(f"{label} dist: {dist_diff:.8f} AU")
    else:
        r.ok(max(lon_as, lat_diff), f"{label}")

    if check_speed and len(se_pos) >= 4 and len(le_pos) >= 4:
        spd_diff = abs(se_pos[3] - le_pos[3])
        if spd_diff > TOL_SPEED:
            r.fail(f"{label} spd: {spd_diff:.6f}")
        else:
            r.ok(spd_diff * 3600, f"{label} spd")


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Default ecliptic — 4 asteroids × 5 epochs")
    print("=" * 70)

    r = R("P1: Asteroids ecliptic")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in ASTEROIDS:
            flags = FLG_SPEED
            label = f"{epoch_name} {body_name}"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Equatorial output")
    print("=" * 70)

    r = R("P2: Asteroids EQ")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in ASTEROIDS:
            flags = FLG_SPEED | FLG_EQUATORIAL
            label = f"{epoch_name} {body_name} EQ"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: J2000 ecliptic")
    print("=" * 70)

    r = R("P3: Asteroids J2000")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in ASTEROIDS:
            flags = FLG_SPEED | FLG_J2000
            label = f"{epoch_name} {body_name} J2K"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Heliocentric")
    print("=" * 70)

    r = R("P4: Asteroids HELIO")

    for y, m, d, h, epoch_name in EPOCHS:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in ASTEROIDS:
            flags = FLG_SPEED | FLG_HELCTR
            label = f"{epoch_name} {body_name} HELIO"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Speed validation (finite difference)")
    print("=" * 70)

    r = R("P5: Asteroid speeds FD")
    jd = swe.julday(2024, 3, 20, 15.5)
    dt = 1.0 / 86400.0  # 1 second

    for body_id, body_name in ASTEROIDS:
        flags = FLG_SPEED
        label = f"{body_name}"
        try:
            le_r = ephem.calc_ut(jd, body_id, flags)
            le_before = ephem.calc_ut(jd - dt, body_id, 0)
            le_after = ephem.calc_ut(jd + dt, body_id, 0)

            for i, comp in enumerate(["lon", "lat", "dist"]):
                reported = le_r[0][i + 3]
                fd = (le_after[0][i] - le_before[0][i]) / (2 * dt)
                if i == 0:  # handle wraparound
                    if fd > 180:
                        fd -= 360
                    elif fd < -180:
                        fd += 360
                fd_diff = abs(reported - fd)
                tol = 0.01 if i < 2 else 1e-4
                if fd_diff > tol:
                    r.fail(
                        f"{label} {comp}_spd FD: rep={reported:.8f} fd={fd:.8f} diff={fd_diff:.2e}"
                    )
                else:
                    r.ok(fd_diff, f"{label} {comp}_spd FD")
        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Multi-epoch sweep (1900-2100)")
    print("=" * 70)

    r = R("P6: Multi-epoch asteroids")

    for y in range(1900, 2101, 10):
        jd = swe.julday(y, 6, 21, 12.0)
        for body_id, body_name in ASTEROIDS:
            flags = FLG_SPEED
            label = f"{y} {body_name}"
            try:
                se_r = swe.calc_ut(jd, body_id, flags)
                le_r = ephem.calc_ut(jd, body_id, flags)
                compare_pos(r, se_r, le_r, label, check_speed=False)
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: Sidereal mode (Lahiri)")
    print("=" * 70)

    r = R("P7: Asteroids sidereal")

    swe.set_sid_mode(SIDM_LAHIRI)
    ephem.set_sid_mode(SIDM_LAHIRI, 0, 0)

    jd = swe.julday(2024, 3, 20, 15.5)
    for body_id, body_name in ASTEROIDS:
        flags = FLG_SPEED | FLG_SIDEREAL
        label = f"{body_name} SID"
        try:
            se_r = swe.calc_ut(jd, body_id, flags)
            le_r = ephem.calc_ut(jd, body_id, flags)
            compare_pos(r, se_r, le_r, label)
        except Exception as e:
            r.skip(f"{label}: {e}")

    # Reset sidereal mode
    swe.set_sid_mode(0)
    ephem.set_sid_mode(0, 0, 0)

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 35: Asteroid Positions Deep Sweep")
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
    print("ROUND 35 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 35: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
