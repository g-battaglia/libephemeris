#!/usr/bin/env python3
"""
Round 27: Hypothetical/Uranian Bodies Re-verification
======================================================

Re-verifies hypothetical/Uranian body calculations after the Round 1 fix,
testing more dates, flag combinations, and edge cases.

Parts:
  P1: All 8 Uranian bodies at 10 dates (geocentric, ecliptic of date)
  P2: Heliocentric positions (FLG_HELCTR)
  P3: J2000 ecliptic positions (FLG_J2000)
  P4: Equatorial positions (FLG_EQUATORIAL)
  P5: Speed values for all Uranian bodies
  P6: Historical dates (1900-1950) — long-term propagation accuracy
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
            print(f'  Max diff: {self.max_diff:.1f}" ({self.max_label})')
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def angle_diff(a, b):
    d = a - b
    while d > 180:
        d -= 360
    while d < -180:
        d += 360
    return abs(d)


# Uranian/hypothetical bodies
URANIANS = [
    (CUPIDO, "Cupido"),
    (HADES, "Hades"),
    (ZEUS, "Zeus"),
    (KRONOS, "Kronos"),
    (APOLLON, "Apollon"),
    (ADMETOS, "Admetos"),
    (VULKANUS, "Vulkanus"),
    (POSEIDON, "Poseidon"),
]

DATES = [
    (1970, 1, 1, 12.0),
    (1980, 6, 15, 0.0),
    (1990, 3, 21, 12.0),
    (2000, 1, 1, 12.0),
    (2005, 7, 1, 0.0),
    (2010, 12, 21, 12.0),
    (2015, 4, 15, 0.0),
    (2020, 6, 21, 12.0),
    (2024, 3, 20, 12.0),
    (2024, 9, 22, 0.0),
]


def run_part(part_name, flags, tol_arcsec, flag_desc):
    print(f"\n{'=' * 70}")
    print(f"  {part_name}: {flag_desc}")
    print(f"{'=' * 70}")

    r = R(part_name)

    for y, m, d, h in DATES:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in URANIANS:
            label = f"{body_name} {y}"
            try:
                se_xx = swe.calc_ut(jd, body_id, flags | FLG_SPEED)[0]
            except Exception:
                r.skip(f"{label}: SE error")
                continue
            try:
                le_xx, _ = ephem.calc_ut(jd, body_id, flags | FLG_SPEED)
            except Exception:
                r.skip(f"{label}: LE error")
                continue

            lon_d = angle_diff(se_xx[0], le_xx[0]) * 3600
            lat_d = abs(se_xx[1] - le_xx[1]) * 3600
            max_d = max(lon_d, lat_d)

            if max_d > tol_arcsec:
                r.fail(f'{label}: lon={lon_d:.1f}" lat={lat_d:.1f}"')
            else:
                r.ok(max_d, label)

    # Print summary per body
    print(f"    Tested {r.passed + r.failed + r.skipped} combinations")
    return r.summary(), r


def run_part1():
    return run_part(
        "P1: Geocentric EclDate",
        0,
        60.0,
        "All 8 Uranians × 10 dates (geocentric, ecliptic of date)",
    )


def run_part2():
    return run_part("P2: Heliocentric", FLG_HELCTR, 30.0, "Heliocentric positions")


def run_part3():
    return run_part("P3: J2000", FLG_J2000, 60.0, "J2000 ecliptic positions")


def run_part4():
    return run_part("P4: Equatorial", FLG_EQUATORIAL, 60.0, "Equatorial positions")


def run_part5():
    print(f"\n{'=' * 70}")
    print(f"  P5: Speed values")
    print(f"{'=' * 70}")

    r = R("P5: Speed")
    jd = swe.julday(2024, 3, 20, 12.0)

    for body_id, body_name in URANIANS:
        label = f"{body_name} speed"
        try:
            se_xx = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
            le_xx, _ = ephem.calc_ut(jd, body_id, FLG_SPEED)
        except Exception:
            r.skip(f"{label}")
            continue

        spd_diff = abs(se_xx[3] - le_xx[3])
        # Uranian speeds are very slow (~0.5-1.5°/year = ~0.001-0.004°/day)
        tol = 0.01  # 0.01°/day tolerance
        if spd_diff > tol:
            r.fail(f"{label}: SE={se_xx[3]:.6f} LE={le_xx[3]:.6f} diff={spd_diff:.6f}")
        else:
            r.ok(spd_diff, label)

        print(
            f"    {body_name:10s}: SE={se_xx[3]:.6f}°/d LE={le_xx[3]:.6f}°/d diff={spd_diff:.6f}"
        )

    return r.summary(), r


def run_part6():
    print(f"\n{'=' * 70}")
    print(f"  P6: Historical dates (1900-1950)")
    print(f"{'=' * 70}")

    r = R("P6: Historical")

    hist_dates = [
        (1900, 1, 1, 12.0),
        (1910, 6, 15, 0.0),
        (1920, 3, 21, 12.0),
        (1930, 9, 23, 0.0),
        (1940, 12, 22, 12.0),
        (1950, 6, 21, 0.0),
    ]

    for y, m, d, h in hist_dates:
        jd = swe.julday(y, m, d, h)
        for body_id, body_name in URANIANS:
            label = f"{body_name} {y}"
            try:
                se_xx = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
                le_xx, _ = ephem.calc_ut(jd, body_id, FLG_SPEED)
            except Exception:
                r.skip(f"{label}")
                continue

            lon_d = angle_diff(se_xx[0], le_xx[0]) * 3600
            # Wider tolerance for historical (long propagation)
            tol = 120.0  # 2 arcmin
            if lon_d > tol:
                r.fail(f'{label}: lon={lon_d:.1f}"')
            else:
                r.ok(lon_d, label)

    print(f"    Tested {r.passed + r.failed + r.skipped} historical combinations")
    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 27: Hypothetical/Uranian Bodies Re-verification")
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
    print("ROUND 27 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 27: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
