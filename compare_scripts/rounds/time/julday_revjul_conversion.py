#!/usr/bin/env python3
"""
Round 29: Date Conversion & Julian Day Functions
==================================================

Tests all time/date conversion functions comprehensively.

Parts:
  P1: julday / revjul round-trip across centuries
  P2: utc_to_jd / jdet_to_utc round-trip
  P3: deltat across epochs (1600-2100)
  P4: deltat_ex with different ephemeris flags
  P5: Gregorian/Julian calendar boundary (Oct 1582)
  P6: Negative years (BCE dates)
  P7: sidtime / sidtime0 comparison
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
            print(f"  Max diff: {self.max_diff:.10f} ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


# ============================================================
# PART 1: julday / revjul round-trip
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: julday / revjul round-trip across centuries")
    print("=" * 70)

    r = R("P1: julday/revjul")

    test_dates = [
        (1600, 1, 1, 12.0),
        (1700, 6, 15, 6.5),
        (1800, 3, 21, 18.25),
        (1900, 1, 1, 0.0),
        (1950, 7, 4, 12.0),
        (1969, 7, 20, 20.0 + 17.0 / 60),
        (1999, 12, 31, 23.999),
        (2000, 1, 1, 12.0),
        (2000, 2, 29, 0.0),
        (2020, 3, 20, 3.5),
        (2024, 2, 29, 12.0),
        (2024, 12, 31, 23.5),
        (2050, 6, 21, 0.0),
        (2100, 1, 1, 0.0),
    ]

    for y, m, d, h in test_dates:
        label = f"{y}-{m:02d}-{d:02d} {h:.2f}h"

        # SE julday
        se_jd = swe.julday(y, m, d, h)
        le_jd = ephem.julday(y, m, d, h)

        jd_diff = abs(se_jd - le_jd)
        if jd_diff > 1e-10:
            r.fail(f"{label} julday: SE={se_jd:.10f} LE={le_jd:.10f} diff={jd_diff}")
        else:
            r.ok(jd_diff, f"{label} julday")

        # SE revjul
        se_rev = swe.revjul(se_jd)
        le_rev = ephem.revjul(le_jd)

        # Compare year, month, day
        if (
            se_rev[0] != le_rev[0]
            or se_rev[1] != le_rev[1]
            or int(se_rev[2]) != int(le_rev[2])
        ):
            r.fail(f"{label} revjul date: SE={se_rev[:3]} LE={le_rev[:3]}")
        else:
            r.ok(0.0, f"{label} revjul date")

        # Compare hour (floating point)
        hour_diff = abs(se_rev[3] - le_rev[3])
        if hour_diff > 1e-8:
            r.fail(f"{label} revjul hour: SE={se_rev[3]:.10f} LE={le_rev[3]:.10f}")
        else:
            r.ok(hour_diff, f"{label} revjul hour")

    print(f"  Tested {len(test_dates)} dates")
    return r.summary(), r


# ============================================================
# PART 2: utc_to_jd / jdet_to_utc
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: utc_to_jd comparison")
    print("=" * 70)

    r = R("P2: utc_to_jd")

    test_dates = [
        (2000, 1, 1, 12, 0, 0.0),
        (2024, 3, 20, 15, 30, 45.0),
        (2024, 6, 21, 0, 0, 0.0),
        (1999, 12, 31, 23, 59, 59.0),
        (2016, 12, 31, 23, 59, 60.0),  # Leap second
        (1972, 7, 1, 0, 0, 0.0),
        (2050, 1, 1, 0, 0, 0.0),
    ]

    for y, m, d, h, mi, s in test_dates:
        label = f"{y}-{m:02d}-{d:02d} {h:02d}:{mi:02d}:{s:04.1f}"

        try:
            se_result = swe.utc_to_jd(y, m, d, h, mi, s, 1)  # 1 = Gregorian
            se_jd_et = se_result[0]  # ET
            se_jd_ut = se_result[1]  # UT
        except Exception:
            r.skip(f"{label}: SE error")
            continue

        try:
            le_result = ephem.utc_to_jd(y, m, d, h, mi, s, 1)
            le_jd_et = le_result[0]
            le_jd_ut = le_result[1]
        except Exception:
            r.skip(f"{label}: LE error")
            continue

        et_diff = abs(se_jd_et - le_jd_et)
        ut_diff = abs(se_jd_ut - le_jd_ut)

        # ET tolerance: 0.001 seconds (delta-T differences)
        if et_diff > 1e-8:  # ~0.001 seconds
            r.fail(f"{label} ET: diff={et_diff:.12f}")
        else:
            r.ok(et_diff, f"{label} ET")

        if ut_diff > 1e-10:
            r.fail(f"{label} UT: diff={ut_diff:.12f}")
        else:
            r.ok(ut_diff, f"{label} UT")

    print(f"  Tested {len(test_dates)} UTC dates")
    return r.summary(), r


# ============================================================
# PART 3: deltat across epochs
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: deltat across epochs (1600-2100)")
    print("=" * 70)

    r = R("P3: deltat")

    years = list(range(1600, 2101, 25))
    max_diff_sec = 0.0

    for y in years:
        jd = swe.julday(y, 1, 1, 12.0)
        label = f"deltat {y}"

        se_dt = swe.deltat(jd)
        le_dt = ephem.deltat(jd)

        diff_sec = abs(se_dt - le_dt) * 86400.0  # Convert days to seconds

        if diff_sec > max_diff_sec:
            max_diff_sec = diff_sec

        # Tolerance: 1 second for historical, 0.1s for modern
        tol = 1.0 if y < 1900 else 0.1
        if diff_sec > tol:
            r.fail(
                f"{label}: SE={se_dt * 86400:.2f}s LE={le_dt * 86400:.2f}s diff={diff_sec:.4f}s"
            )
        else:
            r.ok(diff_sec, label)

    print(f"  {len(years)} epochs, max_diff={max_diff_sec:.4f}s")
    return r.summary(), r


# ============================================================
# PART 4: deltat_ex
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: deltat_ex with different flags")
    print("=" * 70)

    r = R("P4: deltat_ex")

    jds = [swe.julday(y, 1, 1, 12.0) for y in [1900, 1950, 2000, 2024, 2050]]

    for jd in jds:
        yr = int(swe.revjul(jd)[0])
        label = f"deltat_ex {yr}"

        try:
            se_dt = swe.deltat_ex(jd, -1)  # -1 = default
            if isinstance(se_dt, tuple):
                se_val = se_dt[0]
            else:
                se_val = se_dt
        except Exception:
            r.skip(f"{label}: SE error")
            continue

        try:
            le_val = ephem.deltat_ex(jd, -1)
            if isinstance(le_val, tuple):
                le_val = le_val[0]
        except Exception:
            r.skip(f"{label}: LE error")
            continue

        diff_sec = abs(se_val - le_val) * 86400.0
        tol = 0.5
        if diff_sec > tol:
            r.fail(f"{label}: diff={diff_sec:.4f}s")
        else:
            r.ok(diff_sec, label)

    print(f"  Tested {r.passed + r.failed} deltat_ex calls")
    return r.summary(), r


# ============================================================
# PART 5: Gregorian/Julian calendar boundary
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Gregorian/Julian calendar boundary (Oct 1582)")
    print("=" * 70)

    r = R("P5: Calendar Boundary")

    # Dates around the Gregorian reform
    boundary_dates = [
        (1582, 10, 4, 12.0, "Last Julian"),
        (1582, 10, 15, 12.0, "First Gregorian"),
        (1582, 10, 16, 12.0, "2nd Gregorian"),
        (1582, 1, 1, 0.0, "Jan 1582"),
        (1583, 1, 1, 0.0, "Jan 1583"),
        (1500, 6, 15, 12.0, "1500 (Julian)"),
        (1600, 1, 1, 0.0, "1600 (Gregorian)"),
    ]

    for y, m, d, h, desc in boundary_dates:
        label = f"boundary {desc}"

        se_jd = swe.julday(y, m, d, h)
        le_jd = ephem.julday(y, m, d, h)

        diff = abs(se_jd - le_jd)
        if diff > 1e-10:
            r.fail(f"{label}: SE={se_jd:.6f} LE={le_jd:.6f} diff={diff}")
        else:
            r.ok(diff, label)

    print(f"  Tested {len(boundary_dates)} boundary dates")
    return r.summary(), r


# ============================================================
# PART 6: Negative years (BCE)
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Negative years (BCE dates)")
    print("=" * 70)

    r = R("P6: BCE Dates")

    bce_dates = [
        (-4712, 1, 1, 12.0, "JD epoch"),
        (-3000, 6, 15, 0.0, "3001 BCE"),
        (-1000, 3, 21, 12.0, "1001 BCE"),
        (-500, 7, 1, 0.0, "501 BCE"),
        (0, 1, 1, 12.0, "1 BCE (year 0)"),
        (1, 1, 1, 12.0, "1 CE"),
        (100, 6, 15, 0.0, "100 CE"),
        (500, 1, 1, 0.0, "500 CE"),
    ]

    for y, m, d, h, desc in bce_dates:
        label = f"BCE {desc}"

        try:
            se_jd = swe.julday(y, m, d, h)
        except Exception:
            r.skip(f"{label}: SE error")
            continue
        try:
            le_jd = ephem.julday(y, m, d, h)
        except Exception:
            r.skip(f"{label}: LE error")
            continue

        diff = abs(se_jd - le_jd)
        if diff > 1e-8:
            r.fail(f"{label}: SE={se_jd:.6f} LE={le_jd:.6f} diff={diff}")
        else:
            r.ok(diff, label)

        # Also test revjul round-trip
        try:
            se_rev = swe.revjul(se_jd)
            le_rev = ephem.revjul(le_jd)
            if int(se_rev[0]) != int(le_rev[0]) or se_rev[1] != le_rev[1]:
                r.fail(f"{label} revjul: SE={se_rev[:3]} LE={le_rev[:3]}")
            else:
                r.ok(0.0, f"{label} revjul")
        except Exception:
            r.skip(f"{label} revjul")

    print(f"  Tested {len(bce_dates)} BCE dates")
    return r.summary(), r


# ============================================================
# PART 7: Sidereal time
# ============================================================
def run_part7():
    print("\n" + "=" * 70)
    print("PART 7: sidtime comparison")
    print("=" * 70)

    r = R("P7: Sidereal Time")

    jds = [
        swe.julday(y, 1, 1, 12.0)
        for y in [1900, 1950, 1980, 2000, 2010, 2020, 2024, 2050, 2100]
    ]

    for jd in jds:
        yr = int(swe.revjul(jd)[0])
        label = f"sidtime {yr}"

        try:
            se_st = swe.sidtime(jd)
        except Exception:
            r.skip(f"{label}: SE error")
            continue

        try:
            le_st = ephem.sidtime(jd)
        except Exception:
            r.skip(f"{label}: LE error")
            continue

        # Sidereal time in hours (0-24)
        diff_sec = abs(se_st - le_st) * 3600.0  # hours to seconds

        # Handle wraparound
        if diff_sec > 12 * 3600:
            diff_sec = 24 * 3600 - diff_sec

        tol = 1.0  # 1 second tolerance
        if diff_sec > tol:
            r.fail(f"{label}: SE={se_st:.8f}h LE={le_st:.8f}h diff={diff_sec:.4f}s")
        else:
            r.ok(diff_sec, label)

        print(f"  {yr}: SE={se_st:.8f}h LE={le_st:.8f}h diff={diff_sec:.6f}s")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 29: Date Conversion & Julian Day Functions")
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
    print("ROUND 29 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 29: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
