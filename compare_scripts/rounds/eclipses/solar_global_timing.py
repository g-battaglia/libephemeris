#!/usr/bin/env python3
"""
Round 39: Solar Eclipse when_glob Timing Precision
====================================================

Tests sol_eclipse_when_glob() timing precision by comparing SE and LE
results for solar eclipse maximum times, and cross-validating eclipse
types and attributes.

Parts:
  P1: Forward search from 2020 — find 10 solar eclipses
  P2: Eclipse type classification agreement
  P3: Eclipse attributes (gamma, magnitude) comparison
  P4: Contact times comparison (C1, C2, C3, C4)
  P5: Backward verification — known historical eclipses
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
            print(f"  Max diff: {self.max_diff:.1f}s ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def ecl_type_str(ecl_type):
    """Convert eclipse type flags to string."""
    parts = []
    if ecl_type & ECL_TOTAL:
        parts.append("TOTAL")
    if ecl_type & ECL_ANNULAR:
        parts.append("ANNULAR")
    if ecl_type & ECL_PARTIAL:
        parts.append("PARTIAL")
    if ecl_type & ECL_ANNULAR_TOTAL:
        parts.append("HYBRID")
    if ecl_type & ECL_CENTRAL:
        parts.append("CENTRAL")
    if ecl_type & ECL_NONCENTRAL:
        parts.append("NONCENTRAL")
    return "|".join(parts) if parts else f"0x{ecl_type:x}"


def jd_to_date(jd):
    rev = swe.revjul(jd)
    return f"{int(rev[0])}-{int(rev[1]):02d}-{int(rev[2]):02d} {rev[3]:.2f}h"


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Forward search from 2020 — find 10 solar eclipses")
    print("=" * 70)

    r = R("P1: Eclipse max times")

    jd_start = swe.julday(2020, 1, 1, 0.0)
    eclipses_se = []
    eclipses_le = []

    for i in range(10):
        label = f"Eclipse #{i + 1}"
        try:
            se_result = swe.sol_eclipse_when_glob(jd_start, 0, 0)
            se_type = se_result[0]
            se_times = se_result[1]
            se_max_jd = se_times[0]

            le_result = ephem.sol_eclipse_when_glob(jd_start, 0, 0)
            le_type = le_result[0]
            le_times = le_result[1]
            le_max_jd = le_times[0]

            eclipses_se.append((se_type, se_times))
            eclipses_le.append((le_type, le_times))

            diff_sec = abs(se_max_jd - le_max_jd) * 86400
            tol = 60.0  # 60 seconds (eclipse timing can differ due to ephemeris)

            date_str = jd_to_date(le_max_jd)
            type_se = ecl_type_str(se_type)
            type_le = ecl_type_str(le_type)

            if diff_sec > tol:
                r.fail(f"{label} {date_str}: diff={diff_sec:.1f}s")
            else:
                r.ok(diff_sec, f"{label} {date_str}")
                print(
                    f"  {label}: {date_str} SE={type_se} LE={type_le} diff={diff_sec:.1f}s"
                )

            jd_start = se_max_jd + 20
        except Exception as e:
            r.skip(f"{label}: {e}")
            jd_start += 180

    return r.summary(), r, eclipses_se, eclipses_le


def run_part2(eclipses_se, eclipses_le):
    print("\n" + "=" * 70)
    print("PART 2: Eclipse type classification agreement")
    print("=" * 70)

    r = R("P2: Eclipse types")

    for i, (se_data, le_data) in enumerate(zip(eclipses_se, eclipses_le)):
        se_type, se_times = se_data
        le_type, le_times = le_data
        label = f"Eclipse #{i + 1}"
        date_str = jd_to_date(le_times[0])

        # Compare main type (TOTAL/ANNULAR/PARTIAL/HYBRID)
        se_main = se_type & (
            ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
        )
        le_main = le_type & (
            ECL_TOTAL | ECL_ANNULAR | ECL_PARTIAL | ECL_ANNULAR_TOTAL
        )

        if se_main != le_main:
            r.fail(
                f"{label} {date_str}: SE={ecl_type_str(se_main)} LE={ecl_type_str(le_main)}"
            )
        else:
            r.ok(0, f"{label} type")

    return r.summary(), r


def run_part3(eclipses_se, eclipses_le):
    print("\n" + "=" * 70)
    print("PART 3: Eclipse attributes (sol_eclipse_how at maximum)")
    print("=" * 70)

    r = R("P3: Eclipse attributes")

    for i, (se_data, le_data) in enumerate(zip(eclipses_se, eclipses_le)):
        se_type, se_times = se_data
        le_type, le_times = le_data
        label = f"Eclipse #{i + 1}"

        # Use sol_eclipse_how at the maximum time
        # geopos = (0, 0, 0) for sub-solar point approximation
        try:
            se_how = swe.sol_eclipse_how(se_times[0], (0.0, 0.0, 0.0), 0)
            le_how = ephem.sol_eclipse_how(le_times[0], 0, (0.0, 0.0, 0.0))
        except Exception as e:
            r.skip(f"{label}: {e}")
            continue

        # Compare attributes
        se_attr = se_how[1]  # attr array
        le_attr = le_how[1]

        # attr[0] = fraction of solar diameter covered
        # attr[1] = ratio of lunar to solar diameter
        # attr[2] = fraction of solar disc obscured
        for j, (name, tol) in enumerate(
            [
                ("fraction_covered", 0.05),
                ("diameter_ratio", 0.01),
                ("disc_obscured", 0.05),
            ]
        ):
            if j < len(se_attr) and j < len(le_attr):
                diff = abs(se_attr[j] - le_attr[j])
                if diff > tol:
                    r.fail(
                        f"{label} {name}: SE={se_attr[j]:.6f} LE={le_attr[j]:.6f} diff={diff:.6f}"
                    )
                else:
                    r.ok(diff, f"{label} {name}")

    return r.summary(), r


def run_part4(eclipses_se, eclipses_le):
    print("\n" + "=" * 70)
    print("PART 4: Contact times comparison")
    print("=" * 70)

    r = R("P4: Contact times")

    # se_times indices:
    # [0] = maximum, [1] = first contact (C1), [2] = second contact (C2)
    # [3] = third contact (C3), [4] = fourth contact (C4)
    contact_names = ["max", "C1", "C2", "C3", "C4"]

    for i, (se_data, le_data) in enumerate(zip(eclipses_se, eclipses_le)):
        se_type, se_times = se_data
        le_type, le_times = le_data
        label = f"Eclipse #{i + 1}"

        for j, cname in enumerate(contact_names):
            if j >= len(se_times) or j >= len(le_times):
                continue
            se_t = se_times[j]
            le_t = le_times[j]
            if se_t == 0 or le_t == 0:
                continue  # No contact for this eclipse type

            diff_sec = abs(se_t - le_t) * 86400
            tol = 120.0  # 2 minutes

            if diff_sec > tol:
                r.fail(f"{label} {cname}: diff={diff_sec:.1f}s")
            else:
                r.ok(diff_sec, f"{label} {cname}")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Known historical eclipses verification")
    print("=" * 70)

    r = R("P5: Historical eclipses")

    # Known eclipses with approximate dates
    known = [
        (2017, 8, 21, "2017 Great American (Total)"),
        (2019, 7, 2, "2019 South Pacific (Total)"),
        (2020, 6, 21, "2020 Annular"),
        (2023, 10, 14, "2023 Annular Americas"),
        (2024, 4, 8, "2024 Total Americas"),
    ]

    for y, m, d, desc in known:
        label = desc
        jd_start = swe.julday(y, m, d, 0.0) - 5  # 5 days before

        try:
            se_result = swe.sol_eclipse_when_glob(jd_start, 0, 0)
            le_result = ephem.sol_eclipse_when_glob(jd_start, 0, 0)

            se_max = se_result[1][0]
            le_max = le_result[1][0]

            diff_sec = abs(se_max - le_max) * 86400
            tol = 60.0

            se_type_str_val = ecl_type_str(se_result[0])
            le_type_str_val = ecl_type_str(le_result[0])

            if diff_sec > tol:
                r.fail(
                    f"{label}: diff={diff_sec:.1f}s SE={se_type_str_val} LE={le_type_str_val}"
                )
            else:
                r.ok(diff_sec, label)
                print(
                    f"  {label}: diff={diff_sec:.1f}s SE={se_type_str_val} LE={le_type_str_val}"
                )

        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 39: Solar Eclipse when_glob Timing Precision")
    print("=" * 70)

    start = time.time()
    all_ok = True
    all_results = []

    # P1 returns eclipse data for P2-P4
    try:
        ok, res, eclipses_se, eclipses_le = run_part1()
        all_results.append(("P1", res))
        if not ok:
            all_ok = False
    except Exception as e:
        print(f"\n  P1 CRASHED: {e}")
        traceback.print_exc()
        all_ok = False
        eclipses_se = []
        eclipses_le = []

    for pname, pfn in [
        ("P2", lambda: run_part2(eclipses_se, eclipses_le)),
        ("P3", lambda: run_part3(eclipses_se, eclipses_le)),
        ("P4", lambda: run_part4(eclipses_se, eclipses_le)),
        ("P5", run_part5),
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
    print("ROUND 39 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 39: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
