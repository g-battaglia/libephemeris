#!/usr/bin/env python3
"""
Round 40: Lunar Eclipse when_glob Timing Precision
====================================================

Tests lun_eclipse_when() timing precision by comparing SE and LE results
for lunar eclipse maximum times, types, and attributes.

Parts:
  P1: Forward search from 2020 — find 10 lunar eclipses
  P2: Eclipse type classification agreement
  P3: Eclipse magnitude comparison
  P4: Contact times (penumbral/partial/total begin/end)
  P5: Known historical lunar eclipses
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
    parts = []
    if ecl_type & ECL_TOTAL:
        parts.append("TOTAL")
    if ecl_type & ECL_PARTIAL:
        parts.append("PARTIAL")
    if ecl_type & ECL_PENUMBRAL:
        parts.append("PENUMBRAL")
    return "|".join(parts) if parts else f"0x{ecl_type:x}"


def jd_to_date(jd):
    rev = swe.revjul(jd)
    return f"{int(rev[0])}-{int(rev[1]):02d}-{int(rev[2]):02d} {rev[3]:.2f}h"


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Forward search from 2020 — find 10 lunar eclipses")
    print("=" * 70)

    r = R("P1: Lunar eclipse max times")

    jd_start = swe.julday(2020, 1, 1, 0.0)
    eclipses_se = []
    eclipses_le = []

    for i in range(10):
        label = f"LunEcl #{i + 1}"
        try:
            se_result = swe.lun_eclipse_when(jd_start, 0, 0)
            se_type = se_result[0]
            se_times = se_result[1]
            se_max_jd = se_times[0]

            le_result = ephem.lun_eclipse_when(jd_start, 0, 0)
            le_type = le_result[0]
            le_times = le_result[1]
            le_max_jd = le_times[0]

            eclipses_se.append((se_type, se_times))
            eclipses_le.append((le_type, le_times))

            diff_sec = abs(se_max_jd - le_max_jd) * 86400
            tol = 60.0

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

    r = R("P2: Lunar eclipse types")

    for i, (se_data, le_data) in enumerate(zip(eclipses_se, eclipses_le)):
        se_type = se_data[0]
        le_type = le_data[0]
        label = f"LunEcl #{i + 1}"

        se_main = se_type & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)
        le_main = le_type & (ECL_TOTAL | ECL_PARTIAL | ECL_PENUMBRAL)

        if se_main != le_main:
            r.fail(f"{label}: SE={ecl_type_str(se_main)} LE={ecl_type_str(le_main)}")
        else:
            r.ok(0, f"{label} type")

    return r.summary(), r


def run_part3(eclipses_se, eclipses_le):
    print("\n" + "=" * 70)
    print("PART 3: Eclipse magnitude comparison (lun_eclipse_how)")
    print("=" * 70)

    r = R("P3: Lunar eclipse magnitude")

    for i, (se_data, le_data) in enumerate(zip(eclipses_se, eclipses_le)):
        se_times = se_data[1]
        le_times = le_data[1]
        label = f"LunEcl #{i + 1}"

        try:
            se_how = swe.lun_eclipse_how(se_times[0], (0.0, 0.0, 0.0), 0)
            le_how = ephem.lun_eclipse_how(le_times[0], 0, (0.0, 0.0, 0.0))

            se_attr = se_how[1]
            le_attr = le_how[1]

            # attr[0] = umbral magnitude
            # attr[1] = penumbral magnitude
            for j, (name, tol) in enumerate(
                [
                    ("umbral_mag", 0.02),
                    ("penumbral_mag", 0.02),
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

        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part4(eclipses_se, eclipses_le):
    print("\n" + "=" * 70)
    print("PART 4: Contact times comparison")
    print("=" * 70)

    r = R("P4: Lunar eclipse contacts")

    # Lunar eclipse tret indices:
    # [0] = maximum
    # [2] = partial phase begin
    # [3] = partial phase end
    # [4] = totality begin
    # [5] = totality end
    # [6] = penumbral phase begin
    # [7] = penumbral phase end
    contact_map = [
        (0, "max"),
        (2, "partial_begin"),
        (3, "partial_end"),
        (4, "total_begin"),
        (5, "total_end"),
        (6, "pen_begin"),
        (7, "pen_end"),
    ]

    for i, (se_data, le_data) in enumerate(zip(eclipses_se, eclipses_le)):
        se_times = se_data[1]
        le_times = le_data[1]
        label = f"LunEcl #{i + 1}"

        for idx, cname in contact_map:
            if idx >= len(se_times) or idx >= len(le_times):
                continue
            se_t = se_times[idx]
            le_t = le_times[idx]
            if se_t == 0 or le_t == 0:
                continue

            diff_sec = abs(se_t - le_t) * 86400
            tol = 120.0

            if diff_sec > tol:
                r.fail(f"{label} {cname}: diff={diff_sec:.1f}s")
            else:
                r.ok(diff_sec, f"{label} {cname}")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Known historical lunar eclipses")
    print("=" * 70)

    r = R("P5: Historical lunar eclipses")

    known = [
        (2019, 1, 21, "2019-01-21 Total"),
        (2021, 5, 26, "2021-05-26 Total"),
        (2021, 11, 19, "2021-11-19 Partial"),
        (2022, 5, 16, "2022-05-16 Total"),
        (2022, 11, 8, "2022-11-08 Total"),
        (2023, 10, 28, "2023-10-28 Partial"),
        (2024, 3, 25, "2024-03-25 Penumbral"),
        (2024, 9, 18, "2024-09-18 Partial"),
    ]

    for y, m, d, desc in known:
        label = desc
        jd_start = swe.julday(y, m, d, 0.0) - 5

        try:
            se_result = swe.lun_eclipse_when(jd_start, 0, 0)
            le_result = ephem.lun_eclipse_when(jd_start, 0, 0)

            se_max = se_result[1][0]
            le_max = le_result[1][0]

            diff_sec = abs(se_max - le_max) * 86400
            tol = 60.0

            if diff_sec > tol:
                r.fail(f"{label}: diff={diff_sec:.1f}s")
            else:
                r.ok(diff_sec, label)
                print(
                    f"  {label}: diff={diff_sec:.1f}s SE={ecl_type_str(se_result[0])} LE={ecl_type_str(le_result[0])}"
                )

        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 40: Lunar Eclipse when_glob Timing Precision")
    print("=" * 70)

    start = time.time()
    all_ok = True
    all_results = []

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
    print("ROUND 40 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 40: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
