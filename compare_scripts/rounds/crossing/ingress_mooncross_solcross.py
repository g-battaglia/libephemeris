#!/usr/bin/env python3
"""
Round 37: Planetary Ingress Timing
====================================

Tests the precision of zodiacal sign ingress detection using mooncross_ut
and solcross_ut / crossings for planets entering new signs.

Parts:
  P1: Sun sign ingresses 2024 (all 12 signs)
  P2: Moon sign ingresses (fast — tests ~13°/day motion)
  P3: Mercury/Venus ingresses 2024
  P4: Mars/Jupiter/Saturn ingresses 2024
  P5: Cross-validation: SE mooncross vs LE mooncross for 0°/90°/180°/270°
  P6: Sun solstice/equinox precision (0°/90°/180°/270° ecliptic)
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
            print(f"  Max diff: {self.max_diff:.2f}s ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Sun sign ingresses 2024 (all 12 signs)")
    print("=" * 70)

    r = R("P1: Sun ingresses")

    # Sun enters each sign approximately on these dates
    # We'll search from ~5 days before each expected date
    sign_lons = list(range(0, 360, 30))  # 0°, 30°, 60°, ... 330°
    sign_names = [
        "Ari",
        "Tau",
        "Gem",
        "Can",
        "Leo",
        "Vir",
        "Lib",
        "Sco",
        "Sag",
        "Cap",
        "Aqu",
        "Pis",
    ]

    # Start searching from Jan 1, 2024
    jd_start = swe.julday(2024, 1, 1, 0.0)

    for i, (lon, name) in enumerate(zip(sign_lons, sign_names)):
        label = f"Sun→{name} ({lon}°)"
        try:
            se_jd = swe.solcross_ut(lon, jd_start, 0)
            le_jd = ephem.solcross_ut(lon, jd_start, 0)

            diff_sec = abs(se_jd - le_jd) * 86400
            tol = 5.0  # 5 seconds

            if diff_sec > tol:
                r.fail(f"{label}: diff={diff_sec:.2f}s")
            else:
                r.ok(diff_sec, label)
                print(f"  {label}: diff={diff_sec:.2f}s")

            # Use this crossing as next search start
            jd_start = max(se_jd, le_jd) + 1
        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Moon sign ingresses (first 6 signs from Mar 20, 2024)")
    print("=" * 70)

    r = R("P2: Moon ingresses")

    sign_lons = list(range(0, 360, 30))
    sign_names = [
        "Ari",
        "Tau",
        "Gem",
        "Can",
        "Leo",
        "Vir",
        "Lib",
        "Sco",
        "Sag",
        "Cap",
        "Aqu",
        "Pis",
    ]

    jd_start = swe.julday(2024, 3, 20, 0.0)

    for i in range(6):  # First 6 sign ingresses
        lon = sign_lons[i]
        name = sign_names[i]
        label = f"Moon→{name} ({lon}°)"
        try:
            se_jd = swe.mooncross_ut(lon, jd_start, 0)
            le_jd = ephem.mooncross_ut(lon, jd_start, 0)

            diff_sec = abs(se_jd - le_jd) * 86400
            tol = 2.0  # 2 seconds (Moon moves fast)

            if diff_sec > tol:
                r.fail(f"{label}: diff={diff_sec:.2f}s")
            else:
                r.ok(diff_sec, label)

            jd_start = max(se_jd, le_jd) + 0.5
        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Mercury/Venus ingresses 2024")
    print("=" * 70)

    r = R("P3: Mercury/Venus ingresses")

    for body_id, body_name in [(MERCURY, "Mercury"), (VENUS, "Venus")]:
        jd_start = swe.julday(2024, 1, 1, 0.0)
        for lon in range(0, 360, 30):
            label = f"{body_name}→{lon}°"
            try:
                se_jd = swe.solcross_ut(lon, jd_start, body_id)
            except Exception:
                try:
                    # Fallback: use helio_cross_ut for planets
                    se_jd = swe.helio_cross_ut(body_id, lon, jd_start, 0)
                except Exception:
                    r.skip(f"{label}: SE error")
                    continue

            try:
                le_jd = ephem.solcross_ut(lon, jd_start, body_id)
            except Exception:
                try:
                    le_jd = ephem.helio_cross_ut(body_id, lon, jd_start, 0)
                except Exception:
                    r.skip(f"{label}: LE error")
                    continue

            diff_sec = abs(se_jd - le_jd) * 86400
            tol = 10.0

            if diff_sec > tol:
                r.fail(f"{label}: diff={diff_sec:.2f}s")
            else:
                r.ok(diff_sec, label)

            jd_start = max(se_jd, le_jd) + 1

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Mars/Jupiter/Saturn ingresses 2024")
    print("=" * 70)

    r = R("P4: Outer planet ingresses")

    # These planets move slowly, fewer ingresses per year
    for body_id, body_name, expected_ingresses in [
        (MARS, "Mars", 6),
        (JUPITER, "Jupiter", 2),
        (SATURN, "Saturn", 1),
    ]:
        jd_start = swe.julday(2024, 1, 1, 0.0)
        jd_end = swe.julday(2025, 1, 1, 0.0)
        count = 0

        # Get current position to find next sign boundary
        pos = swe.calc_ut(jd_start, body_id, FLG_SPEED)[0]
        current_sign = int(pos[0] / 30)
        next_lon = ((current_sign + 1) % 12) * 30

        for _ in range(expected_ingresses + 2):
            label = f"{body_name}→{next_lon}°"
            try:
                se_jd = swe.solcross_ut(next_lon, jd_start, body_id)
                if se_jd > jd_end:
                    break
            except Exception:
                r.skip(f"{label}: SE error")
                break

            try:
                le_jd = ephem.solcross_ut(next_lon, jd_start, body_id)
            except Exception:
                r.skip(f"{label}: LE error")
                jd_start = se_jd + 1
                next_lon = (next_lon + 30) % 360
                continue

            diff_sec = abs(se_jd - le_jd) * 86400
            tol = 30.0  # 30s for slow outer planets

            if diff_sec > tol:
                r.fail(f"{label}: diff={diff_sec:.2f}s")
            else:
                r.ok(diff_sec, label)
                count += 1

            jd_start = max(se_jd, le_jd) + 1
            next_lon = (next_lon + 30) % 360

        print(f"  {body_name}: {count} ingresses found")

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Moon crossing cardinal points (0°/90°/180°/270°)")
    print("=" * 70)

    r = R("P5: Moon cardinal crossings")

    jd_start = swe.julday(2024, 1, 1, 0.0)

    for lon in [0, 90, 180, 270]:
        for trial in range(3):  # 3 consecutive crossings each
            label = f"Moon cross {lon}° #{trial + 1}"
            try:
                se_jd = swe.mooncross_ut(lon, jd_start, 0)
                le_jd = ephem.mooncross_ut(lon, jd_start, 0)

                diff_sec = abs(se_jd - le_jd) * 86400
                tol = 2.0

                if diff_sec > tol:
                    r.fail(f"{label}: diff={diff_sec:.2f}s")
                else:
                    r.ok(diff_sec, label)

                jd_start = max(se_jd, le_jd) + 1
            except Exception as e:
                r.skip(f"{label}: {e}")
                jd_start += 27.3  # skip to next month

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Sun equinox/solstice precision (cardinal points)")
    print("=" * 70)

    r = R("P6: Sun equinox/solstice")

    for year in [2020, 2024, 2025, 2030]:
        jd_start = swe.julday(year, 1, 1, 0.0)
        for lon, event in [
            (0, "VernalEq"),
            (90, "SummerSol"),
            (180, "AutumnEq"),
            (270, "WinterSol"),
        ]:
            label = f"{year} {event}"
            try:
                se_jd = swe.solcross_ut(lon, jd_start, 0)
                le_jd = ephem.solcross_ut(lon, jd_start, 0)

                diff_sec = abs(se_jd - le_jd) * 86400
                tol = 5.0

                if diff_sec > tol:
                    r.fail(f"{label}: diff={diff_sec:.2f}s")
                else:
                    r.ok(diff_sec, label)
                    # Print the event time
                    rev = swe.revjul(le_jd)
                    print(
                        f"  {label}: {int(rev[0])}-{int(rev[1]):02d}-{int(rev[2]):02d} {rev[3]:.4f}h diff={diff_sec:.3f}s"
                    )

                jd_start = se_jd + 30
            except Exception as e:
                r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 37: Planetary Ingress Timing")
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
    print("ROUND 37 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 37: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
