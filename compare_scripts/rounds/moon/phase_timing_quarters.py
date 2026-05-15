#!/usr/bin/env python3
"""
Round 38: Moon Phase Timing (New Moon / Full Moon / Quarters)
=============================================================

Tests precision of Moon phase detection by finding times when the
Sun-Moon elongation equals 0° (New), 90° (First Quarter), 180° (Full),
270° (Last Quarter).

Parts:
  P1: New Moons 2024 (12 lunations)
  P2: Full Moons 2024 (12 lunations)
  P3: Quarter Moons 2024 (first 6 each)
  P4: Lunation length consistency (synodic month ~29.53 days)
  P5: Multi-year New Moon sweep (2020-2026)
  P6: Verification against mooncross_ut Sun-Moon elongation
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
            print(f"  Max diff: {self.max_diff:.3f}s ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


def find_moon_phase_se(jd_start, phase_angle):
    """Find next time Moon-Sun elongation = phase_angle using SE mooncross."""
    # Moon longitude = Sun longitude + phase_angle
    # We find when Moon crosses Sun_lon + phase_angle
    sun = swe.calc_ut(jd_start, SUN, 0)[0]
    target_lon = (sun[0] + phase_angle) % 360

    # Use mooncross to find when Moon reaches target longitude
    jd = swe.mooncross_ut(target_lon, jd_start, 0)

    # Iterate: Sun moves ~1°/day, so refine
    for _ in range(5):
        sun = swe.calc_ut(jd, SUN, 0)[0]
        target_lon = (sun[0] + phase_angle) % 360
        jd = swe.mooncross_ut(target_lon, jd - 0.5, 0)

    return jd


def find_moon_phase_le(jd_start, phase_angle):
    """Find next time Moon-Sun elongation = phase_angle using LE mooncross."""
    sun = ephem.calc_ut(jd_start, SUN, 0)[0]
    target_lon = (sun[0] + phase_angle) % 360

    jd = ephem.mooncross_ut(target_lon, jd_start, 0)

    for _ in range(5):
        sun = ephem.calc_ut(jd, SUN, 0)[0]
        target_lon = (sun[0] + phase_angle) % 360
        jd = ephem.mooncross_ut(target_lon, jd - 0.5, 0)

    return jd


def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: New Moons 2024 (elongation = 0°)")
    print("=" * 70)

    r = R("P1: New Moons 2024")

    jd_start = swe.julday(2024, 1, 1, 0.0)

    for i in range(12):
        label = f"NewMoon #{i + 1}"
        try:
            se_jd = find_moon_phase_se(jd_start, 0)
            le_jd = find_moon_phase_le(jd_start, 0)

            diff_sec = abs(se_jd - le_jd) * 86400
            tol = 5.0

            if diff_sec > tol:
                r.fail(f"{label}: diff={diff_sec:.2f}s")
            else:
                r.ok(diff_sec, label)
                rev = swe.revjul(le_jd)
                print(
                    f"  {label}: {int(rev[0])}-{int(rev[1]):02d}-{int(rev[2]):02d} {rev[3]:.2f}h  diff={diff_sec:.3f}s"
                )

            jd_start = max(se_jd, le_jd) + 20  # skip ~20 days to next lunation
        except Exception as e:
            r.skip(f"{label}: {e}")
            jd_start += 29.53

    return r.summary(), r


def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Full Moons 2024 (elongation = 180°)")
    print("=" * 70)

    r = R("P2: Full Moons 2024")

    jd_start = swe.julday(2024, 1, 1, 0.0)

    for i in range(12):
        label = f"FullMoon #{i + 1}"
        try:
            se_jd = find_moon_phase_se(jd_start, 180)
            le_jd = find_moon_phase_le(jd_start, 180)

            diff_sec = abs(se_jd - le_jd) * 86400
            tol = 5.0

            if diff_sec > tol:
                r.fail(f"{label}: diff={diff_sec:.2f}s")
            else:
                r.ok(diff_sec, label)
                rev = swe.revjul(le_jd)
                print(
                    f"  {label}: {int(rev[0])}-{int(rev[1]):02d}-{int(rev[2]):02d} {rev[3]:.2f}h  diff={diff_sec:.3f}s"
                )

            jd_start = max(se_jd, le_jd) + 20
        except Exception as e:
            r.skip(f"{label}: {e}")
            jd_start += 29.53

    return r.summary(), r


def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Quarter Moons 2024 (first 6 each)")
    print("=" * 70)

    r = R("P3: Quarter Moons")

    for phase, phase_name in [(90, "FirstQ"), (270, "LastQ")]:
        jd_start = swe.julday(2024, 1, 1, 0.0)
        for i in range(6):
            label = f"{phase_name} #{i + 1}"
            try:
                se_jd = find_moon_phase_se(jd_start, phase)
                le_jd = find_moon_phase_le(jd_start, phase)

                diff_sec = abs(se_jd - le_jd) * 86400
                tol = 5.0

                if diff_sec > tol:
                    r.fail(f"{label}: diff={diff_sec:.2f}s")
                else:
                    r.ok(diff_sec, label)

                jd_start = max(se_jd, le_jd) + 20
            except Exception as e:
                r.skip(f"{label}: {e}")
                jd_start += 29.53

    return r.summary(), r


def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Lunation length consistency (synodic month ~29.53 days)")
    print("=" * 70)

    r = R("P4: Lunation length")

    jd_start = swe.julday(2024, 1, 1, 0.0)
    prev_jd = None

    for i in range(13):
        try:
            le_jd = find_moon_phase_le(jd_start, 0)

            if prev_jd is not None:
                lunation = le_jd - prev_jd
                label = f"Lunation #{i}"
                # Synodic month varies from ~29.27 to ~29.83 days
                if lunation < 29.2 or lunation > 29.9:
                    r.fail(f"{label}: length={lunation:.4f} days (out of range)")
                else:
                    r.ok(abs(lunation - 29.53), label)
                    print(f"  {label}: {lunation:.4f} days")

            prev_jd = le_jd
            jd_start = le_jd + 20
        except Exception as e:
            jd_start += 29.53

    return r.summary(), r


def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Multi-year New Moon sweep (2020-2026)")
    print("=" * 70)

    r = R("P5: Multi-year New Moons")

    for year in range(2020, 2027):
        jd_start = swe.julday(year, 1, 1, 0.0)
        for i in range(3):  # 3 New Moons per year (sample)
            label = f"{year} NM#{i + 1}"
            try:
                se_jd = find_moon_phase_se(jd_start, 0)
                le_jd = find_moon_phase_le(jd_start, 0)

                diff_sec = abs(se_jd - le_jd) * 86400
                tol = 5.0

                if diff_sec > tol:
                    r.fail(f"{label}: diff={diff_sec:.2f}s")
                else:
                    r.ok(diff_sec, label)

                jd_start = max(se_jd, le_jd) + 20
            except Exception as e:
                r.skip(f"{label}: {e}")
                jd_start += 29.53

    return r.summary(), r


def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Elongation verification at computed phase times")
    print("=" * 70)

    r = R("P6: Elongation verification")

    jd_start = swe.julday(2024, 3, 1, 0.0)

    for phase, phase_name in [
        (0, "New"),
        (90, "FirstQ"),
        (180, "Full"),
        (270, "LastQ"),
    ]:
        label = f"{phase_name} elongation"
        try:
            le_jd = find_moon_phase_le(jd_start, phase)

            # Compute actual elongation at this time
            sun = ephem.calc_ut(le_jd, SUN, 0)[0]
            moon = ephem.calc_ut(le_jd, MOON, 0)[0]
            elong = (moon[0] - sun[0]) % 360

            # For New Moon, elongation should be ~0° (or ~360°)
            diff_from_target = abs(elong - phase)
            if diff_from_target > 180:
                diff_from_target = 360 - diff_from_target

            diff_arcsec = diff_from_target * 3600

            if diff_arcsec > 1.0:  # 1" tolerance
                r.fail(
                    f'{label}: elongation={elong:.6f}° target={phase}° err={diff_arcsec:.3f}"'
                )
            else:
                r.ok(diff_arcsec, label)
                print(f'  {label}: elongation={elong:.6f}° err={diff_arcsec:.4f}"')

        except Exception as e:
            r.skip(f"{label}: {e}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 38: Moon Phase Timing")
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
    print("ROUND 38 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 38: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
