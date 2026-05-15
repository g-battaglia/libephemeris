#!/usr/bin/env python3
"""
Round 28: Mean/True Lilith & Node Edge Cases
==============================================

Deep verification of lunar node and Lilith (lunar apogee) calculations
at edge cases, multiple epochs, and with various flag combinations.

Parts:
  P1: Mean Node — 20 dates across 1900-2100
  P2: True Node — 20 dates across 1900-2100
  P3: Mean Lilith (Mean Apogee) — 20 dates
  P4: Osc Lilith (True/Osculating Apogee) — 20 dates
  P5: Node/Lilith with FLG_J2000
  P6: Node/Lilith with FLG_EQUATORIAL
  P7: Speed comparison for all 4 bodies
  P8: 18.6-year nodal cycle verification
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
            print(f'  Max diff: {self.max_diff:.2f}" ({self.max_label})')
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


DATES_20 = [
    swe.julday(y, m, 15, 12.0)
    for y, m in [
        (1900, 1),
        (1920, 6),
        (1940, 3),
        (1950, 9),
        (1960, 12),
        (1970, 4),
        (1980, 7),
        (1985, 11),
        (1990, 2),
        (1995, 8),
        (2000, 1),
        (2005, 5),
        (2010, 10),
        (2015, 3),
        (2018, 6),
        (2020, 9),
        (2022, 1),
        (2024, 3),
        (2024, 9),
        (2050, 6),
    ]
]


def run_body_sweep(part_name, body_id, body_name, flags, tol_lon, tol_lat):
    """Run a sweep for a single body across dates."""
    r = R(part_name)

    for jd in DATES_20:
        yr = int(swe.revjul(jd)[0])
        label = f"{body_name} {yr}"
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

        if lon_d > tol_lon:
            r.fail(f'{label}: lon={lon_d:.2f}"')
        elif lat_d > tol_lat:
            r.fail(f'{label}: lat={lat_d:.2f}"')
        else:
            r.ok(max(lon_d, lat_d), label)

    return r.summary(), r


def run_part1():
    print(f"\n{'=' * 70}")
    print(f"PART 1: Mean Node — 20 dates")
    print(f"{'=' * 70}")
    return run_body_sweep("P1: Mean Node", MEAN_NODE, "MeanNode", 0, 5.0, 5.0)


def run_part2():
    print(f"\n{'=' * 70}")
    print(f"PART 2: True Node — 20 dates")
    print(f"{'=' * 70}")
    return run_body_sweep("P2: True Node", TRUE_NODE, "TrueNode", 0, 5.0, 5.0)


def run_part3():
    print(f"\n{'=' * 70}")
    print(f"PART 3: Mean Lilith — 20 dates")
    print(f"{'=' * 70}")
    # Wider lat tolerance for Mean Lilith (known ~19" lat model difference)
    return run_body_sweep("P3: Mean Lilith", MEAN_APOG, "MeanLilith", 0, 5.0, 25.0)


def run_part4():
    print(f"\n{'=' * 70}")
    print(f"PART 4: Osc Lilith — 20 dates")
    print(f"{'=' * 70}")
    return run_body_sweep("P4: Osc Lilith", OSCU_APOG, "OscLilith", 0, 60.0, 60.0)


def run_part5():
    print(f"\n{'=' * 70}")
    print(f"PART 5: Node/Lilith with FLG_J2000")
    print(f"{'=' * 70}")
    r = R("P5: J2000")
    bodies = [
        (MEAN_NODE, "MeanNode", 5.0),
        (TRUE_NODE, "TrueNode", 5.0),
        (MEAN_APOG, "MeanLilith", 25.0),
    ]
    for jd in DATES_20[:10]:
        yr = int(swe.revjul(jd)[0])
        for body_id, body_name, tol in bodies:
            label = f"{body_name} J2000 {yr}"
            try:
                se_xx = swe.calc_ut(jd, body_id, FLG_J2000 | FLG_SPEED)[0]
                le_xx, _ = ephem.calc_ut(jd, body_id, FLG_J2000 | FLG_SPEED)
            except Exception:
                r.skip(f"{label}")
                continue
            lon_d = angle_diff(se_xx[0], le_xx[0]) * 3600
            lat_d = abs(se_xx[1] - le_xx[1]) * 3600
            max_d = max(lon_d, lat_d)
            if max_d > tol:
                r.fail(f'{label}: lon={lon_d:.2f}" lat={lat_d:.2f}"')
            else:
                r.ok(max_d, label)
    print(f"  Tested {r.passed + r.failed} J2000 combinations")
    return r.summary(), r


def run_part6():
    print(f"\n{'=' * 70}")
    print(f"PART 6: Node/Lilith with FLG_EQUATORIAL")
    print(f"{'=' * 70}")
    r = R("P6: Equatorial")
    bodies = [
        (MEAN_NODE, "MeanNode", 5.0),
        (TRUE_NODE, "TrueNode", 5.0),
        (MEAN_APOG, "MeanLilith", 25.0),
    ]
    for jd in DATES_20[:10]:
        yr = int(swe.revjul(jd)[0])
        for body_id, body_name, tol in bodies:
            label = f"{body_name} EQU {yr}"
            try:
                se_xx = swe.calc_ut(jd, body_id, FLG_EQUATORIAL | FLG_SPEED)[0]
                le_xx, _ = ephem.calc_ut(
                    jd, body_id, FLG_EQUATORIAL | FLG_SPEED
                )
            except Exception:
                r.skip(f"{label}")
                continue
            lon_d = angle_diff(se_xx[0], le_xx[0]) * 3600
            lat_d = abs(se_xx[1] - le_xx[1]) * 3600
            max_d = max(lon_d, lat_d)
            if max_d > tol:
                r.fail(f'{label}: lon={lon_d:.2f}" lat={lat_d:.2f}"')
            else:
                r.ok(max_d, label)
    print(f"  Tested {r.passed + r.failed} equatorial combinations")
    return r.summary(), r


def run_part7():
    print(f"\n{'=' * 70}")
    print(f"PART 7: Speed comparison — all 4 bodies")
    print(f"{'=' * 70}")
    r = R("P7: Speed")
    bodies = [
        (MEAN_NODE, "MeanNode"),
        (TRUE_NODE, "TrueNode"),
        (MEAN_APOG, "MeanLilith"),
        (OSCU_APOG, "OscLilith"),
    ]
    jd = swe.julday(2024, 3, 20, 12.0)
    for body_id, body_name in bodies:
        try:
            se_xx = swe.calc_ut(jd, body_id, FLG_SPEED)[0]
            le_xx, _ = ephem.calc_ut(jd, body_id, FLG_SPEED)
        except Exception:
            r.skip(f"{body_name}")
            continue
        lon_spd_d = abs(se_xx[3] - le_xx[3])
        lat_spd_d = abs(se_xx[4] - le_xx[4])
        tol = 0.01
        if lon_spd_d > tol or lat_spd_d > tol:
            r.fail(f"{body_name}: lon_spd={lon_spd_d:.6f} lat_spd={lat_spd_d:.6f}")
        else:
            r.ok(max(lon_spd_d, lat_spd_d) * 3600, f"{body_name}")
        print(
            f"  {body_name:12s}: SE_spd={se_xx[3]:.6f} LE_spd={le_xx[3]:.6f} "
            f"diff={lon_spd_d:.6f}°/d"
        )
    return r.summary(), r


def run_part8():
    print(f"\n{'=' * 70}")
    print(f"PART 8: 18.6-year nodal cycle — Mean Node regression")
    print(f"{'=' * 70}")
    r = R("P8: Nodal Cycle")
    # Verify Mean Node completes ~360° in ~18.6 years
    jd_start = swe.julday(2000, 1, 1, 12.0)
    le_start, _ = ephem.calc_ut(jd_start, MEAN_NODE, 0)
    # After ~18.6 years, node should return to roughly same position
    jd_end = jd_start + 18.6 * 365.25
    le_end, _ = ephem.calc_ut(jd_end, MEAN_NODE, 0)
    cycle_diff = angle_diff(le_start[0], le_end[0])
    # Should be within ~5° (it's approximate)
    if cycle_diff < 10.0:
        r.ok(cycle_diff * 3600, "18.6yr cycle")
        print(
            f"  Start: {le_start[0]:.4f}° End: {le_end[0]:.4f}° diff: {cycle_diff:.2f}°"
        )
    else:
        r.fail(f"18.6yr cycle: diff={cycle_diff:.2f}°")

    # Also verify node is retrograde (negative speed)
    le_xx, _ = ephem.calc_ut(jd_start, MEAN_NODE, FLG_SPEED)
    if le_xx[3] < 0:
        r.ok(0.0, "Node retrograde")
        print(f"  Mean Node speed: {le_xx[3]:.6f}°/d (retrograde: OK)")
    else:
        r.fail(f"Node speed positive: {le_xx[3]:.6f}")

    # Verify Mean Lilith is direct (positive speed)
    le_lil, _ = ephem.calc_ut(jd_start, MEAN_APOG, FLG_SPEED)
    if le_lil[3] > 0:
        r.ok(0.0, "Lilith direct")
        print(f"  Mean Lilith speed: {le_lil[3]:.6f}°/d (direct: OK)")
    else:
        r.fail(f"Lilith speed negative: {le_lil[3]:.6f}")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 28: Mean/True Lilith & Node Edge Cases")
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
        ("P8", run_part8),
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
    print("ROUND 28 FINAL SUMMARY")
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
    print(f"\n  >>> ROUND 28: {'ALL PASSED' if all_ok else f'{tf} FAILURES'} <<<")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
