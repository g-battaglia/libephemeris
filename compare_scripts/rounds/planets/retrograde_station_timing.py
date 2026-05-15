#!/usr/bin/env python3
"""
Round 22: Retrograde Detection & Station Timing Deep Verification
==================================================================

Tests retrograde/direct motion detection and station (stationary point) timing
by comparing planetary speeds and crossing functions between SE and LE.

Parts:
  P1: Speed sign agreement (retrograde detection) — all planets × 20 dates
  P2: Station timing via speed zero-crossing — Mercury through Pluto
  P3: Mercury retrograde periods 2020-2025 (well-known dates)
  P4: Speed magnitude comparison across full orbits
  P5: Outer planet stations — Jupiter through Pluto (slow-moving)
  P6: Venus retrograde (rare, ~40-day events every 19 months)
"""

from __future__ import annotations

import os
import sys
import time
import traceback
import math

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
            print(f"  Max diff: {self.max_diff:.6f} ({self.max_label})")
        if self.failures:
            for f in self.failures[:20]:
                print(f"    - {f}")
            if len(self.failures) > 20:
                print(f"    ... and {len(self.failures) - 20} more")
        print(f"{'=' * 70}")
        return self.failed == 0


PLANETS = [
    (MERCURY, "Mercury"),
    (VENUS, "Venus"),
    (MARS, "Mars"),
    (JUPITER, "Jupiter"),
    (SATURN, "Saturn"),
    (URANUS, "Uranus"),
    (NEPTUNE, "Neptune"),
    (PLUTO, "Pluto"),
]


def get_speed(body, jd):
    """Get ecliptic longitude speed from both SE and LE."""
    try:
        se_xx = swe.calc_ut(jd, body, FLG_SPEED)[0]
        se_speed = se_xx[3]  # lon speed in deg/day
    except Exception:
        se_speed = None

    try:
        le_xx, _ = ephem.calc_ut(jd, body, FLG_SPEED)
        le_speed = le_xx[3]
    except Exception:
        le_speed = None

    return se_speed, le_speed


# ============================================================
# PART 1: Speed sign agreement (retrograde detection)
# ============================================================
def run_part1():
    print("\n" + "=" * 70)
    print("PART 1: Speed sign agreement — all planets × 20 dates (2020-2025)")
    print("=" * 70)

    r = R("P1: Speed Sign")

    # 20 dates spanning 2020-2025
    dates = []
    for y in range(2020, 2026):
        for m in [1, 4, 7, 10]:
            if y == 2025 and m > 1:
                break
            dates.append((y, m, 15, 12.0))

    for body_id, body_name in PLANETS:
        agree = disagree = 0
        for y, m, d, h in dates:
            jd = swe.julday(y, m, d, h)
            se_spd, le_spd = get_speed(body_id, jd)

            if se_spd is None or le_spd is None:
                r.skip(f"{body_name} {y}-{m:02d}-{d:02d}")
                continue

            label = f"{body_name} {y}-{m:02d}-{d:02d}"

            # Check sign agreement
            se_retro = se_spd < 0
            le_retro = le_spd < 0

            if se_retro == le_retro:
                # Also check speed magnitude
                speed_diff = abs(se_spd - le_spd)
                # Relative tolerance: 1% of speed magnitude, min 0.001°/day
                tol = max(abs(se_spd) * 0.01, 0.001)
                if speed_diff > tol:
                    r.fail(
                        f"{label}: speed diff={speed_diff:.6f}°/d (SE={se_spd:.6f} LE={le_spd:.6f})"
                    )
                else:
                    r.ok(speed_diff, label)
                agree += 1
            else:
                # Near-zero speed (station) — both within 0.01°/day of zero
                if abs(se_spd) < 0.01 and abs(le_spd) < 0.01:
                    r.ok(0.0, f"{label} (near-station)")
                else:
                    r.fail(f"{label}: SIGN DISAGREE SE={se_spd:.6f} LE={le_spd:.6f}")
                    disagree += 1

        status = "OK" if disagree == 0 else f"{disagree} DISAGREE"
        print(f"  {body_name:10s}: {agree} agree, {disagree} disagree — {status}")

    return r.summary(), r


# ============================================================
# PART 2: Station timing via speed zero-crossing
# ============================================================
def run_part2():
    print("\n" + "=" * 70)
    print("PART 2: Station timing — find speed=0 crossings (2024)")
    print("=" * 70)

    r = R("P2: Station Timing")

    # For each planet, scan 2024 looking for speed sign changes
    jd_start = swe.julday(2024, 1, 1, 0.0)
    jd_end = swe.julday(2025, 1, 1, 0.0)
    step = 5.0  # 5-day steps for initial scan

    for body_id, body_name in PLANETS:
        # Skip Sun and Moon (no retrograde)
        stations_found = 0

        jd = jd_start
        se_spd_prev, le_spd_prev = get_speed(body_id, jd)

        while jd < jd_end:
            jd += step
            se_spd, le_spd = get_speed(body_id, jd)

            if se_spd is None or le_spd is None:
                se_spd_prev, le_spd_prev = se_spd, le_spd
                continue

            # Check for SE speed sign change
            if se_spd_prev is not None and se_spd_prev * se_spd < 0:
                # Station detected in SE — bisect to find exact time
                se_station = _bisect_speed_zero(body_id, jd - step, jd, use_se=True)
                le_station = _bisect_speed_zero(body_id, jd - step, jd, use_se=False)

                if se_station is not None and le_station is not None:
                    diff_hours = abs(se_station - le_station) * 24.0
                    station_type = "Rx" if se_spd < 0 else "D"
                    label = f"{body_name} {station_type} JD={se_station:.2f}"

                    # Tolerance: 2 hours for inner planets, 6 hours for outer
                    tol_h = 2.0 if body_id <= MARS else 6.0
                    if diff_hours > tol_h:
                        r.fail(f"{label}: diff={diff_hours:.2f}h")
                    else:
                        r.ok(diff_hours, label)

                    stations_found += 1
                    # Convert JD to date for display
                    yr, mo, da, hr = swe.revjul(se_station)
                    print(
                        f"  {body_name:10s} {station_type}: {yr}-{mo:02d}-{int(da):02d} "
                        f"{hr:.1f}h — diff={diff_hours:.2f}h"
                    )

            se_spd_prev = se_spd
            le_spd_prev = le_spd

        if stations_found == 0:
            print(f"  {body_name:10s}: no stations in 2024")

    return r.summary(), r


def _bisect_speed_zero(body_id, jd_a, jd_b, use_se=True, max_iter=50):
    """Bisect to find speed=0 crossing."""
    for _ in range(max_iter):
        jd_mid = (jd_a + jd_b) / 2
        if use_se:
            try:
                xx = swe.calc_ut(jd_mid, body_id, FLG_SPEED)[0]
                spd = xx[3]
            except Exception:
                return None
        else:
            try:
                xx, _ = ephem.calc_ut(jd_mid, body_id, FLG_SPEED)
                spd = xx[3]
            except Exception:
                return None

        # Get speed at jd_a for comparison
        if use_se:
            try:
                xa = swe.calc_ut(jd_a, body_id, FLG_SPEED)[0]
                spd_a = xa[3]
            except Exception:
                return None
        else:
            try:
                xa, _ = ephem.calc_ut(jd_a, body_id, FLG_SPEED)
                spd_a = xa[3]
            except Exception:
                return None

        if spd_a * spd < 0:
            jd_b = jd_mid
        else:
            jd_a = jd_mid

        if abs(jd_b - jd_a) < 1.0 / 86400.0:  # 1 second precision
            return (jd_a + jd_b) / 2

    return (jd_a + jd_b) / 2


# ============================================================
# PART 3: Mercury retrograde periods 2020-2025
# ============================================================
def run_part3():
    print("\n" + "=" * 70)
    print("PART 3: Mercury retrograde periods 2020-2025")
    print("=" * 70)

    r = R("P3: Mercury Retro")

    # Known Mercury retrograde approximate start dates (Rx station)
    # Source: widely published astrology/astronomy tables
    mercury_rx_dates = [
        # (year, month, day) — approximate Rx station dates
        (2020, 2, 16),
        (2020, 6, 18),
        (2020, 10, 14),
        (2021, 1, 30),
        (2021, 5, 29),
        (2021, 9, 27),
        (2022, 1, 14),
        (2022, 5, 10),
        (2022, 9, 10),
        (2022, 12, 29),
        (2023, 4, 21),
        (2023, 8, 23),
        (2023, 12, 13),
        (2024, 4, 1),
        (2024, 8, 5),
        (2024, 11, 26),
    ]

    for y, m, d in mercury_rx_dates:
        # Search ±5 days around the expected date
        jd_center = swe.julday(y, m, d, 12.0)
        jd_search_start = jd_center - 5.0
        jd_search_end = jd_center + 5.0

        se_station = None
        le_station = None

        # Scan for speed sign change (positive → negative = Rx)
        step = 0.25  # 6-hour steps
        jd = jd_search_start
        se_prev = swe.calc_ut(jd, MERCURY, FLG_SPEED)[0][3]
        le_prev_xx, _ = ephem.calc_ut(jd, MERCURY, FLG_SPEED)
        le_prev = le_prev_xx[3]

        while jd < jd_search_end:
            jd += step
            se_spd = swe.calc_ut(jd, MERCURY, FLG_SPEED)[0][3]
            le_xx, _ = ephem.calc_ut(jd, MERCURY, FLG_SPEED)
            le_spd = le_xx[3]

            if se_prev > 0 and se_spd <= 0 and se_station is None:
                se_station = _bisect_speed_zero(MERCURY, jd - step, jd, use_se=True)

            if le_prev > 0 and le_spd <= 0 and le_station is None:
                le_station = _bisect_speed_zero(MERCURY, jd - step, jd, use_se=False)

            if se_station is not None and le_station is not None:
                break

            se_prev = se_spd
            le_prev = le_spd

        label = f"Mercury Rx {y}-{m:02d}-{d:02d}"

        if se_station is None or le_station is None:
            r.fail(f"{label}: not found (SE={se_station} LE={le_station})")
            continue

        diff_hours = abs(se_station - le_station) * 24.0

        if diff_hours > 2.0:
            r.fail(f"{label}: diff={diff_hours:.2f}h")
        else:
            r.ok(diff_hours, label)

        yr, mo, da, hr = swe.revjul(se_station)
        print(
            f"  {label}: found {yr}-{mo:02d}-{int(da):02d} {hr:.1f}h, diff={diff_hours:.3f}h"
        )

    return r.summary(), r


# ============================================================
# PART 4: Speed magnitude comparison across full orbits
# ============================================================
def run_part4():
    print("\n" + "=" * 70)
    print("PART 4: Speed magnitude — full 2024 sweep, weekly samples")
    print("=" * 70)

    r = R("P4: Speed Magnitude")

    jd_start = swe.julday(2024, 1, 1, 0.0)
    n_weeks = 52

    for body_id, body_name in PLANETS:
        max_diff = 0.0
        max_pct = 0.0

        for w in range(n_weeks):
            jd = jd_start + w * 7.0
            se_spd, le_spd = get_speed(body_id, jd)

            if se_spd is None or le_spd is None:
                r.skip(f"{body_name} week {w}")
                continue

            diff = abs(se_spd - le_spd)
            label = f"{body_name} week {w}"

            # Speed tolerance: 0.5% of speed + 0.0001°/day baseline
            tol = abs(se_spd) * 0.005 + 0.0001
            if diff > tol:
                r.fail(f"{label}: SE={se_spd:.6f} LE={le_spd:.6f} diff={diff:.6f}")
            else:
                r.ok(diff, label)

            if diff > max_diff:
                max_diff = diff
            if abs(se_spd) > 0.001:
                pct = diff / abs(se_spd) * 100
                if pct > max_pct:
                    max_pct = pct

        print(f"  {body_name:10s}: max_diff={max_diff:.6f}°/d ({max_pct:.4f}%)")

    return r.summary(), r


# ============================================================
# PART 5: Outer planet stations — Jupiter through Pluto
# ============================================================
def run_part5():
    print("\n" + "=" * 70)
    print("PART 5: Outer planet stations 2020-2025")
    print("=" * 70)

    r = R("P5: Outer Stations")

    outer_planets = [
        (JUPITER, "Jupiter"),
        (SATURN, "Saturn"),
        (URANUS, "Uranus"),
        (NEPTUNE, "Neptune"),
        (PLUTO, "Pluto"),
    ]

    jd_start = swe.julday(2020, 1, 1, 0.0)
    jd_end = swe.julday(2026, 1, 1, 0.0)
    step = 10.0  # 10-day steps for outer planets

    for body_id, body_name in outer_planets:
        stations_found = 0
        jd = jd_start
        se_spd_prev, _ = get_speed(body_id, jd)

        while jd < jd_end:
            jd += step
            se_spd, le_spd = get_speed(body_id, jd)

            if se_spd is None or se_spd_prev is None:
                se_spd_prev = se_spd
                continue

            if se_spd_prev * se_spd < 0:
                se_station = _bisect_speed_zero(body_id, jd - step, jd, use_se=True)
                le_station = _bisect_speed_zero(body_id, jd - step, jd, use_se=False)

                if se_station is not None and le_station is not None:
                    diff_hours = abs(se_station - le_station) * 24.0
                    station_type = "Rx" if se_spd < 0 else "D"
                    label = f"{body_name} {station_type}"

                    # Outer planets: tolerance 6 hours
                    if diff_hours > 6.0:
                        r.fail(f"{label}: diff={diff_hours:.2f}h")
                    else:
                        r.ok(diff_hours, label)

                    stations_found += 1
                    yr, mo, da, hr = swe.revjul(se_station)
                    print(
                        f"  {body_name:10s} {station_type}: {yr}-{mo:02d}-{int(da):02d} "
                        f"diff={diff_hours:.3f}h"
                    )

            se_spd_prev = se_spd

        if stations_found == 0:
            print(f"  {body_name:10s}: no stations found")

    return r.summary(), r


# ============================================================
# PART 6: Venus retrograde
# ============================================================
def run_part6():
    print("\n" + "=" * 70)
    print("PART 6: Venus retrograde events 2020-2025")
    print("=" * 70)

    r = R("P6: Venus Retro")

    jd_start = swe.julday(2020, 1, 1, 0.0)
    jd_end = swe.julday(2026, 1, 1, 0.0)
    step = 5.0

    stations_found = 0
    jd = jd_start
    se_spd_prev, _ = get_speed(VENUS, jd)

    while jd < jd_end:
        jd += step
        se_spd, le_spd = get_speed(VENUS, jd)

        if se_spd is None or se_spd_prev is None:
            se_spd_prev = se_spd
            continue

        if se_spd_prev * se_spd < 0:
            se_station = _bisect_speed_zero(VENUS, jd - step, jd, use_se=True)
            le_station = _bisect_speed_zero(VENUS, jd - step, jd, use_se=False)

            if se_station is not None and le_station is not None:
                diff_hours = abs(se_station - le_station) * 24.0
                station_type = "Rx" if se_spd < 0 else "D"
                label = f"Venus {station_type}"

                if diff_hours > 2.0:
                    r.fail(f"{label}: diff={diff_hours:.2f}h")
                else:
                    r.ok(diff_hours, label)

                stations_found += 1
                yr, mo, da, hr = swe.revjul(se_station)
                print(
                    f"  Venus {station_type}: {yr}-{mo:02d}-{int(da):02d} "
                    f"diff={diff_hours:.3f}h"
                )

        se_spd_prev = se_spd

    if stations_found == 0:
        print("  No Venus stations found")

    return r.summary(), r


def main():
    print("=" * 70)
    print("ROUND 22: Retrograde Detection & Station Timing")
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
    print("ROUND 22 FINAL SUMMARY")
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
    if all_ok:
        print(f"\n  >>> ROUND 22: ALL PASSED <<<")
    else:
        print(f"\n  >>> ROUND 22: {tf} FAILURES <<<")

    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
