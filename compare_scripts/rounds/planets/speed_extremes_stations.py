#!/usr/bin/env python3
"""Round 69: Planetary Speed Extremes & Station Detection"""

from __future__ import annotations
import sys, os
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

# Reference ephemeris data path (set via REF_EPHE_PATH env var)
_REF_EPHE_PATH = os.environ.get("REF_EPHE_PATH", "./ephe")

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256  # FLG_SPEED

NAMES = {
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}

print("=" * 70)
print("ROUND 69: Planetary Speed Extremes & Stations")
print("=" * 70)

# P1: All planet speeds daily for 5 years
print("\n=== P1: All planet speeds (5-year daily sweep) ===")
jd0 = 2451545.0
for body in range(2, 10):
    body_pass = body_fail = 0
    for i in range(0, 1825, 1):
        jd = jd0 + i
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff = abs(se[0][3] - le[0][3])
            tol = 0.0001 if body <= 6 else 0.001
            if diff < tol:
                passed += 1
                body_pass += 1
            else:
                failed += 1
                body_fail += 1
        except:
            errors += 1
    if body_fail > 0:
        print(
            f"  {NAMES[body]}: {body_pass}/{body_pass + body_fail} ({100 * body_pass / (body_pass + body_fail):.1f}%)"
        )
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Detect stations (speed near zero) and compare timing
print("\n=== P2: Station detection ===")
for body in [4, 5, 6, 7]:  # Mars-Uranus
    prev_se_speed = None
    stations_found = 0
    for i in range(0, 730):
        jd = jd0 + i
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            se_speed = se[0][3]
            le_speed = le[0][3]
            # Check if station (speed crosses zero)
            if prev_se_speed is not None and prev_se_speed * se_speed < 0:
                stations_found += 1
                # At station, speeds should both be near zero
                if abs(se_speed) < 0.1 and abs(le_speed) < 0.1:
                    passed += 1
                else:
                    failed += 1
                    print(
                        f"  FAIL P2 {NAMES[body]} station at jd={jd:.1f}: SE_spd={se_speed:.6f} LE_spd={le_speed:.6f}"
                    )
            prev_se_speed = se_speed
        except:
            errors += 1
    if stations_found > 0:
        pass  # good
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Latitude speed consistency
print("\n=== P3: Latitude speed ===")
for body in range(0, 10):
    for i in range(0, 365, 5):
        jd = jd0 + i
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff = abs(se[0][4] - le[0][4])  # lat speed
            if diff < 0.001:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Distance speed consistency
print("\n=== P4: Distance speed ===")
for body in range(0, 10):
    for i in range(0, 365, 5):
        jd = jd0 + i
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            diff = abs(se[0][5] - le[0][5])  # dist speed
            if diff < 0.0001:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 69 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
