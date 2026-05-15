#!/usr/bin/env python3
"""Round 68: Moon Distance (Perigee/Apogee Cycle) + Sun Distance"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256  # FLG_SPEED

print("=" * 70)
print("ROUND 68: Moon & Sun Distance Precision")
print("=" * 70)

# P1: Moon distance every 6 hours for 2 years (covers ~26 perigee/apogee cycles)
print("\n=== P1: Moon distance (2-year sweep, 6h steps) ===")
jd_start = 2451545.0
for i in range(2920):  # ~2 years at 6h steps
    jd = jd_start + i * 0.25
    try:
        se = swe.calc_ut(jd, 1, FLAGS)
        le = ephem.calc_ut(jd, 1, FLAGS)
        se_dist = se[0][2]
        le_dist = le[0][2]
        diff_km = abs(se_dist - le_dist) * 149597870.7  # AU to km
        if diff_km < 0.5:  # sub-500m
            passed += 1
        else:
            failed += 1
            if failed <= 10:
                print(
                    f"  FAIL P1 jd={jd:.2f}: SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_km:.1f}km"
                )
    except Exception as e:
        errors += 1
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Sun distance every day for 2 years
print("\n=== P2: Sun distance (2-year sweep, daily) ===")
for i in range(730):
    jd = jd_start + i
    try:
        se = swe.calc_ut(jd, 0, FLAGS)
        le = ephem.calc_ut(jd, 0, FLAGS)
        se_dist = se[0][2]
        le_dist = le[0][2]
        diff_km = abs(se_dist - le_dist) * 149597870.7
        if diff_km < 50:  # sub-50km
            passed += 1
        else:
            failed += 1
            if failed <= 10:
                print(
                    f"  FAIL P2 jd={jd:.1f}: SE={se_dist:.10f} LE={le_dist:.10f} diff={diff_km:.1f}km"
                )
    except Exception as e:
        errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Moon speed (deg/day) — correlates with distance
print("\n=== P3: Moon speed consistency ===")
for i in range(365):
    jd = jd_start + i
    try:
        se = swe.calc_ut(jd, 1, FLAGS)
        le = ephem.calc_ut(jd, 1, FLAGS)
        se_speed = se[0][3]  # lon speed
        le_speed = le[0][3]
        diff = abs(se_speed - le_speed)
        if diff < 0.001:  # 0.001 deg/day = 3.6 arcsec/day
            passed += 1
        else:
            failed += 1
            if failed <= 10:
                print(
                    f"  FAIL P3 jd={jd:.1f}: SE={se_speed:.6f} LE={le_speed:.6f} diff={diff:.6f}"
                )
    except Exception as e:
        errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Moon at known perigee/apogee dates (2020)
print("\n=== P4: Moon at known perigee/apogee 2020 ===")
# Known perigees and apogees of 2020 (approximate JD)
perigee_dates = [
    2458854.5,  # 2020-01-13
    2458882.0,  # 2020-02-10
    2458908.5,  # 2020-03-10
    2458936.5,  # 2020-04-07
    2458965.0,  # 2020-05-06
]
apogee_dates = [
    2458868.0,  # 2020-01-27
    2458895.5,  # 2020-02-24
    2458923.0,  # 2020-03-24
    2458950.0,  # 2020-04-20
    2458978.0,  # 2020-05-18
]
for jd in perigee_dates + apogee_dates:
    try:
        se = swe.calc_ut(jd, 1, FLAGS)
        le = ephem.calc_ut(jd, 1, FLAGS)
        diff_km = abs(se[0][2] - le[0][2]) * 149597870.7
        if diff_km < 0.5:
            passed += 1
        else:
            failed += 1
            print(f"  FAIL P4 jd={jd:.1f}: diff={diff_km:.1f}km")
    except Exception as e:
        errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Planet distances (Mars, Jupiter, Saturn) over 5 years
print("\n=== P5: Outer planet distances ===")
for body in [4, 5, 6]:  # Mars, Jupiter, Saturn
    names = {4: "Mars", 5: "Jupiter", 6: "Saturn"}
    for i in range(0, 1825, 10):  # 5 years, every 10 days
        jd = jd_start + i
        try:
            se = swe.calc_ut(jd, body, FLAGS)
            le = ephem.calc_ut(jd, body, FLAGS)
            se_dist = se[0][2]
            le_dist = le[0][2]
            if se_dist > 0:
                ratio = le_dist / se_dist
                if 0.99999 < ratio < 1.00001:
                    passed += 1
                else:
                    failed += 1
                    if failed <= 15:
                        print(f"  FAIL P5 {names[body]} jd={jd:.0f}: ratio={ratio:.8f}")
        except Exception as e:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 68 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
