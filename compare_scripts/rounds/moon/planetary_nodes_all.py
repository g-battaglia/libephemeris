#!/usr/bin/env python3
"""Round 74: Planetary Nodes All Planets"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256  # FLG_SPEED

NAMES = {
    0: "Sun",
    1: "Moon",
    2: "Mercury",
    3: "Venus",
    4: "Mars",
    5: "Jupiter",
    6: "Saturn",
    7: "Uranus",
    8: "Neptune",
    9: "Pluto",
}
# TRUE_NODE=11, MEAN_NODE=10, MEAN_APOG=12, OSCU_APOG=13
DATES = [2451545.0 + i * 365.25 for i in range(25)]

print("=" * 70)
print("ROUND 74: Planetary Nodes All Planets")
print("=" * 70)

# P1: True Node longitude
print("\n=== P1: True Node ===")
for jd in DATES:
    try:
        se = swe.calc_ut(jd, 11, FLAGS)
        le = ephem.calc_ut(jd, 11, FLAGS)
        diff = abs(se[0][0] - le[0][0])
        if diff > 180:
            diff = 360 - diff
        if diff * 3600 < 1.0:
            passed += 1
        else:
            failed += 1
            print(
                f'  FAIL P1 jd={jd:.0f}: SE={se[0][0]:.6f} LE={le[0][0]:.6f} d={diff * 3600:.2f}"'
            )
    except Exception as e:
        errors += 1
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Mean Node longitude
print("\n=== P2: Mean Node ===")
for jd in DATES:
    try:
        se = swe.calc_ut(jd, 10, FLAGS)
        le = ephem.calc_ut(jd, 10, FLAGS)
        diff = abs(se[0][0] - le[0][0])
        if diff > 180:
            diff = 360 - diff
        if diff * 3600 < 1.0:
            passed += 1
        else:
            failed += 1
    except:
        errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Mean Apogee (Mean Lilith)
print("\n=== P3: Mean Apogee (Lilith) ===")
for jd in DATES:
    try:
        se = swe.calc_ut(jd, 12, FLAGS)
        le = ephem.calc_ut(jd, 12, FLAGS)
        diff = abs(se[0][0] - le[0][0])
        if diff > 180:
            diff = 360 - diff
        if diff * 3600 < 1.0:
            passed += 1
        else:
            failed += 1
    except:
        errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Osculating Apogee (True Lilith)
print("\n=== P4: Osculating Apogee (True Lilith) ===")
for jd in DATES:
    try:
        se = swe.calc_ut(jd, 13, FLAGS)
        le = ephem.calc_ut(jd, 13, FLAGS)
        diff = abs(se[0][0] - le[0][0])
        if diff > 180:
            diff = 360 - diff
        if diff * 3600 < 30.0:  # Wider tolerance for osculating
            passed += 1
        else:
            failed += 1
            if failed <= 5:
                print(
                    f'  FAIL P4 jd={jd:.0f}: SE={se[0][0]:.4f} LE={le[0][0]:.4f} d={diff * 3600:.1f}"'
                )
    except:
        errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: Node/Apsides via nod_aps_ut for all planets
print("\n=== P5: nod_aps_ut for planets ===")
for body in range(0, 10):
    for jd in DATES[:10]:
        try:
            se = swe.nod_aps_ut(jd, body, FLAGS, 0)  # method 0 = osculating
            le = ephem.nod_aps_ut(jd, body, FLAGS, 0)
            # Compare ascending node longitude
            se_asc = se[0][0]
            le_asc = le[0][0]
            diff = abs(se_asc - le_asc)
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 5.0:
                passed += 1
            else:
                failed += 1
                if failed <= 5:
                    print(
                        f'  FAIL P5 {NAMES[body]} jd={jd:.0f}: SE_asc={se_asc:.4f} LE_asc={le_asc:.4f} d={diff * 3600:.1f}"'
                    )
        except Exception as e:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# P6: Mean elements nod_aps_ut (method 1)
print("\n=== P6: nod_aps_ut mean elements ===")
for body in [1, 2, 3, 4, 5, 6]:
    for jd in DATES[:10]:
        try:
            se = swe.nod_aps_ut(jd, body, FLAGS, 1)  # method 1 = mean
            le = ephem.nod_aps_ut(jd, body, FLAGS, 1)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 5.0:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

# P7: True Node at dense time steps (1 year daily)
print("\n=== P7: True Node dense sweep ===")
for i in range(365):
    jd = 2451545.0 + i
    try:
        se = swe.calc_ut(jd, 11, FLAGS)
        le = ephem.calc_ut(jd, 11, FLAGS)
        diff = abs(se[0][0] - le[0][0])
        if diff > 180:
            diff = 360 - diff
        if diff * 3600 < 1.0:
            passed += 1
        else:
            failed += 1
    except:
        errors += 1
print(f"  After P7: {passed} passed, {failed} failed, {errors} errors")

# P8: Node speed comparison
print("\n=== P8: True Node speed ===")
for jd in DATES:
    try:
        se = swe.calc_ut(jd, 11, FLAGS)
        le = ephem.calc_ut(jd, 11, FLAGS)
        diff = abs(se[0][3] - le[0][3])
        if diff < 0.001:
            passed += 1
        else:
            failed += 1
    except:
        errors += 1
print(f"  After P8: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 74 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
