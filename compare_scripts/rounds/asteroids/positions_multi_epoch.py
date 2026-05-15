#!/usr/bin/env python3
"""Round 78: Asteroid Positions Deep (Ceres/Pallas/Juno/Vesta)"""

from __future__ import annotations
import sys, os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import swisseph as swe
import libephemeris as ephem

swe.set_ephe_path(_REF_EPHE_PATH)
passed = failed = errors = 0
FLAGS = 256  # FLG_SPEED

# libephemeris uses CERES=17, PALLAS=18, JUNO=19, VESTA=20
# pyswisseph uses AST_OFFSET + N (10000 + N)
LE_BODIES = {17: "Ceres", 18: "Pallas", 19: "Juno", 20: "Vesta"}
SE_BODIES = {10001: "Ceres", 10002: "Pallas", 10003: "Juno", 10004: "Vesta"}

DATES = [2451545.0 + i * 30 for i in range(240)]  # 20 years monthly

print("=" * 70)
print("ROUND 78: Asteroid Positions Deep")
print("=" * 70)

# P1: Longitude comparison
print("\n=== P1: Asteroid longitude (20-year monthly) ===")
for le_id, name in LE_BODIES.items():
    se_id = 10000 + (le_id - 16)  # Map LE IDs to SE AST_OFFSET
    body_pass = body_fail = 0
    for jd in DATES:
        try:
            se = swe.calc_ut(jd, se_id, FLAGS)
            le = ephem.calc_ut(jd, le_id, FLAGS)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 1.0:
                passed += 1
                body_pass += 1
            else:
                failed += 1
                body_fail += 1
                if body_fail <= 2:
                    print(
                        f'  FAIL {name} jd={jd:.0f}: SE={se[0][0]:.6f} LE={le[0][0]:.6f} d={diff * 3600:.2f}"'
                    )
        except Exception as e:
            errors += 1
    t = body_pass + body_fail
    print(
        f"  {name}: {body_pass}/{t} ({100 * body_pass / t:.1f}%)"
        if t > 0
        else f"  {name}: no data"
    )
print(f"  After P1: {passed} passed, {failed} failed, {errors} errors")

# P2: Latitude comparison
print("\n=== P2: Asteroid latitude ===")
for le_id, name in LE_BODIES.items():
    se_id = 10000 + (le_id - 16)
    for jd in DATES[:120]:
        try:
            se = swe.calc_ut(jd, se_id, FLAGS)
            le = ephem.calc_ut(jd, le_id, FLAGS)
            diff = abs(se[0][1] - le[0][1]) * 3600
            if diff < 1.0:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P2: {passed} passed, {failed} failed, {errors} errors")

# P3: Distance comparison
print("\n=== P3: Asteroid distance ===")
for le_id, name in LE_BODIES.items():
    se_id = 10000 + (le_id - 16)
    for jd in DATES[:120]:
        try:
            se = swe.calc_ut(jd, se_id, FLAGS)
            le = ephem.calc_ut(jd, le_id, FLAGS)
            if se[0][2] > 0 and le[0][2] > 0:
                ratio = le[0][2] / se[0][2]
                if 0.9999 < ratio < 1.0001:
                    passed += 1
                else:
                    failed += 1
            else:
                passed += 1
        except:
            errors += 1
print(f"  After P3: {passed} passed, {failed} failed, {errors} errors")

# P4: Speed comparison
print("\n=== P4: Asteroid speed ===")
for le_id, name in LE_BODIES.items():
    se_id = 10000 + (le_id - 16)
    for jd in DATES[:120]:
        try:
            se = swe.calc_ut(jd, se_id, FLAGS)
            le = ephem.calc_ut(jd, le_id, FLAGS)
            diff = abs(se[0][3] - le[0][3])
            if diff < 0.001:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P4: {passed} passed, {failed} failed, {errors} errors")

# P5: J2000 mode
print("\n=== P5: Asteroids in J2000 ===")
J2000 = 256 | 32
for le_id, name in LE_BODIES.items():
    se_id = 10000 + (le_id - 16)
    for jd in DATES[:60]:
        try:
            se = swe.calc_ut(jd, se_id, J2000)
            le = ephem.calc_ut(jd, le_id, J2000)
            diff = abs(se[0][0] - le[0][0])
            if diff > 180:
                diff = 360 - diff
            if diff * 3600 < 1.0:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P5: {passed} passed, {failed} failed, {errors} errors")

# P6: Equatorial mode
print("\n=== P6: Asteroids equatorial ===")
EQ = 256 | 2048
for le_id, name in LE_BODIES.items():
    se_id = 10000 + (le_id - 16)
    for jd in DATES[:60]:
        try:
            se = swe.calc_ut(jd, se_id, EQ)
            le = ephem.calc_ut(jd, le_id, EQ)
            diff_ra = abs(se[0][0] - le[0][0])
            if diff_ra > 180:
                diff_ra = 360 - diff_ra
            diff_dec = abs(se[0][1] - le[0][1])
            if diff_ra * 3600 < 1.0 and diff_dec * 3600 < 1.0:
                passed += 1
            else:
                failed += 1
        except:
            errors += 1
print(f"  After P6: {passed} passed, {failed} failed, {errors} errors")

print(f"\n{'=' * 70}")
total = passed + failed
pct = 100 * passed / total if total > 0 else 0
print(f"ROUND 78 FINAL: {passed}/{total} passed ({pct:.1f}%)")
print(f"  Passed:  {passed}")
print(f"  Failed:  {failed}")
print(f"  Errors:  {errors}")
print(f"{'=' * 70}")
